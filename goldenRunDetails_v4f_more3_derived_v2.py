# goldenRunDetails_v4f_more3_derived.py
# Phase-1 "derived-only" runner with optional params JSON, summary export, and clean prints.

import time, json, numpy as np
from math import sqrt
from scipy.integrate import solve_bvp
import math

import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from goldenRunDetails_v4f_more3 import (
    initial_guess_tanh,
    bounce_odes_rescaled,
    boundary_conditions_rescaled,
    bounce_odes_jac,
    export_background_csv,
    analyze_solution,
    wall_radius_from_solution,
    calculate_action,
)

# ==== Your existing imports/utilities go here (unchanged) ====
# from your_module import (
#     initial_guess_tanh, bounce_odes_rescaled, boundary_conditions_rescaled,
#     bounce_odes_jac, export_background_csv, analyze_solution,
#     wall_radius_from_solution, calculate_action
# )
# Keep the names you already use.

# ---------- Helpers we add ----------
def enforce_manifold_consistency(P):
    """
    Make (phi_FV, bias_t, g_portal_t) consistent with (lambda_t, mu2_t, v_phi_t, lam_phi),
    using a calibrated m_Phi if provided, or calibrating it once from the JSON's g_portal_t.
    Returns a modified copy of P.
    """
    P = dict(P)  # copy

    v_t     = float(P["v_phi_t"])
    lam_phi = float(P["lam_phi"])
    lam_t   = float(P["lambda_t"])
    mu2_t   = float(P["mu2_t"])

    # 1) Infer phi_FV from (lambda_t, mu2_t)
    phi2  = 1.0 + 2.0 * (mu2_t / lam_t)
    phiFV = math.sqrt(phi2)

    # 2) SFV-anchored bias (general form; reduces to -(9/4) lam_t/v_t^2 when phi2=3)
    bias = - (phi2**2) / 4.0 * lam_t / (v_t**2)

    # 3) Choose m_Phi:
    #    - Prefer a calibrated m_Phi in JSON
    #    - Else calibrate it ONCE from the JSON's g_portal_t (if present), then store it
    if "m_Phi_used" in P:
        mPhi = float(P["m_Phi_used"])
    else:
        g_json = float(P.get("g_portal_t", 0.0))
        mPhi2  = 2.0*phi2*g_json - lam_phi*(v_t**2) + 2.0*bias
        mPhi   = math.sqrt(max(mPhi2, 0.0))
        P["m_Phi_used"] = mPhi

    # 4) Recompute the required g along the manifold
    g_req = (mPhi**2 + lam_phi*(v_t**2) - 2.0*bias) / (2.0*phi2)

    # 5) Save back (force consistency)
    P["phi_FV"]      = -phiFV     # negative branch (matches your BC convention)
    P["bias_t"]      = bias
    P["g_portal_t"]  = g_req

    # Optional: quick sanity print for humans
    P["_manifold_debug"] = {
        "phiFV_from_mu2_over_lambda": phiFV,
        "bias_rederived": bias,
        "g_required": g_req,
        "m_Phi_used": mPhi,
    }
    return P

def load_params_json(path):
    with open(path, "r") as f:
        P = json.load(f)
    required = ["v_phi_t","lam_phi","lambda_t","mu2_t","bias_t","g_portal_t"]
    missing = [k for k in required if k not in P]
    if missing:
        raise ValueError(f"Missing required keys in params JSON: {missing}")
    return P

def export_run_summary(summary_path, params, g_hist, action_hist, R_star, S_E, elapsed_s):
    data = {
        "params": params,
        "g_history": g_hist,
        "S_E_history": action_hist,
        "R_star": R_star,
        "S_E_final": S_E,
        "elapsed_s": elapsed_s,
    }
    with open(summary_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"[save] summary -> {summary_path}")

# ---------- Main solver ----------

def run_simulation(v_phys, v_phi_phys, lam_phi_phys, lam_phys, g_portal_max, R0, w,
                   derived_params=None, out_prefix=""):
    start_time = time.time()

    print("\n" + "="*60)
    print(f"RUNNING: v_phi={v_phi_phys:.2e}, lam_phi={lam_phi_phys:.3f}, lam={lam_phys:.6f}, g_max={g_portal_max:.2f}")
    print("="*60)

    # ===== Derived parameter set =====
    if derived_params is None:
        # Built-in Option A defaults (if no JSON provided)
        v_phi_t_derived = (v_phi_phys / v_phys)
        lam_phi_derived = 0.1
        lambda_t_derived = 2.3382436740743693e-4
        mu2_t_derived    = 2.3382436740743693e-4
        bias_t_derived   = -1.1457394002964407e-4
        g_portal_t_derived = 1.972261917118631
    else:
        v_phi_t_derived    = float(derived_params["v_phi_t"])
        lam_phi_derived    = float(derived_params["lam_phi"])
        lambda_t_derived   = float(derived_params["lambda_t"])
        mu2_t_derived      = float(derived_params["mu2_t"])
        bias_t_derived     = float(derived_params["bias_t"])
        g_portal_t_derived = float(derived_params["g_portal_t"])
        # For a true derived run, target g_max := g_portal_t exactly
        g_portal_max = g_portal_t_derived

    params_tilde = {
        'v_phi_t':  v_phi_t_derived,
        'lam_phi':  lam_phi_derived,
        'lam_t':    lambda_t_derived,
        'mu2_t':    mu2_t_derived,
        'bias_t':   bias_t_derived,
        'g_portal_t': g_portal_t_derived,
    }

    # FV/TV values: phi_FV will be ~ +sqrt(3) here; BC flips it to the negative branch.
    params_tilde['Phi_FV_t'] = 0.0
    params_tilde['phi_FV_t'] = float(np.sqrt(1.0 + 2.0*params_tilde['mu2_t']/params_tilde['lam_t']))
    params_tilde['Phi_TV_t'] = params_tilde['v_phi_t']
    params_tilde['phi_TV_t'] = 0.0

    print("[derived] v_phi_t={:.6f}, lam_phi={:.3f}, lam_t={:.6e}, mu2_t={:.6e}, bias_t={:.6e}, g={:.6f}, phi_FV≈{:.6f}"
          .format(params_tilde['v_phi_t'], params_tilde['lam_phi'], params_tilde['lam_t'],
                  params_tilde['mu2_t'], params_tilde['bias_t'], params_tilde['g_portal_t'],
                  float(params_tilde['phi_FV_t'])))

    # ===== Build mesh and initial decoupled guess =====
    r_max = 100.0
    mesh = np.unique(np.concatenate([
        np.linspace(1e-9, max(1e-9, R0 - 6*w), 300),
        R0 + w*np.tanh(np.linspace(-4, 4, 2500)),
        np.linspace(R0 + 6*w, r_max, 300),
    ]))

    params_decoupled = params_tilde.copy()
    params_decoupled['g_portal_t'] = 0.0
    guess = initial_guess_tanh(mesh, params_decoupled, R0, w)

    solution = solve_bvp(
        lambda r, y: bounce_odes_rescaled(r, y, params_decoupled),
        lambda ya, yb: boundary_conditions_rescaled(ya, yb, params_decoupled),
        mesh, guess,
        fun_jac=lambda r, y: bounce_odes_jac(r, y, params_decoupled),
        verbose=1, max_nodes=500000, tol=1e-5,
    )
    if not solution.success:
        print("!!! FAILURE: Decoupled solution (g=0) not viable for these parameters !!!")
        return False, None, None, None, None

    print("Decoupled solution OK. Starting adaptive continuation...")

    # Export the g=0 background (optional)
    export_background_csv(solution, out_path=(out_prefix + "background_profile_g0.csv"))

    # Diagnostics at g=0
    R0_calc = wall_radius_from_solution(params_decoupled, solution)
    S0_calc = calculate_action(params_decoupled, solution)
    print(f"[cont] OK   g=0.000  rmax={r_max:.1f}  R≈{R0_calc:.3f}  S_E≈{S0_calc:.6g}")

    g_values = [0.0]
    action_values = [S0_calc]

    # ===== Continuation to the derived portal coupling =====
    params_now = params_decoupled.copy()
    g_step = max(g_portal_max/50.0, 1e-4)
    g_val = 0.0

    while g_val < g_portal_max:
        g_val_next = min(g_val + g_step, g_portal_max)
        params_tilde['g_portal_t'] = g_val_next
        guessY = solution.y.copy()
        sol_next = solve_bvp(
            lambda r, y: bounce_odes_rescaled(r, y, params_tilde),
            lambda ya, yb: boundary_conditions_rescaled(ya, yb, params_tilde),
            solution.x, guessY,
            fun_jac=lambda r, y: bounce_odes_jac(r, y, params_tilde),
            verbose=0, max_nodes=500000, tol=1e-5,
        )
        if sol_next.success:
            solution = sol_next
            g_val = g_val_next
            params_now = params_tilde.copy()
            current_action = calculate_action(params_tilde, solution)
            R_now = wall_radius_from_solution(params_tilde, solution)
            g_values.append(g_val)
            action_values.append(current_action)
            print(f" -> Reached g_portal={g_val:.6f}, R≈{R_now:.3f}, S_E={current_action:.6f}")
        else:
            g_step *= 0.5
            if g_step < 1e-4:
                print(f" -> Continuation stalled near g ~ {g_val:.6f}.")
                break

    end_time = time.time()
    final_action = action_values[-1]
    R_final = wall_radius_from_solution(params_tilde, solution)
    print(f"\nMax stable g_portal = {g_val:.6f}")
    print(f"Final: R≈{R_final:.3f}, S_E={final_action:.6f}")
    print(f"Calculation took {end_time-start_time:.2f}s")

    # Export final background
    export_background_csv(solution, out_path=(out_prefix + "background_profile.csv"))

    return True, solution, params_now, g_values, action_values

# ---------- CLI ----------

if __name__=="__main__":
    import argparse, os
    ap = argparse.ArgumentParser()
    ap.add_argument("--v", type=float, default=4.2e-5, dest="v_phys")
    ap.add_argument("--vphi", type=float, default=9.0e-5, dest="v_phi_phys")
    ap.add_argument("--lamphi", type=float, default=0.1, dest="lam_phi_phys")
    ap.add_argument("--lam", type=float, default=0.00013, dest="lam_phys")
    ap.add_argument("--gmax", type=float, default=2.0, dest="g_portal_max")
    ap.add_argument("--R0", type=float, default=7.0)
    ap.add_argument("--w", type=float, default=3.0)
    ap.add_argument("--params", type=str, default=None, help="Path to derived params JSON")
    ap.add_argument("--summary-out", type=str, default="run_summary.json", help="Where to save summary JSON")
    ap.add_argument("--out-prefix", type=str, default="", help="Prefix for output CSV/figs")
    args = ap.parse_args()
    
    P = load_params_json(args.params) if args.params else None
    if P is not None:
        P = enforce_manifold_consistency(P)

    print("--- Starting Single 'Golden Run' Simulation ---")
    P = load_params_json(args.params) if args.params else None

    ok, sol, final_params, g_hist, action_hist = run_simulation(
    args.v_phys, args.v_phi_phys, args.lam_phi_phys, args.lam_phys,
    args.g_portal_max, args.R0, args.w, derived_params=P, out_prefix=args.out_prefix


    )

    if ok:
        # enrich params for analysis/summary
        final_params = final_params or {}
        final_params.update({'v_phys': args.v_phys, 'v_phi_phys': args.v_phi_phys, 'lam_phys': args.lam_phys})

        # figures (uses your existing function)
        analyze_solution(final_params, sol, title_prefix="Final", g_hist=g_hist, action_hist=action_hist)

        # compute S_E and R*
        S_E = calculate_action(final_params, sol)
        R_star = wall_radius_from_solution(final_params, sol)

        print("\n[done] CSVs ready:", args.out_prefix + "background_profile_g0.csv",
              "and", args.out_prefix + "background_profile.csv")
        export_run_summary(args.summary_out, final_params, g_hist, action_hist, R_star, S_E,
                           elapsed_s=None)
    else:
        print("\n[abort] No CSV written (solve unsuccessful).")
