# derive_params_from_SFV.py  (corrected)
import json, math, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--v-phys", type=float, default=4.2e-5)
ap.add_argument("--vphi-phys", type=float, required=True)
ap.add_argument("--lamphi", type=float, default=0.1, help="Option A uses 0.1")
ap.add_argument("--m-Phi", type=float, default=3.37255966313)

# Choose ONE: set phiFV & m-phi   OR   set phiFV & lambda-t
ap.add_argument("--phiFV", type=float, default=None, help="SFV-set φ_FV (>0). If omitted, √3 is used.")
grp = ap.add_mutually_exclusive_group(required=True)
grp.add_argument("--m-phi", type=float, default=None, help="Brane tail mass. Then λ_t = m_φ²/(2 φ_FV²).")
grp.add_argument("--lambda-t", type=float, default=None, help="Directly set λ_t; then m_φ = √(2 φ_FV² λ_t).")

ap.add_argument("-o","--out", default="params_derived.json")
args = ap.parse_args()

v_t   = args.vphi_phys / args.v_phys
phiFV = float(args.phiFV) if args.phiFV is not None else math.sqrt(3.0)

if args.lambda_t is not None:
    lam_t = float(args.lambda_t)
    m_phi = math.sqrt(2.0*(phiFV**2)*lam_t)
else:
    m_phi = float(args.m_phi)
    lam_t = (m_phi*m_phi)/(2.0*(phiFV**2))

mu2_t  = 0.5*lam_t*(phiFV**2 - 1.0)
bias_t = - (phiFV**4)/4.0 * lam_t / (v_t**2)
g_t    = (args.m_Phi**2 + args.lamphi*(v_t**2) - 2.0*bias_t) / (2.0*(phiFV**2))

P = {
  "notes": "Derived (no tuning) from SFV properties + O(4) tail mass or λ_t.",
  "Phi_FV": 0.0,
  "phi_FV": -phiFV,  # negative-branch FV (sign matches your BC convention)
  "v_phys": args.v_phys,
  "v_phi_phys": args.vphi_phys,
  "v_phi_t": v_t,
  "lam_phi": args.lamphi,
  "lambda_t": lam_t,
  "mu2_t": mu2_t,
  "bias_t": bias_t,
  "g_portal_t": g_t,
  "m_phi_used": m_phi,
  "m_Phi_used": args.m_Phi
}
with open(args.out,"w") as f:
    json.dump(P, f, indent=2)

print(f"[ok] wrote {args.out}")
print(f"     v_phi_t={v_t:.9f}  lam_phi={args.lamphi:g}  phi_FV={phiFV:.6f}  "
      f"lam_t={lam_t:.9e}  mu2_t={mu2_t:.9e}  bias_t={bias_t:.9e}  g_portal_t={g_t:.9f}  "
      f"m_phi={m_phi:.9f}")
