
# Phase-1: Derived Two-Field Bounce in the SFV dSB Model

This repo reproduces the *derivation-first* bounce used in the preprint. Every parameter is derived from SFV properties and FV observables; **no scans** are used. A single, transparent calibration of the FV curvature in the SFV sector (`m_Phi_used`) anchors the narrow `g(λ_t)` corridor.

## 0) Environment

- Python 3.10+
- `numpy`, `scipy`, `pandas`, `matplotlib` (for optional plots)

## 1) Files

- `goldenRunDetails_v4f_more3_derived_v2.py` — solver/BVP with manifold enforcement.
- `derive_params_from_SFV.py` — builds a params JSON from SFV properties.
- `params_corridor.json` — example JSON used to achieve `S_E ≈ 1620` (see below).
- Outputs: `background_profile_g0.csv`, `background_profile.csv`, `run_summary.json`.

## 2) Quickstart (Option A: λ_φ = 0.1)

**A. Build a corridor-consistent JSON** (set φ_FV from μ²/λ, calibrate `m_Phi_used`, fix λ_t):

```bash
python derive_params_from_SFV.py ^
  --v-phys 4.2e-5 --vphi-phys 8.9e-5 --lamphi 0.1 ^
  --phiFV 1.79984323175 ^
  --lambda-t 2.0824e-4 ^
  --m-Phi 3.5089 ^
  -o params_corridor.json
```

> Tip: You can specify `--m-phi` instead of `--lambda-t` if you prefer to anchor via tail mass.

**B. Run the bounce**

```bash
python goldenRunDetails_v4f_more3_derived_v2.py --params params_corridor.json
```

The solver enforces the FV manifold at startup:
- recomputes `phi_FV` from `mu2_t/lambda_t`,
- sets `bias_t = -(phi_FV^4/4) * lambda_t / v_t^2`,
- uses `m_Phi_used` to compute the required `g` from the FV relation.

**Expected outcome:** continuation in `g` converges cleanly on the corridor; recent runs give `S_E ≈ 1620`.

## 3) Measuring from a g=0 profile (future work)

We include an optional mode to read `background_profile_g0.csv` and fit `(phi_FV, m_phi, m_Phi)` from the O(4) tail, then derive all parameters. In practice, tail-window robustness is delicate; if the auto-fit fails quality checks, prefer the corridor-calibrated JSON above.

```bash
python derive_params_from_SFV.py --profile background_profile_g0.csv ^
  --v-phys 4.2e-5 --vphi-phys 8.9e-5 --lamphi 0.1 ^
  -o params_from_g0.json
```

or corridor-locked:

```bash
python derive_params_from_SFV.py --profile background_profile_g0.csv ^
  --v-phys 4.2e-5 --vphi-phys 8.9e-5 --lamphi 0.1 ^
  --target-g 1.968 ^
  -o params_from_g0_corridor.json
```

## 4) Repro checklist

- [ ] Commit the exact JSON used (`params_corridor.json`).
- [ ] Commit `run_summary.json` with `S_E` and `R_*` from the solver.
- [ ] Record `m_Phi_used` in the JSON; keep it fixed across λ_t nudges.
- [ ] Keep `v_phi_t = v_phi_phys / v_phys` in the JSON (don’t edit `v_phi_phys` alone).

## 5) Notes on the solvability corridor

At fixed `(m_Phi, λ_φ, v_t, φ_FV)` the FV consistency enforces
\[
g(λ_t) = \frac{m_\Phi^2 + λ_φ v_t^2}{2 φ_FV^2} + \frac{φ_FV^2}{4 v_t^2} λ_t.
\]

Small deviations from this line lead to trivial solutions or exploding action. The solver corrects minor JSON inconsistencies before solving, ensuring it stays on the corridor.

## 6) Citing the preprint

Please cite the preprint included in this repo (see `paper.tex`).
