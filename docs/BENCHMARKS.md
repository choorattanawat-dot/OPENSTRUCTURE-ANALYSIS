# Benchmarks

Each test below names the source of the expected value. The arithmetic is in the test, not only in this note.

## Analysis

Classical Euler-Bernoulli beam theory, as tabulated in structural-analysis textbooks (Hibbeler, *Structural Analysis*, beam-deflection tables, and the same closed forms in any strength-of-materials reference). Shear deformation is not included, and OpenSees `elasticBeamColumn` does not include it either.

| Case | Result | Test |
|---|---|---|
| Simply supported beam, uniform load `w` | `δ = 5 w L^4 / 384 EI`, `M = w L^2 / 8`, `R = w L / 2` | `tests/test_analysis.py` |
| Fixed beam, uniform load | `M_end = w L^2 / 12`, `M_mid = w L^2 / 24`, `R = w L / 2` | `tests/test_analysis.py` |
| Simply supported beam, midspan point load `P` | `δ = P L^3 / 48 EI`, `M = P L / 4` | `tests/test_analysis.py` |
| Fixed portal, flexible girder | Independent direct-stiffness solution in the test file (axial + bending). Same element assumptions as OpenSees `elasticBeamColumn` with a linear transformation. | `tests/test_analysis.py` |
| Fixed portal, very stiff girder, horizontal load `P` | `Δ = P h^3 / 24 EI_column`, `M_base = P h / 4`. This closed form ignores axial deformation. The test allows 5% for that difference. | `tests/test_analysis.py` |

Section forces for the fixed beam and the simple-beam point load are also checked from equilibrium alone, without OpenSees, in `tests/test_section_forces.py`. Global equilibrium of every example combination (reactions + applied loads) is checked in `tests/test_example_and_cli.py`.

## Design

Hand arithmetic is in `tests/test_design_hand_calcs.py`.

| Check | What is verified |
|---|---|
| ACI singly reinforced beam, Mu = 200 kN·m, 300×500, d = 450, f'c = 28 MPa, fy = 400 MPa | `ρ` from `Rn = Mu / (φ b d²)` with φ = 0.90, and `Mu / φMn ≈ 1` from an independent Whitney-block solution |
| ACI doubly reinforced beam, Mu = 400 kN·m, d' = 60 mm | Compression steel is required at the εt = 0.005 limit and the provided steel satisfies the section check |
| ACI shear, Vu = 200 kN, 2-leg 9 mm, fyt = 400 MPa | `Vc = 0.17 √f'c bw d`, φ = 0.75, spacing rounded down from `Av fyt d / Vs` |
| EIT shear | `Vc = 0.53 √f'c b d` in ksc, cm, kgf, and φ = 0.85 |
| Column Po | `Po = 0.85 f'c (Ag − Ast) + fy Ast` for four 500 mm² bars in a 400 mm square |
| φPn,max | ACI tied `0.65 × 0.80 Po`; EIT tied `0.70 × 0.80 Po` |
| Balanced strain point | Explicit Whitney arithmetic at `c = 0.003/(0.003+0.002) × 340 mm` matches `nominal_at_c` |
| Utilization | `Pu = 0.5 φPn,max` gives about 0.5; `Pu = 1.2 φPn,max` is not acceptable |

Code-factor sources and the items still marked as assumptions are in `docs/CODE_ASSUMPTIONS.md`.
