# Code assumptions for the reviewing engineer

v0.1 implements strength design for reinforced-concrete beams (flexure and shear) and a short-column uniaxial P–M check. This file is the list to confirm. Anything marked **confirmed** is used as stated. Anything marked **assumption** is implemented so the library can run, and the calculation report repeats it, but it is not a quotation of a clause read from the printed standard.

SI units are used in the API. EIT textbook formulas that are published in ksc (kgf/cm²) and centimetres are converted with 1 kgf/cm² = 0.0980665 MPa.

## Confirmed and used as stated

### ACI 318-19

- **§5.3.1** load combinations (a) through (g). The live-load factor in (c), (d), and (e) is kept at 1.0. **§5.3.3** allows 0.5 for some occupancies; that reduction is not applied automatically.
- **§21.2.1** strength-reduction factors: 0.90 tension-controlled flexure and axial tension, 0.75 shear, 0.65 tied compression-controlled, 0.75 spiral compression-controlled.
- **§21.2.2** linear transition of φ with net tensile strain εt from εy to 0.005.
- **§22.2.2.4** Whitney block: 0.85 f'c, a = β1 c, εcu = 0.003. β1 = 0.85 for f'c ≤ 28 MPa, then reduced by 0.05 for each 7 MPa, but not below 0.65.
- **§9.6.1.2** As,min = max(0.25 √f'c / fy, 1.4 / fy) bw d.
- **§9.3.3.1** a flexural member needs εt ≥ 0.004. Design targets the tension-controlled limit εt = 0.005 so that φ = 0.90.
- **§22.5.5.1 (a)** Vc = 0.17 λ √f'c bw d for Nu = 0. Equation (b), which can give a higher Vc, is not used.
- **§22.5.1.2** Vs ≤ 0.66 √f'c bw d.
- **§22.4.2** Po = 0.85 f'c (Ag − Ast) + fy Ast. Pn,max = 0.80 Po for tied reinforcement and 0.85 Po for spiral reinforcement. φPn,max uses the compression-controlled φ.
- **§10.6.1.1** 0.01 Ag ≤ Ast ≤ 0.08 Ag.
- **§20.2.2.2** Es = 200 000 MPa.
- **§19.2.2.1** if a material omits Ec, the loader uses 4700 √f'c MPa. The example file sets Ec explicitly.

### วสท. 1008 (EIT), from published Thai sources rather than from memory of a clause number

These are used for the good-workmanship set. No EIT clause number is printed in the report, because the clause was not read from the book in this repository.

- φ = 0.90 for flexure, 0.85 for shear, 0.70 for a tied column, 0.75 for a spiral column, when construction and materials are controlled. Source: Thai Concrete Association journal comparing ว.ส.ท. 1008 with ACI 318-05, and standard วสท. 1008 column teaching notes.
- Gravity combination U = 1.4D + 1.7L. Same sources.
- Po and Pn,max = 0.80 Po (tied) or 0.85 Po (spiral). Same column notes.
- 0.01 Ag ≤ Ast ≤ 0.08 Ag. Same column notes.

The ordinary-QC set (φ flexure 0.75, tied 0.60, spiral 0.625, shear 0.70) is **not** the default. It is available as `EIT1008(quality_control="ordinary")` and is an assumption; see below.

## Assumptions — please confirm

1. **EIT wind combinations.** U = 0.75(1.4D + 1.7L + 1.7W) and U = 0.9D + 1.3W are the ACI 318-99 Section 9.2 forms from the generation วสท. 1008 follows. They are not given an EIT clause number here.
2. **EIT earthquake combinations.** W in those expressions is replaced by 1.1E, which gives 0.75(1.4D + 1.7L + 1.87E) and 0.9D + 1.43E. Confirm against the printed standard and the ministerial regulation that actually governs the project. Earthquake load E is whatever case the user defines; vertical earthquake Ev is not added automatically.
3. **EIT φ transition for columns.** For symmetric reinforcement, fy ≤ 420 MPa, and a clear distance between the extreme bar layers of at least 0.70 h, φ moves from the compression value to the flexure value as φPn decreases from 0.10 f'c Ag to zero. Otherwise the ramp starts at min(0.10 f'c Ag, φ Pb). This is ACI 318-99 §9.3.2, not a verified วสท. sentence. ACI 318-19 instead uses the εt transition.
4. **EIT beam ρmax.** Taken as 0.75 ρb (ACI 318-89), with φ staying at the flexure value when the tension steel yields. ACI 318-19 uses εt instead (0.005 for φ = 0.90, and 0.004 as the lower limit for a beam).
5. **EIT β1.** 0.85 up to 280 ksc, then minus 0.05 per 70 ksc, not below 0.65. That is the usual textbook rounding. The ACI 318-89 parent rule uses 4000 psi (about 281 ksc). At f'c = 28 MPa the two codes therefore give slightly different β1 (0.85 for ACI, about 0.846 for this EIT rounding).
6. **EIT As,min.** The ACI 318-89 expression max(3√f'c / fy, 200 / fy) bw d, evaluated in psi and converted. Textbooks often round the second term to 14/fy with fy in ksc.
7. **EIT shear coefficients,** in ksc, cm, and kgf: Vc = 0.53 √f'c b d, Vs,max = 2.12 √f'c b d, spacing break at 1.06 √f'c b d, Av,min/s = max(0.20 √f'c, 3.5) b / fy, and Av,min required when Vu > 0.5 φ Vc. These are the rounded textbook forms of the ACI 318-89 psi equations. Confirm the printed coefficients.
8. **EIT fyt** is not capped at 420 MPa. ACI caps the fyt used for Vs and Av,min at 420 MPa (§20.2.2.4).
9. **Ordinary-QC φ factors** (0.75, 0.60, 0.625, 0.70) come from a university summary of the two sets in the Thai standard, not from a clause transcribed here. Default is the good-QC set.
10. **ACI √f'c cap of 8.3 MPa** is applied to Vc and Vs,max only when minimum shear reinforcement is not provided (§22.5.3.1; confirm the printed exception).
11. **ACI Av,min trigger** is Vu > 0.083 φ λ √f'c bw d. The 0.083 is the SI form of the psi expression in §9.6.3.1. Confirm the printed metric equation. Spacing limits are min(d/2, 600 mm) or min(d/4, 300 mm) about the 0.33 √f'c threshold (§9.7.6.2.2; 0.33 is the SI form of 4√f'c in psi).
12. **Vc for beams** ignores axial force and uses only the simple 0.17 λ √f'c expression (ACI) or 0.53 √f'c (EIT). The more detailed ACI equation (b) is not evaluated.
13. **Critical shear section.** The library does not move Vu to a distance d from the support face. The model file's `face_offset_mm` is added to d and measured from the member node. For a centre-line model, set `face_offset_mm` to half the support depth. Flexural Mu is the extreme moment on the centre-line member, which for hogging is usually at the node rather than at the column face. The diagram can be read at `face_offset_mm` if the face moment is wanted; v0.1 does not reduce Mu to the face automatically.
14. **Column moments** are about the geometric centroid. v0.1 does not shift the origin to the plastic centroid of an unsymmetric bar layout.
15. **Not checked in v0.1:** slenderness and moment magnification, P-delta as a design requirement (the analysis option `PDelta` exists but is not a benchmark), creep, cracked inertia, torsion, two-way shear, development and splices, crack control, deflection, T-beams, and shear in columns. Gross elastic I is used for the frame analysis. `i_factor` on a section is only a manual stiffness modifier.

## Where the two codes deliberately differ

| Item | ACI 318-19 | EIT 1008 in this library |
|---|---|---|
| Gravity combination | 1.2D + 1.6L | 1.4D + 1.7L |
| φ shear | 0.75 | 0.85 |
| φ tied column | 0.65 | 0.70 |
| φ transition | with εt | with φPn (assumption 3) |
| Beam steel limit | εt = 0.005 | 0.75 ρb (assumption 4) |
| β1 at 28 MPa | 0.85 | about 0.846 (assumption 5) |

A singly reinforced beam whose tension steel is below both limits, with the same φ = 0.90, gets the same As from both codes. The 300×500 example beam is in that range. Shear spacing and column φPn,max still differ, although rounding of stirrup spacing can hide a small difference.
