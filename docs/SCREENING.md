# Self-consistent screening (gap 1)

`screening.py` (solver + benchmark), `tc_ceiling.py` (design rule)

## Method
- Fixed-density e–h bilayer, T = 0 mean field. Conventions from Pieri,
  Neilson & Strinati, PRB 75, 113301 (2007) [P07]: one species per layer,
  n = Σ v_k², r_s = 1/√(πn), units a_B* and Ry*.
- Screening in the superfluid state (SS), Eq. (4) of Neilson, Perali & Hamilton,
  PRB 89, 060502(R) (2014), arXiv:1308.0280 [N14]. The calculation starts from
  the unscreened coherent state and iterates the BCS normal plus anomalous
  polarisability to self-consistency. Π(q → 0) = 0 in the superfluid, and Π
  reduces to the 2D Lindhard function when Δ → 0.

## Benchmarks
| quantity | paper | this code |
|---|---|---|
| unscreened Δ_max, d = 1, r_s = 1.5 / 3 / 20 | ≈ 0.66 / 0.68 / 0.18 Ry* [P07 Fig. 1b] | 0.659 / 0.700 / 0.173 |
| μ = 0 crossing, d = 1 | r_s ≈ 3 [P07] | μ = +0.07 Ry* at r_s = 3, −0.04 at 4 |
| unscreened Δ_max, d = 1, r_s = 0.5 | ≈ 0.22 Ry*, peak at k ≈ 4 [P07 Fig. 1] | 0.24–0.25 at k = 4.0 (needs nk ≥ 280; 0.13 at the default nk = 140) |
| SS onset (condensate jump), d = 1 | SS ≈ 3; DQMC ≈ 2.5 [N14 Fig. 1–2] | **2.38**, jump 0.015 → 0.38 |
| SS onset, d = 0.7 | DQMC ≈ 2.1 [N14 Fig. 1a] | **1.97** |
| SS → US at low density | equal for r_s ≳ 5 [N14 Fig. 3] | SS/US = 0.91 (r_s = 5), 0.94 (r_s = 6) |

Paper values are read from figures (about ±0.1). Our onset sits slightly
left of N14's SS line and close to DQMC.

The condensate fraction uses c = Σu²v²/Σv². N14 does not print its
definition, so only the jump location is compared, not c itself.

## Result: screening ceiling on T_c
The SS superfluid exists only for r_s > r_s,c(d/a_B), with
r_s,c = 1.29, 1.97, 2.38, 2.98, 3.53, 4.46 at d/a_B = 0.3, 0.7, 1, 1.5, 2, 3.
With n_s ≤ n, the Nelson–Kosterlitz bound gives **T_c ≤ Ry*/(4 r_s,c²)**.

- MoSe₂/hBN/WSe₂ (m_r = 0.22, ε = 4.5, d = 1.6 nm): **T_c ≤ 49 K**. With the
  thinnest barrier (d = 0.7 nm): ≤ 121 K.
- Over the box m_r 0.05–0.6, ε 2–15, d 0.6–3 nm: **no point reaches 350 K**.
  36 of 90 000 reach 300 K, all at ε = 2 (near vacuum). **None with ε ≥ 4.5
  (hBN) reaches 300 K.**
- Valley/spin degeneracy (TMDs, graphene) adds screening channels and lowers
  this ceiling further. A 2023 DQMC study (PRB 107, L041409) finds that
  valley degeneracy shrinks the condensate region.

Robustness: for d/a_B < 0.3, `tc_ceiling.py` holds r_s,c at its d/a_B = 0.3
value. Extrapolating r_s,c linearly below 0.3 instead (floored at 1.0 or 0.5)
leaves the result unchanged: the maximum is 328 K, nothing reaches 350 K, and
no point with ε ≥ 4.5 reaches 300 K.

## Limitations
- T = 0 SS screening, static RPA, no vertex corrections. N14 argues these are
  small for r_s ≲ 3 and r_s ≳ 5.
- Equal masses in Π.
- The 2D bound uses n_s = n, so the ceiling is generous.
