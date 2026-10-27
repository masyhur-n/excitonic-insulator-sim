# Dirac bands: graphene double layer, the BiSFET platform (gap 2)

`python dirac.py`

## Model
- Two-band Dirac model of Min, Bistritzer, Su & MacDonald, PRB 78, 121401(R)
  (2008) [M08]: ξ = ħvk − E_F, with the chirality form factor (1+cos θ)/2 in
  the interlayer interaction [M08 Eq. 3], and N_sN_v = 4 spin–valley copies.
  The bare-Coulomb angular integral with the chirality factor is analytic
  (elliptic K and E). The log singularity is cell-averaged.
- Screening:
  - **US** (unscreened, as in M08).
  - **NS** (normal-state RPA, Kharitonov & Efetov, PRB 78, 241401(R) (2008)
    [KE08] Eq. 5).
  - **SS** (superfluid-state, the approach validated against DQMC by Neilson,
    Perali & Hamilton 2014): intraband Π from BCS coherence factors with
    chirality; graphene's interband polarisability kept at its normal-state
    value, N_sN_vν − Π_intra,normal.

## Benchmarks
| quantity | paper | this code |
|---|---|---|
| V_us(p_F)/V_NS(p_F), r_s = e²/εħv = 1, d = 0 | ≈ 9 [KE08 p.3] | 8.99 |
| maximum coupling λ_c (α → ∞, d → 0) | 1/16 = 0.0625 [KE08 Eq. 2] | 0.0624 |
| λ_c at α = 1 from Eq. 12 vs slope of our numerical kernel eigenvalue in ln(1/T) | equal [KE08 Eqs. 10–12] | 0.0567 vs 0.0546 (4%), grid-converged |
| T_c (NS), α = 1 | ~10⁻⁷ E_F [KE08 Eq. 3] | ≈ 2×10⁻⁷ E_F (extrapolated) |
| BiSFET density → E_F | 4.9×10¹² cm⁻² → 0.1 E_F/k_B ≈ 300 K [Banerjee 2009] | E_F = 258 meV → 300 K |
| US gap scale | ≈ 0.1 E_F (M08 Fig. 3) | Δ = 31 meV = 0.12 E_F (ε = 3.9, d = 1 nm) |

## Result at the BiSFET operating point (E_F = 258 meV, d = 1 nm, ε = 3.9)
| screening | Δ_max | consequence |
|---|---|---|
| US | 31 meV | T_KT ≤ E_F/8 = 374 K [M08]: the room-temperature premise |
| NS | 0 (T_c ~ 10⁻⁷ E_F) | the KE08 conclusion |
| SS | **0**, collapses within 2 iterations | agrees with N14 p.4: "no evidence of superfluidity … consistent with … the (SS) mean-field approach" |

Even in the most favourable corner (ε = 1, d = 0.3 nm, α = 2.2), the SS gap is
only 0.087 meV, a mean-field scale of about 1 K. **Under both literature-validated
screening treatments, monolayer-graphene BiSFETs do not condense at any
useful temperature.** This settles the Min vs Kharitonov–Efetov disagreement
from REFERENCES.md, with the DQMC-validated SS approach on the KE side for
monolayer graphene.

## Limitations
- Two-band model without intralayer exchange (M08 includes it in Δ^z).
- Π_inter is kept at its normal-state value. It is constant for q ≤ 2k_F, and
  for q > 2k_F the real Π grows, so the screening estimate is if anything too weak.
- T = 0 only. Λ = 8 k_F UV cut-off.
