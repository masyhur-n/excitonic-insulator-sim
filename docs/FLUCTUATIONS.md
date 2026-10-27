# Fluctuation test: does the steep mean-field switching survive?

`python fluctuations.py [T_op] [n_Eg]` → `fluctuations_<T>K.png`

## Question

`device_scoring.py` found sub-CMOS switching only where mean field turns the
gap on abruptly as Eg crosses the EI transition. Mean-field theory has no
fluctuations, and near a 2D transition fluctuations are large: R&T 2014 (p.9)
notes the critical region is of order unity near the crossover. This test adds
Gaussian fluctuations of the order parameter.

## Method

1. **Landau function** Ω(A) along Δ_k = A·φ_k, from the Gibbs–Bogoliubov
   variational free energy. This is a true upper bound for any trial gap, and
   it is stationary exactly at the gap equation (`ei_model.py`).
2. **Ginzburg–Landau coefficients** a, b from a quartic fit of Ω(A). The pair
   field is normalised to the pair density, |ψ|² = Σ_k F_k². The gradient term
   is c = ħ²/2M with M = m_e + m_h. This is exact in the BEC limit. On the BCS
   side it overestimates the stiffness, so it underestimates fluctuations, which
   makes the test *conservative* (it favours steep switching).
3. **Self-consistent Gaussian (Hartree / large-N)** treatment of the quartic
   term [Ullah & Dorsey 1991; Larkin & Varlamov 2005]. It uses Bose occupation
   of the pair modes:
   `m = Σ_Q n_B(a + b·m + cQ²)`.
   Gaussian fluctuations about the normal state are equivalent to NSR
   [R&T 2014, p.9; Nozières & Schmitt-Rink 1985]. The Bose occupation matters:
   the classical form T/(a + cQ²), tested first, produced ~10¹² cm⁻²
   thermal excitons deep in the insulator. A Boltzmann estimate gives about
   10⁻⁵ of that.
4. **Gap with fluctuations:** E_gap = 2 min √(ξ² + ⟨A²⟩φ²). This is the
   standard pseudogap form Δ² = Δ_sc² + Δ_pg² [Chen, Stajic, Tan & Levin 2005,
   ref. [12] in R&T]. It overestimates the gap deep in the insulating region,
   which doesn't affect the minimum swing.

Checks:
- The GL quartic reproduces the mean-field amplitude: A²(GL)/A²(MF) = 0.96–1.04.
- No discontinuities in the cut-off or the fit range across the transition.
- The Eg grid is refined to 121 points.

## Result

| T_op | material | mean-field S_min | with fluctuations (Λ = 0.5 / 1 / 2 × 1/a_B) | CMOS |
|---|---|---|---|---|
| 350 K | MoSe₂/hBN/WSe₂ | **60.4** (sub-CMOS) | 140 / 133 / 102 | 69.4 |
| 350 K | WTe₂ monolayer | 117.5 | 463 / 219 / 147 | 69.4 |
| 300 K | MoSe₂/hBN/WSe₂ | 88.1 | 120 / 120 / 104 | 59.5 |
| 300 K | WTe₂ monolayer | **38.6** (sub-CMOS) | 280 / 169 / 137 | 59.5 |

Robustness: changing the GL fit range (1–3 k_BT) gives S_min = 120–204 mV/dec.
All values stay above CMOS.

**The steep switching does not survive fluctuations.** Both mean-field
sub-CMOS points (60.4 mV/dec at 350 K, 38.6 mV/dec at 300 K) come out at
100–470 mV/dec, i.e. 1.5–8× worse than CMOS. The abrupt onset turns into a
smooth pseudogap crossover. This removes the last route to sub-CMOS switching found
in the mean-field screen.

## Limitations

- Gaussian / Hartree level. It is not a full BKT (vortex) treatment; vortices
  would add more smoothing, not less.
- Fixed pair shape φ_k. Static fluctuations only (classical in time), with no
  dynamics or transport.
- Same model assumptions as before: parabolic, equal-mass bands, estimated
  material parameters and an ideal gate.

## References
Ullah & Dorsey, PRB 44, 262 (1991); Larkin & Varlamov, *Theory of Fluctuations
in Superconductors* (OUP 2005); Randeria & Taylor (2014); Nozières &
Schmitt-Rink, JLTP 59, 195 (1985); Chen, Stajic, Tan & Levin, Phys. Rep. 412, 1
(2005); Gibbs–Bogoliubov–Feynman inequality (see REFERENCES.md).
