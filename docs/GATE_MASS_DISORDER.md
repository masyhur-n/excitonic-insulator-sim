# Real gate (Hartree charging), unequal masses, disorder (gap 4)

## Real gate: `python real_gate.py`
- Electrostatics of a dual-gated bilayer from Zeng & MacDonald, PRB 102,
  085154 (2020) [ZM20], Eqs. (3) and (15). In the paired phase,
  d(eV) = (C + dμ_x/dn) dn with C = 4πe²d_s/ε. Part of the gate voltage goes
  into charging the interlayer capacitor instead of tuning the gap.
- Benchmark:
  - Interlayer binding energy at d = a_B: exact Wannier solution
    **0.866 Ry***. Linear intercept of our μ_x(n): 0.855 Ry*. ZM20's
    mean-field linear fit: 0.77 Ry*. That fit includes intralayer exchange,
    which our model lacks.
  - Exciton–exciton term: ours is g ≈ +32 Ry*a_B*² (no intralayer exchange);
    ZM20 gives g = −12 (exchange-dominated). Both are used as brackets.
- Result (300 K, d_s = a_B, Ry* = 110 meV as in ZM20):
  - With T = 0 compressibility only, ZM20's attractive g gives
    49–56 mV/dec at the lowest densities. This is a negative-compressibility
    gain, like a negative-capacitance FET.
  - But at 300 K the dilute pair gas is non-degenerate
    (k_BT_deg = πRy*·na_B*² ≪ k_BT). Its thermal compressibility k_BT/n > 0
    dominates. **With it, S ≥ 285 mV/dec, about 5× worse than CMOS, and worse
    than the ideal gate (≥ 152 mV/dec).**

## Unequal masses
- Pieri, Neilson & Strinati, PRB 75, 113301 (2007), use GaAs masses
  m_e = 0.07, m_h = 0.30. Their balanced (α = 0) gaps are reproduced to ±3%
  using m_r alone (SCREENING.md). At T = 0 with equal densities, the mass
  asymmetry δξ drops out of the gap equation unless E_k ± δξ_k changes sign
  (Sarma/FFLO; their Fig. 2).
- The BKT bound uses M = m_e + m_h = 4m_r · (m_e+m_h)²/(4m_e m_h) ≥ 4m_r.
  **Equal masses, as used everywhere above, is the most optimistic case.**
  For GaAs the bound drops by a factor 1.63.

## Disorder: `python disorder.py`
- Abrikosov–Gor'kov pair breaking for bilayers: Bistritzer & MacDonald, PRL
  101, 256406 (2008) [BM08], Eq. (9).
- Benchmarks:
  - weak-disorder slope −0.7855 vs −π/4 = −0.7854 [BM08 p.2];
  - critical pair breaking 0.8819 T_c0 vs ≈ 0.88 [BM08].
- Requirement: keeping T_c0 needs τ > ħ/(1.76 k_BT_c0). That is 12 fs at
  350 K and 89 fs at 49 K. **Disorder is not the bottleneck for room-temperature
  operation.** It only removes the low-T_c cases. The mobility conversion in
  the script uses an illustrative m* = 0.44 m₀, which is not sourced.
