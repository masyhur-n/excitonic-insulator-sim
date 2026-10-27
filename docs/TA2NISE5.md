# Ta₂NiSe₅: quasi-1D chains and lattice coupling (gaps 3 and 6)

`python ta2nise5_chain.py`

## Model
- Published three-chain model of Kaneko, Toriyama, Konishi & Ohta, PRB 87,
  035121 (2013) [K13], Eq. (1), with its DFT-fitted parameters:
  t_c = −0.8, t_f = 0.4, D = 0.2 eV; U = 4V; exciton–phonon coupling λ = g²/ω
  (uniform shear).
- K13's mean-field equations are given in its Supplemental Material, which was not accessed here,
  so the equations here are **derived from their Eq. (1)**: Hartree for U and V,
  Fock on both c–f bonds, a static q = 0 phonon, and opposite order on the two
  Ta chains, so the nonbonding c-band decouples (as stated in K13).

## Benchmarks
| quantity | K13 | this code |
|---|---|---|
| distortion onset at V = 0 | λ_c = 0.125 eV (Fig. 3f) | **0.125 eV** |
| EI onset at λ = 0 | V_c = 0.6 eV (Fig. 3e) | 0.692 eV (+15%) |
| k_BT_c at V = 0.55, λ = 0.04 | 0.0783 eV (Fig. 4) | 0.0564 eV (−28%) |

Two alternative Fock decouplings bracket K13's numbers:
- on-site only: V_c ≈ 1.0 eV, no T_c;
- both bonds with the on-site amplitude: V_c ≲ 0.5 eV, k_BT_c = 0.126 eV.

The remaining 15–28% mismatch is attributed to the unavailable supplement.
Closing it requires the K13 Supplemental Material.

## Result: excitonic vs structural (V = 0.55 eV)
| λ (eV) | excitonic hybridisation | lattice δ | excitonic share | gap | T_c (MF) → ÷2.2 |
|---|---|---|---|---|---|
| 0.02 | 0 | 0 | – | 0.200 | no order |
| 0.04 | 0.048 | 0.017 | **74%** | 0.220 | 655 K → 298 K |
| 0.05 | 0.061 | 0.028 | **69%** | 0.236 | 830 K → 377 K |

- K13 estimate λ ≈ 0.02–0.05 eV from the measured 0.5–1° distortion. In that
  range, order appears only cooperatively, and it is about 70% electronic.
- Purely structural order needs λ ≥ 0.125 eV, three times that estimate.
- The ÷2.2 mean-field correction (Randeria & Taylor p.10) gives 298–377 K,
  which brackets the measured T_c = 326 K [Lu et al. 2017]. This is a
  consistency check, not a fit.
- The model gap (0.22–0.24 eV) is above the measured optical gap of 0.16 eV
  [Lu 2017].

## Relation to the symmetry debate
Mazza et al., PRL 124, 197601 (2020) [arXiv:1911.11835] use realistic Wannier
models with U ≈ 2.1–2.4 eV and V ≈ 0.9 eV. They find:
- a purely electronic excitonic order parameter that breaks only **discrete**
  mirror symmetries, because Ta–Ni hopping (≈ 36 meV) is symmetry-allowed;
- an excitonic phase only in a narrow V window near a Lifshitz (zero-gap) point;
- and they conclude: "there is no dissipationless transport or excitonic
  superfluidity" (p.5).

Baldini et al., PNAS 2023 argue that the symmetry breaking is structural.

Consequences for this project:
1. **Ta₂NiSe₅ cannot host a BiSFET-type condensate device.** There is no
   continuous U(1) symmetry, so no counterflow superfluid.
2. Only gap switching is possible, and Ta₂NiSe₅ is a bulk crystal without
   electrostatic gap control (pressure or strain only).
3. Within K13's model the order is mostly electronic but needs the lattice.
   This is compatible with both camps: excitonic correlations dominate the
   gap, while the lattice decides whether order sets in.
