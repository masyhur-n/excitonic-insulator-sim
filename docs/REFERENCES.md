# References and verification status

Every reference used in the code, documentation or manuscript, with how it was checked.

| tag | meaning |
|---|---|
| **FULL** | full text read (publisher or arXiv PDF) |
| **ABS** | bibliographic record and main claim checked against the journal, arXiv, PMC or an institutional repository; full text not read |
| **CIT** | cited via the reference list of a FULL source; original not opened |
| **REPRO** | the quoted number is reproduced by this code (`run_benchmarks.py`) |
| **⚠** | not verified; check before relying on it |

## Benchmarks reproduced by the code (manuscript Table I)

| reference | quantity | literature | this code | status |
|---|---|---|---|---|
| M. Randeria & E. Taylor, "Crossover from Bardeen-Cooper-Schrieffer to Bose-Einstein Condensation and the Unitary Fermi Gas", Annu. Rev. Condens. Matter Phys. 5, 209 (2014), arXiv:1306.5785 | 3D unitarity mean field Δ/E_F, μ/E_F (Table I) | 0.68, 0.59 | 0.686, 0.591 | FULL + REPRO |
| same | BCS-limit T_c = (8γ/πe²)E_F e^{−π/2k_F\|a\|} (p.9) | formula | < 1% | FULL + REPRO |
| D. Jérome, T. M. Rice, W. Kohn, Phys. Rev. 158, 462 (1967) | EI exists only for E_g < E_b | — | yes | ABS + REPRO |
| textbook | 2D exciton E_b = 4 Ry* | 435.4 meV | 432.9 meV | REPRO |
| P. Pieri, D. Neilson, G. C. Strinati, PRB 75, 113301 (2007), arXiv:cond-mat/0610311 | Δ_max at d = a_B, r_s = 1.5, 3, 20 (Fig. 1b, read from figure) | ≈ 0.66, 0.68, 0.18 Ry* | 0.659, 0.700, 0.173 | FULL + REPRO |
| D. Neilson, A. Perali, A. R. Hamilton, PRB 89, 060502(R) (2014), arXiv:1308.0280 | superfluid onset r_s,c at d = a_B (Figs. 1–2) | DQMC ≈ 2.5; SS ≈ 3 | 2.38 | FULL + REPRO |
| same | onset at d = 0.7 a_B | DQMC ≈ 2.1 | 1.97 | FULL + REPRO |
| same | SS gap → US gap for r_s ≳ 5 (Fig. 3) | equal | 0.94 at r_s = 6 | FULL + REPRO |
| M. Yu. Kharitonov & K. B. Efetov, PRB 78, 241401(R) (2008), arXiv:0808.2164 | V_us(p_F)/V(p_F), r_s = 1 (p.3) | ≈ 9 | 8.99 | FULL + REPRO |
| same | maximum coupling λ_c (Eq. 2) | 1/16 | 0.0624 | FULL + REPRO |
| same | kernel log-slope = λ_c (Eqs. 10–12) | 0.0567 | 0.0546 | FULL + REPRO |
| T. Kaneko, T. Toriyama, T. Konishi, Y. Ohta, PRB 87, 035121 (2013), arXiv:1210.2787 | distortion onset λ_c (Fig. 3f) | 0.125 eV | 0.125 eV | FULL + REPRO |
| same | EI onset V_c (Fig. 3e) | 0.6 eV | 0.69 eV (+15%) | FULL; **deviation** (supplement not accessed) |
| same | k_BT_c at V = 0.55, λ = 0.04 eV (Fig. 4) | 0.0783 eV | 0.0564 eV (−28%) | FULL; **deviation** |
| Y. Zeng & A. H. MacDonald, PRB 102, 085154 (2020), arXiv:2006.09937 | interlayer binding at d = a_B (Fig. 2a) | 0.77 Ry* (linear fit) | 0.866 Ry* (exact) | FULL + REPRO |
| R. Bistritzer & A. H. MacDonald, PRL 101, 256406 (2008), arXiv:0808.1310 | Abrikosov–Gor'kov slope; critical pair breaking (Eq. 9, p.2) | −π/4; 0.88 T_c0 | −0.7855; 0.8819 T_c0 | FULL + REPRO |
| S. K. Banerjee, L. F. Register, E. Tutuc, D. Reddy, A. H. MacDonald, IEEE Electron Device Lett. 30, 158 (2009), doi:10.1109/LED.2008.2009362 | switching energy (p.160) | 0.008 aJ (SPICE) | 0.0044 aJ (C·V²) | FULL + REPRO |

## Theory used (not benchmarked numerically)

| reference | used for | status |
|---|---|---|
| H. Min, R. Bistritzer, J.-J. Su, A. H. MacDonald, PRB 78, 121401(R) (2008), arXiv:0802.3462 | Dirac two-band model, chirality factor (Eq. 3), T_KT ≤ E_F/8 (p.3) | FULL |
| G. Mazza et al., PRL 124, 197601 (2020), arXiv:1911.11835 | discrete symmetry breaking in Ta₂NiSe₅; "no dissipationless transport or excitonic superfluidity" (p.5) | FULL |
| D. R. Nelson & J. M. Kosterlitz, PRL 39, 1201 (1977) | universal jump, T_BKT bound | ABS |
| L. V. Keldysh & Yu. V. Kopaev, Sov. Phys. Solid State 6, 2219 (1965) | BCS-like excitonic insulator | ABS (citation only; no full text found) |
| Yu. E. Lozovik & V. I. Yudson, Sov. Phys. JETP 44, 389 (1976) | spatially separated pairing | ABS |
| T. C. Berkelbach, M. S. Hybertsen, D. R. Reichman, PRB 88, 045318 (2013) | Rytova–Keldysh 2D potential | ABS |
| P. Nozières & S. Schmitt-Rink, J. Low Temp. Phys. 59, 195 (1985) | NSR = Gaussian fluctuations | CIT (Randeria & Taylor ref. 16) |
| L. P. Gor'kov & T. K. Melik-Barkhudarov, Zh. Eksp. Teor. Fiz. 40, 1452 (1961) | (4e)^{1/3} ≈ 2.2 mean-field overestimate | CIT (Randeria & Taylor ref. 38, p.10) |
| Q. Chen, J. Stajic, S. Tan, K. Levin, Phys. Rep. 412, 1 (2005) | pseudogap Δ² = Δ_sc² + Δ_pg² | CIT (Randeria & Taylor ref. 12) |
| S. Ullah & A. T. Dorsey, "Effect of fluctuations on the transport properties of type-II superconductors in a magnetic field", PRB 44, 262 (1991), doi:10.1103/PhysRevB.44.262 | self-consistent Hartree treatment of the GL quartic term | ABS; title, volume and pages confirmed via Crossref |
| A. Larkin & A. Varlamov, *Theory of Fluctuations in Superconductors* (OUP, 2005) | GL fluctuation theory | ABS (publisher record) |
| Gibbs–Bogoliubov–Feynman inequality | variational Landau function | textbook |
| A. Perali, D. Neilson, A. R. Hamilton, PRL 110, 146803 (2013) | screening, double-bilayer graphene | ABS |
| S. De Palo, F. Rapisarda, G. Senatore, PRL 88, 206401 (2002), arXiv:cond-mat/0201414 | DQMC phase boundary (values taken as reproduced in Neilson et al. 2014) | FULL (downloaded) |

## Devices and motivation

| reference | used for | status |
|---|---|---|
| International Energy Agency, *Energy and AI* (2025) | data-centre electricity ~415 TWh (2024) → ~945 TWh (2030) | ABS (IEA summary page + secondary coverage) |
| A. M. Ionescu & H. Riel, Nature 479, 329 (2011), doi:10.1038/nature10679 | 60 mV/decade thermionic limit | ABS (one page of the full text seen) |
| D. Unuchek et al., Nature 560, 340 (2018) | room-temperature excitonic transistor without condensate | ABS |
| S. M. Sze, *Physics of Semiconductor Devices* | intrinsic activation n_i ∝ e^{−E_g/2k_BT} | textbook ⚠ |

## Materials (experiment)

| material | reference | fact used | status |
|---|---|---|---|
| Ta₂NiSe₅ | Y. Wakisaka et al., PRL 103, 026402 (2009) | ARPES evidence for an EI | ABS |
| | Y. F. Lu et al., Nat. Commun. 8, 14408 (2017) | T_c = 326 K; optical gap ≈ 0.16 eV; T_c dome vs E_g | ABS (open access) |
| | K. Seki et al., PRB 90, 155116 (2014) | BEC regime | ABS |
| | E. Baldini et al., PNAS 120, e2221688120 (2023) | symmetry breaking is structural | ABS |
| 1T-TiSe₂ | A. Kogar et al., Science 358, 1314 (2017) | exciton condensation signatures (~190 K) | ABS |
| Ta₂Pd₃Te₅ | J. Huang et al., PRX 14, 011046 (2024); P. Zhang et al., PRX 14, 011047 (2024) | EI candidate | ABS ⚠ titles/authors not checked on APS |
| WTe₂ monolayer | Y. Jia et al., Nat. Phys. 18, 87 (2022); B. Sun et al., Nat. Phys. 18, 94 (2022) | 2D EI | ABS |
| MoSe₂/WSe₂ | Z. Wang et al., Nature 574, 76 (2019); L. Ma et al., Nature 598, 585 (2021) | interlayer exciton condensation / EI | ABS |
| InAs/GaSb | L. Du et al., Nat. Commun. 8, 1971 (2017) | gap ≈ 2 meV, T_c ≈ 10 K | ABS |
| double bilayer graphene | J. I. A. Li et al., Nat. Phys. 13, 751 (2017); X. Liu et al., Nat. Phys. 13, 746 (2017) | condensate in the quantum Hall regime | ABS |

## Software
F. Pedregosa et al., *Scikit-learn*, J. Mach. Learn. Res. 12, 2825 (2011).

## Not sourced (stated as such in the manuscript)
- The illustrative model parameters in `materials.py` (m_r, ε, r0, d, E_g) for
  the first-generation screen are order-of-magnitude estimates, not literature values. The
  manuscript's quantitative results use instead the published Kaneko et al.
  (2013) and Zeng & MacDonald (2020) parameters, or scans over wide parameter
  boxes.
- The "buildable double layer" limits (d ≥ 0.6 nm, |E_g| ≤ 300 meV) are estimates.
- ⚠ Not verified at all: Comte & Nozières, J. Physique 43, 1069 (1982), Part I;
  Keldysh, JETP Lett. 29, 658 (1979) (journal name uncertain). Neither is cited
  in the manuscript.

## Bibliographic check
All journal entries in `paper/refs.bib` were checked against Crossref (title,
volume, first page) on 2026-10-27, and DOIs were added for 30 of them. The entries
without a Crossref DOI are the Soviet journals (Keldysh–Kopaev 1965,
Lozovik–Yudson 1976, Gor'kov–Melik-Barkhudarov 1961), J. Low Temp. Phys. 1985,
JMLR 2011, the IEA report and the Larkin–Varlamov book.
