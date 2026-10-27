"""
Excitonic-insulator candidate materials, with ILLUSTRATIVE model parameters.

Values are order-of-magnitude literature estimates mapped onto the 2D isotropic
model in ei_model.py. Bulk / quasi-1D crystals (Ta2NiSe5, TiSe2, Ta2Pd3Te5) are
strongly anisotropic and 3D-screened, so treat them as qualitative. Calibrate
m_r, eps, r0 and Eg against ARPES / STM / optical data before quantitative use.

Every entry in `refs` was checked by web search (Oct 2026) against journal or
arXiv records. See docs/REFERENCES.md for the verification status of each one.

Eg: single-particle gap without pairing (meV); negative = band overlap (semimetal).
For gate-tunable bilayers, Eg_range gives the tuning window.
"""

MATERIALS = {
    "Ta2NiSe5": dict(
        refs=["Y. Wakisaka et al., PRL 103, 026402 (2009)",
              "Y. F. Lu et al., Nat. Commun. 8, 14408 (2017): E_op ~ 0.16 eV, T_c = 326 K, "
              "dome of T_c vs Eg centred at Eg ~ 0",
              "K. Seki et al., PRB 90, 155116 (2014): BEC regime",
              "E. Baldini et al., PNAS 120, e2221688120 (2023): symmetry breaking is structural"],
        m_r=0.25, eps=12.0, r0=2.0, d=0.0, Eg=-40.0, Eg_range=(-150, 50),
        kind="bulk, quasi-1D chalcogenide",
        status="leading candidate; T_c = 326 K, flat valence-band top in ARPES. Contested: "
               "Baldini et al. (2023) argue the symmetry breaking is structural.",
        side="near the crossover (Eg around 0, contested)"),
    "1T-TiSe2": dict(
        refs=["A. Kogar et al., Science 358, 1314 (2017): soft plasmon at ~190 K, exciton condensation"],
        m_r=0.20, eps=15.0, r0=3.0, d=0.0, Eg=-80.0, Eg_range=(-150, 20),
        kind="layered TMD semimetal",
        status="CDW at ~190-200 K with an excitonic component (plasmon softening). Strong phonon coupling.",
        side="BCS (small band overlap)"),
    "Ta2Pd3Te5": dict(
        refs=["J. Huang et al., PRX 14, 011046 (2024)",
              "P. Zhang et al., PRX 14, 011047 (2024)"],
        m_r=0.15, eps=10.0, r0=2.0, d=0.0, Eg=-20.0, Eg_range=(-100, 50),
        kind="bulk quasi-1D telluride",
        status="newer candidate (2024); gap opens with minimal lattice distortion",
        side="near the crossover"),
    "WTe2 monolayer": dict(
        refs=["Y. Jia et al., Nat. Phys. 18, 87 (2022)",
              "B. Sun et al., Nat. Phys. 18, 94 (2022): calculated E_b > 100 meV"],
        m_r=0.15, eps=4.5, r0=3.0, d=0.0, Eg=10.0, Eg_range=(-50, 60),
        kind="2D TMD, hBN-encapsulated",
        status="evidence for a 2D excitonic insulator at charge neutrality (2022); gate-tunable",
        side="BEC/crossover"),
    "MoSe2/hBN/WSe2": dict(
        refs=["Z. Wang et al., Nature 574, 76 (2019): EL threshold persisting above 100 K",
              "L. Ma et al., Nature 598, 585 (2021): exciton-compressible, charge-incompressible"],
        m_r=0.22, eps=4.5, r0=0.0, d=1.6, Eg=60.0, Eg_range=(-150, 300),
        kind="TMD double layer with hBN spacer (indirect excitons)",
        status="best platform for a TUNABLE BEC-BCS crossover: Eg is set by an electric "
               "field, d by the number of hBN layers. Evidence, not yet proof, of condensation.",
        side="tunable from BEC to BCS"),
    "InAs/GaSb QW": dict(
        refs=["L. Du et al., Nat. Commun. 8, 1971 (2017): gap ~ 2 meV, T_c ~ 10 K"],
        m_r=0.022, eps=15.0, r0=0.0, d=10.0, Eg=-30.0, Eg_range=(-60, 20),
        kind="III-V broken-gap quantum well",
        status="weak-coupling (BCS) EI: gap ~2 meV, T_c ~10 K, with quantized edge conductance",
        side="BCS"),
    "double bilayer graphene": dict(
        refs=["J. I. A. Li et al., Nat. Phys. 13, 751 (2017)",
              "X. Liu et al., Nat. Phys. 13, 746 (2017)",
              "A. Perali, D. Neilson, A. R. Hamilton, PRL 110, 146803 (2013) (theory, screening)"],
        m_r=0.02, eps=4.5, r0=0.0, d=2.0, Eg=-20.0, Eg_range=(-40, 20),
        kind="graphene bilayers + hBN, magnetic-field-assisted",
        status="interlayer exciton condensate seen in the quantum Hall regime (zero-field case needs more)",
        side="BCS to BEC with filling (in a B field)"),
}
