# Mean-field device feasibility screen

`python device_scoring.py [T_op]` (default 350 K) → `results/device_scoring.png`

This is the mean-field stage. Fluctuations, gate electrostatics and screening are added in FLUCTUATIONS.md, GATE_MASS_DISORDER.md and SCREENING.md. The manuscript quotes the mean-field swing from `fluctuations.py` (121-point E_g grid: 60.4 mV/dec for MoSe₂/hBN/WSe₂ at 350 K). The 81-point grid used here gives 66.0 mV/dec.

**This is a feasibility screen, not a device simulation.** It asks whether the
*physics* (pairing and condensation at the operating temperature, and how
sharply the gap responds to a gate) could beat CMOS. Contacts, transport,
disorder, parasitics and fabrication are not modelled.

## Why

- AI is the main driver of data-centre electricity use, which the IEA projects
  will grow from about 415 TWh (2024) to about 945 TWh (2030) (IEA, *Energy and
  AI*, 2025).
- A CMOS transistor needs at least ln(10)·k_BT/e of gate voltage per decade of
  current: 60 mV/dec at 300 K. That floor limits how far the supply voltage,
  and so the switching energy, can be lowered (Ionescu & Riel, Nature 479, 329, 2011).
- Exciton condensates have been proposed as a way around that floor. The
  BiSFET (Banerjee et al., IEEE EDL 30, 158, 2009) would work at or below
  k_BT/q ≈ 26 mV, but it needs a room-temperature excitonic condensate.

## Choices and assumptions

| item | value | why / source |
|---|---|---|
| operating temperature T_op | **350 K** | chip hot-spot temperature; `python device_scoring.py 300` for room temperature |
| gate model | **ideal ΔEg = e·ΔV** | best case for interlayer bias in a double layer; a real lever arm α < 1 makes S worse by 1/α (see GATE_MASS_DISORDER.md) |
| conduction | G ∝ exp(−E_gap / 2k_BT) | intrinsic activation, n_i ∝ e^(−Eg/2kT) (textbook, e.g. Sze, *Physics of Semiconductor Devices*; not re-checked) |
| swing | S = 2 ln(10) k_BT / \|dE_gap/dEg\| | follows from the two lines above. An ordinary gap-tuned device (slope 1) gives 2× the CMOS limit |
| pairing needed for gap switching | T* > T_op | a pairing pseudogap is enough, coherence is not required (Randeria & Taylor 2014, p.10) |
| mean-field correction | T*/2.2 | MFT overestimates T_c by (4e)^{1/3} ≈ 2.2 even in the BCS limit (R&T p.10; Gor'kov & Melik-Barkhudarov 1961) |
| BiSFET condensate | T_c,2D = min(T*, T_BKT) ≥ T_op | Nelson & Kosterlitz, PRL 39, 1201 (1977) |
| buildable double layer | d ≥ 0.6 nm, \|Eg\| ≤ 300 meV | ⚠ estimate of realistic heterostructures; not sourced |

## Results at 350 K

**Gap-switching transistor** (gate-tunable materials):

| material | T*max/2.2 | mean-field S_min | verdict |
|---|---|---|---|
| WTe₂ monolayer | 164 K | 118.8 mV/dec | worse than CMOS; no pairing at 350 K after the correction |
| MoSe₂/hBN/WSe₂ | 191 K | 66.0 mV/dec | just below CMOS (69.4), but only from the abrupt mean-field onset; no pairing after the correction, so an artefact |
| InAs/GaSb | 6 K | (metal) | gap < k_BT across the whole window: no switching |
| double bilayer graphene | 73 K | (metal) | same |

At 300 K, WTe₂ reaches about 40 mV/dec in mean field at the onset (see
`results/device_scoring.png` after `python device_scoring.py 300`). This is the same
artefact: T*/2.2 = 164 K is below 300 K.

Gap points below k_BT are blanked, because the activated-conduction formula
does not apply to a (semi)metal.

The bulk crystals (Ta₂NiSe₅, TiSe₂, Ta₂Pd₃Te₅) have no electrostatic Eg
control, so they are not scored. Real Ta₂NiSe₅ also fails the 350 K target:
its measured T_c is 326 K (Lu et al., Nat. Commun. 8, 14408, 2017).

**BiSFET**: every double layer **FAILS**. The 2D condensation temperature is
2–6 % of 350 K. Across all 1173 buildable double layers in the dataset, the
best is T_c,2D = 208 K, which is ≤ 95 K after the mean-field correction.
**None reaches 350 K, and none reaches 300 K.**

At 77 K, 93 of the 1173 condense. But the CMOS limit there is only 15 mV/dec,
and no candidate beats it.

## Conclusion

Within this model, no candidate beats CMOS at chip temperatures. The
bottleneck is **condensation (phase stiffness), not pairing**. Pairing scales
T* reach 300 K, but the 2D BKT bound keeps T_c at tens of kelvin, because pair
densities are around 10¹² cm⁻² and the exciton mass is around 1 m₀. This
matches the BEC-side picture in Randeria & Taylor, Fig. 3, where T_c stays far
below T*.

What would change the answer:
1. **Higher pair density with low exciton mass.** T_BKT ∝ n/M. The model shows
   no route that keeps the pairs bound at the same time.
2. **Steep switching near the EI transition.** It appears only in mean field.
   Fluctuations smear it out (R&T p.9 discusses the critical region), so a
   beyond-mean-field calculation or a measurement is needed.
3. **Measured inputs.** The material parameters are estimates (see
   REFERENCES.md), and the model underestimates the measured T_c of Ta₂NiSe₅
   by about 2×. Measured masses and gaps (ARPES, STM, optics) could
   move these numbers.

## References for this file
All are in REFERENCES.md with their verification status: IEA 2025, Ionescu & Riel 2011,
Banerjee et al. 2009, Nelson & Kosterlitz 1977, Randeria & Taylor 2014,
Lu et al. 2017, Unuchek et al. 2018 (room-temperature excitonic transistor, no condensate).
