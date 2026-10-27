"""
Energy per switching event, E = C V^2.

Benchmark: Banerjee et al., IEEE EDL 30, 158 (2009) [PDF read in full]:
  gate capacitance C_g = 3.5e-6 F/cm^2 (p.159), gate length L = 10 nm (p.159),
  gate width W = 20 nm (p.160), clock swing 25 mV (p.160);
  SPICE result 0.008 aJ per clock cycle per BiSFET at 100 GHz (p.160).
A plain C_g W L V^2 estimate should land within a small factor of that
(SPICE also includes quantum capacitance, two gates and the NDR waveform).

For a gap-switching EI transistor the fair comparison is at equal geometry:
switching through N decades of current needs V = N S, so
    E_EI / E_CMOS = (S_EI / S_CMOS)^2.

Usage:  python energy.py
"""
import numpy as np
from ei_model import KB

LN10 = np.log(10)
C_G = 3.5e-6 * 1e4          # F/m^2   (3.5e-6 F/cm^2)
L, W = 10e-9, 20e-9         # m
C = C_G * L * W             # F

# measured swings: fluctuation-corrected minima from docs/FLUCTUATIONS.md (mV/dec)
S_FLUCT = {("MoSe2/hBN/WSe2", 350): (102.0, 140.2), ("WTe2 monolayer", 350): (146.7, 463.3),
           ("MoSe2/hBN/WSe2", 300): (103.5, 120.1), ("WTe2 monolayer", 300): (136.8, 279.5)}


def main():
    E_bisfet = C * 0.025 ** 2 / 1e-18
    print(f"device capacitance C = C_g W L = {C * 1e18:.1f} aF")
    print(f"BiSFET check: C V^2 at 25 mV = {E_bisfet:.4f} aJ   (Banerjee 2009 SPICE: 0.008 aJ)  "
          f"ratio {E_bisfet / 0.008:.2f}")
    decades = 4
    print(f"\nequal-geometry switching energy for {decades} decades of current, relative to ideal CMOS:")
    for (name, T), (lo, hi) in S_FLUCT.items():
        S_c = LN10 * KB * T
        r_lo, r_hi = (lo / S_c) ** 2, (hi / S_c) ** 2
        Vc = decades * S_c / 1e3
        print(f"  {name:16s} {T} K: E_EI/E_CMOS = {r_lo:4.1f} - {r_hi:5.1f}x   "
              f"(CMOS: {decades} x {S_c:.0f} mV = {Vc:.2f} V, {C * Vc ** 2 / 1e-18:.2f} aJ in this geometry)")


if __name__ == "__main__":
    main()
