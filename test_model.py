"""Sanity checks against known limits."""
from ei_model import ExcitonicInsulator2D
import numpy as np

# 1) 2D hydrogen: E_b = 4 Ry*
m = ExcitonicInsulator2D(m_r=0.2, eps=5.0)
Eb = m.exciton_binding_energy()
print(f"2D hydrogen  E_b = {Eb:.1f} meV   (exact 4 Ry* = {4*m.Ry:.1f})")
assert abs(Eb / (4 * m.Ry) - 1) < 0.03

# 2) condensate vanishes exactly at Eg = E_b (Keldysh-Kopaev criterion)
for Eg in (0.9 * Eb, 1.1 * Eb):
    D = m.solve_gap(Eg)
    print(f"Eg = {Eg:6.1f}  max Delta = {D.max():8.3f}")
assert m.solve_gap(0.9 * Eb).max() > 0 and m.solve_gap(1.1 * Eb).max() == 0

# 3) Keldysh screening and layer separation reduce E_b
m2 = ExcitonicInsulator2D(m_r=0.2, eps=5.0, r0=3.0)
m3 = ExcitonicInsulator2D(m_r=0.2, eps=5.0, d=1.0)
print(f"r0=3 nm: E_b = {m2.exciton_binding_energy():.1f}   d=1 nm: E_b = {m3.exciton_binding_energy():.1f}")
assert m2.exciton_binding_energy() < Eb and m3.exciton_binding_energy() < Eb
print("all checks passed")
