"""
Real-gate correction to the switching swing (gap 4: Hartree / interlayer charging).

Electrostatics of a dual-gated e-h bilayer, Zeng & MacDonald, PRB 102, 085154
(2020), arXiv:2006.09937 [ZM20], Eqs. (3), (15): with gates far away (d_e, d_h >> d_s)
    e V_bias - E_g = (phi_e + phi_h) + mu_x,   phi_e + phi_h = 4 pi e^2 n d_s / eps
so in the paired phase  d(eV) = (C + dmu_x/dn) dn,  C = 4 pi e^2 d_s/eps = 8 pi Ry* a_B* d_s.
The switching swing then uses dE_gap/d(eV) = (dE_gap/dn) / (C + dmu_x/dn)
instead of the ideal-gate dE_gap/dEg used in device_scoring.py.
Exciton-exciton term dmu_x/dn: this code (no intralayer exchange) and the ZM20
mean-field fit g = -12 Ry* a_B*^2 (with intralayer exchange) bracket it.
Benchmark: interlayer E_b at d = a_B vs ZM20 Fig. 2a (eps_b = 0.77 Ry* by linear fit).
At 300 K the pair gas is non-degenerate (k_B T_deg = pi Ry* n a_B^2 << k_B T for
n a_B^2 < 0.07), so the T = 0 compressibility must be supplemented by the thermal
(Boltzmann) term k_B T / n, which is always positive.
Units in Ry*, a_B*; device numbers for ZM20's MoSe2/WSe2 values Ry* = 110 meV.
"""
import numpy as np
from screening import Bilayer
from ei_model import KB

RY_MEV = 110.0     # ZM20: m* = 0.2 m_e, eps = 5 -> Ry* ~ 0.11 eV


def gap_and_mu(b, nx):
    rs = 1 / np.sqrt(np.pi * nx)
    mu, D = b.solve_fixed_n(rs)
    xi = b.xi(mu)
    return 2 * np.min(np.sqrt(xi ** 2 + D ** 2)) / b.Ry, 2 * mu / b.Ry


if __name__ == "__main__":
    b = Bilayer(d_aB=1.0)
    ns = np.geomspace(0.002, 0.04, 25)
    gaps, mux = np.array([gap_and_mu(b, n) for n in ns]).T
    dgap = np.gradient(gaps, ns); dmu = np.gradient(mux, ns)
    C = 8 * np.pi * 1.0                       # d_s = a_B
    kT = KB * 300 / RY_MEV
    S_cmos = np.log(10) * KB * 300
    print(f"d_s = a_B: C = {C:.1f} Ry* a_B^2;  g(this code) = {np.median(dmu):.0f}, g(ZM20) = -12")
    print(f"  n a_B^2   E_gap/Ry*   S_ideal   S_real T=0 compressibility (code / ZM20)   + thermal k_BT/n (code / ZM20)   [mV/dec, 300 K]")
    for n, g_, dg, dm in list(zip(ns, gaps, dgap, dmu))[::4]:
        dEg_ideal = -dg / dm                  # ideal gate: E_g <-> -mu_x
        S_ideal = 2 * np.log(10) * kT / max(abs(dEg_ideal), 1e-9) * RY_MEV
        S_code = 2 * np.log(10) * kT / max(abs(dg / (C + dm)), 1e-9) * RY_MEV
        S_zm = 2 * np.log(10) * kT / max(abs(dg / (C - 12)), 1e-9) * RY_MEV
        S_code_T = 2 * np.log(10) * kT / max(abs(dg / (C + dm + kT / n)), 1e-9) * RY_MEV
        S_zm_T = 2 * np.log(10) * kT / max(abs(dg / (C - 12 + kT / n)), 1e-9) * RY_MEV
        print(f"  {n:.4f}    {g_:.3f}     {S_ideal:7.0f}        {S_code:7.0f} / {S_zm:5.0f}                       {S_code_T:7.0f} / {S_zm_T:6.0f}")
    print(f"  CMOS limit at 300 K: {S_cmos:.1f} mV/dec")
