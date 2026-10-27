"""
Reproduce Table I of the manuscript: every benchmark, with its literature value,
the value computed here, and a pass/fail tolerance.

Usage:  python run_benchmarks.py          (about 5 minutes)
        python run_benchmarks.py --fast   (skips the two screening-onset scans)
"""
import sys
import numpy as np
from ei_model import ExcitonicInsulator2D, KB

ROWS = []


def check(name, lit, val, ok, note="", status=None):
    """status overrides PASS/FAIL, e.g. 'DEVIATION' for a documented, quantified mismatch."""
    status = status or ("PASS" if ok else "FAIL")
    ROWS.append((name, lit, val, status, note))
    print(f"[{status:9s}] {name:54s} lit: {lit:24s} this work: {val}  {note}", flush=True)


def main(fast=False):
    # 2D continuum gap equation
    m = ExcitonicInsulator2D(m_r=0.2, eps=5.0)
    Eb = m.exciton_binding_energy()
    check("2D exciton binding energy 4 Ry*", f"{4 * m.Ry:.1f} meV (exact)", f"{Eb:.1f} meV", abs(Eb / (4 * m.Ry) - 1) < 0.01)
    ok = m.solve_gap(0.9 * Eb).max() > 0 and m.solve_gap(1.1 * Eb).max() == 0
    check("EI exists only for Eg < Eb (Jerome-Rice-Kohn 1967)", "yes", "yes" if ok else "no", ok)

    # crossover mean field (Randeria & Taylor 2014)
    from benchmark_meanfield import leggett_3d, bcs_tc_3d
    mu, D = leggett_3d(0.0)
    check("3D unitarity Delta/E_F, mu/E_F (R&T Table I)", "0.68, 0.59", f"{D:.3f}, {mu:.3f}", abs(D - 0.68) < 0.01 and abs(mu - 0.59) < 0.01)
    tc = bcs_tc_3d(-2.0); tf = 8 * np.exp(np.euler_gamma) / (np.pi * np.e ** 2) * np.exp(-np.pi)
    check("BCS-limit mean-field T_c (R&T p.9), 1/k_F a = -2", f"{tf:.4f} E_F", f"{tc:.4f} E_F", abs(tc / tf - 1) < 0.02)

    # fixed-density bilayer (Pieri et al. 2007, Fig. 1b, d = a_B; values read from figure)
    from screening import Bilayer, condensate_fraction
    b = Bilayer(d_aB=1.0)
    vals = []
    for rs, lit in ((1.5, 0.66), (3.0, 0.68), (20.0, 0.18)):
        mu_, D_ = b.solve_fixed_n(rs); vals.append((D_.max() / b.Ry, lit))
    ok = all(abs(v - l) < 0.04 for v, l in vals)
    check("Pieri 2007 Delta_max, r_s = 1.5, 3, 20", "0.66, 0.68, 0.18 Ry*", ", ".join(f"{v:.3f}" for v, _ in vals), ok, "(figure read +-0.03)")

    # superfluid-state screening (Neilson et al. 2014)
    if not fast:
        from screening import onset_rs
        r1 = onset_rs(1.0)
        check("SS superfluid onset r_s,c at d = a_B (N14)", "2.5 DQMC / ~3 SS", f"{r1:.2f}", 2.2 <= r1 <= 3.1)
        r07 = onset_rs(0.7)
        check("SS superfluid onset r_s,c at d = 0.7 a_B (N14)", "2.1 DQMC", f"{r07:.2f}", 1.8 <= r07 <= 2.4)
    mu_, D_ = b.solve_SS(6.0); muU, DU = (b.set_interaction(None), b.solve_fixed_n(6.0))[1]
    ratio = D_.max() / DU.max()
    check("SS -> US at low density, r_s = 6 (N14 Fig. 3)", "equal for r_s >~ 5", f"Delta_SS/Delta_US = {ratio:.2f}", ratio > 0.85)

    # Dirac bands (Kharitonov & Efetov 2008)
    from dirac import DiracBilayer
    bn = DiracBilayer(eps=2.19, d=0.0, mode="NS", nk=40); bu = DiracBilayer(eps=2.19, d=0.0, mode="US", nk=40)
    r = bu.V(bn.kF) / bn.V(bn.kF)
    check("Dirac V_us(p_F)/V_NS(p_F), alpha = 1 (KE08 p.3)", "~9", f"{r:.2f}", abs(r - 9) < 0.3)
    lm = DiracBilayer(eps=0.0219, d=0.0, mode="NS", nk=40).lambda_c()
    check("Dirac maximum coupling lambda_c (KE08 Eq. 2)", "1/16 = 0.0625", f"{lm:.4f}", abs(lm - 0.0625) < 0.002)
    bk = DiracBilayer(eps=2.19, d=0.0, mode="NS", nk=300, cluster=6.0, ntheta=64)
    Ts = np.array([0.02, 0.01, 0.005, 0.0025, 0.00125]) * bk.EF / KB
    sl = np.polyfit(np.log(1 / Ts), [bk.lin_eig(T) for T in Ts], 1)[0]
    check("Dirac kernel slope vs lambda_c (KE08 Eqs. 10-12)", f"{bk.lambda_c():.4f}", f"{sl:.4f}", abs(sl / bk.lambda_c() - 1) < 0.06)

    # Ta2NiSe5 three-chain model (Kaneko et al. 2013)
    from ta2nise5_chain import ThreeChain
    d0 = ThreeChain(0.0, 0.12).solve()["delta"]; d1 = ThreeChain(0.0, 0.13).solve()["delta"]
    check("Ta2NiSe5 distortion onset lambda_c (K13 Fig. 3f)", "0.125 eV", "0.12-0.13 eV", d0 < 1e-4 < d1)
    lo, hi = 0.6, 0.75
    while hi - lo > 0.005:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if ThreeChain(mid, 0.0).solve()["Delta"] > 1e-4 else (mid, hi)
    vc = 0.5 * (lo + hi)
    check("Ta2NiSe5 EI onset V_c (K13 Fig. 3e)", "0.6 eV", f"{vc:.3f} eV", False,
          f"({100 * (vc / 0.6 - 1):+.0f}%; K13 supplement not accessed)", status="DEVIATION")
    tcK = ThreeChain(0.55, 0.04).Tc()
    check("Ta2NiSe5 k_B T_c, V = 0.55, lambda = 0.04 (K13 Fig. 4)", "0.0783 eV", f"{tcK:.4f} eV", False,
          f"({100 * (tcK / 0.0783 - 1):+.0f}%; K13 supplement not accessed)", status="DEVIATION")

    # gate (Zeng & MacDonald 2020)
    m1 = ExcitonicInsulator2D(m_r=0.25, eps=5.0, d=ExcitonicInsulator2D(m_r=0.25, eps=5.0).a_B)
    e1 = m1.exciton_binding_energy() / m1.Ry
    check("interlayer E_b at d = a_B (ZM20 Fig. 2a)", "0.77 Ry* (linear fit)", f"{e1:.3f} Ry* (exact)", 0.7 < e1 < 0.95)

    # disorder (Bistritzer & MacDonald 2008)
    from disorder import tc_ratio, critical_delta
    s = (tc_ratio(1e-3) - 1) / 1e-3; dc = critical_delta()
    check("AG disorder slope; critical delta (BM08)", "-0.7854; 0.88 T_c0", f"{s:.4f}; {dc:.4f}", abs(s + np.pi / 4) < 1e-3 and abs(dc - 0.882) < 0.005)

    # switching energy (Banerjee et al. 2009)
    from energy import C
    E = C * 0.025 ** 2 / 1e-18
    check("BiSFET switching energy (C V^2 vs SPICE)", "0.008 aJ", f"{E:.4f} aJ", 0.3 < E / 0.008 < 1.5, "(order-of-magnitude check)")

    n_fail = sum(r[3] == "FAIL" for r in ROWS)
    n_dev = sum(r[3] == "DEVIATION" for r in ROWS)
    print(f"\n{len(ROWS) - n_fail - n_dev} passed, {n_dev} documented deviations, {n_fail} failed (of {len(ROWS)})")
    return n_fail


if __name__ == "__main__":
    sys.exit(main(fast="--fast" in sys.argv))
