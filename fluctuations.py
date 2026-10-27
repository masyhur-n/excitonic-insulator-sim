"""
Gaussian order-parameter fluctuations beyond mean field: does the steep
mean-field gap onset (the only sub-CMOS switching found in device_scoring.py)
survive?

Method
------
1. Landau function from the Gibbs-Bogoliubov variational free energy
   (a true upper bound for any trial gap Delta_k; stationary at the gap equation):
       Omega(Delta) = sum_k [xi - E - 2T ln(1 + e^{-E/T}) + 2 Delta F]
                      - sum_kk' F_k V_kk' F_k',      F_k = Delta_k tanh(E/2T) / 2E
   evaluated along Delta_k = A phi_k (phi = mean-field shape, or the leading
   eigenvector of the linearised gap equation in the normal state).
2. Ginzburg-Landau (GL) form for the pair field psi, |psi|^2 = s A^2
   (s from the pair amplitude sum_k F_k^2 at small A):
       f = a|psi|^2 + (b/2)|psi|^4 + c|grad psi|^2,     c = hbar^2 / 2M,  M = m_e + m_h
   (c is exact in the BEC limit; on the BCS side it overestimates the stiffness,
   which UNDERESTIMATES fluctuations, so this test is conservative).
3. Self-consistent Gaussian (Hartree / large-N) treatment of the quartic term
   [Ullah & Dorsey, PRB 44, 262 (1991); Larkin & Varlamov, Theory of Fluctuations
   in Superconductors (OUP 2005)], with Bose occupation of the pair modes
   (quantum, so it stays correct when a >> k_B T, i.e. a dilute thermal exciton gas):
       m = <|psi|^2> = sum_{Q<Lambda} n_B(a_eff + c Q^2)
         = (T / 4 pi c) [ln(1 - e^{-(a_eff + c Lambda^2)/T}) - ln(1 - e^{-a_eff/T})],
       a_eff = a + b m
   Near the transition (a_eff << T) this is the classical GL log; at T -> 0 it
   gives mean field (m = -a/b). Gaussian fluctuations about the
   normal state are equivalent to NSR (Randeria & Taylor 2014, p.9).
   Cut-off Lambda = 1/a_B (pair size scale; GL valid for Q xi < 1), constant per
   material so it introduces no jumps; sensitivity checked x0.5, x2.
4. Spectral (pseudo)gap  E_gap = 2 min_k sqrt(xi_k^2 + <A^2> phi_k^2),
   then the same swing S = 2 ln10 k_B T / |dE_gap/dEg| as device_scoring.py.

Usage:  python fluctuations.py [T_op]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.linalg import eig
from ei_model import ExcitonicInsulator2D, KB, HB2_2M0
from materials import MATERIALS

LN10 = np.log(10)


def pair_shape(model, Eg, T):
    """Mean-field gap shape if ordered, else leading eigenvector of the
    linearised gap equation (normal-state pair channel)."""
    D = model.solve_gap(Eg)
    if D.max() > 0:
        DT = model.solve_gap(Eg, T=T, delta0=D)
        if DT.max() > 0:
            return DT / DT.max(), DT.max()
    xi = model.xi(Eg)
    b = 1.0 / (KB * T)
    f = np.where(np.abs(xi) * b < 1e-8, b / 4, np.tanh(0.5 * b * xi) / (2 * xi))
    w, v = eig(model.M * f[None, :])
    phi = np.abs(np.real(v[:, np.argmax(np.real(w))]))
    return phi / phi.max(), 0.0


def omega(model, Eg, T, Delta):
    """Gibbs-Bogoliubov variational free energy per area [meV nm^-2], minus normal state."""
    xi = model.xi(Eg)
    meas = model.k * model.w / (2 * np.pi)
    kT = KB * T

    def om(D):
        E = np.sqrt(xi ** 2 + D ** 2)
        F = D * np.tanh(E / (2 * kT)) / (2 * E)
        single = xi - E - 2 * kT * np.logaddexp(0, -E / kT)
        return np.sum(meas * (single + 2 * D * F)) - np.sum(meas * F * (model.M @ F))

    return om(Delta) - om(np.zeros_like(Delta))


FIT_KT = 1.0   # thermal part of the fit range, in units of k_B T


def gl_coefficients(model, Eg, T, phi, A_mf):
    xi = model.xi(Eg)
    meas = model.k * model.w / (2 * np.pi)
    kT = KB * T
    A_s = np.hypot(1.1 * A_mf, FIT_KT * kT)           # continuous across the MF transition
    As = np.linspace(0.05, 1.0, 14) * A_s
    dO = np.array([omega(model, Eg, T, A * phi) for A in As])
    # quartic GL fit up to just past the mean-field minimum
    X = np.column_stack([As ** 2, As ** 4 / 2])
    alpha, beta = np.linalg.lstsq(X, dO, rcond=None)[0]
    E0 = np.maximum(np.abs(xi), 1e-9)
    s = np.sum(meas * (phi * np.tanh(E0 / (2 * kT)) / (2 * E0)) ** 2)   # |psi|^2 = s A^2
    return alpha / s, beta / s ** 2, s


def fluctuation_gap(model, Eg, T, cut=1.0):
    phi, A_mf = pair_shape(model, Eg, T)
    a, b, s = gl_coefficients(model, Eg, T, phi, A_mf)
    c = HB2_2M0 / (4 * model.m_r)                    # hbar^2 / 2M, M = 4 m_r
    Lam = cut / model.a_B
    kT = KB * T
    b = max(b, 1e-12)                                 # quartic must stabilise

    def g(m):
        aeff = a + b * m
        return m - kT / (4 * np.pi * c) * (np.log(-np.expm1(-(aeff + c * Lam ** 2) / kT))
                                            - np.log(-np.expm1(-aeff / kT)))

    lo = max(0.0, -a / b) + 1e-14                     # a_eff > 0
    hi = lo + 1.0
    while g(hi) < 0:
        hi *= 2
    m = brentq(g, lo + 1e-16, hi)
    A2 = m / s
    xi = model.xi(Eg)
    gap_fl = 2 * np.min(np.sqrt(xi ** 2 + A2 * phi ** 2))
    gap_mf = 2 * np.min(np.sqrt(xi ** 2 + (A_mf * phi) ** 2))
    return gap_fl, gap_mf, -a / b / s if a < 0 else 0.0, A_mf ** 2


def swing(gaps, Egs, T):
    S = 2 * LN10 * KB * T / np.maximum(np.abs(np.gradient(gaps, Egs)), 1e-9)
    return np.where(gaps < KB * T, np.nan, S)


def main(T_op=350.0):
    S_cmos = LN10 * KB * T_op
    targets = ["MoSe2/hBN/WSe2", "WTe2 monolayer"]
    print(f"T_op = {T_op:.0f} K, CMOS limit {S_cmos:.1f} mV/dec\n")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for name in targets:
        p = MATERIALS[name]
        model = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"], nk=120)
        Egs = np.linspace(*p["Eg_range"], int(sys.argv[2]) if len(sys.argv) > 2 else 121)
        res = {cut: np.array([fluctuation_gap(model, e, T_op, cut) for e in Egs]) for cut in (0.5, 1.0, 2.0)}
        r = res[1.0]
        # consistency: GL quartic reproduces the mean-field amplitude
        ok = r[:, 3] > 0
        ratio = np.median(r[ok, 2] / r[ok, 3]) if ok.any() else np.nan
        S_mf, S_fl = swing(r[:, 1], Egs, T_op), swing(r[:, 0], Egs, T_op)
        print(f"{name}:  GL check  A^2(GL)/A^2(MF) = {ratio:.2f} (median, ordered points; 1 = consistent)")
        print(f"  mean field      S_min = {np.nanmin(S_mf):6.1f} mV/dec")
        for cut in (0.5, 1.0, 2.0):
            Sf = swing(res[cut][:, 0], Egs, T_op)
            print(f"  + fluctuations  S_min = {np.nanmin(Sf):6.1f} mV/dec   (cut-off {cut} / a_B)")
        ax[0].plot(Egs, r[:, 1], "--", label=f"{name} mean field")
        ax[0].plot(Egs, r[:, 0], label=f"{name} + fluctuations")
        ax[1].semilogy(Egs, S_mf, "--", label=f"{name} mean field")
        ax[1].semilogy(Egs, S_fl, label=f"{name} + fluctuations")
    ax[1].axhline(S_cmos, color="k", ls=":", label=f"CMOS {S_cmos:.0f} mV/dec")
    ax[0].set(xlabel="Eg (= gate bias, meV)", ylabel="spectral gap [meV]", title=f"gap at {T_op:.0f} K")
    ax[1].set(xlabel="Eg (= gate bias, meV)", ylabel="S [mV/dec]", ylim=(10, 1e4), title="subthreshold swing")
    for a in ax:
        a.legend(fontsize=7)
    fig.tight_layout(); os.makedirs("results", exist_ok=True); fig.savefig(f"results/fluctuations_{T_op:.0f}K.png", dpi=130)
    print(f"\nsaved results/fluctuations_{T_op:.0f}K.png")


if __name__ == "__main__":
    main(*map(float, sys.argv[1:2]))
