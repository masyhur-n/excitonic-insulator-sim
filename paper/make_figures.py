"""Publication figures for the manuscript (run from the repo root: python paper/make_figures.py)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ei_model import ExcitonicInsulator2D, KB
from tc_ceiling import X as XS, RSC, ceiling
import fluctuations as fl
from materials import MATERIALS
from ta2nise5_chain import ThreeChain

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5, "figure.dpi": 300,
                     "lines.linewidth": 1.2, "axes.linewidth": 0.6})
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(OUT, exist_ok=True)


def fig1():
    p = MATERIALS["MoSe2/hBN/WSe2"]
    m = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"])
    Eb = m.exciton_binding_energy()
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.2))
    for x, c in ((-2.5, "C0"), (0.0, "C1"), (0.8, "C2")):
        D = m.solve_gap(x * Eb); xi = m.xi(x * Eb); E = np.sqrt(xi ** 2 + D ** 2)
        ax[0].plot(m.k * m.a_B, D / Eb, color=c, label=f"$E_g/E_b={x:+.1f}$")
        ax[1].plot(m.k * m.a_B, 0.5 * (1 - xi / E), color=c)
    Egs = np.linspace(-3, 1.15, 40) * Eb
    Ts, gaps = [], []
    Dp = None
    for Eg in Egs[::-1]:
        D = m.solve_gap(Eg, delta0=Dp if Dp is not None and Dp.max() > 0 else None); Dp = D
        Ts.append(m.critical_temperature(Eg) if D.max() > 0 else 0.0)
        gaps.append(m.observables(Eg, D)["T_BKT"])
    Ts, Tb = np.array(Ts[::-1]), np.array(gaps[::-1])
    ax[2].semilogy(Egs / Eb, np.where(Ts > 0, Ts, np.nan), "k", label="$T^*$ (mean field)")
    ax[2].semilogy(Egs / Eb, np.where(Ts > 0, np.minimum(Ts, Tb), np.nan), "C3", label=r"min$(T^*,T_{\rm BKT})$")
    ax[2].axhline(350, color="gray", ls=":", lw=0.8); ax[2].text(0.1, 250, "350 K", fontsize=6, color="gray")
    for a in ax[:2]:
        a.set_xlim(0, 5); a.set_xlabel("$k a_B$")
    ax[0].set_ylabel(r"$\Delta_k / E_b$"); ax[1].set_ylabel("$v_k^2$"); ax[0].legend(frameon=False)
    ax[2].set_xlabel("$E_g/E_b$"); ax[2].set_ylabel("T (K)"); ax[2].legend(frameon=False, loc="lower left")
    for a, t in zip(ax, "abc"):
        a.text(0.03, 0.92, f"({t})", transform=a.transAxes, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig1_crossover.pdf"); plt.close(fig)


def fig2():
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax[0].plot(RSC, XS, "o-", color="C3", ms=3, label="this work (SS screening)")
    ax[0].plot([2.1, 2.5], [0.7, 1.0], "ks", ms=4, mfc="none", label="DQMC [Neilson 2014, Fig. 1a]")
    ax[0].plot([3.0], [1.0], "^", color="C0", ms=4, mfc="none", label="SS [Neilson 2014]")
    ax[0].set_xlabel("$r_s$"); ax[0].set_ylabel("$d/a_B$"); ax[0].legend(frameon=False)
    ax[0].text(1.4, 2.5, "normal", fontsize=7); ax[0].text(3.9, 1.0, "superfluid", fontsize=7)
    e = np.linspace(2, 12, 120); d = np.linspace(0.6, 3.0, 120)
    E, Dd = np.meshgrid(e, d)
    T, _ = ceiling(0.22, E, Dd)
    cs = ax[1].contourf(E, Dd, T, levels=[0, 25, 50, 100, 150, 200, 300, 350], cmap="viridis")
    ax[1].contour(E, Dd, T, levels=[300], colors="w", linewidths=0.8)
    ax[1].axvline(4.5, color="w", ls="--", lw=0.8); ax[1].text(4.6, 2.7, "hBN", color="w", fontsize=7)
    ax[1].plot(4.5, 1.6, "r*", ms=7); ax[1].text(4.8, 1.55, "MoSe$_2$/hBN/WSe$_2$", color="r", fontsize=6)
    ax[1].set_xlabel(r"$\epsilon$"); ax[1].set_ylabel("$d$ (nm)")
    plt.colorbar(cs, ax=ax[1], label=r"$T_c$ ceiling (K), $m_r=0.22$")
    for a, t in zip(ax, "ab"):
        a.text(0.03, 0.92, f"({t})", transform=a.transAxes, fontweight="bold", color="k" if t == "a" else "w")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig2_screening_ceiling.pdf"); plt.close(fig)


def fig3():
    T = 350.0
    p = MATERIALS["MoSe2/hBN/WSe2"]
    model = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"], nk=120)
    Egs = np.linspace(*p["Eg_range"], 121)
    r = np.array([fl.fluctuation_gap(model, e, T) for e in Egs])
    S_mf, S_fl = fl.swing(r[:, 1], Egs, T), fl.swing(r[:, 0], Egs, T)
    S_c = np.log(10) * KB * T
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.4))
    ax[0].plot(Egs, r[:, 1], "--", color="C0", label="mean field")
    ax[0].plot(Egs, r[:, 0], color="C3", label="+ Gaussian fluctuations")
    ax[0].set_xlabel("$E_g$ (meV) = gate bias"); ax[0].set_ylabel("spectral gap (meV)"); ax[0].legend(frameon=False)
    ax[1].semilogy(Egs, S_mf, "--", color="C0", label="mean field")
    ax[1].semilogy(Egs, S_fl, color="C3", label="+ fluctuations")
    ax[1].axhline(S_c, color="k", ls=":", label=f"CMOS limit {S_c:.0f} mV/dec")
    ax[1].set_ylim(20, 1e4); ax[1].set_xlabel("$E_g$ (meV)"); ax[1].set_ylabel("$S$ (mV/decade)")
    ax[1].legend(frameon=False, loc="upper right")
    for a, t in zip(ax, "ab"):
        a.text(0.03, 0.92, f"({t})", transform=a.transAxes, fontweight="bold")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig3_switching.pdf"); plt.close(fig)


def fig4():
    m = ThreeChain(0.55, 0.04)
    Ts = np.linspace(0.002, 0.065, 22)
    D, dl, st = [], [], None
    for t in Ts:
        r = m.solve(t, state=st); st = r["state"]
        D.append(0.55 * abs(st[0] + st[1])); dl.append(r["delta"])
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.4))
    ax.plot(Ts / 8.617e-5, D, "C3", label=r"excitonic $|V(\phi_0+\phi_1)|$")
    ax.plot(Ts / 8.617e-5, dl, "C0", label=r"lattice $|\delta|$")
    ax.axvline(326 * 2.2, color="gray", ls=":", lw=0.8)
    ax.text(326 * 2.2 + 10, 0.04, r"$2.2\times326$ K", fontsize=6, color="gray")
    ax.set_xlabel("T (K, mean field)"); ax.set_ylabel("hybridisation (eV)"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig4_ta2nise5.pdf"); plt.close(fig)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4):
        f(); print("done", f.__name__, flush=True)
