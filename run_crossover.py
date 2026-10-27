"""Sweep Eg across the BEC-BCS crossover and plot the main observables.

Usage:  python run_crossover.py [material]     (default: MoSe2/hBN/WSe2)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ei_model import ExcitonicInsulator2D
from materials import MATERIALS


def sweep(p, Egs, q_s=0.0):
    model = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"], q_s=q_s)
    rows, D_prev, profiles = [], None, {}
    for Eg in Egs[::-1]:                    # start on the insulator side, reuse guesses
        D = model.solve_gap(Eg, delta0=D_prev if D_prev is not None and D_prev.max() > 0 else None)
        D_prev = D
        obs = model.observables(Eg, D)
        obs.update(Eg=Eg, regime=model.regime(Eg, D),
                   Tc=model.critical_temperature(Eg) if D.max() > 0 else 0.0)
        rows.append(obs)
        profiles[Eg] = (D.copy(), model.xi(Eg))
    return model, rows[::-1], profiles


def main(name="MoSe2/hBN/WSe2"):
    p = MATERIALS[name]
    ref = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"])
    Eb = ref.exciton_binding_energy()
    print(f"{name}: E_b = {Eb:.1f} meV, a_B = {ref.a_B:.2f} nm")
    Egs = np.linspace(-3 * Eb, 1.15 * Eb, 56)

    model, bare, prof = sweep(p, Egs)
    # free carriers exist only in the semimetal: screen only for Eg < 0
    qs = 0.5 / ref.a_B
    _, scr, _ = sweep(p, Egs[Egs < 0], q_s=qs)

    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    for rows, ls, tag in ((bare, "-", "unscreened"), (scr, "--", f"screened q_s={qs:.2f}/nm")):
        g = lambda key: np.array([r[key] for r in rows])
        x = g("Eg") / Eb
        ax[0, 0].plot(x, g("qp_gap"), ls, label=f"gap 2·min E ({tag})")
        ax[0, 1].plot(x, g("Tc"), ls, label=f"mean-field T* ({tag})")
        ax[0, 1].plot(x, np.minimum(g("Tc"), g("T_BKT")), ls, lw=2.5,
                      label=f"T_c est. = min(T*, T_BKT) ({tag})")
        ax[0, 2].semilogy(x, np.where(g("size_ratio") > 0, g("size_ratio"), np.nan), ls, label=tag)
        ax[1, 0].semilogy(x, np.where(g("density") > 0, g("density"), np.nan), ls, label=tag)
    ax[0, 0].plot(np.array([r["Eg"] for r in bare]) / Eb, [r["Delta0"] for r in bare], ":", label="max Δ(k) (unscreened)")
    ax[0, 0].set(ylabel="meV", title="spectral gap & order parameter")
    ax[0, 1].set(ylabel="K", title="pairing scale $T^*$ vs condensation $T_c$", yscale="log")
    ax[0, 2].set(title="pair size / spacing  ξ/r$_s$   (>1 BCS, <1 BEC)")
    ax[0, 2].axhline(1, color="gray", lw=0.8)
    ax[1, 0].set(ylabel="nm$^{-2}$", title="e-h pair density")
    for a in ax.flat[:4]:
        a.axvline(0, color="k", lw=0.8, ls=":"); a.axvline(1, color="r", lw=0.8, ls=":")
        a.set_xlabel("$E_g / E_b$   (BCS ← 0 → BEC,  1 = insulator edge)")
        a.legend(fontsize=7)
    for target in (-2.5, 0.0, 0.8):
        Eg = min(prof, key=lambda e: abs(e / Eb - target))
        D, xi = prof[Eg]
        E = np.sqrt(xi ** 2 + D ** 2)
        ax[1, 1].plot(model.k * model.a_B, D, label=f"Eg = {Eg/Eb:+.1f} E_b")
        ax[1, 2].plot(model.k * model.a_B, 0.5 * (1 - xi / E), label=f"Eg = {Eg/Eb:+.1f} E_b")
    ax[1, 1].set(xlim=(0, 6), xlabel="k a$_B$", ylabel="meV", title="Δ(k)  (unscreened)"); ax[1, 1].legend()
    ax[1, 2].set(xlim=(0, 6), xlabel="k a$_B$", title="pair occupation $v_k^2$: Fermi edge (BCS) vs smooth (BEC)")
    ax[1, 2].legend()
    fig.suptitle(f"{name}: 2D excitonic insulator BEC–BCS crossover  "
                 f"(m_r={p['m_r']}, ε={p['eps']}, r0={p['r0']} nm, d={p['d']} nm, E_b={Eb:.0f} meV)")
    fig.tight_layout()
    os.makedirs("results", exist_ok=True)
    out = "results/crossover_" + name.replace("/", "-").replace(" ", "_") + ".png"
    fig.savefig(out, dpi=130)
    print("saved", out)
    for r in bare[::5]:
        print(f"Eg={r['Eg']:7.1f}  {r['regime']:6s} Δ0={r['Delta0']:7.2f}  gap={r['qp_gap']:7.2f}  "
              f"T*={r['Tc']:7.1f}  T_BKT<={r['T_BKT']:7.1f}  ξ/rs={r['size_ratio']:.2f}")


if __name__ == "__main__":
    main(*sys.argv[1:])
