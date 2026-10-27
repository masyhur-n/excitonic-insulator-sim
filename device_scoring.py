"""
Low-power device feasibility screen for excitonic insulators (EI).

This is a FEASIBILITY SCREEN, not a device simulation: no contacts, transport,
disorder or parasitics. Every assumption is listed in docs/DEVICE_SCORING.md.

Operating point: T_op = 350 K (chip hot-spot temperature, user choice).
Gate model: ideal, dEg = e dV (best case for a double layer).

1) Gap-switching transistor
   Conduction across the gap is thermally activated, G ~ exp(-E_gap / 2 k_B T)
   (intrinsic carriers, n_i ~ exp(-Eg/2kT): textbook, e.g. Sze, "Physics of
   Semiconductor Devices"). With dEg = e dV the subthreshold swing is
       S = 2 ln(10) k_B T / |dE_gap/dEg|          [mV/decade]
   compared with the thermionic limit S_CMOS = ln(10) k_B T / e
   (60 mV/dec at 300 K; Ionescu & Riel, Nature 479, 329 (2011)).
   E_gap is the mean-field spectral gap at T_op. It needs pairing (T* > T_op),
   not phase coherence: a pairing pseudogap between T_c and T* is enough
   (Randeria & Taylor 2014, p.10).

2) BiSFET (Banerjee, Register, Tutuc, Reddy, MacDonald, IEEE EDL 30, 158 (2009))
   Needs an interlayer exciton CONDENSATE at T_op. In 2D: T_c ~ min(T*, T_BKT)
   (Nelson & Kosterlitz, PRL 39, 1201 (1977)). Pass only if T_c >= T_op.
   Only defined for double layers (d > 0).

Usage:  python device_scoring.py [T_op]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ei_model import ExcitonicInsulator2D, KB
from materials import MATERIALS

LN10 = np.log(10)
MFT_OVERESTIMATE = 2.2   # (4e)^(1/3): MFT overestimates T_c even in the BCS limit (R&T 2014 p.10)


def gap_vs_Eg(model, Egs, T):
    """Mean-field spectral gap 2 min E(k) at temperature T along an Eg sweep."""
    gaps, Ds, D_prev = [], [], None
    for Eg in Egs[::-1]:
        D0 = model.solve_gap(Eg)                       # T = 0 solution as a seed
        D = model.solve_gap(Eg, T=T, delta0=D0 if D0.max() > 0 else None) if D0.max() > 0 else D0
        xi = model.xi(Eg)
        gaps.append(2 * np.min(np.sqrt(xi ** 2 + D ** 2)))
        Ds.append(D.max())
    return np.array(gaps[::-1]), np.array(Ds[::-1])


def score_material(name, p, T_op, n_eg=81):
    model = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"], nk=120)
    Egs = np.linspace(*p["Eg_range"], n_eg)
    gaps, Dmax = gap_vs_Eg(model, Egs, T_op)
    slope = np.abs(np.gradient(gaps, Egs))
    S = 2 * LN10 * KB * T_op / np.maximum(slope, 1e-9)
    # gap < k_B T: (semi)metallic, activated-conduction formula does not apply
    S = np.where(gaps < KB * T_op, np.nan, S)
    i = int(np.nanargmin(S)) if np.isfinite(S).any() else 0
    # resolution limit: a slope from one grid step can be a mean-field jump
    dEg = Egs[1] - Egs[0]
    # BiSFET: condensate at the nominal operating Eg
    D0 = model.solve_gap(p["Eg"])
    Tstar = model.critical_temperature(p["Eg"]) if D0.max() > 0 else 0.0
    Tstar_max = max(model.critical_temperature(e) for e in Egs[::8])  # best point in the gate window
    Tbkt = model.observables(p["Eg"], D0)["T_BKT"]
    Tc2d = min(Tstar, Tbkt)
    return dict(name=name, Egs=Egs, gaps=gaps, S=S, Dmax=Dmax,
                S_min=S[i] if np.isfinite(S).any() else np.nan, Eg_at_Smin=Egs[i], paired_at_Top=bool(Dmax.max() > 0),
                gateable=p["d"] > 0 or "monolayer" in name,
                bilayer=p["d"] > 0, Tstar=Tstar, Tstar_max=Tstar_max, Tbkt=Tbkt, Tc2d=Tc2d,
                grid_mV=dEg)


def bisfet_premise():
    """Re-derive the BiSFET paper's room-temperature numbers and compare with
    the screening critique. Numbers quoted from Banerjee et al., IEEE EDL 30,
    158 (2009), p.158: coherence onset ~0.1 E_F/k_B (their ref. [5], Min et al.,
    PRB 78, 121401 (2008)), reached at n = p = 4.9e12 cm^-2 with
    n = k_F^2/pi (spin x valley = 4) and v = 1e8 cm/s."""
    hbar_v = 658.212 * 1.0          # hbar * 1e6 m/s  [meV nm]
    n = 4.9e12 * 1e-14              # cm^-2 -> nm^-2
    kF = np.sqrt(np.pi * n)         # n = kF^2 / pi
    EF = hbar_v * kF
    print("\n3) BiSFET premise (monolayer graphene double layer, Dirac bands; not the parabolic model above)")
    print(f"  Banerjee 2009: n = 4.9e12 cm^-2 -> E_F = {EF:.0f} meV;  0.1 E_F/k_B = {0.1 * EF / KB:.0f} K"
          f"  (paper: ~300 K)")
    print(f"  Kharitonov & Efetov, PRB 78, 241401 (2008), with carrier screening: T_c ~ 1e-7 E_F"
          f" (as summarised in the reply by Bistritzer et al., arXiv:0810.0331)"
          f"\n    ="
          f" {1e-7 * EF / KB * 1e3:.3f} mK")
    print("  -> the room-temperature premise depends on unscreened mean field; contested in the literature")


def main(T_op=350.0):
    S_cmos = LN10 * KB * T_op * 1.0            # meV per decade == mV/dec for e*V
    print(f"T_op = {T_op:.0f} K   thermionic limit S_CMOS = {S_cmos:.1f} mV/dec"
          f"   (= 59.5 mV/dec at 300 K)\n")
    rows = [score_material(n, p, T_op) for n, p in MATERIALS.items()]

    print("1) gap-switching transistor (ideal gate, mean field)")
    print(f"  {'material':24s} {'gate?':6s} {'paired@T_op':11s} {'T*max/2.2':>9s} {'S_min':>9s} {'at Eg':>8s}  verdict")
    for r in rows:
        robust = r["Tstar_max"] / MFT_OVERESTIMATE
        if not r["gateable"]:
            print(f"  {r['name']:24s} {'no':6s} {'-':11s} {robust:8.0f}K {'-':>9s} {'-':>8s}"
                  f"  bulk crystal: no electrostatic Eg control (pressure/strain only)")
            continue
        if robust < T_op:
            verdict = (f"no pairing at T_op once the MFT overestimate is removed, so it acts as an "
                       f"ordinary gap-tuned device" + (" (MFT S_min is an artefact)" if r["S_min"] < S_cmos else ""))
        elif r["S_min"] < S_cmos:
            verdict = (f"beats CMOS on paper; grid step {r['grid_mV']:.1f} mV, so a steep "
                       f"mean-field onset (smeared by fluctuations, R&T p.9) is likely")
        else:
            verdict = "no gain over CMOS"
        if not np.isfinite(r["S_min"]):
            verdict = "gap < k_B T over the whole window: (semi)metal, no switching"
            s_txt, eg_txt = f"{'metal':>9s}", f"{'-':>8s}"
        else:
            s_txt, eg_txt = f"{r['S_min']:7.1f}mV", f"{r['Eg_at_Smin']:8.1f}"
        print(f"  {r['name']:24s} {'yes':6s} {str(r['paired_at_Top']):11s} {robust:8.0f}K"
              f" {s_txt} {eg_txt}  {verdict}")

    print(f"\n2) BiSFET: needs a 2D condensate at T_op = {T_op:.0f} K")
    print(f"  {'material':24s} {'T*':>7s} {'T_BKT<=':>8s} {'T_c,2D':>7s} {'margin':>7s}  verdict")
    for r in rows:
        if not r["bilayer"]:
            print(f"  {r['name']:24s} {'-':>7s} {'-':>8s} {'-':>7s} {'-':>7s}  not a double layer")
            continue
        m = r["Tc2d"] / T_op
        print(f"  {r['name']:24s} {r['Tstar']:6.0f}K {r['Tbkt']:7.0f}K {r['Tc2d']:6.0f}K {m:6.2f}x"
              f"  {'PASS' if m >= 1 else 'FAIL'}")

    # does ANY double layer in the model's parameter space reach T_op?
    data = np.genfromtxt("ei_dataset.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    Tbkt_all = np.pi * 38.0998 * data["density"] / (4 * data["m_r"] * KB)
    Tc_all = np.minimum(data["Tc"], Tbkt_all)
    # buildable double layers only: d >= 0.6 nm (closest TMD heterobilayer spacing, an
    # estimate), |Eg| <= 300 meV, and drop solver-ceiling values (T* = 2000 K marker)
    dl = (data["d"] >= 0.6) & (np.abs(data["Eg"]) <= 300) & (data["Tc"] < 2000)
    j = np.argmax(np.where(dl, Tc_all, -1))
    print(f"\n  best buildable double layer in the dataset: T_c,2D = {Tc_all[j]:.0f} K, "
          f"after the MFT correction <= {Tc_all[j] / MFT_OVERESTIMATE:.0f} K "
          f"(m_r={data['m_r'][j]:.2f}, eps={data['eps'][j]:.1f}, d={data['d'][j]:.2f} nm, Eg={data['Eg'][j]:.0f} meV)"
          f"  -> {'reaches' if Tc_all[j] >= T_op else 'below'} T_op")
    print(f"  buildable double layers with T_c,2D >= {T_op:.0f} K: {(Tc_all[dl] >= T_op).sum()} of {dl.sum()}")

    bisfet_premise()

    # plot
    gate = [r for r in rows if r["gateable"]]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for r in gate:
        ax[0].plot(r["Egs"], r["gaps"], label=r["name"])
        ax[1].semilogy(r["Egs"], r["S"], label=r["name"])
    ax[0].set(xlabel="Eg (= gate bias, meV)", ylabel="spectral gap at T_op [meV]",
              title=f"gap vs gate at {T_op:.0f} K (mean field)")
    ax[1].axhline(S_cmos, color="k", ls="--", label=f"CMOS limit {S_cmos:.0f} mV/dec")
    ax[1].set(xlabel="Eg (= gate bias, meV)", ylabel="S [mV/decade]", ylim=(10, 1e4),
              title="subthreshold swing of a gap-switching EI transistor")
    for a in ax:
        a.legend(fontsize=7)
    fig.tight_layout(); os.makedirs("results", exist_ok=True); fig.savefig("results/device_scoring.png", dpi=130)
    print("\nsaved results/device_scoring.png")


if __name__ == "__main__":
    main(*map(float, sys.argv[1:]))
