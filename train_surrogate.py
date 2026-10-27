"""Train ML surrogates on the solver dataset, then use them for inverse design.

1. Classifier: parameters -> phase (normal / BCS-EI / BEC-EI)
2. Regressors: parameters -> T_c, spectral gap.
   T_c target = min(T*, T_BKT): the mean-field T* is only a pairing scale on the
   BEC side (Randeria & Taylor 2014, p.9), and in 2D condensation is capped by
   the Nelson-Kosterlitz bound computed from the pair density.
3. Inverse design: search parameter space with the surrogate (~1e5 candidates
   in about a second) for a high T_c, then re-check the best ones with the full solver.

Usage:  python train_surrogate.py [ei_dataset.csv]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, r2_score, mean_absolute_error
from ei_model import ExcitonicInsulator2D, HB2_2M0, KB
from materials import MATERIALS

FEATURES = ["m_r", "eps", "r0", "d", "q_s", "Eg"]


def load(path):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    X = np.column_stack([data[f].astype(float) for f in FEATURES])
    return data, X


def feats(X):
    """Physics-informed features: hydrogenic scales help the model generalise."""
    m_r, eps, r0, d, q_s, Eg = X.T
    Ry = 13605.7 * m_r / eps ** 2
    aB = eps * 0.0529177 / m_r
    return np.column_stack([X, np.log(m_r), Ry, Eg / Ry, r0 / aB, d / aB, q_s * aB])


def main(path="ei_dataset.csv"):
    data, X = load(path)
    F = feats(X)
    y_cls = data["regime"]
    T_bkt = np.pi * HB2_2M0 * data["density"].astype(float) / (4 * data["m_r"].astype(float) * KB)
    y_tc = np.minimum(data["Tc"].astype(float), T_bkt)
    y_gap = data["qp_gap"].astype(float)
    idx_tr, idx_te = train_test_split(np.arange(len(X)), test_size=0.2, random_state=0)

    clf = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.08).fit(F[idx_tr], y_cls[idx_tr])
    acc = accuracy_score(y_cls[idx_te], clf.predict(F[idx_te]))
    reg_tc = HistGradientBoostingRegressor(max_iter=600, learning_rate=0.06).fit(F[idx_tr], np.log1p(y_tc[idx_tr]))
    reg_gap = HistGradientBoostingRegressor(max_iter=600, learning_rate=0.06).fit(F[idx_tr], np.log1p(y_gap[idx_tr]))
    tc_pred = np.expm1(reg_tc.predict(F[idx_te]))
    gap_pred = np.expm1(reg_gap.predict(F[idx_te]))
    print(f"phase classifier accuracy : {acc:.3f}")
    print(f"T_c   R2 = {r2_score(y_tc[idx_te], tc_pred):.3f}  MAE = {mean_absolute_error(y_tc[idx_te], tc_pred):.1f} K")
    print(f"gap   R2 = {r2_score(y_gap[idx_te], gap_pred):.3f}  MAE = {mean_absolute_error(y_gap[idx_te], gap_pred):.1f} meV")

    # --- screen the candidate materials: surrogate vs full solver vs experiment --
    # measured T_c only where the cited paper states it (see docs/REFERENCES.md)
    measured = {"Ta2NiSe5": "326 K (Lu 2017)", "1T-TiSe2": "~190 K (Kogar 2017)",
                "InAs/GaSb QW": "~10 K (Du 2017)", "MoSe2/hBN/WSe2": "signals >100 K (Wang 2019)"}
    lo, hi = X.min(axis=0), X.max(axis=0)
    print("\ncandidate materials at their nominal Eg:")
    print(f"  {'material':24s} {'phase':7s} {'surrogate':>9s} {'solver T*':>9s} {'T_BKT<=':>8s} {'T_c est':>8s}  measured")
    for name, p in MATERIALS.items():
        x = np.array([[p["m_r"], p["eps"], p["r0"], p["d"], 0.0, p["Eg"]]])
        outside = [f for f, v, a, b in zip(FEATURES, x[0], lo, hi) if not a <= v <= b]
        sur = float(np.expm1(reg_tc.predict(feats(x)))[0])
        mdl = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"])
        D = mdl.solve_gap(p["Eg"])
        ts = mdl.critical_temperature(p["Eg"]) if D.max() > 0 else 0.0
        tb = mdl.observables(p["Eg"], D)["T_BKT"]
        flag = f"  [surrogate extrapolating in {','.join(outside)}]" if outside else ""
        print(f"  {name:24s} {mdl.regime(p['Eg'], D):7s} {sur:8.0f}K {ts:8.0f}K {tb:7.0f}K {min(ts, tb):7.0f}K"
              f"  {measured.get(name, '-')}{flag}")
    print("  (bulk crystals Ta2NiSe5, TiSe2, Ta2Pd3Te5 are 3D: the 2D BKT cap does not apply, compare with T*)")

    # --- inverse design: highest-T_c hBN-encapsulated double layer ------------
    rng = np.random.default_rng(1)
    N = 200_000
    cand = np.column_stack([np.exp(rng.uniform(np.log(0.05), np.log(0.6), N)),  # m_r
                            np.full(N, 4.5),                                    # hBN
                            np.zeros(N),                                        # r0
                            rng.uniform(0.7, 3.0, N),                           # d >= 2 hBN layers
                            np.zeros(N),
                            rng.uniform(-50, 300, N)])                          # gate-tunable Eg
    tc = np.expm1(reg_tc.predict(feats(cand)))
    ph = clf.predict(feats(cand))
    ok = ph != "normal"
    best = np.argsort(np.where(ok, tc, -1))[::-1][:5]
    print("\ninverse design (hBN double layer, d >= 0.7 nm): top candidates, verified with full solver")
    for i in best:
        m_r, eps, r0, d, q_s, Eg = cand[i]
        mdl = ExcitonicInsulator2D(m_r=m_r, eps=eps, r0=r0, d=d)
        D = mdl.solve_gap(Eg)
        tc_true = min(mdl.critical_temperature(Eg), mdl.observables(Eg, D)["T_BKT"]) if D.max() > 0 else 0.0
        print(f"  m_r={m_r:.3f} d={d:.2f} nm Eg={Eg:6.1f} meV | surrogate T_c={tc[i]:6.0f} K"
              f" | solver T_c={tc_true:6.0f} K ({mdl.regime(Eg, D)})")

    # --- parity plot ----------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
    ax[0].scatter(y_tc[idx_te], tc_pred, s=6, alpha=0.5)
    lim = [0, max(y_tc[idx_te].max(), tc_pred.max())]
    ax[0].plot(lim, lim, "k--", lw=0.8); ax[0].set(xlabel="solver T_c [K]", ylabel="surrogate T_c [K]", title="T_c surrogate (test set)")
    ax[1].scatter(y_gap[idx_te], gap_pred, s=6, alpha=0.5, color="C1")
    lim = [0, max(y_gap[idx_te].max(), gap_pred.max())]
    ax[1].plot(lim, lim, "k--", lw=0.8); ax[1].set(xlabel="solver gap [meV]", ylabel="surrogate gap [meV]", title="gap surrogate (test set)")
    fig.tight_layout(); os.makedirs("results", exist_ok=True); fig.savefig("results/surrogate_parity.png", dpi=130)
    print("\nsaved results/surrogate_parity.png")


if __name__ == "__main__":
    main(*sys.argv[1:])
