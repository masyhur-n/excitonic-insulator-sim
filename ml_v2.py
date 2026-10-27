"""
ML surrogate v2: retrained on the upgraded physics, with uncertainty.

Label (2D double layers):
  T_c = min( T*_MF / 2.2,                          # mean-field pairing scale, R&T 2014 p.10 correction
             T_BKT bound with M = m_e + m_h,       # Nelson-Kosterlitz, n from the solver
             SS-screening ceiling Ry*/(4 r_s,c^2) * 4 m_r / M )   # tc_ceiling.py (Neilson 2014 benchmark)
Features: m_r, mass ratio m_h/m_e, eps, r0, d, Eg (wider box than v1).
Uncertainty: bootstrap ensemble of gradient-boosted trees (mean +- std); inputs
outside the training box are flagged.
Usage:  python ml_v2.py [n_samples]
"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np
from multiprocessing import Pool
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error
from ei_model import ExcitonicInsulator2D, HB2_2M0, KB
from tc_ceiling import rs_c

FEATS = ["m_r", "mass_ratio", "eps", "r0", "d", "Eg"]


def label(p):
    model = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"], nk=100, ntheta=32)
    D = model.solve_gap(p["Eg"])
    if D.max() == 0:
        return 0.0
    Ts = model.critical_temperature(p["Eg"])
    Ts = 2000.0 if not np.isfinite(Ts) else Ts
    n = model.observables(p["Eg"], D)["density"]
    q = p["mass_ratio"]                                   # m_h / m_e
    M_over_4mr = (1 + q) ** 2 / (4 * q)                   # (m_e+m_h)^2 / (4 m_e m_h)
    T_bkt = np.pi * HB2_2M0 * n / (4 * p["m_r"] * M_over_4mr * KB)
    x = max(p["d"] / model.a_B, 0.3)
    T_scr = model.Ry / (4 * float(rs_c(x)) ** 2) / KB / M_over_4mr
    return float(min(Ts / 2.2, T_bkt, T_scr))


def sample(seed):
    rng = np.random.default_rng(seed)
    p = dict(m_r=float(np.exp(rng.uniform(np.log(0.01), np.log(1.0)))),
             mass_ratio=float(np.exp(rng.uniform(np.log(0.2), np.log(5.0)))),
             eps=float(rng.uniform(1.5, 15)),
             r0=float(rng.choice([0.0, rng.uniform(0, 5)])),
             d=float(rng.uniform(0.3, 10.0)))
    Eb = ExcitonicInsulator2D(m_r=p["m_r"], eps=p["eps"], r0=p["r0"], d=p["d"], nk=100, ntheta=32).exciton_binding_energy()
    p["Eg"] = float(rng.uniform(-3, 1.2) * Eb)
    p["Tc"] = label(p)
    return p


def feats(X):
    m_r, q, eps, r0, d, Eg = X.T
    Ry = 13605.7 * m_r / eps ** 2; aB = eps * 0.0529177 / m_r
    return np.column_stack([X, np.log(m_r), np.log(q), Ry, Eg / Ry, d / aB, r0 / aB])


class Ensemble:
    def __init__(self, n_models=10):
        self.models = [HistGradientBoostingRegressor(max_iter=500, learning_rate=0.06, random_state=i)
                       for i in range(n_models)]

    def fit(self, X, y):
        rng = np.random.default_rng(0)
        self.lo, self.hi = X.min(0), X.max(0)
        for m in self.models:
            idx = rng.integers(0, len(X), len(X))
            m.fit(feats(X[idx]), np.log1p(y[idx]))
        return self

    def predict(self, X):
        P = np.expm1(np.array([m.predict(feats(X)) for m in self.models]))
        outside = np.any((X < self.lo) | (X > self.hi), axis=1)
        return P.mean(0), P.std(0), outside


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    with Pool() as pool:
        rows = pool.map(sample, range(10000, 10000 + n), chunksize=8)
    X = np.array([[r[f] for f in FEATS] for r in rows]); y = np.array([r["Tc"] for r in rows])
    np.savetxt("ei_dataset_v2.csv", np.column_stack([X, y]), delimiter=",", header=",".join(FEATS + ["Tc"]), comments="")
    tr, te = train_test_split(np.arange(n), test_size=0.2, random_state=0)
    ens = Ensemble().fit(X[tr], y[tr])
    mu, sd, _ = ens.predict(X[te])
    cover = np.mean(np.abs(mu - y[te]) <= 2 * sd + 1.0)
    print(f"{n} samples; T_c > 0 in {np.mean(y > 0):.2f}; max label {y.max():.0f} K")
    print(f"test: R2 = {r2_score(y[te], mu):.3f}, MAE = {mean_absolute_error(y[te], mu):.1f} K, "
          f"fraction within 2 sigma (+1 K) = {cover:.2f}")
    # inverse design with uncertainty: hBN-like (eps >= 4.5) vs anything (eps >= 1.5)
    rng = np.random.default_rng(3)
    N = 200_000
    C = np.column_stack([np.exp(rng.uniform(np.log(0.01), np.log(1.0), N)), np.exp(rng.uniform(np.log(0.2), np.log(5), N)),
                         rng.uniform(1.5, 15, N), np.zeros(N), rng.uniform(0.6, 10, N), rng.uniform(-300, 300, N)])
    m_, s_, out = ens.predict(C)
    for name, mask in (("eps >= 4.5 (hBN)", C[:, 2] >= 4.5), ("any eps >= 1.5", C[:, 2] >= 1.5)):
        ucb = np.where(mask & ~out, m_ + 2 * s_, -1)
        i = np.argsort(ucb)[::-1][:3]
        print(f"\ninverse design, {name}: top by mean + 2 sigma, verified with the physics label")
        for j in i:
            p = dict(zip(FEATS, C[j]))
            print(f"  m_r={p['m_r']:.2f} m_h/m_e={p['mass_ratio']:.2f} eps={p['eps']:.1f} d={p['d']:.2f} nm Eg={p['Eg']:6.0f} meV"
                  f" | surrogate {m_[j]:4.0f} +- {s_[j]:3.0f} K | physics {label(p):4.0f} K")
    print(f"\nlabels >= 300 K: {(y >= 300).sum()} of {n};  >= 350 K: {(y >= 350).sum()}")
