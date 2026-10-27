"""Generate a labelled dataset (material parameters -> EI properties) for ML.

Usage:  python make_dataset.py [n_samples] [out.csv]
"""
import os, sys, csv
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")     # one BLAS thread per worker: avoids fork deadlocks
import numpy as np
from multiprocessing import Pool
from ei_model import ExcitonicInsulator2D

FEATURES = ["m_r", "eps", "r0", "d", "q_s", "Eg"]
TARGETS = ["E_b", "Delta0", "qp_gap", "Tc", "size_ratio", "density", "regime"]


def sample(seed):
    rng = np.random.default_rng(seed)
    p = dict(m_r=float(np.exp(rng.uniform(np.log(0.02), np.log(0.6)))),
             eps=float(rng.uniform(2, 15)),
             r0=float(rng.choice([0.0, rng.uniform(0, 5)])),
             d=float(rng.choice([0.0, rng.uniform(0, 3)])),
             q_s=0.0)
    model = ExcitonicInsulator2D(**p, nk=100, ntheta=32)
    Eb = model.exciton_binding_energy()
    Eg = float(rng.uniform(-3, 1.3) * Eb)
    if Eg < 0 and rng.random() < 0.5:        # semimetal: free-carrier screening
        p["q_s"] = float(rng.uniform(0, 1) / model.a_B)
        model = ExcitonicInsulator2D(**p, nk=100, ntheta=32)
    D = model.solve_gap(Eg)
    obs = model.observables(Eg, D)
    Tc = model.critical_temperature(Eg) if D.max() > 0 else 0.0
    return {**p, "Eg": Eg, "E_b": Eb, "Delta0": obs["Delta0"], "qp_gap": obs["qp_gap"],
            "Tc": Tc if np.isfinite(Tc) else 2000.0, "size_ratio": obs["size_ratio"],
            "density": obs["density"], "regime": model.regime(Eg, D)}


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    out = sys.argv[2] if len(sys.argv) > 2 else "ei_dataset.csv"
    with Pool() as pool:
        rows = pool.map(sample, range(n), chunksize=8)
    with open(out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=FEATURES + TARGETS)
        wr.writeheader(); wr.writerows(rows)
    regs = [r["regime"] for r in rows]
    print(f"wrote {n} samples to {out}:", {k: regs.count(k) for k in set(regs)})
