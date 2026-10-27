"""
Screening-limited ceiling on the 2D exciton-superfluid T_c.

Inputs (benchmarked in screening.py against Neilson, Perali & Hamilton 2014):
SS-screened superfluid exists only for r_s > r_s,c(d/a_B). With n_s <= n,
the Nelson-Kosterlitz bound gives T_c <= Ry*/(4 r_s,c^2).
Single-component layers (no spin/valley degeneracy) as in the benchmark;
extra degeneracy adds screening channels and lowers the ceiling further.

Usage:  python tc_ceiling.py
"""
import numpy as np
from ei_model import KB

# r_s,c(d/a_B) from `python screening.py` (SS onset, condensate fraction > 0.1)
X = np.array([0.3, 0.7, 1.0, 1.5, 2.0, 3.0])
RSC = np.array([1.29, 1.97, 2.38, 2.98, 3.53, 4.46])


def rs_c(x):
    slope = (RSC[-1] - RSC[-2]) / (X[-1] - X[-2])
    return np.where(x <= X[-1], np.interp(x, X, RSC), RSC[-1] + slope * (x - X[-1]))


def ceiling(m_r, eps, d_nm):
    Ry = 13605.7 * m_r / eps ** 2
    aB = eps * 0.0529177 / m_r
    x = d_nm / aB
    return Ry / (4 * rs_c(np.maximum(x, X[0])) ** 2) / KB, x


if __name__ == "__main__":
    m, e, d = np.meshgrid(np.geomspace(0.05, 0.6, 60), np.linspace(2, 15, 60), np.linspace(0.6, 3.0, 25), indexing="ij")
    T, x = ceiling(m, e, d)
    print("examples (T_c ceiling):")
    for args in ((0.22, 4.5, 1.6), (0.22, 4.5, 0.7), (0.30, 3.0, 0.7), (0.5, 2.0, 0.6)):
        t, xx = ceiling(*args)
        print(f"  m_r={args[0]}, eps={args[1]}, d={args[2]} nm  (d/a_B={xx:.2f}): T_c <= {t:.0f} K")
    i = np.unravel_index(np.argmax(T), T.shape)
    print(f"\nbest in box: m_r={m[i]:.2f}, eps={e[i]:.1f}, d={d[i]:.1f} nm -> T_c <= {T[i]:.0f} K (d/a_B={x[i]:.2f})")
    for Top in (300, 350):
        ok = T >= Top
        print(f"grid points with ceiling >= {Top} K: {ok.sum()} of {ok.size};"
              f" with eps >= 4.5 (hBN): {(ok & (e >= 4.5)).sum()}")
        if ok.any():
            print(f"   needs eps <= {e[ok].max():.1f}, m_r >= {m[ok].min():.2f}")
    if (x[T >= 300] < X[0]).any():
        print("   note: some passing points have d/a_B < 0.3 (outside benchmarked range)")
