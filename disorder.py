"""
Disorder pair-breaking (gap 4c), Abrikosov-Gor'kov form as applied to e-h bilayers by
R. Bistritzer & A. H. MacDonald, PRL 101, 256406 (2008), arXiv:0808.1310 [BM08], Eq. (9):
    ln(T_c0 / T_c) = psi(1/2 + alpha/pi) - psi(1/2),   alpha = delta / (2 T_c),
    delta = 1/(2 tau_S) + 1/(2 tau_D)   (pair-breaking rate, hbar = k_B = 1).
Benchmarks [BM08 p.2]: weak disorder T_c = T_c0 - pi delta / 4; condensate destroyed
at delta ~ 0.88 T_c0.
"""
import numpy as np
from scipy.special import digamma
from scipy.optimize import brentq

HBAR_MEV_PS = 0.6582
KB = 0.0861733


def tc_ratio(delta_over_tc0):
    """T_c / T_c0 from the AG equation; 0 if the condensate is destroyed."""
    x = delta_over_tc0
    if x == 0:
        return 1.0
    f = lambda t: np.log(1 / t) - (digamma(0.5 + x / (2 * t) / np.pi) - digamma(0.5))
    try:
        return brentq(f, 1e-8, 1.0)
    except ValueError:
        return 0.0


def critical_delta():
    lo, hi = 0.5, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if tc_ratio(mid) > 1e-6 else (lo, mid)
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    x = 1e-3
    print(f"weak-disorder slope d(T_c/T_c0)/d(delta/T_c0) = {(tc_ratio(x) - 1) / x:.4f}   [BM08: -pi/4 = {-np.pi / 4:.4f}]")
    print(f"critical delta = {critical_delta():.4f} T_c0   [BM08: ~0.88]")
    print("\nrequired momentum-scattering time to keep a target T_c0 (delta = hbar/2tau < 0.88 k_B T_c0):")
    for T in (49, 121, 300, 350):
        tau = HBAR_MEV_PS / (2 * 0.88 * KB * T)
        print(f"  T_c0 = {T:3d} K: tau > {tau * 1e3:.0f} fs  (illustrative m* = 0.44 m0: mobility > {tau * 1e-12 * 1.602e-19 / (0.44 * 9.109e-31) * 1e4:.0f} cm^2/Vs)")
