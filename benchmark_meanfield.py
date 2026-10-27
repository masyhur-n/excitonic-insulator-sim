"""
Benchmarks of the crossover mean-field theory (MFT) against published values.

Reference [RT14]: M. Randeria and E. Taylor, "BCS-BEC Crossover and the Unitary
Fermi Gas", Annu. Rev. Condens. Matter Phys. 5, 209 (2014), arXiv:1306.5785.

  (a) 3D contact interaction, Eqs. (4)-(5) of [RT14], at unitarity 1/k_F a = 0:
      MFT gives Delta/E_F = 0.68 and xi_s = mu/E_F = 0.59 (Table I, p.13).
  (b) BCS-limit MFT T_c = (8 gamma / pi e^2) E_F exp(-pi / 2 k_F |a|) ([RT14] p.9).
  (c) In the BEC limit mu -> -E_b/2 ([RT14] p.9). In the excitonic model this is
      built in: mu_eff = -Eg/2 equals -E_b/2 exactly at the insulator edge Eg = E_b.

Units in (a, b): hbar = 1, 2m = 1 (so eps_k = k^2), k_F = 1, E_F = 1.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve, brentq
from ei_model import ExcitonicInsulator2D


def leggett_3d(inv_kfa):
    # gap eq. (renormalised with Eq. 5):  -m/(4 pi a) = int d^3k/(2pi)^3 [1/2E - 1/2eps]
    # number eq.:  n = k_F^3/(3 pi^2) = int d^3k/(2pi)^3 (1 - xi/E)
    m = 0.5

    def eqs(v):
        mu, D = v
        E = lambda k: np.sqrt((k * k - mu) ** 2 + D * D)
        gap = sum(quad(lambda k: k * k * (1 / (2 * E(k)) - 1 / (2 * k * k)), a, b, limit=400)[0]
                  for a, b in ((0, 4), (4, np.inf))) / (2 * np.pi ** 2)
        num = sum(quad(lambda k: k * k * D * D / (E(k) * (E(k) + k * k - mu)), a, b, limit=400)[0]
                  for a, b in ((0, 4), (4, np.inf))) / (2 * np.pi ** 2)
        return [gap + m * inv_kfa / (4 * np.pi), num - 1 / (3 * np.pi ** 2)]

    mu, D = fsolve(eqs, [0.6, 0.7])
    return mu, abs(D)


def bcs_tc_3d(inv_kfa):
    # linearised gap equation at mu = E_F (valid for 1/k_F a << -1)
    m = 0.5

    def f(T):
        h = lambda k: k * k * (np.tanh((k * k - 1) / (2 * T)) / (2 * (k * k - 1)) - 1 / (2 * k * k))
        g = sum(quad(h, a, b, limit=800)[0] for a, b in ((0, 1), (1, 3), (3, np.inf))) / (2 * np.pi ** 2)
        return g + m * inv_kfa / (4 * np.pi)

    return brentq(f, 1e-5, 1.0)


if __name__ == "__main__":
    mu, D = leggett_3d(0.0)
    print("(a) 3D MFT at unitarity   [RT14 Table I: Delta/E_F = 0.68, xi_s = 0.59]")
    print(f"    this code:            Delta/E_F = {D:.3f}, mu/E_F = {mu:.3f}")
    assert abs(D - 0.68) < 0.01 and abs(mu - 0.59) < 0.01

    inv = -2.0
    tc = bcs_tc_3d(inv)
    tc_formula = 8 * np.exp(np.euler_gamma) / (np.pi * np.e ** 2) * np.exp(np.pi / 2 * inv)  # gamma = e^C = 1.781
    print(f"(b) BCS-limit T_c at 1/k_F a = {inv}   [RT14 p.9: (8 e^gamma / pi e^2) E_F exp(-pi/2k_F|a|)]")
    print(f"    this code: T_c/E_F = {tc:.4f}   formula: {tc_formula:.4f}")
    assert abs(tc / tc_formula - 1) < 0.05

    m = ExcitonicInsulator2D(m_r=0.2, eps=5.0)
    Eb = m.exciton_binding_energy()
    print(f"(c) EI model: at the insulator edge Eg = E_b = {Eb:.1f} meV, mu_eff = -Eg/2 = {-Eb/2:.1f} meV = -E_b/2  [RT14 p.9]")
    print("all benchmarks passed")
