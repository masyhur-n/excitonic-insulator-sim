"""
Fixed-density electron-hole bilayer with self-consistent screening in the
superfluid state (SS).

Model and conventions follow
  [P07] P. Pieri, D. Neilson, G. C. Strinati, PRB 75, 113301 (2007), Eqs. (2)-(4):
        one species per layer (no spin factor), n = sum_k v_k^2,
        r_s = 1/sqrt(pi n) in units of a_B*, energies in Ry*;
        xi_k = k^2/(4 m_r) - mu  (equal-mass average dispersion; our Eg = -2 mu).
  [N14] D. Neilson, A. Perali, A. R. Hamilton, PRB 89, 060502(R) (2014), arXiv:1308.0280:
        Eq. (2) unscreened (US), Eq. (4) screened in the superfluid state (SS):
          V_SS(q) = V0(q) / [1 + 2 v_q Pi(q) + (v_q Pi(q))^2 (1 - e^{-2qd})],
          v_q = 2 pi e^2 / (eps q),  V0 = v_q e^{-qd},  Pi = Pi_n + Pi_a.
        Equal masses are assumed in Pi (m_e = m_h = 2 m_r).
  T = 0 superfluid-state polarizabilities (BCS coherence factors; they reduce to
  the static 2D Lindhard function when Delta -> 0, and Pi(q -> 0) = 0 in the
  superfluid, i.e. "superfluidity kills screening" [N14 p.3]):
          Pi_n(q) =  sum_k (u_k^2 v_k'^2 + v_k^2 u_k'^2) / (E_k + E_k'),
          Pi_a(q) = -sum_k 2 u_k v_k u_k' v_k' / (E_k + E_k'),   k' = k + q.
"""
import numpy as np
from scipy.optimize import brentq
from ei_model import ExcitonicInsulator2D, E2


class Bilayer:
    def __init__(self, d_aB=1.0, m_r=0.25, eps=5.0, nk=140, ntheta=96):
        self.base = ExcitonicInsulator2D(m_r=m_r, eps=eps, d=0.0, nk=nk)
        self.aB, self.Ry = self.base.a_B, self.base.Ry
        self.d = d_aB * self.aB
        self.k, self.w = self.base.k, self.base.w
        self.meas = self.k * self.w / (2 * np.pi)
        self.ntheta = ntheta
        self.M_coul = self.base.M.copy()          # bare 1/q part with analytic angle integral
        s, ws = np.polynomial.legendre.leggauss(ntheta)
        s = 0.5 * (s + 1); ws = 0.5 * ws
        self.th = np.pi * s ** 2                  # cluster near theta = 0
        self.wth = 2 * np.pi * s * ws * 2         # d theta, x2 for [pi, 2 pi]
        self.Pi = None
        self.set_interaction(None)

    # ---- interaction -------------------------------------------------------
    def vq(self, q):
        return 2 * np.pi * E2 / (self.base.eps * np.maximum(q, 1e-12))

    def V(self, q):
        v = self.vq(q)
        V0 = v * np.exp(-q * self.d)
        if self.Pi is None:
            return V0
        P = np.interp(q, self.qgrid, self.Pi, right=self.Pi[-1])
        D = 1 + 2 * v * P + (v * P) ** 2 * (-np.expm1(-2 * q * self.d))
        return V0 / D

    def set_interaction(self, Pi_grid):
        self.Pi = Pi_grid
        K, Kp = np.meshgrid(self.k, self.k, indexing="ij")
        q = np.sqrt(K[..., None] ** 2 + Kp[..., None] ** 2 - 2 * K[..., None] * Kp[..., None] * np.cos(self.th))
        R = self.V(q) - self.vq(q)                # regular remainder (singular 1/q kept analytic)
        W = np.sum(R * self.wth, axis=-1)
        self.M = self.M_coul + W * (self.k * self.w)[None, :] / (2 * np.pi) ** 2

    # ---- polarizability in the superfluid state ----------------------------
    def polarizability(self, mu, D):
        xi = self.xi(mu)
        E = np.sqrt(xi ** 2 + D ** 2)
        u = np.sqrt(0.5 * (1 + xi / E)); v = np.sqrt(0.5 * (1 - xi / E))
        self.qgrid = np.concatenate([[0.0], np.geomspace(1e-3, 40, 90) / self.aB])
        phi = np.linspace(0, 2 * np.pi, 129)[:-1]
        dphi = phi[1] - phi[0]
        k = self.k[:, None]
        out = []
        for q in self.qgrid:
            kp = np.sqrt(k ** 2 + q ** 2 + 2 * k * q * np.cos(phi))
            Dp = np.interp(kp, self.k, D, right=0.0)
            xip = HB2_4MR(self) * kp ** 2 - mu
            Ep = np.sqrt(xip ** 2 + Dp ** 2)
            up = np.sqrt(0.5 * (1 + xip / Ep)); vp = np.sqrt(0.5 * (1 - xip / Ep))
            num = (u[:, None] ** 2 * vp ** 2 + v[:, None] ** 2 * up ** 2) - 2 * u[:, None] * v[:, None] * up * vp
            integ = np.sum(num / (E[:, None] + Ep), axis=1) * dphi
            out.append(np.sum(self.k * self.w * integ) / (2 * np.pi) ** 2)
        return np.array(out)

    # ---- mean-field at fixed density ---------------------------------------
    def xi(self, mu):
        return HB2_4MR(self) * self.k ** 2 - mu

    def gap(self, mu, D0=None):
        self.base.M = self.M
        return self.base.solve_gap(-2 * mu, delta0=D0, tol=1e-8, maxiter=6000)

    def density(self, mu, D):
        xi = self.xi(mu)
        E = np.sqrt(xi ** 2 + D ** 2)
        return np.sum(self.meas * 0.5 * (1 - xi / np.where(E > 0, E, 1)))

    def solve_fixed_n(self, rs, D0=None):
        n = 1 / (np.pi * (rs * self.aB) ** 2)
        state = {"D": D0}

        def f(mu):
            D = self.gap(mu, state["D"] if state["D"] is not None and state["D"].max() > 0 else None)
            state["D"] = D
            return self.density(mu, D) - n

        lo, hi = -3 * self.Ry, 1.0
        while f(hi) < 0:
            hi *= 2
        mu = brentq(f, lo, hi, xtol=1e-6 * self.Ry)
        D = self.gap(mu, state["D"])
        return mu, D

    def solve_SS(self, rs, n_iter=30, mix=0.5, verbose=False):
        """Start from the unscreened (US) coherent state, then screen with the
        superfluid-state polarizability until self-consistent [N14]."""
        self.set_interaction(None)
        mu, D = self.solve_fixed_n(rs)
        Pi_old = None
        for it in range(n_iter):
            Pi = self.polarizability(mu, D)
            if Pi_old is not None:
                Pi = mix * Pi + (1 - mix) * Pi_old
            self.set_interaction(Pi)
            mu, Dn = self.solve_fixed_n(rs, D)
            ch = np.max(np.abs(Dn - D)) / max(D.max(), 1e-12)
            D, Pi_old = Dn, Pi
            if verbose:
                print(f"   iter {it}: Dmax = {D.max() / self.Ry:.4f} Ry*, mu = {mu / self.Ry:+.3f} Ry*")
            if ch < 1e-3 or D.max() < 1e-6 * self.Ry:
                break
        return mu, D


def HB2_4MR(b):
    from ei_model import HB2_2M0
    return HB2_2M0 / (2 * b.base.m_r)


def condensate_fraction(b, mu, D):
    """c = sum u^2 v^2 / sum v^2 (standard BCS definition)."""
    xi = b.xi(mu); E = np.sqrt(xi ** 2 + D ** 2)
    v2 = 0.5 * (1 - xi / E); u2 = 1 - v2
    return np.sum(b.meas * u2 * v2) / np.sum(b.meas * v2)


def onset_rs(d_aB, c_min=0.1, lo=1.0, hi=8.0, tol=0.05, **kw):
    """Smallest r_s with a large-gap SS superfluid (condensate fraction > c_min)."""
    b = Bilayer(d_aB=d_aB, **kw)
    def sf(rs):
        mu, D = b.solve_SS(rs)
        return condensate_fraction(b, mu, D) > c_min
    if sf(lo):
        return lo
    if not sf(hi):
        return np.inf
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if sf(mid) else (mid, hi)
    return 0.5 * (lo + hi)


def tc_ceiling(Ry_meV, rs_c):
    """T_BKT <= pi hbar^2 n / (2 M k_B) with n = 1/(pi r_s^2 a_B^2), M = 4 m_r
       = Ry* / (4 r_s^2)   [Nelson & Kosterlitz 1977; n_s <= n]. Returns kelvin."""
    from ei_model import KB
    return Ry_meV / (4 * rs_c ** 2) / KB


if __name__ == "__main__":
    print("onset of the SS superfluid vs layer separation (benchmark: N14 Fig. 1a, DQMC r_s ~ 2.5 at d = 1)")
    for d in (0.3, 0.7, 1.0, 1.5, 2.0, 3.0):
        print(f"  d = {d:.1f} a_B : r_s,c = {onset_rs(d):.2f}", flush=True)
