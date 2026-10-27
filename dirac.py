"""
Graphene double layer (massless Dirac bands) exciton pairing, with three
treatments of screening.

Model (two-band, as in Min, Bistritzer, Su, MacDonald, PRB 78, 121401(R) (2008) [M08]):
  n-layer conduction band and p-layer valence band, xi_k = hbar v k - E_F,
  interlayer interaction with chirality form factor [M08 Eq. (3)]:
      V_kk' = V(|k - k'|) (1 + cos theta)/2,   (1+cos)/2 = ((k+k')^2 - q^2)/(4 k k').
  N_s N_v = 4 independent spin-valley copies.
Screening
  US  unscreened  V = 2 pi e*^2 e^{-qd}/q                       [M08]
  NS  normal-state RPA with graphene's static Pi = N_s N_v nu (q <= 2 p_F),
      Kharitonov & Efetov, PRB 78, 241401(R) (2008) [KE08] Eq. (5):
      V = 2 pi e*^2 e^{-qd} / (q + 2 kappa + kappa^2 (1 - e^{-qd})/q),
      kappa = 2 pi N_s N_v e*^2 nu,  nu = E_F / (2 pi (hbar v)^2)
  SS  screened by the paired state (Neilson, Perali, Hamilton 2014 approach):
      intraband Pi from BCS coherence factors with the chirality factor,
      interband part kept at its normal-state value N_s N_v nu - Pi_intra,normal(q).
Benchmarks (benchmark_dirac in __main__):
  KE08 Eq. (12): lambda_c = nu int dtheta/2pi V(2 p_F sin theta/2) (1+cos theta)/2
                 and its maximum 1/16 [Eq. (2)];  V_us(p_F)/V(p_F) ~ 9 at r_s = 1, d = 0 [p.3].
Units: meV, nm, K.  hbar v = 658.2 meV nm (v = 1e6 m/s).
"""
import numpy as np
from scipy.special import ellipk, ellipe
from ei_model import E2, KB

HBAR_V = 658.212


class DiracBilayer:
    def __init__(self, eps=3.9, d=1.0, E_F=258.0, mode="US", Ns=2, Nv=2, cut=8.0, nk=220, ntheta=96, cluster=4.0):
        self.eps, self.d, self.EF, self.mode = eps, d, E_F, mode
        self.N4 = Ns * Nv
        self.kF = E_F / HBAR_V
        self.cluster = cluster
        self.e2 = E2 / eps                       # e*^2  [meV nm]
        self.alpha = self.e2 / HBAR_V            # r_s = e^2/(eps hbar v)
        self.nu = E_F / (2 * np.pi * HBAR_V ** 2)     # DOS per spin-valley [meV^-1 nm^-2]
        self.kappa = 2 * np.pi * self.N4 * self.e2 * self.nu
        # radial grid: dense near k_F, cut-off at cut * k_F
        x = (np.arange(nk) + 0.5) / nk
        self.x_edges = np.arange(nk + 1) / nk
        self._Lam = cut * self.kF
        self.k = self._kmap(x, self._Lam)
        self.w = np.diff(self._kmap(self.x_edges, self._Lam))
        s, ws = np.polynomial.legendre.leggauss(ntheta)
        s = 0.5 * (s + 1); ws = 0.5 * ws
        self.th = np.pi * s ** 2
        self.wth = 2 * np.pi * s * ws * 2
        self.Pi_q = None
        self.build_kernel()

    def _kmap(self, x, Lam):
        # piecewise: [0, 2kF] clustered at kF (tanh), then [2kF, Lam] geometric-like
        x = np.asarray(x, float)
        a = 0.6
        c = self.cluster
        k1 = self.kF * (1 + np.sinh(c * (2 * x / a - 1)) / np.sinh(c))     # dense at k_F
        t = np.clip((x - a) / (1 - a), 0, 1)
        k2 = 2 * self.kF * (Lam / (2 * self.kF)) ** t
        return np.where(x <= a, k1, k2)

    # ---- interactions -------------------------------------------------------
    def V(self, q):
        q = np.maximum(q, 1e-12)
        ed = np.exp(-q * self.d)
        if self.mode == "US":
            return 2 * np.pi * self.e2 * ed / q
        if self.mode == "NS":
            ka = self.kappa
            return 2 * np.pi * self.e2 * ed / (q + 2 * ka + ka ** 2 * (-np.expm1(-q * self.d)) / q)
        # SS
        v = 2 * np.pi * self.e2 / q
        P = np.interp(q, self.qgrid, self.Pi_q, right=self.Pi_q[-1]) if self.Pi_q is not None else self.N4 * self.nu
        return v * ed / (1 + 2 * v * P + (v * P) ** 2 * (-np.expm1(-2 * q * self.d)))

    def _coul_chiral(self, k, kp):
        """int_0^{2pi} dtheta (2 pi e*^2/q) (1+cos)/2, analytic via elliptic K, E."""
        s = k + kp
        m = np.minimum(4 * k * kp / s ** 2, 1 - 1e-16)
        I1 = 4 * ellipk(m) / s          # int dtheta / q
        I2 = 4 * s * ellipe(m)          # int dtheta q
        return 2 * np.pi * self.e2 * (s ** 2 * I1 - I2) / (4 * k * kp)

    def build_kernel(self):
        k, w = self.k, self.w
        K, Kp = np.meshgrid(k, k, indexing="ij")
        q = np.sqrt(K[..., None] ** 2 + Kp[..., None] ** 2 - 2 * K[..., None] * Kp[..., None] * np.cos(self.th))
        F = 0.5 * (1 + np.cos(self.th))
        if self.mode == "NS":
            W = np.sum(self.V(q) * F * self.wth, axis=-1)            # finite at q -> 0
        else:
            W = self._coul_chiral(K, Kp)
            xs, ws = np.polynomial.legendre.leggauss(16)
            for i in range(len(k)):                                   # cell-average log singularity
                a, b = self.x_edges[i], self.x_edges[i + 1]
                xx = 0.5 * (b - a) * xs + 0.5 * (a + b)
                kk = self._kmap(xx, self._Lam)
                jac = np.gradient(kk, xx)
                W[i, i] = 0.5 * (b - a) * np.sum(ws * self._coul_chiral(k[i], kk) * kk * jac) / (k[i] * w[i])
            R = self.V(q) - 2 * np.pi * self.e2 / np.maximum(q, 1e-12)
            W = W + np.sum(R * F * self.wth, axis=-1)
        self.M = W * (k * w)[None, :] / (2 * np.pi) ** 2

    # ---- gap equation -------------------------------------------------------
    def xi(self):
        return HBAR_V * self.k - self.EF

    def lin_eig(self, T):
        xi = self.xi(); b = 1 / (KB * T)
        f = np.where(np.abs(xi) * b < 1e-8, b / 4, np.tanh(0.5 * b * xi) / (2 * xi))
        return np.max(np.real(np.linalg.eigvals(self.M * f[None, :])))

    def solve_gap(self, T=0.0, D0=None, maxiter=8000, tol=1e-8, mix=0.5):
        xi = self.xi()
        D = np.full_like(xi, 0.1 * self.EF) if D0 is None else D0.copy()
        for _ in range(maxiter):
            E = np.sqrt(xi ** 2 + D ** 2)
            occ = np.tanh(E / (2 * KB * T)) if T > 0 else 1.0
            Dn = self.M @ (D / (2 * E) * occ)
            if np.max(np.abs(Dn - D)) < tol * max(1, D.max()):
                D = Dn; break
            D = mix * Dn + (1 - mix) * D
        D = np.abs(D)
        return np.where(D.max() < 1e-6 * self.EF, 0.0, D)

    # ---- KE08 analytic coupling ---------------------------------------------
    def lambda_c(self):
        th = np.linspace(-np.pi, np.pi, 4001)[1:-1]
        q = 2 * self.kF * np.abs(np.sin(th / 2))
        integ = self.V(q) * 0.5 * (1 + np.cos(th))
        trap = getattr(np, "trapezoid", None) or np.trapz      # numpy < 2.0 compatibility
        return self.nu * trap(integ, th) / (2 * np.pi)


def _intraband_pi(b, D, qgrid, nphi=128):
    """N_s N_v sum_k (1+cos)/2 [u^2 v'^2 + v^2 u'^2 - 2 u v u' v'] / (E + E')  (T = 0)."""
    xi = b.xi()
    E = np.sqrt(xi ** 2 + D ** 2)
    u = np.sqrt(0.5 * (1 + xi / E)); v = np.sqrt(0.5 * (1 - xi / E))
    phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    k = b.k[:, None]
    out = []
    for q in qgrid:
        kp = np.sqrt(k ** 2 + q ** 2 + 2 * k * q * np.cos(phi))
        cos_kk = (k + q * np.cos(phi)) / np.maximum(kp, 1e-12)       # angle between k and k+q
        F = 0.5 * (1 + cos_kk)
        Dp = np.interp(kp, b.k, D, right=0.0)
        xip = HBAR_V * kp - b.EF
        Ep = np.sqrt(xip ** 2 + Dp ** 2)
        up = np.sqrt(0.5 * (1 + xip / Ep)); vp = np.sqrt(0.5 * (1 - xip / Ep))
        num = u[:, None] ** 2 * vp ** 2 + v[:, None] ** 2 * up ** 2 - 2 * u[:, None] * v[:, None] * up * vp
        integ = np.sum(F * num / (E[:, None] + Ep), axis=1) * (2 * np.pi / nphi)
        out.append(b.N4 * np.sum(b.k * b.w * integ) / (2 * np.pi) ** 2)
    return np.array(out)


def solve_SS(b, n_iter=25, mix=0.5, verbose=False):
    """Start from the US paired state, screen with Pi_SS = N nu + [Pi_intra,SS - Pi_intra,normal]."""
    b.mode = "US"; b.Pi_q = None; b.build_kernel()
    D = b.solve_gap()
    b.qgrid = np.concatenate([[0.0], np.geomspace(1e-3, 6, 70) * b.kF])
    Pi_n = _intraband_pi(b, np.zeros_like(D), b.qgrid)
    Pi_old = None
    b.mode = "SS"
    for it in range(n_iter):
        Pi = b.N4 * b.nu + _intraband_pi(b, D, b.qgrid) - Pi_n
        Pi = np.maximum(Pi, 0.0)
        if Pi_old is not None:
            Pi = mix * Pi + (1 - mix) * Pi_old
        b.Pi_q = Pi; b.build_kernel()
        Dn = b.solve_gap(D0=np.maximum(D, 1e-3 * b.EF))
        if verbose:
            print(f"   SS iter {it}: Dmax = {Dn.max():.3f} meV")
        conv = np.max(np.abs(Dn - D)) < 1e-3 * max(D.max(), 1e-9)
        D, Pi_old = Dn, Pi
        if conv or D.max() == 0:
            break
    return D


if __name__ == "__main__":
    print("Dirac benchmarks vs Kharitonov & Efetov 2008 [KE08]")
    b = DiracBilayer(eps=2.19, d=0.0, mode="NS", nk=300, cluster=6.0, ntheta=64)
    bus = DiracBilayer(eps=2.19, d=0.0, mode="US", nk=40)
    print(f"  V_us(pF)/V_NS(pF) = {bus.V(b.kF) / b.V(b.kF):.2f}   [KE08 p.3: ~9]")
    print(f"  lambda_c (alpha = 1)  = {b.lambda_c():.4f};  alpha -> 100: "
          f"{DiracBilayer(eps=0.0219, d=0.0, mode='NS', nk=40).lambda_c():.4f}   [KE08 Eq. 2: max 1/16]")
    Ts = np.array([0.02, 0.01, 0.005, 0.0025, 0.00125]) * b.EF / KB
    sl = np.polyfit(np.log(1 / Ts), [b.lin_eig(T) for T in Ts], 1)[0]
    print(f"  numerical kernel: d eig / d ln(1/T) = {sl:.4f}   [= lambda_c, KE08 Eq. 10-12]")

    print("\nBiSFET operating point (Banerjee 2009: n = 4.9e12 cm^-2 -> E_F = 258 meV; ~1 nm SiO2-like barrier)")
    for eps in (3.9,):
        for mode in ("US", "NS"):
            bb = DiracBilayer(eps=eps, d=1.0, E_F=258.0, mode=mode, nk=220, cluster=6.0, ntheta=64)
            D = bb.solve_gap()
            print(f"  eps={eps} d=1 nm {mode}: Delta_max = {D.max():7.3f} meV")
        bb = DiracBilayer(eps=eps, d=1.0, E_F=258.0, nk=220, cluster=6.0, ntheta=64)
        D = solve_SS(bb, verbose=False)
        print(f"  eps={eps} d=1 nm SS: Delta_max = {D.max():7.3f} meV")
    print(f"  M08 bound T_KT <= E_F/8 = {258 / 8 / KB:.0f} K holds only if a condensate exists")
