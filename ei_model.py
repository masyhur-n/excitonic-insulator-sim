"""
Mean-field (BCS-like) model of a 2D excitonic insulator.

Two bands (conduction c, valence v), electron-hole symmetric pairing:
    eps_c(k) = +Eg/2 + hbar^2 k^2 / 2 m_e
    eps_v(k) = -Eg/2 - hbar^2 k^2 / 2 m_h
    xi(k)    = (eps_c - eps_v)/2 = Eg/2 + hbar^2 k^2 / (4 m_r)

Gap equation (Keldysh / Jerome-Rice-Kohn / Comte-Nozieres):
    Delta(k) = sum_k' V(|k-k'|) Delta(k') / (2 E(k')) * tanh(E(k') / 2 kT)
    E(k)     = sqrt(xi(k)^2 + Delta(k)^2)

Interaction: screened 2D Coulomb with Rytova-Keldysh length r0 and
optional layer separation d (spatially indirect excitons in bilayers):
    V(q) = 2 pi e^2 / (eps sqrt(q^2 + q_s^2) (1 + r0 q)) * exp(-q d)
q_s is an optional static (Thomas-Fermi-like) screening wavevector from free
carriers; without it the BCS-side gap is overestimated at high density.

BEC-BCS crossover:
    Eg < 0  (semimetal, band overlap)  -> pairs form at a Fermi surface: BCS side
    0 < Eg < Eb (semiconductor)        -> preformed excitons condense: BEC side
    Eg > Eb                             -> ordinary band insulator (no condensate)

Units: energy meV, length nm, temperature K.
"""
import numpy as np
from scipy.special import ellipk
from scipy.linalg import eigh, eigvals

HB2_2M0 = 38.0998   # hbar^2 / (2 m0)  [meV nm^2]
E2 = 1439.964       # e^2 / (4 pi eps0) [meV nm]
KB = 0.0861733      # k_B [meV / K]


class ExcitonicInsulator2D:
    def __init__(self, m_r=0.2, eps=5.0, r0=0.0, d=0.0, q_s=0.0, nk=160, ntheta=48, g=1.0):
        """
        m_r   : reduced e-h mass (units of m0)
        eps   : background dielectric constant
        r0    : Rytova-Keldysh screening length [nm] (0 = plain Coulomb)
        d     : electron-hole layer separation [nm] (0 = same layer)
        q_s   : static screening wavevector [nm^-1] (0 = unscreened)
        g     : degeneracy (spin x valley) used only for the pair density
        """
        self.m_r, self.eps, self.r0, self.d, self.q_s, self.g = m_r, eps, r0, d, q_s, g
        # natural scales of the 2D hydrogen problem
        self.a_B = eps * 0.0529177 / m_r                 # exciton Bohr radius [nm]
        self.Ry = 13605.7 * m_r / eps**2                 # effective Rydberg [meV]
        self._build_grid(nk)
        self._build_kernel(ntheta)

    # ---------- grid & interaction kernel -----------------------------------
    def _build_grid(self, nk):
        # midpoint grid in x in (0,1), mapped k = s * tan(pi x / 2): dense near
        # k ~ 1/a_B, long tail for the 1/k^3 decay of Delta(k)
        s = 1.0 / self.a_B
        x = (np.arange(nk) + 0.5) / nk
        self.x_edges = np.arange(nk + 1) / nk
        self.kmap = lambda x: s * np.tan(0.5 * np.pi * x * 0.985)
        self.dkdx = lambda x: s * 0.5 * np.pi * 0.985 / np.cos(0.5 * np.pi * x * 0.985) ** 2
        self.k = self.kmap(x)
        self.w = self.dkdx(x) / nk                       # dk weights
        self.nk = nk

    def _V_regular(self, q):
        """V(q) minus its 2 pi e^2 / (eps sqrt(q^2+q_s^2)) part -- smooth at q -> 0."""
        pref = 2 * np.pi * E2 / self.eps
        qq = np.sqrt(q ** 2 + self.q_s ** 2) + 1e-12
        return pref * (np.exp(-q * self.d) / (1 + self.r0 * q) - 1.0) / qq

    def _coulomb_ang(self, k, kp):
        """int_0^{2pi} dtheta 2 pi e^2/(eps sqrt(|k-k'|^2+q_s^2))
           = (2 pi e^2/eps) 4 K(m) / sqrt((k+k')^2 + q_s^2)."""
        s2 = (k + kp) ** 2 + self.q_s ** 2
        m = np.minimum(4 * k * kp / s2, 1 - 1e-16)
        return 2 * np.pi * E2 / self.eps * 4 * ellipk(m) / np.sqrt(s2)

    def _build_kernel(self, ntheta):
        k, w = self.k, self.w
        K, Kp = np.meshgrid(k, k, indexing="ij")
        # bare Coulomb part (analytic angular integral, log-singular on diagonal)
        W = self._coulomb_ang(K, Kp)
        # diagonal: replace point value by its cell average (log singularity is integrable)
        xs, ws = np.polynomial.legendre.leggauss(16)
        for i in range(self.nk):
            a, b = self.x_edges[i], self.x_edges[i + 1]
            xx = 0.5 * (b - a) * xs + 0.5 * (a + b)
            kk = self.kmap(xx)
            vals = self._coulomb_ang(k[i], kk) * kk * self.dkdx(xx)
            W[i, i] = (0.5 * (b - a) * np.sum(ws * vals)) / (k[i] * w[i])
        # regular remainder (screening + layer separation), numeric angular integral
        if self.r0 > 0 or self.d > 0:
            th, wt = np.polynomial.legendre.leggauss(ntheta)
            th = 0.5 * np.pi * (th + 1)       # [0, pi]
            wt = 0.5 * np.pi * wt * 2         # x2 for [pi, 2pi]
            q = np.sqrt(K[..., None] ** 2 + Kp[..., None] ** 2
                        - 2 * K[..., None] * Kp[..., None] * np.cos(th))
            W = W + np.sum(self._V_regular(q) * wt, axis=-1)
        # Nystrom matrix:  sum_k'  ->  int k' dk' dtheta / (2pi)^2
        self.M = W * (k * w)[None, :] / (2 * np.pi) ** 2

    # ---------- single-particle -------------------------------------------
    def xi(self, Eg):
        return 0.5 * Eg + HB2_2M0 * self.k ** 2 / (2 * self.m_r)

    def exciton_binding_energy(self):
        """Lowest eigenvalue of the Wannier equation, -E_b (no condensate)."""
        H = np.diag(HB2_2M0 * self.k ** 2 / self.m_r) - self.M
        # symmetrise with sqrt(k w) so eigh can be used
        s = np.sqrt(self.k * self.w)
        Hs = (s[:, None] * H) / s[None, :]
        Hs = 0.5 * (Hs + Hs.T)
        return -eigh(Hs, eigvals_only=True, subset_by_index=[0, 0])[0]

    # ---------- gap equation ------------------------------------------------
    def solve_gap(self, Eg, T=0.0, delta0=None, tol=1e-7, maxiter=4000, mix=0.5):
        xi = self.xi(Eg)
        D = np.full(self.nk, max(self.Ry, 1.0)) if delta0 is None else delta0.copy()
        D = D / (1 + (self.k * self.a_B) ** 2) ** 1.5 if delta0 is None else D
        for it in range(maxiter):
            E = np.sqrt(xi ** 2 + D ** 2)
            occ = np.tanh(E / (2 * KB * T)) if T > 0 else 1.0
            Dn = self.M @ (D / (2 * E) * occ)
            err = np.max(np.abs(Dn - D))
            D = mix * Dn + (1 - mix) * D
            if err < tol * max(1.0, np.max(np.abs(D))):
                break
        D = np.abs(D) if D[0] < 0 else D
        if np.max(np.abs(D)) < 1e-4:
            D[:] = 0.0
        return D

    def critical_temperature(self, Eg, Tmax=2000.0):
        """Mean-field T_c from the linearised gap equation (largest eigenvalue = 1).
        On the BEC side this is the pair-formation scale T*, an upper bound on the
        true (BKT) condensation temperature."""
        xi = self.xi(Eg)

        def lam(T):
            b = 1.0 / (KB * T)
            f = np.where(np.abs(xi) * b < 1e-8, b / 4, np.tanh(0.5 * b * xi) / (2 * xi))
            return np.max(np.real(eigvals(self.M * f[None, :])))

        lo, hi = 0.1, Tmax
        if lam(lo) < 1:
            return 0.0
        if lam(hi) > 1:
            return np.inf
        for _ in range(50):
            mid = np.sqrt(lo * hi)
            lo, hi = (mid, hi) if lam(mid) > 1 else (lo, mid)
            if hi / lo < 1.001:
                break
        return np.sqrt(lo * hi)

    # ---------- observables -------------------------------------------------
    def observables(self, Eg, D):
        k, w = self.k, self.w
        xi = self.xi(Eg)
        E = np.sqrt(xi ** 2 + D ** 2)
        v2 = 0.5 * (1 - xi / E)                    # e-h pair occupation
        phi = D / (2 * E)                          # pair wavefunction in k
        meas = k * w / (2 * np.pi)                 # d^2k/(2pi)^2 after angle int.
        n = self.g * np.sum(v2 * meas)             # pair density [nm^-2]
        dphi = np.gradient(phi, k)
        norm = np.sum(phi ** 2 * meas)
        xi_pair = np.sqrt(np.sum(dphi ** 2 * meas) / norm) if norm > 0 else np.nan
        r_s = 1 / np.sqrt(np.pi * n) if n > 0 else np.inf
        return dict(
            Delta0=float(np.max(D)),               # max |Delta(k)|
            qp_gap=float(2 * np.min(E)),           # optical/tunnelling gap (2 min E)
            density=float(n),
            xi_pair=float(xi_pair),                # pair (coherence) size [nm]
            r_s=float(r_s),                        # interparticle spacing [nm]
            size_ratio=float(xi_pair / r_s) if np.isfinite(r_s) else 0.0,
            mu_eff=float(-Eg / 2),                 # >0 BCS, <0 BEC
            T_BKT=float(self.bkt_upper_bound(n)),
        )

    def bkt_upper_bound(self, n):
        """Upper bound on the 2D condensation temperature from the Nelson-Kosterlitz
        universal jump, n_s(T_BKT) = 2 M k_B T / (pi hbar^2) [PRL 39, 1201 (1977)],
        taking the superfluid density n_s <= n (the T = 0 pair density).
        Assumes equal e and h masses, so the exciton mass is M = m_e + m_h = 4 m_r."""
        M = 4 * self.m_r
        return np.pi * HB2_2M0 * n / (M * KB)

    def regime(self, Eg, D):
        if np.max(D) == 0:
            return "normal"
        return "BCS" if Eg < 0 else "BEC"
