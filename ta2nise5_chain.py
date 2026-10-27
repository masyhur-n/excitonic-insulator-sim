"""
Quasi-1D three-chain model of Ta2NiSe5 with exciton-phonon coupling.

Hamiltonian and published parameters: T. Kaneko, T. Toriyama, T. Konishi,
Y. Ohta, PRB 87, 035121 (2013), arXiv:1210.2787 [K13], Eq. (1):
  eps_c(k) = 2 t_c (cos k - 1) + D/2,  eps_f(k) = 2 t_f (cos k - 1) - D/2,
  t_c = -0.8, t_f = 0.4, D = 0.2 eV (fitted to DFT); two Ta c-chains, one Ni-Se f-chain;
  on-site U (U_c = U_f), interchain V (n^c_i + n^c_{i+1}) n^f_i with V = U/4;
  phonon (uniform shear, q = 0) with lambda = g^2/omega; filling 2 e per unit cell.
K13's mean-field equations are in their Supplemental Material, which was not
accessed here; the equations below are derived here from their Eq. (1) (Hartree for U, V;
Fock for V on both c-f bonds; static q = 0 phonon). Chain order parameters are
opposite (delta_1 = -delta_2, monoclinic shear), so f hybridises only with
c_- = (c_1 - c_2)/sqrt2 and c_+ stays the "nonbonding" band (as stated in K13).
Benchmarks [K13 Fig. 3(e,f), Fig. 4]:
  lambda = 0: EI onset at V_c = 0.6 eV;  V = 0: distortion onset at lambda_c = 0.125 eV;
  V = 0.55 eV, lambda = 0.04 eV: k_B T_c = 0.0783 eV.
Units: eV.
"""
import numpy as np
from scipy.optimize import brentq

TC, TF, D = -0.8, 0.4, 0.2


class ThreeChain:
    def __init__(self, V, lam, U=None, nk=512, fock='both'):
        self.V, self.lam, self.fock = V, lam, fock
        self.U = 4 * V if U is None else U
        self.k = -np.pi + 2 * np.pi * (np.arange(nk) + 0.5) / nk
        self.ec = 2 * TC * (np.cos(self.k) - 1) + D / 2
        self.ef = 2 * TF * (np.cos(self.k) - 1) - D / 2

    def _bands(self, phi0, phi1, delta, n, nf):
        """per-spin bands: c_+ (free) and the c_- / f 2x2 block."""
        Sc = self.U * n + 4 * self.V * nf
        Sf = self.U * nf + 8 * self.V * n
        ec, ef = self.ec + Sc, self.ef + Sf
        # Fock (both bonds) + phonon; sqrt2 from c_- projection
        if self.fock == "both":        # both c-f bonds, own amplitudes
            hf = -self.V * (phi0 + phi1 * np.exp(-1j * self.k))
        elif self.fock == "onsite":    # only the on-site bond
            hf = -self.V * phi0 * np.ones_like(self.k)
        else:                          # "onsite2": both bonds with the on-site amplitude (K13 text)
            hf = -2 * self.V * phi0 * np.ones_like(self.k)
        h = np.sqrt(2) * (hf + delta)
        avg, dif = 0.5 * (ec + ef), 0.5 * (ec - ef)
        R = np.sqrt(dif ** 2 + np.abs(h) ** 2)
        return ec, avg, dif, R, h

    def _fermi(self, E, mu, T):
        x = (E - mu) / max(T, 1e-9)
        return 0.5 * (1 - np.tanh(0.5 * x))

    def iterate(self, T, state, mix=0.5):
        phi0, phi1, delta, n, nf = state
        ec, avg, dif, R, h = self._bands(phi0, phi1, delta, n, nf)
        Ep, Em = avg + R, avg - R
        nk = len(self.k)

        def N(mu):    # electrons per spin per unit cell (target 1)
            return (np.sum(self._fermi(ec, mu, T)) + np.sum(self._fermi(Ep, mu, T)
                                                            + self._fermi(Em, mu, T))) / nk - 1.0
        mu = brentq(N, -6, 6, xtol=1e-12)
        fp, fm = self._fermi(Ep, mu, T), self._fermi(Em, mu, T)
        Rs = np.where(R > 1e-14, R, 1e-14)
        # <c_-^dag f> = -h* (f_+ - f_-)/(2R) ... : off-diagonal of 2x2 projector
        cmf = (np.conj(h) / (2 * Rs)) * (fp - fm)         # <f^dag c_-> per k
        ncm = 0.5 * (fp + fm) + dif / (2 * Rs) * (fp - fm)
        nff = 0.5 * (fp + fm) - dif / (2 * Rs) * (fp - fm)
        ncp = self._fermi(ec, mu, T)
        # chain-1 quantities: c_1 = (c_+ + c_-)/sqrt2
        new_phi0 = np.mean(cmf) / np.sqrt(2)
        new_phi1 = np.mean(cmf * np.exp(1j * self.k)) / np.sqrt(2)
        n_new = 0.5 * (np.mean(ncp) + np.mean(ncm))
        nf_new = np.mean(nff)
        # phonon: delta_1 = -2 lambda (1/L) sum_{k sigma} <c_1^dag f + h.c.> = -8 lambda Re phi0  (2 spins)
        new_delta = -8 * self.lam * np.real(new_phi0)
        new = (new_phi0, new_phi1, new_delta, n_new, nf_new)
        return tuple(mix * a + (1 - mix) * b for a, b in zip(new, state)), mu

    def solve(self, T=1e-4, state=None, tol=1e-10, maxiter=20000):
        st = state or (0.05, 0.05, -0.02, 0.05, 0.9)
        for _ in range(maxiter):
            new, mu = self.iterate(T, st)
            if max(abs(a - b) for a, b in zip(new, st)) < tol:
                st = new; break
            st = new
        phi0, phi1, delta, n, nf = st
        return dict(Delta=abs(self.V * phi0), phi1=abs(phi1), delta=abs(delta), n=n, nf=nf, mu=mu, state=st)

    def gap(self, res):
        ec, avg, dif, R, h = self._bands(*res["state"])
        Ep, Em = avg + R, avg - R
        lower = np.max(Em); upper = min(np.min(Ep), np.min(ec))
        return upper - lower

    def Tc(self, Tmax=0.3):
        def ordered(T):
            return self.solve(T)["Delta"] + self.solve(T)["delta"] > 1e-5
        if not ordered(1e-3):
            return 0.0
        lo, hi = 1e-3, Tmax
        while hi - lo > 1e-4:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if ordered(mid) else (lo, mid)
        return 0.5 * (lo + hi)


def lattice_vs_exciton(V=0.55, lams=(0.0, 0.02, 0.04, 0.05)):
    """Share of the c-f hybridisation that is electronic (Fock) vs lattice (delta),
    and T_c, across K13's estimate lambda ~ 0.02-0.05 eV from the measured distortion."""
    rows = []
    for lam in lams:
        m = ThreeChain(V, lam)
        r = m.solve()
        h_exc = V * abs(r["state"][0] + r["state"][1])
        share = h_exc / (h_exc + r["delta"]) if h_exc + r["delta"] > 0 else np.nan
        rows.append((lam, h_exc, r["delta"], share, m.gap(r), m.Tc()))
    return rows


if __name__ == "__main__":
    print("K13 benchmarks")
    for V in (0.5, 0.58, 0.6, 0.62, 0.7):
        r = ThreeChain(V, 0.0).solve()
        print(f"  lambda=0, V={V}: Delta = {r['Delta']:.4f} eV   [K13: onset V_c = 0.6]")
    for lam in (0.10, 0.12, 0.125, 0.13, 0.15):
        r = ThreeChain(0.0, lam).solve()
        print(f"  V=0, lambda={lam}: delta = {r['delta']:.4f} eV   [K13: onset lambda_c = 0.125]")
    m = ThreeChain(0.55, 0.04)
    r = m.solve()
    print(f"  V=0.55, lambda=0.04, T=0: Delta = {r['Delta']:.4f}, delta = {r['delta']:.4f} eV, n = {r['n']:.3f}, gap = {m.gap(r):.3f} eV")
    print(f"  k_B T_c = {m.Tc():.4f} eV   [K13 Fig. 4: 0.0783 eV]")

    print("\nexcitonic vs lattice (V = 0.55 eV, U = 4V)")
    for lam, he, de, sh, g, tc in lattice_vs_exciton():
        print(f"  lambda={lam:.2f}: exc. hyb = {he:.4f}, lattice = {de:.4f} eV, excitonic share = {sh:.2f},"
              f" gap = {g:.3f} eV, k_B T_c = {tc:.4f} eV ({tc / 8.617e-5:.0f} K; /2.2 MFT = {tc / 8.617e-5 / 2.2:.0f} K)")
    r = ThreeChain(0.0, 0.04).solve()
    print(f"  pure lattice (V = 0) at lambda = 0.04: delta = {r['delta']:.4f} eV -> needs lambda >= 0.125 (3x the K13 estimate)")
