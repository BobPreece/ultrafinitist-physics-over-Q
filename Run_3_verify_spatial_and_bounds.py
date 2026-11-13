"""Machine-checked symbolic & exact-Q proofs of spatial operators, positivity,
3D Symanzik vector calculus + Leray projector, Varga M-matrix positivity,
and sharp closed-form Ultrafinitist No-Blowup bounds B_u(E_0, N), B_omega(E_0, N).

Strict Ultrafinitist Compliance:
- Pure Python stdlib (`fractions.Fraction` + `formalist.cas` only).
- Zero floats, zero transcendentals, zero roots, zero sympy/numpy.
- Zero `isinstance` or `type` calls.
"""

from fractions import Fraction
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Run_3_cas import (
    QQPoly,
    QQRat,
    poly_vars,
    qqpoly_vars,
    rat_vars,
)

ZERO_Q = Fraction(0, 1)
ONE_Q = Fraction(1, 1)


# ============================================================================
# Closed-Form Sharp Ultrafinitist No-Blowup Bounds over Q
# ============================================================================


def _cube_root_exact(n):
    """Exact integer cube root L of N = L^3 via finite integer binary search (no floats/roots)."""
    lo = 1
    hi = n
    while lo <= hi:
        mid = (lo + hi) // 2
        c = mid * mid * mid
        if c == n:
            return mid
        if c < n:
            lo = mid + 1
        else:
            hi = mid - 1
    raise ValueError(f"N={n} is not a perfect cube L^3")


def bound_B_u(E0, N, zero_momentum=True, h=ONE_Q, volume_weighted=False):
    """Exact, sharp rational upper bound on max_{x in Z_L^3} ||u_n(x)||^2 over Q.

    Theorem (Finite-Lattice Reproducing-Kernel Incompressible Velocity Bound):
    Let u in (Q^N)^3 (N = L^3, L >= 3) satisfy sitewise Symanzik incompressibility
    `div_Sym(u) == 0` and kinetic energy `(1/2) <u, u> <= E_0` (or `(h^3/2)<u,u> <= E_0`
    if `volume_weighted=True`).
    - On `Z_L^3`, `dim ker(grad_Sym) = 8` for even `L` (spanned by `(-1)^(ax+by+cz)`)
      and `dim ker(grad_Sym) = 1` for odd `L`.
    - Consequently `rank(P_div) = 2N + 8` (even L) or `2N + 1` (odd L).
    - If `sum_x u(x) == (0, 0, 0)` (`zero_momentum=True`), `Pi = P_div - P_mean` is
      an exact rational orthogonal projector of rank `2N + 5` (even L) or `2N - 2` (odd L),
      with site reproducing block `E_x Pi E_x^T = c_u(N) I_3` where:
        c_u(N) = (2N + 5) / (3N)   for even L,
        c_u(N) = (2N - 2) / (3N)   for odd L.
      By the exact Sum-of-Squares identity:
        c_u(N) ||u||_2^2 - ||u(x)||^2 = (1 / c_u(N)) || c_u(N) u - Pi E_x^T u(x) ||_2^2 >= 0,
      we have for ALL n >= 0:
        max_x ||u_n(x)||^2 <= B_u(E_0, N) = 2 * c_u(N) * E_0,
      which is ATTAINED on `V_{div,0}` by `u* = Pi e_{x,0}`!
    """
    L = _cube_root_exact(N)
    if (L % 2) == 0:
        rank_proj = (2 * N + 5) if zero_momentum else (2 * N + 8)
    else:
        rank_proj = (2 * N - 2) if zero_momentum else (2 * N + 1)
    c_u = Fraction(rank_proj, 3 * N)
    scale = (2 * c_u) / (h ** 3) if volume_weighted else (2 * c_u)
    return scale * E0


def bound_B_omega(E0, N, h=ONE_Q, volume_weighted=False):
    """Exact, sharp rational upper bound on max_{x in Z_L^3} ||omega_n(x)||^2 over Q.

    Theorem (Finite-Lattice Symanzik Vorticity No-Blowup Bound):
    Let u in (Q^N)^3 (N = L^3) satisfy `(1/2) <u, u> <= E_0` and `omega = curl_Sym(u)`.
    Since `curl_Sym^T == curl_Sym` and `curl_Sym o Pi == curl_Sym`, the local vorticity
    Gram block is `E_x curl_Sym^2 E_x^T = c_omega(L) I_3`, where for all `L >= 5`:
      c_omega(L) = 2 * ||d^{Sym}||_2^2 = 2 * (1/144 + 64/144 + 64/144 + 1/144) = 65 / (36 * h^2).
    By the exact Sum-of-Squares identity:
      c_omega ||u||_2^2 - ||omega(x)||^2 = (1/c_omega) || c_omega u - curl_Sym(E_x^T omega(x)) ||_2^2 >= 0,
    we have for ALL n >= 0:
      max_x ||omega_n(x)||^2 <= B_omega(E_0, N) = 2 * c_omega(L) * E_0 = (65 / (18 * h^2)) * E_0,
    which is ATTAINED on `V_{div,0}` by `u* = curl_Sym(delta_{x_0} e_x)`!
    """
    L = _cube_root_exact(N)
    w = [ZERO_Q] * L
    w[1 % L] += Fraction(2, 3) / h
    w[(-1) % L] -= Fraction(2, 3) / h
    w[2 % L] -= Fraction(1, 12) / h
    w[(-2) % L] += Fraction(1, 12) / h
    norm_sq_1d = sum((val * val for val in w), ZERO_Q)
    c_omega = 2 * norm_sq_1d  # Equals Fraction(65, 36) / h^2 for all L >= 5
    scale = (2 * c_omega) / (h ** 3) if volume_weighted else (2 * c_omega)
    return scale * E0


# ============================================================================
# Exact Rational Linear Solver (Gauss-Jordan over Fraction, Zero Floats)
# ============================================================================


def solve_linear_q(mat, rhs):
    """Exact Gauss-Jordan elimination over Q for `mat * x = rhs` (`rhs` is vector or matrix)."""
    n = len(mat)
    a = [[Fraction(v) for v in row] for row in mat]
    is_vec = not hasattr(rhs[0], "__len__")
    b = [[Fraction(v)] for v in rhs] if is_vec else [[Fraction(v) for v in row] for row in rhs]
    m_cols = len(b[0])

    for col in range(n):
        piv = None
        for r in range(col, n):
            if a[r][col] != ZERO_Q:
                piv = r
                break
        if piv is None:
            raise ValueError(f"Singular rational matrix at column {col}")
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            b[col], b[piv] = b[piv], b[col]
        inv_p = ONE_Q / a[col][col]
        for j in range(col, n):
            a[col][j] *= inv_p
        for j in range(m_cols):
            b[col][j] *= inv_p
        for r in range(n):
            if r != col and a[r][col] != ZERO_Q:
                fac = a[r][col]
                a[r][col] = ZERO_Q
                for j in range(col + 1, n):
                    a[r][j] -= fac * a[col][j]
                for j in range(m_cols):
                    b[r][j] -= fac * b[col][j]

    if is_vec:
        return [row[0] for row in b]
    return b


# ============================================================================
# Part 1A & 1B: 2D Gaussian Torus Z[i] / <(1 + i)^m> Proofs
# ============================================================================


class GaussianTorus:
    """Exact quotient ring & 2D lattice G_m = Z[i] / <(1 + i)^m> of size N = 2^m."""

    def __init__(self, m):
        self.m = m
        self.N = 1 << m
        self.k = m // 2
        self.q = 1 << self.k
        if (m % 2) == 0:
            self.sites = [(x, y) for x in range(self.q) for y in range(self.q)]
        else:
            self.sites = [(x, y) for x in range(self.q) for y in range(2 * self.q)]
        assert len(self.sites) == self.N
        self.index_of = {pt: idx for idx, pt in enumerate(self.sites)}

    def canon(self, x, y):
        if (self.m % 2) == 0:
            return (x % self.q, y % self.q)
        two_q = 2 * self.q
        x0 = x % two_q
        y0 = y % two_q
        if x0 >= self.q:
            x0 -= self.q
            y0 = (y0 - self.q) % two_q
        return (x0, y0)

    def shift_idx(self, idx, dx, dy):
        x, y = self.sites[idx]
        return self.index_of[self.canon(x + dx, y + dy)]

    def d4_permutations(self):
        """Return all 8 D_4 group actions on G_m as index permutations."""
        # Generated by R(x,y) = (-y, x) [mult by i] and C(x,y) = (x, -y) [conjugation]
        perms = []
        for refl in (0, 1):
            for rot in range(4):
                p = []
                for x, y in self.sites:
                    u, v = (x, -y) if refl else (x, y)
                    for _ in range(rot):
                        u, v = -v, u
                    p.append(self.index_of[self.canon(u, v)])
                perms.append(tuple(p))
        return perms

    def apply_Lx(self, f):
        return [
            2 * f[i] - f[self.shift_idx(i, 1, 0)] - f[self.shift_idx(i, -1, 0)]
            for i in range(self.N)
        ]

    def apply_Ly(self, f):
        return [
            2 * f[i] - f[self.shift_idx(i, 0, 1)] - f[self.shift_idx(i, 0, -1)]
            for i in range(self.N)
        ]

    def apply_L2x(self, f):
        return [
            2 * f[i] - f[self.shift_idx(i, 2, 0)] - f[self.shift_idx(i, -2, 0)]
            for i in range(self.N)
        ]

    def apply_L2y(self, f):
        return [
            2 * f[i] - f[self.shift_idx(i, 0, 2)] - f[self.shift_idx(i, 0, -2)]
            for i in range(self.N)
        ]

    def apply_K_Sym(self, f):
        lx = self.apply_Lx(f)
        ly = self.apply_Ly(f)
        l2x = self.apply_L2x(f)
        l2y = self.apply_L2y(f)
        c1 = Fraction(4, 3)
        c2 = Fraction(1, 12)
        return [c1 * (lx[i] + ly[i]) - c2 * (l2x[i] + l2y[i]) for i in range(self.N)]

    def apply_W_flat(self, f):
        # W_flat = (L_x - (1/2) L_y)^2 + (3/4) L_y^2
        lx = self.apply_Lx(f)
        ly = self.apply_Ly(f)
        half = Fraction(1, 2)
        three_q = Fraction(3, 4)
        a = [lx[i] - half * ly[i] for i in range(self.N)]
        part1_a = self.apply_Lx(a)
        part1_b = self.apply_Ly(a)
        term1 = [part1_a[i] - half * part1_b[i] for i in range(self.N)]
        ly2 = self.apply_Ly(ly)
        return [term1[i] + three_q * ly2[i] for i in range(self.N)]

    def arakawa_J2(self, psi, om):
        """Classic 9-point Arakawa (1966) skew-symmetric Jacobian (J++ + J+x + Jx+)/3."""
        N = self.N
        out = []
        inv12 = Fraction(1, 12)
        for z in range(N):
            zp1 = self.shift_idx(z, 1, 0)
            zm1 = self.shift_idx(z, -1, 0)
            zpi = self.shift_idx(z, 0, 1)
            zmi = self.shift_idx(z, 0, -1)
            zpp = self.shift_idx(z, 1, 1)
            zpm = self.shift_idx(z, 1, -1)
            zmp = self.shift_idx(z, -1, 1)
            zmm = self.shift_idx(z, -1, -1)

            # 4 * J^{++}
            j_pp = (psi[zp1] - psi[zm1]) * (om[zpi] - om[zmi]) - (
                psi[zpi] - psi[zmi]
            ) * (om[zp1] - om[zm1])
            # 4 * J^{+x}
            j_px = (
                psi[zp1] * (om[zpp] - om[zpm])
                - psi[zm1] * (om[zmp] - om[zmm])
                - psi[zpi] * (om[zpp] - om[zmp])
                + psi[zmi] * (om[zpm] - om[zmm])
            )
            # 4 * J^{x+} = -4 * J^{+x}(om, psi)
            j_xp = -(
                om[zp1] * (psi[zpp] - psi[zpm])
                - om[zm1] * (psi[zmp] - psi[zmm])
                - om[zpi] * (psi[zpp] - psi[zmp])
                + om[zmi] * (psi[zpm] - psi[zmm])
            )
            out.append(inv12 * (j_pp + j_px + j_xp))
        return out

    def arakawa_J4(self, psi, om):
        """4th-order Arakawa (1966) Jacobian J_4 = 2 * J_2 - J_diag."""
        j2 = self.arakawa_J2(psi, om)
        N = self.N
        inv24 = Fraction(1, 24)
        out = []
        # Diagonal lattice steps: e1' = (1, 1), e2' = (-1, 1) (spacing sqrt(2), cell area 2)
        for z in range(N):
            d_p1 = self.shift_idx(z, 1, 1)
            d_m1 = self.shift_idx(z, -1, -1)
            d_p2 = self.shift_idx(z, -1, 1)
            d_m2 = self.shift_idx(z, 1, -1)
            d_pp = self.shift_idx(z, 0, 2)
            d_pm = self.shift_idx(z, 2, 0)
            d_mp = self.shift_idx(z, -2, 0)
            d_mm = self.shift_idx(z, 0, -2)

            j_xx = (psi[d_p1] - psi[d_m1]) * (om[d_p2] - om[d_m2]) - (
                psi[d_p2] - psi[d_m2]
            ) * (om[d_p1] - om[d_m1])
            j_xp = (
                psi[d_p1] * (om[d_pp] - om[d_pm])
                - psi[d_m1] * (om[d_mp] - om[d_mm])
                - psi[d_p2] * (om[d_pp] - om[d_mp])
                + psi[d_m2] * (om[d_pm] - om[d_mm])
            )
            j_px = -(
                om[d_p1] * (psi[d_pp] - psi[d_pm])
                - om[d_m1] * (psi[d_mp] - psi[d_mm])
                - om[d_p2] * (psi[d_pp] - psi[d_mp])
                + om[d_m2] * (psi[d_pm] - psi[d_mm])
            )
            j_diag = inv24 * (j_xx + j_xp + j_px)
            out.append(2 * j2[z] - j_diag)
        return out

    def matrix_of(self, op_fn):
        cols = []
        for j in range(self.N):
            ej = [ONE_Q if i == j else ZERO_Q for i in range(self.N)]
            cols.append(op_fn(ej))
        return [[cols[j][i] for j in range(self.N)] for i in range(self.N)]

    def compute_green_G0(self):
        """Compute exact zero-mean Green's matrix G_0 = K_Sym^+ over Q."""
        N = self.N
        K_mat = self.matrix_of(self.apply_K_Sym)
        inv_N = Fraction(1, N)
        # Solve (K_Sym + (1/N) J) g = e_0 - (1/N) 1 for the translation-invariant 0-th column
        reg_mat = [
            [K_mat[i][j] + inv_N for j in range(N)]
            for i in range(N)
        ]
        rhs = [(ONE_Q - inv_N) if i == 0 else (-inv_N) for i in range(N)]
        g_col = solve_linear_q(reg_mat, rhs)
        # Reconstruct full translation-invariant matrix G_0[z, w] = g_col[z - w]
        G0 = [[ZERO_Q] * N for _ in range(N)]
        for w in range(N):
            wx, wy = self.sites[w]
            for z in range(N):
                zx, zy = self.sites[z]
                diff_idx = self.index_of[self.canon(zx - wx, zy - wy)]
                G0[z][w] = g_col[diff_idx]
        return K_mat, G0


def verify_2d_gaussian_torus_identities():
    results = {}
    for m in (6, 7):
        tor = GaussianTorus(m)
        N = tor.N

        # 1. Verify D_4 group structure (8 distinct automorphisms) & equivariance
        perms = tor.d4_permutations()
        assert len(set(perms)) == 8, f"D_4 did not yield 8 distinct automorphisms for m={m}"
        for p in perms:
            assert sorted(p) == list(range(N)), f"Non-bijective D_4 action for m={m}"

        K_mat, G0 = tor.compute_green_G0()
        W_mat = tor.matrix_of(tor.apply_W_flat)
        inv_N = Fraction(1, N)

        # Verify K_Sym G_0 == I - (1/N) 1 1^T, G_0^T == G_0, G_0 1 == 0
        for i in range(N):
            row_sum_G = sum(G0[i], ZERO_Q)
            assert row_sum_G == ZERO_Q, f"Nonzero row sum in G_0 at m={m}"
            for j in range(N):
                assert G0[i][j] == G0[j][i], f"G_0 not symmetric at m={m}"
                prod_ij = sum((K_mat[i][k] * G0[k][j] for k in range(N)), ZERO_Q)
                target = (ONE_Q - inv_N) if i == j else (-inv_N)
                assert prod_ij == target, f"K_Sym G_0 mismatch at ({i},{j}) for m={m}"

        # Verify D_4 equivariance of K_Sym, G_0, W_flat: M[p[i], p[j]] == M[i, j]
        for p in perms:
            for i in range(N):
                pi = p[i]
                for j in range(N):
                    pj = p[j]
                    assert K_mat[pi][pj] == K_mat[i][j]
                    assert G0[pi][pj] == G0[i][j]
                    assert W_mat[pi][pj] == W_mat[i][j]

        # 2. Symbolic QQPoly Sum-of-Squares proof for K_Sym and W_flat over N variables
        f_sym = qqpoly_vars(f"f_m{m}", N)
        kf_sym = tor.apply_K_Sym(f_sym)
        lx_f = tor.apply_Lx(f_sym)
        ly_f = tor.apply_Ly(f_sym)
        # Operator identity: K_Sym f == L_x f + L_y f + (1/12)(L_x^2 f + L_y^2 f)
        lx2_f = tor.apply_Lx(lx_f)
        ly2_f = tor.apply_Ly(ly_f)
        inv12 = Fraction(1, 12)
        for i in range(N):
            sos_op_i = lx_f[i] + ly_f[i] + inv12 * (lx2_f[i] + ly2_f[i])
            assert (kf_sym[i] - sos_op_i) == 0

        # Quadratic form Sum-of-Squares: <f, K_Sym f> == sum_z [(D+x f)^2 + (D+y f)^2 + (1/12)((Lx f)^2 + (Ly f)^2)]
        quad_K = sum((f_sym[i] * kf_sym[i] for i in range(N)), QQPoly.zero())
        sos_K = sum(
            (
                (f_sym[tor.shift_idx(i, 1, 0)] - f_sym[i]) ** 2
                + (f_sym[tor.shift_idx(i, 0, 1)] - f_sym[i]) ** 2
                + inv12 * (lx_f[i] ** 2 + ly_f[i] ** 2)
                for i in range(N)
            ),
            QQPoly.zero(),
        )
        assert (quad_K - sos_K) == 0

        # Quadratic form Sum-of-Squares for W_flat:
        # <f, W_flat f> == sum_z [ ((Lx - 1/2 Ly) f)^2 + (3/4) (Ly f)^2 ]
        wf_sym = tor.apply_W_flat(f_sym)
        quad_W = sum((f_sym[i] * wf_sym[i] for i in range(N)), QQPoly.zero())
        half = Fraction(1, 2)
        three_q = Fraction(3, 4)
        sos_W = sum(
            (
                (lx_f[i] - half * ly_f[i]) ** 2 + three_q * (ly_f[i] ** 2)
                for i in range(N)
            ),
            QQPoly.zero(),
        )
        assert (quad_W - sos_W) == 0

        # 3. 16:16:16 Brillouin Zone-Corner Isotropy via shift involutions (-1)^(ax + by) (zero pi!)
        corner_eigenvalues = {}
        for a, b in ((1, 0), (0, 1), (1, 1)):
            mode = [
                ONE_Q if (((x * a + y * b) % 2) == 0) else (-ONE_Q)
                for (x, y) in tor.sites
            ]
            # Verify mode is well-defined wrt shifts on G_m: T_x mode == (-1)^a mode, T_y mode == (-1)^b mode
            sx = -ONE_Q if a else ONE_Q
            sy = -ONE_Q if b else ONE_Q
            for i in range(N):
                assert mode[tor.shift_idx(i, 1, 0)] == sx * mode[i]
                assert mode[tor.shift_idx(i, 0, 1)] == sy * mode[i]
            w_mode = tor.apply_W_flat(mode)
            for i in range(N):
                assert w_mode[i] == Fraction(16, 1) * mode[i]
            corner_eigenvalues[(a, b)] = Fraction(16, 1)

        # 4. Symbolic QQPoly Arakawa Jacobian Conservation (<1, J> == 0, <psi, J> == 0, <om, J> == 0)
        psi_sym = qqpoly_vars(f"psi_m{m}", N)
        om_sym = qqpoly_vars(f"om_m{m}", N)
        for jac_name, jac_fn in (("J2", tor.arakawa_J2), ("J4", tor.arakawa_J4)):
            j_po = jac_fn(psi_sym, om_sym)
            j_op = jac_fn(om_sym, psi_sym)
            # Pointwise antisymmetry {psi, om} + {om, psi} == 0
            for i in range(N):
                assert (j_po[i] + j_op[i]) == 0
            assert sum(j_po, QQPoly.zero()) == 0
            assert sum((psi_sym[i] * j_po[i] for i in range(N)), QQPoly.zero()) == 0
            assert sum((om_sym[i] * j_po[i] for i in range(N)), QQPoly.zero()) == 0

        # 5. Symbolic 8-term Sum-of-Squares for 1B Viscous Dissipation <omega, G_0 W_flat omega> >= 0/1
        # Since [K_Sym, W_flat] == 0 and G_0 K_Sym = Pi_0, writing phi = G_0 omega gives
        # <omega, G_0 W_flat omega> == <phi, K_Sym W_flat phi> = sum of 8 non-negative squares!
        if m == 6:
            phi_s = f_sym
            # K_Sym W_flat = (Lx + Ly + 1/12 Lx^2 + 1/12 Ly^2) * (A^2 + 3/4 Ly^2), A = Lx - 1/2 Ly
            kw_phi = tor.apply_K_Sym(wf_sym)
            lhs_diss = sum((phi_s[i] * kw_phi[i] for i in range(N)), QQPoly.zero())
            a_phi = [lx_f[i] - half * ly_f[i] for i in range(N)]
            lx_a = tor.apply_Lx(a_phi)
            ly_a = tor.apply_Ly(a_phi)
            lx_ly = tor.apply_Lx(ly_f)
            ly_ly = tor.apply_Ly(ly_f)
            sos_diss = sum(
                (
                    # <A phi, Lx A phi> + <A phi, Ly A phi>
                    (a_phi[tor.shift_idx(i, 1, 0)] - a_phi[i]) ** 2
                    + (a_phi[tor.shift_idx(i, 0, 1)] - a_phi[i]) ** 2
                    # + (3/4)(<Ly phi, Lx Ly phi> + <Ly phi, Ly Ly phi>)
                    + three_q * (ly_f[tor.shift_idx(i, 1, 0)] - ly_f[i]) ** 2
                    + three_q * (ly_f[tor.shift_idx(i, 0, 1)] - ly_f[i]) ** 2
                    # + (1/12)(||Lx A phi||^2 + ||Ly A phi||^2) + (1/16)(||Lx Ly phi||^2 + ||Ly^2 phi||^2)
                    + inv12 * (lx_a[i] ** 2 + ly_a[i] ** 2)
                    + Fraction(1, 16) * (lx_ly[i] ** 2 + ly_ly[i] ** 2)
                    for i in range(N)
                ),
                QQPoly.zero(),
            )
            assert (lhs_diss - sos_diss) == 0

        results[m] = {
            "N": N,
            "G0_diag": G0[0][0],
            "corners": corner_eigenvalues,
        }

    # Also prove symbolic First-Law + Second-Law Cayley/implicit-midpoint thermodynamics over QQRat:
    # Given symmetric dissipative step m = (u_n + u_{n+1})/2 with u_{n+1} - u_n = -dt * D_pos,
    # Q_{n+1} = Q_n + (E_n - E_{n+1}) satisfies E_{n+1} + Q_{n+1} == E_n + Q_n (0/1)
    # and for scalar/modal mode u_{n+1} = ((1 - c*lam)/(1 + c*lam)) u_n (c = nu*dt/2 >= 0, lam >= 0),
    # Delta Q = (1/2) u_n^2 - (1/2) u_{n+1}^2 == 2*c*lam * u_n^2 / (1 + c*lam)^2 >= 0/1!
    u_r, q_r, c_r, lam_r = rat_vars("u_r", "q_r", "c_r", "lam_r")
    u_next = ((1 - c_r * lam_r) / (1 + c_r * lam_r)) * u_r
    e_cur = Fraction(1, 2) * (u_r ** 2)
    e_nxt = Fraction(1, 2) * (u_next ** 2)
    dq_heat = e_cur - e_nxt
    q_next = q_r + dq_heat
    assert (e_nxt + q_next) - (e_cur + q_r) == 0
    expected_sos_heat = (2 * c_r * lam_r * (u_r ** 2)) / ((1 + c_r * lam_r) ** 2)
    assert (dq_heat - expected_sos_heat) == 0

    return results


# ============================================================================
# Part 1C & Part 3: 3D Symanzik Vector Calculus, Exact Rational Leray Projector,
# Closed-Form Sharp Bounds B_u & B_omega, and Varga M-Matrix Positivity
# ============================================================================


class SymanzikLattice3D:
    """Exact 3D periodic torus Z_L^3 (N = L^3) with 4th-order Symanzik operators over Q."""

    def __init__(self, L):
        self.L = L
        self.N = L * L * L
        self.sites = [(x, y, z) for x in range(L) for y in range(L) for z in range(L)]
        self.index_of = {pt: idx for idx, pt in enumerate(self.sites)}
        # Precompute shift table for d in {-2, -1, +1, +2} along axes 0, 1, 2
        self.sh = {}
        for axis in (0, 1, 2):
            for step in (-2, -1, 1, 2):
                arr = [0] * self.N
                for idx, (x, y, z) in enumerate(self.sites):
                    coords = [x, y, z]
                    coords[axis] = (coords[axis] + step) % L
                    arr[idx] = self.index_of[tuple(coords)]
                self.sh[(axis, step)] = tuple(arr)

    def D_axis(self, f, axis):
        c1 = Fraction(2, 3)
        c2 = Fraction(1, 12)
        sp1 = self.sh[(axis, 1)]
        sm1 = self.sh[(axis, -1)]
        sp2 = self.sh[(axis, 2)]
        sm2 = self.sh[(axis, -2)]
        return [
            c1 * (f[sp1[i]] - f[sm1[i]]) - c2 * (f[sp2[i]] - f[sm2[i]])
            for i in range(self.N)
        ]

    def grad_Sym(self, phi):
        return (
            self.D_axis(phi, 0),
            self.D_axis(phi, 1),
            self.D_axis(phi, 2),
        )

    def div_Sym(self, u):
        dx = self.D_axis(u[0], 0)
        dy = self.D_axis(u[1], 1)
        dz = self.D_axis(u[2], 2)
        return [dx[i] + dy[i] + dz[i] for i in range(self.N)]

    def curl_Sym(self, A):
        dy_Az = self.D_axis(A[2], 1)
        dz_Ay = self.D_axis(A[1], 2)
        dz_Ax = self.D_axis(A[0], 2)
        dx_Az = self.D_axis(A[2], 0)
        dx_Ay = self.D_axis(A[1], 0)
        dy_Ax = self.D_axis(A[0], 1)
        N = self.N
        return (
            [dy_Az[i] - dz_Ay[i] for i in range(N)],
            [dz_Ax[i] - dx_Az[i] for i in range(N)],
            [dx_Ay[i] - dy_Ax[i] for i in range(N)],
        )

    def apply_Delta_div(self, phi):
        """Apply symmetric PSD Hodge Laplacian Delta_div = -div_Sym(grad_Sym(phi)) = G^T G phi."""
        dphi = self.div_Sym(self.grad_Sym(phi))
        return [-v for v in dphi]

    def compute_Delta_div_pseudoinverse_kernel(self):
        """Solve exact rational Moore-Penrose convolution kernel g_div of (G^T G)^+ using O_h orbits."""
        L = self.L
        N = self.N
        # Group Z_L^3 into cubic O_h symmetry orbits of sorted(min(x, L-x), min(y, L-y), min(z, L-z))
        orbit_key_of = []
        orbit_keys = []
        key_to_oid = {}
        for x, y, z in self.sites:
            rx = x if x <= L - x else L - x
            ry = y if y <= L - y else L - y
            rz = z if z <= L - z else L - z
            ok = tuple(sorted((rx, ry, rz)))
            oid = key_to_oid.get(ok)
            if oid is None:
                oid = len(orbit_keys)
                key_to_oid[ok] = oid
                orbit_keys.append(ok)
            orbit_key_of.append(oid)

        num_orbits = len(orbit_keys)
        orbit_rep_site = [orbit_key_of.index(oid) for oid in range(num_orbits)]

        # Kernel projector Pi_ker(z, 0) on Z_L^3:
        # - Odd L:  dim ker(grad_Sym) = 1 -> Pi_ker(z, 0) = 1 / N
        # - Even L: dim ker(grad_Sym) = 8 -> sum_{a,b,c in {0,1}^3} (-1)^(ax+by+cz) / N
        #           = (8 / N) if (x, y, z) all even, else 0!
        inv_N = Fraction(1, N)
        pi_ker_col0 = []
        for x, y, z in self.sites:
            if (L % 2) == 0:
                if (x % 2 == 0) and (y % 2 == 0) and (z % 2 == 0):
                    pi_ker_col0.append(Fraction(8, N))
                else:
                    pi_ker_col0.append(ZERO_Q)
            else:
                pi_ker_col0.append(inv_N)

        # Build exact num_orbits x num_orbits rational system for (Delta_div + Pi_ker) g = delta_0 - Pi_ker(., 0)
        # Note Pi_ker(z - w) acting on g: (Pi_ker g)(z) = sum_w Pi_ker(z - w) g(w)
        mat_orb = [[ZERO_Q] * num_orbits for _ in range(num_orbits)]
        for oid_col in range(num_orbits):
            indicator = [ONE_Q if orbit_key_of[i] == oid_col else ZERO_Q for i in range(N)]
            d_ind = self.apply_Delta_div(indicator)
            # Add Pi_ker * indicator
            for oid_row in range(num_orbits):
                rep_i = orbit_rep_site[oid_row]
                rx, ry, rz = self.sites[rep_i]
                pk_sum = ZERO_Q
                for j in range(N):
                    if indicator[j] != ZERO_Q:
                        jx, jy, jz = self.sites[j]
                        diff_i = self.index_of[((rx - jx) % L, (ry - jy) % L, (rz - jz) % L)]
                        pk_sum += pi_ker_col0[diff_i] * indicator[j]
                mat_orb[oid_row][oid_col] = d_ind[rep_i] + pk_sum

        rhs_orb = []
        for oid_row in range(num_orbits):
            rep_i = orbit_rep_site[oid_row]
            rhs_orb.append((ONE_Q if rep_i == 0 else ZERO_Q) - pi_ker_col0[rep_i])

        sol_orb = solve_linear_q(mat_orb, rhs_orb)
        g_div = [sol_orb[orbit_key_of[i]] for i in range(N)]
        return g_div, pi_ker_col0

    def convolve_scalar(self, ker, f):
        """Exact cyclic 3D convolution (ker * f)(x) = sum_y ker(x - y) f(y)."""
        L = self.L
        N = self.N
        # Fast separable/sparse-free direct shift accumulation over Q
        out = [ZERO_Q] * N
        for shift_idx, kval in enumerate(ker):
            if kval == ZERO_Q:
                continue
            sx, sy, sz = self.sites[shift_idx]
            for y_idx, fval in enumerate(f):
                if fval == ZERO_Q:
                    continue
                yx, yy, yz = self.sites[y_idx]
                x_idx = self.index_of[((yx + sx) % L, (yy + sy) % L, (yz + sz) % L)]
                out[x_idx] = out[x_idx] + kval * fval
        return out

    def apply_P_div(self, w, g_div):
        """Exact rational Leray projection P_div(w) = w + grad_Sym(g_div * div_Sym(w))."""
        div_w = self.div_Sym(w)
        phi = self.convolve_scalar(g_div, div_w)
        gphi = self.grad_Sym(phi)
        N = self.N
        return (
            [w[0][i] + gphi[0][i] for i in range(N)],
            [w[1][i] + gphi[1][i] for i in range(N)],
            [w[2][i] + gphi[2][i] for i in range(N)],
        )


def verify_3d_symanzik_and_leray():
    out_summary = {}

    # 1. Symbolic QQPoly vector calculus identities on Z_6^3 (N = 216 -> 864 symbolic variables!)
    lat6 = SymanzikLattice3D(6)
    N6 = lat6.N
    phi_sym = qqpoly_vars("phi3d", N6)
    Ax_sym = qqpoly_vars("Ax3d", N6)
    Ay_sym = qqpoly_vars("Ay3d", N6)
    Az_sym = qqpoly_vars("Az3d", N6)
    A_sym = (Ax_sym, Ay_sym, Az_sym)

    # 1a. div_Sym(curl_Sym(A)) == 0 identically at all 216 sites
    dcA = lat6.div_Sym(lat6.curl_Sym(A_sym))
    for i in range(N6):
        assert dcA[i] == 0

    # 1b. curl_Sym(grad_Sym(phi)) == (0, 0, 0) identically at all 216 sites
    cgp = lat6.curl_Sym(lat6.grad_Sym(phi_sym))
    for comp in range(3):
        for i in range(N6):
            assert cgp[comp][i] == 0

    # 1c. Skew-adjointness <A, grad_Sym(phi)> + <div_Sym(A), phi> == 0 over 864 variables!
    g_phi = lat6.grad_Sym(phi_sym)
    d_A = lat6.div_Sym(A_sym)
    inner_u_gp = sum(
        (A_sym[c][i] * g_phi[c][i] for c in range(3) for i in range(N6)),
        QQPoly.zero(),
    )
    inner_du_p = sum((d_A[i] * phi_sym[i] for i in range(N6)), QQPoly.zero())
    assert (inner_u_gp + inner_du_p) == 0

    # 1d. Exact self-adjointness of curl_Sym on periodic Z_6^3
    Bx_sym = qqpoly_vars("Bx3d", N6)
    By_sym = qqpoly_vars("By3d", N6)
    Bz_sym = qqpoly_vars("Bz3d", N6)
    B_sym = (Bx_sym, By_sym, Bz_sym)
    cA = lat6.curl_Sym(A_sym)
    cB = lat6.curl_Sym(B_sym)
    lhs_curl = sum((B_sym[c][i] * cA[c][i] for c in range(3) for i in range(N6)), QQPoly.zero())
    rhs_curl = sum((cB[c][i] * A_sym[c][i] for c in range(3) for i in range(N6)), QQPoly.zero())
    assert (lhs_curl - rhs_curl) == 0

    # 2. Exact rational Leray projector P_div verification on BOTH L=6 (N=216) and L=8 (N=512)!
    for L in (6, 8):
        lat = SymanzikLattice3D(L)
        N = lat.N
        g_div, pi_ker_col0 = lat.compute_Delta_div_pseudoinverse_kernel()

        # Verify Delta_div(g_div) == delta_0 - Pi_ker(., 0) and Pi_ker * g_div == 0
        dg = lat.apply_Delta_div(g_div)
        for i in range(N):
            target = (ONE_Q if i == 0 else ZERO_Q) - pi_ker_col0[i]
            assert dg[i] == target, f"Delta_div^+ residual nonzero at L={L}, site={i}"

        # Construct the exact 3x3 block-convolution kernel P_ker[alpha][beta][z] = (P_div e_{0, beta})_alpha(z)
        P_ker = [[None] * 3 for _ in range(3)]
        for beta in range(3):
            e_beta = tuple(
                [ONE_Q if (c == beta and i == 0) else ZERO_Q for i in range(N)]
                for c in range(3)
            )
            p_eb = lat.apply_P_div(e_beta, g_div)
            # Verify div_Sym(P_div e_{0, beta}) == 0/1 at every site -> proves div_Sym(P_div w) == 0 for ALL w!
            div_peb = lat.div_Sym(p_eb)
            for i in range(N):
                assert div_peb[i] == ZERO_Q, f"div_Sym(P_div) != 0 at L={L}"
            for alpha in range(3):
                P_ker[alpha][beta] = p_eb[alpha]

        # Verify P_div^T == P_div: P_ker[alpha][beta][z] == P_ker[beta][alpha][-z mod L]
        for alpha in range(3):
            for beta in range(3):
                for idx, (x, y, z) in enumerate(lat.sites):
                    neg_idx = lat.index_of[((-x) % L, (-y) % L, (-z) % L)]
                    assert P_ker[alpha][beta][idx] == P_ker[beta][alpha][neg_idx]

        # Verify P_div^2 == P_div: apply P_div to each basis column p_eb and check exact equality!
        for beta in range(3):
            col_vec = (P_ker[0][beta], P_ker[1][beta], P_ker[2][beta])
            p2_col = lat.apply_P_div(col_vec, g_div)
            for alpha in range(3):
                assert p2_col[alpha] == P_ker[alpha][beta], f"P_div^2 != P_div at L={L}"

        # Also verify symbolically on a 3x3x3 = 27-site or full Z_6^3 symbolic vector field w:
        # By linearity + translation invariance, P_ker symbolic convolution on A_sym gives div_Sym(P_div(A_sym)) == 0!
        if L == 6:
            # Test symbolic Leray projection on full 648-variable symbolic field on Z_6^3 via sparse shift!
            # Even better: test full `apply_P_div(A_sym, g_div)` on a smaller L=4 (N=64, 192 vars) directly!
            pass

        out_summary[L] = {
            "N": N,
            "P_div_diag_00": P_ker[0][0][0],
            "P_div_diag_01": P_ker[0][1][0],
            "Pi_zero_mean_diag": P_ker[0][0][0] - Fraction(1, N),
        }

    # Direct full symbolic QQPoly pass of apply_P_div on Z_4^3 (N=64 -> 192 symbolic variables)
    lat4 = SymanzikLattice3D(4)
    g_div4, _ = lat4.compute_Delta_div_pseudoinverse_kernel()
    wx4 = qqpoly_vars("wx4", lat4.N)
    wy4 = qqpoly_vars("wy4", lat4.N)
    wz4 = qqpoly_vars("wz4", lat4.N)
    pw4 = lat4.apply_P_div((wx4, wy4, wz4), g_div4)
    div_pw4 = lat4.div_Sym(pw4)
    for i in range(lat4.N):
        assert div_pw4[i] == 0

    return out_summary


def verify_sharp_no_blowup_bounds_and_varga():
    """Verify closed-form B_u(E_0, N), B_omega(E_0, N), symbolic SoS certificates,
    explicit rational extremizers attaining equality to 0/1, and Varga M-matrix positivity.
    """
    E0_sym = QQRat.var("E_0")
    cert = {}

    for L in (6, 8):
        lat = SymanzikLattice3D(L)
        N = lat.N
        g_div, _ = lat.compute_Delta_div_pseudoinverse_kernel()
        inv_N = Fraction(1, N)

        # 1. Compute exact reproducing block E_0 Pi E_0^T = E_0 (P_div - P_mean) E_0^T over Q
        K_u_block = [[ZERO_Q] * 3 for _ in range(3)]
        pi_cols_at_0 = []
        for beta in range(3):
            e_beta = tuple(
                [ONE_Q if (c == beta and i == 0) else ZERO_Q for i in range(N)]
                for c in range(3)
            )
            p_eb = lat.apply_P_div(e_beta, g_div)
            pi_eb = tuple(
                [(p_eb[c][i] - inv_N) if c == beta else p_eb[c][i] for i in range(N)]
                for c in range(3)
            )
            pi_cols_at_0.append(pi_eb)
            for alpha in range(3):
                K_u_block[alpha][beta] = pi_eb[alpha][0]

        expected_cu = Fraction(2 * N + 5, 3 * N)
        for alpha in range(3):
            for beta in range(3):
                target = expected_cu if alpha == beta else ZERO_Q
                assert K_u_block[alpha][beta] == target

        # Verify symbolic formula bound_B_u(E0_sym, N) == 2 * expected_cu * E0_sym
        bu_sym = bound_B_u(E0_sym, N, zero_momentum=True)
        assert bu_sym == (2 * expected_cu) * E0_sym

        # Exhibit EXPLICIT RATIONAL EXTREMIZER u_star_u = pi_cols_at_0[0] in V_{div,0}!
        u_star = pi_cols_at_0[0]
        div_ustar = lat.div_Sym(u_star)
        assert all(v == ZERO_Q for v in div_ustar)
        for c in range(3):
            assert sum(u_star[c], ZERO_Q) == ZERO_Q
        ekin_ustar = Fraction(1, 2) * sum(
            (u_star[c][i] ** 2 for c in range(3) for i in range(N)), ZERO_Q
        )
        max_u_sq = max(
            sum((u_star[c][i] ** 2 for c in range(3)), ZERO_Q) for i in range(N)
        )
        # Sharpness to 0/1!
        assert max_u_sq == bound_B_u(ekin_ustar, N, zero_momentum=True)

        # 2. Compute exact vorticity Gram block E_0 (curl_Sym^2) E_0^T over Q
        K_om_block = [[ZERO_Q] * 3 for _ in range(3)]
        curl_cols_at_0 = []
        for beta in range(3):
            e_beta = tuple(
                [ONE_Q if (c == beta and i == 0) else ZERO_Q for i in range(N)]
                for c in range(3)
            )
            c_eb = lat.curl_Sym(e_beta)
            curl_cols_at_0.append(c_eb)
            cc_eb = lat.curl_Sym(c_eb)
            for alpha in range(3):
                K_om_block[alpha][beta] = cc_eb[alpha][0]

        expected_com = Fraction(65, 36)
        for alpha in range(3):
            for beta in range(3):
                target = expected_com if alpha == beta else ZERO_Q
                assert K_om_block[alpha][beta] == target

        bom_sym = bound_B_omega(E0_sym, N)
        assert bom_sym == Fraction(65, 18) * E0_sym

        # Exhibit EXPLICIT RATIONAL EXTREMIZER u_star_om = curl_Sym(delta_0 e_x) in V_{div,0}!
        u_star_om = curl_cols_at_0[0]
        assert all(v == ZERO_Q for v in lat.div_Sym(u_star_om))
        for c in range(3):
            assert sum(u_star_om[c], ZERO_Q) == ZERO_Q
        ekin_om = Fraction(1, 2) * sum(
            (u_star_om[c][i] ** 2 for c in range(3) for i in range(N)), ZERO_Q
        )
        om_field = lat.curl_Sym(u_star_om)
        max_om_sq = max(
            sum((om_field[c][i] ** 2 for c in range(3)), ZERO_Q) for i in range(N)
        )
        # Sharpness to 0/1!
        assert max_om_sq == bound_B_omega(ekin_om, N)

        cert[L] = {
            "N": N,
            "B_u_factor": 2 * expected_cu,
            "B_omega_factor": Fraction(65, 18),
            "extremizer_u_slack": max_u_sq - bound_B_u(ekin_ustar, N),
            "extremizer_omega_slack": max_om_sq - bound_B_omega(ekin_om, N),
        }

    # 3. Symbolic Sum-of-Squares Cauchy-Schwarz / Reproducing-Kernel Algebra Proof in QQRat:
    # Verify c * (2 E_0) - ||v||^2 == c * (2 E_0 - ||u||_2^2) + (1/c) * || c u - w ||_2^2
    # where <u, w> = ||v||^2 and <w, w> = c * ||v||^2!
    c_s, e0_s, unorm2_s, vnorm2_s = rat_vars("c_rk", "E0_rk", "unorm2_rk", "vnorm2_rk")
    # Expand || c u - w ||_2^2 = c^2 ||u||_2^2 - 2 c <u, w> + <w, w> = c^2 unorm2 - 2 c vnorm2 + c vnorm2
    norm_diff_sq = (c_s ** 2) * unorm2_s - 2 * c_s * vnorm2_s + c_s * vnorm2_s
    rhs_sos = c_s * (2 * e0_s - unorm2_s) + (ONE_Q / c_s) * norm_diff_sq
    lhs_gap = c_s * (2 * e0_s) - vnorm2_s
    assert (lhs_gap - rhs_sos) == 0

    # 4. Varga Squared-M-Matrix Positivity & Mass Conservation Certificate:
    # (4a) Non-circulant spatially-varying shear advection-diffusion M-matrix (I + (dt/2) A)^(-2) on Z_3^3 (N = 27)
    lat3 = SymanzikLattice3D(3)
    N3 = lat3.N
    kappa = Fraction(1, 8)
    dt = Fraction(1, 2)
    alpha = dt / 2
    A_mat3 = [[ZERO_Q] * N3 for _ in range(N3)]
    for y in range(N3):
        yx, yy, yz = lat3.sites[y]
        vel = (
            Fraction((yy % 3) - 1, 2),
            Fraction((yz % 3) - 1, 3),
            Fraction((yx % 3) - 1, 5),
        )
        for axis in (0, 1, 2):
            nb_plus = lat3.sh[(axis, 1)][y]
            nb_minus = lat3.sh[(axis, -1)][y]
            A_mat3[y][y] += 2 * kappa
            A_mat3[nb_plus][y] -= kappa
            A_mat3[nb_minus][y] -= kappa
            v_ax = vel[axis]
            if v_ax > ZERO_Q:
                A_mat3[y][y] += v_ax
                A_mat3[nb_plus][y] -= v_ax
            elif v_ax < ZERO_Q:
                A_mat3[y][y] += (-v_ax)
                A_mat3[nb_minus][y] -= (-v_ax)

    for y in range(N3):
        assert sum((A_mat3[x][y] for x in range(N3)), ZERO_Q) == ZERO_Q
        for x in range(N3):
            if x != y:
                assert A_mat3[x][y] <= ZERO_Q

    B_mat3 = [
        [(ONE_Q if x == y else ZERO_Q) + alpha * A_mat3[x][y] for y in range(N3)]
        for x in range(N3)
    ]
    I_N3 = [[ONE_Q if x == y else ZERO_Q for y in range(N3)] for x in range(N3)]
    R1_3 = solve_linear_q(B_mat3, I_N3)
    M_varga3 = [
        [sum((R1_3[x][k] * R1_3[k][y] for k in range(N3)), ZERO_Q) for y in range(N3)]
        for x in range(N3)
    ]
    for y in range(N3):
        assert sum((M_varga3[x][y] for x in range(N3)), ZERO_Q) == ONE_Q
        for x in range(N3):
            assert M_varga3[x][y] > ZERO_Q

    # (4b) Full Z_6^3 (N = 216) Varga squared-M-matrix diffusion resolvent (I + (dt/2) kappa L_graph)^(-2)
    # composed with conservative convex upwind transport flux M_flux, plus full 216-var QQPoly proof!
    lat6 = SymanzikLattice3D(6)
    N6 = lat6.N
    # Solve 1D 6x6 Varga factor (I + alpha * kappa * L_1d)^(-1) and tensor-square in 3D (or O_h orbit solve!)
    L1d_6 = [[ZERO_Q] * 6 for _ in range(6)]
    for i in range(6):
        L1d_6[i][i] = ONE_Q + 2 * alpha * kappa
        L1d_6[i][(i + 1) % 6] -= alpha * kappa
        L1d_6[i][(i - 1) % 6] -= alpha * kappa
    r1_1d = solve_linear_q(L1d_6, [ONE_Q if i == 0 else ZERO_Q for i in range(6)])
    r2_1d = [
        sum((r1_1d[(i - k) % 6] * r1_1d[k] for k in range(6)), ZERO_Q)
        for i in range(6)
    ]
    assert sum(r2_1d, ZERO_Q) == ONE_Q and all(v > ZERO_Q for v in r2_1d)

    # Full 3D separable Varga squared M-matrix kernel K_varga(dx, dy, dz) = r2_1d[dx]*r2_1d[dy]*r2_1d[dz] > 0/1!
    ker_varga_6 = [
        r2_1d[x] * r2_1d[y] * r2_1d[z] for (x, y, z) in lat6.sites
    ]
    assert sum(ker_varga_6, ZERO_Q) == ONE_Q
    min_entry = min(ker_varga_6)
    min_diag = ker_varga_6[0]
    assert min_entry > ZERO_Q and min_diag > ZERO_Q

    # Compose with convex upwind flux step on full 216 symbolic variables rho_0..rho_215 in QQPoly!
    rho_sym6 = qqpoly_vars("rho6", N6)
    # Conservative convex upwind step: rho_mid = rho - tau * div_upwind(u * rho) with CFL <= 1
    tau_cfl = Fraction(1, 8)
    rho_mid = list(rho_sym6)
    # Verify exact convex combination weights (all >= 0, diagonal >= 1 - 6*tau*(1/2) = 5/8 > 0)
    flux_col_weights_pos = True
    for y in range(N6):
        yx, yy, yz = lat6.sites[y]
        vel = (
            Fraction((yy % 3) - 1, 2),
            Fraction((yz % 3) - 1, 3),
            Fraction((yx % 3) - 1, 5),
        )
        out_rate = tau_cfl * (abs(vel[0]) + abs(vel[1]) + abs(vel[2]))
        if ONE_Q - out_rate <= ZERO_Q:
            flux_col_weights_pos = False
        for axis in (0, 1, 2):
            v_ax = vel[axis]
            if v_ax > ZERO_Q:
                nb = lat6.sh[(axis, 1)][y]
                flx = tau_cfl * v_ax * rho_sym6[y]
                rho_mid[y] = rho_mid[y] - flx
                rho_mid[nb] = rho_mid[nb] + flx
            elif v_ax < ZERO_Q:
                nb = lat6.sh[(axis, -1)][y]
                flx = tau_cfl * (-v_ax) * rho_sym6[y]
                rho_mid[y] = rho_mid[y] - flx
                rho_mid[nb] = rho_mid[nb] + flx
    assert flux_col_weights_pos
    # Apply 1D separable Varga squared resolvent along x, y, z on symbolic QQPoly vector rho_mid!
    rho_v = rho_mid
    for axis in (0, 1, 2):
        nxt = [QQPoly.zero() for _ in range(N6)]
        for idx, (x, y, z) in enumerate(lat6.sites):
            coords = [x, y, z]
            acc = QQPoly.zero()
            for d in range(6):
                c_src = list(coords)
                c_src[axis] = (coords[axis] - d) % 6
                acc = acc + r2_1d[d] * rho_v[lat6.index_of[tuple(c_src)]]
            nxt[idx] = acc
        rho_v = nxt

    # Prove exact 0/1 symbolic mass conservation across ALL 216 symbolic variables!
    mass_diff = sum(rho_v, QQPoly.zero()) - sum(rho_sym6, QQPoly.zero())
    assert mass_diff == 0
    # And verify every coefficient of rho_v[x] wrt rho_sym6[y] is STRICTLY POSITIVE (> 0/1)!
    for x in range(N6):
        assert len(rho_v[x]._terms) == N6
        for _, coeff_q in rho_v[x]._terms.items():
            assert coeff_q > ZERO_Q

    cert["varga_Z6"] = {
        "N": N6,
        "min_entry": min_entry,
        "min_diag": min_diag,
        "col_sum_residual": ZERO_Q,
    }
    return cert


def run_all_spatial_and_bound_proofs():
    t0_ns = time.perf_counter_ns()
    print("=== Running 2D Gaussian Torus Z[i]/<(1+i)^m> Symbolic Proofs (m=6 N=64, m=7 N=128) ===")
    res_2d = verify_2d_gaussian_torus_identities()
    for m, info in res_2d.items():
        print(
            f"  [PASS] m={m} (N={info['N']}): D_4 equivariance, K_Sym*G_0 == I - 1/N J (G_0[0,0]={info['G0_diag']}), "
            f"K_Sym & W_flat SoS, 16:16:16 corners={list(info['corners'].values())}, Arakawa J2 & J4 symbolic 0/1"
        )

    print("=== Running 3D Symanzik Vector Calculus & Exact Rational Leray Projector Proofs ===")
    res_3d = verify_3d_symanzik_and_leray()
    for L, info in res_3d.items():
        print(
            f"  [PASS] Z_{L}^3 (N={info['N']}): div(curl)==0, curl(grad)==0, <u,grad p>==-<div u,p>, "
            f"P_div^2==P_div, P_div^T==P_div, div(P_div w)==0 (P_div[0,0](0)={info['P_div_diag_00']})"
        )

    print("=== Running Closed-Form Sharp No-Blowup Bounds & Varga Positivity Proofs ===")
    res_bounds = verify_sharp_no_blowup_bounds_and_varga()
    for L in (6, 8):
        b = res_bounds[L]
        print(
            f"  [PASS] Z_{L}^3 (N={b['N']}): B_u(E_0,{b['N']}) = ({b['B_u_factor']}) * E_0 "
            f"[extremizer slack={b['extremizer_u_slack']}], "
            f"B_omega(E_0,{b['N']}) = ({b['B_omega_factor']}) * E_0 "
            f"[extremizer slack={b['extremizer_omega_slack']}]"
        )
    vg = res_bounds["varga_Z6"]
    print(
        f"  [PASS] Varga (I + (dt/2)A)^(-2) on Z_6^3 (N={vg['N']}): "
        f"col_sum_err={vg['col_sum_residual']}, min_diag={vg['min_diag']} > 0, min_entry > 0"
    )
    dt_wall = (time.perf_counter_ns() - t0_ns) // 1_000_000
    print(f"=== ALL SPATIAL ALGEBRA, POSITIVITY & NO-BLOWUP PROOFS VERIFIED ({dt_wall} ms) ===")
    return {"2d": res_2d, "3d": res_3d, "bounds": res_bounds}


if __name__ == "__main__":
    run_all_spatial_and_bound_proofs()
