from __future__ import annotations

# =============================================================================
# PART I: EXACT ULTRAFINITIST Q-ENGINE (PURE RATIONAL ARITHMETIC OVER Q)
# =============================================================================

from fractions import Fraction

# --- Module: quantizer.py ---
ZERO_Q = Fraction(0, 1)
ONE_Q = Fraction(1, 1)
HALF_Q = Fraction(1, 2)

MANTISSA_BITS_DEFAULT = 56
MAX_SHIFT_DEFAULT = 256
CARRIER_MAX_SHIFT = 60
SNAP_MANTISSA_BITS = 60
SNAP_MAX_SHIFT = 60


def _lcm_pos_ints(a_int, b_int):
    """Exact LCM of positive integers `a_int, b_int` without `math` or `while`."""
    if a_int == b_int or (a_int % b_int) == 0:
        return a_int
    if (b_int % a_int) == 0:
        return b_int
    return a_int * Fraction(a_int, b_int).denominator


def pade_22(r):
    """Exact [2/2] Pade rational guide weight P_{2,2}(r) = (5r^2 + 10r + 1) / (r^2 + 10r + 5)."""
    r2 = r * r
    ten_r = 10 * r
    return (5 * r2 + ten_r + 1) / (r2 + ten_r + 5)


def floor_log2_positive_rational(anum, den):
    """Exact integer floor(log2(anum / den)) for integers anum > 0, den > 0 in O(1) bit ops."""
    eb = anum.bit_length() - den.bit_length()
    if eb >= 0:
        return (eb - 1) if (anum < (den << eb)) else eb
    return (eb - 1) if ((anum << (-eb)) < den) else eb


def identity_scalar(val):
    """Exact identity scalar pass-through for algebraic / symbolic evaluation."""
    return val


def identity_vector(vec):
    """Exact identity vector pass-through for algebraic / symbolic evaluation."""
    return list(vec)


def identity_vec3(vec3):
    """Exact identity 3-component vector pass-through for algebraic / symbolic evaluation."""
    return [list(vec3[0]), list(vec3[1]), list(vec3[2])]


def snap_factored_identity(vec_num, scalar_den=ONE_Q):
    """Exact identity factored pass-through `(vec_num, scalar_den)` without per-site division."""
    return vec_num, scalar_den


def clear_vec_identity(vec):
    """Exact identity pass-through `(list(vec), 1/1)` for symbolic evaluation."""
    return list(vec), ONE_Q


def clear_vec3_identity(vec3):
    """Exact identity 3-component pass-through `(vec3, 1/1)` for symbolic evaluation."""
    return [list(vec3[0]), list(vec3[1]), list(vec3[2])], ONE_Q


def clear_vec_to_int_or_jet(vec):
    """Scale `vec` by exact positive integer LCM of `.denominator`.

    Uses pairwise binary-tree LCM + GCD-free integer scaling (`d_lcm // q0`) so
    large unquantized rational vectors clear in milliseconds while preserving all
    higher Taylor jet coefficients `(c_1, c_2)` to `0/1` exactness on `QQJet`/`AuditJet`.
    """
    lvl = [val.denominator for val in vec]
    for lvl_idx in range(15):
        n_lvl = len(lvl)
        if n_lvl <= 1:
            break
        nxt = []
        n_pairs = n_lvl // 2
        for k_pair in range(n_pairs):
            nxt.append(_lcm_pos_ints(lvl[2 * k_pair], lvl[2 * k_pair + 1]))
        if n_lvl & 1:
            nxt.append(lvl[n_lvl - 1])
        lvl = nxt
    d_lcm = lvl[0] if lvl else 1
    big_den_cutoff = 1 << 160
    out = []
    for val in vec:
        q0 = val.denominator
        n0 = val.numerator
        mult = d_lcm // q0
        if q0 > big_den_cutoff:
            out.append(n0 * mult)
        else:
            prod_q0 = val * q0
            out.append(n0 * mult if prod_q0 == n0 else prod_q0 * mult)
    return out, Fraction(d_lcm, 1)


def clear_vec3_to_int_or_jet(vec3):
    """Scale `[vx, vy, vz]` by exact positive integer LCM of `.denominator`."""
    n0 = len(vec3[0])
    flat = list(vec3[0]) + list(vec3[1]) + list(vec3[2])
    c_flat, d_q = clear_vec_to_int_or_jet(flat)
    return [c_flat[0:n0], c_flat[n0 : 2 * n0], c_flat[2 * n0 : 3 * n0]], d_q


def clear_factored_vec_to_int_or_jet(vec_num, scalar_den=ONE_Q):
    """Exact 0/1 denominator clearing of `vec_num / scalar_den` onto `list[int]` (or exact Jet)."""
    c_vec, d_q = clear_vec_to_int_or_jet(vec_num)
    return c_vec, d_q * scalar_den


def clear_factored_vec3_to_int_or_jet(vec3_num, scalar_den=ONE_Q):
    """Exact 0/1 denominator clearing of `[vx, vy, vz] / scalar_den` onto `list[list[int]]` (or exact Jet)."""
    c_vec3, d_q = clear_vec3_to_int_or_jet(vec3_num)
    return c_vec3, d_q * scalar_den


def wrap_mod4_dyadic(angle_val):
    """Map `angle_val` modulo 4 into `[-2, 2)` via exact integer shift `angle_val - 4 * k`.

    Preserves all higher derivatives when `angle_val` is a Taylor jet (`QQJet` / `AuditJet`).
    """
    shifted = angle_val + 2
    k_int = shifted.numerator // (4 * shifted.denominator)
    return angle_val - 4 * k_int


_WRAP_MOD4_DISPATCH = (identity_scalar, wrap_mod4_dyadic)


def wrap_mod4_symmetric(angle_val, quantize=True):
    """Map `angle_val` modulo 4 into `[-2, 2)` via dispatch table."""
    fn = _WRAP_MOD4_DISPATCH[bool(quantize)]
    return fn(angle_val)


def quantize_dyadic_scalar(val, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=MAX_SHIFT_DEFAULT):
    """Odd-symmetric floating-dyadic scalar quantizer Q(-val) == -Q(val)."""
    num = val.numerator
    if num == 0:
        return ZERO_Q
    neg = num < 0
    anum = -num if neg else num
    den = val.denominator
    ef = floor_log2_positive_rational(anum, den)
    shift = mantissa_bits - ef
    if shift > max_shift:
        shift = max_shift
    if shift >= 0:
        scaled = (anum << (shift + 1)) + den
        q_int = scaled // (den << 1)
        signed_q = -q_int if neg else q_int
        return Fraction(signed_q, 1 << shift)
    rshift = -shift
    den_shifted = den << rshift
    scaled = (anum << 1) + den_shifted
    q_int = (scaled // (den_shifted << 1)) << rshift
    signed_q = -q_int if neg else q_int
    return Fraction(signed_q, 1)


def quantize_dyadic_vector(vec, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=MAX_SHIFT_DEFAULT):
    """Odd-symmetric shared-exponent floating-dyadic vector quantizer.

    Anchors all vector entries to a single dyadic grid `2^(e_max - mantissa_bits)`
    determined by `e_max = max_i floor(log2 |vec[i]|)`. Preserves antisymmetry
    `Q(-vec) == -Q(vec)` on all of `Q^n` (including half-way ties) while bounding
    every output denominator to a single power of 2 (`<= mantissa_bits + max_shift` bits).
    """
    n_items = len(vec)
    emax = None
    for idx in range(n_items):
        q = vec[idx]
        num = q.numerator
        if num != 0:
            anum = -num if num < 0 else num
            ef = floor_log2_positive_rational(anum, q.denominator)
            if emax is None or ef > emax:
                emax = ef
    if emax is None:
        return [ZERO_Q for idx in range(n_items)]

    shift = mantissa_bits - emax
    if shift > max_shift:
        shift = max_shift

    out = []
    if shift >= 0:
        scale_den = 1 << shift
        for idx in range(n_items):
            q = vec[idx]
            num = q.numerator
            if num == 0:
                out.append(ZERO_Q)
            else:
                neg = num < 0
                anum = -num if neg else num
                den = q.denominator
                scaled = (anum << (shift + 1)) + den
                q_int = scaled // (den << 1)
                out.append(Fraction(-q_int if neg else q_int, scale_den))
    else:
        rshift = -shift
        for idx in range(n_items):
            q = vec[idx]
            num = q.numerator
            if num == 0:
                out.append(ZERO_Q)
            else:
                neg = num < 0
                anum = -num if neg else num
                den_shifted = q.denominator << rshift
                scaled = (anum << 1) + den_shifted
                q_int = (scaled // (den_shifted << 1)) << rshift
                out.append(Fraction(-q_int if neg else q_int, 1))
    return out


def quantize_dyadic_vec3(vec3, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=MAX_SHIFT_DEFAULT):
    """Odd-symmetric shared-exponent floating-dyadic quantizer on `[ux, uy, uz]`."""
    n0 = len(vec3[0])
    flat = list(vec3[0]) + list(vec3[1]) + list(vec3[2])
    qflat = quantize_dyadic_vector(flat, mantissa_bits=mantissa_bits, max_shift=max_shift)
    return [qflat[0:n0], qflat[n0 : 2 * n0], qflat[2 * n0 : 3 * n0]]


def quantize_dyadic_floor_nonneg(val, grid_bits=MANTISSA_BITS_DEFAULT):
    """Conservative downward dyadic quantizer `0 <= Q_down(val) <= val` on fixed grid `2^-grid_bits`."""
    num = val.numerator
    if num <= 0:
        return ZERO_Q
    den = val.denominator
    q_int = (num << grid_bits) // den
    return Fraction(q_int, 1 << grid_bits)


def snap_factored_vec_dyadic(
    vec_num,
    scalar_den=ONE_Q,
    mantissa_bits=SNAP_MANTISSA_BITS,
    max_shift=SNAP_MAX_SHIFT,
):
    """GCD-free odd-symmetric shared-exponent integer snap of `vec_num / scalar_den`.

    Returns `(int_list, scale_den_q)` where `int_list` is a `list[int]` of bounded
    (`<= mantissa_bits + 1` bit) Python integers and `scale_den_q` is `Fraction(2^shift, 1)`,
    so `int_list[i] / scale_den_q` approximates `vec_num[i] / scalar_den` to `2^-mantissa_bits`
    while enabling pure C `PyLong` execution in downstream spatial stencils!
    """
    n_items = len(vec_num)
    sd_num_raw = scalar_den.numerator
    sd_den = scalar_den.denominator
    neg_sd = sd_num_raw < 0
    sd_num = -sd_num_raw if neg_sd else sd_num_raw

    emax = None
    for idx in range(n_items):
        q = vec_num[idx]
        num = q.numerator
        if num != 0:
            anum = (-num if num < 0 else num) * sd_den
            den = q.denominator * sd_num
            ef = floor_log2_positive_rational(anum, den)
            if emax is None or ef > emax:
                emax = ef
    if emax is None:
        return [0 for idx in range(n_items)], ONE_Q

    shift = mantissa_bits - emax
    if shift > max_shift:
        shift = max_shift

    out = []
    if shift >= 0:
        for idx in range(n_items):
            q = vec_num[idx]
            num = q.numerator
            if num == 0:
                out.append(0)
            else:
                neg = (num < 0) != neg_sd
                anum = (-num if num < 0 else num) * sd_den
                den = q.denominator * sd_num
                q_int = ((anum << (shift + 1)) + den) // (den << 1)
                out.append(-q_int if neg else q_int)
        return out, Fraction(1 << shift, 1)

    rshift = -shift
    for idx in range(n_items):
        q = vec_num[idx]
        num = q.numerator
        if num == 0:
            out.append(0)
        else:
            neg = (num < 0) != neg_sd
            anum = (-num if num < 0 else num) * sd_den
            den_shifted = (q.denominator * sd_num) << rshift
            q_int = (((anum << 1) + den_shifted) // (den_shifted << 1)) << rshift
            out.append(-q_int if neg else q_int)
    return out, ONE_Q


def snap_factored_vec3_dyadic(
    vec3_num,
    scalar_den=ONE_Q,
    mantissa_bits=SNAP_MANTISSA_BITS,
    max_shift=SNAP_MAX_SHIFT,
):
    """GCD-free odd-symmetric shared-exponent integer snap of `[vx, vy, vz] / scalar_den`."""
    n0 = len(vec3_num[0])
    flat = list(vec3_num[0]) + list(vec3_num[1]) + list(vec3_num[2])
    q_int_flat, scale_den_q = snap_factored_vec_dyadic(
        flat, scalar_den=scalar_den, mantissa_bits=mantissa_bits, max_shift=max_shift
    )
    return (
        [q_int_flat[0:n0], q_int_flat[n0 : 2 * n0], q_int_flat[2 * n0 : 3 * n0]],
        scale_den_q,
    )


def _raw_num_scale(q_val):
    """Extract raw integer numerator of an integer-valued scale `Fraction(k, 1)`."""
    return q_val.numerator


EXACT_IDENTITY_BUNDLE = {
    "wrap_mod4": identity_scalar,
    "vec": identity_vector,
    "vec3": identity_vec3,
    "nonneg": identity_scalar,
    "snap_vec": snap_factored_identity,
    "snap_vec3": snap_factored_identity,
    "clear_vec": clear_vec_identity,
    "clear_vec3": clear_vec3_identity,
    "to_raw_scale": identity_scalar,
}

EXACT_FACTORED_BUNDLE = {
    "wrap_mod4": identity_scalar,
    "vec": identity_vector,
    "vec3": identity_vec3,
    "nonneg": identity_scalar,
    "snap_vec": clear_factored_vec_to_int_or_jet,
    "snap_vec3": clear_factored_vec3_to_int_or_jet,
    "clear_vec": clear_vec_to_int_or_jet,
    "clear_vec3": clear_vec3_to_int_or_jet,
    "to_raw_scale": _raw_num_scale,
}

DYADIC_QUANTIZER_BUNDLE = {
    "wrap_mod4": wrap_mod4_dyadic,
    "vec": lambda v: quantize_dyadic_vector(
        v, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=CARRIER_MAX_SHIFT
    ),
    "vec3": lambda v3: quantize_dyadic_vec3(
        v3, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=CARRIER_MAX_SHIFT
    ),
    "nonneg": lambda val: quantize_dyadic_floor_nonneg(val, grid_bits=MANTISSA_BITS_DEFAULT),
    "snap_vec": snap_factored_vec_dyadic,
    "snap_vec3": snap_factored_vec3_dyadic,
    "clear_vec": clear_vec_to_int_or_jet,
    "clear_vec3": clear_vec3_to_int_or_jet,
    "to_raw_scale": _raw_num_scale,
}

_BUNDLE_DISPATCH = (EXACT_IDENTITY_BUNDLE, DYADIC_QUANTIZER_BUNDLE)


def resolve_quantizer_bundle(mode_flag=True, bundle=None):
    """Resolve quantizer callable bundle via table lookup without mode branching."""
    if bundle is not None:
        return bundle
    return _BUNDLE_DISPATCH[bool(mode_flag)]


# --- Module: spatial_2d.py ---
ZERO_Q = Fraction(0, 1)
ONE_Q = Fraction(1, 1)
HALF_Q = Fraction(1, 2)
THIRD_Q = Fraction(1, 3)
TWO_THIRDS_Q = Fraction(2, 3)
FOUR_THIRDS_Q = Fraction(4, 3)
TWELFTH_Q = Fraction(1, 12)
THREE_QUARTERS_Q = Fraction(3, 4)


def _lcm_positive_ints_q(a_int, b_int):
    """Exact LCM of positive integers `a_int, b_int` via `Fraction` reduction (zero `math` / `while`)."""
    return a_int * Fraction(a_int, b_int).denominator


def solve_linear_system_q(mat, rhs_vec):
    """Exact finite-step Gauss-Jordan elimination over Q (`mat * x = rhs_vec`)."""
    n_dim = len(mat)
    a_rows = [[Fraction(val) for val in row] for row in mat]
    b_vec = [Fraction(val) for val in rhs_vec]

    for col in range(n_dim):
        piv = None
        for row_idx in range(col, n_dim):
            if a_rows[row_idx][col] != ZERO_Q:
                piv = row_idx
                break
        if piv is None:
            raise ValueError(f"Singular matrix over Q at column {col}")
        if piv != col:
            a_rows[col], a_rows[piv] = a_rows[piv], a_rows[col]
            b_vec[col], b_vec[piv] = b_vec[piv], b_vec[col]

        inv_piv = ONE_Q / a_rows[col][col]
        for j_idx in range(col, n_dim):
            a_rows[col][j_idx] *= inv_piv
        b_vec[col] *= inv_piv

        for row_idx in range(n_dim):
            if row_idx != col and a_rows[row_idx][col] != ZERO_Q:
                factor = a_rows[row_idx][col]
                a_rows[row_idx][col] = ZERO_Q
                for j_idx in range(col + 1, n_dim):
                    a_rows[row_idx][j_idx] -= factor * a_rows[col][j_idx]
                b_vec[row_idx] -= factor * b_vec[col]

    return b_vec


def _gaussian_pow_1_plus_i(r_exp):
    """Compute `(1 + i)^r_exp` as exact integer pair `(re, im)`."""
    re_val = 1
    im_val = 0
    for step_idx in range(r_exp):
        next_re = re_val - im_val
        next_im = re_val + im_val
        re_val = next_re
        im_val = next_im
    return re_val, im_val


class GaussianTorusGrid2D:
    """Exact 2D Gaussian integer quotient lattice `Z[i] / <(1 + i)^m>` (`N = 2^m`)."""

    def __bool__(self):
        return True

    def __len__(self):
        return 0

    def __init__(self, m):
        if m < 4:
            raise ValueError("GaussianTorusGrid2D requires m >= 4 for range-2 stencils")
        self.m = m
        self.N = 1 << m
        if (m % 2) == 0:
            self.nx = 1 << (m // 2)
            self.ny = 1 << (m // 2)
            self.twist = 0
        else:
            k_half = (m - 1) // 2
            self.nx = 1 << (k_half + 1)
            self.ny = 1 << k_half
            self.twist = 1 << k_half

        self.coords = tuple(
            (x_coord, y_coord)
            for y_coord in range(self.ny)
            for x_coord in range(self.nx)
        )
        self.coord_to_idx = {pt: idx for idx, pt in enumerate(self.coords)}

        # Precompute range-1 and range-2 shift permutation tables
        self.xp1 = tuple(self.shift_index(i, 1, 0) for i in range(self.N))
        self.xm1 = tuple(self.shift_index(i, -1, 0) for i in range(self.N))
        self.xp2 = tuple(self.shift_index(i, 2, 0) for i in range(self.N))
        self.xm2 = tuple(self.shift_index(i, -2, 0) for i in range(self.N))
        self.yp1 = tuple(self.shift_index(i, 0, 1) for i in range(self.N))
        self.ym1 = tuple(self.shift_index(i, 0, -1) for i in range(self.N))
        self.yp2 = tuple(self.shift_index(i, 0, 2) for i in range(self.N))
        self.ym2 = tuple(self.shift_index(i, 0, -2) for i in range(self.N))

        # D_4 group generators: 90-degree rotation z -> i*z and conjugation z -> conj(z)
        self.rot90_perm = tuple(
            self.coord_to_idx[self.reduce_coord(-self.coords[i][1], self.coords[i][0])]
            for i in range(self.N)
        )
        self.conj_perm = tuple(
            self.coord_to_idx[self.reduce_coord(self.coords[i][0], -self.coords[i][1])]
            for i in range(self.N)
        )

        # Precompute exact translation-invariant zero-mean Symanzik Green kernel G_0 = K_Sym^+
        (
            self.g0_first_row,
            self.g0_matrix,
            self.g0_orbit_stencils,
            self.g0_int_orbit_stencils,
            self.g0_inv_common_den,
        ) = self._precompute_g0_kernel()
        self.g0_common_den = Fraction((ONE_Q / self.g0_inv_common_den).numerator, 1)

        # Precompute exact (1+i)-adic subgroup coset tables for roughness shells tilde_R = 0..m
        self.coset_tables = self._precompute_adic_cosets()

    def reduce_coord(self, x_val, y_val):
        """Reduce `x_val + i*y_val` modulo `<(1 + i)^m>` to canonical `(rx, ry)`."""
        if self.twist == 0:
            return (x_val % self.nx, y_val % self.ny)
        q_y, ry = divmod(y_val, self.ny)
        rx = (x_val + q_y * self.twist) % self.nx
        return (rx, ry)

    def shift_index(self, idx, dx, dy):
        x_val, y_val = self.coords[idx]
        return self.coord_to_idx[self.reduce_coord(x_val + dx, y_val + dy)]

    def _d4_orbit_data(self):
        """Partition sites into exact D_4 symmetry orbits around origin 0."""
        n_sites = self.N
        id_p = tuple(range(n_sites))
        r1 = self.rot90_perm
        r2 = tuple(r1[r1[i]] for i in range(n_sites))
        r3 = tuple(r1[r2[i]] for i in range(n_sites))
        c0 = self.conj_perm
        c1 = tuple(r1[c0[i]] for i in range(n_sites))
        c2 = tuple(r2[c0[i]] for i in range(n_sites))
        c3 = tuple(r3[c0[i]] for i in range(n_sites))
        d4_perms = (id_p, r1, r2, r3, c0, c1, c2, c3)

        site_to_oid = [-1] * n_sites
        orbit_members = []
        for i_idx in range(n_sites):
            if site_to_oid[i_idx] != -1:
                continue
            oid = len(orbit_members)
            members = sorted({perm[i_idx] for perm in d4_perms})
            for m_idx in members:
                site_to_oid[m_idx] = oid
            orbit_members.append(tuple(members))
        return tuple(site_to_oid), tuple(orbit_members)

    def _precompute_g0_kernel(self):
        """Solve `(K_Sym + (1/N) 1 1^T) g = e_0 - (1/N) 1` exactly over Q via D_4 orbits.

        Factors out the single exact LCM `D_common` of all orbit denominators so `apply_G0`
        accumulates pure integer weights `w_int * orb_sum` and divides ONCE per site.
        """
        n_sites = self.N
        inv_n = Fraction(1, n_sites)
        site_to_oid, orbit_members = self._d4_orbit_data()
        n_orbits = len(orbit_members)

        mat_cols = []
        for col_oid in range(n_orbits):
            basis_orb = [ZERO_Q] * n_sites
            for m_idx in orbit_members[col_oid]:
                basis_orb[m_idx] = ONE_Q
            k_col = self.apply_K_sym(basis_orb)
            orb_card_div_n = Fraction(len(orbit_members[col_oid]), n_sites)
            mat_cols.append(
                [k_col[orbit_members[row_oid][0]] + orb_card_div_n for row_oid in range(n_orbits)]
            )

        mat_rows = [
            [mat_cols[col_oid][row_oid] for col_oid in range(n_orbits)]
            for row_oid in range(n_orbits)
        ]
        rhs = [
            (ONE_Q - inv_n) if (orbit_members[row_oid][0] == 0) else (-inv_n)
            for row_oid in range(n_orbits)
        ]
        sol_orbits = solve_linear_system_q(mat_rows, rhs)
        g0_row0 = tuple(sol_orbits[site_to_oid[i_idx]] for i_idx in range(n_sites))

        d_common = 1
        for q_val in sol_orbits:
            if q_val != ZERO_Q:
                d_common = _lcm_positive_ints_q(d_common, q_val.denominator)
        inv_common_den = Fraction(1, d_common)
        int_orbits = [
            (sol_orbits[oid] * d_common).numerator for oid in range(n_orbits)
        ]

        full_rows = []
        orbit_stencils = []
        int_orbit_stencils = []
        for i_idx in range(n_sites):
            xi, yi = self.coords[i_idx]
            row_vals = []
            by_oid = [[] for oid_i in range(n_orbits)]
            for j_idx in range(n_sites):
                xj, yj = self.coords[j_idx]
                diff_idx = self.coord_to_idx[self.reduce_coord(xi - xj, yi - yj)]
                row_vals.append(g0_row0[diff_idx])
                by_oid[site_to_oid[diff_idx]].append(j_idx)
            full_rows.append(tuple(row_vals))
            orbit_stencils.append(
                tuple(
                    (sol_orbits[oid], tuple(by_oid[oid]))
                    for oid in range(n_orbits)
                    if sol_orbits[oid] != ZERO_Q
                )
            )
            int_orbit_stencils.append(
                tuple(
                    (int_orbits[oid], tuple(by_oid[oid]))
                    for oid in range(n_orbits)
                    if int_orbits[oid] != 0
                )
            )
        return (
            g0_row0,
            tuple(full_rows),
            tuple(orbit_stencils),
            tuple(int_orbit_stencils),
            inv_common_den,
        )

    def _precompute_adic_cosets(self):
        """Build exact disjoint subgroup cosets `z + <(1+i)^r>` for `r = 0..m`."""
        n_sites = self.N
        tables = []
        disjoint_levels = []
        for r_lvl in range(self.m + 1):
            wr_x, wr_y = _gaussian_pow_1_plus_i(r_lvl)
            subgroup_set = set()
            bound = self.nx
            for a_idx in range(bound):
                for b_idx in range(bound):
                    gx = a_idx * wr_x - b_idx * wr_y
                    gy = a_idx * wr_y + b_idx * wr_x
                    subgroup_set.add(self.reduce_coord(gx, gy))
            subgroup_pts = sorted(subgroup_set)
            site_cosets = []
            visited = [False] * n_sites
            coset_groups = []
            site_to_cid = [0] * n_sites
            for i_idx in range(n_sites):
                xi, yi = self.coords[i_idx]
                coset_indices = tuple(
                    self.coord_to_idx[self.reduce_coord(xi + hx, yi + hy)]
                    for hx, hy in subgroup_pts
                )
                site_cosets.append(coset_indices)
                if not visited[i_idx]:
                    cid = len(coset_groups)
                    coset_groups.append(coset_indices)
                    for m_idx in coset_indices:
                        visited[m_idx] = True
                        site_to_cid[m_idx] = cid
            tables.append(tuple(site_cosets))
            disjoint_levels.append(
                (Fraction(1, len(subgroup_pts)), tuple(coset_groups), tuple(site_to_cid))
            )
        self.disjoint_coset_levels = tuple(disjoint_levels)
        return tuple(tables)

    # --- Spatial Operators (100% duck-typed over Fraction, QQPoly, QQRat, QQJet, AuditJet) ---

    def apply_Lx(self, field, range_step=1):
        """Positive semi-definite 1D graph Laplacian `(L_x f)(z) = 2 f(z) - f(z+s) - f(z-s)`."""
        n_sites = self.N
        if range_step == 1:
            p_tab, m_tab = self.xp1, self.xm1
        elif range_step == 2:
            p_tab, m_tab = self.xp2, self.xm2
        else:
            p_tab = tuple(self.shift_index(i, range_step, 0) for i in range(n_sites))
            m_tab = tuple(self.shift_index(i, -range_step, 0) for i in range(n_sites))
        return [
            2 * field[i] - field[p_tab[i]] - field[m_tab[i]]
            for i in range(n_sites)
        ]

    def apply_Ly(self, field, range_step=1):
        """Positive semi-definite 1D graph Laplacian `(L_y f)(z) = 2 f(z) - f(z+i*s) - f(z-i*s)`."""
        n_sites = self.N
        if range_step == 1:
            p_tab, m_tab = self.yp1, self.ym1
        elif range_step == 2:
            p_tab, m_tab = self.yp2, self.ym2
        else:
            p_tab = tuple(self.shift_index(i, 0, range_step) for i in range(n_sites))
            m_tab = tuple(self.shift_index(i, 0, -range_step) for i in range(n_sites))
        return [
            2 * field[i] - field[p_tab[i]] - field[m_tab[i]]
            for i in range(n_sites)
        ]

    def apply_K_sym(self, field):
        """4th-order Symanzik Laplacian `K_Sym = (4/3)(L_x + L_y) - (1/12)(L_{2x} + L_{2y})`."""
        n_sites = self.N
        xp1, xm1, xp2, xm2 = self.xp1, self.xm1, self.xp2, self.xm2
        yp1, ym1, yp2, ym2 = self.yp1, self.ym1, self.yp2, self.ym2
        inv_12 = TWELFTH_Q
        out = []
        for i in range(n_sites):
            fi = field[i]
            two_fi = 2 * fi
            l1 = (two_fi - field[xp1[i]] - field[xm1[i]]) + (two_fi - field[yp1[i]] - field[ym1[i]])
            l2 = (two_fi - field[xp2[i]] - field[xm2[i]]) + (two_fi - field[yp2[i]] - field[ym2[i]])
            out.append(inv_12 * (16 * l1 - l2))
        return out

    def apply_G0_int(self, omega):
        """Integer-weighted `g0_common_den * (G_0 omega)` (pure C PyLong on int, duck-typed on Q/CAS)."""
        n_sites = self.N
        stencils = self.g0_int_orbit_stencils
        out = []
        for i in range(n_sites):
            acc = 0
            for w_int, nbrs in stencils[i]:
                orb_sum = omega[nbrs[0]]
                for k_idx in range(1, len(nbrs)):
                    orb_sum = orb_sum + omega[nbrs[k_idx]]
                acc = acc + w_int * orb_sum
            out.append(acc)
        return out

    def apply_G0(self, omega):
        """Exact zero-mean streamfunction `psi = G_0 omega` (`K_Sym G_0 == I - (1/N) 1 1^T`)."""
        inv_d = self.g0_inv_common_den
        return [inv_d * acc for acc in self.apply_G0_int(omega)]

    def apply_Dx_12(self, field):
        """Integer-scaled `12 * D_x^{Sym}` (`8*(T_{+x}-T_{-x}) - (T_{+2x}-T_{-2x})`, zero divisions)."""
        n_sites = self.N
        xp1, xm1, xp2, xm2 = self.xp1, self.xm1, self.xp2, self.xm2
        return [
            8 * (field[xp1[i]] - field[xm1[i]]) - (field[xp2[i]] - field[xm2[i]])
            for i in range(n_sites)
        ]

    def apply_Dy_12(self, field):
        """Integer-scaled `12 * D_y^{Sym}` (`8*(T_{+y}-T_{-y}) - (T_{+2y}-T_{-2y})`, zero divisions)."""
        n_sites = self.N
        yp1, ym1, yp2, ym2 = self.yp1, self.ym1, self.yp2, self.ym2
        return [
            8 * (field[yp1[i]] - field[ym1[i]]) - (field[yp2[i]] - field[ym2[i]])
            for i in range(n_sites)
        ]

    def apply_Dx_sym(self, field):
        """4th-order skew-adjoint derivative `D_x^{Sym} = (2/3)(T_{+x} - T_{-x}) - (1/12)(T_{+2x} - T_{-2x})`."""
        inv_12 = TWELFTH_Q
        return [inv_12 * val for val in self.apply_Dx_12(field)]

    def apply_Dy_sym(self, field):
        """4th-order skew-adjoint derivative `D_y^{Sym} = (2/3)(T_{+y} - T_{-y}) - (1/12)(T_{+2y} - T_{-2y})`."""
        inv_12 = TWELFTH_Q
        return [inv_12 * val for val in self.apply_Dy_12(field)]

    def velocity_from_psi(self, psi):
        """Incompressible 2D velocity `(u_x, u_y) = (+D_y^{Sym} psi, -D_x^{Sym} psi)`."""
        ux = self.apply_Dy_sym(psi)
        dx_psi = self.apply_Dx_sym(psi)
        uy = [-val for val in dx_psi]
        return [ux, uy]

    def arakawa_bracket_432(self, psi, omega):
        """Division-free integer-scaled `432 * {psi, omega}` (preserves `den == 1` with zero GCDs)."""
        n_sites = self.N
        ax_psi = self.apply_Dx_12(psi)
        by_psi = self.apply_Dy_12(psi)
        ax_om = self.apply_Dx_12(omega)
        by_om = self.apply_Dy_12(omega)

        j1 = [ax_psi[i] * by_om[i] - by_psi[i] * ax_om[i] for i in range(n_sites)]

        psi_by_om = [psi[i] * by_om[i] for i in range(n_sites)]
        psi_ax_om = [psi[i] * ax_om[i] for i in range(n_sites)]
        d_psi_by_om = self.apply_Dx_12(psi_by_om)
        d_psi_ax_om = self.apply_Dy_12(psi_ax_om)
        j2 = [d_psi_by_om[i] - d_psi_ax_om[i] for i in range(n_sites)]

        ax_psi_om = [ax_psi[i] * omega[i] for i in range(n_sites)]
        by_psi_om = [by_psi[i] * omega[i] for i in range(n_sites)]
        d_ax_psi_om = self.apply_Dy_12(ax_psi_om)
        d_by_psi_om = self.apply_Dx_12(by_psi_om)
        j3 = [d_ax_psi_om[i] - d_by_psi_om[i] for i in range(n_sites)]

        return [j1[i] + j2[i] + j3[i] for i in range(n_sites)]

    def arakawa_bracket(self, psi, omega):
        """4th-order Symanzik Arakawa skew-symmetric Jacobian `{psi, omega}`.

        Identically satisfies over `Q` and `QQPoly`:
          - `sum_z {psi, omega}(z) == 0`,
          - `<psi, {psi, omega}> == 0`,
          - `<omega, {psi, omega}> == 0`,
          - `{omega, psi} == -{psi, omega}`.
        """
        inv_432 = Fraction(1, 432)
        return [inv_432 * val for val in self.arakawa_bracket_432(psi, omega)]

    def apply_W_flat(self, field):
        """Sum-of-Squares equal-corner (`16:16:16`) biharmonic `W_flat = (L_x - (1/2)L_y)^2 + (3/4)L_y^2`."""
        n_sites = self.N
        lx_f = self.apply_Lx(field, 1)
        ly_f = self.apply_Ly(field, 1)
        a_f = [lx_f[i] - HALF_Q * ly_f[i] for i in range(n_sites)]

        lx_af = self.apply_Lx(a_f, 1)
        ly_af = self.apply_Ly(a_f, 1)
        a2_f = [lx_af[i] - HALF_Q * ly_af[i] for i in range(n_sites)]

        ly2_f = self.apply_Ly(ly_f, 1)
        return [a2_f[i] + THREE_QUARTERS_Q * ly2_f[i] for i in range(n_sites)]

    def project_roughness_shells(self, field):
        """Decompose `field` into `m + 1` orthogonal `(1+i)`-adic roughness shells `tilde_R = 0..m`.

        Sums each disjoint coset once per level (`(m+1)*N` additions instead of `2*N^2`),
        100% duck-typed over `Fraction`, `QQPoly`, `QQRat`, `QQJet`, and `AuditJet`.
        """
        n_sites = self.N
        avgs = []
        for r_lvl in range(self.m + 1):
            inv_card, coset_groups, site_to_cid = self.disjoint_coset_levels[r_lvl]
            group_means = []
            for c_members in coset_groups:
                acc = field[c_members[0]]
                for k_idx in range(1, len(c_members)):
                    acc = acc + field[c_members[k_idx]]
                group_means.append(inv_card * acc)
            avgs.append([group_means[site_to_cid[i]] for i in range(n_sites)])

        shells = [avgs[0]]
        for tilde_r in range(1, self.m + 1):
            fine_lvl = self.m + 1 - tilde_r
            coarse_lvl = self.m - tilde_r
            fine_a = avgs[fine_lvl]
            coarse_a = avgs[coarse_lvl]
            shells.append([fine_a[i] - coarse_a[i] for i in range(n_sites)])
        return shells

    def energy_spectrum_shells(self, omega):
        """Exact rational roughness-shell kinetic energy spectrum `E(tilde_R)` (`sum == E_kin` to `0/1`).

        Uses translation-invariance `[G_0, P_{tilde_R}] == 0` so `G_0` is evaluated ONCE on `omega`
        rather than `m + 1` times (`0/1` identical over `Q` and `QQPoly`, `(m + 1)x` faster).
        """
        om_shells = self.project_roughness_shells(omega)
        psi = self.apply_G0(omega)
        psi_shells = self.project_roughness_shells(psi)
        spec = []
        for r_lvl in range(self.m + 1):
            sh = om_shells[r_lvl]
            psi_sh = psi_shells[r_lvl]
            acc = ZERO_Q
            for i in range(self.N):
                acc = acc + sh[i] * psi_sh[i]
            spec.append(HALF_Q * acc)
        return spec

    def enstrophy_spectrum_shells(self, omega):
        """Exact rational roughness-shell enstrophy spectrum `Z(tilde_R)` (`sum == Z` to `0/1`)."""
        om_shells = self.project_roughness_shells(omega)
        spec = []
        for r_lvl in range(self.m + 1):
            sh = om_shells[r_lvl]
            acc = ZERO_Q
            for i in range(self.N):
                acc = acc + sh[i] * sh[i]
            spec.append(HALF_Q * acc)
        return spec


_GRID_2D_CACHE = {}


def get_grid_2d(m):
    """Cached factory for `GaussianTorusGrid2D(m)`."""
    if m not in _GRID_2D_CACHE:
        _GRID_2D_CACHE[m] = GaussianTorusGrid2D(m)
    return _GRID_2D_CACHE[m]


# --- Module: birational_cores.py ---
ZERO_Q = Fraction(0, 1)
ONE_Q = Fraction(1, 1)
HALF_Q = Fraction(1, 2)
QUARTER_Q = Fraction(1, 4)
EIGHTH_Q = Fraction(1, 8)
TWELFTH_Q = Fraction(1, 12)
NUM_DOUBLE_CAYLEY_PLANES_2D = 2
DEFAULT_DENSITY_SOUND_SPEED = Fraction(6, 1)


class Jet2:
    """Exact truncated 2-jet `c0 + c1*eps + c2*eps^2` with pure left-scalar algebraic dispatch."""

    __slots__ = ("c0", "c1", "c2")

    def __init__(self, c0, c1=ZERO_Q, c2=ZERO_Q):
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2

    def __pos__(self):
        return self

    def __neg__(self):
        return Jet2(-self.c0, -self.c1, -self.c2)

    def __add__(self, other):
        return Jet2(self.c0 + other.c0, self.c1 + other.c1, self.c2 + other.c2)

    def __radd__(self, scalar):
        return Jet2(scalar + self.c0, self.c1, self.c2)

    def __sub__(self, other):
        return Jet2(self.c0 - other.c0, self.c1 - other.c1, self.c2 - other.c2)

    def __rsub__(self, scalar):
        return Jet2(scalar - self.c0, -self.c1, -self.c2)

    def __mul__(self, other):
        a0, a1, a2 = self.c0, self.c1, self.c2
        b0, b1, b2 = other.c0, other.c1, other.c2
        return Jet2(
            a0 * b0,
            a0 * b1 + a1 * b0,
            a0 * b2 + a1 * b1 + a2 * b0,
        )

    def __rmul__(self, scalar):
        return Jet2(scalar * self.c0, scalar * self.c1, scalar * self.c2)

    def inv(self):
        inv0 = ONE_Q / self.c0
        u1 = self.c1 * inv0
        u2 = self.c2 * inv0
        return Jet2(inv0, -u1 * inv0, (u1 * u1 - u2) * inv0)


def _dot_1d_raw(u_vec, v_vec):
    """Raw inner product initialized with `0` for C `PyLong` speed + CAS duck-typing."""
    acc = 0
    for i in range(len(u_vec)):
        acc = acc + u_vec[i] * v_vec[i]
    return acc


def _dot_vec3_raw(u_vec3, v_vec3):
    """Raw 3D vector inner product initialized with `0` for C `PyLong` speed + CAS duck-typing."""
    n_sites = len(u_vec3[0])
    acc = 0
    for comp in range(3):
        uc = u_vec3[comp]
        vc = v_vec3[comp]
        for i in range(n_sites):
            acc = acc + uc[i] * vc[i]
    return acc


def dot_1d(u_vec, v_vec):
    """Exact inner product `sum_i u_vec[i] * v_vec[i]` over `Q` or CAS ring."""
    return ZERO_Q + _dot_1d_raw(u_vec, v_vec)


def dot_vec3(u_vec3, v_vec3):
    """Exact 3D vector field inner product `sum_{c=0..2} sum_x u_c(x) * v_c(x)`."""
    return ZERO_Q + _dot_vec3_raw(u_vec3, v_vec3)


def clear_vec_denominators_q(vec):
    """Scale `vec` by positive integer LCM of `.denominator` so `Fraction` entries have `den == 1`."""
    d_lcm = 1
    for val in vec:
        d_lcm = _lcm_positive_ints_q(d_lcm, val.denominator)
    return [val * d_lcm for val in vec], Fraction(d_lcm, 1)


def solve_3x3_cramer(g_mat, rhs):
    """Exact closed-form 3x3 Cramer solver (branch-free, 100% duck-typed over Q / QQRat / Jet2)."""
    a00, a01, a02 = g_mat[0][0], g_mat[0][1], g_mat[0][2]
    a10, a11, a12 = g_mat[1][0], g_mat[1][1], g_mat[1][2]
    a20, a21, a22 = g_mat[2][0], g_mat[2][1], g_mat[2][2]
    b0, b1, b2 = rhs[0], rhs[1], rhs[2]

    m00 = a11 * a22 - a12 * a21
    m01 = a10 * a22 - a12 * a20
    m02 = a10 * a21 - a11 * a20

    det = a00 * m00 - a01 * m01 + a02 * m02
    inv_det = ONE_Q / det

    x0 = (b0 * m00 - a01 * (b1 * a22 - a12 * b2) + a02 * (b1 * a21 - a11 * b2)) * inv_det
    x1 = (a00 * (b1 * a22 - a12 * b2) - b0 * m01 + a02 * (a10 * b2 - b1 * a20)) * inv_det
    x2 = (a00 * (a11 * b2 - b1 * a21) - a01 * (a10 * b2 - b1 * a20) + b0 * m02) * inv_det
    return (x0, x1, x2)


def invert_4x4_q(mat4):
    """Exact Gauss-Jordan inverse of a 4x4 matrix over Q or QQRat (fixed 4x4 loop, zero `while`)."""
    aug = [
        [mat4[r][c] for c in range(4)] + [ONE_Q if r == c else ZERO_Q for c in range(4)]
        for r in range(4)
    ]
    for col in range(4):
        piv = col
        if aug[piv][col] == ZERO_Q:
            for r_try in range(col + 1, 4):
                if aug[r_try][col] != ZERO_Q:
                    piv = r_try
                    break
        if piv != col:
            aug[col], aug[piv] = aug[piv], aug[col]
        inv_p = ONE_Q / aug[col][col]
        aug[col] = [val * inv_p for val in aug[col]]
        for r in range(4):
            if r != col:
                fac = aug[r][col]
                if fac != ZERO_Q:
                    aug[r] = [aug[r][k] - fac * aug[col][k] for k in range(8)]
    return [row[4:8] for row in aug]


# =============================================================================
# 1. Galilean-Covariant 1-Quadric Stereographic Readouts & Pullbacks
# =============================================================================


def _split_mean_and_fluct_3d(n_sites, alpha_vec3):
    """Split `alpha_vec3 in V_div` into spatial mean `u_mean[3]` and zero-mean `alpha_0[3][N]`."""
    inv_n = Fraction(1, n_sites)
    u_mean = [sum(alpha_vec3[c], ZERO_Q) * inv_n for c in range(3)]
    alpha_0 = [
        [alpha_vec3[c][i] - u_mean[c] for i in range(n_sites)]
        for c in range(3)
    ]
    return u_mean, alpha_0


def build_stereographic_anchor_3d(grid, alpha_vec3, q_bundle=None):
    """Precompute Galilean anchor data parameterized by `q_bundle` (C `PyLong` & `QQRat` unified)."""
    qb = resolve_quantizer_bundle(True, bundle=q_bundle)
    n_sites = grid.N
    u_mean, alpha_0 = _split_mean_and_fluct_3d(n_sites, alpha_vec3)
    a_sq = dot_vec3(alpha_0, alpha_0)
    combined = list(u_mean) + list(alpha_0[0]) + list(alpha_0[1]) + list(alpha_0[2])
    comb_c, d_a_q = qb["clear_vec"](combined)
    d_a_raw = qb["to_raw_scale"](d_a_q)
    u_mean_int = comb_c[0:3]
    alpha_0_int = [
        comb_c[3 : 3 + n_sites],
        comb_c[3 + n_sites : 3 + 2 * n_sites],
        comb_c[3 + 2 * n_sites : 3 + 3 * n_sites],
    ]
    a_sq_int = _dot_vec3_raw(alpha_0_int, alpha_0_int)
    return {
        "u_mean": u_mean,
        "alpha_0": alpha_0,
        "a_sq": a_sq,
        "inv_a": ONE_Q / a_sq,
        "d_a": d_a_q,
        "d_a_raw": d_a_raw,
        "u_mean_int": u_mean_int,
        "alpha_0_int": alpha_0_int,
        "a_sq_int": a_sq_int,
    }


def stereographic_readout_3d(grid, alpha_vec3, eta_raw_vec3, amp_scale=ONE_Q):
    """Galilean-covariant one-shot birational stereographic readout onto `V_div`.

    Splits `alpha_vec3 = u_mean + alpha_0` (`u_mean = P_0 / N`, `alpha_0 in V_{div,0}`):
      `u = u_mean + amp_scale * ((4A - B)*alpha_0 + 4A*eta) / (4A + B)`,
    where `A = <alpha_0, alpha_0>`, `eta in V_{div,0} cap alpha_0^perp`, `B = <eta, eta>`.
    Identically satisfies:
      - `div_Sym(u)(x) == 0/1` sitewise,
      - `sum_x u(x) == sum_x alpha_vec3(x) == P_0` (`0/1` for arbitrary nonzero `P_0`!),
      - `1/2 <u, u> == (N/2)||u_mean||^2 + amp_scale^2 * (1/2 <alpha_0, alpha_0>)`.
    """
    n_sites = grid.N
    u_mean, alpha_0 = _split_mean_and_fluct_3d(n_sites, alpha_vec3)
    eta_div = grid.project_leray_zero_mean(eta_raw_vec3)
    a_sq = dot_vec3(alpha_0, alpha_0)
    inv_a = ONE_Q / a_sq
    dot_ae = dot_vec3(alpha_0, eta_div) * inv_a
    eta = [
        [eta_div[c][i] - dot_ae * alpha_0[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    b_sq = dot_vec3(eta, eta)
    four_a = 4 * a_sq
    denom = four_a + b_sq
    inv_den = ONE_Q / denom
    c_alpha = amp_scale * ((four_a - b_sq) * inv_den)
    c_eta = amp_scale * (four_a * inv_den)
    return [
        [u_mean[c] + c_alpha * alpha_0[c][i] + c_eta * eta[c][i] for i in range(n_sites)]
        for c in range(3)
    ]


def stereographic_readout_anchored_3d(grid, anc3, eta_raw_vec3, amp_scale=ONE_Q, q_bundle=None):
    """Exact `0/1` anchored stereographic readout executing in C `PyLong` on concrete `Fraction`."""
    qb = resolve_quantizer_bundle(True, bundle=q_bundle)
    n_sites = grid.N
    um_int = anc3["u_mean_int"]
    a0_int = anc3["alpha_0_int"]
    a_sq_int = anc3["a_sq_int"]
    d_a = anc3["d_a_raw"]

    eta_i, d_eta_q = qb["clear_vec3"](eta_raw_vec3)
    d_eta = qb["to_raw_scale"](d_eta_q)
    ed_int, s_ed_q = grid.project_leray_zero_mean_int(eta_i)
    s_ed = qb["to_raw_scale"](s_ed_q)

    dot_ae_s = _dot_vec3_raw(a0_int, ed_int)
    eta_s = [
        [a_sq_int * ed_int[c][i] - dot_ae_s * a0_int[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    k_raw = a_sq_int * d_eta * s_ed
    da_sq = d_a * d_a
    b_s = da_sq * _dot_vec3_raw(eta_s, eta_s)
    fa_s = 4 * a_sq_int * (k_raw * k_raw)
    d_s = fa_s + b_s
    ca_s = fa_s - b_s
    ce_s = 4 * a_sq_int * d_a * k_raw

    amp_list, d_amp_q = qb["clear_vec"]([amp_scale])
    amp_i = amp_list[0]
    d_amp = qb["to_raw_scale"](d_amp_q)
    um_coeff = d_amp * d_s
    inv_tot = ONE_Q / (d_amp * d_a * d_s)
    return [
        [
            (um_coeff * um_int[c] + amp_i * (ca_s * a0_int[c][i] + ce_s * eta_s[c][i])) * inv_tot
            for i in range(n_sites)
        ]
        for c in range(3)
    ]


def stereographic_pullback_3d(grid, alpha_vec3, eta_raw_vec3, x_vec3):
    """Exact differential pullback `d(eta)_{u_fluct} [X_tan]` of `stereographic_readout_3d` (`amp_scale = 1`)."""
    n_sites = grid.N
    _u_mean, alpha_0 = _split_mean_and_fluct_3d(n_sites, alpha_vec3)
    eta_div = grid.project_leray_zero_mean(eta_raw_vec3)
    a_sq = dot_vec3(alpha_0, alpha_0)
    inv_a = ONE_Q / a_sq
    dot_ae = dot_vec3(alpha_0, eta_div) * inv_a
    eta = [
        [eta_div[c][i] - dot_ae * alpha_0[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    b_sq = dot_vec3(eta, eta)
    four_a = 4 * a_sq
    inv_den = ONE_Q / (four_a + b_sq)
    u_fluct = [
        [((four_a - b_sq) * inv_den) * alpha_0[c][i] + (four_a * inv_den) * eta[c][i] for i in range(n_sites)]
        for c in range(3)
    ]

    x_div = grid.project_leray_zero_mean(x_vec3)
    dot_ux = dot_vec3(u_fluct, x_div) * inv_a
    x_tan = [
        [x_div[c][i] - dot_ux * u_fluct[c][i] for i in range(n_sites)]
        for c in range(3)
    ]

    scale_main = (four_a + b_sq) * (QUARTER_Q * inv_a)
    dot_ax = dot_vec3(alpha_0, x_tan)
    scale_corr = (HALF_Q * inv_a) * scale_main * dot_ax
    return [
        [
            scale_main * x_tan[c][i] - scale_corr * (2 * alpha_0[c][i] + eta[c][i])
            for i in range(n_sites)
        ]
        for c in range(3)
    ]


def stereographic_inviscid_pullback_and_metrics_3d(grid, anc3, eta_raw_vec3, q_bundle=None):
    """Unified branch-free factored inviscid 3D Euler stereographic pullback.

    Executes 100% in C `PyLong` integer arrays on concrete `DYADIC_QUANTIZER_BUNDLE` runs
    (`< 0.12 s/step` on `L=8` prime-perturbed Galilean-boosted ICs!) while evaluating the
    exact `0/1` rational identity on `EXACT_IDENTITY_BUNDLE` (`QQRat` / `QQJet` / `AuditJet`).
    """
    qb = resolve_quantizer_bundle(True, bundle=q_bundle)
    n_sites = grid.N
    um_int = anc3["u_mean_int"]
    a0_int = anc3["alpha_0_int"]
    a_sq_int = anc3["a_sq_int"]
    d_a = anc3["d_a_raw"]

    eta_i, d_eta_q = qb["clear_vec3"](eta_raw_vec3)
    d_eta = qb["to_raw_scale"](d_eta_q)
    ed_int, s_ed_q = grid.project_leray_zero_mean_int(eta_i)
    s_ed = qb["to_raw_scale"](s_ed_q)

    dot_ae_s = _dot_vec3_raw(a0_int, ed_int)
    eta_s = [
        [a_sq_int * ed_int[c][i] - dot_ae_s * a0_int[c][i] for i in range(n_sites)]
        for c in range(3)
    ]

    k_raw = a_sq_int * d_eta * s_ed
    da_sq = d_a * d_a
    b_s = da_sq * _dot_vec3_raw(eta_s, eta_s)
    fa_s = 4 * a_sq_int * (k_raw * k_raw)
    d_s = fa_s + b_s
    ca_s = fa_s - b_s
    ce_s = 4 * a_sq_int * d_a * k_raw

    uf_num = [
        [ca_s * a0_int[c][i] + ce_s * eta_s[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    u_num = [
        [d_s * um_int[c] + uf_num[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    u_scale_q = ONE_Q * (d_a * d_s)
    u_snap, d_u_q = qb["snap_vec3"](u_num, u_scale_q)

    adv_24h = grid.skew_advection_24h(u_snap)
    x_div_int, s_xdiv_q = grid.project_leray_zero_mean_int(adv_24h)
    s_xdiv = qb["to_raw_scale"](s_xdiv_q)
    x_tot_scale = 24 * grid.h * (d_u_q * d_u_q * s_xdiv)

    uf_den2 = d_s * d_s * a_sq_int
    dot_ux_n = _dot_vec3_raw(uf_num, x_div_int)
    x_tan_num = [
        [uf_den2 * x_div_int[c][i] - dot_ux_n * uf_num[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    dot_ax_n = _dot_vec3_raw(a0_int, x_tan_num)
    two_ak = 2 * a_sq_int * k_raw
    deta_num = [
        [
            two_ak * x_tan_num[c][i] - dot_ax_n * (2 * k_raw * a0_int[c][i] + d_a * eta_s[c][i])
            for i in range(n_sites)
        ]
        for c in range(3)
    ]
    deta_scale_q = ONE_Q * (two_ak * fa_s * d_s * a_sq_int) * x_tot_scale
    dr_snap, d_dr_q = qb["snap_vec3"](deta_num, deta_scale_q)
    return deta_num, deta_scale_q, dr_snap, d_dr_q, x_div_int, x_tot_scale


def fast_stereographic_inviscid_pullback_3d(grid, anc3, eta_raw_vec3, q_bundle=None):
    """Exact inviscid Euler stereographic pullback wrapper returning `deta_raw`."""
    deta_num, deta_scale_q, _dr_s, _ddr, _xi, _xs = stereographic_inviscid_pullback_and_metrics_3d(
        grid, anc3, eta_raw_vec3, q_bundle=q_bundle
    )
    inv_scale = ONE_Q / deta_scale_q
    n_sites = grid.N
    return [
        [deta_num[c][i] * inv_scale for i in range(n_sites)]
        for c in range(3)
    ]


def stereographic_chart_jet2_accel_3d(alpha_vec3, eta_vec3, deta_1_vec3):
    """Exact 2nd-order chart curvature `(1/2) d^2 u_fluct(alpha_0, eta)[deta_1, deta_1]` in O(N) ops."""
    n_sites = len(alpha_vec3[0])
    _u_mean, alpha_0 = _split_mean_and_fluct_3d(n_sites, alpha_vec3)
    a_sq = dot_vec3(alpha_0, alpha_0)
    four_a = 4 * a_sq

    b0 = dot_vec3(eta_vec3, eta_vec3)
    b1 = 2 * dot_vec3(eta_vec3, deta_1_vec3)
    b2 = dot_vec3(deta_1_vec3, deta_1_vec3)
    b_jet = Jet2(b0, b1, b2)
    den_jet = four_a + b_jet
    inv_den = den_jet.inv()
    ca_jet = (four_a - b_jet) * inv_den
    ce_jet = four_a * inv_den

    ca2 = ca_jet.c2
    ce1 = ce_jet.c1
    ce2 = ce_jet.c2
    return [
        [
            ca2 * alpha_0[c][i] + ce2 * eta_vec3[c][i] + ce1 * deta_1_vec3[c][i]
            for i in range(n_sites)
        ]
        for c in range(3)
    ]


def stereographic_readout_2d_energy(grid, alpha, eta_raw, amp_scale=ONE_Q):
    """One-shot birational stereographic readout onto 2D energy sphere `<om, G_0 om> == 2 E_0`."""
    n_sites = grid.N
    inv_n = Fraction(1, n_sites)
    mean_e = sum(eta_raw, ZERO_Q) * inv_n
    eta_z = [eta_raw[i] - mean_e for i in range(n_sites)]

    g0_alpha = grid.apply_G0(alpha)
    a_sq = dot_1d(alpha, g0_alpha)
    inv_a = ONE_Q / a_sq
    dot_ae = dot_1d(g0_alpha, eta_z) * inv_a
    eta = [eta_z[i] - dot_ae * alpha[i] for i in range(n_sites)]

    g0_eta = grid.apply_G0(eta)
    b_sq = dot_1d(eta, g0_eta)
    four_a = 4 * a_sq
    inv_den = ONE_Q / (four_a + b_sq)
    c_alpha = amp_scale * ((four_a - b_sq) * inv_den)
    c_eta = amp_scale * (four_a * inv_den)
    return [c_alpha * alpha[i] + c_eta * eta[i] for i in range(n_sites)]


def stereographic_pullback_2d_energy(grid, alpha, eta_raw, x_vec):
    """Exact differential pullback `d(eta)_om [X_tan]` for `<., G_0 .>` sphere (`amp_scale = 1`)."""
    n_sites = grid.N
    inv_n = Fraction(1, n_sites)
    mean_e = sum(eta_raw, ZERO_Q) * inv_n
    eta_z = [eta_raw[i] - mean_e for i in range(n_sites)]

    g0_alpha = grid.apply_G0(alpha)
    a_sq = dot_1d(alpha, g0_alpha)
    inv_a = ONE_Q / a_sq
    dot_ae = dot_1d(g0_alpha, eta_z) * inv_a
    eta = [eta_z[i] - dot_ae * alpha[i] for i in range(n_sites)]

    g0_eta = grid.apply_G0(eta)
    b_sq = dot_1d(eta, g0_eta)
    four_a = 4 * a_sq
    inv_den = ONE_Q / (four_a + b_sq)
    c_a0 = (four_a - b_sq) * inv_den
    c_e0 = four_a * inv_den
    om_cur = [c_a0 * alpha[i] + c_e0 * eta[i] for i in range(n_sites)]
    g0_om = [c_a0 * g0_alpha[i] + c_e0 * g0_eta[i] for i in range(n_sites)]

    mean_x = sum(x_vec, ZERO_Q) * inv_n
    x_z = [x_vec[i] - mean_x for i in range(n_sites)]
    dot_ox = dot_1d(g0_om, x_z) * inv_a
    x_tan = [x_z[i] - dot_ox * om_cur[i] for i in range(n_sites)]

    scale_main = (four_a + b_sq) * (QUARTER_Q * inv_a)
    dot_ax = dot_1d(g0_alpha, x_tan)
    scale_corr = (HALF_Q * inv_a) * scale_main * dot_ax
    return [
        scale_main * x_tan[i] - scale_corr * (2 * alpha[i] + eta[i])
        for i in range(n_sites)
    ]


def stereographic_chart_jet2_accel_2d_energy(grid, alpha, eta_raw, deta_1):
    """Exact 2nd-order chart curvature `(1/2) d^2 om(alpha, eta)[deta_1, deta_1]` for 2D energy sphere."""
    n_sites = grid.N
    inv_n = Fraction(1, n_sites)
    mean_e = sum(eta_raw, ZERO_Q) * inv_n
    eta_z = [eta_raw[i] - mean_e for i in range(n_sites)]
    g0_alpha = grid.apply_G0(alpha)
    a_sq = dot_1d(alpha, g0_alpha)
    inv_a = ONE_Q / a_sq
    dot_ae = dot_1d(g0_alpha, eta_z) * inv_a
    eta = [eta_z[i] - dot_ae * alpha[i] for i in range(n_sites)]

    g0_eta = grid.apply_G0(eta)
    g0_deta1 = grid.apply_G0(deta_1)

    four_a = 4 * a_sq
    b0 = dot_1d(eta, g0_eta)
    b1 = 2 * dot_1d(g0_eta, deta_1)
    b2 = dot_1d(deta_1, g0_deta1)
    b_jet = Jet2(b0, b1, b2)
    inv_den = (four_a + b_jet).inv()
    ca_jet = (four_a - b_jet) * inv_den
    ce_jet = four_a * inv_den

    ca2 = ca_jet.c2
    ce1 = ce_jet.c1
    ce2 = ce_jet.c2
    return [
        ca2 * alpha[i] + ce2 * eta[i] + ce1 * deta_1[i]
        for i in range(n_sites)
    ]


# =============================================================================
# 2. 2-Quadric Sparse-Pivot Bezout-Vieta + Disjoint-Parity Double-Cayley Core (1A)
# =============================================================================


_CAYLEY_PLANES_CACHE = {}


def get_2d_commuting_cayley_planes(grid):
    """2 disjoint-parity `D_4`-symmetric mutually orthogonal 2D integer planes in `K_Sym = 14/3`.

    - Plane 0 `(p1, q1)` is supported exclusively on even-parity sites `(x + y) % 2 == 0`.
    - Plane 1 `(p2, q2)` is supported exclusively on odd-parity sites `(x + y) % 2 == 1`.
    At every site `(x, y)`, exactly ONE of `{p1, q1, p2, q2}` is `+/- 1` and the other three are `0`!
    """
    m_key = grid.m
    if m_key not in _CAYLEY_PLANES_CACHE:
        c4 = (1, 0, -1, 0)
        s4 = (0, 1, 0, -1)
        n_sites = grid.N
        p1 = tuple(c4[x % 4] * c4[y % 4] for (x, y) in grid.coords)
        q1 = tuple(s4[x % 4] * s4[y % 4] for (x, y) in grid.coords)
        p2 = tuple(c4[x % 4] * s4[y % 4] for (x, y) in grid.coords)
        q2 = tuple(s4[x % 4] * c4[y % 4] for (x, y) in grid.coords)
        _CAYLEY_PLANES_CACHE[m_key] = (
            (p1, q1, Fraction(4, n_sites)),
            (p2, q2, Fraction(4, n_sites)),
        )
    return _CAYLEY_PLANES_CACHE[m_key]


def _apply_double_cayley_2d_parity_split(grid, rot_vec, field_num, quantize=True, q_bundle=None):
    """Return `(out_split, even_mask, d0_n, d1_n)` avoiding cross-parity `d0 * d1` multiplication."""
    qb = resolve_quantizer_bundle(quantize, bundle=q_bundle)
    planes = get_2d_commuting_cayley_planes(grid)
    n_sites = grid.N

    p0, q0, _inv_n0 = planes[0]
    p1, q1, _inv_n1 = planes[1]
    sig_raw = [qb["wrap_mod4"](rot_vec[0]), qb["wrap_mod4"](rot_vec[1])]
    sig_c, d_sig_q = qb["clear_vec"](sig_raw)
    d_sig = qb["to_raw_scale"](d_sig_q)
    s0_i = sig_c[0]
    s1_i = sig_c[1]
    ds_sq = d_sig * d_sig
    s0_sq = s0_i * s0_i
    s1_sq = s1_i * s1_i

    d0 = (4 * ds_sq + s0_sq) * (4 * ds_sq + s0_sq)
    c0_m1 = -32 * ds_sq * s0_sq
    s0_v = 8 * d_sig * s0_i * (4 * ds_sq - s0_sq)
    dp0_4 = 4 * _dot_1d_raw(p0, field_num)
    dq0_4 = 4 * _dot_1d_raw(q0, field_num)
    cp0_n = c0_m1 * dp0_4 + s0_v * dq0_4
    cq0_n = c0_m1 * dq0_4 - s0_v * dp0_4

    d1 = (4 * ds_sq + s1_sq) * (4 * ds_sq + s1_sq)
    c1_m1 = -32 * ds_sq * s1_sq
    s1_v = 8 * d_sig * s1_i * (4 * ds_sq - s1_sq)
    dp1_4 = 4 * _dot_1d_raw(p1, field_num)
    dq1_4 = 4 * _dot_1d_raw(q1, field_num)
    cp1_n = c1_m1 * dp1_4 + s1_v * dq1_4
    cq1_n = c1_m1 * dq1_4 - s1_v * dp1_4

    d0_n = n_sites * d0
    d1_n = n_sites * d1
    out_split = []
    even_mask = []
    for i_idx in range(n_sites):
        if p0[i_idx] != 0:
            out_split.append(d0_n * field_num[i_idx] + p0[i_idx] * cp0_n)
            even_mask.append(True)
        elif q0[i_idx] != 0:
            out_split.append(d0_n * field_num[i_idx] + q0[i_idx] * cq0_n)
            even_mask.append(True)
        elif p1[i_idx] != 0:
            out_split.append(d1_n * field_num[i_idx] + p1[i_idx] * cp1_n)
            even_mask.append(False)
        else:
            out_split.append(d1_n * field_num[i_idx] + q1[i_idx] * cq1_n)
            even_mask.append(False)
    return out_split, even_mask, d0, d1


def apply_double_cayley_2d_factored(grid, rot_vec, field_num, quantize=True, q_bundle=None):
    """Unified branch-free factored Double-Cayley rotation `(out_num, dc_den)`."""
    out_split, even_mask, d0, d1 = _apply_double_cayley_2d_parity_split(
        grid, rot_vec, field_num, quantize=quantize, q_bundle=q_bundle
    )
    n_sites = grid.N
    dc_den = ONE_Q * (n_sites * d0 * d1)
    out_num = [
        (d1 * out_split[i_idx]) if even_mask[i_idx] else (d0 * out_split[i_idx])
        for i_idx in range(n_sites)
    ]
    return out_num, dc_den


def apply_double_cayley_2d(grid, rot_vec, field, quantize=True, q_bundle=None):
    """Apply pole-free `C^1`-periodic Double-Cayley rotation on `sigma_k = rot_k mod 4 in [-2, 2)`."""
    out_split, even_mask, d0, d1 = _apply_double_cayley_2d_parity_split(
        grid, rot_vec, field, quantize=quantize, q_bundle=q_bundle
    )
    n_sites = grid.N
    inv_d0 = ONE_Q / (n_sites * d0)
    inv_d1 = ONE_Q / (n_sites * d1)
    return [
        (out_split[i_idx] * inv_d0) if even_mask[i_idx] else (out_split[i_idx] * inv_d1)
        for i_idx in range(n_sites)
    ]


def apply_4plane_cayley_2d(grid, s_vec, field):
    """Legacy single-Cayley helper retained for symbolic compatibility."""
    planes = get_2d_commuting_cayley_planes(grid)
    n_sites = grid.N
    out = list(field)
    for k_idx in range(len(s_vec)):
        pk, qk, inv_nsq = planes[k_idx % len(planes)]
        sk = s_vec[k_idx]
        dp = inv_nsq * dot_1d(pk, field)
        dq = inv_nsq * dot_1d(qk, field)
        denom = ONE_Q + QUARTER_Q * (sk * sk)
        inv_d = ONE_Q / denom
        c_sin = sk * inv_d
        c_cos_m1 = -(HALF_Q * (sk * sk)) * inv_d
        coef_p = c_cos_m1 * dp + c_sin * dq
        coef_q = c_cos_m1 * dq - c_sin * dp
        for i_idx in range(n_sites):
            out[i_idx] = out[i_idx] + pk[i_idx] * coef_p + qk[i_idx] * coef_q
    return out


def build_bezout_vieta_anchor_2d(grid, alpha):
    """Build exact sparse-4-pivot Bezout-Vieta chart constants for `alpha` and `beta = T_{(1,0)} alpha`."""
    n_sites = grid.N
    beta = [alpha[grid.xp1[i]] for i in range(n_sites)]
    g0_alpha = grid.apply_G0(alpha)
    g0_beta = [g0_alpha[grid.xp1[i]] for i in range(n_sites)]
    delta = [alpha[i] - beta[i] for i in range(n_sites)]
    g0_delta = [g0_alpha[i] - g0_beta[i] for i in range(n_sites)]
    a1 = dot_1d(delta, delta)
    a2 = dot_1d(delta, g0_delta)
    ell_neg = [2 * (a1 * g0_delta[i] - a2 * delta[i]) for i in range(n_sites)]
    e1 = [ONE_Q for i in range(n_sites)]

    cand_coords = ((0, 0), (1, 0), (0, 1), (2, 1))
    pivs = tuple(grid.coord_to_idx[grid.reduce_coord(cx, cy)] for cx, cy in cand_coords)
    piv_set = set(pivs)
    free_idx = tuple(i for i in range(n_sites) if i not in piv_set)

    rows4 = (e1, beta, g0_beta, ell_neg)
    mat4 = [[rows4[r][pivs[c]] for c in range(4)] for r in range(4)]
    inv_mat4 = invert_4x4_q(mat4)
    coup = [
        [
            -sum(inv_mat4[r][k] * rows4[k][j] for k in range(4))
            for j in free_idx
        ]
        for r in range(4)
    ]

    # Sparse 4-site reference vector `r_dir` in `T_beta = span(1, beta, g0_beta)^perp`:
    mat3 = [[rows4[r][pivs[c]] for c in range(3)] for r in range(3)]
    rhs3 = [-rows4[r][pivs[3]] for r in range(3)]
    r0, r1, r2 = solve_3x3_cramer(mat3, rhs3)
    r_dir = [ZERO_Q for i in range(n_sites)]
    r_dir[pivs[0]] = r0
    r_dir[pivs[1]] = r1
    r_dir[pivs[2]] = r2
    r_dir[pivs[3]] = ONE_Q

    l_r = dot_1d(ell_neg, r_dir)
    inv_lr = ONE_Q / l_r
    g0_r = grid.apply_G0(r_dir)
    q_rr = a2 * dot_1d(r_dir, r_dir) - a1 * dot_1d(r_dir, g0_r)
    nu0 = [a1 * g0_beta[i] - a2 * beta[i] for i in range(n_sites)]

    # Precompute 4x4 pivot submatrix of G_0 and anchor inner products for pure-int chart/pullback:
    g0_piv_4x4 = tuple(
        tuple(grid.g0_matrix[pivs[p_i]][pivs[q_i]] for q_i in range(4))
        for p_i in range(4)
    )
    g0_piv_int_4x4 = tuple(
        tuple((grid.g0_matrix[pivs[p_i]][pivs[q_i]] * grid.g0_common_den).numerator for q_i in range(4))
        for p_i in range(4)
    )
    coup_flat = [coup[r_p][k_f] for r_p in range(4) for k_f in range(len(free_idx))]
    r_dir_piv = [r_dir[pivs[r_p]] for r_p in range(4)]

    return {
        "alpha": alpha,
        "beta": beta,
        "g0_beta": g0_beta,
        "delta": delta,
        "a1": a1,
        "inv_a1": ONE_Q / a1,
        "a2": a2,
        "ell_neg": ell_neg,
        "pivs": pivs,
        "free_idx": free_idx,
        "coup": coup,
        "coup_flat": coup_flat,
        "r_dir": r_dir,
        "r_dir_piv": r_dir_piv,
        "l_r": l_r,
        "inv_lr": inv_lr,
        "g0_r": g0_r,
        "q_rr": q_rr,
        "nu0": nu0,
        "g0_piv_4x4": g0_piv_4x4,
        "g0_piv_int_4x4": g0_piv_int_4x4,
        "dot_beta_alpha": dot_1d(beta, alpha),
        "dot_ell_alpha": dot_1d(ell_neg, alpha),
        "dot_ell_delta": dot_1d(ell_neg, delta),
        "dot_g0r_alpha": dot_1d(g0_r, alpha),
        "dot_g0r_delta": dot_1d(g0_r, delta),
        "dot_g0r_r": dot_1d(g0_r, r_dir),
        "dot_r_alpha": sum(r_dir[pivs[p_i]] * alpha[pivs[p_i]] for p_i in range(4)),
        "dot_r_delta": sum(r_dir[pivs[p_i]] * delta[pivs[p_i]] for p_i in range(4)),
        "dot_r_r": sum(r_dir[pivs[p_i]] * r_dir[pivs[p_i]] for p_i in range(4)),
    }


def _bv_chart_factored_2d(grid, anc, u_free, q_bundle=None):
    """Return `(v_num, den_v)` with `bezout_vieta_single_chart_2d == v_num / den_v` identically."""
    qb = resolve_quantizer_bundle(False, bundle=q_bundle)
    n_sites = grid.N
    pivs = anc["pivs"]
    free_idx = anc["free_idx"]
    g0_piv_int_4x4 = anc["g0_piv_int_4x4"]
    n_free = len(free_idx)

    if len(u_free) == n_sites:
        u_eff = [u_free[j] for j in free_idx]
    else:
        u_eff = u_free

    w_free_vec = [ZERO_Q for i in range(n_sites)]
    for k_f in range(n_free):
        w_free_vec[free_idx[k_f]] = u_eff[k_f]

    wf_c, d_wf_q = qb["clear_vec"](w_free_vec)
    wf_free_c = [wf_c[j] for j in free_idx]

    al_c, d_al_q = qb["clear_vec"](anc["alpha"])
    de_c = [al_c[i] - al_c[grid.xp1[i]] for i in range(n_sites)]
    cp_flat_c, d_cp_q = qb["clear_vec"](anc["coup_flat"])
    d_cp_raw = qb["to_raw_scale"](d_cp_q)
    rd_piv_c, d_rd_q = qb["clear_vec"](anc["r_dir_piv"])
    g0r_c, d_g0r_q = qb["clear_vec"](anc["g0_r"])

    w_piv_cp = []
    for r_p in range(4):
        base_r = r_p * n_free
        w_piv_cp.append(_dot_1d_raw(cp_flat_c[base_r : base_r + n_free], wf_free_c))

    wp_cp = [d_cp_raw * wf_c[i] for i in range(n_sites)]
    for r_p in range(4):
        wp_cp[pivs[r_p]] = w_piv_cp[r_p]

    g0_wf_int = grid.apply_G0_int(wf_c)
    piv_quad_int = 0
    for p_i in range(4):
        row_p = g0_piv_int_4x4[p_i]
        wp_i = w_piv_cp[p_i]
        for q_i in range(4):
            piv_quad_int = piv_quad_int + wp_i * row_p[q_i] * w_piv_cp[q_i]
    dot_wp_g0wp_int = (
        (d_cp_raw * d_cp_raw) * _dot_1d_raw(wf_c, g0_wf_int)
        + (2 * d_cp_raw) * sum(w_piv_cp[p_i] * g0_wf_int[pivs[p_i]] for p_i in range(4))
        + piv_quad_int
    )
    dot_wp_wp_int = _dot_1d_raw(wp_cp, wp_cp)
    dr_wp_int = sum(rd_piv_c[p_i] * w_piv_cp[p_i] for p_i in range(4))
    dg0r_wp_int = _dot_1d_raw(g0r_c, wp_cp)

    a1, a2 = anc["a1"], anc["a2"]
    r_dir, inv_lr, q_rr = anc["r_dir"], anc["inv_lr"], anc["q_rr"]
    inv_dcp = ONE_Q / d_cp_q
    inv_dcp_sq = inv_dcp * inv_dcp
    inv_g0 = ONE_Q / grid.g0_common_den

    q_ww_c = a2 * (dot_wp_wp_int * inv_dcp_sq) - a1 * (dot_wp_g0wp_int * inv_dcp_sq * inv_g0)
    q_wr_c = a2 * (dr_wp_int / (d_rd_q * d_cp_q)) - a1 * (dg0r_wp_int / (d_g0r_q * d_cp_q))

    scale_r_c = inv_lr * q_ww_c
    coef_wp_zeta = d_wf_q * inv_dcp

    d_wf_sq = d_wf_q * d_wf_q
    phi_c = d_wf_sq + (2 * inv_lr * d_wf_q) * q_wr_c + (q_rr * inv_lr) * scale_r_c
    dot_del_zeta = coef_wp_zeta * (_dot_1d_raw(de_c, wp_cp) / d_al_q) + scale_r_c * anc["dot_r_delta"]
    dot_zeta_zeta = (
        (coef_wp_zeta * coef_wp_zeta) * dot_wp_wp_int
        + (2 * coef_wp_zeta * scale_r_c) * (dr_wp_int / d_rd_q)
        + (scale_r_c * scale_r_c) * anc["dot_r_r"]
    )
    b1_c = 2 * dot_del_zeta
    c1_c = dot_zeta_zeta

    den_v = a1 * (phi_c * phi_c) + b1_c * phi_c + c1_c
    c_del_num = -(b1_c * phi_c + c1_c)
    c_zet_num = a1 * phi_c

    inv_dal = ONE_Q / d_al_q
    raw_weights = [
        den_v * inv_dal,
        c_del_num * inv_dal,
        c_zet_num * coef_wp_zeta,
    ] + [c_zet_num * scale_r_c * r_dir[pivs[p_i]] for p_i in range(4)]
    w_c, d_wc_q = qb["clear_vec"](raw_weights)
    k_al, k_de, k_wp = w_c[0], w_c[1], w_c[2]

    v_num = [
        k_al * al_c[i] + k_de * de_c[i] + k_wp * wp_cp[i]
        for i in range(n_sites)
    ]
    for p_i in range(4):
        idx_p = pivs[p_i]
        v_num[idx_p] = v_num[idx_p] + w_c[3 + p_i]
    return v_num, den_v * d_wc_q


def bezout_vieta_single_chart_2d(grid, anc, u_free):
    """Desingularized sparse-4-pivot Bezout-Vieta readout onto `(Gamma, E_kin, Z)` intersection."""
    v_num, den_v = _bv_chart_factored_2d(grid, anc, u_free)
    inv_den = ONE_Q / den_v
    return [val * inv_den for val in v_num]


def bezout_vieta_readout_2d(grid, alpha, w_vec, rot_vec, anc=None, quantize=True, q_bundle=None):
    """Full `(N - 3)`-dimensional birational readout `DC_2(rot_vec)(B_BV(alpha; w_vec))`."""
    qb = resolve_quantizer_bundle(quantize, bundle=q_bundle)
    if anc is None:
        anc = build_bezout_vieta_anchor_2d(grid, alpha)
    v_num, den_v = _bv_chart_factored_2d(grid, anc, w_vec, q_bundle=qb)
    om_split, even_mask, d0, d1 = _apply_double_cayley_2d_parity_split(
        grid, rot_vec, v_num, quantize=quantize, q_bundle=qb
    )
    n_sites = grid.N
    inv_even = ONE_Q / (den_v * n_sites * d0)
    inv_odd = ONE_Q / (den_v * n_sites * d1)
    return [
        (om_split[i_idx] * inv_even) if even_mask[i_idx] else (om_split[i_idx] * inv_odd)
        for i_idx in range(n_sites)
    ]


def bezout_vieta_pullback_factored_2d(grid, anc, u_free, rot_vec, quantize=True, q_bundle=None):
    """Factored branch-free differential lift `(du_num, du_den, drot_num, drot_den, y_0)`.

    Executes all `N`-site and 2-plane inner products in cleared integer/jet rings via `q_bundle`,
    eliminating intermediate `Fraction` GCD reductions on both `DYADIC_QUANTIZER_BUNDLE` and
    `EXACT_FACTORED_BUNDLE` while preserving `0/1` symbolic identity on `EXACT_IDENTITY_BUNDLE`.
    """
    qb = resolve_quantizer_bundle(quantize, bundle=q_bundle)
    n_sites = grid.N
    planes = get_2d_commuting_cayley_planes(grid)

    v_raw, den_v0 = _bv_chart_factored_2d(grid, anc, u_free, q_bundle=qb)
    v_int, den_v = qb["snap_vec"](v_raw, den_v0)

    om_int, dc_den = apply_double_cayley_2d_factored(
        grid, rot_vec, v_int, quantize=quantize, q_bundle=qb
    )
    om_scale = den_v * dc_den

    psi_int = grid.apply_G0_int(om_int)
    x_om_432 = [-val for val in grid.arakawa_bracket_432(psi_int, om_int)]

    neg_rot = [-rot_vec[0], -rot_vec[1]]
    x_tan_num, dc_den_inv = apply_double_cayley_2d_factored(
        grid, neg_rot, x_om_432, quantize=quantize, q_bundle=qb
    )
    xt_den = 432 * grid.g0_common_den * (om_scale * om_scale) * dc_den_inv

    inv_den_v = ONE_Q / den_v

    nu0_c, _d_nu0 = qb["clear_vec"](anc["nu0"])
    need_c = _dot_1d_raw(nu0_c, x_tan_num)
    ak_ints = []
    gk_cs = []
    sos_c = 0
    for k_idx in range(2):
        pk, qk, _inv_nsq = planes[k_idx]
        dp_i = _dot_1d_raw(pk, v_int)
        dq_i = _dot_1d_raw(qk, v_int)
        ak_i = [pk[i] * dq_i - qk[i] * dp_i for i in range(n_sites)]
        ak_ints.append(ak_i)
        gk_c = _dot_1d_raw(nu0_c, ak_i)
        gk_cs.append(gk_c)
        sos_c = sos_c + gk_c * gk_c

    sig_w = [qb["wrap_mod4"](rot_vec[0]), qb["wrap_mod4"](rot_vec[1])]
    sig_c, d_sig_q = qb["clear_vec"](sig_w)
    d_sig_raw = qb["to_raw_scale"](d_sig_q)
    ds_sq = d_sig_raw * d_sig_raw
    drot_num = [
        (4 * ds_sq + sig_c[k_idx] * sig_c[k_idx]) * (n_sites * need_c * gk_cs[k_idx]) * den_v
        for k_idx in range(2)
    ]
    drot_den = (32 * ds_sq) * sos_c * xt_den

    ng0 = need_c * gk_cs[0]
    ng1 = need_c * gk_cs[1]
    ak0, ak1 = ak_ints[0], ak_ints[1]
    dv_sos = [
        sos_c * x_tan_num[i] - ng0 * ak0[i] - ng1 * ak1[i]
        for i in range(n_sites)
    ]

    a1, inv_a1, a2 = anc["a1"], anc["inv_a1"], anc["a2"]
    r_dir, inv_lr, q_rr = anc["r_dir"], anc["inv_lr"], anc["q_rr"]
    pivs = anc["pivs"]
    free_idx = anc["free_idx"]

    al_c, d_al_q = qb["clear_vec"](anc["alpha"])
    be_c = [al_c[grid.xp1[i]] for i in range(n_sites)]
    de_c = [al_c[i] - be_c[i] for i in range(n_sites)]
    el_c, d_el_q = qb["clear_vec"](anc["ell_neg"])
    g0r_c, d_g0r_q = qb["clear_vec"](anc["g0_r"])
    inv_dal = ONE_Q / d_al_q
    inv_del = ONE_Q / d_el_q
    inv_dg0r = ONE_Q / d_g0r_q

    dot_beta_vma = _dot_1d_raw(be_c, v_int) * inv_dal - den_v * anc["dot_beta_alpha"]
    dot_beta_dvsos = _dot_1d_raw(be_c, dv_sos) * inv_dal
    xm1_0_n = -2 * inv_a1 * dot_beta_vma
    xm1_1_sos = -2 * inv_a1 * dot_beta_dvsos
    x_0_n = den_v + xm1_0_n
    x_1_sos = xm1_1_sos

    dlu_0_n = _dot_1d_raw(el_c, v_int) * inv_del - den_v * anc["dot_ell_alpha"] - xm1_0_n * anc["dot_ell_delta"]
    dlu_1_sos = _dot_1d_raw(el_c, dv_sos) * inv_del - xm1_1_sos * anc["dot_ell_delta"]

    c_r0 = inv_lr * dlu_0_n
    c_r1 = inv_lr * dlu_1_sos
    dr_v = sum(r_dir[pivs[p_i]] * v_int[pivs[p_i]] for p_i in range(4))
    dr_dv = sum(r_dir[pivs[p_i]] * dv_sos[pivs[p_i]] for p_i in range(4))

    dr_up0 = dr_v - den_v * anc["dot_r_alpha"] - xm1_0_n * anc["dot_r_delta"] - c_r0 * anc["dot_r_r"]
    dr_up1_sos = dr_dv - xm1_1_sos * anc["dot_r_delta"] - c_r1 * anc["dot_r_r"]
    dg0r_up0 = (
        _dot_1d_raw(g0r_c, v_int) * inv_dg0r
        - den_v * anc["dot_g0r_alpha"]
        - xm1_0_n * anc["dot_g0r_delta"]
        - c_r0 * anc["dot_g0r_r"]
    )
    dg0r_up1_sos = (
        _dot_1d_raw(g0r_c, dv_sos) * inv_dg0r
        - xm1_1_sos * anc["dot_g0r_delta"]
        - c_r1 * anc["dot_g0r_r"]
    )

    q_upr_0_n = a2 * dr_up0 - a1 * dg0r_up0
    q_upr_1_sos = a2 * dr_up1_sos - a1 * dg0r_up1_sos

    y_0_n = x_0_n - (2 * inv_lr) * q_upr_0_n - (q_rr * inv_lr * inv_lr) * dlu_0_n
    y_1_sos = x_1_sos - (2 * inv_lr) * q_upr_1_sos - (q_rr * inv_lr * inv_lr) * dlu_1_sos

    y_0 = y_0_n * inv_den_v
    r_y = y_1_sos / y_0_n
    w4_raw = [
        ONE_Q,
        -r_y,
        (r_y * den_v) * inv_dal,
        (r_y * xm1_0_n - xm1_1_sos) * inv_dal,
    ]
    w4_c, d_w4_q = qb["clear_vec"](w4_raw)
    c_dv, c_vi, c_al, c_de = w4_c[0], w4_c[1], w4_c[2], w4_c[3]
    du_num = [
        c_dv * dv_sos[j] + c_vi * v_int[j] + c_al * al_c[j] + c_de * de_c[j]
        for j in free_idx
    ]
    u_fac = d_w4_q * y_0
    r_fac = 32 * ds_sq
    s_common = sos_c * xt_den
    return du_num, u_fac, drot_num, r_fac, s_common, y_0


def bezout_vieta_pullback_with_contour_2d(grid, anc, u_free, rot_vec, quantize=True, q_bundle=None):
    """Unified branch-free differential lift `(du_free, drot, y_0)` of `bezout_vieta_readout_2d`."""
    du_num, u_fac, drot_num, r_fac, s_common, y_0 = bezout_vieta_pullback_factored_2d(
        grid, anc, u_free, rot_vec, quantize=quantize, q_bundle=q_bundle
    )
    inv_du_den = ONE_Q / (u_fac * s_common)
    inv_dr_den = ONE_Q / (r_fac * s_common)
    du_free = [val * inv_du_den for val in du_num]
    drot = [val * inv_dr_den for val in drot_num]
    return du_free, drot, y_0


def bezout_vieta_pullback_2d(grid, anc, u_free, rot_vec, quantize=True, q_bundle=None):
    """Exact analytic differential lift `(du_free, drot)` of `bezout_vieta_readout_2d`."""
    du_free, drot, _y0 = bezout_vieta_pullback_with_contour_2d(
        grid, anc, u_free, rot_vec, quantize=quantize, q_bundle=q_bundle
    )
    return du_free, drot


# =============================================================================
# 3. Un-Dilated Low-Diffusion [2/2] Pade SOS + Varga Positive Density Core (B3)
# =============================================================================


def varga_sos_departure_rates_3d(grid, u_vec3, kappa, c_s=DEFAULT_DENSITY_SOUND_SPEED):
    """Exact branch-free low-diffusion [2/2] Pade SOS outgoing rates `r_plus[3][N], r_minus[3][N]`.

    Uses local non-Heron rational scale `sigma(u_y) = (4*u_y^2 + 9)/12` and `[2/2]` Pade weight
    `c_a(y) = sigma(u_y) * P_{2,2}(u_y^2 / sigma(u_y)^2)` satisfying exact rational square identity:
      `c_a(y)^2 - u_y^2 == ((4*u_y^2 - 9)^5 / (15360 * (80*u_y^8 + 1088*u_y^6 + 4806*u_y^4 + 4860*u_y^2 + 32805)))^2 >= 0/1`,
    so `r_a^pm(y) = kappa + (1/2)*(c_a(y) +/- u_y) >= kappa > 0/1` and `r_a^+ - r_a^- == u_y` identically.
    """
    n_sites = grid.N
    r_plus = []
    r_minus = []
    for axis in range(3):
        ua = u_vec3[axis]
        rp_ax = []
        rm_ax = []
        for y_idx in range(n_sites):
            uy = ua[y_idx]
            uy_sq = uy * uy
            four_u2_p9 = 4 * uy_sq + 9
            sig_y = four_u2_p9 * TWELFTH_Q
            r0_y = (144 * uy_sq) / (four_u2_p9 * four_u2_p9)
            ca_y = sig_y * pade_22(r0_y)
            rp_ax.append(kappa + HALF_Q * (ca_y + uy))
            rm_ax.append(kappa + HALF_Q * (ca_y - uy))
        r_plus.append(rp_ax)
        r_minus.append(rm_ax)
    return r_plus, r_minus


def _varga_offdiag_over_dt2_1d(grid, kappa, dt_sq):
    """Exact division-by-`dt`-free off-diagonal weights `r2[d] / dt^2 >= 0` (`d = 1..L-1`) for Varga."""
    L = grid.L
    dt_4 = dt_sq * dt_sq
    gv_coeff = (kappa * dt_sq) / (ONE_Q + dt_4)
    alpha_v = dt_sq * gv_coeff
    if L == 6:
        _p0, p1_m, p3_m, p4_m = grid.shell_projectors_1d
        d1 = ONE_Q + alpha_v
        d3 = ONE_Q + 3 * alpha_v
        d4 = ONE_Q + 4 * alpha_v
        c1 = -(gv_coeff * (2 + alpha_v)) / (d1 * d1)
        c3 = -(3 * gv_coeff * (2 + 3 * alpha_v)) / (d3 * d3)
        c4 = -(4 * gv_coeff * (2 + 4 * alpha_v)) / (d4 * d4)
        return [
            ZERO_Q
            if d == 0
            else (c1 * p1_m[0][d] + c3 * p3_m[0][d] + c4 * p4_m[0][d])
            for d in range(6)
        ]
    if L == 8:
        a2 = alpha_v * alpha_v
        d2 = ONE_Q + 2 * alpha_v
        d4 = ONE_Q + 4 * alpha_v
        d_odd = ONE_Q + 4 * alpha_v + 2 * a2
        inv_d_odd_sq = ONE_Q / (d_odd * d_odd)
        c2 = -(2 * gv_coeff * (2 + 2 * alpha_v)) / (d2 * d2)
        c4 = -(4 * gv_coeff * (2 + 4 * alpha_v)) / (d4 * d4)
        c_id_m1 = -(gv_coeff * alpha_v * (6 + 16 * alpha_v + 4 * a2)) * inv_d_odd_sq
        c_lap_dt2 = -(2 * gv_coeff + 4 * gv_coeff * alpha_v) * inv_d_odd_sq
        p2_row = (QUARTER_Q, ZERO_Q, -QUARTER_Q, ZERO_Q, QUARTER_Q, ZERO_Q, -QUARTER_Q, ZERO_Q)
        p4_row = (EIGHTH_Q, -EIGHTH_Q, EIGHTH_Q, -EIGHTH_Q, EIGHTH_Q, -EIGHTH_Q, EIGHTH_Q, -EIGHTH_Q)
        p_odd_row = (HALF_Q, ZERO_Q, ZERO_Q, ZERO_Q, -HALF_Q, ZERO_Q, ZERO_Q, ZERO_Q)
        l_odd_row = (ONE_Q, -HALF_Q, ZERO_Q, HALF_Q, -ONE_Q, HALF_Q, ZERO_Q, -HALF_Q)
        return [
            ZERO_Q
            if d == 0
            else (
                c2 * p2_row[d]
                + c4 * p4_row[d]
                + c_id_m1 * p_odd_row[d]
                + c_lap_dt2 * l_odd_row[d]
            )
            for d in range(8)
        ]
    return [ZERO_Q for d in range(L)]


def _conservative_sos_net_rate_3d(grid, rho_in, r_plus, r_minus, inv_tame_vec, qb):
    """Unified conservative SOS net rate parameterized by `qb` (`sum == 0/1`, local 7-point stencil)."""
    n_sites = grid.N
    q_nonneg = qb["nonneg"]
    net = [ZERO_Q for i in range(n_sites)]
    for axis in range(3):
        p1_ax = grid.p1[axis]
        m1_ax = grid.m1[axis]
        rp_ax = r_plus[axis]
        rm_ax = r_minus[axis]
        for y_idx in range(n_sites):
            ry_scaled = inv_tame_vec[y_idx] * rho_in[y_idx]
            fp_q = q_nonneg(rp_ax[y_idx] * ry_scaled)
            fm_q = q_nonneg(rm_ax[y_idx] * ry_scaled)
            net[y_idx] = net[y_idx] - fp_q - fm_q
            net[p1_ax[y_idx]] = net[p1_ax[y_idx]] + fp_q
            net[m1_ax[y_idx]] = net[m1_ax[y_idx]] + fm_q
    return net


def _varga_1d_pencils(cur, pencils, L, w_off, dt_sq, qb):
    """Apply conservative 1D off-diagonal Varga flux pass per `L`-site pencil using lag symmetry."""
    out = list(cur)
    q_nonneg = qb["nonneg"]
    half_L = L // 2
    for p_tup in pencils:
        f_pen = []
        for s_pos in range(L):
            c_val = cur[p_tup[s_pos]]
            for d_lag in range(1, half_L + 1):
                f_pen.append(q_nonneg(w_off[d_lag] * c_val))
        f_c, d_f = qb["clear_vec"](f_pen)
        net_p = [0 for s_pos in range(L)]
        pos = 0
        for s_pos in range(L):
            for d_lag in range(1, half_L + 1):
                fc_val = f_c[pos]
                pos = pos + 1
                dst_p = (s_pos + d_lag) % L
                net_p[s_pos] = net_p[s_pos] - fc_val
                net_p[dst_p] = net_p[dst_p] + fc_val
                if 2 * d_lag != L:
                    dst_m = (s_pos - d_lag) % L
                    net_p[s_pos] = net_p[s_pos] - fc_val
                    net_p[dst_m] = net_p[dst_m] + fc_val
        scale_dt = dt_sq / d_f
        for s_pos in range(L):
            idx_g = p_tup[s_pos]
            out[idx_g] = cur[idx_g] + scale_dt * net_p[s_pos]
    return out


def _apply_unified_varga_3d(grid, rho_in, kappa, dt_sq, qb):
    """Unified conservative 3D Varga smoothing parameterized by `qb` (`sum == M_0`, `> 0/1`).

    When `qb["nonneg"]` is `identity_scalar`, matches `grid.apply_varga_3d(rho_in, alpha_varga)`
    identically to `0/1` over `Q` and `QQRat`; when `quantize_dyadic_floor_nonneg`,
    preserves `0/1` mass conservation, strict sitewise positivity, and bounded bit-lengths.
    """
    L = grid.L
    L2 = L * L
    w_off = _varga_offdiag_over_dt2_1d(grid, kappa, dt_sq)

    pencils_x = [tuple(x_i * L2 + yz_i for x_i in range(L)) for yz_i in range(L2)]
    pencils_y = [
        tuple(x_i * L2 + y_i * L + z_i for y_i in range(L))
        for x_i in range(L)
        for z_i in range(L)
    ]
    pencils_z = [tuple(xy_i * L + z_i for z_i in range(L)) for xy_i in range(L2)]

    cur = _varga_1d_pencils(rho_in, pencils_x, L, w_off, dt_sq, qb)
    cur = _varga_1d_pencils(cur, pencils_y, L, w_off, dt_sq, qb)
    return _varga_1d_pencils(cur, pencils_z, L, w_off, dt_sq, qb)


def _sos_rates_and_tame_3d(grid, u_in, kappa, dt_sq, c_s, qb):
    """Evaluate snapped `[2/2]` Pade SOS departure rates and quadratic Courant taming vector."""
    q_nonneg = qb["nonneg"]
    n_sites = grid.N
    u_snap, d_u = qb["snap_vec3"](u_in, ONE_Q)
    inv_du = ONE_Q / d_u
    u_eff = [[u_snap[c][i] * inv_du for i in range(n_sites)] for c in range(3)]
    r_plus, r_minus = varga_sos_departure_rates_3d(grid, u_eff, kappa, c_s=c_s)
    nine_fourths = Fraction(9, 4)
    one_third = Fraction(1, 3)
    inv_tame_vec = []
    for y_idx in range(n_sites):
        ux_y = u_eff[0][y_idx]
        uy_y = u_eff[1][y_idx]
        uz_y = u_eff[2][y_idx]
        usq_y = ux_y * ux_y + uy_y * uy_y + uz_y * uz_y
        db_y = 6 * kappa + nine_fourths + one_third * usq_y
        dy_q = q_nonneg(HALF_Q * db_y * db_y)
        inv_tame_vec.append(ONE_Q / (ONE_Q + dt_sq * dy_q))
    return r_plus, r_minus, inv_tame_vec


def varga_density_step_core_3d(
    grid,
    rho_field,
    u_vec3,
    kappa,
    dt,
    c_s=DEFAULT_DENSITY_SOUND_SPEED,
    q_bundle=None,
    u_vec3_end=None,
):
    """Unified 2nd-order Heun-SOS + Varga squared-M-matrix positive density step (`Q`, `QQRat`, `QQJet`).

    Satisfies the exact algebraic identity:
      `rho_heun == (1/2) * rho_field + (1/2) * (I + dt * L_2)((I + dt * L_1) rho_field) >= (1/2) * rho_field > 0/1`
    for all `dt in Q` and any endpoint velocity pair `(u_vec3, u_vec3_end)`.
    """
    qb = resolve_quantizer_bundle(False, bundle=q_bundle)
    n_sites = grid.N
    dt_sq = dt * dt

    rp1, rm1, it1 = _sos_rates_and_tame_3d(grid, u_vec3, kappa, dt_sq, c_s, qb)
    net_1 = _conservative_sos_net_rate_3d(grid, rho_field, rp1, rm1, it1, qb)
    rho_1 = [rho_field[i] + dt * net_1[i] for i in range(n_sites)]

    if u_vec3_end is None:
        rp2, rm2, it2 = rp1, rm1, it1
    else:
        rp2, rm2, it2 = _sos_rates_and_tame_3d(grid, u_vec3_end, kappa, dt_sq, c_s, qb)
    net_2 = _conservative_sos_net_rate_3d(grid, rho_1, rp2, rm2, it2, qb)
    rho_heun = [rho_field[i] + (HALF_Q * dt) * (net_1[i] + net_2[i]) for i in range(n_sites)]
    return _apply_unified_varga_3d(grid, rho_heun, kappa, dt_sq, qb)


# --- Module: spatial_3d.py ---
ZERO_Q = Fraction(0, 1)
ONE_Q = Fraction(1, 1)
HALF_Q = Fraction(1, 2)
TWO_THIRDS_Q = Fraction(2, 3)
FOUR_THIRDS_Q = Fraction(4, 3)
TWELFTH_Q = Fraction(1, 12)


def _exact_cube_root(n_sites):
    """Exact integer cube root L of N = L^3 via bounded integer search (zero `while` loops)."""
    for cand in range(1, 65):
        if cand * cand * cand == n_sites:
            return cand
    raise ValueError(f"N={n_sites} is not a supported cube L^3 (1 <= L <= 64)")


def bound_B_u(E0, N, zero_momentum=True, h=ONE_Q):
    """Exact, sharp closed-form rational upper bound on `max_x ||u_n(x)||^2` over Q.

    Attained on `V_{div,0}` (or `V_{div}`) by the reproducing-kernel column `P_{div,0} e_{x,0}`:
      `max_x ||u_n(x)||^2 <= 2 * (rank(Proj) / (3 * N)) * E_0`.
    """
    L = _exact_cube_root(N)
    if (L % 2) == 0:
        rank_proj = (2 * N + 5) if zero_momentum else (2 * N + 8)
    else:
        rank_proj = (2 * N - 2) if zero_momentum else (2 * N + 1)
    c_u = Fraction(rank_proj, 3 * N)
    return 2 * c_u * E0


def bound_B_omega(E0, N, h=ONE_Q):
    """Exact, sharp closed-form rational upper bound on `max_x ||omega_n(x)||^2` over Q.

    Attained with equality on `Z_L^3` (`L >= 5`) by `u* = curl_Sym(delta_{x_0} e_x)`:
      `max_x ||omega_n(x)||^2 <= (65 / (18 * h^2)) * E_0`.
    """
    L = _exact_cube_root(N)
    w_vec = [ZERO_Q for idx in range(L)]
    w_vec[1 % L] += TWO_THIRDS_Q / h
    w_vec[(-1) % L] -= TWO_THIRDS_Q / h
    w_vec[2 % L] -= TWELFTH_Q / h
    w_vec[(-2) % L] += TWELFTH_Q / h
    norm_sq_1d = ZERO_Q
    for idx in range(L):
        norm_sq_1d = norm_sq_1d + w_vec[idx] * w_vec[idx]
    c_omega = 2 * norm_sq_1d
    return 2 * c_omega * E0


class SymanzikTorus3D:
    """Fast exact rational 3D Symanzik vector calculus & separable Leray projector on `Z_L^3`."""

    def __bool__(self):
        return True

    def __len__(self):
        return 0

    def __init__(self, L, h=ONE_Q):
        if L < 5:
            raise ValueError("SymanzikTorus3D requires L >= 5 for range-2 stencils")
        self.L = L
        self.N = L * L * L
        self.h = Fraction(h)

        # Precompute 3D shift index tables for axes 0, 1, 2 (steps +1, -1, +2, -2)
        self.p1 = [[], [], []]
        self.m1 = [[], [], []]
        self.p2 = [[], [], []]
        self.m2 = [[], [], []]
        for x_i in range(L):
            for y_i in range(L):
                for z_i in range(L):
                    self.p1[0].append(self.idx(x_i + 1, y_i, z_i))
                    self.m1[0].append(self.idx(x_i - 1, y_i, z_i))
                    self.p2[0].append(self.idx(x_i + 2, y_i, z_i))
                    self.m2[0].append(self.idx(x_i - 2, y_i, z_i))

                    self.p1[1].append(self.idx(x_i, y_i + 1, z_i))
                    self.m1[1].append(self.idx(x_i, y_i - 1, z_i))
                    self.p2[1].append(self.idx(x_i, y_i + 2, z_i))
                    self.m2[1].append(self.idx(x_i, y_i - 2, z_i))

                    self.p1[2].append(self.idx(x_i, y_i, z_i + 1))
                    self.m1[2].append(self.idx(x_i, y_i, z_i - 1))
                    self.p2[2].append(self.idx(x_i, y_i, z_i + 2))
                    self.m2[2].append(self.idx(x_i, y_i, z_i - 2))

        for ax in (0, 1, 2):
            self.p1[ax] = tuple(self.p1[ax])
            self.m1[ax] = tuple(self.m1[ax])
            self.p2[ax] = tuple(self.p2[ax])
            self.m2[ax] = tuple(self.m2[ax])

        # Precompute exact rational 3D Green kernel G_div = (-D_x^2 - D_y^2 - D_z^2)^+ via O_h orbits
        # AND integer-factored separable 1D circulant representation!
        self.g_div_row0 = self._precompute_3d_green_kernel(operator_kind="div")
        self.g_sym_row0 = self._precompute_3d_green_kernel(operator_kind="sym")
        self._init_fast_separable_convolution()
        (
            self.shell_projectors_1d,
            self.shell_int_projectors_1d,
            self.shell_roughness_1d,
        ) = self._precompute_1d_roughness_subspaces()

    def idx(self, x_i, y_i, z_i):
        L = self.L
        return ((x_i % L) * L + (y_i % L)) * L + (z_i % L)

    def coords(self, flat_idx):
        L = self.L
        x_i = flat_idx // (L * L)
        rem = flat_idx % (L * L)
        return (x_i, rem // L, rem % L)

    def _canonical_cubic_orbit_rep(self, x_i, y_i, z_i):
        L = self.L
        ax = x_i % L
        if L - ax < ax:
            ax = L - ax
        ay = y_i % L
        if L - ay < ay:
            ay = L - ay
        az = z_i % L
        if L - az < az:
            az = L - az
        vals = sorted((ax, ay, az))
        return (vals[0], vals[1], vals[2])

    def _precompute_3d_green_kernel(self, operator_kind="div"):
        """Compute exact Moore-Penrose pseudoinverse first row over Q via cubic O_h symmetry reduction."""
        L = self.L
        N = self.N

        ker_vectors = []
        if operator_kind == "div" and (L % 2) == 0:
            for a_bit in (0, 1):
                for b_bit in (0, 1):
                    for c_bit in (0, 1):
                        vec = []
                        for x_i in range(L):
                            for y_i in range(L):
                                for z_i in range(L):
                                    parity = (a_bit * x_i + b_bit * y_i + c_bit * z_i) & 1
                                    vec.append(Fraction(-1 if parity else 1, 1))
                        ker_vectors.append(vec)
        else:
            ker_vectors.append([ONE_Q for idx in range(N)])

        orbit_map = {}
        site_to_orbit = [0] * N
        orbit_members = []
        for flat_i in range(N):
            xi, yi, zi = self.coords(flat_i)
            rep = self._canonical_cubic_orbit_rep(xi, yi, zi)
            if rep not in orbit_map:
                orbit_map[rep] = len(orbit_members)
                orbit_members.append([])
            oid = orbit_map[rep]
            site_to_orbit[flat_i] = oid
            orbit_members[oid].append(flat_i)

        n_orbits = len(orbit_members)
        inv_n = Fraction(1, N)

        mat_orbits = []
        for col_oid in range(n_orbits):
            basis_orb = [ZERO_Q] * N
            for member_idx in orbit_members[col_oid]:
                basis_orb[member_idx] = ONE_Q
            if operator_kind == "div":
                op_col = self.apply_L_div(basis_orb)
            else:
                op_col = self.apply_K_sym_scalar(basis_orb)

            for kv in ker_vectors:
                dot_kv = ZERO_Q
                for member_idx in orbit_members[col_oid]:
                    dot_kv = dot_kv + kv[member_idx]
                coeff_kv = dot_kv * inv_n
                if coeff_kv != ZERO_Q:
                    for site_i in range(N):
                        op_col[site_i] = op_col[site_i] + coeff_kv * kv[site_i]

            mat_orbits.append([op_col[orbit_members[row_oid][0]] for row_oid in range(n_orbits)])

        mat_rows = [
            [mat_orbits[col_oid][row_oid] for col_oid in range(n_orbits)]
            for row_oid in range(n_orbits)
        ]

        p_ker_e0 = [ZERO_Q] * N
        for kv in ker_vectors:
            c0 = kv[0] * inv_n
            for site_i in range(N):
                p_ker_e0[site_i] = p_ker_e0[site_i] + c0 * kv[site_i]

        rhs_orbits = []
        for row_oid in range(n_orbits):
            rep_site = orbit_members[row_oid][0]
            val = (ONE_Q if rep_site == 0 else ZERO_Q) - p_ker_e0[rep_site]
            rhs_orbits.append(val)

        sol_orbits = solve_linear_system_q(mat_rows, rhs_orbits)
        row0 = [sol_orbits[site_to_orbit[site_i]] for site_i in range(N)]
        return tuple(row0)

    def _init_fast_separable_convolution(self):
        """Build fast exact O(N * L) 3D circulant application with factored integer common denominator."""
        L = self.L
        half_m = (L // 2) + 1
        self.half_m = half_m

        d_common = 1
        for q_val in self.g_div_row0:
            if q_val != ZERO_Q:
                d_common = _lcm_positive_ints_q(d_common, q_val.denominator)
        self.g_div_inv_common_den = Fraction(1, d_common)

        w_div = [
            [[ZERO_Q for az in range(half_m)] for ay in range(half_m)]
            for ax in range(half_m)
        ]
        w_div_int = [
            [[0 for az in range(half_m)] for ay in range(half_m)]
            for ax in range(half_m)
        ]
        for ax in range(half_m):
            for ay in range(half_m):
                for az in range(half_m):
                    q_entry = self.g_div_row0[self.idx(ax, ay, az)]
                    w_div[ax][ay][az] = q_entry
                    w_div_int[ax][ay][az] = (q_entry * d_common).numerator
        self.w_div_folded = tuple(
            tuple(tuple(w_div[ax][ay][az] for az in range(half_m)) for ay in range(half_m))
            for ax in range(half_m)
        )
        self.w_div_int_folded = tuple(
            tuple(tuple(w_div_int[ax][ay][az] for az in range(half_m)) for ay in range(half_m))
            for ax in range(half_m)
        )
        self.g_div_common_den_int = d_common
        self.g_div_common_den = Fraction(d_common, 1)

        lag_partners = []
        for x_i in range(L):
            x_lags = []
            for d_lag in range(half_m):
                if d_lag == 0:
                    x_lags.append((x_i,))
                elif (L % 2 == 0) and (d_lag == L // 2):
                    x_lags.append(((x_i + d_lag) % L,))
                else:
                    x_lags.append(((x_i + d_lag) % L, (x_i - d_lag) % L))
            lag_partners.append(tuple(x_lags))
        self.lag_partners_1d = tuple(lag_partners)

    def apply_G_div_int(self, scalar_field):
        """Apply integer-weighted folded `g_div_common_den * G_div(scalar_field)` (C PyLong on int)."""
        return self._apply_int_folded_3d_raw(scalar_field, self.w_div_int_folded)

    def apply_G_div(self, scalar_field):
        """Apply exact `G_div = (-D_x^2 - D_y^2 - D_z^2)^+` via integer-weighted folded lags.

        Zero FFTs, zero floats, 1 division per site (`< 3 ms` at `L=8`), 100% duck-typed
        over `Fraction`, `QQPoly`, `QQRat`, `QQJet`, and `AuditJet`.
        """
        inv_common_den = self.g_div_inv_common_den
        return [inv_common_den * val for val in self.apply_G_div_int(scalar_field)]

    def _apply_int_folded_3d_raw(self, scalar_field, w_int_fold):
        L = self.L
        L2 = L * L
        half_m = self.half_m
        partners = self.lag_partners_1d

        # Fold 1 along z: pure additions (initialized with Python int 0 for C PyLong speed + CAS duck-typing)
        fz = [[0] * self.N for az in range(half_m)]
        for z_out in range(L):
            pz = partners[z_out]
            for az in range(half_m):
                tup = pz[az]
                dst_az = fz[az]
                if len(tup) == 1:
                    z0 = tup[0]
                    for xy in range(L2):
                        base = xy * L
                        dst_az[base + z_out] = scalar_field[base + z0]
                else:
                    z0, z1 = tup[0], tup[1]
                    for xy in range(L2):
                        base = xy * L
                        dst_az[base + z_out] = scalar_field[base + z0] + scalar_field[base + z1]

        # Fold 2 along y: pure additions
        fyz = [[[0] * self.N for az in range(half_m)] for ay in range(half_m)]
        for ay in range(half_m):
            for az in range(half_m):
                src = fz[az]
                dst = fyz[ay][az]
                for y_out in range(L):
                    tup = partners[y_out][ay]
                    if len(tup) == 1:
                        y0 = tup[0]
                        for x_i in range(L):
                            bx = x_i * L2
                            by_out = bx + y_out * L
                            by0 = bx + y0 * L
                            for z_i in range(L):
                                dst[by_out + z_i] = src[by0 + z_i]
                    else:
                        y0, y1 = tup[0], tup[1]
                        for x_i in range(L):
                            bx = x_i * L2
                            by_out = bx + y_out * L
                            by0 = bx + y0 * L
                            by1 = bx + y1 * L
                            for z_i in range(L):
                                dst[by_out + z_i] = src[by0 + z_i] + src[by1 + z_i]

        # Fold 3 along x with exact integer weights
        out_int = [0] * self.N
        for ax in range(half_m):
            for x_out in range(L):
                tup = partners[x_out][ax]
                bx_out = x_out * L2
                if len(tup) == 1:
                    bx0 = tup[0] * L2
                    for ay in range(half_m):
                        w_xy = w_int_fold[ax][ay]
                        for az in range(half_m):
                            w_int = w_xy[az]
                            if w_int == 0:
                                continue
                            src = fyz[ay][az]
                            for yz in range(L2):
                                out_int[bx_out + yz] = out_int[bx_out + yz] + w_int * src[bx0 + yz]
                else:
                    bx0 = tup[0] * L2
                    bx1 = tup[1] * L2
                    for ay in range(half_m):
                        w_xy = w_int_fold[ax][ay]
                        for az in range(half_m):
                            w_int = w_xy[az]
                            if w_int == 0:
                                continue
                            src = fyz[ay][az]
                            for yz in range(L2):
                                out_int[bx_out + yz] = out_int[bx_out + yz] + w_int * (
                                    src[bx0 + yz] + src[bx1 + yz]
                                )
        return out_int

    def _precompute_1d_roughness_subspaces(self):
        """Construct exact Z_L-circulant Q-irreducible cyclotomic projectors `Pi_s` (`[Pi_s, T] == 0/1`, `sum Pi_s == I_L`)."""
        L = self.L
        if L == 6:
            b0 = [[Fraction(1, 1) for x in range(6)]]
            b1 = [
                [Fraction(2, 1), Fraction(1, 1), Fraction(-1, 1), Fraction(-2, 1), Fraction(-1, 1), Fraction(1, 1)],
                [Fraction(0, 1), Fraction(1, 1), Fraction(1, 1), Fraction(0, 1), Fraction(-1, 1), Fraction(-1, 1)],
            ]
            b2 = [
                [Fraction(2, 1), Fraction(-1, 1), Fraction(-1, 1), Fraction(2, 1), Fraction(-1, 1), Fraction(-1, 1)],
                [Fraction(0, 1), Fraction(1, 1), Fraction(-1, 1), Fraction(0, 1), Fraction(1, 1), Fraction(-1, 1)],
            ]
            b3 = [[Fraction(1, 1), Fraction(-1, 1), Fraction(1, 1), Fraction(-1, 1), Fraction(1, 1), Fraction(-1, 1)]]
            groups = [(0, b0), (1, b1), (3, b2), (4, b3)]
        elif L == 8:
            # Exact Z_8-circulant Q-irreducible cyclotomic projectors ([Pi_s, T] == 0/1, sum Pi_s == I_8):
            #   r=0 (Phi_1, dim 1): Pi_0 = (1/8) sum_{k=0..7} T^k
            #   r=1 (Phi_8, dim 4): Pi_odd = (1/2) (I_8 - T^4)  [unsplit 4D odd-frequency shell]
            #   r=2 (Phi_4, dim 2): Pi_2 = (1/4) (I_8 - T^2 + T^4 - T^6)
            #   r=4 (Phi_2, dim 1): Pi_4 = (1/8) sum_{k=0..7} (-1)^k T^k
            circ_specs = (
                (0, Fraction(1, 8), (1, 1, 1, 1, 1, 1, 1, 1)),
                (1, Fraction(1, 2), (1, 0, 0, 0, -1, 0, 0, 0)),
                (2, Fraction(1, 4), (1, 0, -1, 0, 1, 0, -1, 0)),
                (4, Fraction(1, 8), (1, -1, 1, -1, 1, -1, 1, -1)),
            )
            projs = []
            int_projs = []
            r_vals = []
            for r_label, inv_den, int_row0 in circ_specs:
                int_mat = tuple(
                    tuple(int_row0[(r_i - c_j) % 8] for c_j in range(8))
                    for r_i in range(8)
                )
                mat = tuple(
                    tuple(inv_den * int_mat[r_i][c_j] for c_j in range(8))
                    for r_i in range(8)
                )
                projs.append(mat)
                int_projs.append((inv_den, int_mat))
                r_vals.append(r_label)
            return tuple(projs), tuple(int_projs), tuple(r_vals)
        else:
            inv_l = Fraction(1, L)
            p0 = tuple(tuple(inv_l for c_j in range(L)) for r_i in range(L))
            p_rest = tuple(
                tuple((ONE_Q - inv_l) if r_i == c_j else (-inv_l) for c_j in range(L))
                for r_i in range(L)
            )
            int0 = tuple(tuple(1 for c_j in range(L)) for r_i in range(L))
            int_rest = tuple(
                tuple((L - 1) if r_i == c_j else -1 for c_j in range(L))
                for r_i in range(L)
            )
            return (p0, p_rest), ((inv_l, int0), (inv_l, int_rest)), (0, 1)

        projs = []
        int_projs = []
        r_vals = []
        for r_label, basis_list in groups:
            mat = [[ZERO_Q for c_j in range(L)] for r_i in range(L)]
            for vec in basis_list:
                nsq = ZERO_Q
                for x_i in range(L):
                    nsq = nsq + vec[x_i] * vec[x_i]
                inv_nsq = ONE_Q / nsq
                for r_i in range(L):
                    for c_j in range(L):
                        mat[r_i][c_j] = mat[r_i][c_j] + (vec[r_i] * vec[c_j]) * inv_nsq
            d_lcm = 1
            for r_i in range(L):
                for c_j in range(L):
                    if mat[r_i][c_j] != ZERO_Q:
                        d_lcm = _lcm_positive_ints_q(d_lcm, mat[r_i][c_j].denominator)
            int_mat = tuple(
                tuple((mat[r_i][c_j] * d_lcm).numerator for c_j in range(L))
                for r_i in range(L)
            )
            projs.append(tuple(tuple(row) for row in mat))
            int_projs.append((Fraction(1, d_lcm), int_mat))
            r_vals.append(r_label)
        return tuple(projs), tuple(int_projs), tuple(r_vals)

    # --- 3D Vector Calculus Operators ---

    def apply_D_axis_12h(self, field, axis):
        """Division-free `12 * h * D_axis` (`8*(T_{+1}-T_{-1}) - (T_{+2}-T_{-2})`, preserves `den == 1`)."""
        p1, m1 = self.p1[axis], self.m1[axis]
        p2, m2 = self.p2[axis], self.m2[axis]
        return [
            8 * (field[p1[i]] - field[m1[i]]) - (field[p2[i]] - field[m2[i]])
            for i in range(self.N)
        ]

    def apply_D_axis(self, field, axis):
        """4th-order skew-adjoint derivative `(1/h)((2/3)(T_{+1} - T_{-1}) - (1/12)(T_{+2} - T_{-2}))`."""
        inv_12h = TWELFTH_Q / self.h
        return [inv_12h * val for val in self.apply_D_axis_12h(field, axis)]

    def grad_sym_12h(self, phi):
        """Division-free `12 * h * grad_Sym(phi)`."""
        return [
            self.apply_D_axis_12h(phi, 0),
            self.apply_D_axis_12h(phi, 1),
            self.apply_D_axis_12h(phi, 2),
        ]

    def grad_sym(self, phi):
        """4th-order Symanzik gradient `[D_x phi, D_y phi, D_z phi]`."""
        return [
            self.apply_D_axis(phi, 0),
            self.apply_D_axis(phi, 1),
            self.apply_D_axis(phi, 2),
        ]

    def div_sym_12h(self, u_vec):
        """Division-free `12 * h * div_Sym(u_vec)`."""
        dx = self.apply_D_axis_12h(u_vec[0], 0)
        dy = self.apply_D_axis_12h(u_vec[1], 1)
        dz = self.apply_D_axis_12h(u_vec[2], 2)
        return [dx[i] + dy[i] + dz[i] for i in range(self.N)]

    def div_sym(self, u_vec):
        """4th-order Symanzik divergence `D_x u_x + D_y u_y + D_z u_z`."""
        inv_12h = TWELFTH_Q / self.h
        return [inv_12h * val for val in self.div_sym_12h(u_vec)]

    def curl_sym(self, u_vec):
        """4th-order Symanzik curl `nabla_Sym x u`."""
        ux, uy, uz = u_vec[0], u_vec[1], u_vec[2]
        dy_uz = self.apply_D_axis(uz, 1)
        dz_uy = self.apply_D_axis(uy, 2)
        dz_ux = self.apply_D_axis(ux, 2)
        dx_uz = self.apply_D_axis(uz, 0)
        dx_uy = self.apply_D_axis(uy, 0)
        dy_ux = self.apply_D_axis(ux, 1)
        n_sites = self.N
        return [
            [dy_uz[i] - dz_uy[i] for i in range(n_sites)],
            [dz_ux[i] - dx_uz[i] for i in range(n_sites)],
            [dx_uy[i] - dy_ux[i] for i in range(n_sites)],
        ]

    def apply_L_div(self, phi):
        """Apply positive semi-definite divergence-gradient Laplacian `L_div = -div_Sym(grad_Sym(phi))`."""
        g_phi = self.grad_sym(phi)
        d_g = self.div_sym(g_phi)
        return [-val for val in d_g]

    def apply_K_sym_12h2_scalar(self, field):
        """Division-free `12 * h^2 * K_Sym(field)` (preserves `den == 1`)."""
        n_sites = self.N
        p1x, m1x, p2x, m2x = self.p1[0], self.m1[0], self.p2[0], self.m2[0]
        p1y, m1y, p2y, m2y = self.p1[1], self.m1[1], self.p2[1], self.m2[1]
        p1z, m1z, p2z, m2z = self.p1[2], self.m1[2], self.p2[2], self.m2[2]
        out = []
        for i in range(n_sites):
            fi = field[i]
            six_fi = 6 * fi
            l1 = six_fi - (
                field[p1x[i]] + field[m1x[i]]
                + field[p1y[i]] + field[m1y[i]]
                + field[p1z[i]] + field[m1z[i]]
            )
            l2 = six_fi - (
                field[p2x[i]] + field[m2x[i]]
                + field[p2y[i]] + field[m2y[i]]
                + field[p2z[i]] + field[m2z[i]]
            )
            out.append(16 * l1 - l2)
        return out

    def apply_K_sym_12h2_vec3(self, u_vec):
        """Division-free `12 * h^2 * K_Sym` applied componentwise to `[ux, uy, uz]`."""
        return [
            self.apply_K_sym_12h2_scalar(u_vec[0]),
            self.apply_K_sym_12h2_scalar(u_vec[1]),
            self.apply_K_sym_12h2_scalar(u_vec[2]),
        ]

    def apply_K_sym_scalar(self, field):
        """Apply 3D 4th-order Symanzik scalar Laplacian `(1/h^2) sum_{ax} ((4/3)L_1 - (1/12)L_2)`."""
        inv_12h2 = TWELFTH_Q / (self.h * self.h)
        return [inv_12h2 * val for val in self.apply_K_sym_12h2_scalar(field)]

    def apply_K_sym_vec3(self, u_vec):
        """Apply 3D 4th-order Symanzik Laplacian componentwise to `[ux, uy, uz]`."""
        return [
            self.apply_K_sym_scalar(u_vec[0]),
            self.apply_K_sym_scalar(u_vec[1]),
            self.apply_K_sym_scalar(u_vec[2]),
        ]

    def project_leray(self, w_vec):
        """Exact rational Leray projector `P_div(w) = w + grad_Sym(G_div(div_Sym(w)))`.

        Evaluates `div_sym_12h -> apply_G_div_int -> grad_sym_12h` before the single scalar
        division so `denominator == 1` inputs execute zero intermediate GCDs!
        """
        div_12h = self.div_sym_12h(w_vec)
        pot_int = self.apply_G_div_int(div_12h)
        grad_144h2 = self.grad_sym_12h(pot_int)
        scale = self.g_div_inv_common_den / (144 * self.h * self.h)
        n_sites = self.N
        return [
            [w_vec[0][i] + scale * grad_144h2[0][i] for i in range(n_sites)],
            [w_vec[1][i] + scale * grad_144h2[1][i] for i in range(n_sites)],
            [w_vec[2][i] + scale * grad_144h2[2][i] for i in range(n_sites)],
        ]

    def project_leray_zero_mean_int(self, w_vec_int):
        """Division-free `(144 * h^2 * g_div_common_den * N) * P_{div,0}(w_vec_int)` (C PyLong on int)."""
        n_sites = self.N
        div_12h = self.div_sym_12h(w_vec_int)
        pot_int = self.apply_G_div_int(div_12h)
        grad_144h2 = self.grad_sym_12h(pot_int)
        h_num = self.h.numerator
        h_den = self.h.denominator
        h_num_sq = h_num * h_num
        h_den_sq = h_den * h_den
        base_scale_int = 144 * h_num_sq * self.g_div_common_den_int
        total_scale_int = base_scale_int * n_sites
        grad_scale_int = n_sites * h_den_sq
        out = []
        for comp in range(3):
            wc = w_vec_int[comp]
            gc = grad_144h2[comp]
            sum_wc = 0
            for i in range(n_sites):
                sum_wc = sum_wc + wc[i]
            # Note: sum_i gc[i] == 0 identically because grad_Sym is skew-adjoint on 1!
            mean_scaled = base_scale_int * sum_wc
            out.append(
                [
                    total_scale_int * wc[i] + grad_scale_int * gc[i] - mean_scaled
                    for i in range(n_sites)
                ]
            )
        return out, Fraction(total_scale_int, h_den_sq)

    def project_leray_zero_mean(self, w_vec):
        """Exact rational orthogonal projector `P_{div,0} = P_div - P_mean` onto zero-mean incompressible fields."""
        u_div = self.project_leray(w_vec)
        inv_n = Fraction(1, self.N)
        out = []
        for comp in range(3):
            uc = u_div[comp]
            acc = ZERO_Q
            for i in range(self.N):
                acc = acc + uc[i]
            mean_c = acc * inv_n
            out.append([uc[i] - mean_c for i in range(self.N)])
        return out

    def skew_advection_24h(self, u_vec):
        """Division-free `24 * h * skew_advection(u_vec)` (pure C PyLong on int, duck-typed on Q/CAS)."""
        n_sites = self.N
        ux, uy, uz = u_vec[0], u_vec[1], u_vec[2]
        u_list = (ux, uy, uz)
        out = []
        for comp_c in range(3):
            uc = u_list[comp_c]
            acc = [0 for i in range(n_sites)]
            for axis_a in range(3):
                ua = u_list[axis_a]
                d12_uc = self.apply_D_axis_12h(uc, axis_a)
                prod_ac = [ua[i] * uc[i] for i in range(n_sites)]
                d12_prod = self.apply_D_axis_12h(prod_ac, axis_a)
                for i in range(n_sites):
                    acc[i] = acc[i] - (ua[i] * d12_uc[i] + d12_prod[i])
            out.append(acc)
        return out

    def skew_advection(self, u_vec):
        """Energy- and momentum-skew-symmetric 3D advection `-1/2 sum_a (u_a D_a u_c + D_a(u_a u_c))`."""
        inv_24h = Fraction(1, 24) / self.h
        raw = self.skew_advection_24h(u_vec)
        return [
            [inv_24h * val for val in raw[0]],
            [inv_24h * val for val in raw[1]],
            [inv_24h * val for val in raw[2]],
        ]

    def _apply_1d_proj_axis(self, field, proj_1d, axis):
        L = self.L
        L2 = L * L
        out = [ZERO_Q] * self.N
        if axis == 0:
            for x_out in range(L):
                row = proj_1d[x_out]
                bx_out = x_out * L2
                for x_in in range(L):
                    w_val = row[x_in]
                    if w_val == ZERO_Q:
                        continue
                    bx_in = x_in * L2
                    for yz in range(L2):
                        out[bx_out + yz] = out[bx_out + yz] + w_val * field[bx_in + yz]
        elif axis == 1:
            for y_out in range(L):
                row = proj_1d[y_out]
                for y_in in range(L):
                    w_val = row[y_in]
                    if w_val == ZERO_Q:
                        continue
                    for x_i in range(L):
                        bx = x_i * L2
                        by_out = bx + y_out * L
                        by_in = bx + y_in * L
                        for z_i in range(L):
                            out[by_out + z_i] = out[by_out + z_i] + w_val * field[by_in + z_i]
        else:
            for z_out in range(L):
                row = proj_1d[z_out]
                for z_in in range(L):
                    w_val = row[z_in]
                    if w_val == ZERO_Q:
                        continue
                    for xy in range(L2):
                        base = xy * L
                        out[base + z_out] = out[base + z_out] + w_val * field[base + z_in]
        return out

    def _apply_1d_int_proj_axis_raw(self, field, int_mat, axis):
        """Division-free integer-weighted 1D projection along `axis`."""
        L = self.L
        L2 = L * L
        out_int = [ZERO_Q] * self.N
        if axis == 0:
            for x_out in range(L):
                row = int_mat[x_out]
                bx_out = x_out * L2
                for x_in in range(L):
                    w_int = row[x_in]
                    if w_int == 0:
                        continue
                    bx_in = x_in * L2
                    for yz in range(L2):
                        out_int[bx_out + yz] = out_int[bx_out + yz] + w_int * field[bx_in + yz]
        elif axis == 1:
            for y_out in range(L):
                row = int_mat[y_out]
                for y_in in range(L):
                    w_int = row[y_in]
                    if w_int == 0:
                        continue
                    for x_i in range(L):
                        bx = x_i * L2
                        by_out = bx + y_out * L
                        by_in = bx + y_in * L
                        for z_i in range(L):
                            out_int[by_out + z_i] = out_int[by_out + z_i] + w_int * field[by_in + z_i]
        else:
            for z_out in range(L):
                row = int_mat[z_out]
                for z_in in range(L):
                    w_int = row[z_in]
                    if w_int == 0:
                        continue
                    for xy in range(L2):
                        base = xy * L
                        out_int[base + z_out] = out_int[base + z_out] + w_int * field[base + z_in]
        return out_int

    def _apply_1d_int_proj_axis(self, field, inv_den, int_mat, axis):
        """Integer-factored 1D projection along `axis` (`1` division per site)."""
        out_int = self._apply_1d_int_proj_axis_raw(field, int_mat, axis)
        return [inv_den * val for val in out_int]

    def roughness_shell_energy_spectrum(self, u_vec):
        """Exact rational 3D roughness-shell energy spectrum `dict[tilde_R -> Fraction]` (`sum == E_kin` to `0/1`)."""
        int_projs = self.shell_int_projectors_1d
        r_labels = self.shell_roughness_1d
        n_1d = len(int_projs)
        shell_energies = {}
        for ix in range(n_1d):
            for iy in range(n_1d):
                for iz in range(n_1d):
                    r_tot = r_labels[ix] + r_labels[iy] + r_labels[iz]
                    if r_tot not in shell_energies:
                        shell_energies[r_tot] = ZERO_Q

        for comp in range(3):
            fc = u_vec[comp]
            for ix in range(n_1d):
                inv_dx, mat_x = int_projs[ix]
                fx = self._apply_1d_int_proj_axis_raw(fc, mat_x, 0)
                for iy in range(n_1d):
                    inv_dy, mat_y = int_projs[iy]
                    inv_xy = inv_dx * inv_dy
                    fxy = self._apply_1d_int_proj_axis_raw(fx, mat_y, 1)
                    for iz in range(n_1d):
                        inv_dz, mat_z = int_projs[iz]
                        inv_xyz = inv_xy * inv_dz
                        fxyz = self._apply_1d_int_proj_axis_raw(fxy, mat_z, 2)
                        r_tot = r_labels[ix] + r_labels[iy] + r_labels[iz]
                        nsq = ZERO_Q
                        for i in range(self.N):
                            nsq = nsq + fxyz[i] * fxyz[i]
                        shell_energies[r_tot] = (
                            shell_energies[r_tot] + HALF_Q * (inv_xyz * inv_xyz) * nsq
                        )

        return shell_energies

    def varga_1d_kernel_sq(self, step_alpha):
        """Exact closed-form rational 1st row of `(I + step_alpha * L_1d)^(-2)` on `Z_L`."""
        L = self.L
        if L == 6:
            p0, p1_m, p3_m, p4_m = self.shell_projectors_1d
            d1 = ONE_Q + step_alpha
            d3 = ONE_Q + 3 * step_alpha
            d4 = ONE_Q + 4 * step_alpha
            w1 = ONE_Q / (d1 * d1)
            w3 = ONE_Q / (d3 * d3)
            w4 = ONE_Q / (d4 * d4)
            return [
                p0[0][d] + w1 * p1_m[0][d] + w3 * p3_m[0][d] + w4 * p4_m[0][d]
                for d in range(6)
            ]
        if L == 8:
            a2 = step_alpha * step_alpha
            d2 = ONE_Q + 2 * step_alpha
            d4 = ONE_Q + 4 * step_alpha
            d_odd = ONE_Q + 4 * step_alpha + 2 * a2
            inv_d_odd_sq = ONE_Q / (d_odd * d_odd)
            w2 = ONE_Q / (d2 * d2)
            w4 = ONE_Q / (d4 * d4)
            c_id = (ONE_Q + 8 * step_alpha + 14 * a2) * inv_d_odd_sq
            c_lap = -(2 * step_alpha + 4 * a2) * inv_d_odd_sq

            eighth = Fraction(1, 8)
            quarter = Fraction(1, 4)
            p0_row = (eighth, eighth, eighth, eighth, eighth, eighth, eighth, eighth)
            p2_row = (quarter, ZERO_Q, -quarter, ZERO_Q, quarter, ZERO_Q, -quarter, ZERO_Q)
            p4_row = (eighth, -eighth, eighth, -eighth, eighth, -eighth, eighth, -eighth)
            p_odd_row = (HALF_Q, ZERO_Q, ZERO_Q, ZERO_Q, -HALF_Q, ZERO_Q, ZERO_Q, ZERO_Q)
            l_odd_row = (ONE_Q, -HALF_Q, ZERO_Q, HALF_Q, -ONE_Q, HALF_Q, ZERO_Q, -HALF_Q)
            return [
                p0_row[d]
                + w2 * p2_row[d]
                + w4 * p4_row[d]
                + c_id * p_odd_row[d]
                + c_lap * l_odd_row[d]
                for d in range(8)
            ]
        raise ValueError(f"Closed-form Varga kernel implemented for L in (6, 8), got L={L}")

    def apply_varga_3d(self, rho_field, step_alpha):
        """Apply separable 3D Varga squared M-matrix resolvent `(I + step_alpha * L_1d)^(-2)(x)3`."""
        L = self.L
        r2 = self.varga_1d_kernel_sq(step_alpha)
        circ_mat = [[r2[(r_i - c_j) % L] for c_j in range(L)] for r_i in range(L)]
        f_x = self._apply_1d_proj_axis(rho_field, circ_mat, 0)
        f_xy = self._apply_1d_proj_axis(f_x, circ_mat, 1)
        return self._apply_1d_proj_axis(f_xy, circ_mat, 2)


_GRID_3D_CACHE = {}


def get_grid_3d(L, h=ONE_Q):
    """Cached factory for `SymanzikTorus3D(L, h)`."""
    key = (L, Fraction(h))
    if key not in _GRID_3D_CACHE:
        _GRID_3D_CACHE[key] = SymanzikTorus3D(L, h=h)
    return _GRID_3D_CACHE[key]


# --- Module: battles.py ---
# Alias retained for symbolic Sturm verification compatibility:
_b3_varga_offdiag_over_dt2_1d = _varga_offdiag_over_dt2_1d


def _b_q_vec(vec):
    """Dynamic module-level vector quantizer wrapper (supports external tie instrumentation)."""
    return quantize_dyadic_vector(
        vec, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=CARRIER_MAX_SHIFT
    )


def _b_q_vec3(vec3):
    """Dynamic module-level 3-vector quantizer wrapper (supports external tie instrumentation)."""
    return quantize_dyadic_vec3(
        vec3, mantissa_bits=MANTISSA_BITS_DEFAULT, max_shift=CARRIER_MAX_SHIFT
    )


def _b_q_nonneg(val):
    """Dynamic module-level downward floor quantizer wrapper (supports external tie instrumentation)."""
    return quantize_dyadic_floor_nonneg(val, grid_bits=MANTISSA_BITS_DEFAULT)


_BATTLES_DYADIC_BUNDLE = dict(
    DYADIC_QUANTIZER_BUNDLE,
    vec=_b_q_vec,
    vec3=_b_q_vec3,
    nonneg=_b_q_nonneg,
)
_BATTLES_BUNDLE_DISPATCH = (EXACT_FACTORED_BUNDLE, _BATTLES_DYADIC_BUNDLE)


def _resolve_b_bundle(mode_flag=True, q_bundle=None):
    """Resolve callable quantizer bundle via 2-tuple lookup without mode branching."""
    if q_bundle is not None:
        return q_bundle
    return _BATTLES_BUNDLE_DISPATCH[bool(mode_flag)]


_ANC_2D_CACHE = {}
_ANC_3D_CACHE = {}


def _get_cached_anchor_2d(grid, a_vec):
    """Cache Bezout-Vieta anchor dictionary keyed by `(grid.m, tuple(a_vec))`."""
    key = (grid.m, tuple(a_vec))
    if key not in _ANC_2D_CACHE:
        if len(_ANC_2D_CACHE) > 32:
            _ANC_2D_CACHE.clear()
        _ANC_2D_CACHE[key] = build_bezout_vieta_anchor_2d(grid, a_vec)
    return _ANC_2D_CACHE[key]


def _build_factored_anchor_3d(grid, alpha_vec3):
    """Build exact integer-cleared 3D anchor (`EXACT_FACTORED_BUNDLE`) without tuple hashing."""
    return build_stereographic_anchor_3d(grid, alpha_vec3, q_bundle=EXACT_FACTORED_BUNDLE)


def _get_cached_anchor_3d(grid, alpha_vec3):
    """Cache integer-cleared 3D Galilean stereographic anchor keyed by `(grid.L, alpha)`."""
    key = (grid.L, tuple(tuple(comp) for comp in alpha_vec3))
    if key not in _ANC_3D_CACHE:
        if len(_ANC_3D_CACHE) > 32:
            _ANC_3D_CACHE.clear()
        _ANC_3D_CACHE[key] = build_stereographic_anchor_3d(
            grid, alpha_vec3, q_bundle=DYADIC_QUANTIZER_BUNDLE
        )
    return _ANC_3D_CACHE[key]


_ANC2D_DISPATCH = (build_bezout_vieta_anchor_2d, _get_cached_anchor_2d)
_ANC3D_DISPATCH = (_build_factored_anchor_3d, _get_cached_anchor_3d)


# =============================================================================
# BATTLE 1: 2D Gaussian Torus Z[i] / <(1 + i)^m> (1A Inviscid & 1B Viscous)
# =============================================================================


def make_battle1_initial_omega(m=6, prime_perturbation=ZERO_Q):
    """Construct exact rational zero-mean 2D unstable 4-vortex / shear merger vorticity on `G_m`."""
    grid = get_grid_2d(m)
    s2 = (ZERO_Q, ONE_Q, ZERO_Q, -ONE_Q, ZERO_Q, ONE_Q, ZERO_Q, -ONE_Q)
    c2 = (ONE_Q, ZERO_Q, -ONE_Q, ZERO_Q, ONE_Q, ZERO_Q, -ONE_Q, ZERO_Q)
    s1_8 = (
        ZERO_Q,
        Fraction(5, 7),
        ONE_Q,
        Fraction(5, 7),
        ZERO_Q,
        Fraction(-5, 7),
        -ONE_Q,
        Fraction(-5, 7),
    )
    c1_8 = (
        ONE_Q,
        Fraction(5, 7),
        ZERO_Q,
        Fraction(-5, 7),
        -ONE_Q,
        Fraction(-5, 7),
        ZERO_Q,
        Fraction(5, 7),
    )
    c16 = (
        ONE_Q,
        Fraction(12, 13),
        Fraction(5, 7),
        Fraction(5, 13),
        ZERO_Q,
        Fraction(-5, 13),
        Fraction(-5, 7),
        Fraction(-12, 13),
        -ONE_Q,
        Fraction(-12, 13),
        Fraction(-5, 7),
        Fraction(-5, 13),
        ZERO_Q,
        Fraction(5, 13),
        Fraction(5, 7),
        Fraction(12, 13),
    )
    s16 = tuple(c16[(k_i - 4) % 16] for k_i in range(16))
    alt = (ONE_Q, -ONE_Q)
    is_odd_twisted = (m % 2) == 1

    raw_om = []
    for idx_i in range(grid.N):
        x_c, y_c = grid.coords[idx_i]
        v_main = Fraction(6, 5) * (
            s2[x_c % 8] * s2[y_c % 8]
            + Fraction(4, 5) * c2[x_c % 8] * c2[y_c % 8]
            + Fraction(1, 8) * c2[x_c % 8] * s2[y_c % 8]
        )
        v_ir = Fraction(7, 20) * (
            s1_8[x_c % 8] * c1_8[y_c % 8] - Fraction(1, 3) * c1_8[x_c % 8] * s1_8[y_c % 8]
        )
        v_planes = (
            Fraction(1, 16) * c2[x_c % 8]
            + Fraction(1, 20) * c2[y_c % 8]
            + Fraction(1, 24) * c2[x_c % 8] * alt[y_c % 2]
            + Fraction(1, 28) * alt[x_c % 2] * c2[y_c % 8]
        )
        v_twist = ZERO_Q
        if is_odd_twisted:
            v_twist = (
                Fraction(3, 10) * c16[(x_c + y_c) % 16]
                - Fraction(1, 6) * s16[(x_c - y_c) % 16]
                + Fraction(9, 20) * s16[(3 * x_c + y_c) % 16]
                - Fraction(7, 20) * c16[(x_c + 3 * y_c) % 16]
            )
        pert = ZERO_Q
        if prime_perturbation != ZERO_Q:
            if idx_i == 0:
                pert = Fraction(prime_perturbation)
            elif idx_i == 1:
                pert = -Fraction(prime_perturbation)
        raw_om.append(v_main + v_ir + v_planes + v_twist + pert)

    inv_n = Fraction(1, grid.N)
    mean_om = sum(raw_om, ZERO_Q) * inv_n
    return [raw_om[i] - mean_om for i in range(grid.N)]


def continuum_rhs_battle1(arg0, m=6, nu_4=ZERO_Q):
    """Exact spatial semi-discretization RHS `d omega / dt = -{psi, omega} - nu_4 W_flat omega`.

    Supports both `(omega, m=6, nu_4=0)` and `(grid, omega, nu_4=0)` calling conventions.
    """
    if len(arg0) == 0:
        grid = arg0
        omega = m
        nu4_val = nu_4
    else:
        omega = arg0
        grid = get_grid_2d(m)
        nu4_val = nu_4
    psi = grid.apply_G0(omega)
    jac = grid.arakawa_bracket(psi, omega)
    if nu4_val == ZERO_Q:
        return [-jac[i] for i in range(grid.N)]
    w_om = grid.apply_W_flat(omega)
    return [-jac[i] - nu4_val * w_om[i] for i in range(grid.N)]


def init_battle1(m=6, nu_4=ZERO_Q, prime_perturbation=ZERO_Q, omega_0=None):
    """Initialize unified Battle 1 state container over pure `fractions.Fraction`."""
    grid = get_grid_2d(m)
    nu_q = Fraction(nu_4)
    if omega_0 is None:
        om0 = make_battle1_initial_omega(m=m, prime_perturbation=prime_perturbation)
    else:
        inv_n = Fraction(1, grid.N)
        m_val = sum(omega_0, ZERO_Q) * inv_n
        om0 = [Fraction(omega_0[i]) - m_val for i in range(grid.N)]

    if nu_q == ZERO_Q:
        n_free = grid.N - 4
        return {
            "m": m,
            "N": grid.N,
            "alpha": list(om0),
            "w_a": [ZERO_Q for i in range(n_free)],
            "w_b": [ZERO_Q for i in range(n_free)],
            "rot_a": [ZERO_Q for k in range(NUM_DOUBLE_CAYLEY_PLANES_2D)],
            "rot_b": [ZERO_Q for k in range(NUM_DOUBLE_CAYLEY_PLANES_2D)],
        }
    return {
        "m": m,
        "N": grid.N,
        "nu_4": nu_q,
        "alpha": list(om0),
        "eta": [ZERO_Q for i in range(grid.N)],
        "amp_scale": ONE_Q,
        "Q_heat": ZERO_Q,
    }


def extract_omega_battle1(state, quantize=True, q_bundle=None):
    """Extract exact rational vorticity field `omega` (`list[Fraction]` of length `N`)."""
    qb = _resolve_b_bundle(quantize, q_bundle)
    grid = get_grid_2d(state["m"])
    if "w_a" in state:
        w_a = state["w_a"]
        w_b = state["w_b"]
        rot_a = state["rot_a"]
        rot_b = state["rot_b"]
        n_w = len(w_a)
        w_mid = [HALF_Q * (w_a[i] + w_b[i]) for i in range(n_w)]
        rot_mid = [HALF_Q * (rot_a[k] + rot_b[k]) for k in range(len(rot_a))]
        anc = _ANC2D_DISPATCH[bool(quantize)](grid, state["alpha"])
        return bezout_vieta_readout_2d(
            grid, state["alpha"], w_mid, rot_mid, anc=anc, quantize=quantize, q_bundle=qb
        )
    return stereographic_readout_2d_energy(
        grid, state["alpha"], state["eta"], amp_scale=state["amp_scale"]
    )


def extract_fields_battle1(state, quantize=True, q_bundle=None):
    """Extract physical fields `{"omega", "psi", "u", "Q_heat"}` strictly in `Q`."""
    qb = _resolve_b_bundle(quantize, q_bundle)
    grid = get_grid_2d(state["m"])
    om = extract_omega_battle1(state, quantize=quantize, q_bundle=qb)
    om_c, d_om = qb["clear_vec"](om)
    psi_c = grid.apply_G0_int(om_c)
    denom_psi = d_om * grid.g0_common_den
    inv_psi = ONE_Q / denom_psi
    psi = [val * inv_psi for val in psi_c]
    inv_u = ONE_Q / (12 * denom_psi)
    u_xy = [
        [val * inv_u for val in grid.apply_Dy_12(psi_c)],
        [-val * inv_u for val in grid.apply_Dx_12(psi_c)],
    ]
    q_heat = state["Q_heat"] if ("Q_heat" in state) else ZERO_Q
    return {
        "omega": om,
        "psi": psi,
        "u": u_xy,
        "Q_heat": q_heat,
    }


def invariants_battle1(state, quantize=True, q_bundle=None):
    """Compute exact rational invariants `(Gamma, E_kin, Z, Q_heat, E_total)` of Battle 1."""
    qb = _resolve_b_bundle(quantize, q_bundle)
    grid = get_grid_2d(state["m"])
    om = extract_omega_battle1(state, quantize=quantize, q_bundle=qb)
    om_c, d_om = qb["clear_vec"](om)
    psi_c = grid.apply_G0_int(om_c)
    inv_psi = ONE_Q / (d_om * grid.g0_common_den)
    psi = [val * inv_psi for val in psi_c]
    circ = sum(om, ZERO_Q)
    ekin = HALF_Q * dot_1d(om, psi)
    enst = HALF_Q * dot_1d(om, om)
    q_heat = state["Q_heat"] if ("Q_heat" in state) else ZERO_Q
    return {
        "Gamma": circ,
        "E_kin": ekin,
        "Z": enst,
        "Q_heat": q_heat,
        "E_total": ekin + q_heat,
    }


def _b1a_tamed_pullback_rate(grid, anc, w_in, rot_in, dt_sq, quantize, qb):
    """Evaluate even-in-`dt` Courant-tamed Bezout-Vieta + Double-Cayley shadow velocity `(dw, drot)`.

    Zero radial `mu_conf` restoring force (`C2-3`), so carrier split stays `O(dt^2)` and
    off-manifold small-`dt` tangents match `X(z(t))` cleanly while bounding bit-lengths to `O(1)`.
    """
    du_num, u_fac, drot_num, r_fac, s_common, y0_val = bezout_vieta_pullback_factored_2d(
        grid, anc, w_in, rot_in, quantize=quantize, q_bundle=qb
    )
    n_w = len(du_num)
    n_planes = len(drot_num)

    inv_a1 = anc["inv_a1"]
    w_norm_rel = dot_1d(w_in, w_in) * inv_a1
    inv_y0_sq = ONE_Q / (y0_val * y0_val)
    tame_den = ONE_Q + dt_sq * (QUARTER_Q * inv_y0_sq + QUARTER_Q * w_norm_rel)
    inv_dw_total = ONE_Q / (u_fac * s_common * tame_den)
    inv_dr_total = ONE_Q / (r_fac * s_common * tame_den)

    dw_eff = [du_num[i] * inv_dw_total for i in range(n_w)]
    drot_eff = [drot_num[k] * inv_dr_total for k in range(n_planes)]
    combined = qb["vec"](dw_eff + drot_eff)
    return combined[0:n_w], combined[n_w : n_w + n_planes]


def step_battle1(state, dt, quantize=True, q_bundle=None):
    """Unified branch-free explicit finite-op `Q` macro-step for Battle 1 (1A & 1B)."""
    qb = _resolve_b_bundle(quantize, q_bundle)
    m = state["m"]
    grid = get_grid_2d(m)
    n_sites = grid.N
    dt_sq = dt * dt

    if "w_a" in state:
        # --- Battle 1A: Reversible Strang Two-Carrier + Pole-Free Double-Cayley ---
        alpha = state["alpha"]
        w_a = state["w_a"]
        w_b = state["w_b"]
        rot_a = state["rot_a"]
        rot_b = state["rot_b"]
        n_w = len(w_a)
        n_planes = len(rot_a)
        half_dt = HALF_Q * dt

        anc = _ANC2D_DISPATCH[bool(quantize)](grid, alpha)

        dw_1, dr_1 = _b1a_tamed_pullback_rate(grid, anc, w_b, rot_b, dt_sq, quantize, qb)
        w_a_half = [w_a[i] + half_dt * dw_1[i] for i in range(n_w)]
        rot_a_half = [rot_a[k] + half_dt * dr_1[k] for k in range(n_planes)]

        dw_2, dr_2 = _b1a_tamed_pullback_rate(grid, anc, w_a_half, rot_a_half, dt_sq, quantize, qb)
        w_b_next = [w_b[i] + dt * dw_2[i] for i in range(n_w)]
        rot_b_next = [rot_b[k] + dt * dr_2[k] for k in range(n_planes)]

        dw_3, dr_3 = _b1a_tamed_pullback_rate(grid, anc, w_b_next, rot_b_next, dt_sq, quantize, qb)
        w_a_next = [w_a_half[i] + half_dt * dw_3[i] for i in range(n_w)]
        rot_a_next = [rot_a_half[k] + half_dt * dr_3[k] for k in range(n_planes)]

        return {
            "m": m,
            "N": n_sites,
            "alpha": alpha,
            "w_a": w_a_next,
            "w_b": w_b_next,
            "rot_a": rot_a_next,
            "rot_b": rot_b_next,
        }

    # --- Battle 1B: Viscous `W_flat` Biharmonic + Exact SOS Thermal Reservoir `Q_heat` ---
    alpha = state["alpha"]
    eta = state["eta"]
    lam = state["amp_scale"]
    nu_4 = state["nu_4"]
    q_heat_cur = state["Q_heat"]

    om_unit = stereographic_readout_2d_energy(grid, alpha, eta, amp_scale=ONE_Q)
    om_u_c, d_omu = qb["clear_vec"](om_unit)
    psi_u_c = grid.apply_G0_int(om_u_c)
    inv_psi_u = ONE_Q / (d_omu * grid.g0_common_den)
    psi_unit = [val * inv_psi_u for val in psi_u_c]
    jac_u_432 = grid.arakawa_bracket_432(psi_u_c, om_u_c)
    inv_jac_u = ONE_Q / (432 * d_omu * d_omu * grid.g0_common_den)
    jac_unit = [val * inv_jac_u for val in jac_u_432]

    two_e_unit = dot_1d(om_unit, psi_unit)
    inv_2e_unit = ONE_Q / two_e_unit
    e_cur = HALF_Q * (lam * lam) * two_e_unit

    w_unit = grid.apply_W_flat(om_unit)
    x1_div_lam = [-lam * jac_unit[i] - nu_4 * w_unit[i] for i in range(n_sites)]
    x1 = [lam * x1_div_lam[i] for i in range(n_sites)]

    diss_unit = dot_1d(psi_unit, w_unit)
    r1 = nu_4 * diss_unit * inv_2e_unit
    deta_1 = stereographic_pullback_2d_energy(grid, alpha, eta, x1_div_lam)

    x1_c, d_x1 = qb["clear_vec"](x1)
    psi_x1_c = grid.apply_G0_int(x1_c)
    jac2_a_432 = grid.arakawa_bracket_432(psi_x1_c, om_u_c)
    jac2_b_432 = grid.arakawa_bracket_432(psi_u_c, x1_c)
    scale_j2 = lam / (432 * d_x1 * d_omu * grid.g0_common_den)
    jac2_sum = [(jac2_a_432[i] + jac2_b_432[i]) * scale_j2 for i in range(n_sites)]

    w_x1 = grid.apply_W_flat(x1)
    ode_acc = [-HALF_Q * (jac2_sum[i] + nu_4 * w_x1[i]) for i in range(n_sites)]

    u_dot = [x1_div_lam[i] + r1 * om_unit[i] for i in range(n_sites)]
    ud_c, d_ud = qb["clear_vec"](u_dot)
    g0_ud_c = grid.apply_G0_int(ud_c)
    inv_g0_ud = ONE_Q / (d_ud * grid.g0_common_den)
    g0_udot = [val * inv_g0_ud for val in g0_ud_c]

    udot_nsq = dot_1d(u_dot, g0_udot)
    lam_d1_neg = lam * r1
    lam_c2 = dot_1d(psi_unit, ode_acc) * inv_2e_unit + HALF_Q * lam * udot_nsq * inv_2e_unit

    inv_lam = ONE_Q / lam
    u_acc_target = [
        (ode_acc[i] - lam_c2 * om_unit[i] + lam_d1_neg * u_dot[i]) * inv_lam
        for i in range(n_sites)
    ]
    chart_acc_u = stereographic_chart_jet2_accel_2d_energy(grid, alpha, eta, deta_1)
    deta_2 = stereographic_pullback_2d_energy(
        grid, alpha, eta, [u_acc_target[i] - chart_acc_u[i] for i in range(n_sites)]
    )

    q2_coeff = r1 * r1 - lam_c2 * inv_lam
    gamma3_coeff = (q2_coeff * q2_coeff) / (2 * r1)
    num_over_dt = r1 + q2_coeff * dt + gamma3_coeff * dt_sq
    den_lam = ONE_Q + dt * num_over_dt
    decay_rate = lam * (num_over_dt / den_lam)

    inv_a_norm = ONE_Q / dot_1d(alpha, alpha)
    deta_rel_sq = dot_1d(deta_1, deta_1) * inv_a_norm
    inv_tame_eta = ONE_Q / (ONE_Q + dt_sq * (QUARTER_Q * deta_rel_sq))
    eta_rate = [(deta_1[i] + dt * deta_2[i]) * inv_tame_eta for i in range(n_sites)]

    q_eta_rate = qb["vec"](eta_rate)
    q_decay = qb["nonneg"](decay_rate)

    eta_next = [eta[i] + dt * q_eta_rate[i] for i in range(n_sites)]
    lam_next = lam - dt * q_decay

    om_next = stereographic_readout_2d_energy(grid, alpha, eta_next, amp_scale=lam_next)
    omn_c, d_omn = qb["clear_vec"](om_next)
    psin_c = grid.apply_G0_int(omn_c)
    e_next = (HALF_Q * _dot_1d_raw(omn_c, psin_c)) / (d_omn * d_omn * grid.g0_common_den)
    q_heat_next = q_heat_cur + (e_cur - e_next)

    return {
        "m": m,
        "N": n_sites,
        "nu_4": nu_4,
        "alpha": alpha,
        "eta": eta_next,
        "amp_scale": lam_next,
        "Q_heat": q_heat_next,
    }


# =============================================================================
# BATTLE 2: 3D Taylor-Green Vortex Stretching on Z_L^3 (L = 6 & L = 8)
# =============================================================================


def make_battle2_initial_velocity(L=6, amp=None, prime_perturbation=ZERO_Q):
    """Construct exact divergence-free zero-mean 3D Taylor-Green Vortex-Stretching IC over `Q`."""
    grid = get_grid_3d(L)
    if L == 6:
        s1 = (ZERO_Q, ONE_Q, ONE_Q, ZERO_Q, -ONE_Q, -ONE_Q)
        c1 = (ONE_Q, HALF_Q, -HALF_Q, -ONE_Q, -HALF_Q, HALF_Q)
        amp_q = Fraction(4, 1) if amp is None else Fraction(amp)
    elif L == 8:
        s1 = (
            ZERO_Q,
            Fraction(5, 7),
            ONE_Q,
            Fraction(5, 7),
            ZERO_Q,
            Fraction(-5, 7),
            -ONE_Q,
            Fraction(-5, 7),
        )
        c1 = (
            ONE_Q,
            Fraction(5, 7),
            ZERO_Q,
            Fraction(-5, 7),
            -ONE_Q,
            Fraction(-5, 7),
            ZERO_Q,
            Fraction(5, 7),
        )
        amp_q = Fraction(3, 1) if amp is None else Fraction(amp)
    else:
        raise ValueError(f"make_battle2_initial_velocity supports L in (6, 8), got L={L}")

    w0 = [[], [], []]
    for idx_i in range(grid.N):
        x_i, y_i, z_i = grid.coords(idx_i)
        ux_i = amp_q * (s1[x_i] * c1[y_i] + s1[x_i] * c1[z_i])
        uy_i = amp_q * (-c1[x_i] * s1[y_i] + s1[y_i] * c1[z_i])
        uz_i = amp_q * (-c1[x_i] * s1[z_i] - c1[y_i] * s1[z_i])
        if prime_perturbation != ZERO_Q and idx_i == 0:
            ux_i = ux_i + Fraction(prime_perturbation)
        w0[0].append(ux_i)
        w0[1].append(uy_i)
        w0[2].append(uz_i)

    return grid.project_leray_zero_mean(w0)


def continuum_rhs_battle2(arg0, L=6, nu=ZERO_Q):
    """Exact 3D spatial semi-discretization RHS `du/dt = P_{div,0}(skew_adv(u) - nu K_Sym u)`.

    Supports both `(u_vec3, L=6, nu=0)` and `(grid, u_vec3, nu=0)` calling conventions.
    """
    if len(arg0) == 0:
        grid = arg0
        u_vec3 = L
        nu_val = nu
    else:
        u_vec3 = arg0
        grid = get_grid_3d(L)
        nu_val = nu
    adv = grid.skew_advection(u_vec3)
    if nu_val == ZERO_Q:
        return grid.project_leray_zero_mean(adv)
    ku = grid.apply_K_sym_vec3(u_vec3)
    n_sites = grid.N
    comb = [
        [adv[c][i] - nu_val * ku[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    return grid.project_leray_zero_mean(comb)


def init_battle2(
    L=6,
    nu=ZERO_Q,
    amp=None,
    prime_perturbation=ZERO_Q,
    u_0=None,
    galilean_boost=None,
):
    """Initialize unified Battle 2 state container over pure `fractions.Fraction`."""
    grid = get_grid_3d(L)
    nu_q = Fraction(nu)
    if u_0 is None:
        u0 = make_battle2_initial_velocity(L=L, amp=amp, prime_perturbation=prime_perturbation)
    else:
        u0 = grid.project_leray(
            [[Fraction(u_0[c][i]) for i in range(grid.N)] for c in range(3)]
        )
    if galilean_boost is not None:
        bx = Fraction(galilean_boost[0])
        by = Fraction(galilean_boost[1])
        bz = Fraction(galilean_boost[2])
        b_tup = (bx, by, bz)
        u0 = [
            [u0[c][i] + b_tup[c] for i in range(grid.N)]
            for c in range(3)
        ]

    if nu_q == ZERO_Q:
        return {
            "L": L,
            "N": grid.N,
            "alpha": u0,
            "eta_a": [[ZERO_Q for i in range(grid.N)] for c in range(3)],
            "eta_b": [[ZERO_Q for i in range(grid.N)] for c in range(3)],
        }
    return {
        "L": L,
        "N": grid.N,
        "nu": nu_q,
        "alpha": u0,
        "eta": [[ZERO_Q for i in range(grid.N)] for c in range(3)],
        "amp_scale": ONE_Q,
        "Q_heat": ZERO_Q,
    }


def _extract_u_anchored_3d(grid, state, mode_flag=True, q_bundle=None):
    """Extract exact `0/1` divergence-free 3D velocity using C `PyLong` anchored readout."""
    qb = _resolve_b_bundle(mode_flag, q_bundle)
    anc3 = _ANC3D_DISPATCH[bool(mode_flag)](grid, state["alpha"])
    n_sites = grid.N
    if "eta_a" in state:
        eta_a = state["eta_a"]
        eta_b = state["eta_b"]
        eta_mid = [
            [HALF_Q * (eta_a[c][i] + eta_b[c][i]) for i in range(n_sites)]
            for c in range(3)
        ]
        return stereographic_readout_anchored_3d(grid, anc3, eta_mid, amp_scale=ONE_Q, q_bundle=qb)
    lam = state["amp_scale"] if ("amp_scale" in state) else ONE_Q
    return stereographic_readout_anchored_3d(grid, anc3, state["eta"], amp_scale=lam, q_bundle=qb)


def extract_velocity_battle2(state, quantize=True, q_bundle=None):
    """Extract exact divergence-free rational velocity `[ux, uy, uz]` (`div_Sym(u) == 0/1`)."""
    grid = get_grid_3d(state["L"])
    return _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)


extract_u_battle2 = extract_velocity_battle2


def extract_fields_battle2(state, quantize=True, q_bundle=None):
    """Extract physical fields `{"u", "omega", "Q_heat"}` strictly in `Q`."""
    grid = get_grid_3d(state["L"])
    u_vec3 = _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)
    om_vec3 = grid.curl_sym(u_vec3)
    q_heat = state["Q_heat"] if ("Q_heat" in state) else ZERO_Q
    return {
        "u": u_vec3,
        "omega": om_vec3,
        "Q_heat": q_heat,
    }


def invariants_battle2(state, quantize=True, q_bundle=None):
    """Compute exact algebraic 3D invariants `(P, div_norm_sq, E_kin, Z, H_hel, Q_heat, E_total)`."""
    grid = get_grid_3d(state["L"])
    u_vec3 = _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)
    om_vec3 = grid.curl_sym(u_vec3)
    div_u = grid.div_sym(u_vec3)
    mom = tuple(sum(u_vec3[c], ZERO_Q) for c in range(3))
    div_norm_sq = dot_1d(div_u, div_u)
    ekin = HALF_Q * dot_vec3(u_vec3, u_vec3)
    enst = HALF_Q * dot_vec3(om_vec3, om_vec3)
    h_hel = dot_vec3(u_vec3, om_vec3)
    q_heat = state["Q_heat"] if ("Q_heat" in state) else ZERO_Q
    return {
        "P": mom,
        "div_norm_sq": div_norm_sq,
        "E_kin": ekin,
        "Z": enst,
        "H_hel": h_hel,
        "Q_heat": q_heat,
        "E_total": ekin + q_heat,
    }


def _skew_advection_frechet_24h_3d(grid, u_base, du_vec):
    """Division-free `24 * h * D(skew_advection)(u_base)[du_vec]` (C `PyLong` on `int`)."""
    n_sites = grid.N
    out = []
    for comp_c in range(3):
        uc = u_base[comp_c]
        duc = du_vec[comp_c]
        acc = [0] * n_sites
        for axis_a in range(3):
            ua = u_base[axis_a]
            dua = du_vec[axis_a]
            d12_uc = grid.apply_D_axis_12h(uc, axis_a)
            d12_duc = grid.apply_D_axis_12h(duc, axis_a)
            prod_d = [dua[i] * uc[i] + ua[i] * duc[i] for i in range(n_sites)]
            d12_prod = grid.apply_D_axis_12h(prod_d, axis_a)
            for i in range(n_sites):
                acc[i] = acc[i] - (dua[i] * d12_uc[i] + ua[i] * d12_duc[i] + d12_prod[i])
        out.append(acc)
    return out


def _pullback_from_tangent_div0_3d(alpha_0, a_sq, inv_a, eta_tan, u_unit_fluct, b_sq, x_div):
    """Exact differential pullback when `eta_tan in V_{div,0} cap alpha_0^perp` and `x_div in V_{div,0}`."""
    n_sites = len(alpha_0[0])
    four_a = 4 * a_sq
    dot_ux = dot_vec3(u_unit_fluct, x_div) * inv_a
    x_tan = [
        [x_div[c][i] - dot_ux * u_unit_fluct[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    scale_main = (four_a + b_sq) * (QUARTER_Q * inv_a)
    dot_ax = dot_vec3(alpha_0, x_tan)
    scale_corr = (HALF_Q * inv_a) * scale_main * dot_ax
    return [
        [
            scale_main * x_tan[c][i] - scale_corr * (2 * alpha_0[c][i] + eta_tan[c][i])
            for i in range(n_sites)
        ]
        for c in range(3)
    ]


def _b2a_tamed_pullback_rate(grid, anc3, eta_in, dt_sq, qb):
    """Evaluate Galilean-covariant even-in-`dt` Courant-tamed stereographic shadow velocity `deta`.

    - Zero radial `mu_conf` restoring force (`C2-3`), keeping `dt * omega_eff < 2` strictly
      inside the elliptic leapfrog stability disk (`|slope| < 0.02 b/step`, small carrier split,
      crisp off-manifold tangents, and `O(dt^2)` accuracy).
    - Branch-free factored C `PyLong` evaluation via `qb` (`< 0.12 s/step` on `L=8` prime).
    """
    n_sites = grid.N
    deta_num, deta_scale_q, dr_snap, d_dr_q, x_div_int, x_tot_scale = (
        stereographic_inviscid_pullback_and_metrics_3d(grid, anc3, eta_in, q_bundle=qb)
    )
    a_sq = anc3["a_sq"]
    rate_rel = dot_vec3(dr_snap, dr_snap) / (4 * a_sq * d_dr_q * d_dr_q)
    adv_rel = dot_vec3(x_div_int, x_div_int) / (a_sq * x_tot_scale * x_tot_scale)
    om_sq = qb["nonneg"](rate_rel + QUARTER_Q * adv_rel)
    inv_eff_scale = ONE_Q / (deta_scale_q * (ONE_Q + dt_sq * om_sq))
    deta_eff = [
        [deta_num[c][i] * inv_eff_scale for i in range(n_sites)]
        for c in range(3)
    ]
    return qb["vec3"](deta_eff)


def _strang_two_carrier_3d(grid, anc3, eta_a, eta_b, dt, qb):
    """Advance 3D `(eta_a, eta_b)` by symmetric Strang (`0/1` reversible for all `dt`)."""
    n_sites = grid.N
    dt_sq = dt * dt
    half_dt = HALF_Q * dt

    v_1 = _b2a_tamed_pullback_rate(grid, anc3, eta_b, dt_sq, qb)
    eta_a_half = [
        [eta_a[c][i] + half_dt * v_1[c][i] for i in range(n_sites)]
        for c in range(3)
    ]

    v_2 = _b2a_tamed_pullback_rate(grid, anc3, eta_a_half, dt_sq, qb)
    eta_b_next = [
        [eta_b[c][i] + dt * v_2[c][i] for i in range(n_sites)]
        for c in range(3)
    ]

    v_3 = _b2a_tamed_pullback_rate(grid, anc3, eta_b_next, dt_sq, qb)
    eta_a_next = [
        [eta_a_half[c][i] + half_dt * v_3[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    return eta_a_next, eta_b_next


def step_battle2(state, dt, quantize=True, q_bundle=None):
    """Unified branch-free explicit finite-op `Q` macro-step for Battle 2 (2A & 2B)."""
    qb = _resolve_b_bundle(quantize, q_bundle)
    L = state["L"]
    grid = get_grid_3d(L)
    n_sites = grid.N
    alpha = state["alpha"]
    anc3 = _ANC3D_DISPATCH[bool(quantize)](grid, alpha)

    if "nu" not in state:
        # --- Battle 2A: Reversible Strang Two-Carrier on V_div Sphere ---
        if "eta_a" in state:
            eta_a_next, eta_b_next = _strang_two_carrier_3d(
                grid, anc3, state["eta_a"], state["eta_b"], dt, qb
            )
            return {
                "L": L,
                "N": n_sites,
                "alpha": alpha,
                "eta_a": eta_a_next,
                "eta_b": eta_b_next,
            }
        eta_a_next, eta_b_next = _strang_two_carrier_3d(
            grid, anc3, state["eta"], state["eta"], dt, qb
        )
        return {
            "L": L,
            "N": n_sites,
            "alpha": alpha,
            "eta": [
                [HALF_Q * (eta_a_next[c][i] + eta_b_next[c][i]) for i in range(n_sites)]
                for c in range(3)
            ],
        }

    # --- Battle 2B: Unified Galilean-Covariant Viscous Navier-Stokes + SOS `Q_heat` ---
    nu = state["nu"]
    eta = state["eta"]
    lam = state["amp_scale"]
    q_heat_cur = state["Q_heat"]
    dt_sq = dt * dt
    inv_lam = ONE_Q / lam

    u_mean = anc3["u_mean"]
    alpha_0 = anc3["alpha_0"]
    a_sq = anc3["a_sq"]
    inv_a = anc3["inv_a"]
    u_mean_int = anc3["u_mean_int"]
    alpha_0_int = anc3["alpha_0_int"]
    a_sq_int = anc3["a_sq_int"]
    d_a = anc3["d_a_raw"]

    eta_i, d_eta_q = qb["clear_vec3"](eta)
    d_eta = qb["to_raw_scale"](d_eta_q)
    ed_int, s_ed_q = grid.project_leray_zero_mean_int(eta_i)
    s_ed = qb["to_raw_scale"](s_ed_q)
    dot_ae_s = _dot_vec3_raw(alpha_0_int, ed_int)
    eta_s = [
        [a_sq_int * ed_int[c][i] - dot_ae_s * alpha_0_int[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    k_scale = a_sq_int * d_eta * s_ed
    inv_k = ONE_Q / k_scale
    eta_tan = [[eta_s[c][i] * inv_k for i in range(n_sites)] for c in range(3)]

    k_sq = k_scale * k_scale
    da_sq = d_a * d_a
    b_s = da_sq * _dot_vec3_raw(eta_s, eta_s)
    b_sq = (ONE_Q * b_s) / (da_sq * k_sq)
    fa_s = 4 * a_sq_int * k_sq
    d_s = fa_s + b_s
    ca_s = fa_s - b_s
    ce_s = 4 * a_sq_int * d_a * k_scale

    u_unit_fluct_int = [
        [ca_s * alpha_0_int[c][i] + ce_s * eta_s[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    u_unit_scale = d_a * d_s
    inv_u_unit_scale = ONE_Q / u_unit_scale
    u_unit_fluct = [
        [u_unit_fluct_int[c][i] * inv_u_unit_scale for i in range(n_sites)]
        for c in range(3)
    ]

    lam_c_list, d_lam_q = qb["clear_vec"]([lam])
    lam_c = lam_c_list[0]
    d_lam = qb["to_raw_scale"](d_lam_q)
    u_cur_int = [
        [
            d_lam * d_s * u_mean_int[c] + lam_c * u_unit_fluct_int[c][i]
            for i in range(n_sites)
        ]
        for c in range(3)
    ]
    u_cur_scale = d_lam * u_unit_scale
    e_cur = HALF_Q * n_sites * dot_1d(u_mean, u_mean) + HALF_Q * (lam * lam) * a_sq

    inv_2e_unit = inv_a

    adv_24h = grid.skew_advection_24h(u_cur_int)
    inv_adv = ONE_Q / (24 * grid.h * u_cur_scale * u_cur_scale)
    ku_12h2 = grid.apply_K_sym_12h2_vec3(u_unit_fluct_int)
    inv_ku = ONE_Q / (12 * (grid.h * grid.h) * u_unit_scale)

    rhs1_raw = [
        [
            (inv_lam * inv_adv) * adv_24h[c][i] - (nu * inv_ku) * ku_12h2[c][i]
            for i in range(n_sites)
        ]
        for c in range(3)
    ]
    rhs1_c, d_r1 = qb["clear_vec3"](rhs1_raw)
    x1_dl_int, s_x1 = grid.project_leray_zero_mean_int(rhs1_c)
    inv_x1_dl = ONE_Q / (d_r1 * s_x1)
    x1_div_lam = [
        [x1_dl_int[c][i] * inv_x1_dl for i in range(n_sites)]
        for c in range(3)
    ]
    x1 = [[lam * x1_div_lam[c][i] for i in range(n_sites)] for c in range(3)]

    diss_unit = dot_vec3(u_unit_fluct_int, ku_12h2) * inv_u_unit_scale * inv_ku
    r1 = nu * diss_unit * inv_2e_unit
    deta_1 = _pullback_from_tangent_div0_3d(
        alpha_0, a_sq, inv_a, eta_tan, u_unit_fluct, b_sq, x1_div_lam
    )

    x1_c, d_x1_v = qb["clear_vec3"](x1)
    d_adv_24h = _skew_advection_frechet_24h_3d(grid, u_cur_int, x1_c)
    inv_dadv = ONE_Q / (24 * grid.h * u_cur_scale * d_x1_v)
    ku_x1_12h2 = grid.apply_K_sym_12h2_vec3(x1_c)
    inv_kux1 = ONE_Q / (12 * (grid.h * grid.h) * d_x1_v)

    rhs2_raw = [
        [
            HALF_Q * (inv_dadv * d_adv_24h[c][i] - (nu * inv_kux1) * ku_x1_12h2[c][i])
            for i in range(n_sites)
        ]
        for c in range(3)
    ]
    rhs2_c, d_r2 = qb["clear_vec3"](rhs2_raw)
    ode_acc_int, s_oa = grid.project_leray_zero_mean_int(rhs2_c)
    inv_oa = ONE_Q / (d_r2 * s_oa)
    ode_acc = [
        [ode_acc_int[c][i] * inv_oa for i in range(n_sites)]
        for c in range(3)
    ]

    u_dot = [
        [x1_div_lam[c][i] + r1 * u_unit_fluct[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    udot_nsq = dot_vec3(u_dot, u_dot)
    lam_d1_neg = lam * r1
    lam_c2 = (
        dot_vec3(u_unit_fluct, ode_acc) * inv_2e_unit
        + HALF_Q * lam * udot_nsq * inv_2e_unit
    )
    u_acc_target = [
        [
            (ode_acc[c][i] - lam_c2 * u_unit_fluct[c][i] + lam_d1_neg * u_dot[c][i]) * inv_lam
            for i in range(n_sites)
        ]
        for c in range(3)
    ]
    chart_acc = stereographic_chart_jet2_accel_3d(alpha, eta_tan, deta_1)
    deta_2 = _pullback_from_tangent_div0_3d(
        alpha_0,
        a_sq,
        inv_a,
        eta_tan,
        u_unit_fluct,
        b_sq,
        [[u_acc_target[c][i] - chart_acc[c][i] for i in range(n_sites)] for c in range(3)],
    )

    q2_coeff = r1 * r1 - lam_c2 * inv_lam
    gamma3_coeff = (q2_coeff * q2_coeff) / (2 * r1)
    num_over_dt = r1 + q2_coeff * dt + gamma3_coeff * dt_sq
    den_lam = ONE_Q + dt * num_over_dt
    decay_rate = lam * (num_over_dt / den_lam)

    inv_4a = QUARTER_Q * inv_a
    rate_rel = dot_vec3(deta_1, deta_1) * inv_4a
    inv_tame_eta = ONE_Q / (ONE_Q + dt_sq * (HALF_Q * rate_rel))
    eta_rate = [
        [(deta_1[c][i] + dt * deta_2[c][i]) * inv_tame_eta for i in range(n_sites)]
        for c in range(3)
    ]

    q_eta_rate = qb["vec3"](eta_rate)
    q_decay = qb["nonneg"](decay_rate)
    eta_next = [
        [eta[c][i] + dt * q_eta_rate[c][i] for i in range(n_sites)]
        for c in range(3)
    ]
    lam_next = lam - dt * q_decay
    e_next = HALF_Q * n_sites * dot_1d(u_mean, u_mean) + HALF_Q * (lam_next * lam_next) * a_sq
    q_heat_next = q_heat_cur + (e_cur - e_next)

    return {
        "L": L,
        "N": n_sites,
        "nu": nu,
        "alpha": alpha,
        "eta": eta_next,
        "amp_scale": lam_next,
        "Q_heat": q_heat_next,
    }


# =============================================================================
# BATTLE 3: Extreme 3D Anti-Parallel Vortex Collision + [2/2] Pade SOS Density
# =============================================================================


def make_battle3_initial_fields(L=6, amp=Fraction(10, 1), prime_perturbation=ZERO_Q):
    """Construct extreme colliding 3D anti-parallel vortex tubes + sharp `rho_0 > 0/1` blob over `Q`."""
    grid = get_grid_3d(L)
    if L == 6:
        s1 = (ZERO_Q, ONE_Q, ONE_Q, ZERO_Q, -ONE_Q, -ONE_Q)
        c1 = (ONE_Q, HALF_Q, -HALF_Q, -ONE_Q, -HALF_Q, HALF_Q)
    elif L == 8:
        s1 = (
            ZERO_Q,
            Fraction(5, 7),
            ONE_Q,
            Fraction(5, 7),
            ZERO_Q,
            Fraction(-5, 7),
            -ONE_Q,
            Fraction(-5, 7),
        )
        c1 = (
            ONE_Q,
            Fraction(5, 7),
            ZERO_Q,
            Fraction(-5, 7),
            -ONE_Q,
            Fraction(-5, 7),
            ZERO_Q,
            Fraction(5, 7),
        )
    else:
        raise ValueError(f"make_battle3_initial_fields supports L in (6, 8), got L={L}")

    amp_q = Fraction(amp)
    w0 = [[], [], []]
    rho0 = []
    for idx_i in range(grid.N):
        x_i, y_i, z_i = grid.coords(idx_i)
        tube = s1[y_i] * (ONE_Q + c1[x_i])
        ux_i = amp_q * (tube * c1[z_i] + HALF_Q * s1[x_i] * c1[y_i])
        uy_i = amp_q * (-c1[x_i] * s1[y_i] * c1[z_i] + HALF_Q * s1[y_i] * c1[z_i])
        uz_i = amp_q * (-s1[x_i] * s1[z_i] + Fraction(3, 4) * c1[y_i] * s1[z_i])
        if prime_perturbation != ZERO_Q and idx_i == 0:
            ux_i = ux_i + Fraction(prime_perturbation)
        w0[0].append(ux_i)
        w0[1].append(uy_i)
        w0[2].append(uz_i)

        if (x_i, y_i, z_i) in ((1, 1, 1), (1, 2, 1), (2, 1, 1), (1, 1, 2)):
            r_val = Fraction(10, 1)
        else:
            r_val = Fraction(1, 20)
        if prime_perturbation != ZERO_Q and idx_i == 0:
            r_val = r_val + Fraction(prime_perturbation)
        rho0.append(r_val)

    u0 = grid.project_leray_zero_mean(w0)
    return u0, rho0


def continuum_rhs_battle3(
    arg0,
    arg1,
    L=6,
    kappa=Fraction(1, 50),
    centered_density=False,
):
    """Exact semi-discrete continuum RHS `(du/dt, drho/dt)` for Battle 3 at `1:1` physical time.

    Supports both `(u_vec3, rho_field, L=6, kappa=1/50, centered_density=False)` and
    `(grid, u_vec3, rho_field, kappa=1/50, centered_density=False)` calling conventions.
    """
    if len(arg0) == 0:
        grid = arg0
        u_vec3 = arg1
        rho_field = L
        kappa_val = kappa
    else:
        u_vec3 = arg0
        rho_field = arg1
        grid = get_grid_3d(L)
        kappa_val = kappa

    n_sites = grid.N
    du = continuum_rhs_battle2(grid, u_vec3, ZERO_Q)

    if centered_density:
        div_flux = [ZERO_Q for i in range(n_sites)]
        for axis in range(3):
            ur = [u_vec3[axis][i] * rho_field[i] for i in range(n_sites)]
            d_ur = grid.apply_D_axis(ur, axis)
            for i in range(n_sites):
                div_flux[i] = div_flux[i] + d_ur[i]
        k_rho = grid.apply_K_sym_scalar(rho_field)
        drho = [-div_flux[i] - kappa_val * k_rho[i] for i in range(n_sites)]
        return du, drho

    r_plus, r_minus = varga_sos_departure_rates_3d(
        grid, u_vec3, kappa_val, c_s=DEFAULT_DENSITY_SOUND_SPEED
    )
    drho = [ZERO_Q for i in range(n_sites)]
    for axis in range(3):
        p1_ax = grid.p1[axis]
        m1_ax = grid.m1[axis]
        rp = r_plus[axis]
        rm = r_minus[axis]
        for y_idx in range(n_sites):
            fp = rp[y_idx] * rho_field[y_idx]
            fm = rm[y_idx] * rho_field[y_idx]
            drho[y_idx] = drho[y_idx] - fp - fm
            drho[p1_ax[y_idx]] = drho[p1_ax[y_idx]] + fp
            drho[m1_ax[y_idx]] = drho[m1_ax[y_idx]] + fm
    return du, drho


def init_battle3(
    L=6,
    kappa=Fraction(1, 50),
    amp=Fraction(10, 1),
    prime_perturbation=ZERO_Q,
    u_0=None,
    rho_0=None,
    galilean_boost=None,
):
    """Initialize unified Battle 3 state container over pure `fractions.Fraction`."""
    grid = get_grid_3d(L)
    u_def, rho_def = make_battle3_initial_fields(L=L, amp=amp, prime_perturbation=prime_perturbation)
    if u_0 is None:
        u0 = u_def
    else:
        u0 = grid.project_leray(
            [[Fraction(u_0[c][i]) for i in range(grid.N)] for c in range(3)]
        )
    if galilean_boost is not None:
        bx = Fraction(galilean_boost[0])
        by = Fraction(galilean_boost[1])
        bz = Fraction(galilean_boost[2])
        b_tup = (bx, by, bz)
        u0 = [
            [u0[c][i] + b_tup[c] for i in range(grid.N)]
            for c in range(3)
        ]
    rho_init = rho_def if (rho_0 is None) else [Fraction(val) for val in rho_0]
    return {
        "L": L,
        "N": grid.N,
        "kappa": Fraction(kappa),
        "alpha": u0,
        "eta_a": [[ZERO_Q for i in range(grid.N)] for c in range(3)],
        "eta_b": [[ZERO_Q for i in range(grid.N)] for c in range(3)],
        "rho": rho_init,
    }


def extract_u_battle3(state, quantize=True, q_bundle=None):
    """Extract exact divergence-free rational velocity `[ux, uy, uz]` of Battle 3."""
    grid = get_grid_3d(state["L"])
    return _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)


def extract_rho_battle3(state):
    """Extract exact positive rational density `rho` (`list[Fraction]` of length `N`)."""
    return list(state["rho"])


def extract_fields_battle3(state, quantize=True, q_bundle=None):
    """Extract physical fields `{"u", "omega", "rho"}` strictly in `Q`."""
    grid = get_grid_3d(state["L"])
    u_vec3 = _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)
    om_vec3 = grid.curl_sym(u_vec3)
    return {
        "u": u_vec3,
        "omega": om_vec3,
        "rho": list(state["rho"]),
    }


def invariants_battle3(state, quantize=True, q_bundle=None):
    """Compute exact algebraic Battle 3 invariants `(P, div_norm_sq, E_kin, Z, mass, B_u, B_omega)`."""
    grid = get_grid_3d(state["L"])
    n_sites = grid.N
    u_vec3 = _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)
    om_vec3 = grid.curl_sym(u_vec3)
    div_u = grid.div_sym(u_vec3)
    rho = state["rho"]
    mom = tuple(sum(u_vec3[c], ZERO_Q) for c in range(3))
    zero_mom = (mom[0] == ZERO_Q) and (mom[1] == ZERO_Q) and (mom[2] == ZERO_Q)

    ekin = HALF_Q * dot_vec3(u_vec3, u_vec3)
    enst = HALF_Q * dot_vec3(om_vec3, om_vec3)
    mass = sum(rho, ZERO_Q)
    b_u = bound_B_u(ekin, n_sites, zero_momentum=zero_mom)
    b_om = bound_B_omega(ekin, n_sites)
    return {
        "P": mom,
        "div_norm_sq": dot_1d(div_u, div_u),
        "E_kin": ekin,
        "Z": enst,
        "mass": mass,
        "B_u": b_u,
        "B_omega": b_om,
    }


def pointwise_bounds_battle3(state, quantize=True, q_bundle=None):
    """Concrete rational pointwise diagnostics `(min_rho, max_u_sq, B_u, max_omega_sq, B_omega)` over `Q`."""
    grid = get_grid_3d(state["L"])
    n_sites = grid.N
    u_vec3 = _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=q_bundle)
    om_vec3 = grid.curl_sym(u_vec3)
    rho = state["rho"]
    mom = tuple(sum(u_vec3[c], ZERO_Q) for c in range(3))
    zero_mom = (mom[0] == ZERO_Q) and (mom[1] == ZERO_Q) and (mom[2] == ZERO_Q)
    ekin = HALF_Q * dot_vec3(u_vec3, u_vec3)

    min_rho = rho[0]
    max_u_sq = ZERO_Q
    max_om_sq = ZERO_Q
    for i_idx in range(n_sites):
        if rho[i_idx] < min_rho:
            min_rho = rho[i_idx]
        usq = (
            u_vec3[0][i_idx] * u_vec3[0][i_idx]
            + u_vec3[1][i_idx] * u_vec3[1][i_idx]
            + u_vec3[2][i_idx] * u_vec3[2][i_idx]
        )
        omsq = (
            om_vec3[0][i_idx] * om_vec3[0][i_idx]
            + om_vec3[1][i_idx] * om_vec3[1][i_idx]
            + om_vec3[2][i_idx] * om_vec3[2][i_idx]
        )
        if usq > max_u_sq:
            max_u_sq = usq
        if omsq > max_om_sq:
            max_om_sq = omsq

    return {
        "min_rho": min_rho,
        "max_u_sq": max_u_sq,
        "B_u": bound_B_u(ekin, n_sites, zero_momentum=zero_mom),
        "max_omega_sq": max_om_sq,
        "B_omega": bound_B_omega(ekin, n_sites),
    }


def step_battle3(state, dt, quantize=True, q_bundle=None):
    """Unified branch-free explicit finite-op `Q` macro-step for Battle 3 (`C2-4`).

    Couples the `0/1`-reversible Galilean-covariant Strang velocity stepper (`max ||u||^2 <= B_u`,
    `max ||omega||^2 <= B_omega`) with the unified branch-free `[2/2]` Pade Heun-SOS + Varga
    squared-M-matrix positive conservative density core `varga_density_step_core_3d(..., q_bundle=qb)`.
    """
    qb = _resolve_b_bundle(quantize, q_bundle)
    L = state["L"]
    grid = get_grid_3d(L)
    n_sites = grid.N
    kappa = state["kappa"]
    alpha = state["alpha"]
    rho = state["rho"]
    anc3 = _ANC3D_DISPATCH[bool(quantize)](grid, alpha)

    u_cur = _extract_u_anchored_3d(grid, state, mode_flag=quantize, q_bundle=qb)
    if "eta_a" in state:
        eta_a_next, eta_b_next = _strang_two_carrier_3d(
            grid, anc3, state["eta_a"], state["eta_b"], dt, qb
        )
    else:
        eta_a_next, eta_b_next = _strang_two_carrier_3d(
            grid, anc3, state["eta"], state["eta"], dt, qb
        )

    eta_mid_next = [
        [HALF_Q * (eta_a_next[c][i] + eta_b_next[c][i]) for i in range(n_sites)]
        for c in range(3)
    ]

    # Exact O(dt^2) non-autonomous Heun endpoint velocity `u_pred = u_cur + dt * X_u(u_cur)`:
    u_cur_c, d_uc = qb["clear_vec3"](u_cur)
    adv_uc = grid.skew_advection_24h(u_cur_c)
    du_cur_int, s_duc = grid.project_leray_zero_mean_int(adv_uc)
    coef_end = dt / (24 * grid.h * d_uc * d_uc * s_duc)
    u_pred_adv = [
        [u_cur[c][i] + coef_end * du_cur_int[c][i] for i in range(n_sites)]
        for c in range(3)
    ]

    rho_next = varga_density_step_core_3d(
        grid,
        rho,
        u_cur,
        kappa,
        dt,
        c_s=DEFAULT_DENSITY_SOUND_SPEED,
        q_bundle=qb,
        u_vec3_end=u_pred_adv,
    )

    if "eta_a" in state:
        return {
            "L": L,
            "N": n_sites,
            "kappa": kappa,
            "alpha": alpha,
            "eta_a": eta_a_next,
            "eta_b": eta_b_next,
            "rho": rho_next,
        }
    return {
        "L": L,
        "N": n_sites,
        "kappa": kappa,
        "alpha": alpha,
        "eta": eta_mid_next,
        "rho": rho_next,
    }


# =============================================================================
# PART II: ZERO-HANDICAP CONTINUUM FLOAT64 BASELINES
# =============================================================================

from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
import math
import operator
import contextlib
import signal
import time
import os
import html
import struct
import zlib
import json

# --- Component: float64_integrators.py ---
RHSFunc = Callable[[List[float]], List[float]]


def to_float64_vec(seq: Iterable[Any]) -> List[float]:
    """Convert an exact rational sequence into float64 with zero artificial handicap."""
    return [float(x) for x in seq]


def flatten_rational_fields(*fields: Sequence[Any]) -> List[float]:
    """Flatten one or more 1D or nested sequences of rationals into a single float64 state vector."""
    out: List[float] = []
    for f in fields:
        for item in f:
            if isinstance(item, (list, tuple)):
                for sub in item:
                    out.append(float(sub))
            else:
                out.append(float(item))
    return out


def step_forward_euler(y: List[float], dt: float, rhs: RHSFunc) -> List[float]:
    """Single Forward Euler step: y_{n+1} = y_n + dt * f(y_n)."""
    k1 = rhs(y)
    return [yi + dt * k1i for yi, k1i in zip(y, k1)]


def step_midpoint_rk2(y: List[float], dt: float, rhs: RHSFunc) -> List[float]:
    """Single Explicit Midpoint (RK2) step:
    k1 = f(y_n)
    k2 = f(y_n + (dt/2) * k1)
    y_{n+1} = y_n + dt * k2
    """
    half_dt = 0.5 * dt
    k1 = rhs(y)
    y_mid = [yi + half_dt * k1i for yi, k1i in zip(y, k1)]
    k2 = rhs(y_mid)
    return [yi + dt * k2i for yi, k2i in zip(y, k2)]


def step_classical_rk4(y: List[float], dt: float, rhs: RHSFunc) -> List[float]:
    """Single Classical 4-stage Runge-Kutta (RK4) step."""
    half_dt = 0.5 * dt
    sixth_dt = dt / 6.0
    k1 = rhs(y)
    y2 = [yi + half_dt * k1i for yi, k1i in zip(y, k1)]
    k2 = rhs(y2)
    y3 = [yi + half_dt * k2i for yi, k2i in zip(y, k2)]
    k3 = rhs(y3)
    y4 = [yi + dt * k3i for yi, k3i in zip(y, k3)]
    k4 = rhs(y4)
    return [
        yi + sixth_dt * (k1i + 2.0 * k2i + 2.0 * k3i + k4i)
        for yi, k1i, k2i, k3i, k4i in zip(y, k1, k2, k3, k4)
    ]


INTEGRATOR_TABLE: Dict[str, Callable[[List[float], float, RHSFunc], List[float]]] = {
    "euler": step_forward_euler,
    "forward_euler": step_forward_euler,
    "rk2": step_midpoint_rk2,
    "midpoint": step_midpoint_rk2,
    "midpoint_rk2": step_midpoint_rk2,
    "rk4": step_classical_rk4,
    "classical_rk4": step_classical_rk4,
}


def check_finite_and_bounded(y: Sequence[float], blowup_threshold: float = 1e12) -> Tuple[bool, float]:
    """Return ``(is_healthy, max_abs)``. Detects NaN, Inf, or amplitude > blowup_threshold."""
    max_abs = 0.0
    for val in y:
        if math.isnan(val) or math.isinf(val):
            return False, float("inf")
        av = abs(val)
        if av > max_abs:
            max_abs = av
        if av > blowup_threshold:
            return False, max_abs
    return True, max_abs


def vec_error_norms(y_approx: Sequence[float], y_ref: Sequence[float]) -> Dict[str, float]:
    """Compute l2 (root-mean-square), l2_sum, and linf error norms between two vectors."""
    n = len(y_ref)
    if n == 0:
        return {"l2_rms": 0.0, "l2_sum": 0.0, "linf": 0.0, "rel_l2": 0.0}
    sq_err = 0.0
    sq_ref = 0.0
    linf = 0.0
    for a, r in zip(y_approx, y_ref):
        d = abs(a - r)
        if math.isnan(d) or math.isinf(d):
            return {
                "l2_rms": float("inf"),
                "l2_sum": float("inf"),
                "linf": float("inf"),
                "rel_l2": float("inf"),
            }
        sq_err += d * d
        sq_ref += r * r
        if d > linf:
            linf = d
    l2_sum = math.sqrt(sq_err)
    l2_rms = math.sqrt(sq_err / float(n))
    rel_l2 = math.sqrt(sq_err / sq_ref) if sq_ref > 0.0 else l2_sum
    return {
        "l2_rms": l2_rms,
        "l2_sum": l2_sum,
        "linf": linf,
        "rel_l2": rel_l2,
    }


def run_float64_trajectory(
    y0: Sequence[float],
    dt: float,
    n_steps: int,
    method: str,
    rhs: RHSFunc,
    blowup_threshold: float = 1e12,
    snapshot_steps: Optional[Iterable[int]] = None,
    observer: Optional[Callable[[int, float, List[float]], Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Integrate ``y' = rhs(y)`` from ``y0`` for ``n_steps`` of ``dt`` using ``method``.

    Stops cleanly on blow-up (NaN, Inf, OverflowError, or ``|y_i| > blowup_threshold``)
    and records the exact ``blowup_step``.
    """
    key = method.lower()
    if key not in INTEGRATOR_TABLE:
        raise ValueError(f"Unknown float64 integrator method: {method!r}")
    stepper = INTEGRATOR_TABLE[key]

    snap_set = set(snapshot_steps) if snapshot_steps is not None else set()
    y = list(y0)
    snapshots: Dict[int, List[float]] = {}
    if 0 in snap_set:
        snapshots[0] = list(y)

    history: List[Dict[str, Any]] = []
    if observer is not None:
        obs0 = dict(observer(0, 0.0, y))
        obs0["step"] = 0
        obs0["t"] = 0.0
        history.append(obs0)

    blew_up = False
    blowup_step: Optional[int] = None
    completed_steps = 0

    for step in range(1, n_steps + 1):
        try:
            y_next = stepper(y, dt, rhs)
        except (OverflowError, ZeroDivisionError, ValueError):
            blew_up = True
            blowup_step = step
            break

        healthy, _ = check_finite_and_bounded(y_next, blowup_threshold=blowup_threshold)
        if not healthy:
            blew_up = True
            blowup_step = step
            break

        y = y_next
        completed_steps = step
        t = step * dt
        if step in snap_set:
            snapshots[step] = list(y)
        if observer is not None:
            obs = dict(observer(step, t, y))
            obs["step"] = step
            obs["t"] = t
            history.append(obs)

    return {
        "method": method,
        "dt": dt,
        "n_steps_requested": n_steps,
        "completed_steps": completed_steps,
        "blew_up": blew_up,
        "blowup_step": blowup_step,
        "final_t": completed_steps * dt,
        "final_y": y if not blew_up else None,
        "snapshots": snapshots,
        "history": history,
    }


def run_rk4_reference(
    y0: Sequence[float],
    T_final: float,
    rhs: RHSFunc,
    dt_ref: float = 1.0 / 512.0,
    record_times: Optional[Iterable[float]] = None,
    observer: Optional[Callable[[int, float, List[float]], Dict[str, Any]]] = None,
    observer_stride: int = 1,
) -> Dict[str, Any]:
    """Run fine continuum reference Classical RK4 at ``dt_ref = 1/512`` up to ``T_final``.

    Exact dyadic ``record_times`` (multiples of ``dt_ref``) are captured without interpolation.
    """
    n_steps = int(round(T_final / dt_ref))
    if abs(n_steps * dt_ref - T_final) > 1e-12:
        raise ValueError(f"T_final={T_final} is not an exact integer multiple of dt_ref={dt_ref}")

    time_to_step: Dict[int, float] = {}
    if record_times is not None:
        for t_req in record_times:
            s_idx = int(round(t_req / dt_ref))
            if abs(s_idx * dt_ref - t_req) > 1e-12:
                raise ValueError(f"record_time={t_req} is not an exact integer multiple of dt_ref={dt_ref}")
            time_to_step[s_idx] = t_req
    time_to_step[0] = 0.0
    time_to_step[n_steps] = T_final

    def _strided_observer(step: int, t: float, y_cur: List[float]) -> Dict[str, Any]:
        if observer is not None and (step % observer_stride == 0 or step == n_steps or step in time_to_step):
            return observer(step, t, y_cur)
        return {}

    traj = run_float64_trajectory(
        y0=y0,
        dt=dt_ref,
        n_steps=n_steps,
        method="classical_rk4",
        rhs=rhs,
        snapshot_steps=time_to_step.keys(),
        observer=_strided_observer if observer is not None else None,
    )

    time_snapshots: Dict[float, List[float]] = {}
    for s_idx, t_val in time_to_step.items():
        if s_idx in traj["snapshots"]:
            time_snapshots[t_val] = traj["snapshots"][s_idx]

    filtered_history = [h for h in traj["history"] if len(h) > 2] if observer is not None else []

    return {
        "dt_ref": dt_ref,
        "T_final": T_final,
        "n_steps": n_steps,
        "blew_up": traj["blew_up"],
        "blowup_step": traj["blowup_step"],
        "final_y": traj["final_y"],
        "time_snapshots": time_snapshots,
        "history": filtered_history,
    }


def compute_pairwise_orders(dts: Sequence[float], errors: Sequence[float]) -> List[Optional[float]]:
    """Compute empirical convergence orders log(e_i / e_{i+1}) / log(dt_i / dt_{i+1})."""
    orders: List[Optional[float]] = []
    for i in range(len(dts) - 1):
        e1, e2 = errors[i], errors[i + 1]
        dt1, dt2 = dts[i], dts[i + 1]
        if (
            e1 is None
            or e2 is None
            or not math.isfinite(e1)
            or not math.isfinite(e2)
            or e1 <= 0.0
            or e2 <= 1e-15
            or dt1 <= 0.0
            or dt2 <= 0.0
            or dt1 == dt2
        ):
            orders.append(None)
        else:
            orders.append(math.log(e1 / e2) / math.log(dt1 / dt2))
    return orders


def verify_convergence_slopes(
    y0: Sequence[float],
    T_eval: float,
    rhs: RHSFunc,
    dts: Sequence[float] = (1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0, 1.0 / 64.0),
    dt_ref: float = 1.0 / 512.0,
) -> Dict[str, Any]:
    """Verify exact O(dt^1), O(dt^2), O(dt^4) convergence of Forward Euler, Midpoint/RK2, and RK4
    against the fine reference RK4(dt=dt_ref) at T_eval on the exact same spatial RHS.
    """
    ref = run_rk4_reference(y0=y0, T_final=T_eval, rhs=rhs, dt_ref=dt_ref, record_times=[T_eval])
    if ref["blew_up"] or ref["final_y"] is None:
        raise RuntimeError(f"Reference RK4(dt={dt_ref}) blew up before T_eval={T_eval}")
    y_ref = ref["final_y"]

    results: Dict[str, Any] = {
        "T_eval": T_eval,
        "dt_ref": dt_ref,
        "dts": list(dts),
        "methods": {},
    }

    expected_orders = {
        "forward_euler": 1.0,
        "midpoint_rk2": 2.0,
        "classical_rk4": 4.0,
    }

    for method, expected_p in expected_orders.items():
        l2_errs: List[float] = []
        linf_errs: List[float] = []
        rel_l2_errs: List[float] = []
        used_dts: List[float] = []

        for dt in dts:
            n_steps = int(round(T_eval / dt))
            if abs(n_steps * dt - T_eval) > 1e-12:
                raise ValueError(f"T_eval={T_eval} is not an integer multiple of dt={dt}")
            traj = run_float64_trajectory(y0=y0, dt=dt, n_steps=n_steps, method=method, rhs=rhs)
            used_dts.append(dt)
            if traj["blew_up"] or traj["final_y"] is None:
                l2_errs.append(float("inf"))
                linf_errs.append(float("inf"))
                rel_l2_errs.append(float("inf"))
            else:
                norms = vec_error_norms(traj["final_y"], y_ref)
                l2_errs.append(norms["l2_rms"])
                linf_errs.append(norms["linf"])
                rel_l2_errs.append(norms["rel_l2"])

        orders_l2 = compute_pairwise_orders(used_dts, l2_errs)
        orders_linf = compute_pairwise_orders(used_dts, linf_errs)
        valid_orders = [p for p in orders_l2 if p is not None]
        # Take asymptotic (finest valid pair) order for strict Gate 4 check
        measured_order = valid_orders[-1] if valid_orders else None

        results["methods"][method] = {
            "expected_order": expected_p,
            "dts": used_dts,
            "l2_rms_errors": l2_errs,
            "linf_errors": linf_errs,
            "rel_l2_errors": rel_l2_errs,
            "orders_l2": orders_l2,
            "orders_linf": orders_linf,
            "asymptotic_order_l2": measured_order,
            "order_matches": (
                measured_order is not None and abs(measured_order - expected_p) <= 0.25
            ),
        }

    results["all_passed"] = all(m["order_matches"] for m in results["methods"].values())
    return results


# --- Component: float64_spatial.py ---
def matvec_f64(mat: Sequence[Sequence[float]], vec: Sequence[float]) -> List[float]:
    """Fast pure-Python float64 dense matrix-vector product using map(operator.mul, ...)."""
    mul = operator.mul
    return [sum(map(mul, row, vec)) for row in mat]


def convert_q_matrix_to_f64(q_mat: Sequence[Sequence[Any]]) -> List[List[float]]:
    """Convert an exact rational 2D matrix into a float64 list-of-lists."""
    return [[float(entry) for entry in row] for row in q_mat]


def build_circulant_3d_matrix_f64(kernel_3d: Sequence[float], L: int) -> List[List[float]]:
    """Expand a translation-invariant 3D kernel K[dx, dy, dz] on Z_L^3 (flat index (x*L + y)*L + z)
    into an N x N float64 matrix for fast map(mul, row, v) evaluation (N = L^3 <= 512).
    """
    n = L * L * L
    if len(kernel_3d) != n:
        raise ValueError(f"Expected kernel length {n}, got {len(kernel_3d)}")
    mat: List[List[float]] = []
    for ix in range(L):
        for iy in range(L):
            for iz in range(L):
                row = [0.0] * n
                for jx in range(L):
                    dx = (ix - jx) % L
                    for jy in range(L):
                        dy = (iy - jy) % L
                        base_k = (dx * L + dy) * L
                        base_j = (jx * L + jy) * L
                        for jz in range(L):
                            dz = (iz - jz) % L
                            row[base_j + jz] = float(kernel_3d[base_k + dz])
                mat.append(row)
    return mat


def assert_rhs_equivalence(
    q_state_vec: Sequence[Fraction],
    rhs_q_fn: Callable[[Sequence[Fraction]], Sequence[Fraction]],
    rhs_f64_fn: Callable[[List[float]], List[float]],
    tol: float = 1e-13,
    label: str = "spatial_rhs",
) -> Dict[str, float]:
    """Cross-verify that ``rhs_f64_fn([float(x) for x in q_state_vec])`` agrees with
    ``[float(r) for r in rhs_q_fn(q_state_vec)]`` to machine precision ``tol``.
    """
    q_out = rhs_q_fn(q_state_vec)
    f64_in = [float(x) for x in q_state_vec]
    f64_out = rhs_f64_fn(f64_in)

    if len(q_out) != len(f64_out):
        raise AssertionError(f"{label}: length mismatch {len(q_out)} vs {len(f64_out)}")

    max_abs_diff = 0.0
    max_ref = 0.0
    for rq, rf in zip(q_out, f64_out):
        rq_f = float(rq)
        diff = abs(rq_f - rf)
        if diff > max_abs_diff:
            max_abs_diff = diff
        if abs(rq_f) > max_ref:
            max_ref = abs(rq_f)

    rel_diff = max_abs_diff / max(max_ref, 1.0)
    if rel_diff > tol:
        raise AssertionError(
            f"{label}: float64 RHS disagrees with exact Q continuum RHS! "
            f"max_abs_diff={max_abs_diff:.3e}, rel_diff={rel_diff:.3e} > tol={tol:.1e}"
        )
    return {
        "max_abs_diff": max_abs_diff,
        "max_rhs_amp": max_ref,
        "rel_diff": rel_diff,
    }


class Battle1Float64Spatial:
    """Zero-handicap float64 mirror of Battle 1 spatial semi-discretization on Z[i]/<(1+i)^m>."""

    def __init__(self, m: int = 6, nu_4: Fraction = Fraction(0, 1)) -> None:
        self.m = m
        self.nu_4_q = Fraction(nu_4)
        self.nu_4 = float(self.nu_4_q)
        self.grid = get_grid_2d(m)
        self.N = self.grid.N
        self.xp1 = self.grid.xp1
        self.xm1 = self.grid.xm1
        self.xp2 = self.grid.xp2
        self.xm2 = self.grid.xm2
        self.yp1 = self.grid.yp1
        self.ym1 = self.grid.ym1
        self.yp2 = self.grid.yp2
        self.ym2 = self.grid.ym2
        self.g0_mat = convert_q_matrix_to_f64(self.grid.g0_matrix)
        self.coset_tables = self.grid.coset_tables

    def apply_G0(self, omega: Sequence[float]) -> List[float]:
        return matvec_f64(self.g0_mat, omega)

    def apply_Dx_sym(self, f: Sequence[float]) -> List[float]:
        c1 = 2.0 / 3.0
        c2 = 1.0 / 12.0
        xp1, xm1, xp2, xm2 = self.xp1, self.xm1, self.xp2, self.xm2
        return [
            c1 * (f[xp1[i]] - f[xm1[i]]) - c2 * (f[xp2[i]] - f[xm2[i]])
            for i in range(self.N)
        ]

    def apply_Dy_sym(self, f: Sequence[float]) -> List[float]:
        c1 = 2.0 / 3.0
        c2 = 1.0 / 12.0
        yp1, ym1, yp2, ym2 = self.yp1, self.ym1, self.yp2, self.ym2
        return [
            c1 * (f[yp1[i]] - f[ym1[i]]) - c2 * (f[yp2[i]] - f[ym2[i]])
            for i in range(self.N)
        ]

    def arakawa_bracket(self, psi: Sequence[float], omega: Sequence[float]) -> List[float]:
        n = self.N
        ax_psi = self.apply_Dx_sym(psi)
        by_psi = self.apply_Dy_sym(psi)
        ax_om = self.apply_Dx_sym(omega)
        by_om = self.apply_Dy_sym(omega)

        j1 = [ax_psi[i] * by_om[i] - by_psi[i] * ax_om[i] for i in range(n)]

        d_psi_by_om = self.apply_Dx_sym([psi[i] * by_om[i] for i in range(n)])
        d_psi_ax_om = self.apply_Dy_sym([psi[i] * ax_om[i] for i in range(n)])

        d_ax_psi_om = self.apply_Dy_sym([ax_psi[i] * omega[i] for i in range(n)])
        d_by_psi_om = self.apply_Dx_sym([by_psi[i] * omega[i] for i in range(n)])

        third = 1.0 / 3.0
        return [
            third * (j1[i] + (d_psi_by_om[i] - d_psi_ax_om[i]) + (d_ax_psi_om[i] - d_by_psi_om[i]))
            for i in range(n)
        ]

    def apply_Lx1(self, f: Sequence[float]) -> List[float]:
        xp1, xm1 = self.xp1, self.xm1
        return [2.0 * f[i] - f[xp1[i]] - f[xm1[i]] for i in range(self.N)]

    def apply_Ly1(self, f: Sequence[float]) -> List[float]:
        yp1, ym1 = self.yp1, self.ym1
        return [2.0 * f[i] - f[yp1[i]] - f[ym1[i]] for i in range(self.N)]

    def apply_W_flat(self, f: Sequence[float]) -> List[float]:
        n = self.N
        lx_f = self.apply_Lx1(f)
        ly_f = self.apply_Ly1(f)
        a_f = [lx_f[i] - 0.5 * ly_f[i] for i in range(n)]
        lx_af = self.apply_Lx1(a_f)
        ly_af = self.apply_Ly1(a_f)
        ly2_f = self.apply_Ly1(ly_f)
        return [(lx_af[i] - 0.5 * ly_af[i]) + 0.75 * ly2_f[i] for i in range(n)]

    def rhs(self, omega: List[float]) -> List[float]:
        psi = self.apply_G0(omega)
        jac = self.arakawa_bracket(psi, omega)
        if self.nu_4 == 0.0:
            return [-v for v in jac]
        w_om = self.apply_W_flat(omega)
        nu4 = self.nu_4
        return [-jac[i] - nu4 * w_om[i] for i in range(self.N)]

    def invariants(self, omega: Sequence[float]) -> Dict[str, float]:
        psi = self.apply_G0(omega)
        mul = operator.mul
        ekin = 0.5 * sum(map(mul, omega, psi))
        enst = 0.5 * sum(map(mul, omega, omega))
        gamma = sum(omega)
        return {"Gamma": gamma, "E_kin": ekin, "Z": enst}

    def shell_energies(self, omega: Sequence[float]) -> List[float]:
        n = self.N
        avgs: List[List[float]] = []
        for r_lvl in range(self.m + 1):
            cosets = self.coset_tables[r_lvl]
            inv_card = 1.0 / float(len(cosets[0]))
            avgs.append([inv_card * sum(omega[c_idx] for c_idx in cosets[i]) for i in range(n)])
        shells = [avgs[0]]
        for tilde_r in range(1, self.m + 1):
            fa = avgs[self.m + 1 - tilde_r]
            ca = avgs[self.m - tilde_r]
            shells.append([fa[i] - ca[i] for i in range(n)])
        mul = operator.mul
        return [0.5 * sum(map(mul, sh, self.apply_G0(sh))) for sh in shells]


class Battle2Float64Spatial:
    """Zero-handicap float64 mirror of Battle 2 3D Symanzik spatial semi-discretization on Z_L^3.

    Flat vector convention: ``y`` of length ``3 * N`` storing ``ux + uy + uz``.
    """

    def __init__(self, L: int = 6, nu: Fraction = Fraction(0, 1)) -> None:
        self.L = L
        self.nu_q = Fraction(nu)
        self.nu = float(self.nu_q)
        self.grid = get_grid_3d(L)
        self.N = self.grid.N
        self.h = float(self.grid.h)
        self.p1 = self.grid.p1
        self.m1 = self.grid.m1
        self.p2 = self.grid.p2
        self.m2 = self.grid.m2
        self.half_m = self.grid.half_m
        self.lag_partners_1d = self.grid.lag_partners_1d
        self.w_div_folded = tuple(
            tuple(tuple(float(v) for v in row_z) for row_z in row_yz)
            for row_yz in self.grid.w_div_folded
        )

    def _apply_G_div_separable_f64(self, scalar_field: Sequence[float]) -> List[float]:
        """Fast O(N * (L//2 + 1)^3) reflection-folded 3D Green convolution in float64."""
        L = self.L
        L2 = L * L
        half_m = self.half_m
        partners = self.lag_partners_1d

        fz = [[0.0] * self.N for _ in range(half_m)]
        for z_out in range(L):
            pz = partners[z_out]
            for az in range(half_m):
                tup = pz[az]
                dst_az = fz[az]
                if len(tup) == 1:
                    z0 = tup[0]
                    for xy in range(L2):
                        base = xy * L
                        dst_az[base + z_out] = scalar_field[base + z0]
                else:
                    z0, z1 = tup[0], tup[1]
                    for xy in range(L2):
                        base = xy * L
                        dst_az[base + z_out] = scalar_field[base + z0] + scalar_field[base + z1]

        fyz = [[[0.0] * self.N for _ in range(half_m)] for _ in range(half_m)]
        for ay in range(half_m):
            for az in range(half_m):
                src = fz[az]
                dst = fyz[ay][az]
                for y_out in range(L):
                    tup = partners[y_out][ay]
                    if len(tup) == 1:
                        y0 = tup[0]
                        for x_i in range(L):
                            bx = x_i * L2
                            by_out = bx + y_out * L
                            by0 = bx + y0 * L
                            dst[by_out : by_out + L] = src[by0 : by0 + L]
                    else:
                        y0, y1 = tup[0], tup[1]
                        for x_i in range(L):
                            bx = x_i * L2
                            by_out = bx + y_out * L
                            by0 = bx + y0 * L
                            by1 = bx + y1 * L
                            for z_i in range(L):
                                dst[by_out + z_i] = src[by0 + z_i] + src[by1 + z_i]

        out = [0.0] * self.N
        w_fold = self.w_div_folded
        for ax in range(half_m):
            for x_out in range(L):
                tup = partners[x_out][ax]
                bx_out = x_out * L2
                if len(tup) == 1:
                    bx0 = tup[0] * L2
                    for ay in range(half_m):
                        w_xy = w_fold[ax][ay]
                        for az in range(half_m):
                            w_val = w_xy[az]
                            if w_val == 0.0:
                                continue
                            src = fyz[ay][az]
                            for yz in range(L2):
                                out[bx_out + yz] += w_val * src[bx0 + yz]
                else:
                    bx0 = tup[0] * L2
                    bx1 = tup[1] * L2
                    for ay in range(half_m):
                        w_xy = w_fold[ax][ay]
                        for az in range(half_m):
                            w_val = w_xy[az]
                            if w_val == 0.0:
                                continue
                            src = fyz[ay][az]
                            for yz in range(L2):
                                out[bx_out + yz] += w_val * (src[bx0 + yz] + src[bx1 + yz])
        return out

    def unpack_vec3(self, y: Sequence[float]) -> Tuple[List[float], List[float], List[float]]:
        n = self.N
        return list(y[0:n]), list(y[n : 2 * n]), list(y[2 * n : 3 * n])

    def pack_vec3(self, ux: Sequence[float], uy: Sequence[float], uz: Sequence[float]) -> List[float]:
        return list(ux) + list(uy) + list(uz)

    def apply_D_axis(self, f: Sequence[float], axis: int) -> List[float]:
        c1 = (2.0 / 3.0) / self.h
        c2 = (1.0 / 12.0) / self.h
        p1, m1 = self.p1[axis], self.m1[axis]
        p2, m2 = self.p2[axis], self.m2[axis]
        return [
            c1 * (f[p1[i]] - f[m1[i]]) - c2 * (f[p2[i]] - f[m2[i]])
            for i in range(self.N)
        ]

    def grad_sym(self, phi: Sequence[float]) -> Tuple[List[float], List[float], List[float]]:
        return (
            self.apply_D_axis(phi, 0),
            self.apply_D_axis(phi, 1),
            self.apply_D_axis(phi, 2),
        )

    def div_sym(self, ux: Sequence[float], uy: Sequence[float], uz: Sequence[float]) -> List[float]:
        dx = self.apply_D_axis(ux, 0)
        dy = self.apply_D_axis(uy, 1)
        dz = self.apply_D_axis(uz, 2)
        return [dx[i] + dy[i] + dz[i] for i in range(self.N)]

    def curl_sym(
        self, ux: Sequence[float], uy: Sequence[float], uz: Sequence[float]
    ) -> Tuple[List[float], List[float], List[float]]:
        dy_uz = self.apply_D_axis(uz, 1)
        dz_uy = self.apply_D_axis(uy, 2)
        dz_ux = self.apply_D_axis(ux, 2)
        dx_uz = self.apply_D_axis(uz, 0)
        dx_uy = self.apply_D_axis(uy, 0)
        dy_ux = self.apply_D_axis(ux, 1)
        n = self.N
        return (
            [dy_uz[i] - dz_uy[i] for i in range(n)],
            [dz_ux[i] - dx_uz[i] for i in range(n)],
            [dx_uy[i] - dy_ux[i] for i in range(n)],
        )

    def apply_K_sym_scalar(self, f: Sequence[float]) -> List[float]:
        inv_h2 = 1.0 / (self.h * self.h)
        c1 = (4.0 / 3.0) * inv_h2
        c2 = (1.0 / 12.0) * inv_h2
        p1x, m1x, p2x, m2x = self.p1[0], self.m1[0], self.p2[0], self.m2[0]
        p1y, m1y, p2y, m2y = self.p1[1], self.m1[1], self.p2[1], self.m2[1]
        p1z, m1z, p2z, m2z = self.p1[2], self.m1[2], self.p2[2], self.m2[2]
        out: List[float] = []
        for i in range(self.N):
            fi = f[i]
            six_fi = 6.0 * fi
            l1 = six_fi - (f[p1x[i]] + f[m1x[i]] + f[p1y[i]] + f[m1y[i]] + f[p1z[i]] + f[m1z[i]])
            l2 = six_fi - (f[p2x[i]] + f[m2x[i]] + f[p2y[i]] + f[m2y[i]] + f[p2z[i]] + f[m2z[i]])
            out.append(c1 * l1 - c2 * l2)
        return out

    def project_leray_zero_mean(
        self, wx: Sequence[float], wy: Sequence[float], wz: Sequence[float]
    ) -> Tuple[List[float], List[float], List[float]]:
        div_w = self.div_sym(wx, wy, wz)
        pot = self._apply_G_div_separable_f64(div_w)
        gx, gy, gz = self.grad_sym(pot)
        n = self.N
        inv_n = 1.0 / float(n)
        ux = [wx[i] + gx[i] for i in range(n)]
        uy = [wy[i] + gy[i] for i in range(n)]
        uz = [wz[i] + gz[i] for i in range(n)]
        mx = sum(ux) * inv_n
        my = sum(uy) * inv_n
        mz = sum(uz) * inv_n
        return (
            [v - mx for v in ux],
            [v - my for v in uy],
            [v - mz for v in uz],
        )

    def skew_advection(
        self, ux: Sequence[float], uy: Sequence[float], uz: Sequence[float]
    ) -> Tuple[List[float], List[float], List[float]]:
        n = self.N
        u_list = (ux, uy, uz)
        out: List[List[float]] = []
        for comp_c in range(3):
            uc = u_list[comp_c]
            acc = [0.0] * n
            for axis_a in range(3):
                ua = u_list[axis_a]
                d_uc = self.apply_D_axis(uc, axis_a)
                d_prod = self.apply_D_axis([ua[i] * uc[i] for i in range(n)], axis_a)
                for i in range(n):
                    acc[i] -= 0.5 * (ua[i] * d_uc[i] + d_prod[i])
            out.append(acc)
        return out[0], out[1], out[2]

    def rhs_components(
        self, ux: Sequence[float], uy: Sequence[float], uz: Sequence[float]
    ) -> Tuple[List[float], List[float], List[float]]:
        ax, ay, az = self.skew_advection(ux, uy, uz)
        if self.nu != 0.0:
            nu = self.nu
            kx = self.apply_K_sym_scalar(ux)
            ky = self.apply_K_sym_scalar(uy)
            kz = self.apply_K_sym_scalar(uz)
            n = self.N
            ax = [ax[i] - nu * kx[i] for i in range(n)]
            ay = [ay[i] - nu * ky[i] for i in range(n)]
            az = [az[i] - nu * kz[i] for i in range(n)]
        return self.project_leray_zero_mean(ax, ay, az)

    def rhs(self, y: List[float]) -> List[float]:
        ux, uy, uz = self.unpack_vec3(y)
        dux, duy, duz = self.rhs_components(ux, uy, uz)
        return self.pack_vec3(dux, duy, duz)

    def invariants(self, y: Sequence[float]) -> Dict[str, float]:
        ux, uy, uz = self.unpack_vec3(y)
        ox, oy, oz = self.curl_sym(ux, uy, uz)
        div_u = self.div_sym(ux, uy, uz)
        mul = operator.mul
        ekin = 0.5 * (sum(map(mul, ux, ux)) + sum(map(mul, uy, uy)) + sum(map(mul, uz, uz)))
        enst = 0.5 * (sum(map(mul, ox, ox)) + sum(map(mul, oy, oy)) + sum(map(mul, oz, oz)))
        h_hel = sum(map(mul, ux, ox)) + sum(map(mul, uy, oy)) + sum(map(mul, uz, oz))
        max_u_sq = max(ux[i] ** 2 + uy[i] ** 2 + uz[i] ** 2 for i in range(self.N))
        max_om_sq = max(ox[i] ** 2 + oy[i] ** 2 + oz[i] ** 2 for i in range(self.N))
        return {
            "Px": sum(ux),
            "Py": sum(uy),
            "Pz": sum(uz),
            "max_abs_div": max(abs(v) for v in div_u),
            "E_kin": ekin,
            "Z": enst,
            "H_hel": h_hel,
            "max_u_sq": max_u_sq,
            "max_omega_sq": max_om_sq,
        }


class Battle3Float64Spatial:
    """Zero-handicap float64 mirror of Battle 3 (3D Collision + Density Transport).

    Flat vector convention: ``y`` of length ``4 * N`` storing ``ux + uy + uz + rho``.
    Supports both:
      - ``centered_density=True``: 4th-order centered ``-div_Sym(u * rho) - kappa * K_Sym(rho)``
        (the continuum PDE semi-discretization that exhibits negative Gibbs ringing ``min rho < 0``).
      - ``centered_density=False``: CFL-normalized Varga-SOS continuous-time generator.
    """

    def __init__(
        self,
        L: int = 6,
        kappa: Fraction = Fraction(1, 50),
        centered_density: bool = True,
    ) -> None:
        self.L = L
        self.kappa_q = Fraction(kappa)
        self.kappa = float(self.kappa_q)
        self.centered_density = centered_density
        self.vel_op = Battle2Float64Spatial(L=L, nu=Fraction(0, 1))
        self.N = self.vel_op.N
        self.p1 = self.vel_op.p1
        self.m1 = self.vel_op.m1
        # Cycle 2: un-dilated 1:1 low-diffusion Varga-SOS rates r^{+-} = kappa + (c_s +- u)^2/(4 c_s)
        self.c_s = float(DEFAULT_DENSITY_SOUND_SPEED)

    def unpack_state(self, y: Sequence[float]) -> Tuple[List[float], List[float], List[float], List[float]]:
        n = self.N
        return (
            list(y[0:n]),
            list(y[n : 2 * n]),
            list(y[2 * n : 3 * n]),
            list(y[3 * n : 4 * n]),
        )

    def pack_state(
        self, ux: Sequence[float], uy: Sequence[float], uz: Sequence[float], rho: Sequence[float]
    ) -> List[float]:
        return list(ux) + list(uy) + list(uz) + list(rho)

    def rhs(self, y: List[float]) -> List[float]:
        ux, uy, uz, rho = self.unpack_state(y)
        dux, duy, duz = self.vel_op.rhs_components(ux, uy, uz)
        n = self.N
        kap = self.kappa

        if self.centered_density:
            dx_ur = self.vel_op.apply_D_axis([ux[i] * rho[i] for i in range(n)], 0)
            dy_ur = self.vel_op.apply_D_axis([uy[i] * rho[i] for i in range(n)], 1)
            dz_ur = self.vel_op.apply_D_axis([uz[i] * rho[i] for i in range(n)], 2)
            k_rho = self.vel_op.apply_K_sym_scalar(rho)
            drho = [-(dx_ur[i] + dy_ur[i] + dz_ur[i]) - kap * k_rho[i] for i in range(n)]
            return self.pack_state(dux, duy, duz, drho)

        # Cycle 3: un-dilated 1:1 low-diffusion [2/2] Pade SOS rates
        kap = self.kappa
        inv12 = 1.0 / 12.0
        drho = [0.0] * n
        u_list = (ux, uy, uz)
        for axis in range(3):
            ua = u_list[axis]
            p1_ax = self.p1[axis]
            m1_ax = self.m1[axis]
            for i in range(n):
                u_i = ua[i]
                r_i = rho[i]
                u2 = u_i * u_i
                f4p9 = 4.0 * u2 + 9.0
                sig = f4p9 * inv12
                r0 = (144.0 * u2) / (f4p9 * f4p9)
                r0_2 = r0 * r0
                ten_r0 = 10.0 * r0
                ca = sig * ((5.0 * r0_2 + ten_r0 + 1.0) / (r0_2 + ten_r0 + 5.0))
                fp = (kap + 0.5 * (ca + u_i)) * r_i
                fm = (kap + 0.5 * (ca - u_i)) * r_i
                drho[i] -= fp + fm
                drho[p1_ax[i]] += fp
                drho[m1_ax[i]] += fm
        return self.pack_state(dux, duy, duz, drho)

    def invariants(self, y: Sequence[float]) -> Dict[str, float]:
        n = self.N
        inv = self.vel_op.invariants(y[0 : 3 * n])
        rho = y[3 * n : 4 * n]
        inv["mass"] = sum(rho)
        inv["min_rho"] = min(rho)
        inv["max_rho"] = max(rho)
        return inv


def build_verified_battle1_baseline(
    m: int = 6,
    nu_4: Fraction = Fraction(0, 1),
    prime_perturbation: Fraction = Fraction(0, 1),
) -> Tuple[List[Fraction], List[float], Battle1Float64Spatial, Dict[str, float]]:
    """Build and machine-verify zero-handicap Battle 1 float64 initial vector and RHS."""
    q_om0 = make_battle1_initial_omega(m=m, prime_perturbation=prime_perturbation)
    op = Battle1Float64Spatial(m=m, nu_4=nu_4)
    eq_stats = assert_rhs_equivalence(
        q_state_vec=q_om0,
        rhs_q_fn=lambda vec: continuum_rhs_battle1(list(vec), m=m, nu_4=nu_4),
        rhs_f64_fn=op.rhs,
        label=f"Battle1(m={m}, nu_4={nu_4})",
    )
    return q_om0, [float(x) for x in q_om0], op, eq_stats


def build_verified_battle2_baseline(
    L: int = 6,
    nu: Fraction = Fraction(0, 1),
    amp: Any = None,
    prime_perturbation: Fraction = Fraction(0, 1),
) -> Tuple[List[Fraction], List[float], Battle2Float64Spatial, Dict[str, float]]:
    """Build and machine-verify zero-handicap Battle 2 float64 initial vector and RHS."""
    u0_q = make_battle2_initial_velocity(L=L, amp=amp, prime_perturbation=prime_perturbation)
    flat_q = list(u0_q[0]) + list(u0_q[1]) + list(u0_q[2])
    op = Battle2Float64Spatial(L=L, nu=nu)
    n = op.N

    def _q_rhs_flat(vec: Sequence[Fraction]) -> List[Fraction]:
        uv = [list(vec[0:n]), list(vec[n : 2 * n]), list(vec[2 * n : 3 * n])]
        du = continuum_rhs_battle2(uv, L=L, nu=nu)
        return list(du[0]) + list(du[1]) + list(du[2])

    eq_stats = assert_rhs_equivalence(
        q_state_vec=flat_q,
        rhs_q_fn=_q_rhs_flat,
        rhs_f64_fn=op.rhs,
        label=f"Battle2(L={L}, nu={nu})",
    )
    return flat_q, [float(x) for x in flat_q], op, eq_stats


def build_verified_battle3_baseline(
    L: int = 6,
    kappa: Fraction = Fraction(1, 50),
    amp: Fraction = Fraction(10, 1),
    centered_density: bool = True,
    prime_perturbation: Fraction = Fraction(0, 1),
) -> Tuple[List[Fraction], List[float], Battle3Float64Spatial, Dict[str, float]]:
    """Build and machine-verify zero-handicap Battle 3 float64 initial vector and RHS."""
    u0_q, rho0_q = make_battle3_initial_fields(L=L, amp=amp, prime_perturbation=prime_perturbation)
    flat_q = list(u0_q[0]) + list(u0_q[1]) + list(u0_q[2]) + list(rho0_q)
    op = Battle3Float64Spatial(L=L, kappa=kappa, centered_density=centered_density)
    n = op.N

    def _q_rhs_flat(vec: Sequence[Fraction]) -> List[Fraction]:
        uv = [list(vec[0:n]), list(vec[n : 2 * n]), list(vec[2 * n : 3 * n])]
        rv = list(vec[3 * n : 4 * n])
        du, drho = continuum_rhs_battle3(uv, rv, L=L, kappa=kappa, centered_density=centered_density)
        return list(du[0]) + list(du[1]) + list(du[2]) + list(drho)

    eq_stats = assert_rhs_equivalence(
        q_state_vec=flat_q,
        rhs_q_fn=_q_rhs_flat,
        rhs_f64_fn=op.rhs,
        label=f"Battle3(L={L}, centered={centered_density})",
    )
    return flat_q, [float(x) for x in flat_q], op, eq_stats


# =============================================================================
# PART III: SELF-VERIFYING ASSERTION HARNESS (BATTLES 1, 2, AND 3)
# =============================================================================

# --- Component: guards_and_inspector.py ---
class StepStallError(TimeoutError):
    """Raised when a single macro-step exceeds the strict wall-clock limit."""


@contextlib.contextmanager
def sigalrm_step_guard(timeout_sec: int = 5, label: str = "macro-step") -> Iterator[Dict[str, float]]:
    """Enforce a strict SIGALRM per-step timeout and record wall-clock duration.

    Yields a mutable dict ``timing`` populated with ``timing['elapsed_sec']`` upon exit.
    """
    timing: Dict[str, float] = {"elapsed_sec": 0.0}

    def _alarm_handler(signum: int, frame: Any) -> None:
        raise StepStallError(f"SIGALRM timeout: {label} exceeded {timeout_sec}s wall-clock limit")

    prev_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _alarm_handler)
    t0 = time.perf_counter()
    signal.alarm(timeout_sec)
    try:
        yield timing
    finally:
        signal.alarm(0)
        timing["elapsed_sec"] = time.perf_counter() - t0
        signal.signal(signal.SIGALRM, prev_handler)


def run_guarded_step(
    step_fn: Callable[..., Any],
    *args: Any,
    timeout_sec: int = 5,
    label: str = "macro-step",
    **kwargs: Any,
) -> Tuple[Any, float]:
    """Execute ``step_fn(*args, **kwargs)`` inside a 5s SIGALRM guard.

    Returns ``(result, elapsed_sec)``.
    """
    with sigalrm_step_guard(timeout_sec=timeout_sec, label=label) as timing:
        result = step_fn(*args, **kwargs)
    return result, timing["elapsed_sec"]


def _scalar_bits(val: Any) -> Optional[Tuple[int, int]]:
    """Return ``(num_bits, den_bits)`` if ``val`` is a Fraction or int, else None."""
    if isinstance(val, bool):
        return None
    if isinstance(val, Fraction):
        return abs(val.numerator).bit_length(), abs(val.denominator).bit_length()
    if isinstance(val, int):
        return abs(val).bit_length(), 1
    if hasattr(val, "numerator") and hasattr(val, "denominator"):
        num = val.numerator
        den = val.denominator
        if isinstance(num, int) and isinstance(den, int):
            return abs(num).bit_length(), abs(den).bit_length()
    return None


def _walk_scalars(obj: Any, prefix: str, out: List[Tuple[str, int, int]], visited: set) -> None:
    """Recursively collect ``(slot_key, num_bits, den_bits)`` from nested containers."""
    sb = _scalar_bits(obj)
    if sb is not None:
        out.append((prefix or "root", sb[0], sb[1]))
        return
    if obj is None or isinstance(obj, (str, bytes, float)):
        return

    obj_id = id(obj)
    if obj_id in visited:
        return
    visited.add(obj_id)

    if isinstance(obj, dict):
        for k, v in obj.items():
            sub_prefix = f"{prefix}.{k}" if prefix else str(k)
            _walk_scalars(v, sub_prefix, out, visited)
    elif isinstance(obj, tuple) and hasattr(obj, "_fields"):
        # namedtuple
        for k in obj._fields:
            sub_prefix = f"{prefix}.{k}" if prefix else str(k)
            _walk_scalars(getattr(obj, k), sub_prefix, out, visited)
    elif isinstance(obj, (list, tuple)):
        # Group list/tuple elements under their parent slot name to avoid thousands of index keys,
        # unless at the root where top-level tuple slots distinguish carrier components.
        if not prefix:
            for idx, v in enumerate(obj):
                _walk_scalars(v, f"slot_{idx}", out, visited)
        else:
            for v in obj:
                _walk_scalars(v, prefix, out, visited)
    elif hasattr(obj, "__dict__"):
        for k, v in vars(obj).items():
            if k.startswith("_"):
                continue
            sub_prefix = f"{prefix}.{k}" if prefix else str(k)
            _walk_scalars(v, sub_prefix, out, visited)
    elif hasattr(obj, "__slots__"):
        for k in obj.__slots__:
            if k.startswith("_") or not hasattr(obj, k):
                continue
            sub_prefix = f"{prefix}.{k}" if prefix else str(k)
            _walk_scalars(getattr(obj, k), sub_prefix, out, visited)


def inspect_bits(state: Any, physical_fields: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Recursively inspect numerator and denominator bit-lengths across ``state`` and ``physical_fields``.

    Returns a dictionary containing:
      - ``max_bits``: overall peak bit-length across both state and physical fields
      - ``state_max_bits``: peak bit-length in ``state``
      - ``fields_max_bits``: peak bit-length in ``physical_fields``
      - ``max_num_bits``: overall peak numerator bit-length
      - ``max_den_bits``: overall peak denominator bit-length
      - ``mean_bits``: mean scalar max(num_bits, den_bits)
      - ``scalar_count``: total number of rational scalars scanned
      - ``per_slot_max_bits``: dict mapping slot/field names to peak bit-lengths
    """
    state_entries: List[Tuple[str, int, int]] = []
    _walk_scalars(state, "state", state_entries, set())

    field_entries: List[Tuple[str, int, int]] = []
    if physical_fields is not None:
        _walk_scalars(physical_fields, "phys", field_entries, set())

    all_entries = state_entries + field_entries
    if not all_entries:
        return {
            "max_bits": 0,
            "state_max_bits": 0,
            "fields_max_bits": 0,
            "max_num_bits": 0,
            "max_den_bits": 0,
            "mean_bits": 0.0,
            "scalar_count": 0,
            "per_slot_max_bits": {},
        }

    per_slot: Dict[str, int] = {}
    max_num = 0
    max_den = 0
    total_b = 0

    for slot, nb, db in all_entries:
        b = nb if nb >= db else db
        if nb > max_num:
            max_num = nb
        if db > max_den:
            max_den = db
        total_b += b
        prev = per_slot.get(slot, 0)
        if b > prev:
            per_slot[slot] = b

    state_max = max((nb if nb >= db else db for _, nb, db in state_entries), default=0)
    fields_max = max((nb if nb >= db else db for _, nb, db in field_entries), default=0)
    overall_max = max_num if max_num >= max_den else max_den

    return {
        "max_bits": overall_max,
        "state_max_bits": state_max,
        "fields_max_bits": fields_max,
        "max_num_bits": max_num,
        "max_den_bits": max_den,
        "mean_bits": total_b / len(all_entries),
        "scalar_count": len(all_entries),
        "per_slot_max_bits": per_slot,
    }


def check_bit_plateau(
    max_bits_history: List[int],
    bit_cap: int = 10000,
    max_late_slope_bits_per_step: float = 2.0,
) -> Dict[str, Any]:
    """Verify both the hard <= 10,000-bit bound and O(1) late-time plateau (no linear drift).

    Computes ordinary least-squares slope over the second half of the trajectory ``[n//2:]``.
    """
    if not max_bits_history:
        return {"passed": True, "peak_bits": 0, "late_slope": 0.0, "within_cap": True, "is_plateau": True}

    peak_bits = max(max_bits_history)
    within_cap = peak_bits <= bit_cap

    n = len(max_bits_history)
    start = max(1, n // 2)
    tail = max_bits_history[start:]
    m = len(tail)
    if m >= 4:
        mean_x = (m - 1) / 2.0
        mean_y = sum(tail) / float(m)
        num = sum((i - mean_x) * (y - mean_y) for i, y in enumerate(tail))
        den = sum((i - mean_x) ** 2 for i in range(m))
        late_slope = num / den if den > 0.0 else 0.0
    else:
        late_slope = 0.0

    is_plateau = late_slope <= max_late_slope_bits_per_step
    return {
        "passed": within_cap and is_plateau,
        "peak_bits": peak_bits,
        "within_cap": within_cap,
        "late_slope": late_slope,
        "is_plateau": is_plateau,
        "final_bits": max_bits_history[-1],
        "mid_bits": max_bits_history[n // 2],
    }


# --- Component: diagnostics.py ---
ZERO = Fraction(0, 1)
HALF = Fraction(1, 2)


def q_sum(vals: Iterable[Fraction]) -> Fraction:
    """Exact rational sum over Fraction scalars."""
    acc = ZERO
    for v in vals:
        acc += v
    return acc


def q_dot_scalar(a: Sequence[Fraction], b: Sequence[Fraction]) -> Fraction:
    """Exact rational dot product <a, b> = sum_i a_i * b_i."""
    acc = ZERO
    for ai, bi in zip(a, b):
        acc += ai * bi
    return acc


def normalize_vec3_field(u: Any, n_sites: Optional[int] = None) -> Tuple[List[Fraction], List[Fraction], List[Fraction]]:
    """Normalize a 3D vector field into ``(ux, uy, uz)`` component lists of length ``N``.

    Supports:
      - ``(ux, uy, uz)`` tuple/list of 3 sequences each of length ``N`` (when ``N != 3`` or nested sequences)
      - ``[(ux_0, uy_0, uz_0), ...]`` sequence of ``N`` 3-tuples
      - ``{"ux": ..., "uy": ..., "uz": ...}`` or ``{"x": ..., "y": ..., "z": ...}`` dict
      - flat sequence ``[ux_0, uy_0, uz_0, ux_1, ...]`` of length ``3 * N``
    """
    if isinstance(u, dict):
        if "ux" in u and "uy" in u and "uz" in u:
            return list(u["ux"]), list(u["uy"]), list(u["uz"])
        if "x" in u and "y" in u and "z" in u:
            return list(u["x"]), list(u["y"]), list(u["z"])
        raise ValueError(f"Unsupported dict keys for 3D vector field: {list(u.keys())}")

    if isinstance(u, (list, tuple)) and len(u) == 3 and isinstance(u[0], (list, tuple)):
        if n_sites is None or len(u[0]) == n_sites:
            if len(u[0]) != 3 or n_sites == 3:
                return list(u[0]), list(u[1]), list(u[2])

    if isinstance(u, (list, tuple)) and len(u) > 0 and isinstance(u[0], (list, tuple)) and len(u[0]) == 3:
        ux = [ pt[0] for pt in u ]
        uy = [ pt[1] for pt in u ]
        uz = [ pt[2] for pt in u ]
        return ux, uy, uz

    if isinstance(u, (list, tuple)):
        if len(u) % 3 != 0:
            raise ValueError(f"Flat 3D vector field length {len(u)} is not divisible by 3")
        ux = list(u[0::3])
        uy = list(u[1::3])
        uz = list(u[2::3])
        return ux, uy, uz

    raise TypeError(f"Cannot normalize 3D vector field of type {type(u)!r}")


def q_dot_vec3(u: Any, v: Any, n_sites: Optional[int] = None) -> Fraction:
    """Exact rational inner product <u, v> = sum_x (u_x v_x + u_y v_y + u_z v_z)."""
    ux, uy, uz = normalize_vec3_field(u, n_sites=n_sites)
    vx, vy, vz = normalize_vec3_field(v, n_sites=n_sites)
    return q_dot_scalar(ux, vx) + q_dot_scalar(uy, vy) + q_dot_scalar(uz, vz)


def q_kinetic_energy_vec3(u: Any, n_sites: Optional[int] = None) -> Fraction:
    """Exact 3D kinetic energy E_kin = (1/2) <u, u>."""
    return HALF * q_dot_vec3(u, u, n_sites=n_sites)


def q_enstrophy_vec3(omega: Any, n_sites: Optional[int] = None) -> Fraction:
    """Exact 3D enstrophy Z = (1/2) <omega, omega>."""
    return HALF * q_dot_vec3(omega, omega, n_sites=n_sites)


def q_helicity_vec3(u: Any, omega: Any, n_sites: Optional[int] = None) -> Fraction:
    """Exact discrete Moffatt helicity H_hel = <u, omega>."""
    return q_dot_vec3(u, omega, n_sites=n_sites)


def q_momentum_vec3(u: Any, n_sites: Optional[int] = None) -> Tuple[Fraction, Fraction, Fraction]:
    """Exact 3D total linear momentum (P_x, P_y, P_z) = (sum u_x, sum u_y, sum u_z)."""
    ux, uy, uz = normalize_vec3_field(u, n_sites=n_sites)
    return q_sum(ux), q_sum(uy), q_sum(uz)


def q_max_norm_sq_vec3(u: Any, n_sites: Optional[int] = None) -> Fraction:
    """Exact sitewise peak squared norm max_x ||u(x)||^2 in Q (zero square roots)."""
    ux, uy, uz = normalize_vec3_field(u, n_sites=n_sites)
    best = ZERO
    for a, b, c in zip(ux, uy, uz):
        nsq = a * a + b * b + c * c
        if nsq > best:
            best = nsq
    return best


def q_max_abs(vals: Iterable[Fraction]) -> Fraction:
    """Exact peak absolute value max_x |vals(x)| in Q."""
    best = ZERO
    for v in vals:
        av = -v if v < ZERO else v
        if av > best:
            best = av
    return best


def q_min_val(vals: Sequence[Fraction]) -> Fraction:
    """Exact minimum value min_x vals(x) in Q."""
    if not vals:
        raise ValueError("q_min_val called on empty sequence")
    best = vals[0]
    for v in vals[1:]:
        if v < best:
            best = v
    return best


def max_state_diff_exact(a: Any, b: Any) -> Fraction:
    """Recursively compute max |a_i - b_i| across nested Fraction/int state trees."""
    if isinstance(a, (Fraction, int)) and not isinstance(a, bool) and isinstance(b, (Fraction, int)) and not isinstance(b, bool):
        diff = Fraction(a) - Fraction(b)
        return -diff if diff < ZERO else diff

    if type(a) is not type(b) and not (isinstance(a, (list, tuple)) and isinstance(b, (list, tuple))):
        raise TypeError(f"State structure mismatch: {type(a)!r} vs {type(b)!r}")

    if isinstance(a, dict):
        if set(a.keys()) != set(b.keys()):
            raise ValueError(f"State dict key mismatch: {set(a.keys())} vs {set(b.keys())}")
        best = ZERO
        for k in a:
            d = max_state_diff_exact(a[k], b[k])
            if d > best:
                best = d
        return best

    if isinstance(a, tuple) and hasattr(a, "_fields"):
        best = ZERO
        for k in a._fields:
            d = max_state_diff_exact(getattr(a, k), getattr(b, k))
            if d > best:
                best = d
        return best

    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            raise ValueError(f"State sequence length mismatch: {len(a)} vs {len(b)}")
        best = ZERO
        for ai, bi in zip(a, b):
            d = max_state_diff_exact(ai, bi)
            if d > best:
                best = d
        return best

    if hasattr(a, "__dict__") and hasattr(b, "__dict__"):
        return max_state_diff_exact(vars(a), vars(b))

    if a == b:
        return ZERO
    raise TypeError(f"Cannot compare leaf objects of type {type(a)!r}")


def check_exact_reversibility(
    step_fn: Callable[[Any, Fraction], Any],
    state0: Any,
    dt: Fraction,
    n_steps: int = 1,
    step_guard_fn: Optional[Callable[..., Tuple[Any, float]]] = None,
) -> Dict[str, Any]:
    """Verify exact 0/1 reversibility: F_{-dt}^n ( F_{+dt}^n (state0) ) == state0."""
    cur = state0
    max_step_sec = 0.0
    for i in range(n_steps):
        if step_guard_fn is not None:
            cur, sec = step_guard_fn(step_fn, cur, dt, label=f"rev_fwd_{i+1}")
            if sec > max_step_sec:
                max_step_sec = sec
        else:
            cur = step_fn(cur, dt)

    mid_diff = max_state_diff_exact(cur, state0)

    for i in range(n_steps):
        if step_guard_fn is not None:
            cur, sec = step_guard_fn(step_fn, cur, -dt, label=f"rev_bwd_{i+1}")
            if sec > max_step_sec:
                max_step_sec = sec
        else:
            cur = step_fn(cur, -dt)

    rev_residual = max_state_diff_exact(cur, state0)
    return {
        "n_steps": n_steps,
        "dt": str(dt),
        "exact_zero_residual": (rev_residual == ZERO),
        "max_abs_residual": str(rev_residual),
        "forward_displacement_nonzero": (mid_diff > ZERO),
        "forward_max_displacement_float": float(mid_diff),
        "max_step_sec": max_step_sec,
    }


def check_small_dt_tangent(
    step_fn: Callable[[Any, Fraction], Any],
    extract_vec_fn: Callable[[Any], List[Fraction]],
    exact_rhs_vec_fn: Callable[[Any], List[Fraction]],
    state0: Any,
    exponents: Sequence[int] = (4, 10, 20, 30),
) -> Dict[str, Any]:
    """Verify (F_dt(z0) - z0)/dt -> X(z0) smoothly down to dt = 2^-30 with zero freezing."""
    z0 = extract_vec_fn(state0)
    rhs0 = exact_rhs_vec_fn(state0)
    rhs_linf = float(q_max_abs(rhs0))

    dts_float: List[float] = []
    tangent_errs: List[float] = []
    records: List[Dict[str, Any]] = []

    for p in exponents:
        dt = Fraction(1, 1 << p)
        s1 = step_fn(state0, dt)
        z1 = extract_vec_fn(s1)
        inv_dt = Fraction(1 << p, 1)
        diff_quot = [(z1i - z0i) * inv_dt for z1i, z0i in zip(z1, z0)]
        err_linf_q = q_max_abs([dq - r for dq, r in zip(diff_quot, rhs0)])
        err_linf = float(err_linf_q)
        dq_linf = float(q_max_abs(diff_quot))

        dts_float.append(float(dt))
        tangent_errs.append(err_linf)
        records.append({
            "exponent": p,
            "dt": f"2^-{p}",
            "diff_quot_linf": dq_linf,
            "rhs_linf": rhs_linf,
            "tangent_error_linf": err_linf,
        })

    # Compute pairwise slopes of tangent error vs dt
    slopes: List[Optional[float]] = []
    for i in range(len(dts_float) - 1):
        e1, e2 = tangent_errs[i], tangent_errs[i + 1]
        if e1 > 0.0 and e2 > 0.0:
            slopes.append(math.log(e1 / e2) / math.log(dts_float[i] / dts_float[i + 1]))
        elif e2 == 0.0:
            slopes.append(float("inf"))
        else:
            slopes.append(None)

    non_frozen = tangent_errs[-1] < 1e-4 * max(rhs_linf, 1e-12)
    return {
        "rhs_linf": rhs_linf,
        "records": records,
        "slopes": slopes,
        "non_frozen_at_2_minus_30": non_frozen,
    }


def make_zero_mean_prime_perturbation(n_sites: int, amplitude: Fraction = Fraction(1, 1)) -> List[Fraction]:
    """Generate a deterministic ugly-prime zero-mean perturbation over Q of length ``n_sites``.

    Uses primes 10007, 65537, and 104729 and subtracts the exact rational mean so sum == 0/1.
    """
    primes = (10007, 65537, 104729)
    numerators = (1, -3, 5, -2, 7, -11)
    raw: List[Fraction] = []
    for i in range(n_sites):
        p = primes[i % len(primes)]
        num = numerators[(i * 3 + i // 4) % len(numerators)]
        raw.append(amplitude * Fraction(num, p))
    mean_val = q_sum(raw) / Fraction(n_sites, 1)
    pert = [v - mean_val for v in raw]
    assert q_sum(pert) == ZERO
    return pert


def fmt_exact_zero_or_sci(val: Fraction) -> str:
    """Format Fraction as exact '0/1' if zero, else compact fraction or scientific float."""
    if val == ZERO:
        return "0/1"
    if abs(val.numerator).bit_length() <= 40 and abs(val.denominator).bit_length() <= 40:
        return f"{val.numerator}/{val.denominator}"
    return f"{float(val):.6e}"


def make_json_safe(obj: Any) -> Any:
    """Recursively convert Fraction and non-JSON scalars into clean JSON-serializable structures."""
    if isinstance(obj, Fraction):
        if obj == ZERO:
            return "0/1"
        if abs(obj.numerator).bit_length() <= 50 and abs(obj.denominator).bit_length() <= 50:
            return f"{obj.numerator}/{obj.denominator}"
        return float(obj)
    if isinstance(obj, float):
        if math.isnan(obj):
            return "NaN"
        if math.isinf(obj):
            return "Inf" if obj > 0 else "-Inf"
        return obj
    if isinstance(obj, dict):
        return {str(k): make_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [make_json_safe(v) for v in obj]
    return obj


# --- Component: standalone_part3_runner.py ---
def _standalone_run_b1(
    tag: str,
    m: int,
    nu_4: Fraction,
    dt: Fraction,
    n_steps: int,
    prime_p: Fraction = Fraction(0, 1),
) -> Dict[str, Any]:
    zero = Fraction(0, 1)
    half = Fraction(1, 2)
    grid = get_grid_2d(m)
    st0 = init_battle1(m=m, nu_4=nu_4, prime_perturbation=prime_p)
    f0 = extract_fields_battle1(st0)
    om0 = f0["omega"]
    psi0 = f0["psi"]
    gamma0 = sum(om0, zero)
    ekin0 = half * dot_1d(om0, psi0)
    enst0 = half * dot_1d(om0, om0)
    qheat0 = st0.get("Q_heat", zero)
    etot0 = ekin0 + qheat0
    sh0_raw = grid.energy_spectrum_shells(om0)
    sh0 = {str(k): float(v) for k, v in enumerate(sh0_raw)}
    r1_0 = float(sh0_raw[1]) if len(sh0_raw) > 1 else 0.0
    etot_f = float(ekin0) if ekin0 > zero else 1.0
    ir0 = sum(float(sh0_raw[k]) for k in (1, 2) if k < len(sh0_raw)) / etot_f
    uv0 = sum(float(sh0_raw[k]) for k in range(3, len(sh0_raw))) / etot_f
    max_ir = ir0

    b0 = inspect_bits(st0, f0)
    bits_hist = [b0["max_bits"]]
    st_bits = [b0["state_max_bits"]]
    f_bits = [b0["fields_max_bits"]]

    max_g_err = zero
    max_e_err = zero
    max_z_err = zero
    max_fl_err = zero
    max_pars_err = abs(sum(sh0_raw, zero) - ekin0)
    min_dq: Optional[Fraction] = None
    prev_q = qheat0

    ts: List[Dict[str, Any]] = [{
        "step": 0,
        "t": 0.0,
        "E_kin": float(ekin0),
        "Z": float(enst0),
        "Q_heat": float(qheat0),
        "E_total": float(etot0),
        "ir_frac": ir0,
        "uv_frac": uv0,
        "max_bits": b0["max_bits"],
        "state_bits": b0["state_max_bits"],
        "fields_bits": b0["fields_max_bits"],
    }]
    step_times: List[float] = []
    cur = st0
    for s_idx in range(1, n_steps + 1):
        cur, el = run_guarded_step(step_battle1, cur, dt, timeout_sec=5, label=f"{tag}:s{s_idx}")
        assert el < 5.0, f"{tag} step {s_idx} exceeded 5.0s ({el:.3f}s)"
        step_times.append(el)
        qh = cur.get("Q_heat", zero)
        dq = qh - prev_q
        prev_q = qh
        if min_dq is None or dq < min_dq:
            min_dq = dq

        if s_idx <= 12 or s_idx % 8 == 0 or s_idx == n_steps:
            ff, _ = run_guarded_step(extract_fields_battle1, cur, timeout_sec=5, label=f"{tag}:ext{s_idx}")
            om = ff["omega"]
            psi = ff["psi"]
            ek = half * dot_1d(om, psi)
            zn = half * dot_1d(om, om)
            max_g_err = max(max_g_err, abs(sum(om, zero) - gamma0))
            max_e_err = max(max_e_err, abs(ek - ekin0))
            max_z_err = max(max_z_err, abs(zn - enst0))
            max_fl_err = max(max_fl_err, abs((ek + qh) - etot0))
            if s_idx in (1, 4, 12, n_steps // 2, n_steps):
                sh_r = grid.energy_spectrum_shells(om)
                max_pars_err = max(max_pars_err, abs(sum(sh_r, zero) - ek))
                ir_c = sum(float(sh_r[k]) for k in (1, 2) if k < len(sh_r)) / etot_f
                uv_c = sum(float(sh_r[k]) for k in range(3, len(sh_r))) / etot_f
                max_ir = max(max_ir, ir_c)
            else:
                ir_c = ts[-1]["ir_frac"]
                uv_c = ts[-1]["uv_frac"]
            bi = inspect_bits(cur, ff)
            bits_hist.append(bi["max_bits"])
            st_bits.append(bi["state_max_bits"])
            f_bits.append(bi["fields_max_bits"])
            ts.append({
                "step": s_idx,
                "t": float(s_idx * dt),
                "E_kin": float(ek),
                "Z": float(zn),
                "Q_heat": float(qh),
                "E_total": float(ek + qh),
                "ir_frac": ir_c,
                "uv_frac": uv_c,
                "max_bits": bi["max_bits"],
                "state_bits": bi["state_max_bits"],
                "fields_bits": bi["fields_max_bits"],
            })
        else:
            bi = inspect_bits(cur, None)
            sb = bi["state_max_bits"]
            bits_hist.append(max(sb, f_bits[-1]))
            st_bits.append(sb)
            f_bits.append(f_bits[-1])

    ff_end = extract_fields_battle1(cur)
    plat = check_bit_plateau(bits_hist, bit_cap=10000)
    assert plat["within_cap"], f"{tag} exceeded 10,000-bit ceiling: {plat['peak_bits']}"
    assert max_g_err == zero, f"{tag} circulation drift: {max_g_err}"
    assert max_pars_err == zero, f"{tag} Parseval shell drift: {max_pars_err}"

    rev1 = None
    rev10 = None
    if nu_4 == zero:
        assert max_e_err == zero, f"{tag} inviscid energy drift: {max_e_err}"
        assert max_z_err == zero, f"{tag} inviscid enstrophy drift: {max_z_err}"
        rev1 = check_exact_reversibility(step_battle1, st0, dt, n_steps=1, step_guard_fn=run_guarded_step)
        rev10 = check_exact_reversibility(step_battle1, st0, dt, n_steps=10 if m == 6 else 4, step_guard_fn=run_guarded_step)
        assert rev1["exact_zero_residual"], f"{tag} 1-step reversibility failed"
        assert rev10["exact_zero_residual"], f"{tag} multi-step reversibility failed"
    else:
        assert max_fl_err == zero, f"{tag} viscous First Law drift: {max_fl_err}"
        assert min_dq is not None and min_dq >= zero, f"{tag} viscous Second Law violated: {min_dq}"

    return {
        "tag": tag,
        "m": m,
        "N": grid.N,
        "nu_4": str(nu_4),
        "dt": str(dt),
        "completed_steps": n_steps,
        "n_steps_requested": n_steps,
        "stalled": False,
        "mean_step_sec": sum(step_times) / len(step_times),
        "max_step_sec": max(step_times),
        "spectral_transfer": {
            "initial_shells": sh0,
            "shell_tilde_R1_t0": r1_0,
            "ir_frac_t0": ir0,
            "ir_frac_max": max_ir,
        },
        "bits": {
            "peak_bits": plat["peak_bits"],
            "peak_state_bits": max(st_bits),
            "peak_fields_bits": max(f_bits),
            "late_slope_bits_per_step": plat["late_slope"],
            "within_10000_cap": plat["within_cap"],
            "is_plateau": plat["is_plateau"],
        },
        "reversibility_1step": rev1,
        "reversibility_10step": rev10,
        "time_series": ts,
        "final_omega_float": [float(x) for x in ff_end["omega"]],
    }


def _standalone_run_b2(
    tag: str,
    L: int,
    nu: Fraction,
    dt: Fraction,
    n_steps: int,
    prime_p: Fraction = Fraction(0, 1),
    galilean_boost: Optional[Tuple[Fraction, Fraction, Fraction]] = None,
) -> Dict[str, Any]:
    zero = Fraction(0, 1)
    half = Fraction(1, 2)
    grid = get_grid_3d(L)
    st0 = init_battle2(L=L, nu=nu, prime_perturbation=prime_p, galilean_boost=galilean_boost)
    f0 = extract_fields_battle2(st0)
    u0 = f0["u"]
    om0 = f0["omega"]
    p0 = tuple(sum(u0[c], zero) for c in range(3))
    zero_mom = all(p0[c] == zero for c in range(3))
    ekin0 = half * dot_vec3(u0, u0)
    enst0 = half * dot_vec3(om0, om0)
    qheat0 = st0.get("Q_heat", zero)
    etot0 = ekin0 + qheat0
    b_u = bound_B_u(ekin0, grid.N, zero_momentum=zero_mom)
    b_om = bound_B_omega(ekin0, grid.N)

    sh0_q = grid.roughness_shell_energy_spectrum(u0)
    sh0 = {str(k): float(v) for k, v in sorted(sh0_q.items())}
    max_pars_err = abs(sum(sh0_q.values(), zero) - ekin0)
    assert max_pars_err == zero, f"{tag} initial Parseval shell error: {max_pars_err}"

    b0 = inspect_bits(st0, f0)
    bits_hist = [b0["max_bits"]]
    st_bits = [b0["state_max_bits"]]
    f_bits = [b0["fields_max_bits"]]

    max_div = q_max_abs(grid.div_sym(u0))
    max_p_err = [zero, zero, zero]
    max_e_err = zero
    max_fl_err = zero
    min_dq: Optional[Fraction] = None
    prev_q = qheat0
    z_max = enst0
    z_max_step = 0
    z_max_shells = dict(sh0)

    ts: List[Dict[str, Any]] = [{
        "step": 0,
        "t": 0.0,
        "E_kin": float(ekin0),
        "Z": float(enst0),
        "Z_ratio": 1.0,
        "Q_heat": float(qheat0),
        "E_total": float(etot0),
        "max_bits": b0["max_bits"],
        "state_bits": b0["state_max_bits"],
        "fields_bits": b0["fields_max_bits"],
    }]
    step_times: List[float] = []
    cur = st0
    for s_idx in range(1, n_steps + 1):
        cur, el = run_guarded_step(step_battle2, cur, dt, timeout_sec=5, label=f"{tag}:s{s_idx}")
        assert el < 5.0, f"{tag} step {s_idx} exceeded 5.0s ({el:.3f}s)"
        step_times.append(el)
        qh = cur.get("Q_heat", zero)
        dq = qh - prev_q
        prev_q = qh
        if min_dq is None or dq < min_dq:
            min_dq = dq

        t_phys = s_idx * dt
        in_peak = (t_phys <= Fraction(5, 2))
        if in_peak or s_idx % 8 == 0 or s_idx == n_steps:
            ff, _ = run_guarded_step(extract_fields_battle2, cur, timeout_sec=5, label=f"{tag}:ext{s_idx}")
            u = ff["u"]
            om = ff["omega"]
            max_div = max(max_div, q_max_abs(grid.div_sym(u)))
            for c in range(3):
                max_p_err[c] = max(max_p_err[c], abs(sum(u[c], zero) - p0[c]))
            ek = half * dot_vec3(u, u)
            zn = half * dot_vec3(om, om)
            max_e_err = max(max_e_err, abs(ek - ekin0))
            max_fl_err = max(max_fl_err, abs((ek + qh) - etot0))
            for i in range(grid.N):
                assert u[0][i] ** 2 + u[1][i] ** 2 + u[2][i] ** 2 <= b_u, f"{tag} violated B_u"
                assert om[0][i] ** 2 + om[1][i] ** 2 + om[2][i] ** 2 <= b_om, f"{tag} violated B_omega"
            if zn > z_max:
                z_max = zn
                z_max_step = s_idx
            bi = inspect_bits(cur, ff)
            bits_hist.append(bi["max_bits"])
            st_bits.append(bi["state_max_bits"])
            f_bits.append(bi["fields_max_bits"])
            ts.append({
                "step": s_idx,
                "t": float(t_phys),
                "E_kin": float(ek),
                "Z": float(zn),
                "Z_ratio": float(zn / enst0) if enst0 > zero else 1.0,
                "Q_heat": float(qh),
                "E_total": float(ek + qh),
                "max_bits": bi["max_bits"],
                "state_bits": bi["state_max_bits"],
                "fields_bits": bi["fields_max_bits"],
            })
            if zn == z_max and L == 6:
                sh_pk = grid.roughness_shell_energy_spectrum(u)
                z_max_shells = {str(k): float(v) for k, v in sorted(sh_pk.items())}
        else:
            bi = inspect_bits(cur, None)
            sb = bi["state_max_bits"]
            bits_hist.append(max(sb, f_bits[-1]))
            st_bits.append(sb)
            f_bits.append(f_bits[-1])

    ff_fin = extract_fields_battle2(cur)
    sh_fin_q = grid.roughness_shell_energy_spectrum(ff_fin["u"])
    sh_fin = {str(k): float(v) for k, v in sorted(sh_fin_q.items())}
    ek_fin = half * dot_vec3(ff_fin["u"], ff_fin["u"])
    max_pars_err = max(max_pars_err, abs(sum(sh_fin_q.values(), zero) - ek_fin))
    assert max_div == zero, f"{tag} divergence residual: {max_div}"
    assert all(x == zero for x in max_p_err), f"{tag} momentum drift: {max_p_err}"
    assert max_pars_err == zero, f"{tag} Parseval residual: {max_pars_err}"

    plat = check_bit_plateau(bits_hist, bit_cap=10000)
    assert plat["within_cap"], f"{tag} exceeded 10,000-bit cap: {plat['peak_bits']}"

    rev10 = None
    if nu == zero:
        assert max_e_err == zero, f"{tag} inviscid energy drift: {max_e_err}"
        rev10 = check_exact_reversibility(
            step_battle2, st0, dt, n_steps=8 if L == 6 else 4, step_guard_fn=run_guarded_step
        )
        assert rev10["exact_zero_residual"], f"{tag} reversibility failed"
    else:
        assert max_fl_err == zero, f"{tag} viscous First Law drift: {max_fl_err}"
        assert min_dq is not None and min_dq >= zero, f"{tag} viscous Second Law violated: {min_dq}"

    return {
        "tag": tag,
        "L": L,
        "N": grid.N,
        "nu": str(nu),
        "dt": str(dt),
        "completed_steps": n_steps,
        "n_steps_requested": n_steps,
        "stalled": False,
        "mean_step_sec": sum(step_times) / len(step_times),
        "max_step_sec": max(step_times),
        "stretching": {
            "Z_max_over_Z0": float(z_max / enst0) if enst0 > zero else 1.0,
            "peak_step": z_max_step,
            "peak_t": float(z_max_step * dt),
        },
        "spectral_transfer": {
            "initial_shells": sh0,
            "peak_z_shells": z_max_shells,
            "final_shells": sh_fin,
        },
        "bits": {
            "peak_bits": plat["peak_bits"],
            "peak_state_bits": max(st_bits),
            "peak_fields_bits": max(f_bits),
            "late_slope_bits_per_step": plat["late_slope"],
            "within_10000_cap": plat["within_cap"],
            "is_plateau": plat["is_plateau"],
        },
        "reversibility_10step": rev10,
        "time_series": ts,
    }


def _standalone_run_b3(
    tag: str,
    L: int,
    dt: Fraction,
    n_steps: int,
    prime_p: Fraction = Fraction(0, 1),
) -> Dict[str, Any]:
    zero = Fraction(0, 1)
    half = Fraction(1, 2)
    grid = get_grid_3d(L)
    st0 = init_battle3(L=L, kappa=Fraction(1, 50), amp=Fraction(10, 1), prime_perturbation=prime_p)
    f0 = extract_fields_battle3(st0)
    u0 = f0["u"]
    om0 = f0["omega"]
    rho0 = list(f0["rho"])
    p0 = tuple(sum(u0[c], zero) for c in range(3))
    zero_mom = all(p0[c] == zero for c in range(3))
    ekin0 = half * dot_vec3(u0, u0)
    enst0 = half * dot_vec3(om0, om0)
    mass0 = sum(rho0, zero)
    b_u = bound_B_u(ekin0, grid.N, zero_momentum=zero_mom)
    b_om = bound_B_omega(ekin0, grid.N)

    min_rho_all = min(rho0)
    max_rho_disp = zero
    max_u_sq = max(u0[0][i] ** 2 + u0[1][i] ** 2 + u0[2][i] ** 2 for i in range(grid.N))
    max_om_sq = max(om0[0][i] ** 2 + om0[1][i] ** 2 + om0[2][i] ** 2 for i in range(grid.N))
    b0 = inspect_bits(st0, f0)
    bits_hist = [b0["max_bits"]]
    st_bits = [b0["state_max_bits"]]
    f_bits = [b0["fields_max_bits"]]

    ts: List[Dict[str, Any]] = [{
        "step": 0,
        "t": 0.0,
        "E_kin": float(ekin0),
        "Z": float(enst0),
        "min_rho": float(min_rho_all),
        "max_u_sq": float(max_u_sq),
        "max_omega_sq": float(max_om_sq),
        "max_bits": b0["max_bits"],
        "state_bits": b0["state_max_bits"],
        "fields_bits": b0["fields_max_bits"],
    }]
    step_times: List[float] = []
    cur = st0
    rho_max_t1 = float(max(rho0))
    for s_idx in range(1, n_steps + 1):
        cur, el = run_guarded_step(step_battle3, cur, dt, timeout_sec=5, label=f"{tag}:s{s_idx}")
        assert el < 5.0, f"{tag} step {s_idx} exceeded 5.0s ({el:.3f}s)"
        step_times.append(el)
        rc = cur["rho"]
        assert sum(rc, zero) == mass0, f"{tag} mass drift at step {s_idx}"
        mr_i = min(rc)
        assert mr_i > zero, f"{tag} density non-positive at step {s_idx}: {mr_i}"
        min_rho_all = min(min_rho_all, mr_i)
        max_rho_disp = max(max_rho_disp, max(abs(rc[i] - rho0[i]) for i in range(grid.N)))
        if s_idx * dt == Fraction(1, 1):
            rho_max_t1 = float(max(rc))

        if s_idx <= 10 or s_idx % 8 == 0 or s_idx == n_steps:
            ff, _ = run_guarded_step(extract_fields_battle3, cur, timeout_sec=5, label=f"{tag}:ext{s_idx}")
            u = ff["u"]
            om = ff["omega"]
            assert q_max_abs(grid.div_sym(u)) == zero, f"{tag} div_Sym != 0/1"
            for c in range(3):
                assert sum(u[c], zero) == p0[c], f"{tag} momentum drift"
            assert half * dot_vec3(u, u) == ekin0, f"{tag} energy drift"
            musq_i = max(u[0][i] ** 2 + u[1][i] ** 2 + u[2][i] ** 2 for i in range(grid.N))
            momsq_i = max(om[0][i] ** 2 + om[1][i] ** 2 + om[2][i] ** 2 for i in range(grid.N))
            assert musq_i <= b_u, f"{tag} violated B_u: {float(musq_i)} > {float(b_u)}"
            assert momsq_i <= b_om, f"{tag} violated B_omega: {float(momsq_i)} > {float(b_om)}"
            max_u_sq = max(max_u_sq, musq_i)
            max_om_sq = max(max_om_sq, momsq_i)
            bi = inspect_bits(cur, ff)
            bits_hist.append(bi["max_bits"])
            st_bits.append(bi["state_max_bits"])
            f_bits.append(bi["fields_max_bits"])
            ts.append({
                "step": s_idx,
                "t": float(s_idx * dt),
                "E_kin": float(ekin0),
                "Z": float(half * dot_vec3(om, om)),
                "min_rho": float(mr_i),
                "max_u_sq": float(musq_i),
                "max_omega_sq": float(momsq_i),
                "max_bits": bi["max_bits"],
                "state_bits": bi["state_max_bits"],
                "fields_bits": bi["fields_max_bits"],
            })
        else:
            bi = inspect_bits(cur, None)
            sb = bi["state_max_bits"]
            bits_hist.append(max(sb, f_bits[-1]))
            st_bits.append(sb)
            f_bits.append(f_bits[-1])

    # Velocity-sector 4-step reversibility
    fwd = st0
    for k in range(4):
        fwd, _ = run_guarded_step(step_battle3, fwd, dt, timeout_sec=5, label=f"{tag}:rf{k}")
    bwd = fwd
    for k in range(4):
        bwd, _ = run_guarded_step(step_battle3, bwd, -dt, timeout_sec=5, label=f"{tag}:rb{k}")
    for slot in ("alpha", "eta_a", "eta_b"):
        if slot in st0:
            assert bwd[slot] == st0[slot], f"{tag} velocity slot {slot} failed 4-step reversibility"

    plat = check_bit_plateau(bits_hist, bit_cap=10000)
    assert plat["within_cap"], f"{tag} exceeded 10,000-bit cap: {plat['peak_bits']}"
    assert float(max_rho_disp) > 1.0, f"{tag} density frozen (max_disp={float(max_rho_disp)})"

    return {
        "tag": tag,
        "L": L,
        "N": grid.N,
        "dt": str(dt),
        "completed_steps": n_steps,
        "n_steps_requested": n_steps,
        "stalled": False,
        "mean_step_sec": sum(step_times) / len(step_times),
        "max_step_sec": max(step_times),
        "positivity_and_bounds": {
            "min_rho_float": float(min_rho_all),
            "rho_max_t1_float": rho_max_t1,
            "max_rho_displacement_float": float(max_rho_disp),
            "max_u_sq_float": float(max_u_sq),
            "B_u_float": float(b_u),
            "u_bound_utilization": float(max_u_sq / b_u),
            "max_omega_sq_float": float(max_om_sq),
            "B_omega_float": float(b_om),
            "omega_bound_utilization": float(max_om_sq / b_om),
        },
        "bits": {
            "peak_bits": plat["peak_bits"],
            "peak_state_bits": max(st_bits),
            "peak_fields_bits": max(f_bits),
            "late_slope_bits_per_step": plat["late_slope"],
            "within_10000_cap": plat["within_cap"],
            "is_plateau": plat["is_plateau"],
        },
        "time_series": ts,
    }


def run_standalone_verification_and_plots(output_plots_dir: str) -> Dict[str, Any]:
    """Execute full Part III self-verifying assertion suite & Part IV pure-stdlib plot rendering."""
    t0 = time.perf_counter()
    zero = Fraction(0, 1)
    print("=" * 78)
    print("CONSTRUCTIVE ULTRAFINITIST FLUID TURBULENCE & NO-BLOWUP ENGINE OVER Q")
    print("Standalone Zero-Dependency Verification Harness (Parts I - IV)")
    print("=" * 78)

    # 1. Gate 4 Zero-Handicap Continuum Float64 Baselines & Convergence Orders
    print("[1/5] Verifying Gate 4 zero-handicap continuum float64 baselines & slopes...")
    f64_data: Dict[str, Any] = {
        "gate4_slope_verification": {},
        "coarse_baseline_comparisons": {},
        "battle2_rk4_ref_stretching": {},
        "battle3_centered_gibbs_ringing": {},
        "q_engine_convergence": {},
        "small_dt_tangents": {},
    }
    _, y0_b1a6, op_b1a6, eq_b1a6 = build_verified_battle1_baseline(m=6, nu_4=zero)
    _, y0_b1a7, op_b1a7, eq_b1a7 = build_verified_battle1_baseline(m=7, nu_4=zero)
    _, y0_b1b6, op_b1b6, eq_b1b6 = build_verified_battle1_baseline(m=6, nu_4=Fraction(1, 200))
    _, y0_b2a6, op_b2a6, eq_b2a6 = build_verified_battle2_baseline(L=6, nu=zero)
    _, y0_b2a8, op_b2a8, eq_b2a8 = build_verified_battle2_baseline(L=8, nu=zero)
    _, y0_b2b6, op_b2b6, eq_b2b6 = build_verified_battle2_baseline(L=6, nu=Fraction(1, 50))
    _, y0_b3c, op_b3c, eq_b3c = build_verified_battle3_baseline(L=6, centered_density=True)
    _, y0_b3v, op_b3v, eq_b3v = build_verified_battle3_baseline(L=6, centered_density=False)

    for cname, y0_v, op_v, eq_v, t_ev, dts_v in [
        ("B1A_m6", y0_b1a6, op_b1a6, eq_b1a6, 0.5, (1 / 8, 1 / 16, 1 / 32, 1 / 64)),
        ("B1A_m7", y0_b1a7, op_b1a7, eq_b1a7, 0.5, (1 / 8, 1 / 16, 1 / 32, 1 / 64)),
        ("B1B_m6", y0_b1b6, op_b1b6, eq_b1b6, 0.5, (1 / 8, 1 / 16, 1 / 32, 1 / 64)),
        ("B2A_L6", y0_b2a6, op_b2a6, eq_b2a6, 0.5, (1 / 8, 1 / 16, 1 / 32, 1 / 64)),
        ("B2A_L8", y0_b2a8, op_b2a8, eq_b2a8, 0.25, (1 / 8, 1 / 16, 1 / 32)),
        ("B2B_L6", y0_b2b6, op_b2b6, eq_b2b6, 0.5, (1 / 8, 1 / 16, 1 / 32, 1 / 64)),
        ("B3_L6_centered", y0_b3c, op_b3c, eq_b3c, 0.125, (1 / 32, 1 / 64, 1 / 128)),
        ("B3_L6_varga_gen", y0_b3v, op_b3v, eq_b3v, 0.125, (1 / 32, 1 / 64, 1 / 128)),
    ]:
        conv_v = verify_convergence_slopes(y0=y0_v, T_eval=t_ev, rhs=op_v.rhs, dts=dts_v, dt_ref=1.0 / 512.0)
        assert conv_v["all_passed"], f"Gate 4 slope failure on {cname}: {conv_v}"
        f64_data["gate4_slope_verification"][cname] = {"rhs_equivalence": eq_v, "convergence": conv_v}

    for coarse_tag, y0_c, op_c, n_c in [
        ("B1A_m6", y0_b1a6, op_b1a6, 120),
        ("B2A_L6", y0_b2a6, op_b2a6, 80),
        ("B3_L6_centered", y0_b3c, op_b3c, 80),
    ]:
        f64_data["coarse_baseline_comparisons"][coarse_tag] = {}
        for dt_f, dt_lbl in [(0.5, "1/2"), (0.25, "1/4")]:
            f64_data["coarse_baseline_comparisons"][coarse_tag][dt_lbl] = {}
            for meth in ("forward_euler", "midpoint_rk2", "classical_rk4"):
                tr = run_float64_trajectory(
                    y0=y0_c,
                    dt=dt_f,
                    n_steps=n_c,
                    method=meth,
                    rhs=op_c.rhs,
                    observer=lambda s, t, y, _op=op_c: _op.invariants(y) if (s % 4 == 0 or s <= 16) else {},
                )
                f64_data["coarse_baseline_comparisons"][coarse_tag][dt_lbl][meth] = {
                    "blew_up": tr["blew_up"],
                    "blowup_step": tr["blowup_step"],
                    "completed_steps": tr["completed_steps"],
                    "history": [h for h in tr["history"] if "E_kin" in h][:40],
                }
    # Assert coarse explicit RK blowup on extreme Battle 3 collision
    for m_k in ("forward_euler", "midpoint_rk2", "classical_rk4"):
        assert f64_data["coarse_baseline_comparisons"]["B3_L6_centered"]["1/4"][m_k]["blew_up"], (
            f"Expected coarse {m_k} blow-up on Battle 3 centered at dt=1/4"
        )

    # Centered continuum Gibbs ringing check (min_rho < 0)
    ref_b3c = run_rk4_reference(
        y0=y0_b3c, T_final=1.0, rhs=op_b3c.rhs, dt_ref=1.0 / 512.0, record_times=[1.0],
        observer=lambda s, t, y: op_b3c.invariants(y), observer_stride=4,
    )
    min_rho_gibbs = min(h["min_rho"] for h in ref_b3c["history"])
    assert min_rho_gibbs < 0.0, f"Expected negative Gibbs ringing on centered Battle 3, got {min_rho_gibbs}"
    f64_data["battle3_centered_gibbs_ringing"] = {
        "min_rho_continuum_centered": min_rho_gibbs,
        "history": ref_b3c["history"],
    }
    print(f"  [PASS] Gate 4 slopes O(dt^1, dt^2, dt^4) verified; centered Gibbs min_rho = {min_rho_gibbs:.3f} < 0")

    # 2. Battle 2 RK4(dt=1/512) Stretching References & Q-Engine O(dt^2) Convergence
    print("[2/5] Verifying Q-engine O(dt^2) convergence & 3D vortex stretching against RK4(dt=1/512)...")
    for L_val, y0_2, op_2, t_end in [(6, y0_b2a6, op_b2a6, 2.5), (8, y0_b2a8, op_b2a8, 2.0)]:
        z0_2 = op_2.invariants(y0_2)["Z"]
        rec_ts = [1.0, 2.0] if L_val == 6 else [1.0]
        ref_2 = run_rk4_reference(
            y0=y0_2, T_final=t_end, rhs=op_2.rhs, dt_ref=1.0 / 512.0,
            record_times=rec_ts, observer=lambda s, t, y, _op=op_2: _op.invariants(y), observer_stride=8,
        )
        hist_2 = [{"t": h["t"], "Z_ratio": h["Z"] / z0_2, "E_kin": h["E_kin"]} for h in ref_2["history"]]
        f64_data["battle2_rk4_ref_stretching"][f"L={L_val}"] = {
            "Z_max_over_Z0": max(h["Z_ratio"] for h in hist_2),
            "history": hist_2,
        }
        if L_val == 6:
            dts_q = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]
            u_e1, u_e2, zr1, zr2 = [], [], [], []
            z_t1_r = op_2.invariants(ref_2["time_snapshots"][1.0])["Z"] / z0_2
            z_t2_r = op_2.invariants(ref_2["time_snapshots"][2.0])["Z"] / z0_2
            for dq in dts_q:
                st = init_battle2(L=6, nu=zero)
                n1 = int(Fraction(1, 1) / dq)
                for _ in range(n1):
                    st = step_battle2(st, dq)
                u1 = [float(x) for c in range(3) for x in extract_velocity_battle2(st)[c]]
                u_e1.append(vec_error_norms(u1, ref_2["time_snapshots"][1.0])["l2_rms"])
                zr1.append(op_2.invariants(u1)["Z"] / z0_2)
                for _ in range(n1):
                    st = step_battle2(st, dq)
                u2 = [float(x) for c in range(3) for x in extract_velocity_battle2(st)[c]]
                u_e2.append(vec_error_norms(u2, ref_2["time_snapshots"][2.0])["l2_rms"])
                zr2.append(op_2.invariants(u2)["Z"] / z0_2)
            dts_f = [float(d) for d in dts_q]
            ord1 = compute_pairwise_orders(dts_f, u_e1)
            ord2 = compute_pairwise_orders(dts_f, u_e2)
            assert ord1[-1] is not None and ord1[-1] > 1.75, f"Battle 2A T=1 order too low: {ord1}"
            assert ord2[-1] is not None and ord2[-1] > 1.75, f"Battle 2A T=2 order too low: {ord2}"
            f64_data["q_engine_convergence"]["B2A_L6"] = {
                "dts_float": dts_f,
                "T1": {"u_l2_rms_errors": u_e1, "u_orders": ord1, "Z_ratio_values": zr1, "Z_ratio_ref": z_t1_r},
                "T2": {"u_l2_rms_errors": u_e2, "u_orders": ord2, "Z_ratio_values": zr2, "Z_ratio_ref": z_t2_r},
            }

    ref_1 = run_rk4_reference(y0=y0_b1a6, T_final=1.0, rhs=op_b1a6.rhs, dt_ref=1.0 / 512.0, record_times=[1.0])
    dts_q1 = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]
    errs_1: List[float] = []
    for dq in dts_q1:
        st1 = init_battle1(m=6, nu_4=zero)
        for _ in range(int(Fraction(1, 1) / dq)):
            st1 = step_battle1(st1, dq)
        errs_1.append(vec_error_norms([float(x) for x in extract_omega_battle1(st1)], ref_1["time_snapshots"][1.0])["l2_rms"])
    dts_f1 = [float(d) for d in dts_q1]
    ords_1 = compute_pairwise_orders(dts_f1, errs_1)
    assert ords_1[-1] is not None and ords_1[-1] > 1.90, f"Battle 1A order too low: {ords_1}"
    f64_data["q_engine_convergence"]["B1A_m6"] = {
        "dts_float": dts_f1,
        "omega_l2_rms_errors": errs_1,
        "orders": ords_1,
    }
    print(f"  [PASS] Q-engine O(dt^2) orders verified: B1A_m6={ords_1[-1]:.3f}, B2A_L6(T=1)={ord1[-1]:.3f}, (T=2)={ord2[-1]:.3f}")

    # 3. Small-dt Tangents down to 2^-30 (at t=0 AND off-manifold t=2.0) + L=8 Shift Equivariance
    print("[3/5] Verifying small-dt tangents down to 2^-30 (t=0 & t>0) and L=8 shift equivariance...")
    s1_0 = init_battle1(m=6, nu_4=zero)
    s1b_0 = init_battle1(m=6, nu_4=Fraction(1, 200))
    s2_0 = init_battle2(L=6, nu=zero)
    s3_0 = init_battle3(L=6)

    def _ext_b3(st: Any) -> List[Fraction]:
        f = extract_fields_battle3(st)
        return [x for c in range(3) for x in f["u"][c]] + list(f["rho"])

    def _rhs_b3_v(st: Any) -> List[Fraction]:
        f = extract_fields_battle3(st)
        du, drho = continuum_rhs_battle3(f["u"], f["rho"], L=6, kappa=st.get("kappa", Fraction(1, 50)), centered_density=False)
        return [x for c in range(3) for x in du[c]] + list(drho)

    for t_key, s_fn, ex_fn, rh_fn, st_init in [
        ("B1A_m6_step_battle1", step_battle1, extract_omega_battle1, lambda st: continuum_rhs_battle1(extract_omega_battle1(st), m=6, nu_4=zero), s1_0),
        ("B1B_m6_step_battle1", step_battle1, extract_omega_battle1, lambda st: continuum_rhs_battle1(extract_omega_battle1(st), m=6, nu_4=Fraction(1, 200)), s1b_0),
        ("B2A_L6_step_battle2", step_battle2, lambda st: [x for c in range(3) for x in extract_velocity_battle2(st)[c]], lambda st: [x for c in range(3) for x in continuum_rhs_battle2(extract_velocity_battle2(st), L=6, nu=zero)[c]], s2_0),
        ("B3_L6_step_battle3", step_battle3, _ext_b3, _rhs_b3_v, s3_0),
    ]:
        t_res = check_small_dt_tangent(s_fn, ex_fn, rh_fn, st_init, exponents=(4, 10, 20, 30))
        assert t_res["non_frozen_at_2_minus_30"], f"Small-dt tangent frozen on {t_key}: {t_res}"
        f64_data["small_dt_tangents"][t_key] = t_res

    # L=8 shift-equivariance check (C2-2)
    g8 = get_grid_3d(8)
    u8 = extract_velocity_battle2(init_battle2(L=8, nu=zero, prime_perturbation=Fraction(3, 65537)))
    sh8_0 = g8.roughness_shell_energy_spectrum(u8)
    u8_shift = [[zero] * 512 for _ in range(3)]
    for x in range(8):
        for y in range(8):
            for z in range(8):
                src = (x * 8 + y) * 8 + z
                dst = (((x + 1) % 8) * 8 + ((y + 2) % 8)) * 8 + ((z + 3) % 8)
                for c in range(3):
                    u8_shift[c][dst] = u8[c][src]
    sh8_s = g8.roughness_shell_energy_spectrum((u8_shift[0], u8_shift[1], u8_shift[2]))
    for k in sh8_0:
        assert sh8_0[k] == sh8_s[k], f"L=8 shell {k} broke shift equivariance: {sh8_0[k]} != {sh8_s[k]}"
    print("  [PASS] All small-dt tangents unfrozen down to 2^-30 and L=8 shell spectrum 0/1 shift-equivariant")

    # 4. Full Multi-Step Q-Engine Assertions Across Battles 1, 2, and 3
    print("[4/5] Running multi-step exact Q-engine trajectories across Battles 1, 2, and 3...")
    q_runs: Dict[str, Any] = {}

    # Battle 1
    q_runs["B1A_m6_dt1_2"] = _standalone_run_b1("B1A_m6_dt1_2", 6, zero, Fraction(1, 2), 300)
    q_runs["B1A_m6_dt1_4"] = _standalone_run_b1("B1A_m6_dt1_4", 6, zero, Fraction(1, 4), 300)
    q_runs["B1A_m6_dt1_4_prime"] = _standalone_run_b1("B1A_m6_dt1_4_prime", 6, zero, Fraction(1, 4), 120, Fraction(3, 65537))
    q_runs["B1A_m7_dt1_4"] = _standalone_run_b1("B1A_m7_dt1_4", 7, zero, Fraction(1, 4), 80)
    q_runs["B1A_m7_dt1_2_prime"] = _standalone_run_b1("B1A_m7_dt1_2_prime", 7, zero, Fraction(1, 2), 40, Fraction(3, 65537))
    assert q_runs["B1A_m7_dt1_4"]["spectral_transfer"]["shell_tilde_R1_t0"] > 1.0, "Twisted m=7 shell R=1 not excited"
    assert q_runs["B1A_m6_dt1_4"]["spectral_transfer"]["ir_frac_max"] > q_runs["B1A_m6_dt1_4"]["spectral_transfer"]["ir_frac_t0"], "Expected 2D inverse IR transfer"
    assert abs(q_runs["B1A_m6_dt1_4"]["bits"]["late_slope_bits_per_step"]) < 0.05, (
        f"Expected O(1) bit plateau on B1A_m6_dt1_4, got slope {q_runs['B1A_m6_dt1_4']['bits']['late_slope_bits_per_step']}"
    )
    q_runs["B1B_m6_dt1_4"] = _standalone_run_b1("B1B_m6_dt1_4", 6, Fraction(1, 200), Fraction(1, 4), 120)

    # Battle 2
    q_runs["B2A_L6_dt1_2"] = _standalone_run_b2("B2A_L6_dt1_2", 6, zero, Fraction(1, 2), 60)
    q_runs["B2A_L6_dt1_4"] = _standalone_run_b2("B2A_L6_dt1_4", 6, zero, Fraction(1, 4), 80)
    q_runs["B2A_L6_dt1_8"] = _standalone_run_b2("B2A_L6_dt1_8", 6, zero, Fraction(1, 8), 40)
    q_runs["B2A_L6_dt1_4_prime"] = _standalone_run_b2("B2A_L6_dt1_4_prime", 6, zero, Fraction(1, 4), 40, Fraction(3, 65537))
    q_runs["B2A_L6_dt1_4_galilean"] = _standalone_run_b2(
        "B2A_L6_dt1_4_galilean", 6, zero, Fraction(1, 4), 30, zero, (Fraction(1, 3), Fraction(-1, 5), Fraction(2, 7))
    )
    q_runs["B2A_L8_dt1_8"] = _standalone_run_b2("B2A_L8_dt1_8", 8, zero, Fraction(1, 8), 24)
    q_runs["B2A_L8_dt1_4_prime"] = _standalone_run_b2("B2A_L8_dt1_4_prime", 8, zero, Fraction(1, 4), 16, Fraction(3, 65537))
    q_runs["B2B_L6_dt1_4"] = _standalone_run_b2("B2B_L6_dt1_4", 6, Fraction(1, 50), Fraction(1, 4), 60)
    assert q_runs["B2A_L6_dt1_8"]["stretching"]["Z_max_over_Z0"] >= 1.35, "Battle 2 L=6 failed Z_max/Z0 >= 1.35"
    assert q_runs["B2A_L8_dt1_8"]["stretching"]["Z_max_over_Z0"] >= 1.35, "Battle 2 L=8 failed Z_max/Z0 >= 1.35"

    # Battle 3
    q_runs["B3_L6_dt1_2"] = _standalone_run_b3("B3_L6_dt1_2", 6, Fraction(1, 2), 60)
    q_runs["B3_L6_dt1_4"] = _standalone_run_b3("B3_L6_dt1_4", 6, Fraction(1, 4), 80)
    q_runs["B3_L6_dt1_4_prime"] = _standalone_run_b3("B3_L6_dt1_4_prime", 6, Fraction(1, 4), 40, Fraction(1, 10007))
    q_runs["B3_L8_dt1_4"] = _standalone_run_b3("B3_L8_dt1_4", 8, Fraction(1, 4), 16)
    assert q_runs["B3_L6_dt1_4"]["positivity_and_bounds"]["rho_max_t1_float"] > 5.0, (
        f"Battle 3 density plume over-diffused: rho_max(t=1)={q_runs['B3_L6_dt1_4']['positivity_and_bounds']['rho_max_t1_float']}"
    )

    peak_bits_all = max(r["bits"]["peak_bits"] for r in q_runs.values())
    max_sec_all = max(r["max_step_sec"] for r in q_runs.values())
    print(f"  [PASS] All {len(q_runs)} trajectories passed 0/1 invariants & reversibility (peak={peak_bits_all}b, max_step={max_sec_all:.3f}s)")

    # 5. Part IV Pure-Stdlib Diagnostic Plot Generation (.svg and .png)
    print(f"[5/5] Rendering 4 pure-stdlib diagnostic .svg and .png plot pairs in {output_plots_dir}...")
    data_bundle = {
        "cycle": 3,
        "wall_clock_sec": time.perf_counter() - t0,
        "q_runs": q_runs,
        "float64_and_convergence": f64_data,
    }
    rendered = render_all_plots_from_data(data_bundle, output_plots_dir)
    assert len(rendered) == 4, f"Expected 4 plot pairs, got {len(rendered)}"
    for svg_p, png_p in rendered:
        assert os.path.isfile(svg_p) and os.path.getsize(svg_p) > 1024, f"Invalid SVG: {svg_p}"
        assert os.path.isfile(png_p) and os.path.getsize(png_p) > 1024, f"Invalid PNG: {png_p}"
        with open(svg_p, "rb") as fp:
            assert b"<svg" in fp.read(256), f"Missing <svg header in {svg_p}"
        with open(png_p, "rb") as fp:
            assert fp.read(8) == b"\x89PNG\r\n\x1a\n", f"Invalid PNG signature in {png_p}"
        print(f"  [OK] {os.path.basename(svg_p)} & {os.path.basename(png_p)}")

    print("=" * 78)
    print(f"ALL ASSERTIONS PASSED IN {time.perf_counter() - t0:.2f}s (EXIT 0)")
    print("=" * 78)
    return data_bundle


# =============================================================================
# PART IV: PURE-STDLIB DIAGNOSTIC SVG & ZLIB PNG RENDERER
# =============================================================================

# --- Component: plotting.py ---
Color = Tuple[int, int, int]

# Colorblind-safe Okabe-Ito inspired scientific palette
PALETTE: Dict[str, Color] = {
    "navy": (20, 54, 93),
    "blue": (0, 114, 178),
    "orange": (230, 159, 0),
    "vermilion": (213, 94, 0),
    "green": (0, 158, 115),
    "purple": (204, 121, 167),
    "sky": (86, 180, 233),
    "crimson": (180, 28, 46),
    "slate": (90, 102, 118),
    "dark": (28, 33, 40),
    "grid": (222, 227, 234),
    "panel_bg": (255, 255, 255),
    "fig_bg": (246, 248, 251),
    "gold": (196, 142, 12),
}

# 5x7 bitmap font: each glyph is 7 rows of 5-bit integers (bit 4 = leftmost pixel, bit 0 = rightmost)
FONT_5X7: Dict[str, Tuple[int, int, int, int, int, int, int]] = {
    " ": (0, 0, 0, 0, 0, 0, 0),
    "!": (4, 4, 4, 4, 4, 0, 4),
    '"': (10, 10, 0, 0, 0, 0, 0),
    "#": (10, 31, 10, 10, 31, 10, 0),
    "$": (4, 15, 20, 14, 5, 30, 4),
    "%": (25, 25, 2, 4, 8, 19, 19),
    "&": (12, 18, 20, 8, 21, 18, 13),
    "'": (4, 4, 0, 0, 0, 0, 0),
    "(": (2, 4, 8, 8, 8, 4, 2),
    ")": (8, 4, 2, 2, 2, 4, 8),
    "*": (0, 4, 21, 14, 21, 4, 0),
    "+": (0, 4, 4, 31, 4, 4, 0),
    ",": (0, 0, 0, 0, 4, 4, 8),
    "-": (0, 0, 0, 31, 0, 0, 0),
    ".": (0, 0, 0, 0, 0, 6, 6),
    "/": (1, 1, 2, 4, 8, 16, 16),
    "0": (14, 17, 19, 21, 25, 17, 14),
    "1": (4, 12, 4, 4, 4, 4, 14),
    "2": (14, 17, 1, 6, 8, 16, 31),
    "3": (31, 2, 4, 6, 1, 17, 14),
    "4": (2, 6, 10, 18, 31, 2, 2),
    "5": (31, 16, 30, 1, 1, 17, 14),
    "6": (6, 8, 16, 30, 17, 17, 14),
    "7": (31, 1, 2, 4, 8, 8, 8),
    "8": (14, 17, 17, 14, 17, 17, 14),
    "9": (14, 17, 17, 15, 1, 2, 12),
    ":": (0, 6, 6, 0, 6, 6, 0),
    ";": (0, 6, 6, 0, 6, 4, 8),
    "<": (2, 4, 8, 16, 8, 4, 2),
    "=": (0, 0, 31, 0, 31, 0, 0),
    ">": (8, 4, 2, 1, 2, 4, 8),
    "?": (14, 17, 1, 2, 4, 0, 4),
    "@": (14, 17, 23, 21, 23, 16, 14),
    "A": (14, 17, 17, 31, 17, 17, 17),
    "B": (30, 17, 17, 30, 17, 17, 30),
    "C": (14, 17, 16, 16, 16, 17, 14),
    "D": (30, 17, 17, 17, 17, 17, 30),
    "E": (31, 16, 16, 30, 16, 16, 31),
    "F": (31, 16, 16, 30, 16, 16, 16),
    "G": (14, 17, 16, 23, 17, 17, 15),
    "H": (17, 17, 17, 31, 17, 17, 17),
    "I": (14, 4, 4, 4, 4, 4, 14),
    "J": (7, 2, 2, 2, 2, 18, 12),
    "K": (17, 18, 20, 24, 20, 18, 17),
    "L": (16, 16, 16, 16, 16, 16, 31),
    "M": (17, 27, 21, 21, 17, 17, 17),
    "N": (17, 25, 21, 19, 17, 17, 17),
    "O": (14, 17, 17, 17, 17, 17, 14),
    "P": (30, 17, 17, 30, 16, 16, 16),
    "Q": (14, 17, 17, 17, 21, 18, 13),
    "R": (30, 17, 17, 30, 20, 18, 17),
    "S": (15, 16, 16, 14, 1, 1, 30),
    "T": (31, 4, 4, 4, 4, 4, 4),
    "U": (17, 17, 17, 17, 17, 17, 14),
    "V": (17, 17, 17, 17, 17, 10, 4),
    "W": (17, 17, 17, 21, 21, 27, 17),
    "X": (17, 17, 10, 4, 10, 17, 17),
    "Y": (17, 17, 10, 4, 4, 4, 4),
    "Z": (31, 1, 2, 4, 8, 16, 31),
    "[": (14, 8, 8, 8, 8, 8, 14),
    "\\": (16, 16, 8, 4, 2, 1, 1),
    "]": (14, 2, 2, 2, 2, 2, 14),
    "^": (4, 10, 17, 0, 0, 0, 0),
    "_": (0, 0, 0, 0, 0, 0, 31),
    "`": (8, 4, 0, 0, 0, 0, 0),
    "a": (0, 0, 14, 1, 15, 17, 15),
    "b": (16, 16, 30, 17, 17, 17, 30),
    "c": (0, 0, 14, 16, 16, 17, 14),
    "d": (1, 1, 15, 17, 17, 17, 15),
    "e": (0, 0, 14, 17, 31, 16, 14),
    "f": (6, 9, 8, 28, 8, 8, 8),
    "g": (0, 15, 17, 17, 15, 1, 14),
    "h": (16, 16, 22, 25, 17, 17, 17),
    "i": (4, 0, 12, 4, 4, 4, 14),
    "j": (2, 0, 6, 2, 2, 18, 12),
    "k": (16, 16, 18, 20, 24, 20, 18),
    "l": (12, 4, 4, 4, 4, 4, 14),
    "m": (0, 0, 26, 21, 21, 17, 17),
    "n": (0, 0, 22, 25, 17, 17, 17),
    "o": (0, 0, 14, 17, 17, 17, 14),
    "p": (0, 0, 30, 17, 30, 16, 16),
    "q": (0, 0, 15, 17, 15, 1, 1),
    "r": (0, 0, 22, 25, 16, 16, 16),
    "s": (0, 0, 15, 16, 14, 1, 30),
    "t": (8, 8, 28, 8, 8, 9, 6),
    "u": (0, 0, 17, 17, 17, 19, 13),
    "v": (0, 0, 17, 17, 17, 10, 4),
    "w": (0, 0, 17, 17, 21, 21, 10),
    "x": (0, 0, 17, 10, 4, 10, 17),
    "y": (0, 0, 17, 17, 15, 1, 14),
    "z": (0, 0, 31, 2, 4, 8, 31),
    "{": (2, 4, 4, 8, 4, 4, 2),
    "|": (4, 4, 4, 4, 4, 4, 4),
    "}": (8, 4, 4, 2, 4, 4, 8),
    "~": (0, 0, 8, 21, 2, 0, 0),
}


def _rgb_hex(c: Color) -> str:
    return f"#{c[0]:02x}{c[1]:02x}{c[2]:02x}"


def encode_png_rgb(width: int, height: int, rgb_bytes: bytes) -> bytes:
    """Encode a raw ``width * height * 3`` RGB buffer into a valid PNG byte stream."""
    if len(rgb_bytes) != width * height * 3:
        raise ValueError(f"RGB buffer size {len(rgb_bytes)} != {width * height * 3}")

    def _chunk(chunk_type: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(chunk_type)
        crc = zlib.crc32(data, crc) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + chunk_type + data + struct.pack(">I", crc)

    # Prepend filter byte 0 (None) to each scanline
    row_stride = width * 3
    raw_scanlines = bytearray((row_stride + 1) * height)
    dst = 0
    src = 0
    for _ in range(height):
        raw_scanlines[dst] = 0
        raw_scanlines[dst + 1 : dst + 1 + row_stride] = rgb_bytes[src : src + row_stride]
        dst += row_stride + 1
        src += row_stride

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    compressed = zlib.compress(bytes(raw_scanlines), level=6)
    return (
        b"\x89PNG\r\n\x1a\n"
        + _chunk(b"IHDR", ihdr)
        + _chunk(b"IDAT", compressed)
        + _chunk(b"IEND", b"")
    )


def viridis_like_colormap(t: float) -> Color:
    """Map normalized scalar ``t in [0, 1]`` to a clean perceptual scientific RGB color."""
    if math.isnan(t):
        return (180, 180, 180)
    u = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
    stops: Tuple[Tuple[float, Color], ...] = (
        (0.00, (44, 24, 84)),
        (0.25, (49, 104, 142)),
        (0.50, (33, 145, 140)),
        (0.75, (94, 201, 98)),
        (1.00, (253, 231, 37)),
    )
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if u <= t1:
            w = (u - t0) / (t1 - t0) if t1 > t0 else 0.0
            return (
                int(round(c0[0] + w * (c1[0] - c0[0]))),
                int(round(c0[1] + w * (c1[1] - c0[1]))),
                int(round(c0[2] + w * (c1[2] - c0[2]))),
            )
    return stops[-1][1]


def diverging_bwr_colormap(t: float) -> Color:
    """Map ``t in [-1, 1]`` to a crisp blue-white-crimson diverging colormap."""
    if math.isnan(t):
        return (180, 180, 180)
    u = -1.0 if t < -1.0 else (1.0 if t > 1.0 else t)
    if u < 0.0:
        w = u + 1.0  # 0 at -1, 1 at 0
        c0, c1 = (18, 78, 158), (248, 249, 251)
    else:
        w = u  # 0 at 0, 1 at +1
        c0, c1 = (248, 249, 251), (186, 30, 45)
    return (
        int(round(c0[0] + w * (c1[0] - c0[0]))),
        int(round(c0[1] + w * (c1[1] - c0[1]))),
        int(round(c0[2] + w * (c1[2] - c0[2]))),
    )


class DualFigure:
    """Simultaneous vector SVG + RGB raster PNG figure canvas."""

    def __init__(self, width: int = 1360, height: int = 900, bg: Color = PALETTE["fig_bg"]) -> None:
        self.width = width
        self.height = height
        self.bg = bg
        self.pixels = bytearray(bytes(bg) * (width * height))
        self.svg_elements: List[str] = [
            f'<rect x="0" y="0" width="{width}" height="{height}" fill="{_rgb_hex(bg)}"/>'
        ]

    def _set_px(self, x: int, y: int, c: Color) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            idx = (y * self.width + x) * 3
            self.pixels[idx] = c[0]
            self.pixels[idx + 1] = c[1]
            self.pixels[idx + 2] = c[2]

    def _blend_px(self, x: int, y: int, c: Color, alpha: float) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            idx = (y * self.width + x) * 3
            inv = 1.0 - alpha
            self.pixels[idx] = int(self.pixels[idx] * inv + c[0] * alpha)
            self.pixels[idx + 1] = int(self.pixels[idx + 1] * inv + c[1] * alpha)
            self.pixels[idx + 2] = int(self.pixels[idx + 2] * inv + c[2] * alpha)

    def fill_rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: Optional[Color] = None,
        stroke: Optional[Color] = None,
        stroke_width: int = 1,
        alpha: float = 1.0,
    ) -> None:
        ix0 = max(0, int(round(x)))
        iy0 = max(0, int(round(y)))
        ix1 = min(self.width, int(round(x + w)))
        iy1 = min(self.height, int(round(y + h)))

        if fill is not None and ix1 > ix0 and iy1 > iy0:
            if alpha >= 0.999:
                row_bytes = bytes(fill) * (ix1 - ix0)
                w_stride = self.width * 3
                for py in range(iy0, iy1):
                    start = py * w_stride + ix0 * 3
                    self.pixels[start : start + len(row_bytes)] = row_bytes
            else:
                for py in range(iy0, iy1):
                    for px in range(ix0, ix1):
                        self._blend_px(px, py, fill, alpha)

        if stroke is not None and stroke_width > 0:
            self._raster_line(x, y, x + w, y, stroke, stroke_width)
            self._raster_line(x + w, y, x + w, y + h, stroke, stroke_width)
            self._raster_line(x + w, y + h, x, y + h, stroke, stroke_width)
            self._raster_line(x, y + h, x, y, stroke, stroke_width)

        fill_attr = _rgb_hex(fill) if fill is not None else "none"
        op_attr = f' fill-opacity="{alpha:.2f}"' if (fill is not None and alpha < 0.999) else ""
        stroke_attr = (
            f' stroke="{_rgb_hex(stroke)}" stroke-width="{stroke_width}"'
            if stroke is not None and stroke_width > 0
            else ' stroke="none"'
        )
        self.svg_elements.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="{fill_attr}"{op_attr}{stroke_attr}/>'
        )

    def _raster_line(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        color: Color,
        width: int = 1,
        dashed: bool = False,
    ) -> None:
        ix0, iy0 = int(round(x0)), int(round(y0))
        ix1, iy1 = int(round(x1)), int(round(y1))
        dx = abs(ix1 - ix0)
        dy = abs(iy1 - iy0)
        sx = 1 if ix0 < ix1 else -1
        sy = 1 if iy0 < iy1 else -1
        err = dx - dy
        r_half = max(0, (width - 1) // 2)
        r_extra = 1 if (width % 2 == 0 and width >= 2) else 0
        step_count = 0

        while True:
            if not dashed or ((step_count // 6) % 2 == 0):
                if width <= 1:
                    self._set_px(ix0, iy0, color)
                else:
                    for oy in range(-r_half, r_half + 1 + r_extra):
                        for ox in range(-r_half, r_half + 1 + r_extra):
                            self._set_px(ix0 + ox, iy0 + oy, color)
            if ix0 == ix1 and iy0 == iy1:
                break
            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                ix0 += sx
            if e2 < dx:
                err += dx
                iy0 += sy
            step_count += 1

    def draw_line(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        color: Color,
        width: int = 2,
        dashed: bool = False,
    ) -> None:
        self._raster_line(x0, y0, x1, y1, color, width=width, dashed=dashed)
        dash_attr = ' stroke-dasharray="6,4"' if dashed else ""
        self.svg_elements.append(
            f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
            f'stroke="{_rgb_hex(color)}" stroke-width="{width}"{dash_attr} stroke-linecap="round"/>'
        )

    def draw_polyline(
        self,
        points: Sequence[Tuple[float, float]],
        color: Color,
        width: int = 2,
        dashed: bool = False,
    ) -> None:
        if len(points) < 2:
            return
        for i in range(len(points) - 1):
            self._raster_line(
                points[i][0], points[i][1], points[i + 1][0], points[i + 1][1],
                color=color, width=width, dashed=dashed,
            )
        pts_str = " ".join(f"{px:.2f},{py:.2f}" for px, py in points)
        dash_attr = ' stroke-dasharray="6,4"' if dashed else ""
        self.svg_elements.append(
            f'<polyline points="{pts_str}" fill="none" stroke="{_rgb_hex(color)}" '
            f'stroke-width="{width}"{dash_attr} stroke-linejoin="round" stroke-linecap="round"/>'
        )

    def draw_circle(
        self,
        cx: float,
        cy: float,
        radius: int,
        fill: Color,
        stroke: Optional[Color] = None,
    ) -> None:
        icx, icy = int(round(cx)), int(round(cy))
        r2 = radius * radius
        for oy in range(-radius, radius + 1):
            for ox in range(-radius, radius + 1):
                d2 = ox * ox + oy * oy
                if d2 <= r2:
                    if stroke is not None and d2 >= (radius - 1) * (radius - 1) and radius >= 2:
                        self._set_px(icx + ox, icy + oy, stroke)
                    else:
                        self._set_px(icx + ox, icy + oy, fill)
        stroke_attr = f' stroke="{_rgb_hex(stroke)}" stroke-width="1"' if stroke is not None else ""
        self.svg_elements.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{radius}" fill="{_rgb_hex(fill)}"{stroke_attr}/>'
        )

    def draw_cross(self, cx: float, cy: float, size: int, color: Color, width: int = 2) -> None:
        self.draw_line(cx - size, cy - size, cx + size, cy + size, color=color, width=width)
        self.draw_line(cx - size, cy + size, cx + size, cy - size, color=color, width=width)

    def draw_text(
        self,
        x: float,
        y: float,
        text: str,
        color: Color = PALETTE["dark"],
        scale: int = 2,
        anchor: str = "start",
        bold: bool = False,
    ) -> None:
        """Render crisp text at ``(x, y)`` where ``y`` is the top of the glyph box."""
        clean = "".join(ch if ch in FONT_5X7 else "?" for ch in text)
        char_adv = 6 * scale
        text_w = max(0, len(clean) * char_adv - scale)
        text_h = 7 * scale

        if anchor == "middle":
            rx = int(round(x - text_w / 2.0))
            svg_anchor = "middle"
            svg_x = x
        elif anchor == "end":
            rx = int(round(x - text_w))
            svg_anchor = "end"
            svg_x = x
        else:
            rx = int(round(x))
            svg_anchor = "start"
            svg_x = x

        ry = int(round(y))
        for idx, ch in enumerate(clean):
            glyph = FONT_5X7[ch]
            gx0 = rx + idx * char_adv
            for row_i, row_bits in enumerate(glyph):
                if row_bits == 0:
                    continue
                for col_i in range(5):
                    if (row_bits >> (4 - col_i)) & 1:
                        px0 = gx0 + col_i * scale
                        py0 = ry + row_i * scale
                        for sy in range(scale):
                            for sx in range(scale):
                                self._set_px(px0 + sx, py0 + sy, color)
                                if bold and sx == scale - 1:
                                    self._set_px(px0 + sx + 1, py0 + sy, color)

        font_px = max(10, int(round(text_h * 0.92)))
        weight = "bold" if bold else "normal"
        svg_y = y + text_h * 0.86
        self.svg_elements.append(
            f'<text x="{svg_x:.1f}" y="{svg_y:.1f}" fill="{_rgb_hex(color)}" '
            f'font-family="ui-monospace, SFMono-Regular, Menlo, Consolas, monospace" '
            f'font-size="{font_px}px" font-weight="{weight}" text-anchor="{svg_anchor}">'
            f"{html.escape(clean)}</text>"
        )

    def save(self, svg_path: str, png_path: str) -> None:
        """Write both the vector SVG and compressed RGB PNG files to disk."""
        os.makedirs(os.path.dirname(os.path.abspath(svg_path)), exist_ok=True)
        os.makedirs(os.path.dirname(os.path.abspath(png_path)), exist_ok=True)

        svg_doc = (
            f'<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n'
            + "\n".join(self.svg_elements)
            + "\n</svg>\n"
        )
        with open(svg_path, "w", encoding="utf-8") as f_svg:
            f_svg.write(svg_doc)

        png_bytes = encode_png_rgb(self.width, self.height, bytes(self.pixels))
        with open(png_path, "wb") as f_png:
            f_png.write(png_bytes)


def _nice_linear_ticks(vmin: float, vmax: float, n_ticks: int = 5) -> Tuple[float, float, List[float]]:
    if not math.isfinite(vmin) or not math.isfinite(vmax):
        return 0.0, 1.0, [0.0, 0.25, 0.5, 0.75, 1.0]
    if abs(vmax - vmin) < 1e-15:
        pad = 1.0 if abs(vmin) < 1e-12 else 0.1 * abs(vmin)
        vmin -= pad
        vmax += pad
    span = vmax - vmin
    raw_step = span / max(1, n_ticks - 1)
    mag = 10.0 ** math.floor(math.log10(raw_step))
    norm = raw_step / mag
    if norm <= 1.5:
        step = 1.0 * mag
    elif norm <= 3.0:
        step = 2.0 * mag
    elif norm <= 7.0:
        step = 5.0 * mag
    else:
        step = 10.0 * mag
    lo = math.floor(vmin / step + 1e-9) * step
    hi = math.ceil(vmax / step - 1e-9) * step
    ticks: List[float] = []
    cur = lo
    while cur <= hi + 0.5 * step:
        ticks.append(0.0 if abs(cur) < 1e-12 * step else cur)
        cur += step
    return lo, hi, ticks


def _fmt_tick(val: float) -> str:
    if abs(val) == 0.0:
        return "0"
    av = abs(val)
    if av >= 1e4 or av < 1e-3:
        return f"{val:.1e}"
    if abs(val - round(val)) < 1e-9:
        return str(int(round(val)))
    if av >= 10.0:
        return f"{val:.1f}"
    if av >= 1.0:
        return f"{val:.2f}".rstrip("0").rstrip(".")
    return f"{val:.3f}".rstrip("0").rstrip(".")


def draw_line_panel(
    fig: DualFigure,
    box: Tuple[int, int, int, int],
    title: str,
    xlabel: str,
    ylabel: str,
    series_list: Sequence[Dict[str, Any]],
    x_log10: bool = False,
    y_log10: bool = False,
    y_limits: Optional[Tuple[float, float]] = None,
    hlines: Optional[Sequence[Dict[str, Any]]] = None,
) -> None:
    """Render a publication-grade scientific line/scatter panel inside ``box = (x, y, w, h)``.

    Each item in ``series_list`` supports keys:
      - ``x``: sequence of floats
      - ``y``: sequence of floats
      - ``label``: legend string
      - ``color``: RGB tuple
      - ``width``: int (default 2)
      - ``dashed``: bool (default False)
      - ``markers``: bool (default False)
      - ``blowup_x``: optional float where series blew up (renders crimson X marker)
    """
    bx, by, bw, bh = box
    pad_l, pad_r, pad_t, pad_b = 76, 22, 42, 46
    px0 = bx + pad_l
    py0 = by + pad_t
    pw = max(40, bw - pad_l - pad_r)
    ph = max(40, bh - pad_t - pad_b)

    fig.fill_rect(bx, by, bw, bh, fill=PALETTE["panel_bg"], stroke=PALETTE["grid"], stroke_width=1)
    fig.draw_text(bx + 12, by + 10, title, color=PALETTE["navy"], scale=2, bold=True)
    fig.draw_text(bx + 12, by + 27, f"y: {ylabel}", color=PALETTE["slate"], scale=1)

    xs_all: List[float] = []
    ys_all: List[float] = []
    for s in series_list:
        for xv, yv in zip(s.get("x", []), s.get("y", [])):
            if xv is None or yv is None or not math.isfinite(xv) or not math.isfinite(yv):
                continue
            tx = math.log10(max(xv, 1e-30)) if x_log10 else float(xv)
            ty = math.log10(max(yv, 1e-30)) if y_log10 else float(yv)
            xs_all.append(tx)
            ys_all.append(ty)
        bx_val = s.get("blowup_x")
        if bx_val is not None and math.isfinite(bx_val):
            xs_all.append(math.log10(max(bx_val, 1e-30)) if x_log10 else float(bx_val))

    if hlines:
        for hl in hlines:
            hv = hl["y"]
            ys_all.append(math.log10(max(hv, 1e-30)) if y_log10 else float(hv))

    x_min_raw = min(xs_all) if xs_all else 0.0
    x_max_raw = max(xs_all) if xs_all else 1.0
    if y_limits is not None:
        y_min_raw = math.log10(max(y_limits[0], 1e-30)) if y_log10 else y_limits[0]
        y_max_raw = math.log10(max(y_limits[1], 1e-30)) if y_log10 else y_limits[1]
    else:
        y_min_raw = min(ys_all) if ys_all else 0.0
        y_max_raw = max(ys_all) if ys_all else 1.0

    xlo, xhi, xticks = _nice_linear_ticks(x_min_raw, x_max_raw, n_ticks=5)
    ylo, yhi, yticks = _nice_linear_ticks(y_min_raw, y_max_raw, n_ticks=5)
    if y_limits is not None:
        ylo, yhi = y_min_raw, y_max_raw

    def _map_x(xv: float) -> float:
        u = (xv - xlo) / (xhi - xlo) if xhi > xlo else 0.5
        return px0 + max(0.0, min(1.0, u)) * pw

    def _map_y(yv: float) -> float:
        u = (yv - ylo) / (yhi - ylo) if yhi > ylo else 0.5
        return py0 + ph - max(0.0, min(1.0, u)) * ph

    # Gridlines & ticks
    for xt in xticks:
        if xlo - 1e-9 <= xt <= xhi + 1e-9:
            sx = _map_x(xt)
            fig.draw_line(sx, py0, sx, py0 + ph, color=PALETTE["grid"], width=1)
            lbl = f"10^{_fmt_tick(xt)}" if x_log10 else _fmt_tick(xt)
            fig.draw_text(sx, py0 + ph + 7, lbl, color=PALETTE["slate"], scale=1, anchor="middle")

    for yt in yticks:
        if ylo - 1e-9 <= yt <= yhi + 1e-9:
            sy = _map_y(yt)
            fig.draw_line(px0, sy, px0 + pw, sy, color=PALETTE["grid"], width=1)
            lbl = f"10^{_fmt_tick(yt)}" if y_log10 else _fmt_tick(yt)
            fig.draw_text(px0 - 6, sy - 3, lbl, color=PALETTE["slate"], scale=1, anchor="end")

    # Reference horizontal lines
    if hlines:
        for hl in hlines:
            hv_raw = math.log10(max(hl["y"], 1e-30)) if y_log10 else float(hl["y"])
            if ylo <= hv_raw <= yhi:
                sy = _map_y(hv_raw)
                hc = hl.get("color", PALETTE["crimson"])
                fig.draw_line(px0, sy, px0 + pw, sy, color=hc, width=hl.get("width", 1), dashed=True)
                if hl.get("label"):
                    fig.draw_text(px0 + pw - 4, sy - 10, str(hl["label"]), color=hc, scale=1, anchor="end")

    # Plot series
    for s in series_list:
        color = s.get("color", PALETTE["blue"])
        lw = s.get("width", 2)
        dashed = s.get("dashed", False)
        pts: List[Tuple[float, float]] = []
        for xv, yv in zip(s.get("x", []), s.get("y", [])):
            if xv is None or yv is None or not math.isfinite(xv) or not math.isfinite(yv):
                if len(pts) >= 2:
                    fig.draw_polyline(pts, color=color, width=lw, dashed=dashed)
                pts = []
                continue
            tx = math.log10(max(xv, 1e-30)) if x_log10 else float(xv)
            ty = math.log10(max(yv, 1e-30)) if y_log10 else float(yv)
            pts.append((_map_x(tx), _map_y(ty)))

        if len(pts) >= 2:
            fig.draw_polyline(pts, color=color, width=lw, dashed=dashed)
        if s.get("markers", False):
            stride = max(1, len(pts) // 18) if len(pts) > 24 else 1
            for px, py in pts[::stride]:
                fig.draw_circle(px, py, radius=3, fill=color)
        bx_val = s.get("blowup_x")
        if bx_val is not None and math.isfinite(bx_val):
            tx = math.log10(max(bx_val, 1e-30)) if x_log10 else float(bx_val)
            sx = _map_x(tx)
            sy = pts[-1][1] if pts else (py0 + 10)
            fig.draw_cross(sx, sy, size=5, color=PALETTE["crimson"], width=2)

    # Axis border and X label
    fig.fill_rect(px0, py0, pw, ph, fill=None, stroke=PALETTE["slate"], stroke_width=1)
    fig.draw_text(px0 + pw / 2.0, by + bh - 15, xlabel, color=PALETTE["dark"], scale=1, anchor="middle")

    # Compact legend in top-right interior
    labeled = [s for s in series_list if s.get("label")]
    if labeled:
        max_chars = max(len(str(s["label"])) for s in labeled)
        leg_w = min(pw - 12, max_chars * 6 + 34)
        leg_h = len(labeled) * 13 + 8
        lx = px0 + pw - leg_w - 6
        ly = py0 + 6
        fig.fill_rect(lx, ly, leg_w, leg_h, fill=(252, 253, 255), stroke=PALETTE["grid"], stroke_width=1, alpha=0.90)
        for idx, s in enumerate(labeled):
            ry = ly + 6 + idx * 13
            c = s.get("color", PALETTE["blue"])
            fig.draw_line(lx + 6, ry + 3, lx + 22, ry + 3, color=c, width=2, dashed=s.get("dashed", False))
            if s.get("markers", False):
                fig.draw_circle(lx + 14, ry + 3, radius=2, fill=c)
            fig.draw_text(lx + 27, ry, str(s["label"]), color=PALETTE["dark"], scale=1)


def draw_heatmap_panel(
    fig: DualFigure,
    box: Tuple[int, int, int, int],
    title: str,
    subtitle: str,
    grid_2d: Sequence[Sequence[float]],
    diverging: bool = True,
    xlabel: str = "x",
    ylabel: str = "y",
) -> None:
    """Render a 2D lattice slice/field heatmap inside ``box = (x, y, w, h)``."""
    bx, by, bw, bh = box
    pad_l, pad_r, pad_t, pad_b = 58, 66, 42, 42
    px0 = bx + pad_l
    py0 = by + pad_t
    pw = max(40, bw - pad_l - pad_r)
    ph = max(40, bh - pad_t - pad_b)

    fig.fill_rect(bx, by, bw, bh, fill=PALETTE["panel_bg"], stroke=PALETTE["grid"], stroke_width=1)
    fig.draw_text(bx + 12, by + 10, title, color=PALETTE["navy"], scale=2, bold=True)
    fig.draw_text(bx + 12, by + 27, subtitle, color=PALETTE["slate"], scale=1)

    ny = len(grid_2d)
    nx = len(grid_2d[0]) if ny > 0 else 0
    if ny == 0 or nx == 0:
        return

    flat = [float(v) for row in grid_2d for v in row if math.isfinite(v)]
    vmin = min(flat) if flat else -1.0
    vmax = max(flat) if flat else 1.0
    vabs = max(abs(vmin), abs(vmax), 1e-12)

    cell_w = pw / float(nx)
    cell_h = ph / float(ny)

    for iy in range(ny):
        # iy = 0 at bottom of physical domain
        ry = py0 + ph - (iy + 1) * cell_h
        for ix in range(nx):
            rx = px0 + ix * cell_w
            val = float(grid_2d[iy][ix])
            if diverging:
                col = diverging_bwr_colormap(val / vabs)
            else:
                span = (vmax - vmin) if vmax > vmin else 1.0
                col = viridis_like_colormap((val - vmin) / span)
            fig.fill_rect(rx, ry, cell_w + 0.6, cell_h + 0.6, fill=col, stroke=None)

    fig.fill_rect(px0, py0, pw, ph, fill=None, stroke=PALETTE["slate"], stroke_width=1)
    fig.draw_text(px0 + pw / 2.0, by + bh - 15, xlabel, color=PALETTE["dark"], scale=1, anchor="middle")
    fig.draw_text(px0 - 8, py0 + ph / 2.0 - 3, ylabel, color=PALETTE["dark"], scale=1, anchor="end")

    # Colorbar on the right
    cb_x = px0 + pw + 10
    cb_w = 12
    n_cb = 32
    for i in range(n_cb):
        frac = i / float(n_cb - 1)
        cy = py0 + ph - (i + 1) * (ph / float(n_cb))
        col = diverging_bwr_colormap(2.0 * frac - 1.0) if diverging else viridis_like_colormap(frac)
        fig.fill_rect(cb_x, cy, cb_w, ph / float(n_cb) + 0.8, fill=col)
    fig.fill_rect(cb_x, py0, cb_w, ph, fill=None, stroke=PALETTE["slate"], stroke_width=1)

    top_lbl = _fmt_tick(vabs if diverging else vmax)
    mid_lbl = _fmt_tick(0.0 if diverging else 0.5 * (vmin + vmax))
    bot_lbl = _fmt_tick(-vabs if diverging else vmin)
    fig.draw_text(cb_x + cb_w + 4, py0, top_lbl, color=PALETTE["slate"], scale=1)
    fig.draw_text(cb_x + cb_w + 4, py0 + ph / 2.0 - 3, mid_lbl, color=PALETTE["slate"], scale=1)
    fig.draw_text(cb_x + cb_w + 4, py0 + ph - 7, bot_lbl, color=PALETTE["slate"], scale=1)


# --- Component: render_plots.py ---
def _make_2x2_boxes(width: int = 1360, height: int = 920) -> Tuple[Tuple[int, int, int, int], ...]:
    top_y = 56
    margin_x = 20
    gap_x = 20
    gap_y = 18
    bottom_pad = 16
    pw = (width - 2 * margin_x - gap_x) // 2
    ph = (height - top_y - bottom_pad - gap_y) // 2
    return (
        (margin_x, top_y, pw, ph),
        (margin_x + pw + gap_x, top_y, pw, ph),
        (margin_x, top_y + ph + gap_y, pw, ph),
        (margin_x + pw + gap_x, top_y + ph + gap_y, pw, ph),
    )


def _draw_figure_banner(fig: DualFigure, title: str, subtitle: str) -> None:
    fig.fill_rect(0, 0, fig.width, 46, fill=PALETTE["navy"])
    fig.draw_text(20, 9, title, color=(255, 255, 255), scale=2, bold=True)
    fig.draw_text(20, 29, subtitle, color=(205, 220, 240), scale=1)


def render_plot1_battle1(data: Dict[str, Any], plots_dir: str) -> Tuple[str, str]:
    q_runs = data["q_runs"]
    f64 = data["float64_and_convergence"]

    fig = DualFigure(width=1360, height=920)
    _draw_figure_banner(
        fig,
        "Battle 1: 2D Gaussian Torus Z[i]/<(1+i)^m> Turbulence, Thermodynamics & Bit Plateau",
        "Exact 0/1 (Gamma, E_kin, Z) in 1A | Exact 0/1 First Law & SOS dQ >= 0 in 1B | Zero FFTs/Roots",
    )
    b_tl, b_tr, b_bl, b_br = _make_2x2_boxes()

    # Panel A: 1A Energy Conservation vs Coarse Continuum Baselines (m=6, dt=1/4)
    r_1a = q_runs["B1A_m6_dt1_4"]
    e0 = r_1a["time_series"][0]["E_kin"]
    s_panel_a: List[Dict[str, Any]] = [
        {
            "x": [pt["t"] for pt in r_1a["time_series"]],
            "y": [pt["E_kin"] for pt in r_1a["time_series"]],
            "label": "Q-Engine 1A (0/1 exact)",
            "color": PALETTE["green"],
            "width": 3,
        }
    ]
    coarse_b1 = f64["coarse_baseline_comparisons"]["B1A_m6"]["1/4"]
    for meth, lbl, col, dash in [
        ("forward_euler", "Float64 Euler (dt=1/4)", PALETTE["crimson"], True),
        ("midpoint_rk2", "Float64 RK2 (dt=1/4)", PALETTE["orange"], True),
        ("classical_rk4", "Float64 RK4 (dt=1/4)", PALETTE["blue"], False),
    ]:
        c_item = coarse_b1[meth]
        xs = [h["t"] for h in c_item["history"] if h["E_kin"] <= 5.0 * e0]
        ys = [h["E_kin"] for h in c_item["history"] if h["E_kin"] <= 5.0 * e0]
        bx_t = (c_item["blowup_step"] * 0.25) if c_item["blew_up"] and c_item["blowup_step"] else None
        if xs and ys and ys[-1] >= 2.5 * e0 and bx_t is None:
            bx_t = xs[-1]
        s_panel_a.append({
            "x": xs,
            "y": ys,
            "label": lbl + (" [BLOWUP]" if c_item["blew_up"] else ""),
            "color": col,
            "width": 2,
            "dashed": dash,
            "blowup_x": min(bx_t, 75.0) if bx_t is not None else None,
        })

    draw_line_panel(
        fig,
        b_tl,
        title="A. Battle 1A Energy E_kin (N=64, dt=1/4)",
        xlabel="Physical time t = n * dt",
        ylabel="Kinetic energy E_kin(t)",
        series_list=s_panel_a,
        y_limits=(0.0, max(2.2 * e0, 1.0)),
    )

    # Panel B: Battle 1B Viscous First & Second Laws (N=64, nu_4=1/200, dt=1/4)
    r_1b = q_runs["B1B_m6_dt1_4"]
    ts_1b = r_1b["time_series"]
    draw_line_panel(
        fig,
        b_tr,
        title="B. Battle 1B Thermodynamics (nu_4=1/200, dt=1/4)",
        xlabel="Physical time t = n * dt",
        ylabel="Energy / Thermal Reservoir",
        series_list=[
            {
                "x": [p["t"] for p in ts_1b],
                "y": [p["E_total"] for p in ts_1b],
                "label": "E_kin + Q_heat (0/1 == E_0)",
                "color": PALETTE["navy"],
                "width": 3,
            },
            {
                "x": [p["t"] for p in ts_1b],
                "y": [p["E_kin"] for p in ts_1b],
                "label": "E_kin(t) (dissipating)",
                "color": PALETTE["blue"],
                "width": 2,
            },
            {
                "x": [p["t"] for p in ts_1b],
                "y": [p["Q_heat"] for p in ts_1b],
                "label": "Q_heat(t) (dQ >= 0/1)",
                "color": PALETTE["vermilion"],
                "width": 2,
            },
        ],
    )

    # Panel C: Intrinsic (1+i)-adic Roughness-Shell Inverse Transfer
    r_1a_m7 = q_runs["B1A_m7_dt1_4"]
    draw_line_panel(
        fig,
        b_bl,
        title="C. 2D Inverse Transfer to IR Shells (R in {1,2})",
        xlabel="Physical time t = n * dt",
        ylabel="IR / UV Energy Fraction of E_0",
        series_list=[
            {
                "x": [p["t"] for p in r_1a["time_series"]],
                "y": [p["ir_frac"] for p in r_1a["time_series"]],
                "label": "N=64 IR Shell Frac (R=1,2)",
                "color": PALETTE["green"],
                "width": 2,
                "markers": True,
            },
            {
                "x": [p["t"] for p in r_1a["time_series"]],
                "y": [p["uv_frac"] for p in r_1a["time_series"]],
                "label": "N=64 UV Shell Frac (R>=3)",
                "color": PALETTE["purple"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": [p["t"] for p in r_1a_m7["time_series"]],
                "y": [p["ir_frac"] for p in r_1a_m7["time_series"]],
                "label": "N=128 IR Shell Frac (R=1,2)",
                "color": PALETTE["blue"],
                "width": 2,
            },
        ],
        y_limits=(0.0, 1.05),
    )

    # Panel D: Bit-length plateaus across N=64, N=128, dt=1/2 or Rev, and Prime-Perturbed (+3/65537)
    r_1a_p = q_runs["B1A_m6_dt1_4_prime"]
    if "B1A_rev_m6_dt1_4" in q_runs:
        r_1a_alt = q_runs["B1A_rev_m6_dt1_4"]
        alt_lbl_1 = f"1A Rev N=64 (max={r_1a_alt['bits']['peak_bits']}b)"
    else:
        r_1a_alt = q_runs["B1A_m6_dt1_2"]
        alt_lbl_1 = f"1A dt=1/2 N=64 (max={r_1a_alt['bits']['peak_bits']}b)"
    draw_line_panel(
        fig,
        b_br,
        title="D. Ultrafinitist Bit Plateau (n=300 Steps)",
        xlabel="Macro-step n",
        ylabel="Peak numerator/denominator bits",
        series_list=[
            {
                "x": [p["step"] for p in r_1a["time_series"]],
                "y": [p["max_bits"] for p in r_1a["time_series"]],
                "label": f"1A N=64 dt=1/4 (max={r_1a['bits']['peak_bits']}b)",
                "color": PALETTE["green"],
                "width": 2,
            },
            {
                "x": [p["step"] for p in r_1a_m7["time_series"]],
                "y": [p["max_bits"] for p in r_1a_m7["time_series"]],
                "label": f"1A N=128 dt=1/4 (max={r_1a_m7['bits']['peak_bits']}b)",
                "color": PALETTE["blue"],
                "width": 2,
            },
            {
                "x": [p["step"] for p in r_1a_alt["time_series"]],
                "y": [p["max_bits"] for p in r_1a_alt["time_series"]],
                "label": alt_lbl_1,
                "color": PALETTE["orange"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": [p["step"] for p in r_1a_p["time_series"]],
                "y": [p["max_bits"] for p in r_1a_p["time_series"]],
                "label": f"1A Prime +3/65537 (max={r_1a_p['bits']['peak_bits']}b)",
                "color": PALETTE["purple"],
                "width": 2,
            },
        ],
        y_limits=(
            0.0,
            max(
                10800.0,
                1.08 * max(
                    r_1a["bits"]["peak_bits"],
                    r_1a_m7["bits"]["peak_bits"],
                    r_1a_alt["bits"]["peak_bits"],
                    r_1a_p["bits"]["peak_bits"],
                ),
            ),
        ),
        hlines=[{"y": 10000.0, "color": PALETTE["crimson"], "width": 2, "label": "10,000-bit Gate 2 Ceiling"}],
    )

    svg_p = os.path.join(plots_dir, "battle1_2d_turbulence_and_bits.svg")
    png_p = os.path.join(plots_dir, "battle1_2d_turbulence_and_bits.png")
    fig.save(svg_p, png_p)
    return svg_p, png_p


def render_plot2_battle2(data: Dict[str, Any], plots_dir: str) -> Tuple[str, str]:
    q_runs = data["q_runs"]
    f64 = data["float64_and_convergence"]

    fig = DualFigure(width=1360, height=920)
    _draw_figure_banner(
        fig,
        "Battle 2: 3D Taylor-Green Vortex Stretching on Z_L^3, UV Cascade & O(dt^2) Accuracy",
        "Exact 0/1 div_Sym(u) == 0, Momentum P, Energy E_kin | Active 3D Stretching Z_max/Z_0 >= 1.35",
    )
    b_tl, b_tr, b_bl, b_br = _make_2x2_boxes()

    # Panel A: Enstrophy Amplification Z(t)/Z_0 vs RK4(dt=1/512)
    ref_l6 = f64["battle2_rk4_ref_stretching"]["L=6"]["history"]
    ref_l8 = f64["battle2_rk4_ref_stretching"]["L=8"]["history"]
    r_l6_d4 = q_runs["B2A_L6_dt1_4"]
    r_l6_d8 = q_runs["B2A_L6_dt1_8"]
    r_l8_d8 = q_runs.get("B2A_L8_dt1_8", q_runs.get("B2A_L8_dt1_8_T2", r_l6_d8))

    draw_line_panel(
        fig,
        b_tl,
        title="A. 3D Vortex Stretching Z(t)/Z_0 vs RK4(dt=1/512)",
        xlabel="Physical time t (early stretching window)",
        ylabel="Enstrophy Amplification Z(t)/Z_0",
        series_list=[
            {
                "x": [h["t"] for h in ref_l8 if h["t"] <= 2.0],
                "y": [h["Z_ratio"] for h in ref_l8 if h["t"] <= 2.0],
                "label": "L=8 RK4(dt=1/512) Ref",
                "color": PALETTE["navy"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": [p["t"] for p in r_l8_d8["time_series"] if p["t"] <= 2.0],
                "y": [p["Z_ratio"] for p in r_l8_d8["time_series"] if p["t"] <= 2.0],
                "label": f"L=8 Q-Engine dt=1/8 (pk={r_l8_d8['stretching']['Z_max_over_Z0']:.3f})",
                "color": PALETTE["blue"],
                "width": 2,
                "markers": True,
            },
            {
                "x": [h["t"] for h in ref_l6 if h["t"] <= 2.5],
                "y": [h["Z_ratio"] for h in ref_l6 if h["t"] <= 2.5],
                "label": "L=6 RK4(dt=1/512) Ref",
                "color": PALETTE["slate"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": [p["t"] for p in r_l6_d8["time_series"] if p["t"] <= 2.5],
                "y": [p["Z_ratio"] for p in r_l6_d8["time_series"] if p["t"] <= 2.5],
                "label": f"L=6 Q-Engine dt=1/8 (pk={r_l6_d8['stretching']['Z_max_over_Z0']:.3f})",
                "color": PALETTE["green"],
                "width": 2,
                "markers": True,
            },
            {
                "x": [p["t"] for p in r_l6_d4["time_series"] if p["t"] <= 2.5],
                "y": [p["Z_ratio"] for p in r_l6_d4["time_series"] if p["t"] <= 2.5],
                "label": "L=6 Q-Engine dt=1/4",
                "color": PALETTE["orange"],
                "width": 2,
            },
        ],
        hlines=[{"y": 1.35, "color": PALETTE["crimson"], "width": 2, "label": "Z_max/Z_0 >= 1.35 Gate"}],
    )

    # Panel B: Q-Engine O(dt^2) Convergence vs Continuum RK Baselines on Battle 2A (L=6)
    q_conv = f64["q_engine_convergence"]["B2A_L6"]
    g4_b2 = f64["gate4_slope_verification"]["B2A_L6"]["convergence"]["methods"]
    draw_line_panel(
        fig,
        b_tr,
        title="B. Battle 2A Convergence Orders vs RK4(dt=1/512)",
        xlabel="log10(dt)",
        ylabel="log10(L2 RMS Error)",
        series_list=[
            {
                "x": q_conv["dts_float"],
                "y": q_conv["T1"]["u_l2_rms_errors"],
                "label": "Q-Engine u(T=1) O(dt^2)",
                "color": PALETTE["green"],
                "width": 3,
                "markers": True,
            },
            {
                "x": q_conv["dts_float"],
                "y": q_conv["T2"]["u_l2_rms_errors"],
                "label": "Q-Engine u(T=2) O(dt^2)",
                "color": PALETTE["sky"],
                "width": 2,
                "markers": True,
            },
            {
                "x": g4_b2["forward_euler"]["dts"],
                "y": g4_b2["forward_euler"]["l2_rms_errors"],
                "label": f"Euler (p={g4_b2['forward_euler']['asymptotic_order_l2']:.2f})",
                "color": PALETTE["crimson"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": g4_b2["midpoint_rk2"]["dts"],
                "y": g4_b2["midpoint_rk2"]["l2_rms_errors"],
                "label": f"RK2 (p={g4_b2['midpoint_rk2']['asymptotic_order_l2']:.2f})",
                "color": PALETTE["orange"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": g4_b2["classical_rk4"]["dts"],
                "y": g4_b2["classical_rk4"]["l2_rms_errors"],
                "label": f"RK4 (p={g4_b2['classical_rk4']['asymptotic_order_l2']:.2f})",
                "color": PALETTE["blue"],
                "width": 2,
                "dashed": True,
            },
        ],
        x_log10=True,
        y_log10=True,
    )

    # Panel C: 3D Roughness-Shell Forward Energy Cascade E(tilde_R)
    sp = r_l6_d8["spectral_transfer"]
    shells_sorted = sorted(int(k) for k in sp["initial_shells"].keys())
    e_init = [sp["initial_shells"][str(k)] for k in shells_sorted]
    e_peak = [sp["peak_z_shells"][str(k)] for k in shells_sorted]
    e_fin = [sp["final_shells"][str(k)] for k in shells_sorted]
    draw_line_panel(
        fig,
        b_bl,
        title="C. 3D Forward Cascade Across Shells tilde_R (L=6)",
        xlabel="Intrinsic Roughness Shell tilde_R = r_x + r_y + r_z",
        ylabel="Shell Kinetic Energy E(tilde_R)",
        series_list=[
            {"x": shells_sorted, "y": e_init, "label": "Initial t=0 (pure IR)", "color": PALETTE["navy"], "width": 2, "markers": True},
            {"x": shells_sorted, "y": e_peak, "label": "Peak Stretching (UV excited)", "color": PALETTE["vermilion"], "width": 2, "markers": True},
            {"x": shells_sorted, "y": e_fin, "label": "Final step n=250", "color": PALETTE["green"], "width": 2, "dashed": True, "markers": True},
        ],
    )

    # Panel D: Battle 2 Bit-Length Plateaus (2A, 2B Viscous, dt=1/2 or Rev, Prime)
    r_2b = q_runs["B2B_L6_dt1_4"]
    r_2p = q_runs["B2A_L6_dt1_4_prime"]
    if "B2A_rev_L6_dt1_4" in q_runs:
        r_2alt = q_runs["B2A_rev_L6_dt1_4"]
        alt_lbl_2 = f"2A Rev (max={r_2alt['bits']['peak_bits']}b)"
    else:
        r_2alt = q_runs["B2A_L6_dt1_2"]
        alt_lbl_2 = f"2A dt=1/2 (max={r_2alt['bits']['peak_bits']}b)"
    draw_line_panel(
        fig,
        b_br,
        title="D. Battle 2 3D Bit-Length Plateau (L=6, n=250 Steps)",
        xlabel="Macro-step n",
        ylabel="Peak numerator/denominator bits",
        series_list=[
            {
                "x": [p["step"] for p in r_l6_d4["time_series"]],
                "y": [p["max_bits"] for p in r_l6_d4["time_series"]],
                "label": f"2A dt=1/4 (max={r_l6_d4['bits']['peak_bits']}b)",
                "color": PALETTE["green"],
                "width": 2,
            },
            {
                "x": [p["step"] for p in r_2alt["time_series"]],
                "y": [p["max_bits"] for p in r_2alt["time_series"]],
                "label": alt_lbl_2,
                "color": PALETTE["orange"],
                "width": 2,
                "dashed": True,
            },
            {
                "x": [p["step"] for p in r_2b["time_series"]],
                "y": [p["max_bits"] for p in r_2b["time_series"]],
                "label": f"2B Viscous nu=1/50 (max={r_2b['bits']['peak_bits']}b)",
                "color": PALETTE["blue"],
                "width": 2,
            },
            {
                "x": [p["step"] for p in r_2p["time_series"]],
                "y": [p["max_bits"] for p in r_2p["time_series"]],
                "label": f"2A Prime +3/65537 (max={r_2p['bits']['peak_bits']}b)",
                "color": PALETTE["purple"],
                "width": 2,
            },
        ],
        y_limits=(
            0.0,
            max(
                10800.0,
                1.08 * max(
                    r_l6_d4["bits"]["peak_bits"],
                    r_2alt["bits"]["peak_bits"],
                    r_2b["bits"]["peak_bits"],
                    r_2p["bits"]["peak_bits"],
                ),
            ),
        ),
        hlines=[{"y": 10000.0, "color": PALETTE["crimson"], "width": 2, "label": "10,000-bit Ceiling"}],
    )

    svg_p = os.path.join(plots_dir, "battle2_3d_taylor_green_stretching.svg")
    png_p = os.path.join(plots_dir, "battle2_3d_taylor_green_stretching.png")
    fig.save(svg_p, png_p)
    return svg_p, png_p


def render_plot3_battle3(data: Dict[str, Any], plots_dir: str) -> Tuple[str, str]:
    q_runs = data["q_runs"]
    f64 = data["float64_and_convergence"]

    fig = DualFigure(width=1360, height=920)
    _draw_figure_banner(
        fig,
        "Battle 3: Extreme 3D Anti-Parallel Vortex Collision, No-Blowup & Positive Density",
        "Explicit Euler/RK2/RK4 Blow Up | Q-Engine Bounded by B_u, B_omega | Varga-SOS min rho > 0/1",
    )
    b_tl, b_tr, b_bl, b_br = _make_2x2_boxes()

    r3_d2 = q_runs["B3_L6_dt1_2"]
    r3_d4 = q_runs["B3_L6_dt1_4"]
    r3_p = q_runs["B3_L6_dt1_4_prime"]
    b_u_val = r3_d4["positivity_and_bounds"]["B_u_float"]
    b_om_val = r3_d4["positivity_and_bounds"]["B_omega_float"]

    # Panel A: Velocity Bound & Coarse Continuum Blowup
    coarse_b3 = f64["coarse_baseline_comparisons"]["B3_L6_centered"]["1/4"]
    s_vel: List[Dict[str, Any]] = [
        {
            "x": [p["t"] for p in r3_d4["time_series"]],
            "y": [p["max_u_sq"] for p in r3_d4["time_series"]],
            "label": "Q-Engine dt=1/4 max ||u||^2",
            "color": PALETTE["green"],
            "width": 2,
        },
        {
            "x": [p["t"] for p in r3_d2["time_series"]],
            "y": [p["max_u_sq"] for p in r3_d2["time_series"]],
            "label": "Q-Engine dt=1/2 max ||u||^2",
            "color": PALETTE["blue"],
            "width": 2,
            "dashed": True,
        },
    ]
    for meth, lbl, col in [
        ("forward_euler", "Euler dt=1/4", PALETTE["crimson"]),
        ("midpoint_rk2", "RK2 dt=1/4", PALETTE["orange"]),
        ("classical_rk4", "RK4 dt=1/4", PALETTE["purple"]),
    ]:
        cm = coarse_b3[meth]
        hist_clamped = [(h["t"], min(h["max_u_sq"], 1.15 * b_u_val)) for h in cm["history"]]
        bx_t = (cm["blowup_step"] * 0.25) if cm["blew_up"] and cm["blowup_step"] else None
        s_vel.append({
            "x": [pt[0] for pt in hist_clamped],
            "y": [pt[1] for pt in hist_clamped],
            "label": f"{lbl} (blowup step {cm['blowup_step']})",
            "color": col,
            "width": 2,
            "blowup_x": bx_t,
        })

    draw_line_panel(
        fig,
        b_tl,
        title="A. Peak ||u(x)||^2 vs Bound B_u(E_0, N)",
        xlabel="Physical time t = n * dt",
        ylabel="Sitewise max ||u_n(x)||^2",
        series_list=s_vel,
        y_limits=(0.0, 1.22 * b_u_val),
        hlines=[{"y": b_u_val, "color": PALETTE["crimson"], "width": 2, "label": f"A Priori B_u = {b_u_val:.1f}"}],
    )

    # Panel B: Vorticity Bound B_omega(E_0, N)
    draw_line_panel(
        fig,
        b_tr,
        title="B. Peak ||omega(x)||^2 vs Bound B_omega(E_0, N)",
        xlabel="Physical time t = n * dt",
        ylabel="Sitewise max ||omega_n(x)||^2",
        series_list=[
            {
                "x": [p["t"] for p in r3_d4["time_series"]],
                "y": [p["max_omega_sq"] for p in r3_d4["time_series"]],
                "label": "Q-Engine dt=1/4 max ||omega||^2",
                "color": PALETTE["green"],
                "width": 2,
            },
            {
                "x": [p["t"] for p in r3_d2["time_series"]],
                "y": [p["max_omega_sq"] for p in r3_d2["time_series"]],
                "label": "Q-Engine dt=1/2 max ||omega||^2",
                "color": PALETTE["blue"],
                "width": 2,
                "dashed": True,
            },
        ],
        y_limits=(0.0, 1.22 * b_om_val),
        hlines=[{"y": b_om_val, "color": PALETTE["crimson"], "width": 2, "label": f"A Priori B_omega = {b_om_val:.1f}"}],
    )

    # Panel C: Density Positivity min_x rho_n(x) > 0/1 vs Continuum Centered Gibbs Ringing (min rho < 0)
    gibbs = f64["battle3_centered_gibbs_ringing"]
    draw_line_panel(
        fig,
        b_bl,
        title="C. Density min rho > 0/1 vs Centered Gibbs (< 0)",
        xlabel="Early Collision Time t in [0, 1.0]",
        ylabel="Sitewise minimum density min_x rho(x)",
        series_list=[
            {
                "x": [p["t"] for p in r3_d4["time_series"] if p["t"] <= 1.0],
                "y": [p["min_rho"] for p in r3_d4["time_series"] if p["t"] <= 1.0],
                "label": f"Q-Engine Varga dt=1/4 (min={r3_d4['positivity_and_bounds']['min_rho_float']:.4f}>0)",
                "color": PALETTE["green"],
                "width": 3,
                "markers": True,
            },
            {
                "x": [p["t"] for p in r3_d2["time_series"] if p["t"] <= 1.0],
                "y": [p["min_rho"] for p in r3_d2["time_series"] if p["t"] <= 1.0],
                "label": f"Q-Engine Varga dt=1/2 (min={r3_d2['positivity_and_bounds']['min_rho_float']:.4f}>0)",
                "color": PALETTE["blue"],
                "width": 2,
                "markers": True,
            },
            {
                "x": [h["t"] for h in gibbs["history"]],
                "y": [h["min_rho"] for h in gibbs["history"]],
                "label": f"Centered RK4(dt=1/512) (min={gibbs['min_rho_continuum_centered']:.3f}<0)",
                "color": PALETTE["crimson"],
                "width": 2,
            },
        ],
        hlines=[{"y": 0.0, "color": PALETTE["dark"], "width": 2, "label": "Vacuum Positivity Floor rho = 0"}],
    )

    # Panel D: Battle 3 Bit-Length Plateau (Standard & Prime +1/10007)
    draw_line_panel(
        fig,
        b_br,
        title="D. Battle 3 Coupled (u, omega, rho) Bit Plateau",
        xlabel="Macro-step n",
        ylabel="Peak numerator/denominator bits",
        series_list=[
            {
                "x": [p["step"] for p in r3_d2["time_series"]],
                "y": [p["max_bits"] for p in r3_d2["time_series"]],
                "label": f"dt=1/2 (max={r3_d2['bits']['peak_bits']}b)",
                "color": PALETTE["blue"],
                "width": 2,
            },
            {
                "x": [p["step"] for p in r3_d4["time_series"]],
                "y": [p["max_bits"] for p in r3_d4["time_series"]],
                "label": f"dt=1/4 (max={r3_d4['bits']['peak_bits']}b)",
                "color": PALETTE["green"],
                "width": 2,
            },
            {
                "x": [p["step"] for p in r3_p["time_series"]],
                "y": [p["max_bits"] for p in r3_p["time_series"]],
                "label": f"Prime +1/10007 (max={r3_p['bits']['peak_bits']}b)",
                "color": PALETTE["purple"],
                "width": 2,
            },
        ],
        y_limits=(
            0.0,
            max(
                10800.0,
                1.08 * max(
                    r3_d2["bits"]["peak_bits"],
                    r3_d4["bits"]["peak_bits"],
                    r3_p["bits"]["peak_bits"],
                ),
            ),
        ),
        hlines=[{"y": 10000.0, "color": PALETTE["crimson"], "width": 2, "label": "10,000-bit Ceiling"}],
    )

    svg_p = os.path.join(plots_dir, "battle3_vortex_collision_positivity_noblowup.svg")
    png_p = os.path.join(plots_dir, "battle3_vortex_collision_positivity_noblowup.png")
    fig.save(svg_p, png_p)
    return svg_p, png_p


def render_plot4_summary(data: Dict[str, Any], plots_dir: str) -> Tuple[str, str]:
    q_runs = data["q_runs"]
    f64 = data["float64_and_convergence"]

    fig = DualFigure(width=1360, height=920)
    _draw_figure_banner(
        fig,
        "Cross-Battle Verification Summary: Small-dt Tangent, Slopes, Speed & Vortex State",
        "Zero Small-dt Quantization Freezing (down to 2^-30) | Fair-Fight Slopes | <5s SIGALRM Audit",
    )
    b_tl, b_tr, b_bl, b_br = _make_2x2_boxes()

    # Panel A: Small-dt tangent error down to dt = 2^-30
    tang = f64["small_dt_tangents"]
    s_tan: List[Dict[str, Any]] = []
    for key, lbl, col in [
        ("B1A_m6_step_battle1", "Battle 1A (m=6)", PALETTE["green"]),
        ("B1B_m6_step_battle1", "Battle 1B (m=6)", PALETTE["orange"]),
        ("B2A_L6_step_battle2", "Battle 2A (L=6)", PALETTE["blue"]),
        ("B3_L6_step_battle3", "Battle 3 (L=6)", PALETTE["purple"]),
    ]:
        recs = tang[key]["records"]
        s_tan.append({
            "x": [2.0 ** (-r["exponent"]) for r in recs],
            "y": [max(r["tangent_error_linf"], 1e-22) for r in recs],
            "label": lbl,
            "color": col,
            "width": 2,
            "markers": True,
        })
    draw_line_panel(
        fig,
        b_tl,
        title="A. Small-dt Tangent ||(F_dt-id)/dt - X(z_0)||_inf",
        xlabel="log10(dt) down to dt = 2^-30 (~9.3e-10)",
        ylabel="log10(Tangent Error)",
        series_list=s_tan,
        x_log10=True,
        y_log10=True,
    )

    # Panel B: Fair-Fight Continuum Slope Verification (Battle 1A m=6 & Battle 3 L=6)
    g4_b1 = f64["gate4_slope_verification"]["B1A_m6"]["convergence"]["methods"]
    q_c1 = f64["q_engine_convergence"]["B1A_m6"]
    draw_line_panel(
        fig,
        b_tr,
        title="B. Gate 4 Slopes & Q-Engine O(dt^2) (Battle 1A)",
        xlabel="log10(dt)",
        ylabel="log10(L2 RMS Error at T=0.5 / T=1.0)",
        series_list=[
            {
                "x": q_c1["dts_float"],
                "y": q_c1["omega_l2_rms_errors"],
                "label": f"Q-Engine 1A (p={q_c1['orders'][-1]:.2f})",
                "color": PALETTE["green"],
                "width": 3,
                "markers": True,
            },
            {
                "x": g4_b1["forward_euler"]["dts"],
                "y": g4_b1["forward_euler"]["l2_rms_errors"],
                "label": f"Euler (p={g4_b1['forward_euler']['asymptotic_order_l2']:.3f})",
                "color": PALETTE["crimson"],
                "width": 2,
                "dashed": True,
                "markers": True,
            },
            {
                "x": g4_b1["midpoint_rk2"]["dts"],
                "y": g4_b1["midpoint_rk2"]["l2_rms_errors"],
                "label": f"Midpoint RK2 (p={g4_b1['midpoint_rk2']['asymptotic_order_l2']:.3f})",
                "color": PALETTE["orange"],
                "width": 2,
                "dashed": True,
                "markers": True,
            },
            {
                "x": g4_b1["classical_rk4"]["dts"],
                "y": g4_b1["classical_rk4"]["l2_rms_errors"],
                "label": f"Classical RK4 (p={g4_b1['classical_rk4']['asymptotic_order_l2']:.3f})",
                "color": PALETTE["blue"],
                "width": 2,
                "dashed": True,
                "markers": True,
            },
        ],
        x_log10=True,
        y_log10=True,
    )

    # Panel C: Per-Step Wall-Clock Speed Across Representative Configurations vs 5.0s SIGALRM
    b2_alt_tag = "B2A_rev_L6_dt1_4" if "B2A_rev_L6_dt1_4" in q_runs else "B2A_L6_dt1_2"
    b2_l8_tag = "B2A_L8_dt1_8" if "B2A_L8_dt1_8" in q_runs else "B2A_L8_dt1_8_T2"
    rep_tags = [
        t
        for t in [
            "B1A_m6_dt1_4",
            "B1A_m7_dt1_4",
            "B1B_m6_dt1_4",
            "B2A_L6_dt1_4",
            b2_alt_tag,
            "B2B_L6_dt1_4",
            b2_l8_tag,
            "B3_L6_dt1_4",
        ]
        if t in q_runs
    ]
    xs_idx = list(range(1, len(rep_tags) + 1))
    ys_max = [q_runs[t]["max_step_sec"] for t in rep_tags]
    ys_mean = [q_runs[t]["mean_step_sec"] for t in rep_tags]
    draw_line_panel(
        fig,
        b_bl,
        title="C. Macro-Step Wall-Clock Speed vs 5.0s SIGALRM Limit",
        xlabel="Config Index (1:B1A6, 2:B1A7, 3:B1B6, 4:B2A6, 5:B2Rev6, 6:B2B6, 7:B2A8, 8:B3)",
        ylabel="Seconds per macro-step",
        series_list=[
            {"x": xs_idx, "y": ys_max, "label": "Peak Step Time (s)", "color": PALETTE["vermilion"], "width": 2, "markers": True},
            {"x": xs_idx, "y": ys_mean, "label": "Mean Step Time (s)", "color": PALETTE["blue"], "width": 2, "markers": True},
        ],
        y_limits=(0.0, 5.5),
        hlines=[{"y": 5.0, "color": PALETTE["crimson"], "width": 2, "label": "5.0s SIGALRM Hard Limit"}],
    )

    # Panel D: 2D Gaussian Torus (m=6, 8x8) Final Vorticity Snapshot
    om_fin = q_runs["B1A_m6_dt1_4"]["final_omega_float"]
    grid_8x8 = [[om_fin[iy * 8 + ix] for ix in range(8)] for iy in range(8)]
    draw_heatmap_panel(
        fig,
        b_br,
        title="D. Battle 1A Final Vorticity omega_300(x, y) on Z_8^2",
        subtitle="Merged IR dipole structure at step n=300 (dt=1/4, t=75)",
        grid_2d=grid_8x8,
        diverging=True,
        xlabel="Lattice x in 0..7",
        ylabel="y",
    )

    svg_p = os.path.join(plots_dir, "convergence_and_ultrafinitist_summary.svg")
    png_p = os.path.join(plots_dir, "convergence_and_ultrafinitist_summary.png")
    fig.save(svg_p, png_p)
    return svg_p, png_p


def render_all_plots_from_data(data: Dict[str, Any], plots_dir: str) -> List[Tuple[str, str]]:
    os.makedirs(plots_dir, exist_ok=True)
    return [
        render_plot1_battle1(data, plots_dir),
        render_plot2_battle2(data, plots_dir),
        render_plot3_battle3(data, plots_dir),
        render_plot4_summary(data, plots_dir),
    ]


def render_all_plots(results_json_path: str, plots_dir: str) -> List[Tuple[str, str]]:
    with open(results_json_path, "r", encoding="utf-8") as fp:
        data = json.load(fp)
    return render_all_plots_from_data(data, plots_dir)


def _standalone_entrypoint() -> None:
    root_dir = os.path.dirname(os.path.abspath(__file__))
    plots_dir = os.path.join(root_dir, 'plots')
    run_standalone_verification_and_plots(plots_dir)


if __name__ == '__main__':
    _standalone_entrypoint()
