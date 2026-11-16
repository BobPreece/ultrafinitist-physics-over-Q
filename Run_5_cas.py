"""Zero-dependency pure-Q symbolic computer algebra system (CAS).

Implements strictly over Q (Python int and fractions.Fraction only):
- Univariate polynomials QQPoly and rational functions QQRat,
- Multivariate polynomials QQMVPoly and rational functions QQMVRat,
- Truncated Taylor jets QQJet of arbitrary truncation order over Q or symbolic rings,
- Exact signed subresultant Sturm chains (sturm_chain) and bisection root isolation,
- Exact square-root-free Cholesky LDL^T Sum-of-Squares (SOS) and Lukacs certificates,
- p-adic valuation tools across all 21 primes in 101..199,
- IdentityQuantizerBundle for unbranched symbolic execution of Battle 3 ASTs.
"""

from fractions import Fraction

PRIMES_101_199 = (
    101,
    103,
    107,
    109,
    113,
    127,
    131,
    137,
    139,
    149,
    151,
    157,
    163,
    167,
    173,
    179,
    181,
    191,
    193,
    197,
    199,
)

# Tag convention via .numerator and .denominator attributes:
# - int / Fraction: x.numerator is int, x.denominator >= 1 (int)
# - QQPoly:   x.numerator is self, x.denominator == -1
# - QQRat:    x.numerator is self, x.denominator == -2
# - QQJet:    x.numerator is self, x.denominator == -3
# - QQMVPoly: x.numerator is self, x.denominator == -4
# - QQMVRat:  x.numerator is self, x.denominator == -5


def q_abs(x):
    return x if x >= 0 else -x


def q_max(a, b):
    return a if a >= b else b


def q_min(a, b):
    return a if a <= b else b


def to_qq(x):
    return Fraction(x.numerator, x.denominator)


def igcd(a, b):
    ap = a if a >= 0 else -a
    bp = b if b >= 0 else -b
    if bp == 0:
        return ap
    return bp // Fraction(ap, bp).denominator


def ilcm(a, b):
    ap = a if a >= 0 else -a
    bp = b if b >= 0 else -b
    if ap == 0 or bp == 0:
        return 0
    return (ap // igcd(ap, bp)) * bp


def dyadic_valuation(n):
    np = n if n >= 0 else -n
    if np == 0:
        return 0
    low = np & (-np)
    return low.bit_length() - 1


def odd_part(n):
    np = n if n >= 0 else -n
    if np == 0:
        return 0
    low = np & (-np)
    return np >> (low.bit_length() - 1)


def int_padic_valuation(n, p):
    np = n if n >= 0 else -n
    if np == 0 or p <= 1:
        return 0
    cur = np
    v = 0
    bound = np.bit_length() + 1
    for _ in range(bound):
        q, r = divmod(cur, p)
        if r != 0:
            break
        cur = q
        v += 1
    return v


def padic_valuation(q, p):
    num = q.numerator
    den = q.denominator
    if num == 0:
        return 0
    return int_padic_valuation(num, p) - int_padic_valuation(den, p)


def _pow_ring(base, n, one_val):
    if n < 0:
        raise ValueError("Negative exponent forbidden in ring power")
    if n == 0:
        return one_val
    res = one_val
    cur = base
    rem = n
    for _ in range(64):
        if rem == 0:
            break
        if rem & 1:
            res = res * cur
        rem >>= 1
        if rem > 0:
            cur = cur * cur
    return res


def _strip_qq_coeffs(coeffs):
    ln = len(coeffs)
    if ln == 0:
        return (Fraction(0, 1),)
    d = ln - 1
    for _ in range(ln):
        if d > 0 and coeffs[d] == 0:
            d -= 1
        else:
            break
    return tuple(to_qq(coeffs[i]) for i in range(d + 1))


def _poly_to_prim_int(poly_coeffs):
    d = len(poly_coeffs) - 1
    if d == 0 and poly_coeffs[0] == 0:
        return [0]
    lcm_den = 1
    for c in poly_coeffs:
        lcm_den = ilcm(lcm_den, c.denominator)
    ints = [(c.numerator * (lcm_den // c.denominator)) for c in poly_coeffs]
    g = 0
    for val in ints:
        if val != 0:
            g = igcd(g, val)
            if g == 1:
                break
    if g > 1:
        ints = [val // g for val in ints]
    return ints


def _strip_int_coeffs(ints):
    ln = len(ints)
    if ln == 0:
        return [0]
    d = ln - 1
    for _ in range(ln):
        if d > 0 and ints[d] == 0:
            d -= 1
        else:
            break
    return ints[: d + 1]


def _prim_part_pos(ints):
    s = _strip_int_coeffs(ints)
    if len(s) == 1 and s[0] == 0:
        return [0]
    g = 0
    for val in s:
        if val != 0:
            g = igcd(g, val)
            if g == 1:
                break
    if g > 1:
        return [val // g for val in s]
    return s


def _subresultant_step_2adic(a, b, denom):
    n = len(b) - 1
    bn = b[n]
    bn2 = bn * bn
    q1 = a[n + 1] * bn
    q0 = a[n] * bn - a[n + 1] * (b[n - 1] if n >= 1 else 0)
    if denom == 1:
        out = [0] * n
        out[0] = bn2 * a[0] - q0 * b[0]
        for i in range(1, n):
            out[i] = bn2 * a[i] - q1 * b[i - 1] - q0 * b[i]
        return out
    a_bits = 0
    for x in a:
        bl = x.bit_length()
        if bl > a_bits:
            a_bits = bl
    b_bits = 0
    for x in b:
        bl = x.bit_length()
        if bl > b_bits:
            b_bits = bl
    q1_pos = q1 if q1 >= 0 else -q1
    q0_pos = q0 if q0 >= 0 else -q0
    v_bound = (bn2 << a_bits) + ((q1_pos + q0_pos) << b_bits)
    q_bits = (v_bound // denom).bit_length() + 3
    tz = (denom & (-denom)).bit_length() - 1
    d_odd = denom >> tz
    tot_bits = q_bits + tz
    mask_tot = (1 << tot_bits) - 1
    mask_q = (1 << q_bits) - 1
    half_q = 1 << (q_bits - 1)
    full_q = 1 << q_bits

    inv = 1
    cur_prec = 1
    for _ in range(28):
        if cur_prec >= tot_bits:
            break
        cur_prec <<= 1
        ms = (1 << cur_prec) - 1
        inv = (inv * (2 - ((d_odd * inv) & ms))) & ms
    inv &= mask_tot

    ca = (bn2 * inv) & mask_tot
    cb1 = ((-q1) * inv) & mask_tot
    cb0 = ((-q0) * inv) & mask_tot

    out = [0] * n
    v0 = ((ca * (a[0] & mask_tot) + cb0 * (b[0] & mask_tot)) & mask_tot) >> tz
    u0 = v0 & mask_q
    out[0] = u0 - full_q if u0 >= half_q else u0
    for i in range(1, n):
        vi = (
            (
                ca * (a[i] & mask_tot)
                + cb1 * (b[i - 1] & mask_tot)
                + cb0 * (b[i] & mask_tot)
            )
            & mask_tot
        ) >> tz
        ui = vi & mask_q
        out[i] = ui - full_q if ui >= half_q else ui
    return out


def _signed_prem_general(a, b, denom):
    m = len(a) - 1
    n = len(b) - 1
    bn = b[n]
    bn_pos = bn if bn >= 0 else -bn
    bn_sgn = 1 if bn >= 0 else -1
    r = list(a)
    for k in range(m, n - 1, -1):
        rk = r[k]
        if rk != 0:
            scale = rk * bn_sgn
            off = k - n
            for j in range(off):
                r[j] *= bn_pos
            for j in range(n):
                r[off + j] = r[off + j] * bn_pos - scale * b[j]
        else:
            for j in range(k):
                r[j] *= bn_pos
        r[k] = 0
    if denom == 1:
        return r[:n]
    return [r[i] // denom for i in range(n)]


def _gcd_int_polys(pa, pb):
    a = _prim_part_pos(pa)
    b = _prim_part_pos(pb)
    if len(a) == 1 and a[0] == 0:
        if b[-1] < 0:
            return [-x for x in b]
        return b
    if len(b) == 1 and b[0] == 0:
        if a[-1] < 0:
            return [-x for x in a]
        return a
    if len(a) < len(b):
        a, b = b, a
    if len(b) == 1:
        return [1]
    g = 1
    h = 1
    step = 0
    max_iters = len(a) + len(b) + 2
    for _ in range(max_iters):
        m = len(a) - 1
        n = len(b) - 1
        if n == 0:
            return [1]
        step += 1
        delta = m - n
        denom = 1 if step == 1 else g * (h ** delta)
        if delta == 1:
            rem = _subresultant_step_2adic(a, b, denom)
        else:
            rem = _signed_prem_general(a, b, denom)
        rem = _strip_int_coeffs(rem)
        if len(rem) == 1 and rem[0] == 0:
            res = _prim_part_pos(b)
            if res[-1] < 0:
                res = [-x for x in res]
            return res
        bn = b[n]
        bn_pos = bn if bn >= 0 else -bn
        g = bn_pos
        h = (g ** delta) if step == 1 else (g ** delta) // (h ** (delta - 1))
        a = b
        b = rem
    res = _prim_part_pos(b)
    if res[-1] < 0:
        res = [-x for x in res]
    return res


class QQPoly:
    """Exact univariate polynomial over Q with Fraction coefficients."""

    def __init__(self, coeffs, var_name="x"):
        self.coeffs = _strip_qq_coeffs(coeffs)
        self.var_name = var_name
        self.numerator = self
        self.denominator = -1

    @classmethod
    def zero(cls, var_name="x"):
        return cls((Fraction(0, 1),), var_name)

    @classmethod
    def one(cls, var_name="x"):
        return cls((Fraction(1, 1),), var_name)

    @classmethod
    def const(cls, c, var_name="x"):
        return cls((to_qq(c),), var_name)

    @classmethod
    def var(cls, var_name="x"):
        return cls((Fraction(0, 1), Fraction(1, 1)), var_name)

    @classmethod
    def vars(cls, names):
        return QQMVPoly.vars(names)

    def degree(self):
        if len(self.coeffs) == 1 and self.coeffs[0] == 0:
            return 0
        return len(self.coeffs) - 1

    def is_zero(self):
        return len(self.coeffs) == 1 and self.coeffs[0] == 0

    def lc(self):
        return self.coeffs[-1]

    def monic(self):
        if self.is_zero():
            return self
        lead = self.coeffs[-1]
        if lead == 1:
            return self
        inv_l = Fraction(1, 1) / lead
        return QQPoly(tuple(c * inv_l for c in self.coeffs), self.var_name)

    def __pos__(self):
        return self

    def __neg__(self):
        return QQPoly(tuple(-c for c in self.coeffs), self.var_name)

    def __eq__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            return len(self.coeffs) == 1 and self.coeffs[0] == q
        if tag == -1:
            return self.coeffs == val.coeffs
        if tag in (-2, -3, -4, -5):
            return val == self
        return False

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self
            c0 = self.coeffs[0] + q
            return QQPoly((c0,) + self.coeffs[1:], self.var_name)
        if tag == -1:
            la = len(self.coeffs)
            lb = len(val.coeffs)
            ln = la if la >= lb else lb
            out = []
            for i in range(ln):
                ca = self.coeffs[i] if i < la else Fraction(0, 1)
                cb = val.coeffs[i] if i < lb else Fraction(0, 1)
                out.append(ca + cb)
            return QQPoly(out, self.var_name)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self
            c0 = self.coeffs[0] - q
            return QQPoly((c0,) + self.coeffs[1:], self.var_name)
        if tag == -1:
            la = len(self.coeffs)
            lb = len(val.coeffs)
            ln = la if la >= lb else lb
            out = []
            for i in range(ln):
                ca = self.coeffs[i] if i < la else Fraction(0, 1)
                cb = val.coeffs[i] if i < lb else Fraction(0, 1)
                out.append(ca - cb)
            return QQPoly(out, self.var_name)
        return NotImplemented

    def __rsub__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return (-self) + Fraction(val, tag)
        return NotImplemented

    def __mul__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0 or self.is_zero():
                return QQPoly.zero(self.var_name)
            if q == 1:
                return self
            return QQPoly(tuple(c * q for c in self.coeffs), self.var_name)
        if tag == -1:
            if self.is_zero() or val.is_zero():
                return QQPoly.zero(self.var_name)
            la = len(self.coeffs)
            lb = len(val.coeffs)
            out = [Fraction(0, 1) for _ in range(la + lb - 1)]
            for i in range(la):
                ci = self.coeffs[i]
                if ci != 0:
                    for j in range(lb):
                        out[i + j] += ci * val.coeffs[j]
            return QQPoly(out, self.var_name)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                raise ZeroDivisionError("QQPoly division by zero scalar")
            inv_q = Fraction(1, 1) / q
            return self * inv_q
        if tag == -1:
            return QQRat(self, val)
        if tag == -2:
            return QQRat(self, QQPoly.one(self.var_name)) / val
        return NotImplemented

    def __rtruediv__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQRat(QQPoly.const(Fraction(val, tag), self.var_name), self)
        return NotImplemented

    def __floordiv__(self, other):
        q, _ = self.divmod_poly(other)
        return q

    def __mod__(self, other):
        _, r = self.divmod_poly(other)
        return r

    def __pow__(self, n):
        return _pow_ring(self, n, QQPoly.one(self.var_name))

    def divmod_poly(self, divisor):
        val = divisor.numerator
        tag = divisor.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            return (self / q, QQPoly.zero(self.var_name))
        if val.is_zero():
            raise ZeroDivisionError("Polynomial division by zero")
        if self.is_zero():
            return (QQPoly.zero(self.var_name), QQPoly.zero(self.var_name))
        da = self.degree()
        db = val.degree()
        if da < db:
            return (QQPoly.zero(self.var_name), self)
        rem = list(self.coeffs)
        quot = [Fraction(0, 1) for _ in range(da - db + 1)]
        inv_lb = Fraction(1, 1) / val.coeffs[db]
        for k in range(da, db - 1, -1):
            ck = rem[k]
            if ck != 0:
                factor = ck * inv_lb
                quot[k - db] = factor
                off = k - db
                for j in range(db + 1):
                    rem[off + j] -= factor * val.coeffs[j]
        return (
            QQPoly(quot, self.var_name),
            QQPoly(rem[:db] if db > 0 else [Fraction(0, 1)], self.var_name),
        )

    def diff(self, m=1):
        cur = self
        for _ in range(m):
            if len(cur.coeffs) <= 1:
                return QQPoly.zero(self.var_name)
            cur = QQPoly(
                tuple(Fraction(i, 1) * cur.coeffs[i] for i in range(1, len(cur.coeffs))),
                self.var_name,
            )
        return cur

    def eval_at(self, x_val):
        d = len(self.coeffs) - 1
        acc = self.coeffs[d]
        for k in range(d - 1, -1, -1):
            acc = acc * x_val + self.coeffs[k]
        return acc

    def compose(self, inner):
        return self.eval_at(inner)

    def shift(self, a):
        x = QQPoly.var(self.var_name)
        return self.eval_at(x + a)

    def reverse(self, deg=None):
        d = self.degree() if deg is None else deg
        out = [Fraction(0, 1) for _ in range(d + 1)]
        for i in range(len(self.coeffs)):
            if d - i >= 0:
                out[d - i] = self.coeffs[i]
        return QQPoly(out, self.var_name)

    def __repr__(self):
        return f"QQPoly({self.coeffs!r}, var_name={self.var_name!r})"


def gcd_poly(a, b):
    """Exact monic GCD of two univariate QQPoly polynomials via subresultant PRS."""
    if a.is_zero() and b.is_zero():
        return QQPoly.zero(a.var_name)
    if a.is_zero():
        return b.monic()
    if b.is_zero():
        return a.monic()
    ia = _poly_to_prim_int(a.coeffs)
    ib = _poly_to_prim_int(b.coeffs)
    ig = _gcd_int_polys(ia, ib)
    return QQPoly(tuple(Fraction(x, 1) for x in ig), a.var_name).monic()


class QQRat:
    """Exact reduced univariate rational function num(x) / den(x) over Q."""

    def __init__(self, num, den=1, var_name="x", _reduced=False):
        self.numerator = self
        self.denominator = -2

        vn = num.numerator
        tn = num.denominator
        if tn >= 1:
            p_num = QQPoly.const(Fraction(vn, tn), var_name)
        elif tn == -1:
            p_num = vn
            var_name = vn.var_name
        else:
            raise ValueError("QQRat numerator must be scalar or QQPoly")

        vd = den.numerator
        td = den.denominator
        if td >= 1:
            p_den = QQPoly.const(Fraction(vd, td), var_name)
        elif td == -1:
            p_den = vd
        else:
            raise ValueError("QQRat denominator must be scalar or QQPoly")

        if p_den.is_zero():
            raise ZeroDivisionError("QQRat zero denominator")

        if p_num.is_zero():
            self.num = QQPoly.zero(var_name)
            self.den = QQPoly.one(var_name)
            self.var_name = var_name
            return

        if _reduced:
            lead = p_den.lc()
            if lead != 1:
                inv_l = Fraction(1, 1) / lead
                p_num = p_num * inv_l
                p_den = p_den * inv_l
            self.num = p_num
            self.den = p_den
            self.var_name = var_name
            return

        g = gcd_poly(p_num, p_den)
        if g.degree() > 0:
            p_num = p_num // g
            p_den = p_den // g
        lead = p_den.lc()
        if lead != 1:
            inv_l = Fraction(1, 1) / lead
            p_num = p_num * inv_l
            p_den = p_den * inv_l
        self.num = p_num
        self.den = p_den
        self.var_name = var_name

    @classmethod
    def var(cls, var_name="x"):
        return cls(QQPoly.var(var_name), QQPoly.one(var_name), var_name, _reduced=True)

    @classmethod
    def vars(cls, names):
        return QQMVRat.vars(names)

    @classmethod
    def const(cls, c, var_name="x"):
        return cls(QQPoly.const(c, var_name), QQPoly.one(var_name), var_name, _reduced=True)

    def is_zero(self):
        return self.num.is_zero()

    def __pos__(self):
        return self

    def __neg__(self):
        return QQRat(-self.num, self.den, self.var_name, _reduced=True)

    def _coerce_other(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQRat.const(Fraction(val, tag), self.var_name)
        if tag == -1:
            return QQRat(val, QQPoly.one(val.var_name), val.var_name, _reduced=True)
        if tag == -2:
            return val
        return None

    def __eq__(self, other):
        co = self._coerce_other(other)
        if co is not None:
            return self.num == co.num and self.den == co.den
        val = other.numerator
        tag = other.denominator
        if tag in (-3, -4, -5):
            return val == self
        return False

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero():
            return co
        if co.is_zero():
            return self
        d_self = self.den.degree()
        d_co = co.den.degree()
        if d_self == 0 and d_co == 0:
            return QQRat(self.num + co.num, QQPoly.one(self.var_name), self.var_name, _reduced=True)
        if d_co == 0:
            return QQRat(self.num + co.num * self.den, self.den, self.var_name, _reduced=True)
        if d_self == 0:
            return QQRat(self.num * co.den + co.num, co.den, self.var_name, _reduced=True)
        if self.den == co.den:
            return QQRat(self.num + co.num, self.den, self.var_name)
        g = gcd_poly(self.den, co.den)
        if g.degree() == 0:
            new_num = self.num * co.den + co.num * self.den
            new_den = self.den * co.den
            return QQRat(new_num, new_den, self.var_name, _reduced=True)
        u1 = self.den // g
        u2 = co.den // g
        t = self.num * u2 + co.num * u1
        if t.is_zero():
            return QQRat.const(0, self.var_name)
        g2 = gcd_poly(t, g)
        if g2.degree() > 0:
            t = t // g2
            g = g // g2
        return QQRat(t, u1 * (u2 * g), self.var_name, _reduced=True)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self + (-co)

    def __rsub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co + (-self)

    def __mul__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero() or co.is_zero():
            return QQRat.const(0, self.var_name)
        if co.num.degree() == 0 and co.den.degree() == 0:
            return QQRat(self.num * co.num.coeffs[0], self.den, self.var_name, _reduced=True)
        if self.num.degree() == 0 and self.den.degree() == 0:
            return QQRat(co.num * self.num.coeffs[0], co.den, self.var_name, _reduced=True)
        if self.den == co.den:
            return QQRat(self.num * co.num, self.den * co.den, self.var_name, _reduced=True)
        g1 = gcd_poly(self.num, co.den) if co.den.degree() > 0 else QQPoly.one(self.var_name)
        g2 = gcd_poly(co.num, self.den) if self.den.degree() > 0 else QQPoly.one(self.var_name)
        n1 = self.num // g1 if g1.degree() > 0 else self.num
        d2 = co.den // g1 if g1.degree() > 0 else co.den
        n2 = co.num // g2 if g2.degree() > 0 else co.num
        d1 = self.den // g2 if g2.degree() > 0 else self.den
        return QQRat(n1 * n2, d1 * d2, self.var_name, _reduced=True)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if co.is_zero():
            raise ZeroDivisionError("QQRat division by zero")
        if co.num.degree() == 0 and co.den.degree() == 0:
            return QQRat(self.num / co.num.coeffs[0], self.den, self.var_name, _reduced=True)
        if self.den == co.den:
            return QQRat(self.num, co.num, self.var_name)
        inv_co = QQRat(co.den, co.num, self.var_name, _reduced=True)
        return self * inv_co

    def __rtruediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co / self

    def __pow__(self, n):
        if n < 0:
            raise ValueError("Negative exponent forbidden in QQRat.__pow__")
        if n == 0:
            return QQRat.const(1, self.var_name)
        return QQRat(self.num ** n, self.den ** n, self.var_name, _reduced=True)

    def diff(self):
        if self.num.is_zero():
            return QQRat.const(0, self.var_name)
        np = self.num.diff()
        if self.den.degree() == 0:
            return QQRat(np, QQPoly.one(self.var_name), self.var_name, _reduced=True)
        dp = self.den.diff()
        g = gcd_poly(self.den, dp)
        d_sqfree = self.den // g if g.degree() > 0 else self.den
        dp_red = dp // g if g.degree() > 0 else dp
        top = np * d_sqfree - self.num * dp_red
        if top.is_zero():
            return QQRat.const(0, self.var_name)
        g2 = gcd_poly(top, d_sqfree)
        if g2.degree() > 0:
            top = top // g2
            den_out = self.den * (d_sqfree // g2)
        else:
            den_out = self.den * d_sqfree
        return QQRat(top, den_out, self.var_name, _reduced=True)

    def eval_at(self, x_val):
        return self.num.eval_at(x_val) / self.den.eval_at(x_val)

    def compose(self, inner):
        return self.eval_at(inner)

    def __repr__(self):
        return f"QQRat({self.num!r}, {self.den!r})"


def _pad_exp(e, nvars):
    if len(e) == nvars:
        return e
    return e + (0,) * (nvars - len(e))


class QQMVPoly:
    """Sparse exact multivariate polynomial over Q."""

    def __init__(self, nvars, terms, var_names=None):
        self.numerator = self
        self.denominator = -4
        self.nvars = nvars
        clean = {}
        for k, v in terms.items():
            qv = to_qq(v)
            if qv != 0:
                ek = _pad_exp(k, nvars)
                clean[ek] = clean.get(ek, Fraction(0, 1)) + qv
                if clean[ek] == 0:
                    del clean[ek]
        self.terms = clean
        if var_names is None:
            self.var_names = tuple(f"x{i}" for i in range(nvars))
        else:
            self.var_names = tuple(var_names)

    @classmethod
    def zero(cls, nvars, var_names=None):
        return cls(nvars, {}, var_names)

    @classmethod
    def one(cls, nvars, var_names=None):
        return cls(nvars, {(0,) * nvars: Fraction(1, 1)}, var_names)

    @classmethod
    def const(cls, c, nvars, var_names=None):
        qc = to_qq(c)
        if qc == 0:
            return cls(nvars, {}, var_names)
        return cls(nvars, {(0,) * nvars: qc}, var_names)

    @classmethod
    def vars(cls, names_seq):
        names = tuple(names_seq)
        n = len(names)
        out = []
        for i in range(n):
            e = tuple(1 if j == i else 0 for j in range(n))
            out.append(cls(n, {e: Fraction(1, 1)}, names))
        return tuple(out)

    @classmethod
    def vars_n(cls, n, prefix="x"):
        names = tuple(f"{prefix}{i}" for i in range(n))
        return cls.vars(names)

    def is_zero(self):
        return len(self.terms) == 0

    def is_const(self):
        if len(self.terms) == 0:
            return True
        if len(self.terms) > 1:
            return False
        z = (0,) * self.nvars
        return z in self.terms

    def const_val(self):
        return self.terms.get((0,) * self.nvars, Fraction(0, 1))

    def total_degree(self):
        best = 0
        for e in self.terms:
            s = sum(e)
            if s > best:
                best = s
        return best

    def _merge_meta(self, other):
        nv = self.nvars if self.nvars >= other.nvars else other.nvars
        vn = self.var_names if len(self.var_names) >= nv else other.var_names
        return nv, vn

    def __pos__(self):
        return self

    def __neg__(self):
        return QQMVPoly(self.nvars, {k: -v for k, v in self.terms.items()}, self.var_names)

    def __eq__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            return self.is_const() and self.const_val() == q
        if tag == -4:
            nv, _ = self._merge_meta(val)
            t1 = {_pad_exp(k, nv): v for k, v in self.terms.items()}
            t2 = {_pad_exp(k, nv): v for k, v in val.terms.items()}
            return t1 == t2
        if tag in (-3, -5):
            return val == self
        return False

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self
            out = dict(self.terms)
            z = (0,) * self.nvars
            out[z] = out.get(z, Fraction(0, 1)) + q
            return QQMVPoly(self.nvars, out, self.var_names)
        if tag == -4:
            nv, vn = self._merge_meta(val)
            out = {_pad_exp(k, nv): v for k, v in self.terms.items()}
            for k, v in val.terms.items():
                ek = _pad_exp(k, nv)
                out[ek] = out.get(ek, Fraction(0, 1)) + v
            return QQMVPoly(nv, out, vn)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return self + (-Fraction(val, tag))
        if tag == -4:
            return self + (-val)
        return NotImplemented

    def __rsub__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return (-self) + Fraction(val, tag)
        return NotImplemented

    def __mul__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0 or self.is_zero():
                return QQMVPoly.zero(self.nvars, self.var_names)
            if q == 1:
                return self
            return QQMVPoly(self.nvars, {k: v * q for k, v in self.terms.items()}, self.var_names)
        if tag == -4:
            nv, vn = self._merge_meta(val)
            if self.is_zero() or val.is_zero():
                return QQMVPoly.zero(nv, vn)
            out = {}
            t1 = [(_pad_exp(k, nv), v) for k, v in self.terms.items()]
            t2 = [(_pad_exp(k, nv), v) for k, v in val.terms.items()]
            for e1, c1 in t1:
                for e2, c2 in t2:
                    ek = tuple(e1[i] + e2[i] for i in range(nv))
                    out[ek] = out.get(ek, Fraction(0, 1)) + c1 * c2
            return QQMVPoly(nv, out, vn)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                raise ZeroDivisionError("QQMVPoly division by zero scalar")
            return self * (Fraction(1, 1) / q)
        if tag == -4:
            if val.is_const():
                c = val.const_val()
                if c == 0:
                    raise ZeroDivisionError("QQMVPoly division by zero")
                return self * (Fraction(1, 1) / c)
            ex = self.exact_div_or_none(val)
            if ex is not None:
                return ex
            return QQMVRat(self, val)
        if tag == -5:
            return QQMVRat(self, QQMVPoly.one(self.nvars, self.var_names)) / val
        return NotImplemented

    def __rtruediv__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQMVRat(QQMVPoly.const(Fraction(val, tag), self.nvars, self.var_names), self)
        return NotImplemented

    def __pow__(self, n):
        return _pow_ring(self, n, QQMVPoly.one(self.nvars, self.var_names))

    def _lex_lead(self):
        best_e = None
        for e in self.terms:
            if best_e is None or e > best_e:
                best_e = e
        return best_e, self.terms[best_e]

    def exact_div_or_none(self, divisor):
        if divisor.is_zero():
            raise ZeroDivisionError("QQMVPoly exact division by zero")
        if self.is_zero():
            return QQMVPoly.zero(self.nvars, self.var_names)
        if self == divisor:
            return QQMVPoly.one(self.nvars, self.var_names)
        if divisor.is_const():
            return self * (Fraction(1, 1) / divisor.const_val())
        nv, vn = self._merge_meta(divisor)
        rem_terms = {_pad_exp(k, nv): v for k, v in self.terms.items()}
        div_terms = {_pad_exp(k, nv): v for k, v in divisor.terms.items()}
        lead_d = None
        for e in div_terms:
            if lead_d is None or e > lead_d:
                lead_d = e
        lc_d = div_terms[lead_d]
        inv_lc_d = Fraction(1, 1) / lc_d
        quot_terms = {}
        bound = len(rem_terms) * 8 + 32
        for _ in range(bound):
            if len(rem_terms) == 0:
                return QQMVPoly(nv, quot_terms, vn)
            lead_r = None
            for e in rem_terms:
                if lead_r is None or e > lead_r:
                    lead_r = e
            for i in range(nv):
                if lead_r[i] < lead_d[i]:
                    return None
            diff_e = tuple(lead_r[i] - lead_d[i] for i in range(nv))
            scale = rem_terms[lead_r] * inv_lc_d
            quot_terms[diff_e] = quot_terms.get(diff_e, Fraction(0, 1)) + scale
            for ed, cd in div_terms.items():
                target_e = tuple(diff_e[i] + ed[i] for i in range(nv))
                nv_val = rem_terms.get(target_e, Fraction(0, 1)) - scale * cd
                if nv_val == 0:
                    if target_e in rem_terms:
                        del rem_terms[target_e]
                else:
                    rem_terms[target_e] = nv_val
        return None

    def diff(self, var_idx):
        out = {}
        for e, c in self.terms.items():
            deg = e[var_idx]
            if deg > 0:
                ne = tuple(e[i] - 1 if i == var_idx else e[i] for i in range(self.nvars))
                out[ne] = c * Fraction(deg, 1)
        return QQMVPoly(self.nvars, out, self.var_names)

    def eval_at(self, point):
        if self.is_zero():
            return Fraction(0, 1)
        one = point[0] * 0 + 1 if len(point) > 0 else Fraction(1, 1)
        acc = point[0] * 0 if len(point) > 0 else Fraction(0, 1)
        for e, c in self.terms.items():
            term = one * c
            for i in range(self.nvars):
                if e[i] > 0:
                    term = term * (point[i] ** e[i])
            acc = acc + term
        return acc

    def __repr__(self):
        return f"QQMVPoly(nvars={self.nvars}, terms={self.terms!r})"


class QQMVRat:
    """Multivariate rational function num / den over QQMVPoly with exact factor cancellation."""

    def __init__(self, num, den=None):
        self.numerator = self
        self.denominator = -5

        vn = num.numerator
        tn = num.denominator
        if tn >= 1:
            raise ValueError("Use QQMVRat.const for scalar construction")
        if tn != -4:
            raise ValueError("QQMVRat numerator must be QQMVPoly")
        p_num = vn

        if den is None:
            p_den = QQMVPoly.one(p_num.nvars, p_num.var_names)
        else:
            vd = den.numerator
            td = den.denominator
            if td >= 1:
                p_den = QQMVPoly.const(Fraction(vd, td), p_num.nvars, p_num.var_names)
            elif td == -4:
                p_den = vd
            else:
                raise ValueError("QQMVRat denominator must be scalar or QQMVPoly")

        if p_den.is_zero():
            raise ZeroDivisionError("QQMVRat zero denominator")

        nv, vnames = p_num._merge_meta(p_den)
        if p_num.is_zero():
            self.num = QQMVPoly.zero(nv, vnames)
            self.den = QQMVPoly.one(nv, vnames)
            return

        if p_num == p_den:
            self.num = QQMVPoly.one(nv, vnames)
            self.den = QQMVPoly.one(nv, vnames)
            return

        if (-p_num) == p_den:
            self.num = -QQMVPoly.one(nv, vnames)
            self.den = QQMVPoly.one(nv, vnames)
            return

        q_direct = p_num.exact_div_or_none(p_den)
        if q_direct is not None:
            self.num = q_direct
            self.den = QQMVPoly.one(nv, vnames)
            return

        q_inv = p_den.exact_div_or_none(p_num)
        if q_inv is not None:
            _, lc_i = q_inv._lex_lead()
            inv_lc = Fraction(1, 1) / lc_i
            self.num = QQMVPoly.const(inv_lc, nv, vnames)
            self.den = q_inv * inv_lc
            return

        _, lc_d = p_den._lex_lead()
        if lc_d != 1:
            inv_l = Fraction(1, 1) / lc_d
            p_num = p_num * inv_l
            p_den = p_den * inv_l

        self.num = p_num
        self.den = p_den

    @classmethod
    def vars(cls, names):
        polys = QQMVPoly.vars(names)
        return tuple(cls(p) for p in polys)

    @classmethod
    def vars_n(cls, n, prefix="x"):
        polys = QQMVPoly.vars_n(n, prefix)
        return tuple(cls(p) for p in polys)

    @classmethod
    def const(cls, c, nvars, var_names=None):
        return cls(QQMVPoly.const(c, nvars, var_names))

    def is_zero(self):
        return self.num.is_zero()

    def _coerce_other(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQMVRat.const(Fraction(val, tag), self.num.nvars, self.num.var_names)
        if tag == -4:
            return QQMVRat(val)
        if tag == -5:
            return val
        return None

    def __pos__(self):
        return self

    def __neg__(self):
        out = QQMVRat.__new__(QQMVRat)
        out.numerator = out
        out.denominator = -5
        out.num = -self.num
        out.den = self.den
        return out

    def __eq__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self.is_zero()
            return (self.num - self.den * q).is_zero()
        co = self._coerce_other(other)
        if co is not None:
            if self.den == co.den:
                return self.num == co.num
            return (self.num * co.den - co.num * self.den).is_zero()
        if tag == -3:
            return val == self
        return False

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero():
            return co
        if co.is_zero():
            return self
        if self.den == co.den:
            return QQMVRat(self.num + co.num, self.den)
        q21 = co.den.exact_div_or_none(self.den)
        if q21 is not None:
            return QQMVRat(self.num * q21 + co.num, co.den)
        q12 = self.den.exact_div_or_none(co.den)
        if q12 is not None:
            return QQMVRat(self.num + co.num * q12, self.den)
        return QQMVRat(self.num * co.den + co.num * self.den, self.den * co.den)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self + (-co)

    def __rsub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co + (-self)

    def __mul__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero() or co.is_zero():
            return QQMVRat.const(0, self.num.nvars, self.num.var_names)
        n1, d1 = self.num, self.den
        n2, d2 = co.num, co.den
        if n1 == d2:
            n1 = QQMVPoly.one(n1.nvars, n1.var_names)
            d2 = n1
        else:
            q = n1.exact_div_or_none(d2)
            if q is not None:
                n1 = q
                d2 = QQMVPoly.one(n1.nvars, n1.var_names)
            else:
                qi = d2.exact_div_or_none(n1)
                if qi is not None:
                    n1 = QQMVPoly.one(n1.nvars, n1.var_names)
                    d2 = qi
        if n2 == d1:
            n2 = QQMVPoly.one(n2.nvars, n2.var_names)
            d1 = n2
        else:
            q = n2.exact_div_or_none(d1)
            if q is not None:
                n2 = q
                d1 = QQMVPoly.one(n2.nvars, n2.var_names)
            else:
                qi = d1.exact_div_or_none(n2)
                if qi is not None:
                    n2 = QQMVPoly.one(n2.nvars, n2.var_names)
                    d1 = qi
        return QQMVRat(n1 * n2, d1 * d2)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if co.is_zero():
            raise ZeroDivisionError("QQMVRat division by zero")
        inv_co = QQMVRat(co.den, co.num)
        return self * inv_co

    def __rtruediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co / self

    def __pow__(self, n):
        return _pow_ring(self, n, QQMVRat.const(1, self.num.nvars, self.num.var_names))

    def eval_at(self, point):
        return self.num.eval_at(point) / self.den.eval_at(point)

    def __repr__(self):
        return f"QQMVRat({self.num!r}, {self.den!r})"


class QQJet:
    """Truncated Taylor series c_0 + c_1*eps + ... + c_k*eps^k over Q or symbolic rings."""

    def __init__(self, coeffs, order=None):
        self.numerator = self
        self.denominator = -3
        seq = list(coeffs)
        if order is None:
            order = len(seq) - 1 if len(seq) > 0 else 0
        self.order = order
        zero_elem = seq[0] * 0 if len(seq) > 0 else Fraction(0, 1)
        padded = []
        for i in range(order + 1):
            if i < len(seq):
                val = seq[i].numerator
                tag = seq[i].denominator
                padded.append(Fraction(val, tag) if tag >= 1 else val)
            else:
                padded.append(zero_elem)
        self.coeffs = tuple(padded)

    @classmethod
    def const(cls, c, order):
        val = c.numerator
        tag = c.denominator
        c0 = Fraction(val, tag) if tag >= 1 else val
        z = c0 * 0
        return cls((c0,) + (z,) * order, order)

    @classmethod
    def var(cls, c0=0, order=4):
        val = c0.numerator
        tag = c0.denominator
        v0 = Fraction(val, tag) if tag >= 1 else val
        z = v0 * 0
        one = z + 1
        if order == 0:
            return cls((v0,), 0)
        return cls((v0, one) + (z,) * (order - 1), order)

    def coeff(self, m):
        if m < 0 or m > self.order:
            return self.coeffs[0] * 0
        return self.coeffs[m]

    def is_zero(self):
        for c in self.coeffs:
            if c != 0:
                return False
        return True

    def __pos__(self):
        return self

    def __neg__(self):
        return QQJet(tuple(-c for c in self.coeffs), self.order)

    def __eq__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            for i in range(ord_m + 1):
                if self.coeffs[i] != val.coeffs[i]:
                    return False
            return True
        c0 = Fraction(val, tag) if tag >= 1 else val
        if self.coeffs[0] != c0:
            return False
        for i in range(1, self.order + 1):
            if self.coeffs[i] != 0:
                return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            return QQJet(
                tuple(self.coeffs[i] + val.coeffs[i] for i in range(ord_m + 1)),
                ord_m,
            )
        sc = Fraction(val, tag) if tag >= 1 else val
        return QQJet((self.coeffs[0] + sc,) + self.coeffs[1:], self.order)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            return QQJet(
                tuple(self.coeffs[i] - val.coeffs[i] for i in range(ord_m + 1)),
                ord_m,
            )
        sc = Fraction(val, tag) if tag >= 1 else val
        return QQJet((self.coeffs[0] - sc,) + self.coeffs[1:], self.order)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            out = []
            for m in range(ord_m + 1):
                acc = self.coeffs[0] * val.coeffs[m]
                for j in range(1, m + 1):
                    acc = acc + self.coeffs[j] * val.coeffs[m - j]
                out.append(acc)
            return QQJet(out, ord_m)
        sc = Fraction(val, tag) if tag >= 1 else val
        return QQJet(tuple(c * sc for c in self.coeffs), self.order)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            v_shift = 0
            for k in range(ord_m + 1):
                if self.coeffs[k] == 0 and val.coeffs[k] == 0:
                    v_shift += 1
                else:
                    break
            if v_shift > ord_m:
                raise ZeroDivisionError("QQJet 0/0 indeterminate up to full jet order")
            if v_shift > 0:
                new_ord = ord_m - v_shift
                a_sub = QQJet(self.coeffs[v_shift : ord_m + 1], new_ord)
                b_sub = QQJet(val.coeffs[v_shift : ord_m + 1], new_ord)
                return a_sub / b_sub
            b0 = val.coeffs[0]
            if b0 == 0:
                raise ZeroDivisionError("QQJet division by jet with zero constant term")
            out = [self.coeffs[0] / b0]
            for m in range(1, ord_m + 1):
                acc = self.coeffs[m]
                for j in range(1, m + 1):
                    acc = acc - val.coeffs[j] * out[m - j]
                out.append(acc / b0)
            return QQJet(out, ord_m)
        sc = Fraction(val, tag) if tag >= 1 else val
        if sc == 0:
            raise ZeroDivisionError("QQJet division by zero scalar")
        return QQJet(tuple(c / sc for c in self.coeffs), self.order)

    def __rtruediv__(self, other):
        return QQJet.const(other, self.order) / self

    def __pow__(self, n):
        one = self.coeffs[0] * 0 + 1
        return _pow_ring(self, n, QQJet.const(one, self.order))

    def __repr__(self):
        return f"QQJet({self.coeffs!r}, order={self.order})"


class IdentityQuantizerBundle:
    """Zero-branch identity QuantizerBundle passing symbolic CAS types straight through."""

    def __call__(self, x, *args):
        return x

    def quantize(self, x, *args):
        return x

    def quantize_scalar(self, x, *args):
        return x

    def quantize_vec(self, x, *args):
        return x

    def quantize_state(self, x, *args):
        return x

    def quantize_pair(self, a, b, *args):
        return (a, b)

    def quantize_step(self, x, *args):
        return x

    def clear_denominators(self, x, *args):
        return x

    def project(self, x, *args):
        return x


# ============================================================================
# Exact Primitive Subresultant Sturm Chains & Bisection Root Isolation over Q
# ============================================================================


def _build_raw_subres_sturm(p0):
    d0 = len(p0) - 1
    if d0 == 0:
        return [p0]
    p1 = _prim_part_pos([i * p0[i] for i in range(1, d0 + 1)])
    raw_chain = [p0, p1]
    g = 1
    h = 1
    step = 0
    for _ in range(d0 + 2):
        a = raw_chain[-2]
        b = raw_chain[-1]
        m = len(a) - 1
        n = len(b) - 1
        if n == 0:
            break
        step += 1
        delta = m - n
        denom = 1 if step == 1 else g * (h ** delta)
        if delta == 1:
            pos_rem = _subresultant_step_2adic(a, b, denom)
        else:
            pos_rem = _signed_prem_general(a, b, denom)
        neg_rem = _strip_int_coeffs([-x for x in pos_rem])
        if len(neg_rem) == 1 and neg_rem[0] == 0:
            break
        raw_chain.append(neg_rem)
        bn = b[n]
        bn_pos = bn if bn >= 0 else -bn
        g = bn_pos
        h = (g ** delta) if step == 1 else (g ** delta) // (h ** (delta - 1))
        if len(neg_rem) == 1:
            break
    return raw_chain


def _sturm_chain_raw_int(poly):
    if poly.is_zero():
        return [[0]]
    v0 = 0
    for c in poly.coeffs:
        if c == 0:
            v0 += 1
        else:
            break
    poly_work = (
        QQPoly((Fraction(0, 1),) + poly.coeffs[v0:], poly.var_name)
        if v0 >= 2
        else poly
    )
    p0 = _poly_to_prim_int(poly_work.coeffs)
    if len(p0) <= 1:
        return [p0]
    raw_chain = _build_raw_subres_sturm(p0)
    if len(raw_chain[-1]) > 1:
        g_int = _prim_part_pos(raw_chain[-1])
        if g_int[-1] < 0:
            g_int = [-x for x in g_int]
        g_poly = QQPoly(tuple(Fraction(c, 1) for c in g_int), poly.var_name).monic()
        sqfree = poly_work // g_poly
        p0_sq = _poly_to_prim_int(sqfree.coeffs)
        raw_chain = _build_raw_subres_sturm(p0_sq)
    return raw_chain


def _sturm_chain_prim_int(poly):
    return [_prim_part_pos(p) for p in _sturm_chain_raw_int(poly)]


def sturm_chain(poly):
    """Exact primitive Sturm sequence [P_0, P_1, ..., P_r] of square-free part of poly."""
    int_chain = _sturm_chain_prim_int(poly)
    return [
        QQPoly(tuple(Fraction(c, 1) for c in p), poly.var_name) for p in int_chain
    ]


def _eval_int_poly_sign(coeffs_int, x_frac):
    d = len(coeffs_int) - 1
    if d == 0:
        c0 = coeffs_int[0]
        return 1 if c0 > 0 else (-1 if c0 < 0 else 0)
    p = x_frac.numerator
    q = x_frac.denominator
    if p == 0:
        c0 = coeffs_int[0]
        return 1 if c0 > 0 else (-1 if c0 < 0 else 0)
    if q == 1:
        acc = coeffs_int[d]
        for k in range(d - 1, -1, -1):
            acc = acc * p + coeffs_int[k]
        return 1 if acc > 0 else (-1 if acc < 0 else 0)
    acc = coeffs_int[d]
    q_pow = 1
    for k in range(d - 1, -1, -1):
        q_pow *= q
        acc = acc * p + coeffs_int[k] * q_pow
    return 1 if acc > 0 else (-1 if acc < 0 else 0)


def _count_sign_changes(signs):
    changes = 0
    prev = 0
    for s in signs:
        if s != 0:
            if prev != 0 and s != prev:
                changes += 1
            prev = s
    return changes


def _sturm_var_at_int_chain(int_chain, x_frac):
    return _count_sign_changes([_eval_int_poly_sign(p, x_frac) for p in int_chain])


def _sturm_var_pos_inf_int_chain(int_chain):
    return _count_sign_changes(
        [1 if p[-1] > 0 else (-1 if p[-1] < 0 else 0) for p in int_chain]
    )


def _sturm_var_neg_inf_int_chain(int_chain):
    signs = []
    for p in int_chain:
        d = len(p) - 1
        lc = p[-1]
        s = 1 if lc > 0 else (-1 if lc < 0 else 0)
        if d & 1:
            s = -s
        signs.append(s)
    return _count_sign_changes(signs)


def sturm_variations_at(chain, x_val):
    rx = to_qq(x_val)
    int_chain = [_poly_to_prim_int(p.coeffs) for p in chain]
    return _sturm_var_at_int_chain(int_chain, rx)


def sturm_count_interval(chain, a, b):
    ra = to_qq(a)
    rb = to_qq(b)
    int_chain = [_poly_to_prim_int(p.coeffs) for p in chain]
    return _sturm_var_at_int_chain(int_chain, ra) - _sturm_var_at_int_chain(
        int_chain, rb
    )


def _sturm_count_halfline_pos_int_chain(int_chain):
    p0 = int_chain[0]
    if p0[0] != 0:
        v0 = _sturm_var_at_int_chain(int_chain, Fraction(0, 1))
    else:
        k0 = 0
        for idx in range(len(p0)):
            if p0[idx] != 0:
                k0 = idx
                break
        c_low = q_abs(p0[k0])
        c_top = 0
        for idx in range(k0 + 1, len(p0)):
            av = q_abs(p0[idx])
            if av > c_top:
                c_top = av
        if c_top == 0:
            return 0
        eps_pos = Fraction(c_low, 2 * (c_low + c_top))
        v0 = _sturm_var_at_int_chain(int_chain, eps_pos)
    v_inf = _sturm_var_pos_inf_int_chain(int_chain)
    return v0 - v_inf


def sturm_count_halfline_pos(chain):
    """Exact count of distinct real roots strictly in (0, +infinity)."""
    int_chain = [_poly_to_prim_int(p.coeffs) for p in chain]
    return _sturm_count_halfline_pos_int_chain(int_chain)


def sturm_count_halfline_pos_poly(poly):
    """Fast direct exact count of distinct real roots of poly strictly in (0, +infinity)."""
    return _sturm_count_halfline_pos_int_chain(_sturm_chain_raw_int(poly))


def sturm_count_real(chain):
    """Exact total count of distinct real roots in (-infinity, +infinity)."""
    int_chain = [_poly_to_prim_int(p.coeffs) for p in chain]
    return _sturm_var_neg_inf_int_chain(int_chain) - _sturm_var_pos_inf_int_chain(
        int_chain
    )


def _nudge_off_root(p0, x, step_dir, limit_bound):
    if _eval_int_poly_sign(p0, x) != 0:
        return x
    delta = q_abs(limit_bound - x) * Fraction(1, 1 << 20)
    for _ in range(40):
        cand = x + step_dir * delta
        if _eval_int_poly_sign(p0, cand) != 0:
            return cand
        delta = delta * Fraction(1, 2)
    return x


def isolate_roots_sturm(
    poly,
    lo,
    hi,
    target_width=Fraction(1, 1 << 40),
    max_bisections=80,
):
    """Isolate and refine all distinct real roots of poly in (lo, hi) to certified rational intervals."""
    r_lo = to_qq(lo)
    r_hi = to_qq(hi)
    if r_lo >= r_hi:
        return []
    int_chain = _sturm_chain_raw_int(poly)
    p0 = int_chain[0]
    if len(p0) <= 1:
        return []

    a0 = _nudge_off_root(p0, r_lo, 1, r_hi)
    b0 = _nudge_off_root(p0, r_hi, -1, a0)
    va0 = _sturm_var_at_int_chain(int_chain, a0)
    vb0 = _sturm_var_at_int_chain(int_chain, b0)
    total_roots = va0 - vb0
    if total_roots <= 0:
        return []

    active = [(a0, b0, va0, vb0)]
    isolated = []
    max_iso_steps = total_roots * 128 + 64
    for _ in range(max_iso_steps):
        if len(active) == 0:
            break
        nxt_active = []
        for a, b, va, vb in active:
            cnt = va - vb
            if cnt <= 0:
                continue
            if cnt == 1 and (b - a) <= target_width:
                isolated.append((a, b))
                continue
            mid = (a + b) * Fraction(1, 2)
            if _eval_int_poly_sign(p0, mid) == 0:
                eps = q_min((b - a) * Fraction(1, 8), target_width * Fraction(1, 4))
                for __ in range(64):
                    left_m = mid - eps
                    right_m = mid + eps
                    vl = _sturm_var_at_int_chain(int_chain, left_m)
                    vr = _sturm_var_at_int_chain(int_chain, right_m)
                    sl = _eval_int_poly_sign(p0, left_m)
                    sr = _eval_int_poly_sign(p0, right_m)
                    if vl - vr == 1 and sl * sr < 0:
                        if va - vl > 0:
                            nxt_active.append((a, left_m, va, vl))
                        isolated.append((left_m, right_m))
                        if vr - vb > 0:
                            nxt_active.append((right_m, b, vr, vb))
                        break
                    eps = eps * Fraction(1, 2)
                continue
            vm = _sturm_var_at_int_chain(int_chain, mid)
            if va - vm > 0:
                if va - vm == 1:
                    isolated.append((a, mid))
                else:
                    nxt_active.append((a, mid, va, vm))
            if vm - vb > 0:
                if vm - vb == 1:
                    isolated.append((mid, b))
                else:
                    nxt_active.append((mid, b, vm, vb))
        active = nxt_active

    refined = []
    for a_init, b_init in isolated:
        a = a_init
        b = b_init
        sa = _eval_int_poly_sign(p0, a)
        sb = _eval_int_poly_sign(p0, b)
        for _ in range(max_bisections):
            if (b - a) <= target_width:
                break
            mid = (a + b) * Fraction(1, 2)
            sm = _eval_int_poly_sign(p0, mid)
            if sm == 0:
                mid = (a * 3 + b * 5) * Fraction(1, 8)
                sm = _eval_int_poly_sign(p0, mid)
                if sm == 0:
                    break
            if sa * sm < 0:
                b = mid
                sb = sm
            else:
                a = mid
                sa = sm
        refined.append((a, b))

    for i in range(1, len(refined)):
        cur_iv = refined[i]
        j = i - 1
        for _ in range(i):
            if j >= 0 and refined[j][0] > cur_iv[0]:
                refined[j + 1] = refined[j]
                j -= 1
            else:
                break
        refined[j + 1] = cur_iv
    return refined


# ============================================================================
# Exact Rational Cholesky LDL^T SOS & Lukacs Positivity Certificates
# ============================================================================


def cholesky_ldlt_sos(sym_mat):
    """Exact rational square-root-free symmetric pivoting P^T Q P = L D L^T over Q."""
    m = len(sym_mat)
    work = [[to_qq(sym_mat[i][j]) for j in range(m)] for i in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if work[i][j] != work[j][i]:
                return {
                    "is_psd": False,
                    "is_pd": False,
                    "rank": 0,
                    "pivots": (),
                    "L": [],
                    "perm": tuple(range(m)),
                    "max_reconstruction_err": Fraction(1, 1),
                }

    l_mat = [[Fraction(1 if i == j else 0, 1) for j in range(m)] for i in range(m)]
    d_vec = [Fraction(0, 1) for _ in range(m)]
    perm = list(range(m))
    rank = 0

    for k in range(m):
        piv = k
        best = work[k][k]
        for r in range(k + 1, m):
            if work[r][r] > best:
                best = work[r][r]
                piv = r
        if best < 0:
            return {
                "is_psd": False,
                "is_pd": False,
                "rank": rank,
                "pivots": tuple(d_vec[:rank]),
                "L": l_mat,
                "perm": tuple(perm),
                "max_reconstruction_err": Fraction(0, 1),
            }
        if best == 0:
            for r in range(k, m):
                for s in range(k, m):
                    if work[r][s] != 0:
                        return {
                            "is_psd": False,
                            "is_pd": False,
                            "rank": rank,
                            "pivots": tuple(d_vec[:rank]),
                            "L": l_mat,
                            "perm": tuple(perm),
                            "max_reconstruction_err": Fraction(0, 1),
                        }
            break
        if piv != k:
            work[k], work[piv] = work[piv], work[k]
            for row in work:
                row[k], row[piv] = row[piv], row[k]
            perm[k], perm[piv] = perm[piv], perm[k]
            for c in range(k):
                l_mat[k][c], l_mat[piv][c] = l_mat[piv][c], l_mat[k][c]

        dk = work[k][k]
        d_vec[k] = dk
        inv_dk = Fraction(1, 1) / dk
        rank += 1
        for i in range(k + 1, m):
            lik = work[i][k] * inv_dk
            l_mat[i][k] = lik
            if lik != 0:
                for j in range(k + 1, i + 1):
                    val = work[i][j] - lik * work[j][k]
                    work[i][j] = val
                    work[j][i] = val

    max_err = Fraction(0, 1)
    for i in range(m):
        for j in range(m):
            target = to_qq(sym_mat[perm[i]][perm[j]])
            recon = Fraction(0, 1)
            lim = i if i <= j else j
            if lim >= rank:
                lim = rank - 1
            for k in range(lim + 1):
                recon += l_mat[i][k] * d_vec[k] * l_mat[j][k]
            err = q_abs(target - recon)
            if err > max_err:
                max_err = err

    return {
        "is_psd": True,
        "is_pd": (rank == m),
        "rank": rank,
        "pivots": tuple(d_vec[:rank]),
        "L": l_mat,
        "perm": tuple(perm),
        "max_reconstruction_err": max_err,
    }


def sos_quad_univariate(poly):
    """Exact Lagrange completion-of-squares SOS decomposition for quadratic a*x^2 + b*x + c."""
    d = poly.degree()
    if d > 2:
        raise ValueError("sos_quad_univariate expects degree <= 2")
    c0 = poly.coeffs[0] if len(poly.coeffs) > 0 else Fraction(0, 1)
    c1 = poly.coeffs[1] if len(poly.coeffs) > 1 else Fraction(0, 1)
    c2 = poly.coeffs[2] if len(poly.coeffs) > 2 else Fraction(0, 1)
    x = QQPoly.var(poly.var_name)
    if c2 == 0:
        if c1 != 0 or c0 < 0:
            return {"is_sos": False}
        recon = QQPoly.const(c0, poly.var_name)
        return {
            "is_sos": True,
            "strict_positive": c0 > 0,
            "terms": ((c0, QQPoly.one(poly.var_name)),),
            "min_bound": c0,
            "remainder": poly - recon,
        }
    if c2 < 0:
        return {"is_sos": False}
    shift = c1 / (2 * c2)
    w0 = c0 - (c1 * c1) / (4 * c2)
    if w0 < 0:
        return {"is_sos": False}
    lin = x + shift
    recon = c2 * (lin ** 2) + w0
    return {
        "is_sos": True,
        "strict_positive": (c2 > 0 and w0 > 0),
        "terms": ((c2, lin), (w0, QQPoly.one(poly.var_name))),
        "min_bound": w0,
        "remainder": poly - recon,
    }


def lukacs_positivity_halfline(poly, max_polya_exp=64):
    """Certify poly(u) >= 0 on [0, +infinity) via exact Polya-Lukacs SOS."""
    if poly.is_zero():
        return {
            "certified": True,
            "zero_order_at_0": 0,
            "polya_even_power": 0,
            "even_sos_terms": (),
            "odd_sos_terms": (),
            "remainder": QQPoly.zero(poly.var_name),
        }
    v_zero = 0
    for c in poly.coeffs:
        if c == 0:
            v_zero += 1
        else:
            break
    core = QQPoly(poly.coeffs[v_zero:], poly.var_name)
    if core.coeffs[0] < 0 or core.lc() < 0:
        return {"certified": False}

    one_plus_u_sq = QQPoly(
        (Fraction(1, 1), Fraction(2, 1), Fraction(1, 1)), poly.var_name
    )
    cur = core
    u_var = QQPoly.var(poly.var_name)

    for r in range(max_polya_exp + 1):
        all_nonneg = True
        for c in cur.coeffs:
            if c < 0:
                all_nonneg = False
                break
        if all_nonneg:
            even_terms = []
            odd_terms = []
            for k in range(len(cur.coeffs)):
                ck = cur.coeffs[k]
                if ck > 0:
                    half_k = k // 2
                    mon = u_var ** half_k
                    if (k & 1) == 0:
                        even_terms.append((ck, mon))
                    else:
                        odd_terms.append((ck, mon))
            sos_even = QQPoly.zero(poly.var_name)
            for w, g_p in even_terms:
                sos_even = sos_even + w * (g_p ** 2)
            sos_odd = QQPoly.zero(poly.var_name)
            for w, h_p in odd_terms:
                sos_odd = sos_odd + w * (h_p ** 2)
            lhs = (one_plus_u_sq ** r) * poly
            rhs = (u_var ** v_zero) * (sos_even + u_var * sos_odd)
            return {
                "certified": True,
                "zero_order_at_0": v_zero,
                "polya_even_power": 2 * r,
                "even_sos_terms": tuple(even_terms),
                "odd_sos_terms": tuple(odd_terms),
                "remainder": lhs - rhs,
            }
        cur = cur * one_plus_u_sq

    return {"certified": False}


def lukacs_positivity_interval(poly, a, b, max_polya_exp=64):
    """Certify poly(x) >= 0 on compact interval [a, b] via Mobius-Bernstein-Lukacs SOS."""
    ra = to_qq(a)
    rb = to_qq(b)
    if ra >= rb:
        raise ValueError("lukacs_positivity_interval requires a < b")
    d = poly.degree()
    t = QQPoly.var("t")
    num_t = QQPoly.const(ra, "t") + QQPoly.const(rb, "t") * t
    den_t = QQPoly.one("t") + t
    q_poly = QQPoly.zero("t")
    for k in range(len(poly.coeffs)):
        ck = poly.coeffs[k]
        if ck != 0:
            q_poly = q_poly + ck * (num_t ** k) * (den_t ** (d - k))
    half_cert = lukacs_positivity_halfline(q_poly, max_polya_exp=max_polya_exp)
    if not half_cert.get("certified", False):
        return {"certified": False}
    return {
        "certified": True,
        "interval": (ra, rb),
        "degree": d,
        "mobius_poly": q_poly,
        "halfline_certificate": half_cert,
        "remainder": half_cert["remainder"],
    }


def sturm_critical_intervals():
    """Return Sturm-isolated rational critical intervals for L in {4, 6, 8} 2D Ising tori.

    Machine-verified from scratch in formalist/verify_proofs.py via isolate_roots_sturm.
    """
    return {
        "binder_crossings": {
            "4_6": {
                "u_lo": Fraction(36368867567419, 219902325555200),
                "u_hi": Fraction(4546108445951, 27487790694400),
            },
            "6_8": {
                "u_lo": Fraction(9281490212173, 54975581388800),
                "u_hi": Fraction(37125960848881, 219902325555200),
            },
            "4_8": {
                "u_lo": Fraction(7349514087299, 43980465111040),
                "u_hi": Fraction(9186892609171, 54975581388800),
            },
        },
        "cv_peaks": {
            "L4": {
                "u_lo": Fraction(247094293497643787680373736053187, 1267650600228229401496703205376000),
                "u_hi": Fraction(3953508695980332826001882687426107, 20282409603651670423947251286016000),
            },
            "L6": {
                "u_lo": Fraction(3814243667971012347914642663874381, 20282409603651670423947251286016000),
                "u_hi": Fraction(7628487335975501599997111478870371, 40564819207303340847894502572032000),
            },
            "L8": {
                "u_lo": Fraction(186822579766479361689671151290409, 1014120480182583521197362564300800),
                "u_hi": Fraction(14945806381382504544203846510291081, 81129638414606681695789005144064000),
            },
        },
    }
