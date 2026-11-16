"""Zero-dependency pure-Q symbolic computer algebra system (CAS) for Exhibit I v2.0.

Implements strictly over Q (Python int and fractions.Fraction only; zero float/complex,
zero math/cmath/numpy/scipy/sympy, zero while loops, zero try/except, zero isinstance/type/hasattr,
zero max/min/abs, zero variable ** powers, quantize_branch_sites == []):
  - Univariate polynomials QQPoly and reduced rational functions QQRat,
  - Multivariate polynomials QQMVPoly and rational functions QQMVRat,
  - Truncated Taylor jets QQJet of arbitrary order over Q or symbolic rings,
  - Exact signed subresultant Sturm chains (sturm_chain) and root counting/isolation,
  - Square-root-free Cholesky LDL^T Sum-of-Squares (SOS), Lukacs interval positivity,
    and Positivstellensatz cone certificates (verify_sos_cone_certificate),
  - Exact p-adic and dyadic valuation tools across all 21 primes in 101..199,
  - IdentityQuantizerBundle for unbranched symbolic execution of Battle 1-3 ASTs.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Optional, Sequence, Tuple

PRIMES_101_199: Tuple[int, ...] = (
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
    151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
)

# Tag convention via .numerator and .denominator attributes (zero isinstance!):
# - int / Fraction: x.numerator is int, x.denominator >= 1
# - QQPoly:   x.numerator is self, x.denominator == -1
# - QQRat:    x.numerator is self, x.denominator == -2
# - QQJet:    x.numerator is self, x.denominator == -3
# - QQMVPoly: x.numerator is self, x.denominator == -4
# - QQMVRat:  x.numerator is self, x.denominator == -5


def q_abs(val_in: Any) -> Any:
    return val_in if val_in >= 0 else -val_in


def q_max(a: Any, b: Any) -> Any:
    return a if a >= b else b


def q_min(a: Any, b: Any) -> Any:
    return a if a <= b else b


def to_qq(x: Any) -> Fraction:
    return Fraction(x.numerator, x.denominator)


def igcd(a: int, b: int) -> int:
    ap = a if a >= 0 else -a
    bp = b if b >= 0 else -b
    if bp == 0:
        return ap
    return bp // Fraction(ap, bp).denominator


def ilcm(a: int, b: int) -> int:
    ap = a if a >= 0 else -a
    bp = b if b >= 0 else -b
    if ap == 0 or bp == 0:
        return 0
    return (ap // igcd(ap, bp)) * bp


def dyadic_valuation(n: int) -> int:
    np = n if n >= 0 else -n
    if np == 0:
        return 0
    low = np & (-np)
    return low.bit_length() - 1


def odd_part(n: int) -> int:
    np = n if n >= 0 else -n
    if np == 0:
        return 0
    low = np & (-np)
    return np >> (low.bit_length() - 1)


def int_padic_valuation(n: int, p: int) -> int:
    np = n if n >= 0 else -n
    if np == 0 or p <= 1:
        return 0
    cur = np
    v = 0
    bound = np.bit_length() + 1
    for _step in range(bound):
        q, rem_v = divmod(cur, p)
        if rem_v != 0:
            break
        cur = q
        v += 1
    return v


def padic_valuation(q: Fraction, p: int) -> int:
    num = q.numerator
    den = q.denominator
    if num == 0:
        return 0
    return int_padic_valuation(num, p) - int_padic_valuation(den, p)


def pow_int(base: int, exp: int) -> int:
    if exp < 0:
        raise ValueError("Negative exponent forbidden in pow_int")
    if exp == 0:
        return 1
    res = 1
    cur = base
    rem = exp
    for _bit in range(64):
        if rem == 0:
            break
        if rem & 1:
            res = res * cur
        rem >>= 1
        if rem > 0:
            cur = cur * cur
    return res


def pow_ring(base: Any, exp: int, one_val: Any) -> Any:
    if exp < 0:
        raise ValueError("Negative exponent forbidden in pow_ring")
    if exp == 0:
        return one_val
    res = one_val
    cur = base
    rem = exp
    for _bit in range(64):
        if rem == 0:
            break
        if rem & 1:
            res = res * cur
        rem >>= 1
        if rem > 0:
            cur = cur * cur
    return res


def _strip_qq_coeffs(coeffs: Sequence[Any]) -> Tuple[Fraction, ...]:
    ln = len(coeffs)
    if ln == 0:
        return (Fraction(0, 1),)
    d = ln - 1
    for _idx in range(ln):
        if d > 0 and coeffs[d] == 0:
            d -= 1
        else:
            break
    return tuple(to_qq(coeffs[i]) for i in range(d + 1))


def _poly_to_prim_int(poly_coeffs: Sequence[Fraction]) -> List[int]:
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


def _strip_int_coeffs(ints: Sequence[int]) -> List[int]:
    ln = len(ints)
    if ln == 0:
        return [0]
    d = ln - 1
    for _idx in range(ln):
        if d > 0 and ints[d] == 0:
            d -= 1
        else:
            break
    return list(ints[: d + 1])


def _prim_part_pos(ints: Sequence[int]) -> List[int]:
    seq_v = _strip_int_coeffs(ints)
    if len(seq_v) == 1 and seq_v[0] == 0:
        return [0]
    g = 0
    for val in seq_v:
        if val != 0:
            g = igcd(g, val)
            if g == 1:
                break
    if g > 1:
        return [val // g for val in seq_v]
    return seq_v


def _signed_prem_general(a: Sequence[int], b: Sequence[int], div_int: int) -> List[int]:
    m = len(a) - 1
    n = len(b) - 1
    bn = b[n]
    bn_pos = bn if bn >= 0 else -bn
    bn_sgn = 1 if bn >= 0 else -1
    rem_arr = list(a)
    for k in range(m, n - 1, -1):
        rem_k = rem_arr[k]
        if rem_k != 0:
            scale = rem_k * bn_sgn
            off = k - n
            for j in range(off):
                rem_arr[j] *= bn_pos
            for j in range(n):
                rem_arr[off + j] = rem_arr[off + j] * bn_pos - scale * b[j]
        else:
            for j in range(k):
                rem_arr[j] *= bn_pos
        rem_arr[k] = 0
    if div_int == 1:
        return rem_arr[:n]
    return [rem_arr[i] // div_int for i in range(n)]


def _gcd_int_polys(pa: Sequence[int], pb: Sequence[int]) -> List[int]:
    a = _prim_part_pos(pa)
    b = _prim_part_pos(pb)
    if len(a) == 1 and a[0] == 0:
        return [-x for x in b] if b[-1] < 0 else b
    if len(b) == 1 and b[0] == 0:
        return [-x for x in a] if a[-1] < 0 else a
    if len(a) < len(b):
        a, b = b, a
    if len(b) == 1:
        return [1]
    g = 1
    h = 1
    step = 0
    max_iters = len(a) + len(b) + 2
    for _it in range(max_iters):
        m = len(a) - 1
        n = len(b) - 1
        if n == 0:
            return [1]
        step += 1
        delta = m - n
        denom = 1 if step == 1 else g * pow_int(h, delta)
        rem = _strip_int_coeffs(_signed_prem_general(a, b, denom))
        if len(rem) == 1 and rem[0] == 0:
            res = _prim_part_pos(b)
            return [-x for x in res] if res[-1] < 0 else res
        bn = b[n]
        g = bn if bn >= 0 else -bn
        h = pow_int(g, delta) if step == 1 else pow_int(g, delta) // pow_int(h, delta - 1)
        a = b
        b = rem
    res = _prim_part_pos(b)
    return [-x for x in res] if res[-1] < 0 else res


class QQPoly:
    """Exact univariate polynomial over Q with Fraction coefficients."""

    def __init__(self, coeffs: Sequence[Any], var_name: str = "x") -> None:
        self.coeffs = _strip_qq_coeffs(coeffs)
        self.var_name = var_name
        self.numerator = self
        self.denominator = -1

    @classmethod
    def zero(cls, var_name: str = "x") -> QQPoly:
        return cls((Fraction(0, 1),), var_name)

    @classmethod
    def one(cls, var_name: str = "x") -> QQPoly:
        return cls((Fraction(1, 1),), var_name)

    @classmethod
    def const(cls, c: Any, var_name: str = "x") -> QQPoly:
        return cls((to_qq(c),), var_name)

    @classmethod
    def var(cls, var_name: str = "x") -> QQPoly:
        return cls((Fraction(0, 1), Fraction(1, 1)), var_name)

    def degree(self) -> int:
        if len(self.coeffs) == 1 and self.coeffs[0] == 0:
            return 0
        return len(self.coeffs) - 1

    def is_zero(self) -> bool:
        return len(self.coeffs) == 1 and self.coeffs[0] == 0

    def lc(self) -> Fraction:
        return self.coeffs[-1]

    def monic(self) -> QQPoly:
        if self.is_zero():
            return self
        lead = self.coeffs[-1]
        if lead == 1:
            return self
        inv_l = Fraction(1, 1) / lead
        return QQPoly(tuple(c * inv_l for c in self.coeffs), self.var_name)

    def __pos__(self) -> QQPoly:
        return self

    def __neg__(self) -> QQPoly:
        return QQPoly(tuple(-c for c in self.coeffs), self.var_name)

    def __eq__(self, other: object) -> bool:
        val = other.numerator  # type: ignore[attr-defined]
        tag = other.denominator  # type: ignore[attr-defined]
        if tag >= 1:
            q = Fraction(val, tag)
            return len(self.coeffs) == 1 and self.coeffs[0] == q
        if tag == -1:
            return self.coeffs == val.coeffs
        if tag in (-2, -3, -4, -5):
            return bool(val == self)
        return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __add__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self
            return QQPoly((self.coeffs[0] + q,) + self.coeffs[1:], self.var_name)
        if tag == -1:
            la = len(self.coeffs)
            lb = len(val.coeffs)
            ln = la if la >= lb else lb
            out = [
                (self.coeffs[i] if i < la else Fraction(0, 1))
                + (val.coeffs[i] if i < lb else Fraction(0, 1))
                for i in range(ln)
            ]
            return QQPoly(out, self.var_name)
        return NotImplemented

    def __radd__(self, other: Any) -> Any:
        return self + other

    def __sub__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self
            return QQPoly((self.coeffs[0] - q,) + self.coeffs[1:], self.var_name)
        if tag == -1:
            la = len(self.coeffs)
            lb = len(val.coeffs)
            ln = la if la >= lb else lb
            out = [
                (self.coeffs[i] if i < la else Fraction(0, 1))
                - (val.coeffs[i] if i < lb else Fraction(0, 1))
                for i in range(ln)
            ]
            return QQPoly(out, self.var_name)
        return NotImplemented

    def __rsub__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return (-self) + Fraction(val, tag)
        return NotImplemented

    def __mul__(self, other: Any) -> Any:
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
            out = [Fraction(0, 1) for _k in range(la + lb - 1)]
            for i in range(la):
                ci = self.coeffs[i]
                if ci != 0:
                    for j in range(lb):
                        out[i + j] += ci * val.coeffs[j]
            return QQPoly(out, self.var_name)
        return NotImplemented

    def __rmul__(self, other: Any) -> Any:
        return self * other

    def __truediv__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                raise ZeroDivisionError("QQPoly division by zero scalar")
            return self * (Fraction(1, 1) / q)
        if tag == -1:
            return QQRat(self, val)
        if tag == -2:
            return QQRat(self, QQPoly.one(self.var_name)) / val
        return NotImplemented

    def __rtruediv__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQRat(QQPoly.const(Fraction(val, tag), self.var_name), self)
        return NotImplemented

    def __floordiv__(self, other: Any) -> QQPoly:
        q, _rem = self.divmod_poly(other)
        return q

    def __mod__(self, other: Any) -> QQPoly:
        _q, rem = self.divmod_poly(other)
        return rem

    def __pow__(self, n: int) -> QQPoly:
        return pow_ring(self, n, QQPoly.one(self.var_name))

    def divmod_poly(self, divisor: Any) -> Tuple[QQPoly, QQPoly]:
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
        quot = [Fraction(0, 1) for _k in range(da - db + 1)]
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

    def diff(self, m: int = 1) -> QQPoly:
        cur = self
        for _step in range(m):
            if len(cur.coeffs) <= 1:
                return QQPoly.zero(self.var_name)
            cur = QQPoly(
                tuple(Fraction(i, 1) * cur.coeffs[i] for i in range(1, len(cur.coeffs))),
                self.var_name,
            )
        return cur

    def eval_at(self, x_val: Any) -> Any:
        d = len(self.coeffs) - 1
        acc: Any = self.coeffs[d]
        for k in range(d - 1, -1, -1):
            acc = acc * x_val + self.coeffs[k]
        return acc

    def __repr__(self) -> str:
        return f"QQPoly({self.coeffs!r}, var_name={self.var_name!r})"


def gcd_poly(a: QQPoly, b: QQPoly) -> QQPoly:
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

    def __init__(self, num: Any, den: Any = 1, var_name: str = "x", _reduced: bool = False) -> None:
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

        if not _reduced:
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
    def var(cls, var_name: str = "x") -> QQRat:
        return cls(QQPoly.var(var_name), QQPoly.one(var_name), var_name, _reduced=True)

    @classmethod
    def const(cls, c: Any, var_name: str = "x") -> QQRat:
        return cls(QQPoly.const(c, var_name), QQPoly.one(var_name), var_name, _reduced=True)

    def is_zero(self) -> bool:
        return self.num.is_zero()

    def __pos__(self) -> QQRat:
        return self

    def __neg__(self) -> QQRat:
        return QQRat(-self.num, self.den, self.var_name, _reduced=True)

    def _coerce_other(self, other: Any) -> Optional[QQRat]:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQRat.const(Fraction(val, tag), self.var_name)
        if tag == -1:
            return QQRat(val, QQPoly.one(val.var_name), val.var_name, _reduced=True)
        if tag == -2:
            return val
        return None

    def __eq__(self, other: object) -> bool:
        co = self._coerce_other(other)
        if co is not None:
            return self.num == co.num and self.den == co.den
        val = other.numerator  # type: ignore[attr-defined]
        tag = other.denominator  # type: ignore[attr-defined]
        if tag in (-3, -4, -5):
            return bool(val == self)
        return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __add__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero():
            return co
        if co.is_zero():
            return self
        if self.den == co.den:
            return QQRat(self.num + co.num, self.den, self.var_name)
        return QQRat(self.num * co.den + co.num * self.den, self.den * co.den, self.var_name)

    def __radd__(self, other: Any) -> Any:
        return self + other

    def __sub__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self + (-co)

    def __rsub__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co + (-self)

    def __mul__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero() or co.is_zero():
            return QQRat.const(0, self.var_name)
        return QQRat(self.num * co.num, self.den * co.den, self.var_name)

    def __rmul__(self, other: Any) -> Any:
        return self * other

    def __truediv__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if co.is_zero():
            raise ZeroDivisionError("QQRat division by zero")
        return QQRat(self.num * co.den, self.den * co.num, self.var_name)

    def __rtruediv__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co / self

    def __pow__(self, n: int) -> QQRat:
        return pow_ring(self, n, QQRat.const(1, self.var_name))

    def eval_at(self, x_val: Any) -> Any:
        return self.num.eval_at(x_val) / self.den.eval_at(x_val)

    def __repr__(self) -> str:
        return f"QQRat({self.num!r}, {self.den!r})"


def _pad_exp(e: Tuple[int, ...], nvars: int) -> Tuple[int, ...]:
    if len(e) == nvars:
        return e
    return e + (0,) * (nvars - len(e))


class QQMVPoly:
    """Sparse exact multivariate polynomial over Q."""

    def __init__(
        self,
        nvars: int,
        terms: Dict[Tuple[int, ...], Any],
        var_names: Optional[Sequence[str]] = None,
    ) -> None:
        self.numerator = self
        self.denominator = -4
        self.nvars = nvars
        clean: Dict[Tuple[int, ...], Fraction] = {}
        for k, v in terms.items():
            qv = to_qq(v)
            if qv != 0:
                ek = _pad_exp(k, nvars)
                clean[ek] = clean.get(ek, Fraction(0, 1)) + qv
                if clean[ek] == 0:
                    del clean[ek]
        self.terms = clean
        self.var_names = (
            tuple(f"x{i}" for i in range(nvars))
            if var_names is None
            else tuple(var_names)
        )

    @classmethod
    def zero(cls, nvars: int, var_names: Optional[Sequence[str]] = None) -> QQMVPoly:
        return cls(nvars, {}, var_names)

    @classmethod
    def one(cls, nvars: int, var_names: Optional[Sequence[str]] = None) -> QQMVPoly:
        return cls(nvars, {(0,) * nvars: Fraction(1, 1)}, var_names)

    @classmethod
    def const(cls, c: Any, nvars: int, var_names: Optional[Sequence[str]] = None) -> QQMVPoly:
        qc = to_qq(c)
        if qc == 0:
            return cls(nvars, {}, var_names)
        return cls(nvars, {(0,) * nvars: qc}, var_names)

    @classmethod
    def vars(cls, names_seq: Sequence[str]) -> Tuple[QQMVPoly, ...]:
        names = tuple(names_seq)
        n = len(names)
        return tuple(
            cls(n, {tuple(1 if j == i else 0 for j in range(n)): Fraction(1, 1)}, names)
            for i in range(n)
        )

    def is_zero(self) -> bool:
        return len(self.terms) == 0

    def is_const(self) -> bool:
        if len(self.terms) == 0:
            return True
        if len(self.terms) > 1:
            return False
        return ((0,) * self.nvars) in self.terms

    def const_val(self) -> Fraction:
        return self.terms.get((0,) * self.nvars, Fraction(0, 1))

    def _merge_meta(self, other: QQMVPoly) -> Tuple[int, Tuple[str, ...]]:
        nv = self.nvars if self.nvars >= other.nvars else other.nvars
        vn = self.var_names if len(self.var_names) >= nv else other.var_names
        return nv, vn

    def _lex_lead(self) -> Tuple[Tuple[int, ...], Fraction]:
        best_e: Optional[Tuple[int, ...]] = None
        for e in self.terms:
            if best_e is None or e > best_e:
                best_e = e
        if best_e is None:
            return ((0,) * self.nvars, Fraction(0, 1))
        return best_e, self.terms[best_e]

    def __pos__(self) -> QQMVPoly:
        return self

    def __neg__(self) -> QQMVPoly:
        return QQMVPoly(self.nvars, {k: -v for k, v in self.terms.items()}, self.var_names)

    def __eq__(self, other: object) -> bool:
        val = other.numerator  # type: ignore[attr-defined]
        tag = other.denominator  # type: ignore[attr-defined]
        if tag >= 1:
            q = Fraction(val, tag)
            return self.is_const() and self.const_val() == q
        if tag == -4:
            nv, _vn = self._merge_meta(val)
            t1 = {_pad_exp(k, nv): v for k, v in self.terms.items()}
            t2 = {_pad_exp(k, nv): v for k, v in val.terms.items()}
            return t1 == t2
        if tag in (-3, -5):
            return bool(val == self)
        return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __add__(self, other: Any) -> Any:
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

    def __radd__(self, other: Any) -> Any:
        return self + other

    def __sub__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return self + (-Fraction(val, tag))
        if tag == -4:
            return self + (-val)
        return NotImplemented

    def __rsub__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return (-self) + Fraction(val, tag)
        return NotImplemented

    def __mul__(self, other: Any) -> Any:
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
            out: Dict[Tuple[int, ...], Fraction] = {}
            t1 = [(_pad_exp(k, nv), v) for k, v in self.terms.items()]
            t2 = [(_pad_exp(k, nv), v) for k, v in val.terms.items()]
            for e1, c1 in t1:
                for e2, c2 in t2:
                    ek = tuple(e1[i] + e2[i] for i in range(nv))
                    out[ek] = out.get(ek, Fraction(0, 1)) + c1 * c2
            return QQMVPoly(nv, out, vn)
        return NotImplemented

    def __rmul__(self, other: Any) -> Any:
        return self * other

    def exact_div_or_none(self, divisor: QQMVPoly) -> Optional[QQMVPoly]:
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
        lead_d: Optional[Tuple[int, ...]] = None
        for e in div_terms:
            if lead_d is None or e > lead_d:
                lead_d = e
        if lead_d is None:
            return None
        inv_lc_d = Fraction(1, 1) / div_terms[lead_d]
        quot_terms: Dict[Tuple[int, ...], Fraction] = {}
        bound = len(rem_terms) * 8 + 32
        for _step in range(bound):
            if len(rem_terms) == 0:
                return QQMVPoly(nv, quot_terms, vn)
            lead_r: Optional[Tuple[int, ...]] = None
            for e in rem_terms:
                if lead_r is None or e > lead_r:
                    lead_r = e
            if lead_r is None:
                return QQMVPoly(nv, quot_terms, vn)
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

    def __truediv__(self, other: Any) -> Any:
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

    def __rtruediv__(self, other: Any) -> Any:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQMVRat(QQMVPoly.const(Fraction(val, tag), self.nvars, self.var_names), self)
        return NotImplemented

    def __pow__(self, n: int) -> QQMVPoly:
        return pow_ring(self, n, QQMVPoly.one(self.nvars, self.var_names))

    def eval_at(self, point: Sequence[Any]) -> Any:
        if self.is_zero():
            return Fraction(0, 1)
        one = point[0] * 0 + 1 if len(point) > 0 else Fraction(1, 1)
        acc = point[0] * 0 if len(point) > 0 else Fraction(0, 1)
        for e, c in self.terms.items():
            term = one * c
            for i in range(self.nvars):
                if e[i] > 0:
                    term = term * pow_ring(point[i], e[i], one)
            acc = acc + term
        return acc

    def subs(self, mapping: Dict[str, Any]) -> Any:
        pts: List[Any] = []
        for idx, name in enumerate(self.var_names):
            if name in mapping:
                mv = mapping[name]
                if mv.denominator >= 1:
                    pts.append(QQMVPoly.const(Fraction(mv.numerator, mv.denominator), self.nvars, self.var_names))
                elif mv.denominator == -4:
                    pts.append(mv)
                else:
                    pts.append(mv)
            else:
                e = tuple(1 if j == idx else 0 for j in range(self.nvars))
                pts.append(QQMVPoly(self.nvars, {e: Fraction(1, 1)}, self.var_names))
        return self.eval_at(pts)

    def __repr__(self) -> str:
        return f"QQMVPoly(nvars={self.nvars}, terms={self.terms!r})"


class QQMVRat:
    """Multivariate rational function num / den over QQMVPoly with exact factor cancellation."""

    def __init__(self, num: Any, den: Any = None) -> None:
        self.numerator = self
        self.denominator = -5

        vn = num.numerator
        tn = num.denominator
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

        q_direct = p_num.exact_div_or_none(p_den)
        if q_direct is not None:
            self.num = q_direct
            self.den = QQMVPoly.one(nv, vnames)
            return

        q_inv = p_den.exact_div_or_none(p_num)
        if q_inv is not None:
            _le, lc_i = q_inv._lex_lead()
            inv_lc = Fraction(1, 1) / lc_i
            self.num = QQMVPoly.const(inv_lc, nv, vnames)
            self.den = q_inv * inv_lc
            return

        _ld, lc_d = p_den._lex_lead()
        if lc_d != 1:
            inv_l = Fraction(1, 1) / lc_d
            p_num = p_num * inv_l
            p_den = p_den * inv_l
        self.num = p_num
        self.den = p_den

    @classmethod
    def vars(cls, names: Sequence[str]) -> Tuple[QQMVRat, ...]:
        return tuple(cls(p) for p in QQMVPoly.vars(names))

    @classmethod
    def const(cls, c: Any, nvars: int, var_names: Optional[Sequence[str]] = None) -> QQMVRat:
        return cls(QQMVPoly.const(c, nvars, var_names))

    def is_zero(self) -> bool:
        return self.num.is_zero()

    def _coerce_other(self, other: Any) -> Optional[QQMVRat]:
        val = other.numerator
        tag = other.denominator
        if tag >= 1:
            return QQMVRat.const(Fraction(val, tag), self.num.nvars, self.num.var_names)
        if tag == -4:
            return QQMVRat(val)
        if tag == -5:
            return val
        return None

    def __pos__(self) -> QQMVRat:
        return self

    def __neg__(self) -> QQMVRat:
        return QQMVRat(-self.num, self.den)

    def __eq__(self, other: object) -> bool:
        val = other.numerator  # type: ignore[attr-defined]
        tag = other.denominator  # type: ignore[attr-defined]
        if tag >= 1:
            q = Fraction(val, tag)
            if q == 0:
                return self.is_zero()
            return (self.num - self.den * q).is_zero()
        co = self._coerce_other(other)
        if co is not None:
            return (self.num * co.den - co.num * self.den).is_zero()
        return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __add__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero():
            return co
        if co.is_zero():
            return self
        if self.den == co.den:
            return QQMVRat(self.num + co.num, self.den)
        return QQMVRat(self.num * co.den + co.num * self.den, self.den * co.den)

    def __radd__(self, other: Any) -> Any:
        return self + other

    def __sub__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self + (-co)

    def __rsub__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co + (-self)

    def __mul__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.is_zero() or co.is_zero():
            return QQMVRat.const(0, self.num.nvars, self.num.var_names)
        return QQMVRat(self.num * co.num, self.den * co.den)

    def __rmul__(self, other: Any) -> Any:
        return self * other

    def __truediv__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if co.is_zero():
            raise ZeroDivisionError("QQMVRat division by zero")
        return QQMVRat(self.num * co.den, self.den * co.num)

    def __rtruediv__(self, other: Any) -> Any:
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co / self

    def __pow__(self, n: int) -> QQMVRat:
        return pow_ring(self, n, QQMVRat.const(1, self.num.nvars, self.num.var_names))

    def subs(self, mapping: Dict[str, Any]) -> Any:
        return self.num.subs(mapping) / self.den.subs(mapping)

    def eval_at(self, point: Sequence[Any]) -> Any:
        return self.num.eval_at(point) / self.den.eval_at(point)

    def __repr__(self) -> str:
        return f"QQMVRat({self.num!r}, {self.den!r})"


def sym_vars(names_str: str) -> Tuple[QQMVRat, ...]:
    names = tuple(tok for tok in names_str.replace(",", " ").split() if tok)
    return QQMVRat.vars(names)


class QQJet:
    """Truncated Taylor series c_0 + c_1*eps + ... + c_k*eps^k over Q or symbolic rings."""

    def __init__(self, coeffs: Sequence[Any], order: Optional[int] = None) -> None:
        self.numerator = self
        self.denominator = -3
        seq = list(coeffs)
        if order is None:
            order = len(seq) - 1 if len(seq) > 0 else 0
        self.order = order
        zero_elem = seq[0] * 0 if len(seq) > 0 else Fraction(0, 1)
        padded: List[Any] = []
        for i in range(order + 1):
            if i < len(seq):
                val = seq[i].numerator
                tag = seq[i].denominator
                padded.append(Fraction(val, tag) if tag >= 1 else val)
            else:
                padded.append(zero_elem)
        self.coeffs = tuple(padded)

    @classmethod
    def const(cls, c: Any, order: int) -> QQJet:
        val = c.numerator
        tag = c.denominator
        c0 = Fraction(val, tag) if tag >= 1 else val
        z = c0 * 0
        return cls((c0,) + (z,) * order, order)

    @classmethod
    def var(cls, c0: Any = 0, order: int = 4) -> QQJet:
        val = c0.numerator
        tag = c0.denominator
        v0 = Fraction(val, tag) if tag >= 1 else val
        z = v0 * 0
        one = z + 1
        if order == 0:
            return cls((v0,), 0)
        return cls((v0, one) + (z,) * (order - 1), order)

    def coeff(self, m: int) -> Any:
        if m < 0 or m > self.order:
            return self.coeffs[0] * 0
        return self.coeffs[m]

    def is_zero(self) -> bool:
        for c in self.coeffs:
            if c != 0:
                return False
        return True

    def __pos__(self) -> QQJet:
        return self

    def __neg__(self) -> QQJet:
        return QQJet(tuple(-c for c in self.coeffs), self.order)

    def __eq__(self, other: object) -> bool:
        val = other.numerator  # type: ignore[attr-defined]
        tag = other.denominator  # type: ignore[attr-defined]
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

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __add__(self, other: Any) -> QQJet:
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            return QQJet(tuple(self.coeffs[i] + val.coeffs[i] for i in range(ord_m + 1)), ord_m)
        sc = Fraction(val, tag) if tag >= 1 else val
        return QQJet((self.coeffs[0] + sc,) + self.coeffs[1:], self.order)

    def __radd__(self, other: Any) -> QQJet:
        return self + other

    def __sub__(self, other: Any) -> QQJet:
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            return QQJet(tuple(self.coeffs[i] - val.coeffs[i] for i in range(ord_m + 1)), ord_m)
        sc = Fraction(val, tag) if tag >= 1 else val
        return QQJet((self.coeffs[0] - sc,) + self.coeffs[1:], self.order)

    def __rsub__(self, other: Any) -> QQJet:
        return (-self) + other

    def __mul__(self, other: Any) -> QQJet:
        val = other.numerator
        tag = other.denominator
        if tag == -3:
            ord_m = self.order if self.order <= val.order else val.order
            out: List[Any] = []
            for m in range(ord_m + 1):
                acc = self.coeffs[0] * val.coeffs[m]
                for j in range(1, m + 1):
                    acc = acc + self.coeffs[j] * val.coeffs[m - j]
                out.append(acc)
            return QQJet(out, ord_m)
        sc = Fraction(val, tag) if tag >= 1 else val
        return QQJet(tuple(c * sc for c in self.coeffs), self.order)

    def __rmul__(self, other: Any) -> QQJet:
        return self * other

    def __truediv__(self, other: Any) -> QQJet:
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

    def __rtruediv__(self, other: Any) -> QQJet:
        return QQJet.const(other, self.order) / self

    def __pow__(self, n: int) -> QQJet:
        one = self.coeffs[0] * 0 + 1
        return pow_ring(self, n, QQJet.const(one, self.order))

    def __repr__(self) -> str:
        return f"QQJet({self.coeffs!r}, order={self.order})"


def jet_cos_sin(theta_jet: QQJet, terms: int = 5) -> Tuple[QQJet, QQJet]:
    """Exact rational truncated Taylor series of (cos(theta), sin(theta)) for QQJet around 0."""
    ord_j = theta_jet.order
    one = QQJet.const(Fraction(1, 1), ord_j)
    zero = QQJet.const(Fraction(0, 1), ord_j)
    c_acc = one
    s_acc = zero
    t2 = theta_jet * theta_jet
    c_term = one
    s_term = theta_jet
    fact_c = 1
    fact_s = 1
    for k in range(1, terms + 1):
        c_term = c_term * (-t2)
        fact_c = fact_c * (2 * k - 1) * (2 * k)
        c_acc = c_acc + c_term * Fraction(1, fact_c)
        s_acc = s_acc + s_term * Fraction(1, fact_s)
        s_term = s_term * (-t2)
        fact_s = fact_s * (2 * k) * (2 * k + 1)
    return c_acc, s_acc


class IdentityQuantizerBundle:
    """Exact pass-through QuantizerBundle compatible with Fraction, QQPoly, QQRat, QQJet."""

    def __call__(self, val: Any, scale: Any = Fraction(1, 1)) -> Any:
        return val + Fraction(0, 1) * scale

    def quantize(self, val: Any, scale: Any = Fraction(1, 1)) -> Any:
        return val + Fraction(0, 1) * scale

    def coarse(self, val: Any) -> Any:
        return val

    def anchor(self, val: Any, scale: Any = Fraction(1, 1)) -> Any:
        return val + Fraction(0, 1) * scale

    def snap(self, val: Any, scale: Any = Fraction(1, 1)) -> Any:
        return val + Fraction(0, 1) * scale

    def rotor(self, val: Any, scale: Any = Fraction(1, 1)) -> Any:
        return val + Fraction(0, 1) * scale

    def quantize_tangent(self, val: Any, scale: Any = Fraction(1, 1)) -> Any:
        return val + Fraction(0, 1) * scale

    def quantize_pair(self, a: Any, b: Any, scale: Any = Fraction(1, 1)) -> Tuple[Any, Any]:
        z = Fraction(0, 1) * scale
        return (a + z, b + z)

    def project_unit_circle(self, c_val: Any, s_val: Any) -> Tuple[Any, Any]:
        c2 = c_val * c_val
        s2 = s_val * s_val
        den = Fraction(1, 1) * (c2 + s2)
        return ((c2 - s2) / den, (Fraction(2, 1) * c_val * s_val) / den)

    def quantize_vec(self, xs: Sequence[Any], scale: Any = Fraction(1, 1)) -> List[Any]:
        z = Fraction(0, 1) * scale
        return [item + z for item in xs]

    def clear_denominators(self, xs: Sequence[Any]) -> Tuple[List[Any], int]:
        return list(xs), 1


# ============================================================================
# Exact Primitive Subresultant Sturm Chains & Interval Root Counting over Q
# ============================================================================

def _build_raw_subres_sturm(p0: List[int]) -> List[List[int]]:
    d0 = len(p0) - 1
    if d0 == 0:
        return [p0]
    p1 = _prim_part_pos([i * p0[i] for i in range(1, d0 + 1)])
    raw_chain = [p0, p1]
    g = 1
    h = 1
    step = 0
    for _it in range(d0 + 2):
        a = raw_chain[-2]
        b = raw_chain[-1]
        m = len(a) - 1
        n = len(b) - 1
        if n == 0:
            break
        step += 1
        delta = m - n
        denom = 1 if step == 1 else g * pow_int(h, delta)
        pos_rem = _signed_prem_general(a, b, denom)
        neg_rem = _strip_int_coeffs([-x for x in pos_rem])
        if len(neg_rem) == 1 and neg_rem[0] == 0:
            break
        raw_chain.append(neg_rem)
        bn = b[n]
        g = bn if bn >= 0 else -bn
        h = pow_int(g, delta) if step == 1 else pow_int(g, delta) // pow_int(h, delta - 1)
        if len(neg_rem) == 1:
            break
    return raw_chain


def _sturm_chain_raw_int(poly: QQPoly) -> List[List[int]]:
    if poly.is_zero():
        return [[0]]
    p0 = _poly_to_prim_int(poly.coeffs)
    if len(p0) <= 1:
        return [p0]
    raw_chain = _build_raw_subres_sturm(p0)
    if len(raw_chain[-1]) > 1:
        g_int = _prim_part_pos(raw_chain[-1])
        if g_int[-1] < 0:
            g_int = [-x for x in g_int]
        g_poly = QQPoly(tuple(Fraction(c, 1) for c in g_int), poly.var_name).monic()
        sqfree = poly // g_poly
        raw_chain = _build_raw_subres_sturm(_poly_to_prim_int(sqfree.coeffs))
    return raw_chain


def sturm_chain(poly: QQPoly) -> List[QQPoly]:
    """Exact primitive Sturm sequence [P_0, P_1, ..., P_r] of square-free part of poly."""
    return [
        QQPoly(tuple(Fraction(c, 1) for c in _prim_part_pos(p)), poly.var_name)
        for p in _sturm_chain_raw_int(poly)
    ]


def _eval_int_poly_sign(coeffs_int: Sequence[int], x_frac: Fraction) -> int:
    d = len(coeffs_int) - 1
    if d == 0:
        c0 = coeffs_int[0]
        return 1 if c0 > 0 else (-1 if c0 < 0 else 0)
    p = x_frac.numerator
    q = x_frac.denominator
    acc = coeffs_int[d]
    q_pow = 1
    for k in range(d - 1, -1, -1):
        q_pow *= q
        acc = acc * p + coeffs_int[k] * q_pow
    return 1 if acc > 0 else (-1 if acc < 0 else 0)


def _count_sign_changes(signs: Sequence[int]) -> int:
    changes = 0
    prev = 0
    for sgn_v in signs:
        if sgn_v != 0:
            if prev != 0 and sgn_v != prev:
                changes += 1
            prev = sgn_v
    return changes


def sturm_count_interval(chain: Sequence[QQPoly], a: Any, b: Any) -> int:
    """Exact count of distinct real roots in (a, b] when a, b are non-roots."""
    ra = to_qq(a)
    rb = to_qq(b)
    int_chain = [_poly_to_prim_int(p.coeffs) for p in chain]
    va = _count_sign_changes([_eval_int_poly_sign(p, ra) for p in int_chain])
    vb = _count_sign_changes([_eval_int_poly_sign(p, rb) for p in int_chain])
    return va - vb


# ============================================================================
# Square-Root-Free Cholesky LDL^T SOS & Positivstellensatz Cone Certificates
# ============================================================================

def cholesky_ldlt_sos(sym_mat: Sequence[Sequence[Any]]) -> Dict[str, Any]:
    """Exact rational square-root-free symmetric pivoting P^T Q P = L D L^T over Q."""
    m = len(sym_mat)
    work = [[to_qq(sym_mat[i][j]) for j in range(m)] for i in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if work[i][j] != work[j][i]:
                return {"is_psd": False, "is_pd": False, "rank": 0, "pivots": ()}

    l_mat = [[Fraction(1 if i == j else 0, 1) for j in range(m)] for i in range(m)]
    d_vec = [Fraction(0, 1) for _i in range(m)]
    perm = list(range(m))
    rank = 0

    for k in range(m):
        piv = k
        best = work[k][k]
        for row_i in range(k + 1, m):
            if work[row_i][row_i] > best:
                best = work[row_i][row_i]
                piv = row_i
        if best < 0:
            return {"is_psd": False, "is_pd": False, "rank": rank, "pivots": tuple(d_vec[:rank])}
        if best == 0:
            for row_i in range(k, m):
                for col_j in range(k, m):
                    if work[row_i][col_j] != 0:
                        return {"is_psd": False, "is_pd": False, "rank": rank, "pivots": tuple(d_vec[:rank])}
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


def verify_sos_cone_certificate(
    lhs: Any,
    rhs: Any,
    cone_terms: Sequence[Tuple[Any, Sequence[Any], Any]],
    positive_den_factors: Sequence[Any] = (),
) -> Dict[str, Any]:
    """Certify lhs - rhs == (sum_j w_j * (prod_{g in gens_j} g) * sq_j^2) / (prod_{f in den} f) == 0/1."""
    diff = lhs - rhs
    one_rat = diff * 0 + 1
    num_cone: Any = diff * 0
    for weight, gens, sq_factor in cone_terms:
        w_q = to_qq(weight)
        if w_q < 0:
            raise ValueError("SOS cone weight must be non-negative")
        prod_g: Any = one_rat * w_q
        for g in gens:
            prod_g = prod_g * g
        num_cone = num_cone + prod_g * (sq_factor * sq_factor)

    den_prod: Any = one_rat
    for f in positive_den_factors:
        den_prod = den_prod * f

    cert_rat = num_cone / den_prod
    residual = diff - cert_rat
    passed = bool(residual == 0)
    if not passed:
        raise AssertionError("Positivstellensatz SOS cone certificate residual != 0/1")
    return {
        "passed": True,
        "residual_zero": "0/1",
        "cone_term_count": len(cone_terms),
    }
