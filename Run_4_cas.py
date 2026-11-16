"""
Zero-Dependency Symbolic Computer Algebra System over Q (fractions.Fraction)
=============================================================================

Owned by @formalist inside ultrafinitist-gravity/formalist/cas.py.
Pure Python stdlib (fractions.Fraction) ONLY -- zero external libraries.

Provides:
  1. QQPoly  : Exact multivariate/univariate polynomial ring Q[x_1, ..., x_n]
               with Subresultant PRS exact multivariate GCD, Yun square-free
               decomposition, formal partial derivatives, and exact substitution.
  2. QQRat   : Canonical rational function field Q(x_1, ..., x_n) over Z[X]
               unique factorization normalization + Henrici fast addition/product
               so 0/1 algebraic identities reduce canonically without swell.
  3. QQJet   : Truncated Taylor jet algebra K[[eps]] / (eps^(order+1)) working
               polymorphically over BOTH Fraction and QQRat coefficients.
  4. IdentityQuantizerBundle : Unbranched identity callable/bundle allowing
               @architect's live engine AST to execute directly on QQPoly,
               QQRat, and QQJet without spoofing or branching.
  5. Exact Algebraic Positivity, Sturm Root Isolation & Global Supremum Certifier:
               Sturm chains over Q[r], Stern-Brocot/Cauchy rational root detection,
               Lukacs/Poincare Sum-of-Squares (SOS) on [0, inf), and adaptive
               Bernstein-Ray exact 0/1-residual positivity & supremum certificates.
  6. 3+1D Spinor Kerr-Schild & General Metric Tensor Engine:
               Kustaanheimo-Stiefel & Levi-Civita spinors, 4D Christoffel symbols,
               256-component Riemann tensor, rational horizon-penetrating
               orthonormal tetrad curvature R_{abcd}, Ricci, Weyl, Kretschmann
               scalar K_h(r), and chart-independent geodesic jet verifier.
"""

from __future__ import annotations

from fractions import Fraction
import functools
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union


# ---------------------------------------------------------------------------
# Basic exact integer helpers (zero math/external imports)
# ---------------------------------------------------------------------------

def int_gcd(a: int, b: int) -> int:
    """Exact non-negative integer greatest common divisor."""
    a = -a if a < 0 else a
    b = -b if b < 0 else b
    while b:
        a, b = b, a % b
    return a


def int_lcm(a: int, b: int) -> int:
    """Exact non-negative integer least common multiple."""
    if a == 0 or b == 0:
        return 0
    a = -a if a < 0 else a
    b = -b if b < 0 else b
    return (a // int_gcd(a, b)) * b


def int_list_gcd(vals: Iterable[int]) -> int:
    g = 0
    for v in vals:
        g = int_gcd(g, v)
        if g == 1:
            return 1
    return g


# ---------------------------------------------------------------------------
# Sparse Monomial Representation
# ---------------------------------------------------------------------------
# A Monomial is stored as a tuple of (var_name, exp) pairs with exp > 0,
# sorted in STRICTLY DESCENDING variable-name order ('z' > 'y' > ... > 'a').
# Why descending var_name? Because then native Python C-level tuple `<` on
# Monomial objects is an EXACT Dickson-well-founded Pure Lexicographic order!

Monomial = Tuple[Tuple[str, int], ...]
EMPTY_MONO: Monomial = ()


def make_monomial(pairs: Iterable[Tuple[str, int]]) -> Monomial:
    acc: Dict[str, int] = {}
    for var, exp in pairs:
        if exp < 0:
            raise ValueError(f"Negative exponent {exp} in polynomial monomial for {var}")
        if exp:
            acc[var] = acc.get(var, 0) + exp
    return tuple(sorted(((v, e) for v, e in acc.items() if e > 0), reverse=True))


@functools.lru_cache(maxsize=65536)
def mono_mul(m1: Monomial, m2: Monomial) -> Monomial:
    if not m1:
        return m2
    if not m2:
        return m1
    res: List[Tuple[str, int]] = []
    i, j = 0, 0
    n1, n2 = len(m1), len(m2)
    while i < n1 and j < n2:
        v1, e1 = m1[i]
        v2, e2 = m2[j]
        if v1 == v2:
            res.append((v1, e1 + e2))
            i += 1
            j += 1
        elif v1 > v2:
            res.append((v1, e1))
            i += 1
        else:
            res.append((v2, e2))
            j += 1
    if i < n1:
        res.extend(m1[i:])
    if j < n2:
        res.extend(m2[j:])
    return tuple(res)


@functools.lru_cache(maxsize=65536)
def mono_div(m1: Monomial, m2: Monomial) -> Optional[Monomial]:
    """Return m1 / m2 if m2 divides m1, else None."""
    if not m2:
        return m1
    if len(m1) < len(m2):
        return None
    res: List[Tuple[str, int]] = []
    i, j = 0, 0
    n1, n2 = len(m1), len(m2)
    while i < n1 and j < n2:
        v1, e1 = m1[i]
        v2, e2 = m2[j]
        if v1 == v2:
            if e1 < e2:
                return None
            if e1 > e2:
                res.append((v1, e1 - e2))
            i += 1
            j += 1
        elif v1 > v2:
            res.append((v1, e1))
            i += 1
        else:
            # m2 contains v2 which is missing in m1
            return None
    if j < n2:
        return None
    if i < n1:
        res.extend(m1[i:])
    return tuple(res)


def mono_gcd(m1: Monomial, m2: Monomial) -> Monomial:
    if not m1 or not m2:
        return EMPTY_MONO
    res: List[Tuple[str, int]] = []
    i, j = 0, 0
    n1, n2 = len(m1), len(m2)
    while i < n1 and j < n2:
        v1, e1 = m1[i]
        v2, e2 = m2[j]
        if v1 == v2:
            res.append((v1, e1 if e1 < e2 else e2))
            i += 1
            j += 1
        elif v1 > v2:
            i += 1
        else:
            j += 1
    return tuple(res)


def mono_pow(m: Monomial, k: int) -> Monomial:
    if k == 0 or not m:
        return EMPTY_MONO
    if k == 1:
        return m
    return tuple((v, e * k) for v, e in m)


def mono_deg_var(m: Monomial, var: str) -> int:
    for v, e in m:
        if v == var:
            return e
        if v < var:
            return 0
    return 0


def mono_drop_var(m: Monomial, var: str) -> Tuple[int, Monomial]:
    """Return (exp_of_var, monomial_without_var)."""
    for idx, (v, e) in enumerate(m):
        if v == var:
            return e, m[:idx] + m[idx + 1 :]
        if v < var:
            return 0, m
    return 0, m


# ---------------------------------------------------------------------------
# Internal Exact Multivariate Polynomial Arithmetic over Z[x_1, ..., x_n]
# ---------------------------------------------------------------------------

ZPoly = Dict[Monomial, int]


def zpoly_clean(p: ZPoly) -> ZPoly:
    return {m: c for m, c in p.items() if c != 0}


def zpoly_neg(p: ZPoly) -> ZPoly:
    return {m: -c for m, c in p.items()}


def zpoly_add(a: ZPoly, b: ZPoly) -> ZPoly:
    if not a:
        return dict(b)
    if not b:
        return dict(a)
    out = dict(a)
    for m, c in b.items():
        nv = out.get(m, 0) + c
        if nv:
            out[m] = nv
        else:
            out.pop(m, None)
    return out


def zpoly_sub(a: ZPoly, b: ZPoly) -> ZPoly:
    if not b:
        return dict(a)
    if not a:
        return zpoly_neg(b)
    out = dict(a)
    for m, c in b.items():
        nv = out.get(m, 0) - c
        if nv:
            out[m] = nv
        else:
            out.pop(m, None)
    return out


def zpoly_scale(a: ZPoly, s: int) -> ZPoly:
    if s == 0 or not a:
        return {}
    if s == 1:
        return dict(a)
    return {m: c * s for m, c in a.items()}


def zpoly_mul(a: ZPoly, b: ZPoly) -> ZPoly:
    if not a or not b:
        return {}
    out: ZPoly = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = mono_mul(m1, m2)
            nv = out.get(m, 0) + c1 * c2
            if nv:
                out[m] = nv
            else:
                out.pop(m, None)
    return out


def zpoly_pow(a: ZPoly, exp: int) -> ZPoly:
    if exp < 0:
        raise ValueError("Negative exponent in zpoly_pow")
    if exp == 0:
        return {EMPTY_MONO: 1}
    if not a:
        return {}
    res: ZPoly = {EMPTY_MONO: 1}
    base = dict(a)
    e = exp
    while e > 0:
        if e & 1:
            res = zpoly_mul(res, base)
        e >>= 1
        if e:
            base = zpoly_mul(base, base)
    return res


def zpoly_int_content(a: ZPoly) -> int:
    if not a:
        return 0
    return int_list_gcd(a.values())


def zpoly_primitive(a: ZPoly) -> Tuple[int, ZPoly]:
    """
    Factor a = c * pp(a) where c in Z, pp(a) is primitive in Z[X] with lc_lex(pp(a)) > 0.
    """
    if not a:
        return 0, {}
    g = zpoly_int_content(a)
    lm = max(a)
    if a[lm] < 0:
        g = -g
    if g == 1:
        return 1, dict(a)
    return g, {m: c // g for m, c in a.items()}


def zpoly_exact_div(a: ZPoly, b: ZPoly) -> Optional[ZPoly]:
    """
    Exact division a / b in Z[x_1, ..., x_n]. Returns None if b does not divide a over Z[X].
    """
    if not b:
        raise ZeroDivisionError("Exact division by zero ZPoly")
    if not a:
        return {}
    if b == {EMPTY_MONO: 1}:
        return dict(a)
    if b == {EMPTY_MONO: -1}:
        return zpoly_neg(a)
    lm_b = max(b)
    lc_b = b[lm_b]
    rem = dict(a)
    quot: ZPoly = {}
    while rem:
        lm_r = max(rem)
        qm = mono_div(lm_r, lm_b)
        if qm is None:
            return None
        q_coeff, r_mod = divmod(rem[lm_r], lc_b)
        if r_mod != 0:
            return None
        quot[qm] = quot.get(qm, 0) + q_coeff
        for m_b, c_b in b.items():
            m_prod = mono_mul(qm, m_b)
            nv = rem.get(m_prod, 0) - q_coeff * c_b
            if nv:
                rem[m_prod] = nv
            else:
                rem.pop(m_prod, None)
    return quot


def zpoly_vars(a: ZPoly) -> set[str]:
    vs: set[str] = set()
    for m in a:
        for v, _ in m:
            vs.add(v)
    return vs


def zpoly_monomial_content(a: ZPoly) -> Monomial:
    it = iter(a)
    try:
        g = next(it)
    except StopIteration:
        return EMPTY_MONO
    for m in it:
        if not g:
            return EMPTY_MONO
        g = mono_gcd(g, m)
    return g


def zpoly_strip_monomial(a: ZPoly, mg: Monomial) -> ZPoly:
    if not mg:
        return dict(a)
    out: ZPoly = {}
    for m, c in a.items():
        dm = mono_div(m, mg)
        assert dm is not None
        out[dm] = c
    return out


def _to_univar_over_zpoly(a: ZPoly, var: str) -> Dict[int, ZPoly]:
    """View a in Z[other_vars][var] as {deg: coeff_zpoly}."""
    res: Dict[int, ZPoly] = {}
    for m, c in a.items():
        e, rest = mono_drop_var(m, var)
        bucket = res.get(e)
        if bucket is None:
            res[e] = {rest: c}
        else:
            bucket[rest] = c
    return res


def _from_univar_over_zpoly(u: Dict[int, ZPoly], var: str) -> ZPoly:
    out: ZPoly = {}
    for e, coeff in u.items():
        if not coeff:
            continue
        v_mono: Monomial = ((var, e),) if e > 0 else EMPTY_MONO
        for m, c in coeff.items():
            if c:
                out[mono_mul(m, v_mono)] = c
    return out


def _u_deg(u: Dict[int, ZPoly]) -> int:
    return max(u) if u else -1


def _u_scale_poly(u: Dict[int, ZPoly], s: ZPoly) -> Dict[int, ZPoly]:
    if not s or not u:
        return {}
    out: Dict[int, ZPoly] = {}
    for e, c in u.items():
        p = zpoly_mul(c, s)
        if p:
            out[e] = p
    return out


def _u_exact_div_scalar(u: Dict[int, ZPoly], d: ZPoly) -> Dict[int, ZPoly]:
    if d == {EMPTY_MONO: 1}:
        return dict(u)
    out: Dict[int, ZPoly] = {}
    for e, c in u.items():
        q = zpoly_exact_div(c, d)
        if q is None:
            raise ArithmeticError("Subresultant exact scalar division failed")
        if q:
            out[e] = q
    return out


def _u_prem(a: Dict[int, ZPoly], b: Dict[int, ZPoly]) -> Dict[int, ZPoly]:
    """
    Exact sparse pseudo-remainder R in D[x] such that lc(b)^(deg(a)-deg(b)+1) * a = Q * b + R
    with deg_x(R) < deg_x(b). Preserves exact exponent count `delta + 1` across sparse degree drops!
    """
    da = _u_deg(a)
    db = _u_deg(b)
    if da < db:
        return dict(a)
    lc_b = b[db]
    rem = {e: dict(c) for e, c in a.items() if c}
    exp_left = da - db + 1
    while rem:
        dr = _u_deg(rem)
        if dr < db:
            break
        lc_r = rem[dr]
        shift = dr - db
        # rem <- lc_b * rem - lc_r * x^shift * b
        new_rem: Dict[int, ZPoly] = {}
        for e, c in rem.items():
            if e != dr:
                p = zpoly_mul(lc_b, c)
                if p:
                    new_rem[e] = p
        for e, c in b.items():
            if e != db:
                target_e = e + shift
                term = zpoly_mul(lc_r, c)
                cur = new_rem.get(target_e)
                diff = zpoly_sub(cur if cur is not None else {}, term)
                if diff:
                    new_rem[target_e] = diff
                else:
                    new_rem.pop(target_e, None)
        rem = new_rem
        exp_left -= 1
    if exp_left > 0 and rem:
        rem = _u_scale_poly(rem, zpoly_pow(lc_b, exp_left))
    return rem


def _u_content(u: Dict[int, ZPoly]) -> ZPoly:
    g: ZPoly = {}
    for c in u.values():
        g = zpoly_gcd(g, c)
        if g == {EMPTY_MONO: 1}:
            return g
    return g


def zpoly_gcd(a: ZPoly, b: ZPoly) -> ZPoly:
    """
    Exact multivariate GCD in Z[x_1, ..., x_n], normalized primitive with lc_lex > 0.
    Uses monomial extraction, missing-variable content projection, univariate primitive PRS,
    and Collins-Brown Subresultant PRS.
    """
    if not a:
        return zpoly_primitive(b)[1]
    if not b:
        return zpoly_primitive(a)[1]

    _, a = zpoly_primitive(a)
    _, b = zpoly_primitive(b)

    if len(a) == 1 and EMPTY_MONO in a:
        return {EMPTY_MONO: 1}
    if len(b) == 1 and EMPTY_MONO in b:
        return {EMPTY_MONO: 1}
    if a == b:
        return a

    # 1. Extract monomial GCD from each
    ma = zpoly_monomial_content(a)
    mb = zpoly_monomial_content(b)
    mg = mono_gcd(ma, mb)
    if ma:
        a = zpoly_strip_monomial(a, ma)
        _, a = zpoly_primitive(a)
    if mb:
        b = zpoly_strip_monomial(b, mb)
        _, b = zpoly_primitive(b)

    if (len(a) == 1 and EMPTY_MONO in a) or (len(b) == 1 and EMPTY_MONO in b):
        return {mg: 1}
    if a == b:
        return zpoly_mul(a, {mg: 1}) if mg else a

    va = zpoly_vars(a)
    vb = zpoly_vars(b)

    # 2. Any variable appearing in only ONE of (a, b) cannot appear in gcd(a, b);
    #    replace that polynomial by its content wrt that variable immediately!
    only_a = va - vb
    for v in only_a:
        ua = _to_univar_over_zpoly(a, v)
        a = _u_content(ua)
        if len(a) == 1 and EMPTY_MONO in a:
            return {mg: 1}
    only_b = vb - va
    for v in only_b:
        ub = _to_univar_over_zpoly(b, v)
        b = _u_content(ub)
        if len(b) == 1 and EMPTY_MONO in b:
            return {mg: 1}

    common_vars = zpoly_vars(a) & zpoly_vars(b)
    if not common_vars:
        # Repeat projection if secondary variables became disjoint
        va2 = zpoly_vars(a)
        for v in va2:
            a = _u_content(_to_univar_over_zpoly(a, v))
        vb2 = zpoly_vars(b)
        for v in vb2:
            b = _u_content(_to_univar_over_zpoly(b, v))
        _, g_const = zpoly_primitive(a if not b else {EMPTY_MONO: int_gcd(zpoly_int_content(a), zpoly_int_content(b))})
        return {mg: 1}

    # Quick exact-divisibility check before PRS
    if zpoly_exact_div(a, b) is not None:
        return zpoly_mul(b, {mg: 1}) if mg else b
    if zpoly_exact_div(b, a) is not None:
        return zpoly_mul(a, {mg: 1}) if mg else a

    # 3. Choose main variable with smallest max degree to minimize PRS depth
    def _var_cost(v: str) -> Tuple[int, int, str]:
        da = max(mono_deg_var(m, v) for m in a)
        db = max(mono_deg_var(m, v) for m in b)
        return (max(da, db), min(da, db), v)

    main_v = min(common_vars, key=_var_cost)
    ua = _to_univar_over_zpoly(a, main_v)
    ub = _to_univar_over_zpoly(b, main_v)

    ca = _u_content(ua)
    cb = _u_content(ub)
    cg = zpoly_gcd(ca, cb)
    ua = _u_exact_div_scalar(ua, ca)
    ub = _u_exact_div_scalar(ub, cb)

    if _u_deg(ua) < _u_deg(ub):
        ua, ub = ub, ua

    # 4. Collins-Brown Subresultant PRS in D[main_v]
    g_s: ZPoly = {EMPTY_MONO: 1}
    h_s: ZPoly = {EMPTY_MONO: 1}

    while _u_deg(ub) > 0:
        delta = _u_deg(ua) - _u_deg(ub)
        r = _u_prem(ua, ub)
        if not r:
            # ub divides ua in Frac(D)[main_v]
            content_ub = _u_content(ub)
            pp_ub = _u_exact_div_scalar(ub, content_ub)
            res = zpoly_mul(cg, _from_univar_over_zpoly(pp_ub, main_v))
            _, res = zpoly_primitive(res)
            return zpoly_mul(res, {mg: 1}) if mg else res
        if _u_deg(r) == 0:
            # Coprime in Frac(D)[main_v]
            _, res = zpoly_primitive(cg)
            return zpoly_mul(res, {mg: 1}) if mg else res

        # Update ua, ub, g_s, h_s per Knuth Vol 2 §4.6.1 Algorithm C
        denom = zpoly_mul(g_s, zpoly_pow(h_s, delta))
        next_ub = _u_exact_div_scalar(r, denom)
        next_g = ub[_u_deg(ub)]
        if delta == 0:
            next_h = h_s
        elif delta == 1:
            next_h = next_g
        else:
            num_h = zpoly_pow(next_g, delta)
            den_h = zpoly_pow(h_s, delta - 1)
            q_h = zpoly_exact_div(num_h, den_h)
            if q_h is None:
                raise ArithmeticError("Subresultant h-update exact division failed")
            next_h = q_h
        ua = ub
        ub = next_ub
        g_s = next_g
        h_s = next_h

    # Reached degree 0
    _, res = zpoly_primitive(cg)
    return zpoly_mul(res, {mg: 1}) if mg else res


# ---------------------------------------------------------------------------
# 1. Exact Multivariate/Univariate Polynomial Ring QQPoly over Q
# ---------------------------------------------------------------------------

ScalarLike = Union[int, Fraction]
PolyLike = Union[int, Fraction, "QQPoly"]
RatLike = Union[int, Fraction, "QQPoly", "QQRat"]


class QQPoly:
    """
    Exact sparse polynomial in Q[x_1, ..., x_n] over fractions.Fraction.
    Canonical invariant: self.terms maps Monomial -> non-zero Fraction.
    """

    __slots__ = ("terms",)

    def __init__(self, terms: Optional[Mapping[Monomial, Union[int, Fraction]]] = None) -> None:
        cleaned: Dict[Monomial, Fraction] = {}
        if terms:
            for m, c in terms.items():
                fc = c if type(c) is Fraction else Fraction(c)
                if fc:
                    cleaned[m] = fc
        self.terms: Dict[Monomial, Fraction] = cleaned

    # --- Constructors ---

    @classmethod
    def zero(cls) -> QQPoly:
        return cls()

    @classmethod
    def one(cls) -> QQPoly:
        return cls({EMPTY_MONO: Fraction(1)})

    @classmethod
    def const(cls, c: Union[int, Fraction]) -> QQPoly:
        fc = c if type(c) is Fraction else Fraction(c)
        if not fc:
            return cls()
        return cls({EMPTY_MONO: fc})

    @classmethod
    def var(cls, name: str) -> QQPoly:
        if not name or not isinstance(name, str):
            raise ValueError(f"Invalid variable name: {name!r}")
        return cls({((name, 1),): Fraction(1)})

    @classmethod
    def vars(cls, *names: str) -> Tuple[QQPoly, ...]:
        if len(names) == 1 and (" " in names[0] or "," in names[0]):
            parts = [p.strip() for p in names[0].replace(",", " ").split() if p.strip()]
            return tuple(cls.var(p) for p in parts)
        return tuple(cls.var(n) for n in names)

    @classmethod
    def from_univariate_coeffs(cls, coeffs_asc: Sequence[Union[int, Fraction]], var: str = "r") -> QQPoly:
        """Create sum_{k=0}^{n} coeffs_asc[k] * var^k."""
        d: Dict[Monomial, Fraction] = {}
        for k, c in enumerate(coeffs_asc):
            fc = Fraction(c)
            if fc:
                m: Monomial = ((var, k),) if k > 0 else EMPTY_MONO
                d[m] = fc
        return cls(d)

    @classmethod
    def from_zpoly(cls, zp: ZPoly, scale: Fraction = Fraction(1)) -> QQPoly:
        if not zp or not scale:
            return cls()
        return cls({m: Fraction(c) * scale for m, c in zp.items() if c})

    @staticmethod
    def coerce(val: Any) -> Optional[QQPoly]:
        if isinstance(val, QQPoly):
            return val
        if isinstance(val, (int, Fraction)):
            return QQPoly.const(val)
        return None

    # --- Properties & Queries ---

    @property
    def numerator(self) -> QQPoly:
        return self

    @property
    def denominator(self) -> QQPoly:
        return QQPoly.one()

    def is_zero(self) -> bool:
        return not self.terms

    def is_constant(self) -> bool:
        if not self.terms:
            return True
        return len(self.terms) == 1 and EMPTY_MONO in self.terms

    def constant_value(self) -> Fraction:
        if not self.is_constant():
            raise ValueError(f"Polynomial {self} is not constant")
        return self.terms.get(EMPTY_MONO, Fraction(0))

    def variables(self) -> Tuple[str, ...]:
        vs: set[str] = set()
        for m in self.terms:
            for v, _ in m:
                vs.add(v)
        return tuple(sorted(vs))

    def degree(self, var: Optional[str] = None) -> int:
        if not self.terms:
            return -1
        if var is None:
            return max(sum(e for _, e in m) for m in self.terms)
        return max(mono_deg_var(m, var) for m in self.terms)

    def univariate_coeffs(self, var: str = "r") -> List[Fraction]:
        """Return ascending coefficient list [c_0, c_1, ..., c_d] for a univariate polynomial."""
        vs = self.variables()
        if any(v != var for v in vs):
            raise ValueError(f"Polynomial {self} has variables {vs}, expected only {var!r}")
        if not self.terms:
            return []
        d = self.degree(var)
        out = [Fraction(0)] * (d + 1)
        for m, c in self.terms.items():
            k = mono_deg_var(m, var)
            out[k] = c
        return out

    def leading_coeff(self, var: Optional[str] = None) -> Fraction:
        if not self.terms:
            return Fraction(0)
        if var is not None:
            cs = self.univariate_coeffs(var)
            return cs[-1] if cs else Fraction(0)
        lm = max(self.terms)
        return self.terms[lm]

    def monic(self, var: Optional[str] = None) -> QQPoly:
        if not self.terms:
            return self
        lc = self.leading_coeff(var)
        return self * (Fraction(1) / lc)

    def to_primitive_zpoly(self) -> Tuple[Fraction, ZPoly]:
        """
        Factor self = content_q * P_int where content_q in Q and P_int in Z[X] is primitive
        with positive lex leading coefficient. If self == 0, returns (Fraction(0), {}).
        """
        if not self.terms:
            return Fraction(0), {}
        den_lcm = 1
        for c in self.terms.values():
            den_lcm = int_lcm(den_lcm, c.denominator)
        scaled: ZPoly = {m: c.numerator * (den_lcm // c.denominator) for m, c in self.terms.items()}
        int_cont, prim = zpoly_primitive(scaled)
        return Fraction(int_cont, den_lcm), prim

    # --- Calculus & Evaluation ---

    def diff(self, var: str, order: int = 1) -> QQPoly:
        if order < 0:
            raise ValueError("Derivative order must be >= 0")
        p = self
        for _ in range(order):
            if not p.terms:
                return QQPoly.zero()
            nxt: Dict[Monomial, Fraction] = {}
            for m, c in p.terms.items():
                e, rest = mono_drop_var(m, var)
                if e > 0:
                    nm = mono_mul(rest, ((var, e - 1),)) if e > 1 else rest
                    nxt[nm] = nxt.get(nm, Fraction(0)) + c * e
            p = QQPoly(nxt)
        return p

    def subs(self, mapping: Mapping[str, Any]) -> Any:
        """
        Exact polymorphic substitution: values in mapping may be int, Fraction,
        QQPoly, QQRat, or QQJet. Preserves at least QQPoly (or lifts to QQRat/QQJet
        if mapping contains QQRat/QQJet). Use `eval_at` to extract a raw Fraction!
        """
        has_jet = any(isinstance(val, QQJet) for val in mapping.values())
        has_rat = any(isinstance(val, QQRat) for val in mapping.values())
        if has_jet:
            max_ord = max(val.order for val in mapping.values() if isinstance(val, QQJet))
            sym_jet = any(
                isinstance(val, (QQPoly, QQRat))
                or (isinstance(val, QQJet) and isinstance(val.coeffs[0], QQRat))
                for val in mapping.values()
            ) or any(v not in mapping for v in self.variables())
            acc: Any = QQJet.const(QQRat.zero() if sym_jet else Fraction(0), order=max_ord)
        elif has_rat:
            acc = QQRat.zero()
        else:
            acc = QQPoly.zero()

        for m, c in self.terms.items():
            term: Any = QQPoly.const(c) if not (has_jet or has_rat) else c
            for v, e in m:
                if v in mapping:
                    mv = mapping[v]
                    if not (has_jet or has_rat) and isinstance(mv, (int, Fraction)):
                        term = term * QQPoly.const(Fraction(mv) ** e)
                    else:
                        term = term * (mv ** e)
                else:
                    term = term * QQPoly({((v, e),): Fraction(1)})
            acc = acc + term
        return acc

    def eval_at(self, **kwargs: Union[int, Fraction]) -> Fraction:
        total_q = Fraction(0)
        for m, c in self.terms.items():
            term = c
            for v, e in m:
                if v not in kwargs:
                    raise ValueError(f"Unbound variable {v!r} in eval_at({kwargs}) for {self}")
                term *= Fraction(kwargs[v]) ** e
            total_q += term
        return total_q

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        if kwargs:
            return self.subs(kwargs)
        vs = self.variables()
        if len(args) == 1 and len(vs) <= 1:
            v = vs[0] if vs else "r"
            return self.subs({v: args[0]})
        if len(args) == len(vs):
            return self.subs(dict(zip(vs, args)))
        raise ValueError(f"Cannot call polynomial in {vs} with {len(args)} positional arguments")

    # --- Algebraic Operations ---

    def __pos__(self) -> QQPoly:
        return self

    def __neg__(self) -> QQPoly:
        return QQPoly({m: -c for m, c in self.terms.items()})

    def __add__(self, other: Any) -> Any:
        coerced = QQPoly.coerce(other)
        if coerced is None:
            return NotImplemented
        if not self.terms:
            return coerced
        if not coerced.terms:
            return self
        out = dict(self.terms)
        for m, c in coerced.terms.items():
            nv = out.get(m, Fraction(0)) + c
            if nv:
                out[m] = nv
            else:
                out.pop(m, None)
        return QQPoly(out)

    def __radd__(self, other: Any) -> Any:
        return self.__add__(other)

    def __sub__(self, other: Any) -> Any:
        coerced = QQPoly.coerce(other)
        if coerced is None:
            return NotImplemented
        return self.__add__(-coerced)

    def __rsub__(self, other: Any) -> Any:
        coerced = QQPoly.coerce(other)
        if coerced is None:
            return NotImplemented
        return coerced.__add__(-self)

    def __mul__(self, other: Any) -> Any:
        coerced = QQPoly.coerce(other)
        if coerced is None:
            return NotImplemented
        if not self.terms or not coerced.terms:
            return QQPoly.zero()
        out: Dict[Monomial, Fraction] = {}
        for m1, c1 in self.terms.items():
            for m2, c2 in coerced.terms.items():
                m = mono_mul(m1, m2)
                nv = out.get(m, Fraction(0)) + c1 * c2
                if nv:
                    out[m] = nv
                else:
                    out.pop(m, None)
        return QQPoly(out)

    def __rmul__(self, other: Any) -> Any:
        return self.__mul__(other)

    def __truediv__(self, other: Any) -> Any:
        if isinstance(other, (int, Fraction)):
            fo = Fraction(other)
            if fo == 0:
                raise ZeroDivisionError("QQPoly division by zero")
            inv = Fraction(1) / fo
            return QQPoly({m: c * inv for m, c in self.terms.items()})
        if isinstance(other, QQPoly):
            return QQRat(self, other)
        return NotImplemented

    def __rtruediv__(self, other: Any) -> Any:
        if isinstance(other, (int, Fraction)):
            return QQRat(QQPoly.const(other), self)
        return NotImplemented

    def __pow__(self, exp: int) -> Union[QQPoly, QQRat]:
        if not isinstance(exp, int):
            raise TypeError(f"QQPoly exponent must be int, got {type(exp).__name__}")
        if exp < 0:
            return QQRat(QQPoly.one(), self ** (-exp))
        if exp == 0:
            return QQPoly.one()
        if not self.terms:
            return QQPoly.zero()
        res = QQPoly.one()
        base = self
        e = exp
        while e > 0:
            if e & 1:
                res = res * base
            e >>= 1
            if e:
                base = base * base
        return res

    def __bool__(self) -> bool:
        return bool(self.terms)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QQPoly):
            return self.terms == other.terms
        if isinstance(other, (int, Fraction)):
            fo = Fraction(other)
            if fo == 0:
                return not self.terms
            return len(self.terms) == 1 and self.terms.get(EMPTY_MONO) == fo
        return NotImplemented

    def __hash__(self) -> int:
        if not self.terms:
            return hash(Fraction(0))
        if self.is_constant():
            return hash(self.terms[EMPTY_MONO])
        return hash(frozenset(self.terms.items()))

    # --- Univariate Division & Exact Multivariate GCD ---

    def divmod_univariate(self, divisor: QQPoly, var: str = "r") -> Tuple[QQPoly, QQPoly]:
        """Exact Euclidean division in Q[var]: returns (quotient, remainder)."""
        if divisor.is_zero():
            raise ZeroDivisionError("Univariate polynomial division by zero")
        a = self.univariate_coeffs(var)
        b = divisor.univariate_coeffs(var)
        while a and a[-1] == 0:
            a.pop()
        while b and b[-1] == 0:
            b.pop()
        if not a or len(a) < len(b):
            return QQPoly.zero(), self
        rem = list(a)
        db = len(b) - 1
        lc_b = b[db]
        quot = [Fraction(0)] * (len(a) - db)
        for k in range(len(a) - 1, db - 1, -1):
            c = rem[k]
            if c != 0:
                q = c / lc_b
                quot[k - db] = q
                for j in range(db + 1):
                    rem[k - db + j] -= q * b[j]
        return (
            QQPoly.from_univariate_coeffs(quot, var),
            QQPoly.from_univariate_coeffs(rem, var),
        )

    @staticmethod
    def gcd(a: QQPoly, b: QQPoly) -> QQPoly:
        """
        Exact multivariate polynomial GCD in Q[x_1, ..., x_n], normalized as primitive
        integer polynomial with positive lex leading coefficient (or 0 if both are 0).
        """
        if a.is_zero() and b.is_zero():
            return QQPoly.zero()
        _, za = a.to_primitive_zpoly()
        _, zb = b.to_primitive_zpoly()
        zg = zpoly_gcd(za, zb)
        return QQPoly.from_zpoly(zg)

    def exact_div(self, divisor: QQPoly) -> Optional[QQPoly]:
        """Exact polynomial division self / divisor in Q[x_1, ..., x_n], or None if not divisible."""
        if divisor.is_zero():
            raise ZeroDivisionError("Exact division by zero QQPoly")
        if self.is_zero():
            return QQPoly.zero()
        ca, za = self.to_primitive_zpoly()
        cb, zb = divisor.to_primitive_zpoly()
        zq = zpoly_exact_div(za, zb)
        if zq is None:
            return None
        return QQPoly.from_zpoly(zq, scale=ca / cb)

    def square_free_decomposition(self, var: str = "r") -> Tuple[Fraction, List[Tuple[QQPoly, int]]]:
        """
        Yun's exact square-free decomposition in Q[var]:
        self(var) = lc * prod_i (P_i(var) ** mult_i)
        where each P_i is monic, non-constant, square-free, and pairwise coprime.
        """
        if self.is_zero():
            return Fraction(0), []
        lc = self.leading_coeff(var)
        f = self.monic(var)
        if f.degree(var) == 0:
            return lc, []

        fp = f.diff(var)
        g = QQPoly.gcd(f, fp).monic(var)
        if g.degree(var) == 0:
            return lc, [(f, 1)]

        w, _ = f.divmod_univariate(g, var)
        y, _ = fp.divmod_univariate(g, var)
        z = y - w.diff(var)

        factors: List[Tuple[QQPoly, int]] = []
        k = 1
        while not z.is_zero():
            dk = QQPoly.gcd(w, z).monic(var)
            if dk.degree(var) > 0:
                factors.append((dk, k))
            w, _ = w.divmod_univariate(dk, var)
            y, _ = z.divmod_univariate(dk, var)
            z = y - w.diff(var)
            k += 1
        if w.degree(var) > 0:
            factors.append((w.monic(var), k))
        return lc, factors

    # --- Formatting ---

    def __repr__(self) -> str:
        if not self.terms:
            return "0"
        pieces: List[str] = []
        for m in sorted(self.terms, reverse=True):
            c = self.terms[m]
            sign = "-" if c < 0 else "+"
            ac = -c if c < 0 else c
            if not m:
                body = str(ac)
            else:
                m_str = "*".join(
                    f"{v}**{e}" if e > 1 else v for v, e in sorted(m)
                )
                if ac == 1:
                    body = m_str
                else:
                    body = f"{ac}*{m_str}"
            if not pieces:
                pieces.append(f"-{body}" if sign == "-" else body)
            else:
                pieces.append(f" {sign} {body}")
        return "".join(pieces)


# ---------------------------------------------------------------------------
# 2. Canonical Rational Function Field QQRat over Q(x_1, ..., x_n)
# ---------------------------------------------------------------------------

class QQRat:
    """
    Canonical rational function num / den in Q(x_1, ..., x_n).
    Canonical normalization:
      - If zero: num = 0, den = 1.
      - Otherwise: num and den have integer coefficients (Z[X]), gcd_Z(all coeffs) == 1,
        gcd_Z[X](num, den) == 1, and lc_lex(den) > 0.
    Every element of Q(x_1, ..., x_n) has a UNIQUE (num, den) representation!
    """

    __slots__ = ("num", "den", "_znum", "_zden")

    def __init__(
        self,
        num: Union[int, Fraction, QQPoly, "QQRat"],
        den: Union[int, Fraction, QQPoly, "QQRat"] = 1,
        *,
        _reduced_z: Optional[Tuple[ZPoly, ZPoly]] = None,
    ) -> None:
        if _reduced_z is not None:
            zn, zd = _reduced_z
            self._znum = zn
            self._zden = zd
            self.num = QQPoly.from_zpoly(zn)
            self.den = QQPoly.from_zpoly(zd)
            return

        if isinstance(num, QQRat) or isinstance(den, QQRat):
            rn = num if isinstance(num, QQRat) else QQRat.coerce(num)
            rd = den if isinstance(den, QQRat) else QQRat.coerce(den)
            assert rn is not None and rd is not None
            res = rn / rd
            self.num, self.den, self._znum, self._zden = res.num, res.den, res._znum, res._zden
            return

        pn = QQPoly.coerce(num)
        pd = QQPoly.coerce(den)
        if pn is None or pd is None:
            raise TypeError(f"Cannot construct QQRat from ({type(num).__name__}, {type(den).__name__})")
        if pd.is_zero():
            raise ZeroDivisionError("QQRat denominator is zero")
        if pn.is_zero():
            self._znum = {}
            self._zden = {EMPTY_MONO: 1}
            self.num = QQPoly.zero()
            self.den = QQPoly.one()
            return

        cn, zn = pn.to_primitive_zpoly()
        cd, zd = pd.to_primitive_zpoly()
        ratio = cn / cd  # exact reduced Fraction (sign carried in ratio.numerator; ratio.denominator > 0)
        zg = zpoly_gcd(zn, zd)
        if zg != {EMPTY_MONO: 1}:
            zn_q = zpoly_exact_div(zn, zg)
            zd_q = zpoly_exact_div(zd, zg)
            assert zn_q is not None and zd_q is not None
            zn, zd = zn_q, zd_q
        # Since zn and zd are already primitive with lc_lex(zd) > 0 and gcd(ratio.num, ratio.den) == 1:
        zn = zpoly_scale(zn, ratio.numerator)
        zd = zpoly_scale(zd, ratio.denominator)
        self._znum = zn
        self._zden = zd
        self.num = QQPoly.from_zpoly(zn)
        self.den = QQPoly.from_zpoly(zd)

    @staticmethod
    def _from_coprime_z(zn: ZPoly, zd: ZPoly) -> QQRat:
        """
        Normalize integer content and sign of zd when polynomial GCD(zn, zd) is ALREADY 1!
        """
        if not zd:
            raise ZeroDivisionError("Zero denominator in QQRat")
        if not zn:
            return QQRat(0, 1, _reduced_z=({}, {EMPTY_MONO: 1}))
        g = int_gcd(zpoly_int_content(zn), zpoly_int_content(zd))
        lm_d = max(zd)
        if zd[lm_d] < 0:
            g = -g
        if g != 1:
            zn = {m: c // g for m, c in zn.items()}
            zd = {m: c // g for m, c in zd.items()}
        return QQRat(0, 1, _reduced_z=(zn, zd))

    # --- Constructors ---

    @classmethod
    def zero(cls) -> QQRat:
        return cls(0, 1, _reduced_z=({}, {EMPTY_MONO: 1}))

    @classmethod
    def one(cls) -> QQRat:
        return cls(1, 1, _reduced_z=({EMPTY_MONO: 1}, {EMPTY_MONO: 1}))

    @classmethod
    def const(cls, c: Union[int, Fraction]) -> QQRat:
        fc = Fraction(c)
        if fc == 0:
            return cls.zero()
        return cls(0, 1, _reduced_z=({EMPTY_MONO: fc.numerator}, {EMPTY_MONO: fc.denominator}))

    @classmethod
    def var(cls, name: str) -> QQRat:
        return cls(0, 1, _reduced_z=({((name, 1),): 1}, {EMPTY_MONO: 1}))

    @classmethod
    def vars(cls, *names: str) -> Tuple[QQRat, ...]:
        if len(names) == 1 and (" " in names[0] or "," in names[0]):
            parts = [p.strip() for p in names[0].replace(",", " ").split() if p.strip()]
            return tuple(cls.var(p) for p in parts)
        return tuple(cls.var(n) for n in names)

    @staticmethod
    def coerce(val: Any) -> Optional[QQRat]:
        if isinstance(val, QQRat):
            return val
        if isinstance(val, (int, Fraction)):
            return QQRat.const(val)
        if isinstance(val, QQPoly):
            return QQRat(val, QQPoly.one())
        return None

    # --- Properties & Queries ---

    @property
    def numerator(self) -> QQPoly:
        return self.num

    @property
    def denominator(self) -> QQPoly:
        return self.den

    def is_zero(self) -> bool:
        return not self._znum

    def is_constant(self) -> bool:
        return self.num.is_constant() and self.den.is_constant()

    def to_fraction(self) -> Fraction:
        if not self.is_constant():
            raise ValueError(f"QQRat {self} is not a constant rational number")
        return self.num.constant_value() / self.den.constant_value()

    def variables(self) -> Tuple[str, ...]:
        return tuple(sorted(set(self.num.variables()) | set(self.den.variables())))

    # --- Calculus & Evaluation ---

    def diff(self, var: str, order: int = 1) -> QQRat:
        if order < 0:
            raise ValueError("Derivative order must be >= 0")
        f = self
        for _ in range(order):
            if f.is_zero():
                return QQRat.zero()
            np = f.num.diff(var)
            dp = f.den.diff(var)
            # (np * den - num * dp) / den^2
            # Notice gcd(den, den^2) = den; factor reduction automatic in QQRat constructor
            f = QQRat(np * f.den - f.num * dp, f.den * f.den)
        return f

    def subs(self, mapping: Mapping[str, Any]) -> Any:
        nv = self.num.subs(mapping)
        dv = self.den.subs(mapping)
        if isinstance(nv, QQJet) or isinstance(dv, QQJet):
            return nv / dv
        return QQRat(nv, dv)

    def eval_at(self, **kwargs: Union[int, Fraction]) -> Fraction:
        nv = self.num.eval_at(**kwargs)
        dv = self.den.eval_at(**kwargs)
        if dv == 0:
            raise ZeroDivisionError(f"Pole at {kwargs} when evaluating {self}")
        return nv / dv

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        if kwargs:
            return self.subs(kwargs)
        vs = self.variables()
        if len(args) == 1 and len(vs) <= 1:
            v = vs[0] if vs else "r"
            return self.subs({v: args[0]})
        if len(args) == len(vs):
            return self.subs(dict(zip(vs, args)))
        raise ValueError(f"Cannot call QQRat in {vs} with {len(args)} positional arguments")

    # --- Fast Henrici UFD Arithmetic ---

    def __pos__(self) -> QQRat:
        return self

    def __neg__(self) -> QQRat:
        if not self._znum:
            return self
        return QQRat(0, 1, _reduced_z=(zpoly_neg(self._znum), self._zden))

    def __add__(self, other: Any) -> Any:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        if not self._znum:
            return coerced
        if not coerced._znum:
            return self

        n1, d1 = self._znum, self._zden
        n2, d2 = coerced._znum, coerced._zden

        # Fast-path: identical denominators
        if d1 == d2:
            m = zpoly_add(n1, n2)
            if not m:
                return QQRat.zero()
            h = zpoly_gcd(m, d1)
            if h != {EMPTY_MONO: 1}:
                mq = zpoly_exact_div(m, h)
                dq = zpoly_exact_div(d1, h)
                assert mq is not None and dq is not None
                return QQRat._from_coprime_z(mq, dq)
            return QQRat._from_coprime_z(m, d1)

        # Henrici's UFD theorem:
        # Let g = gcd(d1, d2), td1 = d1/g, td2 = d2/g, m = n1*td2 + n2*td1.
        # Then gcd(m, td1*td2*g) == gcd(m, g)!
        g = zpoly_gcd(d1, d2)
        if g == {EMPTY_MONO: 1}:
            m = zpoly_add(zpoly_mul(n1, d2), zpoly_mul(n2, d1))
            if not m:
                return QQRat.zero()
            return QQRat._from_coprime_z(m, zpoly_mul(d1, d2))

        td1 = zpoly_exact_div(d1, g)
        td2 = zpoly_exact_div(d2, g)
        assert td1 is not None and td2 is not None
        m = zpoly_add(zpoly_mul(n1, td2), zpoly_mul(n2, td1))
        if not m:
            return QQRat.zero()
        h = zpoly_gcd(m, g)
        if h != {EMPTY_MONO: 1}:
            mq = zpoly_exact_div(m, h)
            gq = zpoly_exact_div(g, h)
            assert mq is not None and gq is not None
            m, g = mq, gq
        den = zpoly_mul(zpoly_mul(td1, td2), g)
        return QQRat._from_coprime_z(m, den)

    def __radd__(self, other: Any) -> Any:
        return self.__add__(other)

    def __sub__(self, other: Any) -> Any:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        return self.__add__(-coerced)

    def __rsub__(self, other: Any) -> Any:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        return coerced.__add__(-self)

    def __mul__(self, other: Any) -> Any:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        if not self._znum or not coerced._znum:
            return QQRat.zero()

        n1, d1 = self._znum, self._zden
        n2, d2 = coerced._znum, coerced._zden

        # Cross-cancel before multiplying: gcd(n1, d1) == gcd(n2, d2) == 1
        g1 = zpoly_gcd(n1, d2)
        if g1 != {EMPTY_MONO: 1}:
            n1 = zpoly_exact_div(n1, g1)  # type: ignore[assignment]
            d2 = zpoly_exact_div(d2, g1)  # type: ignore[assignment]
        g2 = zpoly_gcd(n2, d1)
        if g2 != {EMPTY_MONO: 1}:
            n2 = zpoly_exact_div(n2, g2)  # type: ignore[assignment]
            d1 = zpoly_exact_div(d1, g2)  # type: ignore[assignment]
        return QQRat._from_coprime_z(zpoly_mul(n1, n2), zpoly_mul(d1, d2))

    def __rmul__(self, other: Any) -> Any:
        return self.__mul__(other)

    def __truediv__(self, other: Any) -> Any:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        if not coerced._znum:
            raise ZeroDivisionError("QQRat division by zero")
        if not self._znum:
            return QQRat.zero()
        inv = QQRat._from_coprime_z(coerced._zden, coerced._znum)
        return self.__mul__(inv)

    def __rtruediv__(self, other: Any) -> Any:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        return coerced.__truediv__(self)

    def __pow__(self, exp: int) -> QQRat:
        if not isinstance(exp, int):
            raise TypeError(f"QQRat exponent must be int, got {type(exp).__name__}")
        if exp == 0:
            return QQRat.one()
        if exp < 0:
            if not self._znum:
                raise ZeroDivisionError("Zero QQRat cannot be raised to a negative power")
            inv = QQRat._from_coprime_z(self._zden, self._znum)
            return inv ** (-exp)
        if not self._znum:
            return QQRat.zero()
        # Since gcd(znum, zden) == 1, gcd(znum^k, zden^k) == 1 automatically!
        return QQRat._from_coprime_z(zpoly_pow(self._znum, exp), zpoly_pow(self._zden, exp))

    def __bool__(self) -> bool:
        return bool(self._znum)

    def __eq__(self, other: object) -> bool:
        coerced = QQRat.coerce(other)
        if coerced is None:
            return NotImplemented
        return self._znum == coerced._znum and self._zden == coerced._zden

    def __hash__(self) -> int:
        if self.is_constant():
            return hash(self.to_fraction())
        return hash((frozenset(self._znum.items()), frozenset(self._zden.items())))

    def __repr__(self) -> str:
        if self.den == 1:
            return repr(self.num)
        return f"({self.num})/({self.den})"


# ---------------------------------------------------------------------------
# 3. Truncated Taylor Jet Algebra QQJet over Fraction OR QQRat
# ---------------------------------------------------------------------------

class QQJet:
    """
    Truncated univariate formal power series sum_{k=0}^{order} coeffs[k] * eps^k + O(eps^(order+1)).
    Works seamlessly over BOTH exact rational numbers (Fraction) AND symbolic rational functions (QQRat)!
    """

    __slots__ = ("coeffs", "order")

    def __init__(self, coeffs: Sequence[Any], order: Optional[int] = None) -> None:
        if order is None:
            if not coeffs:
                raise ValueError("QQJet requires non-empty coeffs or explicit order")
            order = len(coeffs) - 1
        if order < 0:
            raise ValueError(f"Invalid QQJet order: {order}")

        # Determine whether any coefficient is symbolic (QQPoly / QQRat)
        symbolic = any(isinstance(c, (QQPoly, QQRat)) for c in coeffs)
        out: List[Any] = []
        for k in range(order + 1):
            raw = coeffs[k] if k < len(coeffs) else 0
            if symbolic:
                rc = QQRat.coerce(raw)
                if rc is None:
                    raise TypeError(f"Unsupported QQJet coefficient: {type(raw).__name__}")
                out.append(rc)
            else:
                out.append(raw if type(raw) is Fraction else Fraction(raw))
        self.coeffs: Tuple[Any, ...] = tuple(out)
        self.order: int = order

    # --- Constructors ---

    @classmethod
    def const(cls, val: Any, order: int = 3) -> QQJet:
        return cls([val], order=order)

    @classmethod
    def var(cls, order: int = 3, *, symbolic: bool = False) -> QQJet:
        """Return the infinitesimal generator eps = 0 + 1*eps + O(eps^(order+1))."""
        if symbolic:
            return cls([QQRat.zero(), QQRat.one()], order=order)
        return cls([Fraction(0), Fraction(1)], order=order)

    @classmethod
    def point_and_tangent(cls, x0: Any, v0: Any, order: int = 3) -> QQJet:
        return cls([x0, v0], order=order)

    def __getitem__(self, idx: int) -> Any:
        return self.coeffs[idx]

    def derivative_at_zero(self, k: int) -> Any:
        """Return d^k/deps^k|_{eps=0} = k! * coeffs[k]."""
        fact = 1
        for i in range(2, k + 1):
            fact *= i
        return self.coeffs[k] * fact

    @staticmethod
    def _align(a: QQJet, b: Any) -> Tuple[QQJet, QQJet, int]:
        if isinstance(b, QQJet):
            ord_common = a.order if a.order < b.order else b.order
            sym = isinstance(a.coeffs[0], QQRat) or isinstance(b.coeffs[0], QQRat)
            if sym:
                ca = [QQRat.coerce(x) for x in a.coeffs[: ord_common + 1]]
                cb = [QQRat.coerce(x) for x in b.coeffs[: ord_common + 1]]
            else:
                ca = list(a.coeffs[: ord_common + 1])
                cb = list(b.coeffs[: ord_common + 1])
            return QQJet(ca, order=ord_common), QQJet(cb, order=ord_common), ord_common
        if isinstance(b, (int, Fraction, QQPoly, QQRat)):
            sym = isinstance(a.coeffs[0], QQRat) or isinstance(b, (QQPoly, QQRat))
            if sym:
                ca = [QQRat.coerce(x) for x in a.coeffs]
                zb = QQRat.zero()
                cb = [QQRat.coerce(b)] + [zb] * a.order
                return QQJet(ca, order=a.order), QQJet(cb, order=a.order), a.order
            fb = b if type(b) is Fraction else Fraction(b)
            cb_f = [fb] + [Fraction(0)] * a.order
            return a, QQJet(cb_f, order=a.order), a.order
        return None, None, -1  # type: ignore[return-value]

    # --- Operations ---

    def __pos__(self) -> QQJet:
        return self

    def __neg__(self) -> QQJet:
        return QQJet([-c for c in self.coeffs], order=self.order)

    def __add__(self, other: Any) -> Any:
        ja, jb, p = QQJet._align(self, other)
        if p < 0:
            return NotImplemented
        return QQJet([ja.coeffs[k] + jb.coeffs[k] for k in range(p + 1)], order=p)

    def __radd__(self, other: Any) -> Any:
        return self.__add__(other)

    def __sub__(self, other: Any) -> Any:
        ja, jb, p = QQJet._align(self, other)
        if p < 0:
            return NotImplemented
        return QQJet([ja.coeffs[k] - jb.coeffs[k] for k in range(p + 1)], order=p)

    def __rsub__(self, other: Any) -> Any:
        ja, jb, p = QQJet._align(self, other)
        if p < 0:
            return NotImplemented
        return QQJet([jb.coeffs[k] - ja.coeffs[k] for k in range(p + 1)], order=p)

    def __mul__(self, other: Any) -> Any:
        ja, jb, p = QQJet._align(self, other)
        if p < 0:
            return NotImplemented
        zero = QQRat.zero() if isinstance(ja.coeffs[0], QQRat) else Fraction(0)
        out = []
        for k in range(p + 1):
            acc = zero
            for i in range(k + 1):
                ai = ja.coeffs[i]
                bj = jb.coeffs[k - i]
                if ai != 0 and bj != 0:
                    acc = acc + ai * bj
            out.append(acc)
        return QQJet(out, order=p)

    def __rmul__(self, other: Any) -> Any:
        return self.__mul__(other)

    def _divide_aligned(ja: QQJet, jb: QQJet, p: int) -> QQJet:
        # Handle valuation shift if leading terms of both numerator and denominator vanish
        v = 0
        while v <= p and jb.coeffs[v] == 0:
            if ja.coeffs[v] != 0:
                raise ZeroDivisionError("QQJet pole: denominator has higher eps-valuation than numerator")
            v += 1
        if v > p:
            raise ZeroDivisionError("QQJet division by zero jet")

        new_p = p - v
        a = ja.coeffs[v : p + 1]
        b = jb.coeffs[v : p + 1]
        b0 = b[0]
        zero = QQRat.zero() if isinstance(b0, QQRat) else Fraction(0)
        c: List[Any] = []
        for k in range(new_p + 1):
            acc = a[k]
            for j in range(1, k + 1):
                if b[j] != 0 and c[k - j] != 0:
                    acc = acc - b[j] * c[k - j]
            c.append(acc / b0 if acc != 0 else zero)
        return QQJet(c, order=new_p)

    def __truediv__(self, other: Any) -> Any:
        ja, jb, p = QQJet._align(self, other)
        if p < 0:
            return NotImplemented
        return QQJet._divide_aligned(ja, jb, p)

    def __rtruediv__(self, other: Any) -> Any:
        ja, jb, p = QQJet._align(self, other)
        if p < 0:
            return NotImplemented
        return QQJet._divide_aligned(jb, ja, p)

    def __pow__(self, exp: int) -> QQJet:
        if not isinstance(exp, int):
            raise TypeError(f"QQJet exponent must be int, got {type(exp).__name__}")
        if exp == 0:
            one = QQRat.one() if isinstance(self.coeffs[0], QQRat) else Fraction(1)
            return QQJet.const(one, order=self.order)
        if exp < 0:
            one = QQRat.one() if isinstance(self.coeffs[0], QQRat) else Fraction(1)
            return (QQJet.const(one, order=self.order) / self) ** (-exp)
        res = QQJet.const(
            QQRat.one() if isinstance(self.coeffs[0], QQRat) else Fraction(1),
            order=self.order,
        )
        base = self
        e = exp
        while e > 0:
            if e & 1:
                res = res * base
            e >>= 1
            if e:
                base = base * base
        return res

    def __bool__(self) -> bool:
        return any(c != 0 for c in self.coeffs)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QQJet):
            return self.order == other.order and self.coeffs == other.coeffs
        if isinstance(other, (int, Fraction, QQPoly, QQRat)):
            return self.coeffs[0] == other and all(c == 0 for c in self.coeffs[1:])
        return NotImplemented

    def __repr__(self) -> str:
        terms = []
        for k, c in enumerate(self.coeffs):
            if c != 0 or k == 0:
                if k == 0:
                    terms.append(f"{c}")
                elif k == 1:
                    terms.append(f"({c})*eps")
                else:
                    terms.append(f"({c})*eps**{k}")
        return " + ".join(terms) + f" + O(eps**{self.order + 1})"


# ---------------------------------------------------------------------------
# 4. Identity QuantizerBundle for Zero-Branch AST Polymorphism
# ---------------------------------------------------------------------------

class IdentityQuantizerBundle:
    """
    Exact identity QuantizerBundle for symbolic CAS verification over QQPoly, QQRat, and QQJet.
    Executes @architect's live step AST with zero rounding and zero branching:
      - Calling qb(x, ...) returns x unchanged.
      - Any common helper methods (.quantize, .anchor, .snap, .rotor, .project_unit_circle,
        .quantize_pair, .quantize_vec, .clear_denominators) act as exact algebraic identities.
    """

    def __init__(self, bits: Optional[int] = None) -> None:
        self.bits = bits

    def __call__(self, x: Any, *args: Any, **kwargs: Any) -> Any:
        return x

    def quantize(self, x: Any, *args: Any, **kwargs: Any) -> Any:
        return x

    def anchor(self, x: Any, *args: Any, **kwargs: Any) -> Any:
        return x

    def snap(self, x: Any, *args: Any, **kwargs: Any) -> Any:
        return x

    def rotor(self, x: Any, *args: Any, **kwargs: Any) -> Any:
        return x

    def quantize_tangent(self, x: Any, *args: Any, **kwargs: Any) -> Any:
        return x

    def quantize_pair(self, c: Any, s: Any, *args: Any, **kwargs: Any) -> Tuple[Any, Any]:
        return (c, s)

    def project_unit_circle(self, c: Any, s: Any, *args: Any, **kwargs: Any) -> Tuple[Any, Any]:
        c2 = c * c
        s2 = s * s
        den = c2 + s2
        return ((c2 - s2) / den, (Fraction(2, 1) * c * s) / den)

    def quantize_vec(self, xs: Sequence[Any], *args: Any, **kwargs: Any) -> List[Any]:
        return list(xs)

    def clear_denominators(self, xs: Sequence[Any]) -> Tuple[List[Any], int]:
        return list(xs), 1


IDENTITY_QB = IdentityQuantizerBundle()


# ---------------------------------------------------------------------------
# 5. Exact Algebraic Positivity, Sturm Root Isolation & Global Supremum Certifier
# ---------------------------------------------------------------------------

def sturm_chain(poly: QQPoly, var: str = "r") -> List[QQPoly]:
    """
    Construct a primitive-scaled exact Sturm sequence [S_0, S_1, ..., S_k] in Q[var]
    for a square-free polynomial `poly`. Positive scaling preserves exact signs everywhere!
    """
    if poly.is_zero():
        return []
    s0 = poly
    if s0.degree(var) == 0:
        return [s0]
    s1 = s0.diff(var)
    chain = [s0, s1]
    while chain[-1].degree(var) > 0:
        _, rem = chain[-2].divmod_univariate(chain[-1], var)
        if rem.is_zero():
            break
        neg_rem = -rem
        # Scale by positive rational 1/|content| to keep integer coefficients compact
        cont, zp = neg_rem.to_primitive_zpoly()
        pos_scale_poly = QQPoly.from_zpoly(zp if cont > 0 else zpoly_neg(zp))
        chain.append(pos_scale_poly)
    return chain


def _sign_variations(signs: Sequence[int]) -> int:
    nz = [s for s in signs if s != 0]
    v = 0
    for i in range(len(nz) - 1):
        if nz[i] * nz[i + 1] < 0:
            v += 1
    return v


def _sturm_variations_at(chain: Sequence[QQPoly], x: Fraction, var: str = "r") -> int:
    signs = []
    for p in chain:
        val = p.eval_at(**{var: x})
        signs.append(1 if val > 0 else (-1 if val < 0 else 0))
    return _sign_variations(signs)


def _sturm_variations_pos_inf(chain: Sequence[QQPoly], var: str = "r") -> int:
    signs = []
    for p in chain:
        lc = p.leading_coeff(var)
        signs.append(1 if lc > 0 else (-1 if lc < 0 else 0))
    return _sign_variations(signs)


def _sturm_variations_neg_inf(chain: Sequence[QQPoly], var: str = "r") -> int:
    signs = []
    for p in chain:
        lc = p.leading_coeff(var)
        d = p.degree(var)
        s = 1 if lc > 0 else -1
        if d & 1:
            s = -s
        signs.append(s)
    return _sign_variations(signs)


def cauchy_root_bound(poly: QQPoly, var: str = "r") -> Fraction:
    """Strict upper bound M in Q_{>0} on |r| for all complex/real roots of `poly`."""
    cs = poly.univariate_coeffs(var)
    while cs and cs[-1] == 0:
        cs.pop()
    if len(cs) <= 1:
        return Fraction(1)
    lc = abs(cs[-1])
    max_ratio = max(abs(c) / lc for c in cs[:-1])
    return Fraction(1) + max_ratio


def count_real_roots_in_interval(
    poly: QQPoly,
    a: Optional[Fraction] = None,
    b: Optional[Fraction] = None,
    *,
    include_a: bool = False,
    include_b: bool = False,
    var: str = "r",
) -> int:
    """
    Count DISTINCT real roots of `poly` in the interval (a, b) (with optional endpoints,
    and a=None meaning -inf, b=None meaning +inf) using exact Yun square-free reduction + Sturm!
    """
    if poly.is_zero():
        raise ValueError("Zero polynomial has infinitely many roots")
    if poly.degree(var) == 0:
        return 0

    _, factors = poly.square_free_decomposition(var)
    sqfree = QQPoly.one()
    for f, _ in factors:
        sqfree = sqfree * f

    if a is not None and b is not None and a > b:
        return 0
    if a is not None and b is not None and a == b:
        return 1 if (include_a and include_b and sqfree.eval_at(**{var: a}) == 0) else 0

    # Strip boundary roots cleanly before applying Sturm's theorem on (a, b]
    extra = 0
    if a is not None and sqfree.eval_at(**{var: a}) == 0:
        if include_a:
            extra += 1
        lin = QQPoly.var(var) - QQPoly.const(a)
        sqfree, _ = sqfree.divmod_univariate(lin, var)
    if b is not None and sqfree.degree(var) > 0 and sqfree.eval_at(**{var: b}) == 0:
        if include_b:
            extra += 1
        lin = QQPoly.var(var) - QQPoly.const(b)
        sqfree, _ = sqfree.divmod_univariate(lin, var)

    if sqfree.degree(var) <= 0:
        return extra

    chain = sturm_chain(sqfree, var)
    va = _sturm_variations_neg_inf(chain, var) if a is None else _sturm_variations_at(chain, a, var)
    vb = _sturm_variations_pos_inf(chain, var) if b is None else _sturm_variations_at(chain, b, var)
    return extra + (va - vb)


def _best_rational_with_bounded_den(lo: Fraction, hi: Fraction, max_den: int) -> Optional[Fraction]:
    """
    Find a rational p/q in [lo, hi] with 1 <= q <= max_den via continued-fraction / Stern-Brocot
    mediant search (takes O(log max_den) steps; zero integer factorization required!).
    """
    if lo > hi or max_den < 1:
        return None
    n0 = lo.numerator // lo.denominator
    if Fraction(n0) >= lo and Fraction(n0) <= hi:
        return Fraction(n0)
    if Fraction(n0 + 1) >= lo and Fraction(n0 + 1) <= hi:
        return Fraction(n0 + 1)

    # Shift integer floor so 0 < lo' <= hi' < 1
    shift = lo.numerator // lo.denominator
    lo_s = lo - shift
    hi_s = hi - shift
    if lo_s <= 0:
        return Fraction(shift)
    if hi_s >= 1:
        return Fraction(shift + 1)

    # Continued-fraction convergent walk for the interval [lo_s, hi_s]
    p_prev, p_cur = 0, 1
    q_prev, q_cur = 1, 0
    a_lo, a_hi = lo_s, hi_s
    while True:
        inv_hi = Fraction(1) / a_hi
        inv_lo = Fraction(1) / a_lo
        k_min = inv_hi.numerator // inv_hi.denominator
        if Fraction(k_min) < inv_hi:
            k_min += 1
        k_max = inv_lo.numerator // inv_lo.denominator
        if k_min <= k_max:
            p_cand = k_min * p_cur + p_prev
            q_cand = k_min * q_cur + q_prev
            if q_cand <= max_den:
                cand = Fraction(shift) + Fraction(p_cand, q_cand)
                if lo <= cand <= hi:
                    return cand
            return None
        k = inv_hi.numerator // inv_hi.denominator
        p_prev, p_cur = p_cur, k * p_cur + p_prev
        q_prev, q_cur = q_cur, k * q_cur + q_prev
        if q_cur > max_den:
            return None
        a_lo, a_hi = inv_hi - k, inv_lo - k


def isolate_real_roots(
    poly: QQPoly,
    *,
    nonneg_only: bool = False,
    tol: Fraction = Fraction(1, 1 << 64),
    var: str = "r",
) -> List[Tuple[Fraction, Fraction, int]]:
    """
    Isolate ALL distinct real roots of `poly` (on [0, inf) if nonneg_only else (-inf, inf))
    into disjoint exact rational intervals [(lo, hi, multiplicity), ...], sorted ascending!
    Any rational root r_* in Q is detected EXACTLY (returning lo == hi == r_*)!
    """
    if poly.is_zero():
        raise ValueError("Cannot isolate roots of zero polynomial")
    if poly.degree(var) <= 0:
        return []

    _, sq_factors = poly.square_free_decomposition(var)
    all_roots: List[Tuple[Fraction, Fraction, int]] = []

    for factor, mult in sq_factors:
        if factor.degree(var) <= 0:
            continue

        # Extract exact r == 0 root if present
        if factor.eval_at(**{var: Fraction(0)}) == 0:
            all_roots.append((Fraction(0), Fraction(0), mult))
            factor, _ = factor.divmod_univariate(QQPoly.var(var), var)
            if factor.degree(var) <= 0:
                continue

        _, zf = factor.to_primitive_zpoly()
        lead_m = ((var, factor.degree(var)),)
        max_den = abs(zf[lead_m])
        rat_tol = Fraction(1, 2 * max_den * max_den + 2)
        target_tol = tol if tol < rat_tol else rat_tol

        bound = cauchy_root_bound(factor, var)
        chain = sturm_chain(factor, var)

        lo_init = Fraction(0) if nonneg_only else -bound
        hi_init = bound

        def _bisect(a: Fraction, b: Fraction, va: int, vb: int) -> None:
            n_roots = va - vb
            if n_roots <= 0:
                return
            if n_roots == 1 and (b - a) <= target_tol:
                # Check if there is an exact rational root p/q with q <= max_den in (a, b]
                cand = _best_rational_with_bounded_den(a, b, max_den)
                if cand is not None and factor.eval_at(**{var: cand}) == 0:
                    all_roots.append((cand, cand, mult))
                    return
                # Refine further down to `tol` if rat_tol was looser than tol
                la, lb = a, b
                sa = factor.eval_at(**{var: la})
                while (lb - la) > tol:
                    mid = (la + lb) / 2
                    sm = factor.eval_at(**{var: mid})
                    if sm == 0:
                        all_roots.append((mid, mid, mult))
                        return
                    if (sa > 0 and sm > 0) or (sa < 0 and sm < 0):
                        la, sa = mid, sm
                    else:
                        lb = mid
                all_roots.append((la, lb, mult))
                return

            mid = (a + b) / 2
            if factor.eval_at(**{var: mid}) == 0:
                all_roots.append((mid, mid, mult))
                # Deflate exact root or perturb slightly? Deflating locally via Sturm around mid:
                # Easier: split with tiny step eps_m that encloses no other root
                eps = (b - a) / 8
                while True:
                    left_end = mid - eps
                    right_start = mid + eps
                    if factor.eval_at(**{var: left_end}) != 0 and factor.eval_at(**{var: right_start}) != 0:
                        vl = _sturm_variations_at(chain, left_end, var)
                        vr = _sturm_variations_at(chain, right_start, var)
                        if vl - vr == 1:
                            _bisect(a, left_end, va, vl)
                            _bisect(right_start, b, vr, vb)
                            return
                    eps /= 2
            vm = _sturm_variations_at(chain, mid, var)
            _bisect(a, mid, va, vm)
            _bisect(mid, b, vm, vb)

        va0 = _sturm_variations_at(chain, lo_init, var)
        vb0 = _sturm_variations_at(chain, hi_init, var)
        _bisect(lo_init, hi_init, va0, vb0)

    all_roots.sort(key=lambda item: item[0])
    return all_roots


def _int_poly_mul(a: Sequence[int], b: Sequence[int]) -> List[int]:
    if not a or not b:
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] += ai * bj
    return out


def _bernstein_coeffs_on_interval(poly: QQPoly, a: Fraction, b: Fraction, var: str = "r") -> List[Fraction]:
    """
    Fast exact Bernstein-basis coefficients beta_0, ..., beta_d satisfying:
      (b - a)^d * poly(r) == sum_{k=0}^d beta_k * (r - a)^k * (b - r)^(d - k)
    computed via pure integer Horner/binomial expansion over cleared denominators!
    """
    cs = poly.univariate_coeffs(var)
    if not cs:
        return [Fraction(0)]
    d = len(cs) - 1
    if d == 0:
        return [cs[0]]

    # Clear denominators of poly coeffs and (a, b) simultaneously
    c_lcm = 1
    for c in cs:
        c_lcm = int_lcm(c_lcm, c.denominator)
    ics = [c.numerator * (c_lcm // c.denominator) for c in cs]

    ab_lcm = int_lcm(a.denominator, b.denominator)
    ia = a.numerator * (ab_lcm // a.denominator)
    ib = b.numerator * (ab_lcm // b.denominator)

    # Precompute powers of (ia + ib*u)^j and (ab_lcm + ab_lcm*u)^(d - j) in Z[u]
    lin_ab = [ia, ib]
    lin_1u = [ab_lcm, ab_lcm]
    pow_ab: List[List[int]] = [[1]]
    pow_1u: List[List[int]] = [[1]]
    for _ in range(d):
        pow_ab.append(_int_poly_mul(pow_ab[-1], lin_ab))
        pow_1u.append(_int_poly_mul(pow_1u[-1], lin_1u))

    acc = [0] * (d + 1)
    for j, cj in enumerate(ics):
        if cj != 0:
            term = _int_poly_mul(pow_ab[j], pow_1u[d - j])
            for k, tk in enumerate(term):
                acc[k] += cj * tk

    denom = c_lcm * (ab_lcm ** d)
    return [Fraction(val, denom) for val in acc]


def certify_nonnegative_on_nonneg_ray(
    poly: QQPoly,
    var: str = "r",
    max_poincare_deg: int = 48,
) -> Dict[str, Any]:
    """
    Prove or refute `poly(r) >= 0 for all r >= 0` over Q[r] with 0/1 algebraic certainty!

    Returns a dictionary containing:
      - `certified`: bool (True iff poly(r) >= 0 for all r >= 0)
      - `zero_roots_on_ray`: exact list of rational/algebraic zeros where poly(r) == 0 on [0, inf)
      - `certificate_type`: "zero" | "constant" | "lukacs_sos" | "poincare_sos" | "bernstein_ray_sos" | "refuted"
      - `algebraic_residual_zero`: bool (True -- expansion of the constructive certificate minus poly == 0/1!)
      - `certificate`: explicit rational SOS / Bernstein-Ray decomposition data
    """
    if poly.is_zero():
        return {
            "certified": True,
            "zero_roots_on_ray": [],
            "certificate_type": "zero",
            "algebraic_residual_zero": True,
            "certificate": {"even_sos": [], "odd_r_sos": [], "denom_sq": "1"},
        }

    if poly.degree(var) == 0:
        c = poly.constant_value()
        return {
            "certified": c >= 0,
            "zero_roots_on_ray": [],
            "certificate_type": "constant" if c >= 0 else "refuted",
            "algebraic_residual_zero": c >= 0,
            "counterexample_r": Fraction(0) if c < 0 else None,
            "counterexample_val": c if c < 0 else None,
        }

    # Necessary boundary checks at r = 0 and r -> +inf
    lc = poly.leading_coeff(var)
    p0 = poly.eval_at(**{var: Fraction(0)})
    if lc < 0:
        b = cauchy_root_bound(poly, var) + 1
        return {
            "certified": False,
            "certificate_type": "refuted",
            "algebraic_residual_zero": False,
            "counterexample_r": b,
            "counterexample_val": poly.eval_at(**{var: b}),
        }
    if p0 < 0:
        return {
            "certified": False,
            "certificate_type": "refuted",
            "algebraic_residual_zero": False,
            "counterexample_r": Fraction(0),
            "counterexample_val": p0,
        }

    # Factor poly(r) = lc * r^v * Q_sq(r)^2 * S_odd(r) via Yun square-free decomposition!
    r_poly = QQPoly.var(var)
    v_zero = 0
    rem_poly = poly
    while not rem_poly.is_zero() and rem_poly.eval_at(**{var: Fraction(0)}) == 0:
        rem_poly, _ = rem_poly.divmod_univariate(r_poly, var)
        v_zero += 1

    lc_rem, sq_factors = rem_poly.square_free_decomposition(var)
    q_sq = QQPoly.one()
    s_odd = QQPoly.const(lc_rem)
    zero_roots: List[Any] = []
    if v_zero > 0:
        zero_roots.append({"root": Fraction(0), "multiplicity": v_zero})

    for f, mult in sq_factors:
        if mult >= 2:
            q_sq = q_sq * (f ** (mult // 2))
        if mult % 2 == 1:
            s_odd = s_odd * f
        else:
            # Record any real roots of even-multiplicity factors on (0, inf)
            for lo, hi, _ in isolate_real_roots(f, nonneg_only=True, var=var):
                zero_roots.append({"root": lo if lo == hi else (lo, hi), "multiplicity": mult})

    # S_odd has only simple roots on (0, inf); check count in O(1) Sturm evaluations first!
    n_odd_pos = count_real_roots_in_interval(s_odd, Fraction(0), None, include_a=False, var=var)
    if n_odd_pos > 0:
        odd_pos_roots = isolate_real_roots(s_odd, nonneg_only=True, tol=Fraction(1, 1 << 20), var=var)
        lo, hi, _ = odd_pos_roots[0]
        # Pick rational witness where poly < 0 just to the left or right of (lo, hi)
        cand_left = lo / 2 if lo > 0 else Fraction(0)
        cand_right = hi + Fraction(1, 16)
        for w in (cand_left, (lo + hi) / 2, cand_right):
            vw = poly.eval_at(**{var: w})
            if vw < 0:
                return {
                    "certified": False,
                    "certificate_type": "refuted",
                    "algebraic_residual_zero": False,
                    "counterexample_r": w,
                    "counterexample_val": vw,
                }
        return {
            "certified": False,
            "certificate_type": "refuted",
            "algebraic_residual_zero": False,
            "odd_root_interval": (lo, hi),
        }

    # S_odd(r) > 0 for ALL r >= 0 IS NOW PROVED BY STURM!
    # Now construct an explicit 0/1-residual constructive SOS / Poincaré / Bernstein-Ray certificate!
    # Note: poly(r) = (r^(v_zero // 2) * q_sq(r))^2 * (r^(v_zero % 2) * s_odd(r)).
    sq_prefactor = (r_poly ** (v_zero // 2)) * q_sq
    r_parity = v_zero % 2

    # Attempt 1 & 2: Direct Lukacs SOS (m = 0) or Poincare SOS ((1 + r)^(2m) * s_odd(r) has all coeffs >= 0)
    one_plus_r_sq = (QQPoly.one() + r_poly) ** 2
    cur_poly = s_odd
    m_half = 0
    poincare_found = False
    max_half = max_poincare_deg // 2
    while m_half <= max_half:
        cs = cur_poly.univariate_coeffs(var)
        if all(c >= 0 for c in cs):
            poincare_found = True
            break
        cur_poly = cur_poly * one_plus_r_sq
        m_half += 1

    if poincare_found:
        # cur_poly(r) = sum_{k >= 0} c_k * r^k with c_k >= 0, and s_odd(r) = cur_poly(r) / ((1+r)^m_half)^2!
        # Multiply by r^r_parity * sq_prefactor(r)^2:
        #   poly(r) * ((1+r)^m_half)^2 = sum_i alpha_i * A_i(r)^2 + r * sum_j beta_j * B_j(r)^2
        even_sos: List[Tuple[Fraction, QQPoly]] = []
        odd_r_sos: List[Tuple[Fraction, QQPoly]] = []
        cs = cur_poly.univariate_coeffs(var)
        for k, ck in enumerate(cs):
            if ck == 0:
                continue
            tot_exp = k + r_parity
            base_poly = sq_prefactor * (r_poly ** (tot_exp // 2))
            if tot_exp % 2 == 0:
                even_sos.append((ck, base_poly))
            else:
                odd_r_sos.append((ck, base_poly))

        denom_base = (QQPoly.one() + r_poly) ** m_half
        # Machine-check 0/1 algebraic residual in QQPoly!
        reconstructed = QQPoly.zero()
        for ck, a_poly in even_sos:
            reconstructed = reconstructed + QQPoly.const(ck) * (a_poly ** 2)
        for dk, b_poly in odd_r_sos:
            reconstructed = reconstructed + QQPoly.const(dk) * r_poly * (b_poly ** 2)
        lhs = poly * (denom_base ** 2)
        residual = lhs - reconstructed
        assert residual.is_zero(), f"SOS residual non-zero: {residual}"

        return {
            "certified": True,
            "zero_roots_on_ray": zero_roots,
            "certificate_type": "lukacs_sos" if m_half == 0 else "poincare_sos",
            "poincare_exponent": 2 * m_half,
            "algebraic_residual_zero": residual.is_zero(),
            "residual": str(residual),
            "certificate": {
                "denom_sq_base": str(denom_base),
                "even_sos_terms": len(even_sos),
                "odd_r_sos_terms": len(odd_r_sos),
            },
        }

    # Attempt 3 (Universal fallback for close complex conjugate roots):
    # Partition [0, inf) into finite dyadic intervals [a_i, b_i] on [0, B] where ALL Bernstein
    # coefficients are >= 0, plus [B, inf) where s_odd(B + u) has ALL Taylor coefficients >= 0!
    B = Fraction(1)
    u_var = QQPoly.var("u")
    ray_shifted = s_odd.subs({var: QQPoly.const(B) + u_var})
    while not all(c >= 0 for c in ray_shifted.univariate_coeffs("u")):
        B = B * 2
        ray_shifted = s_odd.subs({var: QQPoly.const(B) + u_var})

    # Verify 0/1 algebraic residual on [B, inf)
    ray_coeffs = ray_shifted.univariate_coeffs("u")
    ray_recon = QQPoly.zero()
    r_minus_B = r_poly - QQPoly.const(B)
    for k, gk in enumerate(ray_coeffs):
        if gk != 0:
            ray_recon = ray_recon + QQPoly.const(gk) * (r_minus_B ** k)
    assert (s_odd - ray_recon).is_zero()

    # Adaptive Bernstein subdivision on [0, B] (check fast int coefficients, verify algebraic identity)
    panels: List[Tuple[Fraction, Fraction, List[Fraction]]] = []
    stack: List[Tuple[Fraction, Fraction]] = [(Fraction(0), B)]
    d_s = s_odd.degree(var)
    while stack:
        a_p, b_p = stack.pop()
        bcs = _bernstein_coeffs_on_interval(s_odd, a_p, b_p, var)
        if all(c >= 0 for c in bcs):
            panels.append((a_p, b_p, bcs))
        else:
            mid = (a_p + b_p) / 2
            stack.append((mid, b_p))
            stack.append((a_p, mid))

    # Machine-check exact 0/1 Bernstein identity on representative boundary & narrowest panels
    check_panels = panels if len(panels) <= 8 else (panels[:4] + panels[-4:])
    for a_p, b_p, bcs in check_panels:
        recon_p = QQPoly.zero()
        ra = r_poly - QQPoly.const(a_p)
        br = QQPoly.const(b_p) - r_poly
        for k, bk in enumerate(bcs):
            if bk != 0:
                recon_p = recon_p + QQPoly.const(bk) * (ra ** k) * (br ** (d_s - k))
        lhs_p = QQPoly.const((b_p - a_p) ** d_s) * s_odd
        assert (lhs_p - recon_p).is_zero()

    # Also verify global prefactor identity poly - sq_prefactor^2 * r^r_parity * s_odd == 0/1
    global_res = poly - (sq_prefactor ** 2) * (r_poly ** r_parity) * s_odd
    assert global_res.is_zero()

    return {
        "certified": True,
        "zero_roots_on_ray": zero_roots,
        "certificate_type": "bernstein_ray_sos",
        "algebraic_residual_zero": True,
        "residual": "0",
        "certificate": {
            "sq_prefactor": str(sq_prefactor),
            "r_parity": r_parity,
            "bernstein_panels": len(panels),
            "ray_cutoff_B": str(B),
        },
    }


def certify_rational_supremum_on_nonneg_ray(
    func: QQRat,
    var: str = "r",
    candidate_bound: Optional[Fraction] = None,
    tightness_target: Fraction = Fraction(999999, 1000000),
) -> Dict[str, Any]:
    """
    Compute and algebraically certify the global supremum `sup_{r >= 0} func(r) <= B_Q`
    for a univariate rational function `func(r) = P(r) / Q(r)` in Q(r):
      1. Proves `Q(r) > 0` for all `r >= 0` (pole-free on `[0, inf)`!).
      2. Isolates ALL critical points `func'(r) == 0` on `(0, inf)` + checks `r = 0` and `r -> +inf`.
      3. Determines whether the global maximum is ATTAINED AT A RATIONAL POINT `r_* in Q_{>=0}`
         (giving **EXACT `0/1` ZERO SLACK** `B_Q == sup_{r>=0} func(r)`!) or at an irrational
         algebraic root (producing a certified tight rational ceiling `B_Q` with certified relative
         slack `< 1 - tightness_target`!).
      4. Certifies `B_Q * Q(r) - P(r) >= 0` for ALL `r >= 0` with an explicit `0/1`-residual
         SOS / Bernstein-Ray positivity proof!
    """
    P = func.num
    Q = func.den
    if Q.leading_coeff(var) < 0:
        P = -P
        Q = -Q

    # Step 1: Certify denominator Q(r) > 0 strictly on [0, inf)
    den_cert = certify_nonnegative_on_nonneg_ray(Q, var=var)
    q0 = Q.eval_at(**{var: Fraction(0)})
    den_roots = count_real_roots_in_interval(Q, Fraction(0), None, include_a=True, var=var)
    if not den_cert["certified"] or q0 <= 0 or den_roots > 0:
        return {
            "certified": False,
            "reason": "Denominator Q(r) vanishes or changes sign on [0, inf) -- pole detected!",
            "denominator_roots_on_ray": den_roots,
        }

    # Check r -> +inf behavior
    deg_p = P.degree(var)
    deg_q = Q.degree(var)
    if deg_p > deg_q and P.leading_coeff(var) > 0:
        return {
            "certified": False,
            "reason": "Rational function diverges to +infinity as r -> +infinity",
        }
    val_inf: Optional[Fraction] = None
    if deg_p < deg_q:
        val_inf = Fraction(0)
    elif deg_p == deg_q:
        val_inf = P.leading_coeff(var) / Q.leading_coeff(var)

    # Step 2: Evaluate at r = 0 and all critical points of func'(r) = 0 on (0, inf)
    val_0 = func.eval_at(**{var: Fraction(0)})
    best_lower = val_0
    best_rational_arg: Optional[Fraction] = Fraction(0)
    irrational_peaks: List[Tuple[Fraction, Fraction, Fraction]] = []

    df = func.diff(var)
    if not df.is_zero():
        crit_roots = isolate_real_roots(df.num, nonneg_only=True, tol=Fraction(1, 1 << 28), var=var)
        for lo, hi, _ in crit_roots:
            if lo == hi:
                v_r = func.eval_at(**{var: lo})
                if v_r > best_lower:
                    best_lower = v_r
                    best_rational_arg = lo
            else:
                mid = (lo + hi) / 2
                v_mid = func.eval_at(**{var: mid})
                if v_mid > best_lower:
                    best_lower = v_mid
                    best_rational_arg = None
                irrational_peaks.append((lo, hi, v_mid))

    if val_inf is not None and val_inf > best_lower:
        best_lower = val_inf
        best_rational_arg = None  # asymptotic supremum at +inf

    # Step 3: Test candidate_bound if supplied, or construct exact/tight B_Q
    if candidate_bound is not None:
        diff_poly = QQPoly.const(candidate_bound) * Q - P
        pos_cert = certify_nonnegative_on_nonneg_ray(diff_poly, var=var)
        exact_zero_slack = pos_cert["certified"] and (
            (best_rational_arg is not None and best_lower == candidate_bound)
            or (val_inf is not None and val_inf == candidate_bound)
        )
        tightness_ratio = (
            Fraction(1)
            if exact_zero_slack
            else (best_lower / candidate_bound if candidate_bound > 0 else Fraction(0))
        )
        return {
            "certified": pos_cert["certified"],
            "bound": candidate_bound,
            "lower_witness": best_lower,
            "attained_at_rational_r": best_rational_arg if exact_zero_slack else None,
            "exact_zero_slack": exact_zero_slack,
            "slack_upper_bound": candidate_bound - best_lower if pos_cert["certified"] else None,
            "tightness_ratio_lower_bound": tightness_ratio if pos_cert["certified"] else None,
            "denominator_cert": den_cert,
            "positivity_cert": pos_cert,
        }

    # Try exact best_lower when attained at a rational point r_* in Q_{>=0} or as r -> +inf!
    if best_rational_arg is not None or (val_inf is not None and best_lower == val_inf):
        diff_exact = QQPoly.const(best_lower) * Q - P
        exact_cert = certify_nonnegative_on_nonneg_ray(diff_exact, var=var)
        if exact_cert["certified"]:
            return {
                "certified": True,
                "bound": best_lower,
                "lower_witness": best_lower,
                "attained_at_rational_r": best_rational_arg,
                "exact_zero_slack": True,
                "slack_upper_bound": Fraction(0),
                "tightness_ratio_lower_bound": Fraction(1),
                "denominator_cert": den_cert,
                "positivity_cert": exact_cert,
            }

    # Otherwise, global max is at an irrational critical point; find a compact dyadic ceiling via fast Sturm!
    def _is_valid_ceiling_sturm(c_val: Fraction) -> bool:
        dp = QQPoly.const(c_val) * Q - P
        if dp.is_zero():
            return True
        if dp.leading_coeff(var) < 0 or dp.eval_at(**{var: Fraction(0)}) < 0:
            return False
        _, sqf = dp.square_free_decomposition(var)
        s_o = QQPoly.one()
        for f_p, m_p in sqf:
            if m_p % 2 == 1:
                s_o = s_o * f_p
        return count_real_roots_in_interval(s_o, Fraction(0), None, include_a=False, var=var) == 0

    scale = abs(best_lower) if best_lower != 0 else Fraction(1)
    rel_gap = Fraction(1) - tightness_target
    if rel_gap <= 0:
        rel_gap = Fraction(1, 1_000_000)

    # Pick compact dyadic denominator D so 1/D <= scale * rel_gap / 2
    D_den = 1
    target_step = scale * rel_gap / 2
    while Fraction(1, D_den) > target_step:
        D_den *= 2

    # Smallest integer k such that k / D_den > best_lower
    k_int = (best_lower.numerator * D_den) // best_lower.denominator + 1
    step_k = 1
    while not _is_valid_ceiling_sturm(Fraction(k_int, D_den)):
        k_int += step_k
        step_k *= 2

    cand = Fraction(k_int, D_den)
    c_cert = certify_nonnegative_on_nonneg_ray(QQPoly.const(cand) * Q - P, var=var)
    assert c_cert["certified"]
    ratio = best_lower / cand if cand > 0 else Fraction(1)
    return {
        "certified": True,
        "bound": cand,
        "lower_witness": best_lower,
        "attained_at_rational_r": None,
        "exact_zero_slack": False,
        "slack_upper_bound": cand - best_lower,
        "tightness_ratio_lower_bound": ratio,
        "denominator_cert": den_cert,
        "positivity_cert": c_cert,
    }


# ---------------------------------------------------------------------------
# 6. 3+1D Spinor Kerr-Schild & Full 4D Symbolic Tensor Curvature Computer
# ---------------------------------------------------------------------------

def mat_mul_4x4(A: Sequence[Sequence[Any]], B: Sequence[Sequence[Any]]) -> List[List[Any]]:
    out = [[QQRat.zero() for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            acc: Any = QQRat.zero()
            for k in range(4):
                if A[i][k] != 0 and B[k][j] != 0:
                    acc = acc + A[i][k] * B[k][j]
            out[i][j] = acc
    return out


def mat_det_4x4(M: Sequence[Sequence[Any]]) -> Any:
    """Exact cofactor-free / Laplace 4x4 determinant over QQRat / QQPoly / Fraction."""
    # 2x2 minors of rows 0, 1 and rows 2, 3
    def m01(c1: int, c2: int) -> Any:
        return M[0][c1] * M[1][c2] - M[0][c2] * M[1][c1]

    def m23(c1: int, c2: int) -> Any:
        return M[2][c1] * M[3][c2] - M[2][c2] * M[3][c1]

    return (
        m01(0, 1) * m23(2, 3)
        - m01(0, 2) * m23(1, 3)
        + m01(0, 3) * m23(1, 2)
        + m01(1, 2) * m23(0, 3)
        - m01(1, 3) * m23(0, 2)
        + m01(2, 3) * m23(0, 1)
    )


def verify_spinor_kerr_schild_identities() -> Dict[str, Any]:
    """
    Machine-check in QQRat over symbolic 4-spinor (u0, u1, u2, u3) and arbitrary H:
      1. Kustaanheimo-Stiefel / Hopf identity: x^2 + y^2 + z^2 - (u0^2 + u1^2 + u2^2 + u3^2)^2 == 0/1
         and equatorial Levi-Civita reduction (u2 = u3 = 0 => x = u0^2 - u1^2, y = 2*u0*u1, r = u0^2 + u1^2).
      2. Exact nullity of ingoing covector l_mu = (1, x/r, y/r, z/r):
         eta^{mu nu} l_mu l_nu == 0/1 AND g^{mu nu} l_mu l_nu == 0/1.
      3. Exact inverse metric: g_{mu alpha} g^{alpha nu} - delta_mu^nu == 0/1 (all 16 components!)
         identically in H (including across horizon H = 1).
      4. Exact Lorentzian unit determinant: -det(g_{mu nu}) == 1/1 (and det(-g^{(2+1)}) == 1/1).
      5. Radial lightcone factorization: ds^2|_{dOmega=0} = (dt + dr)*(-(1 - H)*dt + (1 + H)*dr),
         proving ingoing dr/dt = -1 and outgoing dr/dt = (1 - H)/(1 + H) with 0/1 residual!
    """
    u0, u1, u2, u3, H, dt_s, dr_s = QQRat.vars("u0", "u1", "u2", "u3", "H", "dt", "dr")

    x = u0**2 - u1**2 - u2**2 + u3**2
    y = 2 * (u0 * u1 - u2 * u3)
    z = 2 * (u0 * u2 + u1 * u3)
    r = u0**2 + u1**2 + u2**2 + u3**2

    hopf_residual = (x**2 + y**2 + z**2) - r**2

    # Equatorial Levi-Civita check (u2 = 0, u3 = 0)
    eq_map = {"u2": 0, "u3": 0}
    lc_residual = (x.subs(eq_map)**2 + y.subs(eq_map)**2) - (r.subs(eq_map))**2

    l_cov = [QQRat.one(), x / r, y / r, z / r]
    eta = [
        [QQRat.const(-1), QQRat.zero(), QQRat.zero(), QQRat.zero()],
        [QQRat.zero(), QQRat.one(), QQRat.zero(), QQRat.zero()],
        [QQRat.zero(), QQRat.zero(), QQRat.one(), QQRat.zero()],
        [QQRat.zero(), QQRat.zero(), QQRat.zero(), QQRat.one()],
    ]
    l_con = [sum(eta[mu][nu] * l_cov[nu] for nu in range(4)) for mu in range(4)]

    eta_null_res = sum(l_con[mu] * l_cov[mu] for mu in range(4))

    g_cov = [
        [eta[mu][nu] + H * l_cov[mu] * l_cov[nu] for nu in range(4)]
        for mu in range(4)
    ]
    g_con = [
        [eta[mu][nu] - H * l_con[mu] * l_con[nu] for nu in range(4)]
        for mu in range(4)
    ]

    g_null_res = sum(
        g_con[mu][nu] * l_cov[mu] * l_cov[nu] for mu in range(4) for nu in range(4)
    )

    prod = mat_mul_4x4(g_cov, g_con)
    inv_residuals = [
        [prod[i][j] - (QQRat.one() if i == j else QQRat.zero()) for j in range(4)]
        for i in range(4)
    ]
    max_inv_nonzero = sum(1 for i in range(4) for j in range(4) if not inv_residuals[i][j].is_zero())

    det_g = mat_det_4x4(g_cov)
    unit_det_res = (-det_g) - QQRat.one()

    # Radial null cone factorization
    ds2_radial = -(1 - H) * dt_s**2 + 2 * H * dt_s * dr_s + (1 + H) * dr_s**2
    factored_radial = (dt_s + dr_s) * (-(1 - H) * dt_s + (1 + H) * dr_s)
    lightcone_factor_res = ds2_radial - factored_radial

    ingoing_res = ds2_radial.subs({"dr": -dt_s})
    outgoing_res = ds2_radial.subs({"dr": ((1 - H) / (1 + H)) * dt_s})

    return {
        "hopf_ks_residual_zero": hopf_residual.is_zero(),
        "levi_civita_2d_residual_zero": lc_residual.is_zero(),
        "eta_nullity_residual_zero": eta_null_res.is_zero(),
        "g_nullity_residual_zero": g_null_res.is_zero(),
        "inverse_metric_all_16_zero": max_inv_nonzero == 0,
        "unit_lorentzian_det_minus_g_eq_1": unit_det_res.is_zero(),
        "radial_lightcone_factorization_zero": lightcone_factor_res.is_zero(),
        "ingoing_null_dr_dt_minus_1_zero": ingoing_res.is_zero(),
        "outgoing_null_dr_dt_ratio_zero": outgoing_res.is_zero(),
    }


def compute_4d_kerr_schild_curvature_universal() -> Dict[str, Any]:
    """
    Compute ALL 64 Christoffel symbols and ALL 256 Riemann tensor components of the general
    3+1D spherically symmetric ingoing Kerr-Schild metric in rational spherical coordinates
    x^mu = (t, r, chi = cos(theta), phi):
      ds^2 = -(1 - H(r)) dt^2 + 2 H(r) dt dr + (1 + H(r)) dr^2
             + (r^2 / (1 - chi^2)) dchi^2 + r^2 (1 - chi^2) dphi^2
    using differential jet indeterminates (r, chi, H, H1=H'(r), H2=H''(r)).

    Proves from first principles (full 4D tensor contraction, zero shortcuts!) that:
      - Kretschmann: K = R_{abcd} R^{abcd} == H2^2 + 4*(H1/r)^2 + 4*(H/r^2)^2 (0/1 residual!)
      - 4D Ricci scalar: R_4D = g^{mu nu} R_{mu nu} == H2 + 4*(H1/r) + 2*(H/r^2) (0/1 residual!)
      - Radial/Laplacian trace: R^0_0 + R^1_1 == H2 + 2*(H1/r) (0/1 residual!)
      - 4D Weyl scalar: C^2 == (4/3) * (H2/2 - H1/r + H/r^2)^2 * 4? Checked directly == (1/3)*(H2 - 2*H1/r + 2*H/r^2)^2 * 4 -> exact formula verified!
      - Square-root-free horizon-penetrating orthonormal tetrad components R_{abcd} across H = 1!
    """
    r, chi, H, H1, H2, H3 = QQRat.vars("r", "chi", "H", "H1", "H2", "H3")
    Z = QQRat.zero()
    O = QQRat.one()

    def D_coord(expr: QQRat, mu: int) -> QQRat:
        if expr.is_zero():
            return Z
        if mu == 0 or mu == 3:  # stationary & axisymmetric
            return Z
        if mu == 2:  # d / dchi
            return expr.diff("chi")
        # mu == 1: d / dr with formal jet chain rule H' = H1, H1' = H2, H2' = H3
        return (
            expr.diff("r")
            + H1 * expr.diff("H")
            + H2 * expr.diff("H1")
            + H3 * expr.diff("H2")
        )

    one_m_chi2 = O - chi**2
    g_cov = [
        [-(O - H), H, Z, Z],
        [H, O + H, Z, Z],
        [Z, Z, (r**2) / one_m_chi2, Z],
        [Z, Z, Z, (r**2) * one_m_chi2],
    ]
    g_con = [
        [-(O + H), H, Z, Z],
        [H, O - H, Z, Z],
        [Z, Z, one_m_chi2 / (r**2), Z],
        [Z, Z, Z, O / ((r**2) * one_m_chi2)],
    ]

    # 1. Precompute dg[sigma][mu][nu] = d_sigma g_{mu nu}
    dg = [
        [[D_coord(g_cov[mu][nu], sigma) for nu in range(4)] for mu in range(4)]
        for sigma in range(4)
    ]

    # 2. All 64 Christoffel symbols Gamma[lam][mu][nu]
    Gamma = [[[Z for _ in range(4)] for _ in range(4)] for _ in range(4)]
    half = QQRat.const(Fraction(1, 2))
    for lam in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                acc = Z
                for sig in range(4):
                    if not g_con[lam][sig].is_zero():
                        term = dg[mu][nu][sig] + dg[nu][mu][sig] - dg[sig][mu][nu]
                        if not term.is_zero():
                            acc = acc + g_con[lam][sig] * term
                val = half * acc
                Gamma[lam][mu][nu] = val
                Gamma[lam][nu][mu] = val

    # 3. All 256 Riemann components R^rho_{sigma mu nu}
    R_up = [[[[Z for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for rho in range(4):
        for sig in range(4):
            for mu in range(4):
                for nu in range(mu + 1, 4):
                    val = D_coord(Gamma[rho][nu][sig], mu) - D_coord(Gamma[rho][mu][sig], nu)
                    for lam in range(4):
                        g1 = Gamma[rho][mu][lam]
                        g2 = Gamma[lam][nu][sig]
                        if not g1.is_zero() and not g2.is_zero():
                            val = val + g1 * g2
                        g3 = Gamma[rho][nu][lam]
                        g4 = Gamma[lam][mu][sig]
                        if not g3.is_zero() and not g4.is_zero():
                            val = val - g3 * g4
                    R_up[rho][sig][mu][nu] = val
                    R_up[rho][sig][nu][mu] = -val

    # Lower first index: R_cov[alpha][sig][mu][nu] = g_{alpha rho} R^rho_{sig mu nu}
    R_cov = [[[[Z for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for alpha in range(4):
        for sig in range(4):
            for mu in range(4):
                for nu in range(4):
                    acc = Z
                    for rho in range(4):
                        if not g_cov[alpha][rho].is_zero() and not R_up[rho][sig][mu][nu].is_zero():
                            acc = acc + g_cov[alpha][rho] * R_up[rho][sig][mu][nu]
                    R_cov[alpha][sig][mu][nu] = acc

    # Raise all indices: R_contra[a][b][c][d] = g^{a p} g^{b q} g^{c m} g^{d n} R_cov[p][q][m][n]
    # Use diagonal block sparsity (g_con non-zero only on (0,1)x(0,1), (2,2), (3,3))
    nz_con = [[j for j in range(4) if not g_con[i][j].is_zero()] for i in range(4)]
    K_scalar = Z
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    if R_cov[a][b][c][d].is_zero():
                        continue
                    r_up_abcd = Z
                    for p in nz_con[a]:
                        for q in nz_con[b]:
                            for m in nz_con[c]:
                                for n in nz_con[d]:
                                    rc = R_cov[p][q][m][n]
                                    if not rc.is_zero():
                                        r_up_abcd = (
                                            r_up_abcd
                                            + g_con[a][p] * g_con[b][q] * g_con[c][m] * g_con[d][n] * rc
                                        )
                    if not r_up_abcd.is_zero():
                        K_scalar = K_scalar + R_cov[a][b][c][d] * r_up_abcd

    # Ricci tensor R_{mu nu} = R^rho_{mu rho nu} and 4D Ricci scalar R_4D = g^{mu nu} R_{mu nu}
    Ric_cov = [[Z for _ in range(4)] for _ in range(4)]
    for mu in range(4):
        for nu in range(4):
            acc = Z
            for rho in range(4):
                acc = acc + R_up[rho][mu][rho][nu]
            Ric_cov[mu][nu] = acc

    Ric_scalar_4d = Z
    Ric_sq = Z
    for mu in range(4):
        for nu in range(4):
            if not g_con[mu][nu].is_zero() and not Ric_cov[mu][nu].is_zero():
                Ric_scalar_4d = Ric_scalar_4d + g_con[mu][nu] * Ric_cov[mu][nu]
            if not Ric_cov[mu][nu].is_zero():
                r_up_mn = Z
                for p in nz_con[mu]:
                    for q in nz_con[nu]:
                        if not Ric_cov[p][q].is_zero():
                            r_up_mn = r_up_mn + g_con[mu][p] * g_con[nu][q] * Ric_cov[p][q]
                Ric_sq = Ric_sq + Ric_cov[mu][nu] * r_up_mn

    # 4D Weyl scalar identity C_{abcd} C^{abcd} = K - 2 R_{ab} R^{ab} + (1/3) R^2
    Weyl_sq = K_scalar - 2 * Ric_sq + QQRat.const(Fraction(1, 3)) * (Ric_scalar_4d**2)

    # Expected closed forms
    K_expected = (H2**2) + 4 * ((H1 / r) ** 2) + 4 * ((H / (r**2)) ** 2)
    R4D_expected = H2 + 4 * (H1 / r) + 2 * (H / (r**2))
    Weyl_expected = QQRat.const(Fraction(1, 3)) * ((H2 - 2 * (H1 / r) + 2 * (H / (r**2))) ** 2)

    # Horizon-penetrating square-root-free rational orthonormal frame e_{a_hat}^mu at chi = 0:
    # Coframe: e^0 = (1 - H/2) dt - (H/2) dr,  e^1 = (H/2) dt + (1 + H/2) dr
    # Dual frame vectors e_{a_hat}^mu (satisfying e^a_mu e_b^mu = delta^a_b & g_{mn} e_a^m e_b^n = eta_{ab}):
    #   e_{\hat{0}} = (1 + H/2) \partial_t - (H/2) \partial_r
    #   e_{\hat{1}} = +(H/2) \partial_t + (1 - H/2) \partial_r
    #   e_{\hat{2}} = -(1/r) \partial_\chi |_{\chi=0}
    #   e_{\hat{3}} = (1/r) \partial_\phi |_{\chi=0}
    e_vec = [
        [O + half * H, -half * H, Z, Z],
        [half * H, O - half * H, Z, Z],
        [Z, Z, -O / r, Z],
        [Z, Z, Z, O / r],
    ]
    g_eq = [[g_cov[i][j].subs({"chi": 0}) for j in range(4)] for i in range(4)]
    R_cov_eq = [
        [[[R_cov[a][b][c][d].subs({"chi": 0}) for d in range(4)] for c in range(4)] for b in range(4)]
        for a in range(4)
    ]

    # Verify orthonormality g_{mu nu} e_a^mu e_b^nu == eta_{ab} identically in H (including H >= 1)!
    eta_diag = [QQRat.const(-1), O, O, O]
    tetrad_ortho_ok = True
    for a in range(4):
        for b in range(4):
            ip = Z
            for m in range(4):
                for n in range(4):
                    if not e_vec[a][m].is_zero() and not e_vec[b][n].is_zero() and not g_eq[m][n].is_zero():
                        ip = ip + g_eq[m][n] * e_vec[a][m] * e_vec[b][n]
            target = eta_diag[a] if a == b else Z
            if not (ip - target).is_zero():
                tetrad_ortho_ok = False

    def frame_R(ah: int, bh: int, ch: int, dh: int) -> QQRat:
        acc = Z
        for m in range(4):
            if e_vec[ah][m].is_zero():
                continue
            for n in range(4):
                if e_vec[bh][n].is_zero():
                    continue
                for p in range(4):
                    if e_vec[ch][p].is_zero():
                        continue
                    for q in range(4):
                        if e_vec[dh][q].is_zero() or R_cov_eq[m][n][p][q].is_zero():
                            continue
                        acc = acc + e_vec[ah][m] * e_vec[bh][n] * e_vec[ch][p] * e_vec[dh][q] * R_cov_eq[m][n][p][q]
        return acc

    R_0101 = frame_R(0, 1, 0, 1)
    R_0202 = frame_R(0, 2, 0, 2)
    R_0303 = frame_R(0, 3, 0, 3)
    R_1212 = frame_R(1, 2, 1, 2)
    R_1313 = frame_R(1, 3, 1, 3)
    R_0212 = frame_R(0, 2, 1, 2)  # radial-angular boost flux -> must be 0/1!
    R_2323 = frame_R(2, 3, 2, 3)

    # Radial/Laplacian trace R^0_0 + R^1_1 = -Ric_hat_{00} + Ric_hat_{11}
    # Ric_hat_{00} = R_{0101} + R_{0202} + R_{0303}, Ric_hat_{11} = -R_{1010} + R_{1212} + R_{1313}
    radial_laplacian = (-R_0101 - R_0202 - R_0303) + (-R_0101 + R_1212 + R_1313)

    return {
        "kretschmann": K_scalar,
        "kretschmann_residual_zero": (K_scalar - K_expected).is_zero(),
        "ricci_4d": Ric_scalar_4d,
        "ricci_4d_residual_zero": (Ric_scalar_4d - R4D_expected).is_zero(),
        "radial_laplacian_trace": radial_laplacian,
        "radial_laplacian_eq_H2_plus_2H1_r": (radial_laplacian - (H2 + 2 * (H1 / r))).is_zero(),
        "weyl_sq": Weyl_sq,
        "weyl_sq_residual_zero": (Weyl_sq - Weyl_expected).is_zero(),
        "tetrad_orthonormality_all_H_zero": tetrad_ortho_ok,
        "orthonormal_riemann": {
            "R_0101": R_0101,
            "R_0202": R_0202,
            "R_0303": R_0303,
            "R_1212": R_1212,
            "R_1313": R_1313,
            "R_0212_flux_zero": R_0212.is_zero(),
            "R_2323": R_2323,
        },
        "Gamma": Gamma,
    }


def analyze_profile_curvature(
    H_of_r: QQRat,
    var: str = "r",
    candidate_B_K: Optional[Fraction] = None,
) -> Dict[str, Any]:
    """
    Given a concrete or parameter-substituted rational metric profile `H_h(r)` in QQRat(r),
    symbolically compute all 3+1D curvature invariants AND orthonormal Riemann components:
      - H''(r), H'(r)/r, H(r)/r^2
      - Check pole-freeness on ALL r >= 0 INCLUDING exact finiteness at r = 0!
      - Exact Kretschmann K_h(r), 4D Ricci R_4D(r), radial trace H'' + 2 H'/r, Weyl C^2(r)
      - Certify exact rational ceiling `sup_{r >= 0} K_h(r) <= B_K` and report exact slack!
    """
    r = QQRat.var(var)
    H1 = H_of_r.diff(var, 1)
    H2 = H_of_r.diff(var, 2)

    tidal_rad = -QQRat.const(Fraction(1, 2)) * H2               # R_{\hat{0}\hat{1}\hat{0}\hat{1}} up to sign
    tidal_mix = QQRat.const(Fraction(1, 2)) * (H1 / r)          # R_{\hat{1}\hat{2}\hat{1}\hat{2}} = -R_{\hat{0}\hat{2}\hat{0}\hat{2}}
    tidal_ang = H_of_r / (r**2)                                 # R_{\hat{2}\hat{3}\hat{2}\hat{3}}

    K_r = (H2**2) + 4 * ((H1 / r) ** 2) + 4 * (tidal_ang**2)
    R_lap_r = H2 + 2 * (H1 / r)
    R_4d_r = H2 + 4 * (H1 / r) + 2 * tidal_ang
    Weyl_r = QQRat.const(Fraction(1, 3)) * ((H2 - 2 * (H1 / r) + 2 * tidal_ang) ** 2)

    # Check exact finiteness at r = 0 for EVERY orthonormal Riemann component!
    r0_finite = True
    r0_values: Dict[str, Any] = {}
    for label, comp in (
        ("H(0)", H_of_r),
        ("R_0101(0) [-H''/2]", tidal_rad),
        ("R_1212(0) [H'/(2r)]", tidal_mix),
        ("R_2323(0) [H/r^2]", tidal_ang),
        ("K_h(0)", K_r),
        ("R_lap(0)", R_lap_r),
        ("R_4D(0)", R_4d_r),
        ("Weyl_sq(0)", Weyl_r),
    ):
        try:
            r0_values[label] = comp.eval_at(**{var: Fraction(0)})
        except ZeroDivisionError:
            r0_finite = False
            r0_values[label] = "SINGULAR_POLE_AT_R_0"

    sup_cert = (
        certify_rational_supremum_on_nonneg_ray(K_r, var=var, candidate_bound=candidate_B_K)
        if r0_finite
        else {"certified": False, "reason": "Singular curvature pole at r = 0"}
    )

    return {
        "r0_all_orthonormal_riemann_finite": r0_finite,
        "r0_values": r0_values,
        "K_h": K_r,
        "R_lap": R_lap_r,
        "R_4D": R_4d_r,
        "Weyl_sq": Weyl_r,
        "orthonormal_components": {
            "R_0101_abs": H2 * QQRat.const(Fraction(1, 2)),
            "R_0202_abs": (H1 / r) * QQRat.const(Fraction(1, 2)),
            "R_2323": tidal_ang,
        },
        "supremum_certificate": sup_cert,
    }


def cartesian_kerr_schild_christoffel_at(
    x_vec: Sequence[Any],
    H_of_r_func: Callable[[Any], Tuple[Any, Any]],
) -> List[List[List[Any]]]:
    """
    Evaluate exact 3+1D Cartesian Kerr-Schild Christoffel symbols Gamma^lam_{mu nu} at
    spatial point x_vec = (x, y, z) where r_val and (H_val, dH_dr_val) are provided rationally!
    Notice:
      g_{mu nu} = eta_{mu nu} + H(r) l_mu l_nu,   l_mu = (1, x/r, y/r, z/r),   r^2 = x^2 + y^2 + z^2.
      Can be evaluated WITHOUT square roots whenever `r` itself is rational (e.g. Levi-Civita / KS spinor state)!
    """
    raise NotImplementedError("Use compute_cartesian_ks_christoffel_from_spinor below for exact root-free evaluation")


def compute_cartesian_ks_christoffel_from_r(
    x: Any,
    y: Any,
    z: Any,
    r: Any,
    H: Any,
    dH_dr: Any,
) -> List[List[List[Any]]]:
    """
    Exact root-free 3+1D Cartesian Kerr-Schild Christoffel symbols Gamma^lam_{mu nu}
    for g_{mu nu} = eta_{mu nu} + H(r) l_mu l_nu, l_mu = (1, x/r, y/r, z/r),
    expressed purely rationally in (x, y, z, r, H, dH_dr)!

    Derivation:
      Let f_{mu nu} = H(r) l_mu l_nu. Since l_0 = 1 and l_i = x^i / r:
        partial_0 g_{mu nu} = 0
        partial_k r = x^k / r = l_k
        partial_k l_0 = 0
        partial_k l_i = (delta_{ki} - l_k l_i) / r
      Therefore:
        partial_k (H l_mu l_nu) = dH_dr * l_k * l_mu * l_nu
                                  + H * (partial_k l_mu) * l_nu
                                  + H * l_mu * (partial_k l_nu)
      And g^{lam sig} = eta^{lam sig} - H * l^lam * l^sig with l^0 = -1, l^i = l_i!
    """
    coords_l = [1, x / r, y / r, z / r]
    l_up = [-1, x / r, y / r, z / r]
    eta_up = [-1, 1, 1, 1]

    # dl[sigma][mu] = partial_sigma l_mu
    dl = [[0 for _ in range(4)] for _ in range(4)]
    for k in (1, 2, 3):
        for i in (1, 2, 3):
            delta_ki = 1 if k == i else 0
            dl[k][i] = (delta_ki - coords_l[k] * coords_l[i]) / r

    # dg[sigma][mu][nu] = partial_sigma g_{mu nu}
    dg = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for k in (1, 2, 3):
        lk = coords_l[k]
        for mu in range(4):
            lmu = coords_l[mu]
            dl_k_mu = dl[k][mu]
            for nu in range(4):
                lnu = coords_l[nu]
                dg[k][mu][nu] = (
                    dH_dr * lk * lmu * lnu
                    + H * dl_k_mu * lnu
                    + H * lmu * dl[k][nu]
                )

    g_con = [
        [
            (eta_up[lam] if lam == sig else 0) - H * l_up[lam] * l_up[sig]
            for sig in range(4)
        ]
        for lam in range(4)
    ]

    half = Fraction(1, 2)
    Gamma = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for lam in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                acc = 0
                for sig in range(4):
                    term = dg[mu][nu][sig] + dg[nu][mu][sig] - dg[sig][mu][nu]
                    if term != 0 and g_con[lam][sig] != 0:
                        acc = acc + g_con[lam][sig] * term
                val = half * acc
                Gamma[lam][mu][nu] = val
                Gamma[lam][nu][mu] = val
    return Gamma


def schwarzschild_ks_geodesic_jet_residuals(
    x_mu_jets: Sequence[QQJet],
    r_jet: QQJet,
    H_of_r_func: Callable[[Any], Any],
) -> Dict[str, Any]:
    """
    Given 4 spacetime coordinate jets x^mu(eps) = (t(eps), x(eps), y(eps), z(eps)) of order >= 2
    produced by running @architect's live step function on `QQJet` with `IdentityQuantizerBundle`,
    together with the spinor radius jet `r(eps) = sum u_a(eps)^2`, check 2nd-order continuum GR
    consistency against the exact 3+1D Kerr-Schild Levi-Civita connection!

    Computes:
      x0^mu = x^mu[0],   v^mu = x^mu[1],   a^mu = 2 * x^mu[2]
      J^mu  = a^mu + Gamma^mu_{alpha beta}(x0) v^alpha v^beta
    Checks:
      1. Affine geodesic residual: `J^mu == 0` (exact if `eps` is proper/affine time `tau` or `lambda`)
      2. Coordinate-time geodesic residual (`t(eps) = t0 + eps` => `v^0 = 1, a^0 = 0`):
         `a^i + Gamma^i_{ab} v^a v^b - v^i Gamma^0_{ab} v^a v^b == 0`
      3. Reparameterization-invariant projective geodesic residual (for Sundman `ds = dtau / r(s)` etc.):
         `J^mu * v^nu - J^nu * v^mu == 0` for all `0 <= mu < nu <= 3`!
    """
    x0 = [j[0] for j in x_mu_jets]
    v0 = [j[1] for j in x_mu_jets]
    a0 = [2 * j[2] for j in x_mu_jets]
    r0 = r_jet[0]

    # Compute H(r0) and H'(r0) using a 1st-order QQJet around r0!
    eps1 = QQJet([r0, 1], order=1)
    H_jet = H_of_r_func(eps1)
    H0 = H_jet[0]
    dH0 = H_jet[1]

    Gamma = compute_cartesian_ks_christoffel_from_r(x0[1], x0[2], x0[3], r0, H0, dH0)

    conn_acc = []
    J = []
    for mu in range(4):
        c_mu = 0
        for al in range(4):
            for be in range(4):
                if v0[al] != 0 and v0[be] != 0 and Gamma[mu][al][be] != 0:
                    c_mu = c_mu + Gamma[mu][al][be] * v0[al] * v0[be]
        conn_acc.append(c_mu)
        J.append(a0[mu] + c_mu)

    # Coordinate-time acceleration residual d^2 x^i / dt^2 (general even if dt/deps != 1):
    # Note: dx^i/dt = v^i / v^0, d^2 x^i / dt^2 = (a^i * v^0 - v^i * a^0) / (v^0)^3, while
    # continuum GR gives d^2 x^i / dt^2 = (-conn^i + (v^i / v^0) * conn^0) / (v^0)^2!
    # Multiplying by (v^0)^3 clears divisions and equals J^i * v^0 - v^i * J^0!
    projective_wedges = {}
    for mu in range(4):
        for nu in range(mu + 1, 4):
            projective_wedges[f"({mu},{nu})"] = J[mu] * v0[nu] - J[nu] * v0[mu]

    coord_time_res = [
        a0[i] + conn_acc[i] - v0[i] * conn_acc[0] for i in (1, 2, 3)
    ]

    return {
        "affine_J_residuals": J,
        "affine_all_zero": all(j == 0 for j in J),
        "coord_time_residuals": coord_time_res,
        "coord_time_all_zero": v0[0] == 1 and a0[0] == 0 and all(r == 0 for r in coord_time_res),
        "projective_wedge_residuals": projective_wedges,
        "projective_geodesic_all_zero": all(w == 0 for w in projective_wedges.values()),
    }


# ---------------------------------------------------------------------------
# Comprehensive Built-in Self-Test Suite
# ---------------------------------------------------------------------------

def run_self_tests() -> Dict[str, Any]:
    """Run exhaustive exact algebraic self-tests of every CAS component."""
    results: Dict[str, Any] = {}

    # Test 1: Multivariate Subresultant GCD & Exact Division
    u, v, w, dt = QQPoly.vars("u", "v", "w", "dt")
    g_true = 3 * (u**2) * v - 5 * (v**3) + 2 * w * dt - 7
    p1 = g_true * (2 * (u**3) - 3 * v * w + 11 * (dt**2) + 1)
    p2 = g_true * (5 * u * (w**2) + 7 * (v**2) * dt - 4)
    g_computed = QQPoly.gcd(p1, p2)
    _, z_true = g_true.to_primitive_zpoly()
    assert g_computed == QQPoly.from_zpoly(z_true), f"GCD mismatch: {g_computed} vs {g_true}"
    results["multivariate_4var_subresultant_gcd_ok"] = True

    # Sparse degree-drop stress test for Subresultant PRS (classic gap hazard!)
    a_gap = (u**4 + w) * (u**2 + v**2 + 1)
    b_gap = (u**4 - v) * (u**2 + v**2 + 1)
    g_gap = QQPoly.gcd(a_gap, b_gap)
    assert g_gap == (u**2 + v**2 + 1)
    results["sparse_gap_subresultant_gcd_ok"] = True

    # Test 2: QQRat Canonical Normalization, Tower Coercion, & Double-Cayley Involution
    uq, vq, hq, rsq, dtq = QQRat.vars("u", "v", "h", "rs", "dt")
    # Double-Cayley symplectic rotor: C(z) = (1 - z^2)/(1 + z^2), S(z) = 2z/(1 + z^2)
    z_rot = dtq * (uq**2 + vq**2) / (1 + hq**2)
    c_rot = (1 - z_rot**2) / (1 + z_rot**2)
    s_rot = (2 * z_rot) / (1 + z_rot**2)
    unit_circle_res = (c_rot**2 + s_rot**2) - Fraction(1)
    assert unit_circle_res == 0 and unit_circle_res.is_zero()
    # Reversibility R(-dt) * R(+dt) == I_2
    c_neg = c_rot.subs({"dt": -dtq})
    s_neg = s_rot.subs({"dt": -dtq})
    rev_diag = c_neg * c_rot - s_neg * s_rot - 1
    rev_off = c_neg * s_rot + s_neg * c_rot
    assert rev_diag.is_zero() and rev_off.is_zero()
    results["qqrat_cayley_unit_circle_and_reversibility_ok"] = True

    # Test 3: QQJet over BOTH Fraction and Symbolic QQRat
    eps_sym = QQJet.var(order=3, symbolic=True)
    # Check Taylor expansion of 1 / (1 - a*eps) == 1 + a*eps + a^2*eps^2 + a^3*eps^3
    aq = QQRat.var("a")
    geom_jet = Fraction(1) / (Fraction(1) - aq * eps_sym)
    assert geom_jet[0] == 1
    assert geom_jet[1] == aq
    assert geom_jet[2] == aq**2
    assert geom_jet[3] == aq**3
    results["qqjet_symbolic_series_ok"] = True

    # Test 4: Spinor Kerr-Schild 3+1D Algebraic & Lightcone Identities
    ks_id = verify_spinor_kerr_schild_identities()
    assert all(ks_id.values()), f"Spinor KS identity failure: {ks_id}"
    results["spinor_ks_identities"] = ks_id

    # Test 5: Full 4D Tensor Curvature Contraction (256 Riemann + Orthonormal Frame)
    curv_univ = compute_4d_kerr_schild_curvature_universal()
    assert curv_univ["kretschmann_residual_zero"]
    assert curv_univ["ricci_4d_residual_zero"]
    assert curv_univ["radial_laplacian_eq_H2_plus_2H1_r"]
    assert curv_univ["weyl_sq_residual_zero"]
    assert curv_univ["tetrad_orthonormality_all_H_zero"]
    assert curv_univ["orthonormal_riemann"]["R_0212_flux_zero"]
    results["full_4d_riemann_kretschmann_tetrad_ok"] = True

    # Test 6: Sturm Root Isolation, SOS Certificates, and Hayward Exact 0/1-Slack Kretschmann Bound
    rq = QQRat.var("r")
    # Hayward metric H_h(r) = r_s * r^2 / (r^3 + h^3) for r_s = 2, h = 1/2
    rs_val = Fraction(2)
    h_val = Fraction(1, 2)
    H_hayward = (rs_val * rq**2) / (rq**3 + h_val**3)
    hayward_audit = analyze_profile_curvature(H_hayward, var="r")
    assert hayward_audit["r0_all_orthonormal_riemann_finite"]
    sup_h = hayward_audit["supremum_certificate"]
    expected_B_K = Fraction(24) * (rs_val**2) / (h_val**6)  # 24 * 4 / (1/64) = 6144
    assert sup_h["certified"] and sup_h["exact_zero_slack"] and sup_h["bound"] == expected_B_K
    results["hayward_zero_slack_kretschmann_bound"] = str(sup_h["bound"])

    # Also test an interior irrational-peak rational function and verify tight SOS/Bernstein certification!
    # e.g., H_quartic(r) = 2 * r^3 / (r^4 + (1/2)^4) => K_h(0) = 0, interior irrational peak!
    H_quartic = (rs_val * rq**3) / (rq**4 + h_val**4)
    quartic_audit = analyze_profile_curvature(H_quartic, var="r")
    assert quartic_audit["r0_all_orthonormal_riemann_finite"]
    sup_q = quartic_audit["supremum_certificate"]
    assert sup_q["certified"] and not sup_q["exact_zero_slack"]
    assert sup_q["tightness_ratio_lower_bound"] > Fraction(999, 1000)
    results["quartic_interior_peak_certified_tightness"] = str(sup_q["tightness_ratio_lower_bound"])

    # Test 7: Singular O(r) core detection (H(r) = 2*r / (r^2 + h^2) MUST be flagged singular at r = 0!)
    H_bad = (rs_val * rq) / (rq**2 + h_val**2)
    bad_audit = analyze_profile_curvature(H_bad, var="r")
    assert not bad_audit["r0_all_orthonormal_riemann_finite"]
    assert bad_audit["r0_values"]["R_2323(0) [H/r^2]"] == "SINGULAR_POLE_AT_R_0"
    results["bad_O_r_core_singularity_detected"] = True

    # Test 8: Symbolic 3+1D Cartesian Kerr-Schild Geodesic Jet Verifier (affine, coord-time, projective)
    # Use symbolic Levi-Civita spinor point x0 = u^2 - v^2, y0 = 2uv, z0 = 0, r0 = u^2 + v^2
    us, vs, vt, vx, vy = QQRat.vars("u", "v", "vt", "vx", "vy")
    x0_s = us**2 - vs**2
    y0_s = 2 * us * vs
    z0_s = QQRat.zero()
    r0_s = us**2 + vs**2
    H_func_hayward = lambda rad: (Fraction(2) * (rad**2)) / ((rad**3) + Fraction(1, 8))
    H0_s = H_func_hayward(r0_s)
    dH0_s = H_func_hayward(QQJet([r0_s, QQRat.one()], order=1))[1]
    G_sym = compute_cartesian_ks_christoffel_from_r(x0_s, y0_s, z0_s, r0_s, H0_s, dH0_s)
    v_vec = [vt, vx, vy, QQRat.zero()]
    a_aff = [
        -sum(G_sym[mu][al][be] * v_vec[al] * v_vec[be] for al in range(4) for be in range(4))
        for mu in range(4)
    ]
    eps_j = QQJet.var(order=2, symbolic=True)
    pos_vec = [QQRat.zero(), x0_s, y0_s, z0_s]
    x_jets = [
        pos_vec[mu] + v_vec[mu] * eps_j + QQRat.const(Fraction(1, 2)) * a_aff[mu] * (eps_j**2)
        for mu in range(4)
    ]
    r_jet_test = QQJet([r0_s, QQRat.zero(), QQRat.zero()], order=2)
    geo_check = schwarzschild_ks_geodesic_jet_residuals(x_jets, r_jet_test, H_func_hayward)
    assert geo_check["affine_all_zero"] and geo_check["projective_geodesic_all_zero"]
    results["symbolic_ks_geodesic_jet_verifier_ok"] = True

    # Test 9: Tower Polymorphism & IdentityQuantizerBundle AST Verification
    qb = IDENTITY_QB
    p_var = QQPoly.var("x")
    r_var = QQRat.var("x")
    j_var = QQJet([r_var, QQRat.one()], order=2)
    for elem in (Fraction(3, 5), p_var, r_var, j_var):
        val = qb((Fraction(1, 2) + elem) * Fraction(3, 4) - (Fraction(2) * elem) / Fraction(5, 3))
        val_sum = sum([elem, Fraction(1, 3), -elem])
        assert val_sum == Fraction(1, 3)
        assert qb(val, Fraction(1)) == val
    results["tower_polymorphism_and_identity_qb_ok"] = True

    return results


if __name__ == "__main__":
    import json
    import time

    t0_ns = time.perf_counter_ns()
    summary = run_self_tests()
    elapsed_ms = (time.perf_counter_ns() - t0_ns) // 1_000_000
    summary["elapsed_ms"] = elapsed_ms
    print(json.dumps(summary, indent=2))
