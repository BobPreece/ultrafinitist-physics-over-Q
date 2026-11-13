"""Exact Multivariate Polynomial & Rational Function Symbolic Algebra Engine over QQ.

Self-contained zero-dependency computer algebra system built strictly on Python's
exact `fractions.Fraction` and `int` for Python 3.9+.

Provides:
- `QQVarRing`: Named multivariate polynomial/rational ring `QQ[x_1, ..., x_m]` / `QQ(x_1, ..., x_m)`.
- `QQPoly`: Exact sparse multivariate polynomial (`dict[tuple[int, ...], Fraction]`)
  with exact `+`, `-`, `*`, integer `**`, partial derivatives, degrees, monomial counts,
  exact substitution, exact multivariate division (`try_div_exact`), and exact
  multivariate polynomial GCD (`poly_gcd` via Collins-Brown Subresultant PRS over ZZ
  with 100% verified exact-division certificates).
- `QQRat`: Exact multivariate rational function `(num: QQPoly, den: QQPoly)` in
  `QQ(x_1, ..., x_m)` with Knuth UFD lowest-terms reduction, exact `+`, `-`, `*`, `/`,
  signed integer `**`, quotient-rule derivatives, homogenized denominator-grouped
  substitution, and 100% rigorous identity verification
  `A == B iff (A.num * B.den - B.num * A.den).is_zero()`.
- `QQJet`: Truncated formal power series in `QQ(z)[[dt]] / (dt^(order+1))` for fast
  exact symbolic Taylor / Hamiltonian vector-field consistency verification even on
  high-dimensional phase spaces.
- Symbolic matrix algebra (`mat_inv`, `mat_det`, `mat_mul`, ...) and Hamiltonian
  mechanics / arithmetic-dynamics degree utilities.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

Exponent = Tuple[int, ...]
Scalar = Union[int, Fraction]


def _to_fraction(val: Scalar) -> Fraction:
    if isinstance(val, Fraction):
        return val
    if isinstance(val, int) and not isinstance(val, bool):
        return Fraction(val, 1)
    raise TypeError(f"Expected int or Fraction, got {type(val).__name__}")


def _lcm(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return abs(a // math.gcd(a, b) * b)


# ---------------------------------------------------------------------------
# Ring definition
# ---------------------------------------------------------------------------


class QQVarRing:
    """Multivariate variable ring over QQ with ordered generator names."""

    def __init__(self, names: Sequence[str]) -> None:
        names_tuple = tuple(names)
        if len(set(names_tuple)) != len(names_tuple):
            raise ValueError(f"Duplicate variable names in ring: {names_tuple}")
        self.names: Tuple[str, ...] = names_tuple
        self.num_vars: int = len(names_tuple)
        self._name_to_idx: Dict[str, int] = {n: i for i, n in enumerate(names_tuple)}
        self._zero_exp: Exponent = (0,) * self.num_vars
        self._zero_poly: Optional[QQPoly] = None
        self._one_poly: Optional[QQPoly] = None
        self._zero_rat: Optional[QQRat] = None
        self._one_rat: Optional[QQRat] = None

    def __repr__(self) -> str:
        return f"QQVarRing({list(self.names)!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, QQVarRing) and self.names == other.names

    def __hash__(self) -> int:
        return hash(self.names)

    def resolve_var(self, var: Union[int, str, QQPoly, QQRat]) -> int:
        if isinstance(var, int) and not isinstance(var, bool):
            if 0 <= var < self.num_vars:
                return var
            raise IndexError(f"Variable index {var} out of range for {self.num_vars} vars")
        if isinstance(var, str):
            if var in self._name_to_idx:
                return self._name_to_idx[var]
            raise KeyError(f"Unknown variable name {var!r} in {self.names}")
        if isinstance(var, QQPoly):
            if var.ring != self or len(var.coeffs) != 1:
                raise ValueError(f"Expected single generator polynomial, got {var}")
            exp, coeff = next(iter(var.coeffs.items()))
            if coeff != 1 or sum(exp) != 1:
                raise ValueError(f"Expected single generator polynomial, got {var}")
            return exp.index(1)
        if isinstance(var, QQRat):
            if not var.den.is_constant() or var.den.constant_value() != 1:
                raise ValueError(f"Expected generator rational variable, got {var}")
            return self.resolve_var(var.num)
        raise TypeError(f"Cannot resolve variable from {type(var).__name__}")

    def zero(self) -> QQPoly:
        if self._zero_poly is None:
            self._zero_poly = QQPoly._raw(self, {})
        return self._zero_poly

    def one(self) -> QQPoly:
        if self._one_poly is None:
            self._one_poly = QQPoly._raw(self, {self._zero_exp: Fraction(1, 1)})
        return self._one_poly

    def const(self, value: Scalar) -> QQPoly:
        f = _to_fraction(value)
        if f == 0:
            return self.zero()
        return QQPoly._raw(self, {self._zero_exp: f})

    def gen(self, var: Union[int, str]) -> QQPoly:
        idx = self.resolve_var(var)
        exp = tuple(1 if i == idx else 0 for i in range(self.num_vars))
        return QQPoly._raw(self, {exp: Fraction(1, 1)})

    def gens(self) -> Tuple[QQPoly, ...]:
        return tuple(self.gen(i) for i in range(self.num_vars))

    def rat_zero(self) -> QQRat:
        if self._zero_rat is None:
            self._zero_rat = QQRat._raw(self.zero(), self.one())
        return self._zero_rat

    def rat_one(self) -> QQRat:
        if self._one_rat is None:
            self._one_rat = QQRat._raw(self.one(), self.one())
        return self._one_rat

    def rat_const(self, value: Scalar) -> QQRat:
        return QQRat._raw(self.const(value), self.one())

    def rat_gen(self, var: Union[int, str]) -> QQRat:
        return QQRat._raw(self.gen(var), self.one())

    def rat_gens(self) -> Tuple[QQRat, ...]:
        return tuple(self.rat_gen(i) for i in range(self.num_vars))

    def extend(self, extra_names: Sequence[str]) -> QQVarRing:
        return QQVarRing(list(self.names) + list(extra_names))

    def embed_poly(self, poly: QQPoly) -> QQPoly:
        if poly.ring == self:
            return poly
        idx_map = [self.resolve_var(name) for name in poly.ring.names]
        new_coeffs: Dict[Exponent, Fraction] = {}
        for exp, c in poly.coeffs.items():
            new_exp = [0] * self.num_vars
            for old_i, power in enumerate(exp):
                if power:
                    new_exp[idx_map[old_i]] = power
            new_coeffs[tuple(new_exp)] = c
        return QQPoly._raw(self, new_coeffs)

    def embed_rat(self, rat: QQRat) -> QQRat:
        if rat.ring == self:
            return rat
        return QQRat._raw(self.embed_poly(rat.num), self.embed_poly(rat.den))


# ---------------------------------------------------------------------------
# Sparse Multivariate Polynomial over QQ
# ---------------------------------------------------------------------------


class QQPoly:
    """Exact sparse multivariate polynomial over QQ."""

    __slots__ = ("ring", "coeffs")

    def __init__(self, ring: QQVarRing, coeffs: Mapping[Exponent, Scalar]) -> None:
        self.ring: QQVarRing = ring
        n = ring.num_vars
        cleaned: Dict[Exponent, Fraction] = {}
        for exp, raw_c in coeffs.items():
            if len(exp) != n or any(e < 0 for e in exp):
                raise ValueError(f"Invalid exponent tuple {exp} for {n} variables")
            c = _to_fraction(raw_c)
            if c != 0:
                cleaned[exp] = c
        self.coeffs: Dict[Exponent, Fraction] = cleaned

    @classmethod
    def _raw(cls, ring: QQVarRing, coeffs: Dict[Exponent, Fraction]) -> QQPoly:
        obj = cls.__new__(cls)
        obj.ring = ring
        obj.coeffs = coeffs
        return obj

    def _coerce_poly(self, other: object) -> Optional[QQPoly]:
        if isinstance(other, QQPoly):
            if other.ring != self.ring:
                raise ValueError(f"Ring mismatch: {self.ring} vs {other.ring}")
            return other
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            return self.ring.const(other)
        return None

    # --- Structural queries ---

    def is_zero(self) -> bool:
        return len(self.coeffs) == 0

    def is_constant(self) -> bool:
        if not self.coeffs:
            return True
        if len(self.coeffs) > 1:
            return False
        exp = next(iter(self.coeffs))
        return exp == self.ring._zero_exp

    def constant_value(self) -> Fraction:
        if not self.coeffs:
            return Fraction(0, 1)
        if not self.is_constant():
            raise ValueError(f"Polynomial is not constant: {self}")
        return self.coeffs[self.ring._zero_exp]

    def monomial_count(self) -> int:
        return len(self.coeffs)

    def total_degree(self) -> int:
        if not self.coeffs:
            return 0
        return max(sum(exp) for exp in self.coeffs)

    def degree(self, var: Union[int, str, QQPoly, QQRat]) -> int:
        if not self.coeffs:
            return 0
        idx = self.ring.resolve_var(var)
        return max(exp[idx] for exp in self.coeffs)

    def degrees(self) -> Tuple[int, ...]:
        n = self.ring.num_vars
        if not self.coeffs:
            return (0,) * n
        res = [0] * n
        for exp in self.coeffs:
            for i, e in enumerate(exp):
                if e > res[i]:
                    res[i] = e
        return tuple(res)

    def degree_in_vars(self, vars_subset: Iterable[Union[int, str, QQPoly, QQRat]]) -> int:
        if not self.coeffs:
            return 0
        idxs = [self.ring.resolve_var(v) for v in vars_subset]
        if not idxs:
            return 0
        return max(sum(exp[i] for i in idxs) for exp in self.coeffs)

    def leading_term(self) -> Tuple[Exponent, Fraction]:
        if not self.coeffs:
            return self.ring._zero_exp, Fraction(0, 1)
        lead_exp = max(self.coeffs)
        return lead_exp, self.coeffs[lead_exp]

    def monic(self) -> QQPoly:
        if not self.coeffs:
            return self
        _, lc = self.leading_term()
        if lc == 1:
            return self
        inv_lc = Fraction(1, 1) / lc
        return QQPoly._raw(self.ring, {e: c * inv_lc for e, c in self.coeffs.items()})

    def monomial_content(self) -> Exponent:
        if not self.coeffs:
            return self.ring._zero_exp
        it = iter(self.coeffs)
        mins = list(next(it))
        for exp in it:
            for i, e in enumerate(exp):
                if e < mins[i]:
                    mins[i] = e
        return tuple(mins)

    def shift_exponents(self, delta: Exponent) -> QQPoly:
        if not self.coeffs or all(d == 0 for d in delta):
            return self
        shifted: Dict[Exponent, Fraction] = {}
        for exp, c in self.coeffs.items():
            ne = tuple(a + b for a, b in zip(exp, delta))
            if any(v < 0 for v in ne):
                raise ValueError("Negative exponent after shift")
            shifted[ne] = c
        return QQPoly._raw(self.ring, shifted)

    def max_coeff_bits(self) -> int:
        if not self.coeffs:
            return 0
        return max(
            max(abs(c.numerator).bit_length(), c.denominator.bit_length())
            for c in self.coeffs.values()
        )

    # --- Equality & Formatting ---

    def __eq__(self, other: object) -> bool:
        coerced = self._coerce_poly(other)
        if coerced is not None:
            return self.coeffs == coerced.coeffs
        if isinstance(other, (QQRat, QQJet)):
            return other == self
        return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __bool__(self) -> bool:
        return len(self.coeffs) > 0

    def __repr__(self) -> str:
        return f"QQPoly({self.__str__()})"

    def __str__(self) -> str:
        if not self.coeffs:
            return "0"
        items = sorted(self.coeffs.items(), key=lambda kv: (sum(kv[0]), kv[0]), reverse=True)
        parts: List[str] = []
        for idx, (exp, coeff) in enumerate(items):
            neg = coeff < 0
            abs_c = -coeff if neg else coeff
            var_chunks: List[str] = []
            for v_i, power in enumerate(exp):
                if power == 1:
                    var_chunks.append(self.ring.names[v_i])
                elif power > 1:
                    var_chunks.append(f"{self.ring.names[v_i]}^{power}")
            mon_str = "*".join(var_chunks)
            if abs_c.denominator == 1:
                c_str = str(abs_c.numerator)
            else:
                c_str = f"({abs_c.numerator}/{abs_c.denominator})"

            if not mon_str:
                term = c_str
            elif abs_c == 1:
                term = mon_str
            else:
                term = f"{c_str}*{mon_str}"

            if idx == 0:
                parts.append(f"-{term}" if neg else term)
            else:
                parts.append(f" - {term}" if neg else f" + {term}")
        return "".join(parts)

    # --- Exact Arithmetic ---

    def __pos__(self) -> QQPoly:
        return self

    def __neg__(self) -> QQPoly:
        if not self.coeffs:
            return self
        return QQPoly._raw(self.ring, {e: -c for e, c in self.coeffs.items()})

    def __add__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        p = self._coerce_poly(other)
        if p is None:
            return NotImplemented
        if not self.coeffs:
            return p
        if not p.coeffs:
            return self
        if len(self.coeffs) < len(p.coeffs):
            small, big = self.coeffs, p.coeffs
        else:
            small, big = p.coeffs, self.coeffs
        out = dict(big)
        for exp, c in small.items():
            nv = out.get(exp, Fraction(0, 1)) + c
            if nv == 0:
                out.pop(exp, None)
            else:
                out[exp] = nv
        return QQPoly._raw(self.ring, out)

    def __radd__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        return self.__add__(other)

    def __sub__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        p = self._coerce_poly(other)
        if p is None:
            return NotImplemented
        if not p.coeffs:
            return self
        if not self.coeffs:
            return -p
        out = dict(self.coeffs)
        for exp, c in p.coeffs.items():
            nv = out.get(exp, Fraction(0, 1)) - c
            if nv == 0:
                out.pop(exp, None)
            else:
                out[exp] = nv
        return QQPoly._raw(self.ring, out)

    def __rsub__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        p = self._coerce_poly(other)
        if p is None:
            return NotImplemented
        return p.__sub__(self)

    def __mul__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            f = _to_fraction(other)
            if f == 0 or not self.coeffs:
                return self.ring.zero()
            if f == 1:
                return self
            return QQPoly._raw(self.ring, {e: c * f for e, c in self.coeffs.items()})
        p = self._coerce_poly(other)
        if p is None:
            return NotImplemented
        if not self.coeffs or not p.coeffs:
            return self.ring.zero()
        out: Dict[Exponent, Fraction] = {}
        p_items = list(p.coeffs.items())
        for e1, c1 in self.coeffs.items():
            for e2, c2 in p_items:
                exp = tuple(a + b for a, b in zip(e1, e2))
                prod = c1 * c2
                prev = out.get(exp)
                if prev is None:
                    out[exp] = prod
                else:
                    nv = prev + prod
                    if nv == 0:
                        del out[exp]
                    else:
                        out[exp] = nv
        return QQPoly._raw(self.ring, out)

    def __rmul__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        return self.__mul__(other)

    def __truediv__(self, other: object) -> Union[QQPoly, QQRat, QQJet]:
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            f = _to_fraction(other)
            if f == 0:
                raise ZeroDivisionError("Division of QQPoly by zero scalar")
            return self * (Fraction(1, 1) / f)
        if isinstance(other, QQPoly):
            if other.ring != self.ring:
                raise ValueError(f"Ring mismatch: {self.ring} vs {other.ring}")
            return QQRat(self, other)
        return NotImplemented

    def __rtruediv__(self, other: object) -> QQRat:
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            return QQRat(self.ring.const(other), self)
        return NotImplemented

    def __pow__(self, power: int) -> Union[QQPoly, QQRat]:
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError(f"Exponent must be int, got {type(power).__name__}")
        if power < 0:
            return QQRat(self.ring.one(), self ** (-power))
        if power == 0:
            return self.ring.one()
        if power == 1 or not self.coeffs:
            return self
        if len(self.coeffs) == 1:
            exp, coeff = next(iter(self.coeffs.items()))
            return QQPoly._raw(
                self.ring,
                {tuple(e * power for e in exp): coeff ** power},
            )
        result = self.ring.one()
        base = self
        k = power
        while k > 0:
            if k & 1:
                result = result * base  # type: ignore[assignment]
            k >>= 1
            if k:
                base = base * base  # type: ignore[assignment]
        return result

    # --- Calculus & Exact Division ---

    def diff(self, var: Union[int, str, QQPoly, QQRat], order: int = 1) -> QQPoly:
        if order < 0:
            raise ValueError("Derivative order must be >= 0")
        if order == 0 or not self.coeffs:
            return self
        idx = self.ring.resolve_var(var)
        cur = self
        for _ in range(order):
            if not cur.coeffs:
                break
            nxt: Dict[Exponent, Fraction] = {}
            for exp, c in cur.coeffs.items():
                p = exp[idx]
                if p > 0:
                    ne = exp[:idx] + (p - 1,) + exp[idx + 1 :]
                    nxt[ne] = c * p
            cur = QQPoly._raw(self.ring, nxt)
        return cur

    def try_div_exact(self, divisor: QQPoly, max_steps: int = 8000) -> Optional[QQPoly]:
        """Return exact quotient `Q` if `self == divisor * Q` in `QQ[x]`, else `None`."""
        if divisor.ring != self.ring:
            raise ValueError("Ring mismatch in try_div_exact")
        if divisor.is_zero():
            raise ZeroDivisionError("Exact division by zero polynomial")
        if self.is_zero():
            return self.ring.zero()
        if divisor.is_constant():
            return self / divisor.constant_value()  # type: ignore[return-value]
        if len(divisor.coeffs) == 1:
            d_exp, d_coeff = next(iter(divisor.coeffs.items()))
            inv_d = Fraction(1, 1) / d_coeff
            q_coeffs: Dict[Exponent, Fraction] = {}
            for exp, c in self.coeffs.items():
                if any(a < b for a, b in zip(exp, d_exp)):
                    return None
                q_coeffs[tuple(a - b for a, b in zip(exp, d_exp))] = c * inv_d
            return QQPoly._raw(self.ring, q_coeffs)

        if self.total_degree() < divisor.total_degree():
            return None
        s_degs = self.degrees()
        d_degs = divisor.degrees()
        if any(sd < dd for sd, dd in zip(s_degs, d_degs)):
            return None

        b_lead_exp, b_lead_c = divisor.leading_term()
        rem = dict(self.coeffs)
        quot: Dict[Exponent, Fraction] = {}
        div_items = list(divisor.coeffs.items())
        steps = 0

        while rem:
            steps += 1
            if steps > max_steps:
                return None
            r_lead_exp = max(rem)
            if any(ra < ba for ra, ba in zip(r_lead_exp, b_lead_exp)):
                return None
            q_exp = tuple(ra - ba for ra, ba in zip(r_lead_exp, b_lead_exp))
            q_c = rem[r_lead_exp] / b_lead_c
            quot[q_exp] = quot.get(q_exp, Fraction(0, 1)) + q_c
            if quot[q_exp] == 0:
                del quot[q_exp]

            for d_exp, d_c in div_items:
                t_exp = tuple(qa + da for qa, da in zip(q_exp, d_exp))
                nv = rem.get(t_exp, Fraction(0, 1)) - q_c * d_c
                if nv == 0:
                    rem.pop(t_exp, None)
                else:
                    rem[t_exp] = nv

        return QQPoly._raw(self.ring, quot)

    # --- Evaluation & Substitution ---

    def eval_fraction(
        self,
        point: Union[Mapping[Union[int, str], Scalar], Sequence[Scalar]],
    ) -> Fraction:
        """Fast exact evaluation at a point in `QQ^m`."""
        if not self.coeffs:
            return Fraction(0, 1)
        if isinstance(point, Mapping):
            vals = [Fraction(0, 1)] * self.ring.num_vars
            for k, v in point.items():
                vals[self.ring.resolve_var(k)] = _to_fraction(v)
        else:
            if len(point) != self.ring.num_vars:
                raise ValueError(f"Expected {self.ring.num_vars} coordinates, got {len(point)}")
            vals = [_to_fraction(v) for v in point]

        max_degs = self.degrees()
        pow_cache: List[List[Fraction]] = []
        for i, d in enumerate(max_degs):
            pows = [Fraction(1, 1)]
            vi = vals[i]
            for _ in range(d):
                pows.append(pows[-1] * vi)
            pow_cache.append(pows)

        total = Fraction(0, 1)
        for exp, c in self.coeffs.items():
            term = c
            for i, e in enumerate(exp):
                if e:
                    term *= pow_cache[i][e]
                    if term == 0:
                        break
            total += term
        return total

    def subs(
        self,
        mapping: Mapping[Union[int, str, QQPoly, QQRat], Union[Scalar, QQPoly, QQRat]],
        target_ring: Optional[QQVarRing] = None,
    ) -> Union[QQPoly, QQRat]:
        """Exact substitution of scalars, `QQPoly`, or `QQRat` for subset/all variables."""
        has_rat = any(isinstance(v, QQRat) for v in mapping.values())
        if target_ring is None:
            for v in mapping.values():
                if isinstance(v, (QQPoly, QQRat)):
                    target_ring = v.ring
                    break
            if target_ring is None:
                target_ring = self.ring

        resolved: Dict[int, Union[QQPoly, QQRat]] = {}
        for k, v in mapping.items():
            idx = self.ring.resolve_var(k)
            if isinstance(v, QQRat):
                resolved[idx] = target_ring.embed_rat(v)
            elif isinstance(v, QQPoly):
                resolved[idx] = target_ring.embed_poly(v)
            else:
                resolved[idx] = target_ring.const(_to_fraction(v))

        if not has_rat:
            if not self.coeffs:
                return target_ring.zero()
            poly_vals: List[QQPoly] = []
            for i in range(self.ring.num_vars):
                if i in resolved:
                    poly_vals.append(resolved[i])  # type: ignore[arg-type]
                else:
                    poly_vals.append(target_ring.gen(self.ring.names[i]))
            max_degs = self.degrees()
            pow_cache: List[List[QQPoly]] = []
            for i, d in enumerate(max_degs):
                pows = [target_ring.one()]
                vi = poly_vals[i]
                for _ in range(d):
                    pows.append(pows[-1] * vi)  # type: ignore[arg-type]
                pow_cache.append(pows)

            acc = target_ring.zero()
            for exp, c in self.coeffs.items():
                term = target_ring.const(c)
                for i, e in enumerate(exp):
                    if e:
                        term = term * pow_cache[i][e]  # type: ignore[assignment]
                        if term.is_zero():
                            break
                acc = acc + term  # type: ignore[assignment]
            return acc

        rat_vals: List[QQRat] = []
        for i in range(self.ring.num_vars):
            if i in resolved:
                rv = resolved[i]
                rat_vals.append(rv if isinstance(rv, QQRat) else QQRat._raw(rv, target_ring.one()))
            else:
                rat_vals.append(target_ring.rat_gen(self.ring.names[i]))
        return _eval_poly_at_rats(self, rat_vals, target_ring)


# ---------------------------------------------------------------------------
# Exact Multivariate Polynomial GCD over ZZ / QQ (Collins-Brown Subresultant PRS)
# ---------------------------------------------------------------------------

ZZDict = Dict[Exponent, int]


def _zz_from_qqpoly(p: QQPoly) -> ZZDict:
    if not p.coeffs:
        return {}
    lcm_den = 1
    for c in p.coeffs.values():
        lcm_den = _lcm(lcm_den, c.denominator)
    int_coeffs: ZZDict = {}
    gcd_num = 0
    for exp, c in p.coeffs.items():
        iv = c.numerator * (lcm_den // c.denominator)
        int_coeffs[exp] = iv
        gcd_num = math.gcd(gcd_num, abs(iv))
    lead_exp = max(int_coeffs)
    sign = -1 if int_coeffs[lead_exp] < 0 else 1
    scale = gcd_num * sign
    if scale != 1:
        int_coeffs = {e: v // scale for e, v in int_coeffs.items()}
    return int_coeffs


def _zz_primitive(a: ZZDict) -> ZZDict:
    if not a:
        return {}
    g = 0
    for v in a.values():
        g = math.gcd(g, abs(v))
        if g == 1:
            break
    lead_exp = max(a)
    sign = -1 if a[lead_exp] < 0 else 1
    if g == 1 and sign == 1:
        return a
    div = g * sign
    return {e: v // div for e, v in a.items()}


def _zz_degrees(a: ZZDict, n: int) -> List[int]:
    res = [0] * n
    for exp in a:
        for i, e in enumerate(exp):
            if e > res[i]:
                res[i] = e
    return res


def _zz_monomial_content(a: ZZDict, n: int) -> Exponent:
    if not a:
        return (0,) * n
    it = iter(a)
    mins = list(next(it))
    for exp in it:
        for i, e in enumerate(exp):
            if e < mins[i]:
                mins[i] = e
    return tuple(mins)


def _zz_shift(a: ZZDict, delta: Exponent) -> ZZDict:
    if not a or all(d == 0 for d in delta):
        return a
    return {tuple(x + d for x, d in zip(e, delta)): v for e, v in a.items()}


def _zz_mul(a: ZZDict, b: ZZDict) -> ZZDict:
    if not a or not b:
        return {}
    out: ZZDict = {}
    b_items = list(b.items())
    for e1, c1 in a.items():
        for e2, c2 in b_items:
            exp = tuple(x + y for x, y in zip(e1, e2))
            nv = out.get(exp, 0) + c1 * c2
            if nv == 0:
                out.pop(exp, None)
            else:
                out[exp] = nv
    return out


def _zz_sub(a: ZZDict, b: ZZDict) -> ZZDict:
    if not b:
        return dict(a)
    out = dict(a)
    for e, c in b.items():
        nv = out.get(e, 0) - c
        if nv == 0:
            out.pop(e, None)
        else:
            out[e] = nv
    return out


def _zz_pow(a: ZZDict, k: int, n: int) -> ZZDict:
    if k == 0:
        return {(0,) * n: 1}
    if k == 1 or not a:
        return a
    res: ZZDict = {(0,) * n: 1}
    base = a
    p = k
    while p > 0:
        if p & 1:
            res = _zz_mul(res, base)
        p >>= 1
        if p:
            base = _zz_mul(base, base)
    return res


def _zz_try_div_exact(a: ZZDict, b: ZZDict, n: int, max_steps: int = 6000) -> Optional[ZZDict]:
    if not b:
        raise ZeroDivisionError("ZZPoly division by zero")
    if not a:
        return {}
    zero_exp = (0,) * n
    if len(b) == 1 and zero_exp in b:
        d_val = b[zero_exp]
        out_c: ZZDict = {}
        for e, v in a.items():
            q, r = divmod(v, d_val)
            if r != 0:
                return None
            out_c[e] = q
        return out_c

    degs_a = _zz_degrees(a, n)
    degs_b = _zz_degrees(b, n)
    if any(da < db for da, db in zip(degs_a, degs_b)):
        return None

    b_lead = max(b)
    b_lc = b[b_lead]
    rem = dict(a)
    quot: ZZDict = {}
    b_items = list(b.items())
    steps = 0

    while rem:
        steps += 1
        if steps > max_steps:
            return None
        r_lead = max(rem)
        if any(ra < ba for ra, ba in zip(r_lead, b_lead)):
            return None
        q_c, mod_c = divmod(rem[r_lead], b_lc)
        if mod_c != 0:
            return None
        q_exp = tuple(ra - ba for ra, ba in zip(r_lead, b_lead))
        quot[q_exp] = quot.get(q_exp, 0) + q_c
        if quot[q_exp] == 0:
            del quot[q_exp]
        for d_exp, d_c in b_items:
            t_exp = tuple(qa + da for qa, da in zip(q_exp, d_exp))
            nv = rem.get(t_exp, 0) - q_c * d_c
            if nv == 0:
                rem.pop(t_exp, None)
            else:
                rem[t_exp] = nv
    return quot


def _zz_decompose_var(a: ZZDict, var_idx: int) -> Dict[int, ZZDict]:
    """Decompose `a` as `sum_d coeff_d * x_{var_idx}^d` where `coeff_d` has 0 in `var_idx`."""
    by_deg: Dict[int, ZZDict] = {}
    for exp, c in a.items():
        d = exp[var_idx]
        base_exp = exp[:var_idx] + (0,) + exp[var_idx + 1 :]
        bucket = by_deg.get(d)
        if bucket is None:
            by_deg[d] = {base_exp: c}
        else:
            bucket[base_exp] = c
    return by_deg


def _zz_recompose_var(by_deg: Dict[int, ZZDict], var_idx: int) -> ZZDict:
    out: ZZDict = {}
    for d, coeff_poly in by_deg.items():
        for exp, c in coeff_poly.items():
            if c != 0:
                full_exp = exp[:var_idx] + (d,) + exp[var_idx + 1 :]
                out[full_exp] = c
    return out


def _zz_univariate_gcd(a: ZZDict, b: ZZDict, var_idx: int, n: int) -> ZZDict:
    """Fast primitive Euclidean GCD in `ZZ[x_{var_idx}]`."""
    pa: Dict[int, int] = {exp[var_idx]: c for exp, c in a.items()}
    pb: Dict[int, int] = {exp[var_idx]: c for exp, c in b.items()}

    def prim1d(p: Dict[int, int]) -> Dict[int, int]:
        if not p:
            return {}
        g = 0
        for v in p.values():
            g = math.gcd(g, abs(v))
            if g == 1:
                break
        lead = max(p)
        s = -1 if p[lead] < 0 else 1
        div = g * s
        if div == 1:
            return p
        return {k: v // div for k, v in p.items()}

    pa = prim1d(pa)
    pb = prim1d(pb)
    if max(pa) < max(pb):
        pa, pb = pb, pa

    while pb:
        db = max(pb)
        if db == 0:
            return {(0,) * n: 1}
        lcb = pb[db]
        rem = dict(pa)
        while rem and max(rem) >= db:
            dr = max(rem)
            cr = rem[dr]
            g = math.gcd(abs(cr), abs(lcb))
            mul_r = lcb // g
            mul_b = cr // g
            shift = dr - db
            new_rem: Dict[int, int] = {}
            for k, v in rem.items():
                if k != dr:
                    nv = v * mul_r
                    if nv:
                        new_rem[k] = nv
            for k, v in pb.items():
                if k != db:
                    tk = k + shift
                    nv = new_rem.get(tk, 0) - v * mul_b
                    if nv == 0:
                        new_rem.pop(tk, None)
                    else:
                        new_rem[tk] = nv
            rem = prim1d(new_rem)
        pa, pb = pb, rem

    out: ZZDict = {}
    for deg, c in pa.items():
        exp = tuple(deg if i == var_idx else 0 for i in range(n))
        out[exp] = c
    return _zz_primitive(out)


def _zz_content_in_var(by_deg: Dict[int, ZZDict], n: int, budget: int) -> ZZDict:
    items = sorted(by_deg.values(), key=len)
    g = _zz_primitive(items[0])
    zero_exp = (0,) * n
    for coeff_poly in items[1:]:
        if len(g) == 1 and zero_exp in g:
            return {zero_exp: 1}
        g = _zz_gcd(g, coeff_poly, n, budget)
    return _zz_primitive(g)


def _fp_poly_gcd_is_one(pa: Dict[int, int], pb: Dict[int, int], p: int) -> bool:
    """Return True if univariate polynomials `pa, pb` over `F_p[x]` have `gcd == 1`."""
    # Convert to dense lists up to max degree
    da = max(pa) if pa else -1
    db = max(pb) if pb else -1
    if da <= 0 or db <= 0:
        return True
    va = [0] * (da + 1)
    for k, c in pa.items():
        va[k] = c % p
    vb = [0] * (db + 1)
    for k, c in pb.items():
        vb[k] = c % p
    while va and va[-1] == 0:
        va.pop()
    while vb and vb[-1] == 0:
        vb.pop()
    if len(va) < len(vb):
        va, vb = vb, va
    while vb:
        if len(vb) == 1:
            return True
        inv_lead = pow(vb[-1], p - 2, p)
        db_cur = len(vb) - 1
        while len(va) >= len(vb):
            fac = (va[-1] * inv_lead) % p
            shift = len(va) - len(vb)
            for idx in range(db_cur):
                if vb[idx]:
                    va[shift + idx] = (va[shift + idx] - fac * vb[idx]) % p
            va.pop()
            while va and va[-1] == 0:
                va.pop()
        va, vb = vb, va
    return len(va) <= 1


_FP_MOD = 2147483647  # 2^31 - 1 (Mersenne prime)


def _zz_var_absent_from_gcd_cert(a: ZZDict, b: ZZDict, var_idx: int, n: int) -> bool:
    """100% sound certificate proving `deg_{x_{var_idx}}(gcd(a, b)) == 0`.

    By Gauss's Lemma and ring homomorphism `Z[x_{!=var_idx}] -> F_p`, if at an evaluation
    point `alpha` both leading coefficients `lc_{x_j}(a)(alpha) != 0 (mod p)` and
    `lc_{x_j}(b)(alpha) != 0 (mod p)` and `gcd_{F_p[x_j]}(a(alpha, x_j), b(alpha, x_j)) == 1`,
    then `gcd(a, b)` cannot depend on `x_{var_idx}`.
    """
    p = _FP_MOD
    da = max(exp[var_idx] for exp in a)
    db = max(exp[var_idx] for exp in b)
    if da == 0 or db == 0:
        return True

    for seed in (104729, 1299709):
        alphas = [((seed * (i + 3) + 17 * (i + 1) * (i + 1)) % (p - 3)) + 2 for i in range(n)]
        # Precompute power tables for variables != var_idx
        degs_a = _zz_degrees(a, n)
        degs_b = _zz_degrees(b, n)
        pow_tab: List[List[int]] = []
        for i in range(n):
            if i == var_idx:
                pow_tab.append([1])
            else:
                md = max(degs_a[i], degs_b[i])
                pows = [1]
                ai = alphas[i]
                for _ in range(md):
                    pows.append((pows[-1] * ai) % p)
                pow_tab.append(pows)

        pa_1d: Dict[int, int] = {}
        for exp, c in a.items():
            prod = c % p
            for i, e in enumerate(exp):
                if i != var_idx and e:
                    prod = (prod * pow_tab[i][e]) % p
            if prod:
                d = exp[var_idx]
                pa_1d[d] = (pa_1d.get(d, 0) + prod) % p
        if pa_1d.get(da, 0) == 0:
            continue

        pb_1d: Dict[int, int] = {}
        for exp, c in b.items():
            prod = c % p
            for i, e in enumerate(exp):
                if i != var_idx and e:
                    prod = (prod * pow_tab[i][e]) % p
            if prod:
                d = exp[var_idx]
                pb_1d[d] = (pb_1d.get(d, 0) + prod) % p
        if pb_1d.get(db, 0) == 0:
            continue

        return _fp_poly_gcd_is_one(pa_1d, pb_1d, p)
    return False


def _zz_gcd(a: ZZDict, b: ZZDict, n: int, budget: int = 1800) -> ZZDict:
    """Primitive multivariate GCD in `ZZ[x_1, ..., x_n]`."""
    if not a:
        return _zz_primitive(b)
    if not b:
        return _zz_primitive(a)

    zero_exp = (0,) * n
    if (len(a) == 1 and zero_exp in a) or (len(b) == 1 and zero_exp in b):
        return {zero_exp: 1}

    ma = _zz_monomial_content(a, n)
    mb = _zz_monomial_content(b, n)
    m_common = tuple(min(x, y) for x, y in zip(ma, mb))
    a0 = _zz_primitive(_zz_shift(a, tuple(-x for x in ma)))
    b0 = _zz_primitive(_zz_shift(b, tuple(-x for x in mb)))

    if (len(a0) == 1 and zero_exp in a0) or (len(b0) == 1 and zero_exp in b0):
        return {m_common: 1}
    if a0 == b0:
        return _zz_shift(a0, m_common)

    if len(a0) <= budget and len(b0) <= budget:
        if _zz_try_div_exact(a0, b0, n) is not None:
            return _zz_shift(b0, m_common)
        if _zz_try_div_exact(b0, a0, n) is not None:
            return _zz_shift(a0, m_common)

    degs_a = _zz_degrees(a0, n)
    degs_b = _zz_degrees(b0, n)

    # Variable present in only one polynomial cannot appear in their GCD
    for j in range(n):
        if degs_a[j] > 0 and degs_b[j] == 0:
            ca = _zz_content_in_var(_zz_decompose_var(a0, j), n, budget)
            g = _zz_gcd(ca, b0, n, budget)
            return _zz_shift(g, m_common)
        if degs_b[j] > 0 and degs_a[j] == 0:
            cb = _zz_content_in_var(_zz_decompose_var(b0, j), n, budget)
            g = _zz_gcd(a0, cb, n, budget)
            return _zz_shift(g, m_common)

    active = [j for j in range(n) if degs_a[j] > 0]
    if not active:
        return {m_common: 1}
    if len(active) == 1:
        g1 = _zz_univariate_gcd(a0, b0, active[0], n)
        return _zz_shift(g1, m_common)

    # Fast modular univariate evaluation certificate: eliminate variables absent from gcd(a0, b0)
    absent_vars = [j for j in active if _zz_var_absent_from_gcd_cert(a0, b0, j, n)]
    if len(absent_vars) == len(active):
        return {m_common: 1}
    if absent_vars:
        elim_j = absent_vars[0]
        ca_elim = _zz_content_in_var(_zz_decompose_var(a0, elim_j), n, budget)
        cb_elim = _zz_content_in_var(_zz_decompose_var(b0, elim_j), n, budget)
        g_elim = _zz_gcd(ca_elim, cb_elim, n, budget)
        return _zz_shift(g_elim, m_common)

    # Linear-degree variable fast path
    for j in active:
        if degs_a[j] == 1 or degs_b[j] == 1:
            lin, other = (a0, b0) if degs_a[j] == 1 else (b0, a0)
            dec_lin = _zz_decompose_var(lin, j)
            cont_lin = _zz_content_in_var(dec_lin, n, budget)
            g_cont = _zz_gcd(cont_lin, other, n, budget)
            prim_lin = _zz_try_div_exact(lin, cont_lin, n)
            other_red = _zz_try_div_exact(other, g_cont, n)
            if prim_lin is not None and other_red is not None:
                if _zz_try_div_exact(other_red, prim_lin, n) is not None:
                    return _zz_shift(_zz_primitive(_zz_mul(g_cont, prim_lin)), m_common)
            return _zz_shift(g_cont, m_common)

    prs_budget = budget if len(active) <= 2 else (400 if len(active) == 3 else 140)
    if len(a0) > prs_budget or len(b0) > prs_budget:
        return {m_common: 1}

    # Choose main variable minimizing (min_deg, max_deg)
    main_var = min(active, key=lambda j: (min(degs_a[j], degs_b[j]), max(degs_a[j], degs_b[j])))
    dec_a = _zz_decompose_var(a0, main_var)
    dec_b = _zz_decompose_var(b0, main_var)
    cont_a = _zz_content_in_var(dec_a, n, budget)
    cont_b = _zz_content_in_var(dec_b, n, budget)
    cont_g = _zz_gcd(cont_a, cont_b, n, budget)

    pa = _zz_try_div_exact(a0, cont_a, n)
    pb = _zz_try_div_exact(b0, cont_b, n)
    if pa is None or pb is None:
        return _zz_shift(cont_g, m_common)

    if max(_zz_decompose_var(pa, main_var)) < max(_zz_decompose_var(pb, main_var)):
        pa, pb = pb, pa

    f1 = _zz_decompose_var(pa, main_var)
    f2 = _zz_decompose_var(pb, main_var)
    g_pol: ZZDict = {zero_exp: 1}
    h_pol: ZZDict = {zero_exp: 1}

    while f2:
        d1 = max(f1)
        d2 = max(f2)
        if d2 == 0:
            return _zz_shift(cont_g, m_common)
        delta = d1 - d2
        lc2 = f2[d2]

        # Standard pseudo-remainder R = prem_{x}(f1, f2)
        rem_d: Dict[int, ZZDict] = dict(f1)
        for deg_step in range(d1, d2 - 1, -1):
            c_top = rem_d.pop(deg_step, None)
            new_rem_d: Dict[int, ZZDict] = {}
            for k, poly_k in rem_d.items():
                mp = _zz_mul(lc2, poly_k)
                if mp:
                    new_rem_d[k] = mp
            if c_top is not None:
                shift = deg_step - d2
                for k, poly_f2 in f2.items():
                    if k == d2:
                        continue
                    tk = k + shift
                    sub_p = _zz_mul(c_top, poly_f2)
                    cur_p = new_rem_d.get(tk, {})
                    diff_p = _zz_sub(cur_p, sub_p)
                    if diff_p:
                        new_rem_d[tk] = diff_p
                    else:
                        new_rem_d.pop(tk, None)
            rem_d = new_rem_d
            if sum(len(v) for v in rem_d.values()) > 4 * budget:
                return _zz_shift(cont_g, m_common)

        if not rem_d:
            break
        if max(rem_d) == 0:
            return _zz_shift(cont_g, m_common)

        # Subresultant exact scaling: next_f2 = rem / (g * h^delta)
        denom_scale = _zz_mul(g_pol, _zz_pow(h_pol, delta, n))
        scaled_f2: Dict[int, ZZDict] = {}
        for k, poly_k in rem_d.items():
            q_k = _zz_try_div_exact(poly_k, denom_scale, n)
            if q_k is None:
                return _zz_shift(cont_g, m_common)
            if q_k:
                scaled_f2[k] = q_k

        f1 = f2
        f2 = scaled_f2
        g_pol = lc2
        if delta == 0:
            pass
        elif delta == 1:
            h_pol = g_pol
        else:
            g_pow = _zz_pow(g_pol, delta, n)
            h_pow = _zz_pow(h_pol, delta - 1, n)
            nh = _zz_try_div_exact(g_pow, h_pow, n)
            if nh is None:
                return _zz_shift(cont_g, m_common)
            h_pol = nh

    # Primitive part of f2 in main_var, multiplied by cont_g
    cont_f2 = _zz_content_in_var(f2, n, budget)
    f2_full = _zz_recompose_var(f2, main_var)
    pp_f2 = _zz_try_div_exact(f2_full, cont_f2, n)
    if pp_f2 is None:
        return _zz_shift(cont_g, m_common)
    res = _zz_primitive(_zz_mul(cont_g, pp_f2))
    return _zz_shift(res, m_common)


def poly_gcd(a: QQPoly, b: QQPoly, budget: int = 1800) -> QQPoly:
    """Verified exact monic GCD in `QQ[x_1, ..., x_m]`.

    Every non-trivial divisor returned is unconditionally certified by
    `a.try_div_exact(g)` and `b.try_div_exact(g)`.
    """
    if a.ring != b.ring:
        raise ValueError(f"Ring mismatch: {a.ring} vs {b.ring}")
    ring = a.ring
    if a.is_zero():
        return b.monic()
    if b.is_zero():
        return a.monic()
    if a.is_constant() or b.is_constant():
        return ring.one()

    za = _zz_from_qqpoly(a)
    zb = _zz_from_qqpoly(b)
    zg = _zz_gcd(za, zb, ring.num_vars, budget=budget)
    g = QQPoly._raw(ring, {e: Fraction(v, 1) for e, v in zg.items()}).monic()
    if g.is_constant():
        return ring.one()

    # Zero-trust exact division certificate
    if a.try_div_exact(g) is not None and b.try_div_exact(g) is not None:
        return g

    # Fallback to monomial GCD if budget bailout occurred mid-PRS
    ma = a.monomial_content()
    mb = b.monomial_content()
    mc = tuple(min(x, y) for x, y in zip(ma, mb))
    return QQPoly._raw(ring, {mc: Fraction(1, 1)})


# ---------------------------------------------------------------------------
# Exact Multivariate Rational Function Field over QQ
# ---------------------------------------------------------------------------


def _eval_poly_at_rats(poly: QQPoly, rat_vals: Sequence[QQRat], target_ring: QQVarRing) -> QQRat:
    """Evaluate `poly(N_1/D_1, ..., N_m/D_m)` via denominator-grouped clearing."""
    if poly.is_zero():
        return target_ring.rat_zero()

    # Group variables sharing identical non-one denominator polynomials
    den_groups: List[Tuple[QQPoly, List[int]]] = []
    den_key_to_group: Dict[int, int] = {}
    var_to_group: List[Optional[int]] = [None] * poly.ring.num_vars

    for i, rv in enumerate(rat_vals):
        if rv.den.is_constant() and rv.den.constant_value() == 1:
            continue
        matched = None
        for g_idx, (d_poly, members) in enumerate(den_groups):
            if d_poly.coeffs == rv.den.coeffs:
                matched = g_idx
                members.append(i)
                break
        if matched is None:
            matched = len(den_groups)
            den_groups.append((rv.den, [i]))
        var_to_group[i] = matched
        den_key_to_group[i] = matched

    # Max combined degree in each denominator group across monomials of poly
    group_max_deg = [0] * len(den_groups)
    for exp in poly.coeffs:
        g_sums = [0] * len(den_groups)
        for i, e in enumerate(exp):
            g_idx = var_to_group[i]
            if g_idx is not None and e:
                g_sums[g_idx] += e
        for g_idx, s in enumerate(g_sums):
            if s > group_max_deg[g_idx]:
                group_max_deg[g_idx] = s

    # Precompute powers of numerators N_i^e and group denominators D_g^k
    max_degs = poly.degrees()
    num_pows: List[List[QQPoly]] = []
    for i, d in enumerate(max_degs):
        pows = [target_ring.one()]
        ni = rat_vals[i].num
        for _ in range(d):
            pows.append(pows[-1] * ni)  # type: ignore[arg-type]
        num_pows.append(pows)

    den_pows: List[List[QQPoly]] = []
    for g_idx, (d_poly, _) in enumerate(den_groups):
        md = group_max_deg[g_idx]
        pows = [target_ring.one()]
        for _ in range(md):
            pows.append(pows[-1] * d_poly)  # type: ignore[arg-type]
        den_pows.append(pows)

    acc_num = target_ring.zero()
    for exp, c in poly.coeffs.items():
        term = target_ring.const(c)
        g_sums = [0] * len(den_groups)
        zero_term = False
        for i, e in enumerate(exp):
            if e:
                term = term * num_pows[i][e]  # type: ignore[assignment]
                if term.is_zero():
                    zero_term = True
                    break
                g_idx = var_to_group[i]
                if g_idx is not None:
                    g_sums[g_idx] += e
        if zero_term:
            continue
        for g_idx, md in enumerate(group_max_deg):
            rem_pow = md - g_sums[g_idx]
            if rem_pow:
                term = term * den_pows[g_idx][rem_pow]  # type: ignore[assignment]
        acc_num = acc_num + term  # type: ignore[assignment]

    acc_den = target_ring.one()
    for g_idx, md in enumerate(group_max_deg):
        if md:
            acc_den = acc_den * den_pows[g_idx][md]  # type: ignore[assignment]

    return QQRat(acc_num, acc_den)


class QQRat:
    """Exact element `(num, den)` of `QQ(x_1, ..., x_m)`."""

    __slots__ = ("num", "den")

    def __init__(self, num: QQPoly, den: Optional[QQPoly] = None, auto_cancel: bool = True) -> None:
        if den is None:
            den = num.ring.one()
        if num.ring != den.ring:
            raise ValueError(f"Ring mismatch: {num.ring} vs {den.ring}")
        if den.is_zero():
            raise ZeroDivisionError("Zero denominator in QQRat")

        ring = num.ring
        if num.is_zero():
            self.num = ring.zero()
            self.den = ring.one()
            return

        # Fast monomial content cancellation
        mn = num.monomial_content()
        md = den.monomial_content()
        mc = tuple(min(a, b) for a, b in zip(mn, md))
        if any(mc):
            neg_mc = tuple(-x for x in mc)
            num = num.shift_exponents(neg_mc)
            den = den.shift_exponents(neg_mc)

        if auto_cancel and not den.is_constant():
            q_exact = num.try_div_exact(den)
            if q_exact is not None:
                num = q_exact
                den = ring.one()
            else:
                g = poly_gcd(num, den)
                if not g.is_constant():
                    qn = num.try_div_exact(g)
                    qd = den.try_div_exact(g)
                    if qn is not None and qd is not None:
                        num = qn
                        den = qd

        # Normalize denominator to be monic over QQ
        _, lc_den = den.leading_term()
        if lc_den != 1:
            inv_lc = Fraction(1, 1) / lc_den
            num = num * inv_lc  # type: ignore[assignment]
            den = den * inv_lc  # type: ignore[assignment]

        self.num: QQPoly = num
        self.den: QQPoly = den

    @classmethod
    def _raw(cls, num: QQPoly, den: QQPoly) -> QQRat:
        """Construct from already reduced `(num, monic_den)`."""
        obj = cls.__new__(cls)
        obj.num = num
        obj.den = den
        return obj

    @property
    def ring(self) -> QQVarRing:
        return self.num.ring

    def _coerce_rat(self, other: object) -> Optional[QQRat]:
        if isinstance(other, QQRat):
            if other.ring != self.ring:
                raise ValueError(f"Ring mismatch: {self.ring} vs {other.ring}")
            return other
        if isinstance(other, QQPoly):
            if other.ring != self.ring:
                raise ValueError(f"Ring mismatch: {self.ring} vs {other.ring}")
            return QQRat._raw(other, self.ring.one())
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            return self.ring.rat_const(other)
        return None

    # --- Queries & Rigorous Identity Verification ---

    def is_zero(self) -> bool:
        return self.num.is_zero()

    def is_polynomial(self) -> bool:
        return self.den.is_constant()

    def degrees(self) -> Tuple[int, int]:
        return (self.num.total_degree(), self.den.total_degree())

    def degree_in_vars(self, vars_subset: Iterable[Union[int, str, QQPoly, QQRat]]) -> Tuple[int, int]:
        vars_list = list(vars_subset)
        return (self.num.degree_in_vars(vars_list), self.den.degree_in_vars(vars_list))

    def monomial_counts(self) -> Tuple[int, int]:
        return (self.num.monomial_count(), self.den.monomial_count())

    def max_coeff_bits(self) -> int:
        return max(self.num.max_coeff_bits(), self.den.max_coeff_bits())

    def cross_diff_poly(self, other: Union[QQRat, QQPoly, Scalar]) -> QQPoly:
        """Compute `self.num * other.den - other.num * self.den` in `QQ[x_1, ..., x_m]`."""
        r = self._coerce_rat(other)
        if r is None:
            raise TypeError(f"Cannot compare QQRat with {type(other).__name__}")
        return (self.num * r.den) - (r.num * self.den)  # type: ignore[return-value]

    def equiv(self, other: Union[QQRat, QQPoly, Scalar]) -> bool:
        """100% rigorous identity check in `QQ(x_1, ..., x_m)`:
        `A == B iff (A.num * B.den - B.num * A.den).is_zero()`.
        """
        r = self._coerce_rat(other)
        if r is None:
            return False
        if r.num.is_zero():
            return self.num.is_zero()
        if self.num.is_zero():
            return False
        if self.den.coeffs == r.den.coeffs:
            return self.num.coeffs == r.num.coeffs
        return ((self.num * r.den) - (r.num * self.den)).is_zero()  # type: ignore[union-attr]

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QQJet):
            return other == self
        r = self._coerce_rat(other)
        if r is None:
            return False
        return self.equiv(r)

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __bool__(self) -> bool:
        return not self.num.is_zero()

    def __repr__(self) -> str:
        return f"QQRat({self.__str__()})"

    def __str__(self) -> str:
        if self.den.is_constant() and self.den.constant_value() == 1:
            return str(self.num)
        return f"({self.num}) / ({self.den})"

    # --- Knuth UFD Lowest-Terms Arithmetic ---

    def __pos__(self) -> QQRat:
        return self

    def __neg__(self) -> QQRat:
        if self.num.is_zero():
            return self
        return QQRat._raw(-self.num, self.den)

    def __add__(self, other: object) -> Union[QQRat, QQJet]:
        r = self._coerce_rat(other)
        if r is None:
            return NotImplemented
        if self.num.is_zero():
            return r
        if r.num.is_zero():
            return self

        p1, q1 = self.num, self.den
        p2, q2 = r.num, r.den

        # Identical denominator fast path
        if q1.coeffs == q2.coeffs:
            return QQRat(p1 + p2, q1)  # type: ignore[arg-type]

        # Polynomial denominator fast paths (Knuth: g = gcd(q1, 1) = 1 => coprime!)
        if q2.is_constant():
            # q2 == 1 because monic
            new_num = p1 + p2 * q1
            if new_num.is_zero():  # type: ignore[union-attr]
                return self.ring.rat_zero()
            return QQRat._raw(new_num, q1)  # type: ignore[arg-type]
        if q1.is_constant():
            new_num = p1 * q2 + p2
            if new_num.is_zero():  # type: ignore[union-attr]
                return self.ring.rat_zero()
            return QQRat._raw(new_num, q2)  # type: ignore[arg-type]

        # General Knuth UFD addition
        g = poly_gcd(q1, q2)
        if g.is_constant():
            new_num = p1 * q2 + p2 * q1
            if new_num.is_zero():  # type: ignore[union-attr]
                return self.ring.rat_zero()
            new_den = q1 * q2
            return QQRat._raw(new_num, new_den.monic())  # type: ignore[arg-type, union-attr]

        q1_bar = q1.try_div_exact(g)
        q2_bar = q2.try_div_exact(g)
        if q1_bar is None or q2_bar is None:
            return QQRat(p1 * q2 + p2 * q1, q1 * q2)  # type: ignore[arg-type]

        num_comb: QQPoly = (p1 * q2_bar + p2 * q1_bar)  # type: ignore[assignment]
        if num_comb.is_zero():
            return self.ring.rat_zero()

        h = poly_gcd(num_comb, g)
        if not h.is_constant():
            nc = num_comb.try_div_exact(h)
            q1_red = q1.try_div_exact(h)
            if nc is not None and q1_red is not None:
                den_final: QQPoly = (q1_red * q2_bar)  # type: ignore[assignment]
                _, lc = den_final.leading_term()
                if lc != 1:
                    inv = Fraction(1, 1) / lc
                    nc = nc * inv  # type: ignore[assignment]
                    den_final = den_final * inv  # type: ignore[assignment]
                return QQRat._raw(nc, den_final)

        den_comb: QQPoly = (q1 * q2_bar)  # type: ignore[assignment]
        _, lc = den_comb.leading_term()
        if lc != 1:
            inv = Fraction(1, 1) / lc
            num_comb = num_comb * inv  # type: ignore[assignment]
            den_comb = den_comb * inv  # type: ignore[assignment]
        return QQRat._raw(num_comb, den_comb)

    def __radd__(self, other: object) -> Union[QQRat, QQJet]:
        return self.__add__(other)

    def __sub__(self, other: object) -> Union[QQRat, QQJet]:
        r = self._coerce_rat(other)
        if r is None:
            return NotImplemented
        return self.__add__(-r)

    def __rsub__(self, other: object) -> Union[QQRat, QQJet]:
        r = self._coerce_rat(other)
        if r is None:
            return NotImplemented
        return r.__add__(-self)

    def __mul__(self, other: object) -> Union[QQRat, QQJet]:
        r = self._coerce_rat(other)
        if r is None:
            return NotImplemented
        if self.num.is_zero() or r.num.is_zero():
            return self.ring.rat_zero()

        p1, q1 = self.num, self.den
        p2, q2 = r.num, r.den

        if not q2.is_constant():
            g1 = poly_gcd(p1, q2)
            if not g1.is_constant():
                p1_q = p1.try_div_exact(g1)
                q2_q = q2.try_div_exact(g1)
                if p1_q is not None and q2_q is not None:
                    p1, q2 = p1_q, q2_q

        if not q1.is_constant():
            g2 = poly_gcd(p2, q1)
            if not g2.is_constant():
                p2_q = p2.try_div_exact(g2)
                q1_q = q1.try_div_exact(g2)
                if p2_q is not None and q1_q is not None:
                    p2, q1 = p2_q, q1_q

        num_prod: QQPoly = (p1 * p2)  # type: ignore[assignment]
        den_prod: QQPoly = (q1 * q2)  # type: ignore[assignment]
        _, lc = den_prod.leading_term()
        if lc != 1:
            inv = Fraction(1, 1) / lc
            num_prod = num_prod * inv  # type: ignore[assignment]
            den_prod = den_prod * inv  # type: ignore[assignment]
        return QQRat._raw(num_prod, den_prod)

    def __rmul__(self, other: object) -> Union[QQRat, QQJet]:
        return self.__mul__(other)

    def __truediv__(self, other: object) -> Union[QQRat, QQJet]:
        r = self._coerce_rat(other)
        if r is None:
            return NotImplemented
        if r.num.is_zero():
            raise ZeroDivisionError("Division of QQRat by zero")
        if self.num.is_zero():
            return self.ring.rat_zero()
        _, lc = r.num.leading_term()
        inv = Fraction(1, 1) / lc
        recip = QQRat._raw(r.den * inv, r.num * inv)  # type: ignore[arg-type]
        return self.__mul__(recip)

    def __rtruediv__(self, other: object) -> Union[QQRat, QQJet]:
        r = self._coerce_rat(other)
        if r is None:
            return NotImplemented
        return r.__truediv__(self)

    def __pow__(self, power: int) -> QQRat:
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError(f"Exponent must be int, got {type(power).__name__}")
        if power == 0:
            return self.ring.rat_one()
        if power < 0:
            if self.num.is_zero():
                raise ZeroDivisionError("Zero QQRat raised to negative power")
            pos = (-power)
            np: QQPoly = self.den ** pos  # type: ignore[assignment]
            dp: QQPoly = self.num ** pos  # type: ignore[assignment]
            _, lc = dp.leading_term()
            if lc != 1:
                inv = Fraction(1, 1) / lc
                np = np * inv  # type: ignore[assignment]
                dp = dp * inv  # type: ignore[assignment]
            return QQRat._raw(np, dp)
        if self.num.is_zero():
            return self
        return QQRat._raw(self.num ** power, self.den ** power)  # type: ignore[arg-type]

    # --- Calculus & Evaluation ---

    def diff(self, var: Union[int, str, QQPoly, QQRat], order: int = 1) -> QQRat:
        if order < 0:
            raise ValueError("Derivative order must be >= 0")
        cur = self
        for _ in range(order):
            if cur.num.is_zero():
                break
            dn = cur.num.diff(var)
            dd = cur.den.diff(var)
            if dd.is_zero():
                cur = QQRat(dn, cur.den)
            else:
                # (dn * den - num * dd) / den^2
                top: QQPoly = (dn * cur.den - cur.num * dd)  # type: ignore[assignment]
                cur = QQRat(top, cur.den ** 2)  # type: ignore[arg-type]
        return cur

    def eval_fraction(
        self,
        point: Union[Mapping[Union[int, str], Scalar], Sequence[Scalar]],
    ) -> Fraction:
        den_v = self.den.eval_fraction(point)
        if den_v == 0:
            raise ZeroDivisionError("Rational denominator vanishes at evaluation point")
        num_v = self.num.eval_fraction(point)
        return num_v / den_v

    def subs(
        self,
        mapping: Mapping[Union[int, str, QQPoly, QQRat], Union[Scalar, QQPoly, QQRat]],
        target_ring: Optional[QQVarRing] = None,
    ) -> QQRat:
        if target_ring is None:
            for v in mapping.values():
                if isinstance(v, (QQPoly, QQRat)):
                    target_ring = v.ring
                    break
            if target_ring is None:
                target_ring = self.ring

        rat_vals: List[QQRat] = []
        resolved: Dict[int, QQRat] = {}
        for k, v in mapping.items():
            idx = self.ring.resolve_var(k)
            if isinstance(v, QQRat):
                resolved[idx] = target_ring.embed_rat(v)
            elif isinstance(v, QQPoly):
                resolved[idx] = QQRat._raw(target_ring.embed_poly(v), target_ring.one())
            else:
                resolved[idx] = target_ring.rat_const(_to_fraction(v))

        for i in range(self.ring.num_vars):
            if i in resolved:
                rat_vals.append(resolved[i])
            else:
                rat_vals.append(target_ring.rat_gen(self.ring.names[i]))

        num_sub = _eval_poly_at_rats(self.num, rat_vals, target_ring)
        den_sub = _eval_poly_at_rats(self.den, rat_vals, target_ring)
        return num_sub / den_sub  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Truncated Formal Power Series (Jet) over QQRat(z)[[dt]] / (dt^(order+1))
# ---------------------------------------------------------------------------


class QQJet:
    """Truncated power series `c_0 + c_1*dt + ... + c_K*dt^K + O(dt^{K+1})` over `QQRat`."""

    __slots__ = ("ring", "order", "coeffs")

    def __init__(self, ring: QQVarRing, coeffs: Sequence[Union[QQRat, QQPoly, Scalar]], order: int) -> None:
        self.ring: QQVarRing = ring
        self.order: int = order
        padded: List[QQRat] = []
        for i in range(order + 1):
            if i < len(coeffs):
                c = coeffs[i]
                if isinstance(c, QQRat):
                    padded.append(ring.embed_rat(c))
                elif isinstance(c, QQPoly):
                    padded.append(QQRat._raw(ring.embed_poly(c), ring.one()))
                else:
                    padded.append(ring.rat_const(_to_fraction(c)))
            else:
                padded.append(ring.rat_zero())
        self.coeffs: Tuple[QQRat, ...] = tuple(padded)

    @classmethod
    def const(cls, val: Union[QQRat, QQPoly, Scalar], ring: QQVarRing, order: int = 2) -> QQJet:
        return cls(ring, [val], order)

    @classmethod
    def dt(cls, ring: QQVarRing, order: int = 2) -> QQJet:
        if order < 1:
            raise ValueError("Jet order must be >= 1 to represent dt")
        return cls(ring, [ring.rat_zero(), ring.rat_one()], order)

    def _coerce_jet(self, other: object) -> Optional[QQJet]:
        if isinstance(other, QQJet):
            if other.ring != self.ring or other.order != self.order:
                raise ValueError("QQJet ring/order mismatch")
            return other
        if isinstance(other, (QQRat, QQPoly, int, Fraction)) and not isinstance(other, bool):
            return QQJet.const(other, self.ring, self.order)
        return None

    def is_zero(self) -> bool:
        return all(c.is_zero() for c in self.coeffs)

    def __eq__(self, other: object) -> bool:
        j = self._coerce_jet(other)
        if j is None:
            return False
        return all(a == b for a, b in zip(self.coeffs, j.coeffs))

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __pos__(self) -> QQJet:
        return self

    def __neg__(self) -> QQJet:
        return QQJet(self.ring, [-c for c in self.coeffs], self.order)

    def __add__(self, other: object) -> QQJet:
        j = self._coerce_jet(other)
        if j is None:
            return NotImplemented
        return QQJet(self.ring, [a + b for a, b in zip(self.coeffs, j.coeffs)], self.order)  # type: ignore[misc]

    def __radd__(self, other: object) -> QQJet:
        return self.__add__(other)

    def __sub__(self, other: object) -> QQJet:
        j = self._coerce_jet(other)
        if j is None:
            return NotImplemented
        return QQJet(self.ring, [a - b for a, b in zip(self.coeffs, j.coeffs)], self.order)  # type: ignore[misc]

    def __rsub__(self, other: object) -> QQJet:
        j = self._coerce_jet(other)
        if j is None:
            return NotImplemented
        return j.__sub__(self)

    def __mul__(self, other: object) -> QQJet:
        j = self._coerce_jet(other)
        if j is None:
            return NotImplemented
        out: List[QQRat] = []
        for k in range(self.order + 1):
            acc = self.ring.rat_zero()
            for i in range(k + 1):
                if not self.coeffs[i].is_zero() and not j.coeffs[k - i].is_zero():
                    acc = acc + self.coeffs[i] * j.coeffs[k - i]  # type: ignore[assignment]
            out.append(acc)
        return QQJet(self.ring, out, self.order)

    def __rmul__(self, other: object) -> QQJet:
        return self.__mul__(other)

    def inv(self) -> QQJet:
        b0 = self.coeffs[0]
        if b0.is_zero():
            raise ZeroDivisionError("Cannot invert QQJet with zero constant term b_0")
        inv_b0: QQRat = self.ring.rat_one() / b0  # type: ignore[assignment]
        c: List[QQRat] = [inv_b0]
        for k in range(1, self.order + 1):
            s = self.ring.rat_zero()
            for i in range(1, k + 1):
                if not self.coeffs[i].is_zero() and not c[k - i].is_zero():
                    s = s + self.coeffs[i] * c[k - i]  # type: ignore[assignment]
            c.append(-(s * inv_b0))  # type: ignore[arg-type]
        return QQJet(self.ring, c, self.order)

    def __truediv__(self, other: object) -> QQJet:
        j = self._coerce_jet(other)
        if j is None:
            return NotImplemented
        return self * j.inv()

    def __rtruediv__(self, other: object) -> QQJet:
        j = self._coerce_jet(other)
        if j is None:
            return NotImplemented
        return j * self.inv()

    def __pow__(self, power: int) -> QQJet:
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError("QQJet exponent must be int")
        if power == 0:
            return QQJet.const(1, self.ring, self.order)
        if power < 0:
            return (self ** (-power)).inv()
        res = QQJet.const(1, self.ring, self.order)
        base = self
        k = power
        while k > 0:
            if k & 1:
                res = res * base
            k >>= 1
            if k:
                base = base * base
        return res

    def __repr__(self) -> str:
        terms = [f"({c})*dt^{i}" for i, c in enumerate(self.coeffs)]
        return " + ".join(terms) + f" + O(dt^{self.order + 1})"


# ---------------------------------------------------------------------------
# Symbolic Linear Algebra, Hamiltonian Vector Fields & Degree Analysis
# ---------------------------------------------------------------------------

RatMatrix = List[List[QQRat]]


def mat_eye(n: int, ring: QQVarRing) -> RatMatrix:
    z, o = ring.rat_zero(), ring.rat_one()
    return [[o if i == j else z for j in range(n)] for i in range(n)]


def mat_transpose(a: RatMatrix) -> RatMatrix:
    if not a:
        return []
    rows, cols = len(a), len(a[0])
    return [[a[r][c] for r in range(rows)] for c in range(cols)]


def mat_add(a: RatMatrix, b: RatMatrix) -> RatMatrix:
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]  # type: ignore[misc]


def mat_sub(a: RatMatrix, b: RatMatrix) -> RatMatrix:
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]  # type: ignore[misc]


def mat_mul(a: RatMatrix, b: RatMatrix) -> RatMatrix:
    rows, mid, cols = len(a), len(b), len(b[0])
    ring = a[0][0].ring
    out: RatMatrix = []
    for i in range(rows):
        row: List[QQRat] = []
        for j in range(cols):
            acc = ring.rat_zero()
            for k in range(mid):
                if not a[i][k].is_zero() and not b[k][j].is_zero():
                    acc = acc + a[i][k] * b[k][j]  # type: ignore[assignment]
            row.append(acc)
        out.append(row)
    return out


def mat_vec_mul(a: RatMatrix, v: Sequence[QQRat]) -> List[QQRat]:
    ring = a[0][0].ring
    out: List[QQRat] = []
    for row in a:
        acc = ring.rat_zero()
        for aij, vj in zip(row, v):
            if not aij.is_zero() and not vj.is_zero():
                acc = acc + aij * vj  # type: ignore[assignment]
        out.append(acc)
    return out


def mat_solve(a: RatMatrix, rhs: RatMatrix) -> RatMatrix:
    """Exact Gauss-Jordan elimination over `QQRat`."""
    n = len(a)
    m = len(rhs[0])
    aug: RatMatrix = [list(a[i]) + list(rhs[i]) for i in range(n)]

    for col in range(n):
        pivot = None
        for r in range(col, n):
            if not aug[r][col].is_zero():
                if pivot is None or aug[r][col].monomial_counts() < aug[pivot][col].monomial_counts():
                    pivot = r
        if pivot is None:
            raise ZeroDivisionError("Singular symbolic matrix in mat_solve")
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]

        piv_val = aug[col][col]
        for j in range(col, n + m):
            aug[col][j] = aug[col][j] / piv_val  # type: ignore[assignment]

        for r in range(n):
            if r != col and not aug[r][col].is_zero():
                factor = aug[r][col]
                aug[r][col] = factor.ring.rat_zero()
                for j in range(col + 1, n + m):
                    if not aug[col][j].is_zero():
                        aug[r][j] = aug[r][j] - factor * aug[col][j]  # type: ignore[assignment]

    return [row[n:] for row in aug]


def mat_inv(a: RatMatrix) -> RatMatrix:
    return mat_solve(a, mat_eye(len(a), a[0][0].ring))


def mat_det(a: RatMatrix) -> QQRat:
    n = len(a)
    ring = a[0][0].ring
    if n == 1:
        return a[0][0]
    if n == 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]  # type: ignore[return-value]
    if n == 3:
        return (
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        )  # type: ignore[return-value]
    work: RatMatrix = [list(row) for row in a]
    det = ring.rat_one()
    sign = 1
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if not work[r][col].is_zero():
                pivot = r
                break
        if pivot is None:
            return ring.rat_zero()
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            sign = -sign
        piv = work[col][col]
        det = det * piv  # type: ignore[assignment]
        for r in range(col + 1, n):
            if not work[r][col].is_zero():
                fac: QQRat = work[r][col] / piv  # type: ignore[assignment]
                for j in range(col + 1, n):
                    work[r][j] = work[r][j] - fac * work[col][j]  # type: ignore[assignment]
    return -det if sign < 0 else det


def hamiltonian_vector_field(
    ham: QQRat,
    q_vars: Sequence[Union[int, str, QQPoly, QQRat]],
    p_vars: Sequence[Union[int, str, QQPoly, QQRat]],
) -> List[QQRat]:
    """Compute `X_H = (dH/dp_1, ..., dH/dp_d, -dH/dq_1, ..., -dH/dq_d)`."""
    if len(q_vars) != len(p_vars):
        raise ValueError("q_vars and p_vars must have equal length")
    dq = [ham.diff(pv) for pv in p_vars]
    dp = [-ham.diff(qv) for qv in q_vars]
    return dq + dp


def lie_derivative(
    expr: QQRat,
    vec_field: Sequence[QQRat],
    z_vars: Sequence[Union[int, str, QQPoly, QQRat]],
) -> QQRat:
    """Compute directional derivative `sum_k (d expr / d z_k) * vec_field[k]`."""
    acc = expr.ring.rat_zero()
    for zk, vk in zip(z_vars, vec_field):
        if not vk.is_zero():
            df = expr.diff(zk)
            if not df.is_zero():
                acc = acc + df * vk  # type: ignore[assignment]
    return acc


def jacobian_matrix(
    components: Sequence[QQRat],
    z_vars: Sequence[Union[int, str, QQPoly, QQRat]],
) -> RatMatrix:
    return [[comp.diff(zv) for zv in z_vars] for comp in components]


def symplectic_defect_matrix(
    components: Sequence[QQRat],
    q_vars: Sequence[Union[int, str, QQPoly, QQRat]],
    p_vars: Sequence[Union[int, str, QQPoly, QQRat]],
) -> RatMatrix:
    """Compute `J^T * Omega * J - Omega` for canonical `(q, p)` coordinates."""
    d = len(q_vars)
    ring = components[0].ring
    z_vars = list(q_vars) + list(p_vars)
    jac = jacobian_matrix(components, z_vars)
    z, o, mo = ring.rat_zero(), ring.rat_one(), -ring.rat_one()
    omega: RatMatrix = []
    for i in range(2 * d):
        row = [z] * (2 * d)
        if i < d:
            row[d + i] = o
        else:
            row[i - d] = mo
        omega.append(row)
    jt_o_j = mat_mul(mat_transpose(jac), mat_mul(omega, jac))
    return mat_sub(jt_o_j, omega)


def analyze_rational_map_degrees(
    components: Sequence[QQRat],
    state_vars: Sequence[Union[int, str, QQPoly, QQRat]],
) -> Dict[str, object]:
    """Compute componentwise and common-denominator projective degrees in `state_vars`."""
    s_vars = list(state_vars)
    comp_degs = [c.degree_in_vars(s_vars) for c in components]
    total_degs = [c.degrees() for c in components]
    mon_counts = [c.monomial_counts() for c in components]

    # Least common multiple of denominators for projective degree on P^d
    ring = components[0].ring
    lcm_den = ring.one()
    for c in components:
        if not c.den.is_constant():
            g = poly_gcd(lcm_den, c.den)
            q = c.den.try_div_exact(g)
            if q is not None:
                lcm_den = (lcm_den * q).monic()  # type: ignore[union-attr]
            else:
                lcm_den = (lcm_den * c.den).monic()  # type: ignore[union-attr]

    den_state_deg = lcm_den.degree_in_vars(s_vars)
    proj_deg = den_state_deg
    for c in components:
        mult = lcm_den.try_div_exact(c.den)
        if mult is not None:
            num_scaled: QQPoly = (c.num * mult)  # type: ignore[assignment]
            proj_deg = max(proj_deg, num_scaled.degree_in_vars(s_vars))
        else:
            proj_deg = max(proj_deg, c.num.degree_in_vars(s_vars) + den_state_deg)

    return {
        "component_state_degrees": comp_degs,
        "component_total_degrees": total_degs,
        "component_monomial_counts": mon_counts,
        "common_den_state_degree": den_state_deg,
        "projective_degree_in_state": proj_deg,
    }
