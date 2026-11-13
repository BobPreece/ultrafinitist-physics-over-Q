"""Zero-dependency exact symbolic CAS over Q (fractions.Fraction only).

Implements:
- QQPoly: Sparse multivariate polynomials over Q with arbitrary variable sets.
- QQRat:  Canonical multivariate rational functions over Q with exact
          multivariate polynomial GCD (Subresultant PRS + Henrici reduction).
- QQJet:  Exact truncated formal Laurent/Taylor jets over Q (or QQPoly/QQRat)
          with tracked precision cutoff O(eps^cutoff) and valuation shift.

Strict Ultrafinitist Compliance:
- Pure Python stdlib (`fractions.Fraction` only).
- Zero floats, zero transcendentals, zero roots, zero sympy/numpy.
- Zero `isinstance` or `type` calls.
"""

from fractions import Fraction

ZERO_Q = Fraction(0, 1)
ONE_Q = Fraction(1, 1)

# Global interned variable registry (preserves deterministic alphabetical / creation order)
_VAR_TO_ID = {}
_ID_TO_VAR = []


def _intern_var(name):
    vid = _VAR_TO_ID.get(name)
    if vid is None:
        vid = len(_ID_TO_VAR)
        _VAR_TO_ID[name] = vid
        _ID_TO_VAR.append(name)
    return vid


def _is_exact_scalar(obj):
    """Detect int or fractions.Fraction without isinstance/type."""
    return (
        getattr(obj, "_cas_tier", 0) == 0
        and hasattr(obj, "numerator")
        and hasattr(obj, "denominator")
    )


def _to_fraction(obj):
    return Fraction(obj.numerator, obj.denominator)


def _igcd(a, b):
    """Exact non-negative integer GCD using C-accelerated Fraction reduction."""
    if a < 0:
        a = -a
    if b < 0:
        b = -b
    if a == 0:
        return b
    if b == 0:
        return a
    q = Fraction(a, b)
    return a // q.numerator


def _ilcm(a, b):
    if a == 0 or b == 0:
        return 0
    g = _igcd(a, b)
    return abs((a // g) * b)


# ============================================================================
# Sparse Monomial Utilities: monomial = tuple of (var_id, exp) with exp > 0
# ============================================================================


def _mono_mul(m1, m2):
    if not m1:
        return m2
    if not m2:
        return m1
    res = []
    i = 0
    j = 0
    n1 = len(m1)
    n2 = len(m2)
    while i < n1 and j < n2:
        v1, e1 = m1[i]
        v2, e2 = m2[j]
        if v1 < v2:
            res.append((v1, e1))
            i += 1
        elif v2 < v1:
            res.append((v2, e2))
            j += 1
        else:
            res.append((v1, e1 + e2))
            i += 1
            j += 1
    if i < n1:
        res.extend(m1[i:])
    if j < n2:
        res.extend(m2[j:])
    return tuple(res)


def _mono_div(m1, m2):
    """Return m1 / m2 if m2 divides m1, else None."""
    if not m2:
        return m1
    if len(m1) < len(m2):
        return None
    res = []
    i = 0
    j = 0
    n1 = len(m1)
    n2 = len(m2)
    while i < n1 and j < n2:
        v1, e1 = m1[i]
        v2, e2 = m2[j]
        if v1 < v2:
            res.append((v1, e1))
            i += 1
        elif v2 < v1:
            return None
        else:
            if e1 < e2:
                return None
            if e1 > e2:
                res.append((v1, e1 - e2))
            i += 1
            j += 1
    if j < n2:
        return None
    if i < n1:
        res.extend(m1[i:])
    return tuple(res)


def _mono_gcd(m1, m2):
    if not m1 or not m2:
        return ()
    res = []
    i = 0
    j = 0
    n1 = len(m1)
    n2 = len(m2)
    while i < n1 and j < n2:
        v1, e1 = m1[i]
        v2, e2 = m2[j]
        if v1 < v2:
            i += 1
        elif v2 < v1:
            j += 1
        else:
            res.append((v1, e1 if e1 <= e2 else e2))
            i += 1
            j += 1
    return tuple(res)


def _mono_deg(m):
    s = 0
    for _, e in m:
        s += e
    return s


def _mono_key(m):
    """Exact graded lexicographic monomial order key (x_0 > x_1 > x_2 > ...)."""
    return (_mono_deg(m), tuple((-v, e) for v, e in m))


# ============================================================================
# Exact Multivariate Polynomial Arithmetic & Subresultant GCD over Z
# Representation of ZZPoly: dict[monomial, int] with nonzero int values only.
# ============================================================================


def _zz_leading_item(p):
    best_m = None
    best_k = None
    for m in p:
        k = _mono_key(m)
        if best_k is None or k > best_k:
            best_k = k
            best_m = m
    return best_m, p[best_m]


def _zz_sign_canonical(p):
    if not p:
        return {}
    _, lc = _zz_leading_item(p)
    if lc < 0:
        return {m: -c for m, c in p.items()}
    return dict(p)


def _zz_int_content(p):
    g = 0
    for c in p.values():
        g = _igcd(g, c)
        if g == 1:
            return 1
    return g


def _zz_primitive(p):
    if not p:
        return 0, {}
    g = _zz_int_content(p)
    _, lc = _zz_leading_item(p)
    if lc < 0:
        g_signed = -g
    else:
        g_signed = g
    return g_signed, {m: c // g_signed for m, c in p.items()}


def _zz_add(a, b):
    if not a:
        return dict(b)
    if not b:
        return dict(a)
    out = dict(a)
    for m, c in b.items():
        nv = out.get(m, 0) + c
        if nv:
            out[m] = nv
        elif m in out:
            del out[m]
    return out


def _zz_sub(a, b):
    if not b:
        return dict(a)
    if not a:
        return {m: -c for m, c in b.items()}
    out = dict(a)
    for m, c in b.items():
        nv = out.get(m, 0) - c
        if nv:
            out[m] = nv
        elif m in out:
            del out[m]
    return out


def _zz_mul(a, b):
    if not a or not b:
        return {}
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = _mono_mul(m1, m2)
            nv = out.get(m, 0) + c1 * c2
            if nv:
                out[m] = nv
            elif m in out:
                del out[m]
    return out


def _zz_scale(a, k):
    if k == 0 or not a:
        return {}
    if k == 1:
        return dict(a)
    return {m: c * k for m, c in a.items()}


def _zz_pow(a, n):
    if n == 0:
        return {(): 1}
    if not a:
        return {}
    res = {(): 1}
    base = a
    e = n
    while e > 0:
        if e & 1:
            res = _zz_mul(res, base)
        e >>= 1
        if e > 0:
            base = _zz_mul(base, base)
    return res


def _zz_vars(p):
    vs = set()
    for m in p:
        for v, _ in m:
            vs.add(v)
    return vs


def _zz_deg_in(p, vid):
    mx = 0
    for m in p:
        for v, e in m:
            if v == vid:
                if e > mx:
                    mx = e
                break
    return mx


def _zz_total_deg(p):
    mx = 0
    for m in p:
        d = _mono_deg(m)
        if d > mx:
            mx = d
    return mx


def _zz_extract_monomial(p):
    """Factor p = mono * q where mono is the monomial GCD of all terms of p."""
    if not p:
        return (), {}
    it = iter(p)
    mg = next(it)
    if not mg:
        return (), dict(p)
    for m in it:
        mg = _mono_gcd(mg, m)
        if not mg:
            return (), dict(p)
    return mg, {_mono_div(m, mg): c for m, c in p.items()}


def _zz_exact_div(a, b):
    """Exact division a / b over Z[vars]; returns quotient dict or None."""
    if not b:
        raise ZeroDivisionError("Polynomial division by zero")
    if not a:
        return {}
    if len(b) == 1:
        bm, bc = next(iter(b.items()))
        q = {}
        for am, ac in a.items():
            qm = _mono_div(am, bm)
            if qm is None or (ac % bc) != 0:
                return None
            q[qm] = ac // bc
        return q
    # Quick degree and variable checks
    if _zz_total_deg(a) < _zz_total_deg(b):
        return None
    if not _zz_vars(b).issubset(_zz_vars(a)):
        return None
    bm_lead, bc_lead = _zz_leading_item(b)
    rem = dict(a)
    quot = {}
    b_deg = _zz_total_deg(b)
    while rem:
        rm_lead, rc_lead = _zz_leading_item(rem)
        qm = _mono_div(rm_lead, bm_lead)
        if qm is None or (rc_lead % bc_lead) != 0:
            return None
        if _mono_deg(rm_lead) < b_deg and qm != ():
            return None
        qc = rc_lead // bc_lead
        quot[qm] = quot.get(qm, 0) + qc
        # Subtract qc * x^qm * b from rem
        for bm, bc in b.items():
            tm = _mono_mul(qm, bm)
            nv = rem.get(tm, 0) - qc * bc
            if nv:
                rem[tm] = nv
            elif tm in rem:
                del rem[tm]
    return quot


def _zz_to_univar(p, vid):
    """Decompose p in Z[vars] as dict[exp_of_vid, coeff_zz_poly_without_vid]."""
    uni = {}
    for m, c in p.items():
        deg_v = 0
        rest = []
        for v, e in m:
            if v == vid:
                deg_v = e
            else:
                rest.append((v, e))
        rm = tuple(rest)
        bucket = uni.get(deg_v)
        if bucket is None:
            uni[deg_v] = {rm: c}
        else:
            bucket[rm] = c
    return uni


def _zz_from_univar(uni, vid):
    out = {}
    for deg_v, coeff_p in uni.items():
        if not coeff_p:
            continue
        if deg_v == 0:
            for rm, c in coeff_p.items():
                if c:
                    out[rm] = out.get(rm, 0) + c
        else:
            v_mono = ((vid, deg_v),)
            for rm, c in coeff_p.items():
                if c:
                    fm = _mono_mul(rm, v_mono)
                    out[fm] = out.get(fm, 0) + c
    return {m: c for m, c in out.items() if c != 0}


def _uni_deg(u):
    return max(u) if u else -1


def _uni_scale_poly(u, poly_c):
    if not u or not poly_c:
        return {}
    out = {}
    for d, coeff in u.items():
        pr = _zz_mul(coeff, poly_c)
        if pr:
            out[d] = pr
    return out


def _uni_exact_div_poly(u, poly_c):
    if not u:
        return {}
    out = {}
    for d, coeff in u.items():
        q = _zz_exact_div(coeff, poly_c)
        if q is None:
            return None
        if q:
            out[d] = q
    return out


def _uni_prem(u_a, u_b):
    """Exact pseudo-remainder of u_a by u_b in R[v], scaled by lc(u_b)^(delta + 1)."""
    deg_a = _uni_deg(u_a)
    deg_b = _uni_deg(u_b)
    delta = deg_a - deg_b
    lc_b = u_b[deg_b]
    rem = dict(u_a)
    steps = 0
    while rem:
        deg_r = _uni_deg(rem)
        if deg_r < deg_b:
            break
        lc_r = rem[deg_r]
        shift = deg_r - deg_b
        # rem = lc_b * rem - lc_r * v^shift * u_b
        new_rem = {}
        for d, coeff in rem.items():
            if d != deg_r:
                pr = _zz_mul(lc_b, coeff)
                if pr:
                    new_rem[d] = pr
        for d, coeff in u_b.items():
            if d != deg_b:
                td = d + shift
                sub_term = _zz_mul(lc_r, coeff)
                cur = new_rem.get(td, {})
                diff = _zz_sub(cur, sub_term)
                if diff:
                    new_rem[td] = diff
                elif td in new_rem:
                    del new_rem[td]
        rem = new_rem
        steps += 1
    extra = (delta + 1) - steps
    if extra > 0 and rem:
        rem = _uni_scale_poly(rem, _zz_pow(lc_b, extra))
    return rem


def _zz_project_out_vars(p, drop_vars):
    """Group p by monomials in drop_vars and return list of coefficient polynomials."""
    buckets = {}
    for m, c in p.items():
        drop_m = []
        keep_m = []
        for v, e in m:
            if v in drop_vars:
                drop_m.append((v, e))
            else:
                keep_m.append((v, e))
        dk = tuple(drop_m)
        km = tuple(keep_m)
        b = buckets.get(dk)
        if b is None:
            buckets[dk] = {km: c}
        else:
            b[km] = c
    return list(buckets.values())


def _zz_gcd_list(polys):
    g = {}
    for p in polys:
        if not p:
            continue
        g = _zz_gcd(g, p)
        if len(g) == 1 and () in g and abs(g[()]) == 1:
            return {(): 1}
    return g if g else {(): 1}


def _zz_gcd(a, b):
    """Exact canonical GCD of two multivariate polynomials over Z."""
    if not a:
        return _zz_sign_canonical(b)
    if not b:
        return _zz_sign_canonical(a)

    # 1. Integer content extraction
    ca, a_prim = _zz_primitive(a)
    cb, b_prim = _zz_primitive(b)
    cg = _igcd(ca, cb)

    if _zz_total_deg(a_prim) == 0 or _zz_total_deg(b_prim) == 0:
        return {(): cg}
    if a_prim == b_prim:
        return _zz_scale(a_prim, cg)

    # 2. Monomial GCD extraction
    ma, a_nom = _zz_extract_monomial(a_prim)
    mb, b_nom = _zz_extract_monomial(b_prim)
    mg = _mono_gcd(ma, mb)

    if _zz_total_deg(a_nom) == 0 or _zz_total_deg(b_nom) == 0:
        return {mg: cg}
    if a_nom == b_nom:
        return {_mono_mul(mg, m): c * cg for m, c in a_nom.items()}

    # 3. Project out non-shared variables immediately!
    va = _zz_vars(a_nom)
    vb = _zz_vars(b_nom)
    shared = va & vb
    if not shared:
        return {mg: cg}
    if va != shared:
        a_proj = _zz_gcd_list(_zz_project_out_vars(a_nom, va - shared))
        core = _zz_gcd(a_proj, b_nom)
        return {_mono_mul(mg, m): c * cg for m, c in core.items()}
    if vb != shared:
        b_proj = _zz_gcd_list(_zz_project_out_vars(b_nom, vb - shared))
        core = _zz_gcd(a_nom, b_proj)
        return {_mono_mul(mg, m): c * cg for m, c in core.items()}

    # 4. Fast exact divisibility & linear irreducibility checks
    deg_a = _zz_total_deg(a_nom)
    deg_b = _zz_total_deg(b_nom)
    if deg_a == 1:
        core = a_nom if _zz_exact_div(b_nom, a_nom) is not None else {(): 1}
        return {_mono_mul(mg, m): c * cg for m, c in core.items()}
    if deg_b == 1:
        core = b_nom if _zz_exact_div(a_nom, b_nom) is not None else {(): 1}
        return {_mono_mul(mg, m): c * cg for m, c in core.items()}
    if deg_a >= deg_b and _zz_exact_div(a_nom, b_nom) is not None:
        return {_mono_mul(mg, m): c * cg for m, c in b_nom.items()}
    if deg_b > deg_a and _zz_exact_div(b_nom, a_nom) is not None:
        return {_mono_mul(mg, m): c * cg for m, c in a_nom.items()}

    # 4b. Fast exact univariate monic Euclidean GCD over Q (Gauss's Lemma) when len(shared) == 1!
    if len(shared) == 1:
        vid = next(iter(shared))
        ra = {}
        for m, c in a_nom.items():
            ra[m[0][1] if m else 0] = Fraction(c, 1)
        rb = {}
        for m, c in b_nom.items():
            rb[m[0][1] if m else 0] = Fraction(c, 1)
        if max(ra) < max(rb):
            ra, rb = rb, ra
        # Normalize rb monic
        db_u = max(rb)
        inv_lcb = Fraction(1, 1) / rb[db_u]
        rb = {d: v * inv_lcb for d, v in rb.items()}
        while rb:
            db_u = max(rb)
            if db_u == 0:
                return {mg: cg}
            while ra and max(ra) >= db_u:
                da_u = max(ra)
                lead_a = ra[da_u]
                shift = da_u - db_u
                del ra[da_u]
                for d, vb in rb.items():
                    if d != db_u:
                        td = d + shift
                        nv = ra.get(td, Fraction(0, 1)) - lead_a * vb
                        if nv:
                            ra[td] = nv
                        elif td in ra:
                            del ra[td]
            if not ra:
                break
            da_u = max(ra)
            inv_lca = Fraction(1, 1) / ra[da_u]
            ra = {d: v * inv_lca for d, v in ra.items()}
            ra, rb = rb, ra
        # Convert monic rb in Q[vid] back to primitive Z[vid]
        lcm_den = 1
        for v in rb.values():
            lcm_den = (lcm_den // _igcd(lcm_den, v.denominator)) * v.denominator
        zz_u = {(((vid, d),) if d > 0 else ()): v.numerator * (lcm_den // v.denominator) for d, v in rb.items()}
        _, prim_u = _zz_primitive(zz_u)
        return {_mono_mul(mg, m): c * cg for m, c in prim_u.items()}

    # 5. Pick main variable v in `shared` minimizing max(deg_v(A), deg_v(B))
    best_v = None
    best_score = None
    for v in shared:
        da = _zz_deg_in(a_nom, v)
        db = _zz_deg_in(b_nom, v)
        hi = da if da >= db else db
        lo = db if da >= db else da
        score = (hi, lo, v)
        if best_score is None or score < best_score:
            best_score = score
            best_v = v

    ua = _zz_to_univar(a_nom, best_v)
    ub = _zz_to_univar(b_nom, best_v)

    cont_a = _zz_gcd_list(list(ua.values()))
    cont_b = _zz_gcd_list(list(ub.values()))
    cont_g = _zz_gcd(cont_a, cont_b)

    ua = _uni_exact_div_poly(ua, cont_a)
    ub = _uni_exact_div_poly(ub, cont_b)

    if _uni_deg(ua) < _uni_deg(ub):
        ua, ub = ub, ua

    # Gauss's Lemma shortcut: primitive degree-1 polynomial in best_v is irreducible!
    if _uni_deg(ub) == 1:
        b_irred = _zz_from_univar(ub, best_v)
        a_red = _zz_from_univar(ua, best_v)
        prim_g = b_irred if _zz_exact_div(a_red, b_irred) is not None else {(): 1}
        full_g = _zz_sign_prod(cont_g, prim_g)
        return {_mono_mul(mg, m): c * cg for m, c in full_g.items()}

    # 6. Collins-Brown-Lazard Subresultant PRS over R[best_v]
    g_pol = {(): 1}
    h_pol = {(): 1}
    while _uni_deg(ub) > 0:
        delta = _uni_deg(ua) - _uni_deg(ub)
        rem = _uni_prem(ua, ub)
        if not rem:
            # ub (made primitive in best_v) is the primitive GCD!
            cb_u = _zz_gcd_list(list(ub.values()))
            ub_prim = _uni_exact_div_poly(ub, cb_u)
            prim_g = _zz_from_univar(ub_prim, best_v)
            full_g = _zz_sign_prod(cont_g, prim_g)
            return {_mono_mul(mg, m): c * cg for m, c in full_g.items()}
        if _uni_deg(rem) == 0:
            # Coprime in best_v!
            full_g = _zz_sign_canonical(cont_g)
            return {_mono_mul(mg, m): c * cg for m, c in full_g.items()}

        denom = _zz_mul(g_pol, _zz_pow(h_pol, delta))
        ua = ub
        ub_next = _uni_exact_div_poly(rem, denom)
        if ub_next is None:
            # Defensive fallback to primitive PRS if ever needed
            cr = _zz_gcd_list(list(rem.values()))
            ub_next = _uni_exact_div_poly(rem, cr)
            g_pol = {(): 1}
            h_pol = {(): 1}
        else:
            g_pol = ua[_uni_deg(ua)]
            if delta == 0:
                pass  # h_pol unchanged
            elif delta == 1:
                h_pol = g_pol
            else:
                h_pol = _zz_exact_div(_zz_pow(g_pol, delta), _zz_pow(h_pol, delta - 1))
        ub = ub_next

    full_g = _zz_sign_canonical(cont_g)
    return {_mono_mul(mg, m): c * cg for m, c in full_g.items()}


def _zz_sign_prod(a, b):
    return _zz_sign_canonical(_zz_mul(a, b))


# ============================================================================
# QQPoly: Multivariate Polynomial over Q
# ============================================================================


class QQPoly:
    """Exact sparse multivariate polynomial over Q (`fractions.Fraction`)."""

    _cas_tier = 1
    __slots__ = ("_terms", "_hash_cache")

    def __init__(self, val=0):
        self._hash_cache = None
        tier = getattr(val, "_cas_tier", 0)
        if tier == 1:
            self._terms = dict(val._terms)
            return
        if _is_exact_scalar(val):
            q = _to_fraction(val)
            self._terms = {(): q} if q != ZERO_Q else {}
            return
        if hasattr(val, "items"):
            cleaned = {}
            for m, c in val.items():
                qc = _to_fraction(c)
                if qc != ZERO_Q:
                    cleaned[m] = qc
            self._terms = cleaned
            return
        raise TypeError(f"Cannot construct QQPoly from {val!r}")

    @classmethod
    def _from_raw(cls, terms):
        obj = cls.__new__(cls)
        obj._terms = terms
        obj._hash_cache = None
        return obj

    @classmethod
    def zero(cls):
        return cls._from_raw({})

    @classmethod
    def one(cls):
        return cls._from_raw({(): ONE_Q})

    @classmethod
    def const(cls, c):
        return cls(c)

    @classmethod
    def var(cls, name):
        vid = _intern_var(str(name))
        return cls._from_raw({((vid, 1),): ONE_Q})

    def is_zero(self):
        return not self._terms

    def is_constant(self):
        if not self._terms:
            return True
        return len(self._terms) == 1 and () in self._terms

    def constant_value(self):
        if not self._terms:
            return ZERO_Q
        if len(self._terms) == 1 and () in self._terms:
            return self._terms[()]
        raise ValueError(f"Polynomial {self} is not constant")

    def degree(self):
        mx = 0
        for m in self._terms:
            d = _mono_deg(m)
            if d > mx:
                mx = d
        return mx

    def degree_in(self, name):
        vid = _VAR_TO_ID.get(str(name))
        if vid is None:
            return 0
        return _zz_deg_in(self._terms, vid)

    def variables(self):
        vids = _zz_vars(self._terms)
        return tuple(sorted(_ID_TO_VAR[v] for v in vids))

    def leading_monomial_and_coeff(self):
        if not self._terms:
            return (), ZERO_Q
        best_m = None
        best_k = None
        for m in self._terms:
            k = _mono_key(m)
            if best_k is None or k > best_k:
                best_k = k
                best_m = m
        return best_m, self._terms[best_m]

    def monic(self):
        if not self._terms:
            return self
        _, lc = self.leading_monomial_and_coeff()
        if lc == ONE_Q:
            return self
        inv = ONE_Q / lc
        return QQPoly._from_raw({m: c * inv for m, c in self._terms.items()})

    def to_zz_primitive(self):
        """Return (scale_q, zz_dict) such that self == scale_q * zz_dict and zz_dict is primitive."""
        if not self._terms:
            return ZERO_Q, {}
        den_lcm = 1
        for c in self._terms.values():
            den_lcm = _ilcm(den_lcm, c.denominator)
        zz_raw = {m: (c.numerator * (den_lcm // c.denominator)) for m, c in self._terms.items()}
        g_signed, zz_prim = _zz_primitive(zz_raw)
        return Fraction(g_signed, den_lcm), zz_prim

    @classmethod
    def from_zz(cls, zz_dict, scale=ONE_Q):
        qscale = _to_fraction(scale)
        if qscale == ZERO_Q or not zz_dict:
            return cls.zero()
        return cls._from_raw({m: Fraction(c, 1) * qscale for m, c in zz_dict.items() if c != 0})

    def terms(self):
        items = []
        for m in sorted(self._terms, key=_mono_key, reverse=True):
            named_m = tuple((_ID_TO_VAR[v], e) for v, e in m)
            items.append((named_m, self._terms[m]))
        return items

    def diff(self, name):
        vid = _VAR_TO_ID.get(str(name))
        if vid is None or not self._terms:
            return QQPoly.zero()
        out = {}
        for m, c in self._terms.items():
            deg_v = 0
            new_m = []
            for v, e in m:
                if v == vid:
                    deg_v = e
                    if e > 1:
                        new_m.append((v, e - 1))
                else:
                    new_m.append((v, e))
            if deg_v > 0:
                km = tuple(new_m)
                nc = c * deg_v
                nv = out.get(km, ZERO_Q) + nc
                if nv != ZERO_Q:
                    out[km] = nv
                elif km in out:
                    del out[km]
        return QQPoly._from_raw(out)

    def eval(self, env):
        if not self._terms:
            return ZERO_Q
        acc = ZERO_Q
        for m, c in self._terms.items():
            term = c
            for v, e in m:
                vname = _ID_TO_VAR[v]
                val = env[vname] if vname in env else QQPoly._from_raw({((v, 1),): ONE_Q})
                term = term * (val ** e)
            acc = acc + term
        return acc

    # --- Arithmetic Dunders ---

    def __bool__(self):
        return bool(self._terms)

    def __pos__(self):
        return self

    def __neg__(self):
        if not self._terms:
            return self
        return QQPoly._from_raw({m: -c for m, c in self._terms.items()})

    def _coerce_other(self, other):
        tier = getattr(other, "_cas_tier", 0)
        if tier == 1:
            return other
        if tier > 1:
            return None  # Let higher tier handle via Reflected op
        if _is_exact_scalar(other):
            q = _to_fraction(other)
            return QQPoly._from_raw({(): q} if q != ZERO_Q else {})
        return None

    def __add__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if not self._terms:
            return co
        if not co._terms:
            return self
        out = dict(self._terms)
        for m, c in co._terms.items():
            nv = out.get(m, ZERO_Q) + c
            if nv != ZERO_Q:
                out[m] = nv
            elif m in out:
                del out[m]
        return QQPoly._from_raw(out)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if not co._terms:
            return self
        if not self._terms:
            return -co
        out = dict(self._terms)
        for m, c in co._terms.items():
            nv = out.get(m, ZERO_Q) - c
            if nv != ZERO_Q:
                out[m] = nv
            elif m in out:
                del out[m]
        return QQPoly._from_raw(out)

    def __rsub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__sub__(self)

    def __mul__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if not self._terms or not co._terms:
            return QQPoly.zero()
        if co.is_constant():
            k = co._terms[()]
            if k == ONE_Q:
                return self
            return QQPoly._from_raw({m: c * k for m, c in self._terms.items()})
        if self.is_constant():
            k = self._terms[()]
            if k == ONE_Q:
                return co
            return QQPoly._from_raw({m: c * k for m, c in co._terms.items()})
        out = {}
        for m1, c1 in self._terms.items():
            for m2, c2 in co._terms.items():
                m = _mono_mul(m1, m2)
                nv = out.get(m, ZERO_Q) + c1 * c2
                if nv != ZERO_Q:
                    out[m] = nv
                elif m in out:
                    del out[m]
        return QQPoly._from_raw(out)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if not co._terms:
            raise ZeroDivisionError("Division by zero QQPoly")
        if co.is_constant():
            inv = ONE_Q / co._terms[()]
            return QQPoly._from_raw({m: c * inv for m, c in self._terms.items()})
        q_exact = poly_exact_div(self, co)
        if q_exact is not None:
            return q_exact
        return QQRat(self, co)

    def __rtruediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__truediv__(self)

    def __pow__(self, exp):
        if not _is_exact_scalar(exp) or _to_fraction(exp).denominator != 1:
            return NotImplemented
        n = int(_to_fraction(exp).numerator)
        if n < 0:
            return QQRat(ONE_Q, self ** (-n))
        if n == 0:
            return QQPoly.one()
        if not self._terms:
            return QQPoly.zero()
        res = QQPoly.one()
        base = self
        e = n
        while e > 0:
            if e & 1:
                res = res * base
            e >>= 1
            if e > 0:
                base = base * base
        return res

    def __eq__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self._terms == co._terms

    def __hash__(self):
        if self._hash_cache is None:
            if self.is_constant():
                self._hash_cache = hash(self.constant_value())
            else:
                self._hash_cache = hash(tuple(sorted(self._terms.items())))
        return self._hash_cache

    def __repr__(self):
        if not self._terms:
            return "0"
        parts = []
        for named_m, c in self.terms():
            if not named_m:
                parts.append(str(c))
            else:
                m_str = "*".join(v if e == 1 else f"{v}^{e}" for v, e in named_m)
                if c == ONE_Q:
                    parts.append(m_str)
                elif c == -ONE_Q:
                    parts.append(f"-{m_str}")
                else:
                    parts.append(f"({c})*{m_str}")
        return " + ".join(parts).replace("+ -", "- ")


def _poly_vars_ids(p):
    vs = set()
    for m in p._terms:
        for v, _ in m:
            vs.add(v)
    return vs


def poly_gcd(a, b):
    """Exact monic multivariate polynomial GCD over Q."""
    pa = a if getattr(a, "_cas_tier", 0) == 1 else QQPoly(a)
    pb = b if getattr(b, "_cas_tier", 0) == 1 else QQPoly(b)
    if pa.is_zero():
        return pb.monic()
    if pb.is_zero():
        return pa.monic()
    if pa.is_constant() or pb.is_constant():
        return QQPoly.one()
    if pa._terms == pb._terms:
        return pa.monic()
    v_all = _poly_vars_ids(pa) | _poly_vars_ids(pb)
    if len(v_all) == 1:
        vid = next(iter(v_all))
        ra = { (m[0][1] if m else 0): c for m, c in pa._terms.items() }
        rb = { (m[0][1] if m else 0): c for m, c in pb._terms.items() }
        if max(ra) < max(rb):
            ra, rb = rb, ra
        inv_lcb = ONE_Q / rb[max(rb)]
        rb = {d: v * inv_lcb for d, v in rb.items()}
        while rb:
            db_u = max(rb)
            if db_u == 0:
                return QQPoly.one()
            while ra and max(ra) >= db_u:
                da_u = max(ra)
                lead_a = ra[da_u]
                shift = da_u - db_u
                del ra[da_u]
                for d, vb in rb.items():
                    if d != db_u:
                        td = d + shift
                        nv = ra.get(td, ZERO_Q) - lead_a * vb
                        if nv != ZERO_Q:
                            ra[td] = nv
                        elif td in ra:
                            del ra[td]
            if not ra:
                break
            inv_lca = ONE_Q / ra[max(ra)]
            ra = {d: v * inv_lca for d, v in ra.items()}
            ra, rb = rb, ra
        return QQPoly._from_raw({ (((vid, d),) if d > 0 else ()): c for d, c in rb.items() })
    _, za = pa.to_zz_primitive()
    _, zb = pb.to_zz_primitive()
    zg = _zz_gcd(za, zb)
    return QQPoly.from_zz(zg).monic()


def poly_exact_div(a, b):
    """Return exact polynomial quotient a / b over Q as QQPoly, or None if not divisible."""
    pa = a if getattr(a, "_cas_tier", 0) == 1 else QQPoly(a)
    pb = b if getattr(b, "_cas_tier", 0) == 1 else QQPoly(b)
    if pb.is_zero():
        raise ZeroDivisionError("Polynomial division by zero")
    if pa.is_zero():
        return QQPoly.zero()
    if pb.is_constant():
        return pa / pb.constant_value()
    if pa._terms == pb._terms:
        return QQPoly.one()
    v_all = _poly_vars_ids(pa) | _poly_vars_ids(pb)
    if len(v_all) == 1:
        vid = next(iter(v_all))
        rem = { (m[0][1] if m else 0): c for m, c in pa._terms.items() }
        div = { (m[0][1] if m else 0): c for m, c in pb._terms.items() }
        db_u = max(div)
        inv_lcb = ONE_Q / div[db_u]
        quot = {}
        while rem:
            dr = max(rem)
            if dr < db_u:
                return None
            q_c = rem[dr] * inv_lcb
            shift = dr - db_u
            quot[shift] = q_c
            del rem[dr]
            for d, vb in div.items():
                if d != db_u:
                    td = d + shift
                    nv = rem.get(td, ZERO_Q) - q_c * vb
                    if nv != ZERO_Q:
                        rem[td] = nv
                    elif td in rem:
                        del rem[td]
        return QQPoly._from_raw({ (((vid, d),) if d > 0 else ()): c for d, c in quot.items() if c != ZERO_Q })
    sa, za = pa.to_zz_primitive()
    sb, zb = pb.to_zz_primitive()
    zq = _zz_exact_div(za, zb)
    if zq is None:
        return None
    return QQPoly.from_zz(zq, scale=sa / sb)


# ============================================================================
# QQRat: Canonical Multivariate Rational Function over Q
# ============================================================================


class QQRat:
    """Exact canonical multivariate rational function num / den over Q.

    Invariants:
    - `num` and `den` are `QQPoly`.
    - `den` is monic (`leading_coeff == 1`), nonzero.
    - `poly_gcd(num, den) == 1` (so any rational identity equal to 0 reduces to `0 / 1`).
    """

    _cas_tier = 2
    __slots__ = ("num", "den", "_hash_cache")

    def __init__(self, num=0, den=1):
        self._hash_cache = None
        t_num = getattr(num, "_cas_tier", 0)
        t_den = getattr(den, "_cas_tier", 0)

        if t_num == 2 and t_den == 0 and _is_exact_scalar(den) and _to_fraction(den) == ONE_Q:
            self.num = num.num
            self.den = num.den
            return

        if t_num == 2 or t_den == 2:
            rn = num if t_num == 2 else QQRat._from_poly_const(num)
            rd = den if t_den == 2 else QQRat._from_poly_const(den)
            res = rn / rd
            self.num = res.num
            self.den = res.den
            return

        p_num = num if t_num == 1 else QQPoly(num)
        p_den = den if t_den == 1 else QQPoly(den)
        n_can, d_can = QQRat._reduce_pair(p_num, p_den)
        self.num = n_can
        self.den = d_can

    @classmethod
    def _from_reduced(cls, num_poly, den_poly):
        """Construct from coprime (num_poly, den_poly), normalizing leading coeff of den_poly."""
        if den_poly.is_zero():
            raise ZeroDivisionError("Zero denominator in QQRat")
        if num_poly.is_zero():
            obj = cls.__new__(cls)
            obj.num = QQPoly.zero()
            obj.den = QQPoly.one()
            obj._hash_cache = None
            return obj
        _, lc = den_poly.leading_monomial_and_coeff()
        if lc != ONE_Q:
            inv = ONE_Q / lc
            num_poly = num_poly * inv
            den_poly = den_poly * inv
        obj = cls.__new__(cls)
        obj.num = num_poly
        obj.den = den_poly
        obj._hash_cache = None
        return obj

    @staticmethod
    def _reduce_pair(p_num, p_den):
        if p_den.is_zero():
            raise ZeroDivisionError("Zero denominator in QQRat")
        if p_num.is_zero():
            return QQPoly.zero(), QQPoly.one()
        if p_den.is_constant():
            inv = ONE_Q / p_den.constant_value()
            return p_num * inv, QQPoly.one()
        if p_num.is_constant():
            _, lc = p_den.leading_monomial_and_coeff()
            inv = ONE_Q / lc
            return p_num * inv, p_den * inv
        g = poly_gcd(p_num, p_den)
        if not g.is_constant():
            p_num = poly_exact_div(p_num, g)
            p_den = poly_exact_div(p_den, g)
        _, lc = p_den.leading_monomial_and_coeff()
        if lc != ONE_Q:
            inv = ONE_Q / lc
            p_num = p_num * inv
            p_den = p_den * inv
        return p_num, p_den

    @classmethod
    def _from_poly_const(cls, val):
        tier = getattr(val, "_cas_tier", 0)
        if tier == 2:
            return val
        p = val if tier == 1 else QQPoly(val)
        return cls._from_reduced(p, QQPoly.one())

    @classmethod
    def zero(cls):
        return cls._from_reduced(QQPoly.zero(), QQPoly.one())

    @classmethod
    def one(cls):
        return cls._from_reduced(QQPoly.one(), QQPoly.one())

    @classmethod
    def const(cls, c):
        return cls._from_poly_const(c)

    @classmethod
    def var(cls, name):
        return cls._from_reduced(QQPoly.var(name), QQPoly.one())

    def is_zero(self):
        return self.num.is_zero()

    def is_polynomial(self):
        return self.den.is_constant()

    def is_constant(self):
        return self.num.is_constant() and self.den.is_constant()

    def constant_value(self):
        if not self.is_constant():
            raise ValueError(f"QQRat {self} is not constant")
        return self.num.constant_value() / self.den.constant_value()

    def variables(self):
        return tuple(sorted(set(self.num.variables()) | set(self.den.variables())))

    def diff(self, name):
        dn = self.num.diff(name)
        dd = self.den.diff(name)
        if dd.is_zero():
            return QQRat(dn, self.den)
        return (QQRat._from_poly_const(dn) * QQRat._from_poly_const(self.den) -
                QQRat._from_poly_const(self.num) * QQRat._from_poly_const(dd)) / (
                    QQRat._from_poly_const(self.den) ** 2
                )

    def eval(self, env):
        ev_n = self.num.eval(env)
        ev_d = self.den.eval(env)
        return ev_n / ev_d

    # --- Henrici-Reduced Arithmetic Dunders ---

    def _coerce_other(self, other):
        tier = getattr(other, "_cas_tier", 0)
        if tier == 2:
            return other
        if tier > 2:
            return None
        if tier == 1 or _is_exact_scalar(other):
            return QQRat._from_poly_const(other)
        return None

    def __bool__(self):
        return not self.num.is_zero()

    def __pos__(self):
        return self

    def __neg__(self):
        if self.num.is_zero():
            return self
        return QQRat._from_reduced(-self.num, self.den)

    def __add__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.num.is_zero():
            return co
        if co.num.is_zero():
            return self
        # Fast path 1: both polynomials (den == 1)
        if self.den.is_constant() and co.den.is_constant():
            return QQRat._from_reduced(self.num + co.num, QQPoly.one())
        # Fast path 2: identical denominators (common in vector norms!)
        if self.den == co.den:
            s = self.num + co.num
            n_can, d_can = QQRat._reduce_pair(s, self.den)
            return QQRat._from_reduced(n_can, d_can)
        # Henrici addition: g = gcd(d1, d2)
        g = poly_gcd(self.den, co.den)
        if g.is_constant():
            num_raw = self.num * co.den + co.num * self.den
            if num_raw.is_zero():
                return QQRat.zero()
            return QQRat._from_reduced(num_raw, self.den * co.den)
        d1_bar = poly_exact_div(self.den, g)
        d2_bar = poly_exact_div(co.den, g)
        num_raw = self.num * d2_bar + co.num * d1_bar
        if num_raw.is_zero():
            return QQRat.zero()
        g2 = poly_gcd(num_raw, g)
        if not g2.is_constant():
            num_raw = poly_exact_div(num_raw, g2)
            g = poly_exact_div(g, g2)
        return QQRat._from_reduced(num_raw, d1_bar * (d2_bar * g))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self.__add__(-co)

    def __rsub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__add__(-self)

    def __mul__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if self.num.is_zero() or co.num.is_zero():
            return QQRat.zero()
        if self.den.is_constant() and co.den.is_constant():
            return QQRat._from_reduced(self.num * co.num, QQPoly.one())
        # Henrici cross-cancellation
        n1, d1 = self.num, self.den
        n2, d2 = co.num, co.den
        if not d2.is_constant() and not n1.is_constant():
            g1 = poly_gcd(n1, d2)
            if not g1.is_constant():
                n1 = poly_exact_div(n1, g1)
                d2 = poly_exact_div(d2, g1)
        if not d1.is_constant() and not n2.is_constant():
            g2 = poly_gcd(n2, d1)
            if not g2.is_constant():
                n2 = poly_exact_div(n2, g2)
                d1 = poly_exact_div(d1, g2)
        return QQRat._from_reduced(n1 * n2, d1 * d2)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if co.num.is_zero():
            raise ZeroDivisionError("Division by zero QQRat")
        if self.num.is_zero():
            return QQRat.zero()
        inv_co = QQRat._from_reduced(co.den, co.num)
        return self.__mul__(inv_co)

    def __rtruediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__truediv__(self)

    def __pow__(self, exp):
        if not _is_exact_scalar(exp) or _to_fraction(exp).denominator != 1:
            return NotImplemented
        n = int(_to_fraction(exp).numerator)
        if n == 0:
            return QQRat.one()
        if n < 0:
            if self.num.is_zero():
                raise ZeroDivisionError("0/1 cannot be raised to a negative power")
            inv = QQRat._from_reduced(self.den, self.num)
            return inv ** (-n)
        if self.num.is_zero():
            return QQRat.zero()
        # Since gcd(num, den) == 1, gcd(num^n, den^n) == 1 automatically!
        return QQRat._from_reduced(self.num ** n, self.den ** n)

    def __eq__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self.num == co.num and self.den == co.den

    def __hash__(self):
        if self._hash_cache is None:
            if self.is_constant():
                self._hash_cache = hash(self.constant_value())
            else:
                self._hash_cache = hash((self.num, self.den))
        return self._hash_cache

    def __repr__(self):
        if self.den.is_constant():
            return repr(self.num)
        return f"({self.num!r})/({self.den!r})"


# ============================================================================
# QQJet: Exact Truncated Formal Laurent/Taylor Jet in `eps`
# ============================================================================


def _coeff_is_zero(c):
    if _is_exact_scalar(c):
        return _to_fraction(c) == ZERO_Q
    return c.is_zero()


class QQJet:
    """Exact truncated formal Laurent/Taylor jet:
       sum_{k=0}^{len(coeffs)-1} coeffs[k] * eps^(val + k) + O(eps^cutoff)

    Supports coefficients in `Fraction`, `QQPoly`, or `QQRat`.
    Handles intermediate Laurent poles (`val < 0`) cleanly so expressions like
    `(<u, dt*g> / <dt*g, dt*g>) * (dt*g)` evaluate without ZeroDivisionError
    and without spurious precision loss.
    """

    _cas_tier = 3
    _audit_jumps = False
    _strict_half_ulp = False
    _nonconst_int_queries = []
    _checked_mantissa_bits = (54,)
    _checked_max_shift = 256
    __slots__ = ("val", "coeffs", "cutoff")

    def __init__(self, val=0, coeffs=(), cutoff=None):
        if coeffs == () and cutoff is None and (
            _is_exact_scalar(val) or getattr(val, "_cas_tier", 0) in (1, 2, 3)
        ):
            tier = getattr(val, "_cas_tier", 0)
            if tier == 3:
                self.val = val.val
                self.coeffs = val.coeffs
                self.cutoff = val.cutoff
                return
            c = _to_fraction(val) if tier == 0 else val
            if _coeff_is_zero(c):
                self.val = 0
                self.coeffs = ()
                self.cutoff = None
            else:
                self.val = 0
                self.coeffs = (c,)
                self.cutoff = None
            return
        v_norm, c_norm, cut_norm = QQJet._normalize(val, coeffs, cutoff)
        self.val = v_norm
        self.coeffs = c_norm
        self.cutoff = cut_norm

    @staticmethod
    def _normalize(val, coeffs, cutoff):
        lst = list(coeffs)
        if cutoff is not None:
            max_len = cutoff - val
            if max_len <= 0:
                return cutoff, (), cutoff
            if len(lst) > max_len:
                lst = lst[:max_len]
        # Strip leading zeros
        idx = 0
        n = len(lst)
        while idx < n and _coeff_is_zero(lst[idx]):
            idx += 1
        if idx == n:
            return (cutoff if cutoff is not None else 0), (), cutoff
        new_val = val + idx
        lst = lst[idx:]
        # Strip trailing zeros cleanly (cutoff tracks exact truncation depth!)
        while lst and _coeff_is_zero(lst[-1]):
            lst.pop()
        return new_val, tuple(lst), cutoff

    @classmethod
    def _raw(cls, val, coeffs, cutoff):
        obj = cls.__new__(cls)
        v_n, c_n, cut_n = cls._normalize(val, coeffs, cutoff)
        obj.val = v_n
        obj.coeffs = c_n
        obj.cutoff = cut_n
        return obj

    @classmethod
    def const(cls, c, order=None):
        val_c = _to_fraction(c) if _is_exact_scalar(c) else c
        cutoff = (order + 1) if order is not None else None
        if _coeff_is_zero(val_c):
            return cls._raw(0, (), cutoff)
        return cls._raw(0, (val_c,), cutoff)

    @classmethod
    def var(cls, order=4, one=ONE_Q):
        """Return `0 + 1*eps + O(eps^(order + 1))`."""
        return cls._raw(1, (one,), order + 1)

    @classmethod
    def dt(cls, order=4, one=ONE_Q):
        """Alias for `QQJet.var(order=order)` representing infinitesimal `dt`."""
        return cls.var(order=order, one=one)

    @property
    def order(self):
        """Highest guaranteed exact exponent k of eps^k (None if exact to all orders)."""
        return (self.cutoff - 1) if self.cutoff is not None else None

    def coeff(self, k):
        """Return exact coefficient of eps^k; raises ValueError if k >= cutoff."""
        if self.cutoff is not None and k >= self.cutoff:
            raise ValueError(
                f"Requested eps^{k} coefficient from QQJet only accurate up to O(eps^{self.cutoff})"
            )
        if not self.coeffs:
            return ZERO_Q
        offset = k - self.val
        if offset < 0 or offset >= len(self.coeffs):
            return ZERO_Q
        return self.coeffs[offset]

    def __getitem__(self, k):
        return self.coeff(k)

    def truncate(self, order):
        new_cut = order + 1
        if self.cutoff is not None and self.cutoff < new_cut:
            new_cut = self.cutoff
        return QQJet._raw(self.val, self.coeffs, new_cut)

    def is_zero_to_order(self, order):
        """Verify cutoff >= order + 1 AND all coefficients up to eps^order are 0."""
        req_cut = order + 1
        if self.cutoff is not None and self.cutoff < req_cut:
            return False
        if not self.coeffs:
            return True
        return self.val >= req_cut

    def has_infinitesimal(self):
        """Return True iff the jet has any non-zero coefficient at exponent k >= 1."""
        if not self.coeffs:
            return False
        for idx, c in enumerate(self.coeffs):
            if self.val + idx >= 1 and not _coeff_is_zero(c):
                return True
        return False

    def _c0_fraction(self):
        """Extract exact Fraction c_0 at eps^0; raises ZeroDivisionError on Laurent poles."""
        if self.val < 0 and self.coeffs:
            raise ZeroDivisionError("QQJet has a Laurent pole at eps=0")
        c0 = self.coeff(0)
        if _is_exact_scalar(c0):
            return _to_fraction(c0)
        tier = getattr(c0, "_cas_tier", 0)
        if tier in (1, 2) and c0.is_constant():
            return c0.constant_value()
        raise TypeError(f"QQJet eps^0 coefficient {c0!r} is symbolic, not a concrete Fraction")

    @classmethod
    def _assert_not_scalar_half_ulp(cls, q0):
        """Assert rational q0 != 0 is not a dyadic half-ulp jump boundary for Q_dyad."""
        num = q0.numerator
        if num == 0:
            return
        den = q0.denominator
        # Half-ulp tie (2m+1)/2^(shift+1) requires denominator to be a power of 2 >= 2
        if den < 2 or (den & (den - 1)) != 0:
            return
        anum = -num if num < 0 else num
        eb = anum.bit_length() - den.bit_length()
        if eb >= 0:
            ef = (eb - 1) if (anum < (den << eb)) else eb
        else:
            ef = (eb - 1) if ((anum << (-eb)) < den) else eb
        for mb in cls._checked_mantissa_bits:
            shift = mb - ef
            if shift > cls._checked_max_shift:
                shift = cls._checked_max_shift
            if shift >= 0:
                if (((anum << (shift + 1)) + den) % (den << 1)) == 0:
                    raise ValueError(
                        f"QQJet eps->0 germ sits on dyadic half-ulp jump point {q0} (mantissa_bits={mb})"
                    )
            else:
                rshift = -shift
                den_shifted = den << rshift
                if (((anum << 1) + den_shifted) % (den_shifted << 1)) == 0:
                    raise ValueError(
                        f"QQJet eps->0 germ sits on dyadic half-ulp jump point {q0} (mantissa_bits={mb})"
                    )

    @property
    def numerator(self):
        """Expose eps->0 constant term numerator for piecewise-constant Q_dyad germ evaluation,
        asserting non-constant jets never sit on a dyadic half-ulp jump boundary.
        """
        q0 = self._c0_fraction()
        if self.has_infinitesimal():
            if QQJet._strict_half_ulp:
                QQJet._assert_not_scalar_half_ulp(q0)
            if QQJet._audit_jumps:
                QQJet._nonconst_int_queries.append(q0)
        return q0.numerator

    @property
    def denominator(self):
        """Expose eps->0 constant term denominator for piecewise-constant Q_dyad germ evaluation."""
        q0 = self._c0_fraction()
        if self.has_infinitesimal() and QQJet._strict_half_ulp:
            QQJet._assert_not_scalar_half_ulp(q0)
        return q0.denominator

    # --- Arithmetic Dunders ---

    def _coerce_other(self, other):
        tier = getattr(other, "_cas_tier", 0)
        if tier == 3:
            return other
        if tier in (1, 2) or _is_exact_scalar(other):
            return QQJet.const(other)
        return None

    def __bool__(self):
        return len(self.coeffs) > 0

    def __pos__(self):
        return self

    def __neg__(self):
        if not self.coeffs:
            return self
        return QQJet._raw(self.val, tuple(-c for c in self.coeffs), self.cutoff)

    @staticmethod
    def _min_cutoff(c1, c2):
        if c1 is None:
            return c2
        if c2 is None:
            return c1
        return c1 if c1 <= c2 else c2

    def __add__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        new_cut = QQJet._min_cutoff(self.cutoff, co.cutoff)
        if not self.coeffs and not co.coeffs:
            return QQJet._raw(0, (), new_cut)
        if not self.coeffs:
            return QQJet._raw(co.val, co.coeffs, new_cut)
        if not co.coeffs:
            return QQJet._raw(self.val, self.coeffs, new_cut)
        v_min = self.val if self.val <= co.val else co.val
        end1 = self.val + len(self.coeffs)
        end2 = co.val + len(co.coeffs)
        v_max = end1 if end1 >= end2 else end2
        if new_cut is not None and v_max > new_cut:
            v_max = new_cut
        if v_max <= v_min:
            return QQJet._raw(new_cut, (), new_cut)
        out = []
        for k in range(v_min, v_max):
            i1 = k - self.val
            c1 = self.coeffs[i1] if (0 <= i1 < len(self.coeffs)) else ZERO_Q
            i2 = k - co.val
            c2 = co.coeffs[i2] if (0 <= i2 < len(co.coeffs)) else ZERO_Q
            out.append(c1 + c2)
        return QQJet._raw(v_min, out, new_cut)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return self.__add__(-co)

    def __rsub__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__add__(-self)

    def __mul__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        # Exact zero handling
        if not self.coeffs and self.cutoff is None:
            return self
        if not co.coeffs and co.cutoff is None:
            return co
        if not self.coeffs and not co.coeffs:
            nc = self.cutoff + co.cutoff
            return QQJet._raw(nc, (), nc)
        if not self.coeffs:
            nc = self.cutoff + co.val
            return QQJet._raw(nc, (), nc)
        if not co.coeffs:
            nc = co.cutoff + self.val
            return QQJet._raw(nc, (), nc)

        # Both have non-empty coeffs starting at self.val and co.val
        c1_cut = (self.cutoff + co.val) if self.cutoff is not None else None
        c2_cut = (co.cutoff + self.val) if co.cutoff is not None else None
        new_cut = QQJet._min_cutoff(c1_cut, c2_cut)
        new_val = self.val + co.val
        n1 = len(self.coeffs)
        n2 = len(self.coeffs)
        prod_len = len(self.coeffs) + len(co.coeffs) - 1
        if new_cut is not None:
            cap = new_cut - new_val
            if cap <= 0:
                return QQJet._raw(new_cut, (), new_cut)
            if prod_len > cap:
                prod_len = cap
        out = [ZERO_Q] * prod_len
        for i, a in enumerate(self.coeffs):
            if i >= prod_len:
                break
            if _coeff_is_zero(a):
                continue
            j_lim = prod_len - i
            for j, b in enumerate(co.coeffs):
                if j >= j_lim:
                    break
                if _coeff_is_zero(b):
                    continue
                out[i + j] = out[i + j] + a * b
        return QQJet._raw(new_val, out, new_cut)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        if not co.coeffs:
            raise ZeroDivisionError("Division by zero QQJet")
        if not self.coeffs:
            if self.cutoff is None:
                return self
            nc = self.cutoff - co.val
            return QQJet._raw(nc, (), nc)

        new_val = self.val - co.val
        # Exact scalar denominator fast path
        if co.cutoff is None and len(co.coeffs) == 1:
            inv_b0 = ONE_Q / co.coeffs[0]
            new_cut = (self.cutoff - co.val) if self.cutoff is not None else None
            return QQJet._raw(new_val, tuple(a * inv_b0 for a in self.coeffs), new_cut)

        p_a = (self.cutoff - self.val) if self.cutoff is not None else None
        p_b = (co.cutoff - co.val) if co.cutoff is not None else None
        rel_prec = QQJet._min_cutoff(p_a, p_b)
        if rel_prec is None:
            rel_prec = 8  # Fallback cap only if someone constructs an un-truncated series
        if rel_prec <= 0:
            nc = new_val + rel_prec
            return QQJet._raw(nc, (), nc)

        b0 = co.coeffs[0]
        inv_b0 = ONE_Q / b0
        q = []
        for k in range(rel_prec):
            acc = self.coeffs[k] if k < len(self.coeffs) else ZERO_Q
            j_max = k if k < len(co.coeffs) - 1 else len(co.coeffs) - 1
            for j in range(1, j_max + 1):
                acc = acc - co.coeffs[j] * q[k - j]
            q.append(acc * inv_b0)
        return QQJet._raw(new_val, q, new_val + rel_prec)

    def __rtruediv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__truediv__(self)

    def __pow__(self, exp):
        if not _is_exact_scalar(exp) or _to_fraction(exp).denominator != 1:
            return NotImplemented
        n = int(_to_fraction(exp).numerator)
        if n == 0:
            return QQJet.const(ONE_Q)
        if n < 0:
            return QQJet.const(ONE_Q) / (self ** (-n))
        res = QQJet.const(ONE_Q)
        base = self
        e = n
        while e > 0:
            if e & 1:
                res = res * base
            e >>= 1
            if e > 0:
                base = base * base
        return res

    def __abs__(self):
        q0 = self._c0_fraction()
        if q0 > ZERO_Q:
            return self
        if q0 < ZERO_Q:
            return -self
        if not self.has_infinitesimal():
            return self
        raise ValueError("abs(QQJet) is non-differentiable at c_0 == 0 with non-zero infinitesimal part")

    def _cmp_diff_c0(self, other):
        co = self._coerce_other(other)
        if co is None:
            return None
        diff = self - co
        d0 = diff._c0_fraction()
        if d0 == ZERO_Q and diff.has_infinitesimal():
            raise ValueError(
                f"QQJet order comparison sits on exact branch-cut tie c_0 == 0 for non-constant diff {diff}"
            )
        return d0

    def __lt__(self, other):
        d0 = self._cmp_diff_c0(other)
        return NotImplemented if d0 is None else (d0 < ZERO_Q)

    def __le__(self, other):
        d0 = self._cmp_diff_c0(other)
        return NotImplemented if d0 is None else (d0 <= ZERO_Q)

    def __gt__(self, other):
        d0 = self._cmp_diff_c0(other)
        return NotImplemented if d0 is None else (d0 > ZERO_Q)

    def __ge__(self, other):
        d0 = self._cmp_diff_c0(other)
        return NotImplemented if d0 is None else (d0 >= ZERO_Q)

    def __floordiv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        a0 = self._c0_fraction()
        b0 = co._c0_fraction()
        if b0 == ZERO_Q:
            raise ZeroDivisionError("QQJet floor division by zero at eps=0")
        ratio = a0 / b0
        if (self.has_infinitesimal() or co.has_infinitesimal()) and ratio.denominator == 1:
            raise ValueError(
                f"QQJet floor/mod sits on exact integer jump discontinuity {ratio} at eps=0"
            )
        return a0 // b0

    def __rfloordiv__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__floordiv__(self)

    def __mod__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        q_int = self.__floordiv__(co)
        return self - q_int * co

    def __rmod__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        return co.__mod__(self)

    def __eq__(self, other):
        co = self._coerce_other(other)
        if co is None:
            return NotImplemented
        diff = self - co
        if diff.cutoff is not None and diff.cutoff <= 0:
            return False
        return len(diff.coeffs) == 0

    def __repr__(self):
        parts = []
        for idx, c in enumerate(self.coeffs):
            if _coeff_is_zero(c):
                continue
            k = self.val + idx
            if k == 0:
                parts.append(f"{c}")
            elif k == 1:
                parts.append(f"({c})*eps")
            else:
                parts.append(f"({c})*eps^{k}")
        if not parts:
            parts.append("0")
        if self.cutoff is not None:
            parts.append(f"O(eps^{self.cutoff})")
        return " + ".join(parts)


# ============================================================================
# Ergonomic Constructor Helpers
# ============================================================================


def poly_vars(*names):
    """Create a tuple of QQPoly variables: `x, y, z = poly_vars('x', 'y', 'z')`."""
    if len(names) == 1 and " " in names[0]:
        names = tuple(names[0].split())
    elif len(names) == 1 and "," in names[0]:
        names = tuple(s.strip() for s in names[0].split(",") if s.strip())
    if len(names) == 1:
        return QQPoly.var(names[0])
    return tuple(QQPoly.var(n) for n in names)


def qqpoly_vars(prefix, count):
    return [QQPoly.var(f"{prefix}_{i}") for i in range(count)]


def rat_vars(*names):
    """Create a tuple of QQRat variables: `x, y, z = rat_vars('x', 'y', 'z')`."""
    if len(names) == 1 and " " in names[0]:
        names = tuple(names[0].split())
    elif len(names) == 1 and "," in names[0]:
        names = tuple(s.strip() for s in names[0].split(",") if s.strip())
    if len(names) == 1:
        return QQRat.var(names[0])
    return tuple(QQRat.var(n) for n in names)


def qqrat_vars(prefix, count):
    return [QQRat.var(f"{prefix}_{i}") for i in range(count)]


def jet_var(order=4):
    return QQJet.var(order=order)
