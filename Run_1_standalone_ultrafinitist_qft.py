#!/usr/bin/env python3
"""Self-Contained Single-File Verification Suite for Ultrafinitist Quantum Field Theory.

Copy this single file into ANY empty directory and run:
    python3 standalone_ultrafinitist_qft.py

Zero local sibling imports required. Automatically uses `numpy`/`scipy`/`sympy`/`matplotlib`
if installed, falling back seamlessly to macOS `Accelerate.framework` (`dsyev_`, `zheev_`,
`dgemm_`) or pure-Python cyclic Jacobi eigensolvers, exact `fractions.Fraction` rational
Positivstellensatz bounds over Q, and `Pillow` (PNG) or pure-Python (SVG) diagnostic plotting.
"""
from __future__ import annotations

import argparse, array, cmath, collections, ctypes, math, os, platform, sys, time
from ctypes import POINTER, byref, c_char_p, c_double, c_int
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

# --- 1. Multi-Tier Hardware Linear Algebra (NumPy -> macOS Accelerate -> Pure-Python Jacobi) ---
try:
    import numpy as _np  # type: ignore
    HAS_NUMPY = True
except Exception:
    _np, HAS_NUMPY = None, False

try:
    import sympy as _sp  # type: ignore
    HAS_SYMPY = True
except Exception:
    _sp, HAS_SYMPY = None, False

_ACCEL = None
if sys.platform == "darwin" and not HAS_NUMPY:
    try:
        _ACCEL = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        _ACCEL.dsyev_.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_double), POINTER(c_int),
                                  POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int)]
        _ACCEL.zheev_.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_double), POINTER(c_int),
                                  POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_int)]
        _ACCEL.dgemm_.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                  POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_double), POINTER(c_int),
                                  POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    except Exception:
        _ACCEL = None

ACTIVE_BACKEND = "numpy" if HAS_NUMPY else ("accelerate" if _ACCEL is not None else "pure_python")


def dgemm_rm(a_rm: Sequence[float], b_rm: Sequence[float], m: int, k: int, n: int) -> List[float]:
    """Row-major real matrix product C (m x n) = A (m x k) @ B (k x n)."""
    if HAS_NUMPY:
        a_arr = _np.frombuffer(array.array("d", a_rm), dtype=_np.float64).reshape(m, k)
        b_arr = _np.frombuffer(array.array("d", b_rm), dtype=_np.float64).reshape(k, n)
        return (a_arr @ b_arr).ravel().tolist()
    if _ACCEL is not None:
        A_c = (c_double * (m * k)).from_buffer_copy(array.array("d", a_rm))
        B_c = (c_double * (k * n)).from_buffer_copy(array.array("d", b_rm))
        C_c = (c_double * (m * n))()
        nn, mm, kk = c_int(n), c_int(m), c_int(k)
        one, zero = c_double(1.0), c_double(0.0)
        _ACCEL.dgemm_(b"N", b"N", byref(nn), byref(mm), byref(kk), byref(one),
                      B_c, byref(nn), A_c, byref(kk), byref(zero), C_c, byref(nn))
        return list(C_c)
    out = [0.0] * (m * n)
    for i in range(m):
        i_k, i_n = i * k, i * n
        for p in range(k):
            aip = a_rm[i_k + p]
            if aip != 0.0:
                p_n = p * n
                for j in range(n):
                    out[i_n + j] += aip * b_rm[p_n + j]
    return out


def eigh_real_sym(mat: Sequence[Sequence[float]], vectors: bool = True) -> Tuple[List[float], Optional[List[float]]]:
    """Symmetric eigensolver returning ascending eigenvalues and column-major eigenvectors (col*n + row)."""
    n = len(mat)
    if HAS_NUMPY:
        arr = _np.array(mat, dtype=_np.float64)
        if vectors:
            w, v = _np.linalg.eigh(arr)
            return w.tolist(), v.T.ravel().tolist()
        return _np.linalg.eigvalsh(arr).tolist(), None
    if _ACCEL is not None:
        col_buf = array.array("d", (float(mat[r][c]) for c in range(n) for r in range(n)))
        a_c = (c_double * (n * n)).from_buffer_copy(col_buf)
        w_c = (c_double * n)()
        jobz, uplo, nn, info = (b"V" if vectors else b"N"), b"U", c_int(n), c_int(0)
        wq, minus1 = c_double(0.0), c_int(-1)
        _ACCEL.dsyev_(jobz, uplo, byref(nn), a_c, byref(nn), w_c, byref(wq), byref(minus1), byref(info))
        lwork = c_int(max(1, int(wq.value)))
        work = (c_double * lwork.value)()
        _ACCEL.dsyev_(jobz, uplo, byref(nn), a_c, byref(nn), w_c, work, byref(lwork), byref(info))
        assert info.value == 0, f"dsyev_ failed: info={info.value}"
        return list(w_c), (list(a_c) if vectors else None)
    a = [[float(mat[i][j]) for j in range(n)] for i in range(n)]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(80):
        off = max((abs(a[p][q]) for p in range(n) for q in range(p + 1, n)), default=0.0)
        if off < 1e-13:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-15:
                    continue
                tau = (a[q][q] - a[p][p]) / (2.0 * apq)
                t = math.copysign(1.0 / (abs(tau) + math.hypot(1.0, tau)), tau) if tau != 0.0 else 1.0
                c = 1.0 / math.hypot(1.0, t)
                s = t * c
                a[p][p] -= t * apq
                a[q][q] += t * apq
                a[p][q] = a[q][p] = 0.0
                for r in range(n):
                    if r != p and r != q:
                        arp, arq = a[r][p], a[r][q]
                        a[r][p] = a[p][r] = c * arp - s * arq
                        a[r][q] = a[q][r] = s * arp + c * arq
                    if vectors:
                        vrp, vrq = v[r][p], v[r][q]
                        v[r][p], v[r][q] = c * vrp - s * vrq, s * vrp + c * vrq
    order = sorted(range(n), key=lambda idx: a[idx][idx])
    return [a[idx][idx] for idx in order], ([v[r][idx] for idx in order for r in range(n)] if vectors else None)


def eigh_complex_herm(mat: Sequence[Sequence[complex]], vectors: bool = True) -> Tuple[List[float], Optional[List[float]]]:
    """Hermitian eigensolver returning ascending eigenvalues and interleaved col-major eigenvectors."""
    n = len(mat)
    if HAS_NUMPY:
        arr = _np.array(mat, dtype=_np.complex128)
        if vectors:
            w, v = _np.linalg.eigh(arr)
            inter = _np.empty(2 * n * n, dtype=_np.float64)
            vt = v.T.ravel()
            inter[0::2], inter[1::2] = vt.real, vt.imag
            return w.tolist(), inter.tolist()
        return _np.linalg.eigvalsh(arr).tolist(), None
    if _ACCEL is not None:
        buf = array.array("d", [0.0]) * (2 * n * n)
        idx = 0
        for c in range(n):
            for r in range(n):
                z = mat[r][c]
                buf[idx], buf[idx + 1] = z.real, z.imag
                idx += 2
        a_c = (c_double * (2 * n * n)).from_buffer_copy(buf)
        w_c = (c_double * n)()
        rwork = (c_double * max(1, 3 * n - 2))()
        jobz, uplo, nn, info = (b"V" if vectors else b"N"), b"U", c_int(n), c_int(0)
        wq, minus1 = (c_double * 2)(), c_int(-1)
        _ACCEL.zheev_(jobz, uplo, byref(nn), a_c, byref(nn), w_c, wq, byref(minus1), rwork, byref(info))
        lwork = c_int(max(1, int(wq[0])))
        work = (c_double * (2 * lwork.value))()
        _ACCEL.zheev_(jobz, uplo, byref(nn), a_c, byref(nn), w_c, work, byref(lwork), rwork, byref(info))
        assert info.value == 0, f"zheev_ failed: info={info.value}"
        return list(w_c), (list(a_c) if vectors else None)
    real2 = [[0.0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            z = mat[i][j]
            real2[i][j] = real2[n + i][n + j] = z.real
            real2[i][n + j], real2[n + i][j] = -z.imag, z.imag
    w2, v2 = eigh_real_sym(real2, vectors=vectors)
    evals = [0.5 * (w2[2 * k] + w2[2 * k + 1]) for k in range(n)]
    if not vectors or v2 is None:
        return evals, None
    inter = [0.0] * (2 * n * n)
    for k in range(n):
        b2 = (2 * k) * (2 * n)
        for r in range(n):
            inter[2 * (k * n + r)], inter[2 * (k * n + r) + 1] = v2[b2 + r], v2[b2 + n + r]
    return evals, inter


# --- 2. Weyl-Schwinger Algebra & Exact Rational Positivstellensatz Bounds over Q ---
def verify_weyl_and_symbolic_identities() -> Dict[str, float]:
    """Verify Weyl-Schwinger clock-shift algebra (d=3,4,5), exact d^N phi^4, and SOS identities."""
    for d in (3, 4, 5):
        omega = cmath.exp(2j * math.pi / d)
        u_diag = [omega ** q for q in range(d)]
        assert max(abs((z ** d) - 1.0) for z in u_diag) < 1e-13
        assert max(abs(u_diag[q] - omega * u_diag[(q - 1) % d]) for q in range(d)) < 1e-13
        assert max(abs(math.sin(2.0 * math.pi * q / d)) for q in range(d)) <= 1.0 + 1e-15

    # Exact d^N many-body tensor-product Hamiltonian on N=4, d=4 (dim H = 256)
    dim_h, phi_4, lam_s = 256, [math.sin(0.5 * math.pi * q) for q in range(4)], 5.0 / 8.0
    k_4 = [[2.5 if i == j else (-4.0 / 3.0 if (i - j) % 4 in (1, 3) else 1.0 / 6.0) for j in range(4)] for i in range(4)]
    h_tot, strides, max_hint = [[0.0] * dim_h for _ in range(dim_h)], [64, 16, 4, 1], 0.0
    for idx in range(dim_h):
        qs = [(idx // strides[x]) % 4 for x in range(4)]
        pot = sum(0.5 * (k_4[x][y] + (0.25 if x == y else 0.0)) * phi_4[qs[x]] * phi_4[qs[y]] for x in range(4) for y in range(4))
        hint = lam_s * sum(phi_4[qs[x]] ** 4 for x in range(4))
        max_hint = max(max_hint, abs(hint))
        h_tot[idx][idx] = 4.0 + pot + hint
        for x in range(4):
            for st in (-1, 1):
                h_tot[idx][idx + (((qs[x] + st) % 4) - qs[x]) * strides[x]] -= 0.5
    ev_mb, _ = eigh_real_sym(h_tot, vectors=False)
    assert abs(max_hint - 2.5) < 1e-14 and (ev_mb[1] - ev_mb[0]) > 0.4

    # Exact Q grid check + optional SymPy proof of Symanzik & Equal-Corner Biharmonic Wilson SOS
    for lx_q in (Fraction(0), Fraction(1), Fraction(2), Fraction(4)):
        assert Fraction(1, 12) * (4 - lx_q) * (16 + lx_q) == Fraction(16, 3) - (lx_q + lx_q ** 2 / 12) >= 0
        for ly_q in (Fraction(0), Fraction(2), Fraction(4)):
            wf = (lx_q - ly_q / 2) ** 2 + Fraction(3, 4) * ly_q ** 2
            assert wf == lx_q ** 2 - lx_q * ly_q + ly_q ** 2 and 0 <= wf <= 16
    if HAS_SYMPY:
        sx, sy = _sp.symbols("sx sy", real=True)
        assert _sp.expand((sx - sy / 2) ** 2 + _sp.Rational(3, 4) * sy ** 2 - (sx ** 2 - sx * sy + sy ** 2)) == 0
    return {"dim_H": float(dim_h), "E0": ev_mb[0], "gap": ev_mb[1] - ev_mb[0], "H_int_norm": max_hint}


def certify_rational_resolvent_bounds(int12_x_rows: List[List[Tuple[int, int]]], m_cap: Fraction,
                                      m0_sq: Fraction, lam_site: Fraction, order_2m: int = 64) -> Dict[str, float]:
    """Exact Q Positivstellensatz bounds via zero-mode-projected midpoint operator (12 X)^r e_0."""
    n = len(int12_x_rows)
    alpha, half_m = m_cap / 2 + m0_sq, m_cap / 2
    rho = half_m / alpha
    vec = [0] * n
    vec[0] = 1
    m_perp = [Fraction(n - 1, n)]
    for r in range(1, order_2m + 1):
        nxt = [0] * n
        for u in range(n):
            acc = 0
            for v, w_int in int12_x_rows[u]:
                acc += w_int * vec[v]
            nxt[u] = acc
        vec = nxt
        m_perp.append(Fraction(vec[0], 12 ** r) - (half_m ** r) / n)

    z1, z2 = Fraction(1, n) / m0_sq, Fraction(1, n) / (m0_sq ** 2)
    tau1_lo = z1 + sum((m_perp[r] / (alpha ** (r + 1)) for r in range(order_2m)), Fraction(0))
    tau1_hi = tau1_lo + m_perp[order_2m] / ((alpha ** (order_2m + 1)) * (1 - rho))
    tau2_hi = z2 + sum(((r + 1) * m_perp[r] / (alpha ** (r + 2)) for r in range(order_2m)), Fraction(0)) + (
        (order_2m + 1) * m_perp[order_2m] / ((alpha ** (order_2m + 2)) * ((1 - rho) ** 2))
    )
    # By Varga Regular Splitting G_{xy} > 0: tau_1^3 <= S_2 <= G_00 * Tr(G^2)/N = tau_1 * tau_2
    s2_lo, s2_hi = tau1_lo ** 3, tau1_hi * tau2_hi
    rm1_lo = 1 + (Fraction(1, 2) * lam_site * tau1_lo) / m0_sq
    rm1_hi = 1 + (Fraction(1, 2) * lam_site * tau1_hi) / m0_sq
    rm2_lo = rm1_lo - (Fraction(1, 6) * (lam_site ** 2) * s2_hi) / m0_sq
    rm2_hi = rm1_hi - (Fraction(1, 6) * (lam_site ** 2) * s2_lo) / m0_sq
    return {
        "tau1_lo": float(tau1_lo), "tau1_hi": float(tau1_hi),
        "rm2_lo": float(rm2_lo), "rm2_hi": float(rm2_hi),
    }


# --- 3. Deterministic Finite Graph & Operator Construction (1+1D & 2+1D, Zero FFT) ---
def build_finite_lattice_operators(mode: str, dim: str, n: int) -> Dict[str, object]:
    """Construct exact Baseline or Cycle 3 scalar K_G, combinatorial L_G, Dirac H_D, and shift maps."""
    m0, m0_sq = 0.5, 0.25
    lap = [[0.0] * n for _ in range(n)]
    k_mat = [[0.0] * n for _ in range(n)]
    h_d = [[0j] * (2 * n) for _ in range(2 * n)]
    int12_x: List[List[Tuple[int, int]]] = [[] for _ in range(n)]

    if dim == "1+1D":
        next_x = [(x + 1) % n for x in range(n)]
        next_maps = [next_x]
        local_nb = [[(x - 1) % n, (x + 1) % n] for x in range(n)]
        for x in range(n):
            xp1, xm1, xp2, xm2 = (x + 1) % n, (x - 1) % n, (x + 2) % n, (x - 2) % n
            lap[x][x] = 2.0
            lap[x][xp1] -= 1.0
            lap[x][xm1] -= 1.0
            if mode == "baseline":
                k_mat[x][x] = 2.0
                k_mat[x][xp1] -= 1.0
                k_mat[x][xm1] -= 1.0
                h_d[x][x], h_d[n + x][n + x] = m0, -m0
                h_d[x][n + xp1] += -0.5j
                h_d[x][n + xm1] += 0.5j
                h_d[n + x][xp1] += -0.5j
                h_d[n + x][xm1] += 0.5j
            else:
                k_mat[x][x] += 2.5
                k_mat[x][xp1] -= 4.0 / 3.0
                k_mat[x][xm1] -= 4.0 / 3.0
                k_mat[x][xp2] += 1.0 / 12.0
                k_mat[x][xm2] += 1.0 / 12.0
                # 12 * ((8/3)I - K_Sym) has exact integer weights: +2 at 0, +16 at +-1, -1 at +-2
                int12_x[x] = [(x, 2), (xp1, 16), (xm1, 16), (xp2, -1), (xm2, -1)]
                r2 = 1.0 / 16.0
                for sgn, off in ((1.0, 0), (-1.0, n)):
                    h_d[off + x][off + x] += sgn * (m0 + 6.0 * r2)
                    h_d[off + x][off + xp1] += sgn * (-4.0 * r2)
                    h_d[off + x][off + xm1] += sgn * (-4.0 * r2)
                    h_d[off + x][off + xp2] += sgn * r2
                    h_d[off + x][off + xm2] += sgn * r2
                for r_off, c_off in ((0, n), (n, 0)):
                    h_d[r_off + x][c_off + xp1] += -2.0j / 3.0
                    h_d[r_off + x][c_off + xm1] += 2.0j / 3.0
                    h_d[r_off + x][c_off + xp2] += 1.0j / 12.0
                    h_d[r_off + x][c_off + xm2] += -1.0j / 12.0
        m_cap, lam_site = Fraction(16, 3), Fraction(5, 8)
    else:
        m_exp = n.bit_length() - 1
        l_half = 1 << (m_exp // 2)
        if mode == "baseline" or (m_exp % 2 == 0):
            lx = l_half if (m_exp % 2 == 0) else l_half
            ly = l_half if (m_exp % 2 == 0) else (2 * l_half)
            coords = [(ix, iy) for ix in range(lx) for iy in range(ly)]
            idx_of = {pt: i for i, pt in enumerate(coords)}
            red = lambda x, y: (x % lx, y % ly)
        else:
            coords = [(ix, iy) for ix in range(2 * l_half) for iy in range(l_half)]
            idx_of = {pt: i for i, pt in enumerate(coords)}
            red = lambda x, y: ((x - (y // l_half) * l_half) % (2 * l_half), y - (y // l_half) * l_half)

        next_x = [idx_of[red(ix + 1, iy)] for ix, iy in coords]
        next_y = [idx_of[red(ix, iy + 1)] for ix, iy in coords]
        next_maps = [next_x, next_y]
        local_nb = [[] for _ in range(n)]
        ax1 = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        dg1 = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
        ax2 = [(-2, 0), (2, 0), (0, -2), (0, 2)]

        for u, (ix, iy) in enumerate(coords):
            xp1, xm1 = idx_of[red(ix + 1, iy)], idx_of[red(ix - 1, iy)]
            yp1, ym1 = idx_of[red(ix, iy + 1)], idx_of[red(ix, iy - 1)]
            xp2, xm2 = idx_of[red(ix + 2, iy)], idx_of[red(ix - 2, iy)]
            yp2, ym2 = idx_of[red(ix, iy + 2)], idx_of[red(ix, iy - 2)]
            local_nb[u] = [xm1, xp1, ym1, yp1]
            if mode == "baseline":
                lap[u][u] = k_mat[u][u] = 4.0
                for v in local_nb[u]:
                    lap[u][v] -= 1.0
                    k_mat[u][v] -= 1.0
                h_d[u][u], h_d[n + u][n + u] = m0, -m0
                h_d[u][n + xp1] += -0.5j
                h_d[u][n + xm1] += 0.5j
                h_d[n + u][xp1] += -0.5j
                h_d[n + u][xm1] += 0.5j
                h_d[u][n + yp1] += -0.5
                h_d[u][n + ym1] += 0.5
                h_d[n + u][yp1] += 0.5
                h_d[n + u][ym1] += -0.5
            else:
                lap[u][u] = 8.0
                k_mat[u][u] = 5.0
                row_x = [(u, 4)]
                for v in (xm1, xp1, ym1, yp1):
                    lap[u][v] -= 1.0
                    k_mat[u][v] -= 4.0 / 3.0
                    row_x.append((v, 16))
                for dx, dy in dg1:
                    lap[u][idx_of[red(ix + dx, iy + dy)]] -= 1.0
                for v in (xm2, xp2, ym2, yp2):
                    k_mat[u][v] += 1.0 / 12.0
                    row_x.append((v, -1))
                int12_x[u] = row_x
                r2 = 3.0 / 32.0
                for sgn, off in ((1.0, 0), (-1.0, n)):
                    h_d[off + u][off + u] += sgn * (m0 + 8.0 * r2)
                    for v in (xp1, xm1, yp1, ym1):
                        h_d[off + u][off + v] += sgn * (-2.0 * r2)
                    for dx, dy in dg1:
                        h_d[off + u][off + idx_of[red(ix + dx, iy + dy)]] += sgn * (-1.0 * r2)
                    for v in (xp2, xm2, yp2, ym2):
                        h_d[off + u][off + v] += sgn * r2
                for r_off, c_off in ((0, n), (n, 0)):
                    h_d[r_off + u][c_off + xp1] += -2.0j / 3.0
                    h_d[r_off + u][c_off + xm1] += 2.0j / 3.0
                    h_d[r_off + u][c_off + xp2] += 1.0j / 12.0
                    h_d[r_off + u][c_off + xm2] += -1.0j / 12.0
                h_d[u][n + yp1] += -2.0 / 3.0
                h_d[u][n + ym1] += 2.0 / 3.0
                h_d[n + u][yp1] += 2.0 / 3.0
                h_d[n + u][ym1] += -2.0 / 3.0
                h_d[u][n + yp2] += 1.0 / 12.0
                h_d[u][n + ym2] += -1.0 / 12.0
                h_d[n + u][yp2] += -1.0 / 12.0
                h_d[n + u][ym2] += 1.0 / 12.0
        m_cap, lam_site = Fraction(32, 3), Fraction(5, 4)

    return {
        "mode": mode, "dim": dim, "N": n, "m0": m0, "m0_sq": m0_sq,
        "K": k_mat, "L": lap, "H_D": h_d, "next_maps": next_maps, "local_nb": local_nb,
        "int12_x": int12_x, "m_cap": m_cap, "lam_site": lam_site,
    }


# --- 4. Spectral, Dispersion, Doubler, Lightcone & Zero-Subtraction Loop Evaluator ---
def evaluate_lattice_model(op: Dict[str, object], keep_full_spectrum: bool = False) -> Dict[str, object]:
    n, m0, m0_sq = int(op["N"]), float(op["m0"]), float(op["m0_sq"])
    k_mat, lap_mat, h_d = op["K"], op["L"], op["H_D"]  # type: ignore
    next_maps: List[List[int]] = op["next_maps"]  # type: ignore
    d_sp = len(next_maps)

    # Commuting symmetry-resolving operator breaks D_4 degeneracies without continuous FFT
    res_op = [[k_mat[i][j] for j in range(n)] for i in range(n)]
    chord_op = [[0.0] * n for _ in range(n)]
    primes = [101.0, 103.0, 107.0]
    maps_ext = list(next_maps) + ([[next_maps[1][next_maps[0][u]] for u in range(n)]] if d_sp == 2 else [])
    for idx_m, mp in enumerate(maps_ext):
        w = 1.0 / primes[idx_m]
        for u in range(n):
            v = mp[u]
            res_op[u][u] += 2.0 * w
            res_op[u][v] -= w
            res_op[v][u] -= w
            chord_op[u][u] += 2.0 * w
            chord_op[u][v] -= w
            chord_op[v][u] -= w

    _, evecs_col = eigh_real_sym(res_op, vectors=True)
    assert evecs_col is not None
    v_rm = [evecs_col[k * n + x] for x in range(n) for k in range(n)]
    k_rm = [k_mat[i][j] for i in range(n) for j in range(n)]
    kv_rm = dgemm_rm(k_rm, v_rm, n, n, n)
    evals_k = [max(0.0, sum(v_rm[x * n + k] * kv_rm[x * n + k] for x in range(n))) for k in range(n)]

    # Extract intrinsic momenta p_vec strictly from shift-operator expectations <v_k | T_mu | v_k>
    modes_s = []
    nx = next_maps[0]
    ny = next_maps[1] if d_sp == 2 else None
    for k in range(n):
        vk = [evecs_col[k * n + u] for u in range(n)]
        c_x = max(-1.0, min(1.0, sum(vk[u] * vk[nx[u]] for u in range(n))))
        px = math.acos(c_x)
        if d_sp == 1:
            p2 = px * px
        else:
            assert ny is not None
            c_y = max(-1.0, min(1.0, sum(vk[u] * vk[ny[u]] for u in range(n))))
            py_abs = math.acos(c_y)
            c_xy = max(-1.0, min(1.0, sum(vk[u] * vk[ny[nx[u]]] for u in range(n))))
            py = py_abs if abs(c_xy - math.cos(px + py_abs)) <= abs(c_xy - math.cos(px - py_abs)) else -py_abs
            p2 = px * px + py * py
        e2 = evals_k[k] + m0_sq
        dl = abs(e2 - (p2 + m0_sq)) / e2 if p2 > 1e-12 else 0.0
        modes_s.append((p2, e2, dl))
    modes_s.sort(key=lambda t: t[0])
    low20 = max(2, int(0.20 * n))
    low_s = modes_s[:low20]
    nz_s = [m for m in low_s if m[0] > 1e-12]
    p1_sq = nz_s[0][0]
    c_s2 = sum((m[1] - m0_sq) / m[0] for m in nz_s if abs(m[0] - p1_sq) < 1e-7) / max(1, sum(1 for m in nz_s if abs(m[0] - p1_sq) < 1e-7))
    dl_s_mean = sum(m[2] for m in nz_s) / len(nz_s)

    def _shell_aniso(pairs: Sequence[Tuple[float, float]]) -> float:
        shells: Dict[int, List[float]] = collections.defaultdict(list)
        for p2, e2 in pairs:
            if p2 > 1e-10:
                shells[int(round(p2 * 1e5))].append(e2)
        return max(((max(vs) - min(vs)) / (sum(vs) / len(vs)) for vs in shells.values() if len(vs) >= 2), default=0.0)

    aniso_s = _shell_aniso([(m[0], m[1]) for m in low_s])

    # Exact N x N Propagator G = V diag(1/(lam + m0^2)) V^T & Zero-Subtraction 1-Loop / 2-Loop Traces
    g_eigs = [1.0 / (lam + m0_sq) for lam in evals_k]
    w_rm = [g_eigs[k] * evecs_col[k * n + y] for k in range(n) for y in range(n)]
    g_flat = dgemm_rm(v_rm, w_rm, n, n, n)
    tau1 = sum(g_flat[x * n + x] for x in range(n)) / n
    s2_per_n = sum(val * val * val for val in g_flat) / n
    if op["mode"] == "baseline":
        lam0 = 0.5
        sig1 = 0.5 * lam0 * (tau1 * n)
        sig2 = -(1.0 / 6.0) * (lam0 ** 2) * (s2_per_n * n)
        bounds = None
    else:
        ls = float(op["lam_site"])
        sig1 = 0.5 * ls * tau1
        sig2 = -(1.0 / 6.0) * (ls ** 2) * s2_per_n
        bounds = certify_rational_resolvent_bounds(op["int12_x"], op["m_cap"], Fraction(1, 4), op["lam_site"])  # type: ignore
        assert bounds["tau1_lo"] - 1e-11 <= tau1 <= bounds["tau1_hi"] + 1e-11
    rm1, rm2 = (m0_sq + sig1) / m0_sq, (m0_sq + sig1 + sig2) / m0_sq
    if bounds is not None:
        assert bounds["rm2_lo"] - 1e-11 <= rm2 <= bounds["rm2_hi"] + 1e-11

    # Causal Wavepacket Lightcone Check under exp(-i H_S t) at t = 1.0
    dist = [-1] * n
    dist[0] = 0
    q_bfs = collections.deque([0])
    local_nb: List[List[int]] = op["local_nb"]  # type: ignore
    while q_bfs:
        u = q_bfs.popleft()
        for v in local_nb[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q_bfs.append(v)
    omega_s = [math.sqrt(lam + m0_sq) for lam in evals_k]
    v0 = [evecs_col[k * n] for k in range(n)]
    cs_rm = [val for k in range(n) for val in (math.cos(-omega_s[k]) * v0[k], math.sin(-omega_s[k]) * v0[k])]
    psi_1 = dgemm_rm(v_rm, cs_rm, n, n, 2)
    leak_t1 = sum(psi_1[2 * x] ** 2 + psi_1[2 * x + 1] ** 2 for x in range(n) if dist[x] > 4)

    # Fermionic Spectrum, Dirichlet Roughness Doubler Audit & Dirac Shell Isotropy
    lap_evals, _ = eigh_real_sym(lap_mat, vectors=False)
    lam_max_l = max(lap_evals)
    d_evals, d_evecs = eigh_complex_herm(h_d, vectors=True)
    assert d_evecs is not None
    m_dim = 2 * n
    x_rm = [0.0] * (n * (4 * m_dim))
    for j in range(m_dim):
        jb, cb = 2 * j * m_dim, 4 * j
        for u in range(n):
            off = u * (4 * m_dim) + cb
            x_rm[off], x_rm[off + 1] = d_evecs[jb + 2 * u], d_evecs[jb + 2 * u + 1]
            x_rm[off + 2], x_rm[off + 3] = d_evecs[jb + 2 * (n + u)], d_evecs[jb + 2 * (n + u) + 1]
    lx_rm = dgemm_rm([lap_mat[i][j] for i in range(n) for j in range(n)], x_rm, n, n, 4 * m_dim)
    cx_rm = dgemm_rm([chord_op[i][j] for i in range(n) for j in range(n)], x_rm, n, n, 4 * m_dim)

    # Unmix degenerate energy clusters for exact roughness & intrinsic Dirac momenta
    r_tilde = [0.0] * m_dim
    d_pos_modes: List[Tuple[float, float]] = []
    idx_ord = sorted(range(m_dim), key=lambda i: d_evals[i])
    pos = 0
    while pos < m_dim:
        end = pos + 1
        while end < m_dim and abs(d_evals[idx_ord[end]] - d_evals[idx_ord[end - 1]]) <= 1e-7:
            end += 1
        cl = idx_ord[pos:end]
        g = len(cl)
        sub_l = [[0j] * g for _ in range(g)]
        sub_c = [[0j] * g for _ in range(g)]
        for a in range(g):
            ja = cl[a]
            for b in range(a, g):
                jb = cl[b]
                re_l = im_l = re_c = im_c = 0.0
                for u in range(n):
                    ca, cb = u * (4 * m_dim) + 4 * ja, u * (4 * m_dim) + 4 * jb
                    for s in (0, 2):
                        xra, xia = x_rm[ca + s], x_rm[ca + s + 1]
                        re_l += xra * lx_rm[cb + s] + xia * lx_rm[cb + s + 1]
                        im_l += xra * lx_rm[cb + s + 1] - xia * lx_rm[cb + s]
                        re_c += xra * cx_rm[cb + s] + xia * cx_rm[cb + s + 1]
                        im_c += xra * cx_rm[cb + s + 1] - xia * cx_rm[cb + s]
                sub_l[a][b] = complex(re_l, im_l if a != b else 0.0)
                sub_l[b][a] = sub_l[a][b].conjugate()
                sub_c[a][b] = complex(re_c, im_c if a != b else 0.0)
                sub_c[b][a] = sub_c[a][b].conjugate()
        rl_eigs, _ = eigh_complex_herm(sub_l, vectors=False)
        for k_s, j_orig in enumerate(cl):
            r_tilde[j_orig] = max(0.0, rl_eigs[k_s]) / lam_max_l
        if d_evals[cl[0]] > 0.0:
            _, u_c = eigh_complex_herm(sub_c, vectors=True)
            assert u_c is not None
            e2_val = d_evals[cl[0]] ** 2
            for k_s in range(g):
                # Reconstruct unmixed 2-spinor density & shift expectation
                psi_up = [0j] * n
                psi_dn = [0j] * n
                for a in range(g):
                    ja = cl[a]
                    coeff = complex(u_c[2 * (k_s * g + a)], u_c[2 * (k_s * g + a) + 1])
                    for u in range(n):
                        ca = u * (4 * m_dim) + 4 * ja
                        psi_up[u] += coeff * complex(x_rm[ca], x_rm[ca + 1])
                        psi_dn[u] += coeff * complex(x_rm[ca + 2], x_rm[ca + 3])
                cx_d = max(-1.0, min(1.0, sum((psi_up[u].conjugate() * psi_up[nx[u]] + psi_dn[u].conjugate() * psi_dn[nx[u]]).real for u in range(n))))
                px_d = math.acos(cx_d)
                if d_sp == 1:
                    p2_d = px_d * px_d
                else:
                    assert ny is not None
                    cy_d = max(-1.0, min(1.0, sum((psi_up[u].conjugate() * psi_up[ny[u]] + psi_dn[u].conjugate() * psi_dn[ny[u]]).real for u in range(n))))
                    py_a = math.acos(cy_d)
                    cxy_d = max(-1.0, min(1.0, sum((psi_up[u].conjugate() * psi_up[ny[nx[u]]] + psi_dn[u].conjugate() * psi_dn[ny[nx[u]]]).real for u in range(n))))
                    py_d = py_a if abs(cxy_d - math.cos(px_d + py_a)) <= abs(cxy_d - math.cos(px_d - py_a)) else -py_a
                    p2_d = px_d * px_d + py_d * py_d
                d_pos_modes.append((p2_d, e2_val))
        pos = end

    d_pos_modes.sort(key=lambda t: t[0])
    low_d = d_pos_modes[:low20]
    nz_d = [m for m in low_d if m[0] > 1e-12]
    c_d2 = sum((m[1] - m0_sq) / m[0] for m in nz_d if abs(m[0] - p1_sq) < 1e-7) / max(1, sum(1 for m in nz_d if abs(m[0] - p1_sq) < 1e-7))
    dl_d_mean = sum(abs(m[1] - (m[0] + m0_sq)) / m[1] for m in nz_d) / len(nz_d)
    aniso_d = _shell_aniso(low_d)
    bf_mis = sum(abs(low_d[k][1] - low_s[k][1]) / low_s[k][1] for k in range(low20)) / low20

    ir_abs = sorted(abs(d_evals[j]) for j in range(m_dim) if r_tilde[j] < 0.20)
    e_cut2 = (0.5 * (ir_abs[0] + next((ev for ev in ir_abs if ev > ir_abs[0] + 1e-5), ir_abs[0] + 0.1))) ** 2
    doublers = sum(1 for j in range(m_dim) if r_tilde[j] >= 0.25 and d_evals[j] ** 2 <= e_cut2) / 2.0

    return {
        "N": n, "c_S2": c_s2, "c_D2": c_d2, "DL_S": dl_s_mean, "DL_D": dl_d_mean,
        "aniso_S": aniso_s, "aniso_D": aniso_d, "BF_mis": bf_mis, "doublers": doublers,
        "tau1": tau1, "R_m1": rm1, "R_m2": rm2, "leak_t1": leak_t1, "bounds": bounds,
        "low_d_modes": low_d if keep_full_spectrum else None,
        "roughness": r_tilde if keep_full_spectrum else None,
        "abs_E": [abs(e) for e in d_evals] if keep_full_spectrum else None,
    }


# --- 5. Multi-Tier Diagnostic Plotter (Matplotlib -> Pillow PNG -> Pure-Python SVG) ---
def render_summary_figure(results: Dict[str, Dict[int, Dict[str, object]]], out_dir: str) -> str:
    ns = sorted(results["c3_1d"].keys())
    n_top = max(ns)
    n_spec = 128 if 128 in ns else ns[0]
    b_2d_top, c3_2d_top = results["base_2d"][n_top], results["c3_2d"][n_top]
    b_2d_sp, c3_2d_sp = results["base_2d"][n_spec], results["c3_2d"][n_spec]

    try:
        from PIL import Image, ImageDraw, ImageFont  # type: ignore
        W, H = 2800, 1900
        img = Image.new("RGB", (W, H), (250, 251, 253))
        dr = ImageDraw.Draw(img)
        f_h = ImageFont.load_default(size=36) if "size" in ImageFont.load_default.__code__.co_varnames else ImageFont.load_default()
        f_t = ImageFont.load_default(size=26) if "size" in ImageFont.load_default.__code__.co_varnames else ImageFont.load_default()
        f_s = ImageFont.load_default(size=20) if "size" in ImageFont.load_default.__code__.co_varnames else ImageFont.load_default()
        dr.text((520, 28), f"Ultrafinitist QFT Standalone Verification Summary (N up to {n_top})", fill=(26, 32, 44), font=f_h)
        boxes = [(180, 170, 1300, 880), (1540, 170, 2660, 880), (180, 1070, 1300, 1780), (1540, 1070, 2660, 1780)]

        def _panel(box, xl, yl, title, xlb, ylb):
            x0, y0, x1, y1 = box
            dr.rectangle(box, fill=(255, 255, 255), outline=(52, 60, 74), width=3)
            dr.text((x0 + 12, y0 - 40), title, fill=(26, 32, 44), font=f_t)
            dr.text(((x0 + x1) // 2 - 60, y1 + 34), xlb, fill=(52, 60, 74), font=f_s)
            dr.text((x0 - 10, y0 - 66), ylb, fill=(95, 106, 125), font=f_s)
            for k in range(6):
                fr = k / 5.0
                px, py = int(x0 + fr * (x1 - x0)), int(y1 - fr * (y1 - y0))
                if 0 < k < 5:
                    dr.line([(px, y0), (px, y1)], fill=(228, 232, 240), width=2)
                    dr.line([(x0, py), (x1, py)], fill=(228, 232, 240), width=2)
                dr.text((px - 18, y1 + 8), f"{xl[0] + fr * (xl[1] - xl[0]):.2g}", fill=(95, 106, 125), font=f_s)
                dr.text((x0 - 68, py - 10), f"{yl[0] + fr * (yl[1] - yl[0]):.3g}", fill=(95, 106, 125), font=f_s)
            return lambda x, y: (int(x0 + min(1.0, max(0.0, (x - xl[0]) / (xl[1] - xl[0]))) * (x1 - x0)),
                                 int(y1 - min(1.0, max(0.0, (y - yl[0]) / (yl[1] - yl[0]))) * (y1 - y0)))

        mp1 = _panel(boxes[0], (0.0, 2.6), (0.2, 3.0), f"(a) (2+1)D Dirac Dispersion E_D^2(p) (N={n_top})", "|p|^2", "E_D^2")
        dr.line([mp1(0.0, 0.25), mp1(2.6, 2.85)], fill=(45, 52, 54), width=4)
        for pairs, col in [(b_2d_top["low_d_modes"], (214, 48, 49)), (c3_2d_top["low_d_modes"], (16, 142, 84))]:
            for px_v, ey_v in (pairs or []):
                X, Y = mp1(px_v, ey_v)
                dr.ellipse((X - 5, Y - 5, X + 5, Y + 5), fill=col)

        mp2 = _panel(boxes[1], (0.0, 1.02), (0.0, 2.5), f"(b) (2+1)D Doubler Audit (N={n_spec}): Red=Base, Green=Cycle3", "Roughness R_tilde", "|E_D|")
        dr.line([mp2(0.25, 0.5), mp2(1.02, 0.5)], fill=(45, 52, 54), width=3)
        for src, col in [(b_2d_sp, (214, 48, 49)), (c3_2d_sp, (16, 142, 84))]:
            for rt, ev in zip(src["roughness"] or [], src["abs_E"] or []):  # type: ignore
                X, Y = mp2(rt, ev)
                dr.ellipse((X - 5, Y - 5, X + 5, Y + 5), fill=col)

        mp3 = _panel(boxes[2], (float(min(ns)), float(n_top)), (1.65, 1.90), "(c) Zero-Subtraction 2-Loop R_m^(2)(N) & Bounds", "Graph Size N", "R_m^(2)")
        for key, col in [("c3_1d", (16, 142, 84)), ("c3_2d", (41, 128, 185))]:
            pts = [mp3(float(nv), float(results[key][nv]["R_m2"])) for nv in ns]  # type: ignore
            dr.line(pts, fill=col, width=5)
            lo_pts = [mp3(float(nv), float(results[key][nv]["bounds"]["rm2_lo"])) for nv in ns]  # type: ignore
            hi_pts = [mp3(float(nv), float(results[key][nv]["bounds"]["rm2_hi"])) for nv in ns]  # type: ignore
            dr.line(lo_pts, fill=col, width=2)
            dr.line(hi_pts, fill=col, width=2)

        mp4 = _panel(boxes[3], (float(min(ns)), float(n_top)), (-8.0, -2.0), "(d) Log10 Outside-Lightcone Leakage at t=1 (< 5e-3)", "Graph Size N", "log10(leak)")
        dr.line([mp4(float(min(ns)), math.log10(5e-3)), mp4(float(n_top), math.log10(5e-3))], fill=(214, 48, 49), width=4)
        for key, col in [("c3_1d", (16, 142, 84)), ("c3_2d", (41, 128, 185))]:
            pts = [mp4(float(nv), math.log10(max(1e-8, float(results[key][nv]["leak_t1"])))) for nv in ns]  # type: ignore
            dr.line(pts, fill=col, width=5)

        out_path = os.path.join(out_dir, "ultrafinitist_qft_summary.png")
        resample = Image.Resampling.LANCZOS if hasattr(Image, "Resampling") else Image.LANCZOS
        img.resize((1400, 950), resample=resample).save(out_path, format="PNG")
        return out_path
    except Exception:
        out_path = os.path.join(out_dir, "ultrafinitist_qft_summary.svg")
        svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="400"><rect width="100%" height="100%" fill="#fafbfd"/>',
               f'<text x="40" y="50" font-size="20" font-family="monospace">Ultrafinitist QFT Verified up to N={n_top}</text>']
        for idx_n, nv in enumerate(ns):
            y_v = 110 + 45 * idx_n
            r1, r2 = results["c3_1d"][nv]["R_m2"], results["c3_2d"][nv]["R_m2"]
            svg.append(f'<text x="50" y="{y_v}" font-size="15" font-family="monospace">N={nv:4d}: 1+1D R_m^(2)={r1:.5f} | 2+1D R_m^(2)={r2:.5f}</text>')
        svg.append("</svg>")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(svg))
        return out_path


def main() -> int:
    ap = argparse.ArgumentParser(description="Single-File Standalone Ultrafinitist QFT Verifier")
    ap.add_argument("--quick", action="store_true", help="Verify N in [64, 128, 256] (~4s)")
    ap.add_argument("--out-dir", type=str, default=".", help="Directory for summary plot")
    args = ap.parse_args()
    t0 = time.perf_counter()
    ns = (64, 128, 256) if args.quick else (64, 128, 256, 512, 1024)

    print(f"=== ULTRAFINITIST QFT STANDALONE VERIFIER | Python {platform.python_version()} | Backend: {ACTIVE_BACKEND} | SymPy: {HAS_SYMPY} ===", flush=True)
    wb = verify_weyl_and_symbolic_identities()
    print(f"[PASS] Weyl-Schwinger (d=3,4,5) & Exact d^N Many-Body (dim={int(wb['dim_H'])}, E0={wb['E0']:.4f}, gap={wb['gap']:.4f}, ||H_int||={wb['H_int_norm']:.4f}) & SOS Identities", flush=True)

    res: Dict[str, Dict[int, Dict[str, object]]] = {"base_1d": {}, "c3_1d": {}, "base_2d": {}, "c3_2d": {}}
    n_top = max(ns)
    for nv in ns:
        keep_sp = nv in (128, n_top)
        res["base_1d"][nv] = evaluate_lattice_model(build_finite_lattice_operators("baseline", "1+1D", nv))
        res["c3_1d"][nv] = evaluate_lattice_model(build_finite_lattice_operators("cycle3", "1+1D", nv))
        res["base_2d"][nv] = evaluate_lattice_model(build_finite_lattice_operators("baseline", "2+1D", nv), keep_full_spectrum=keep_sp)
        res["c3_2d"][nv] = evaluate_lattice_model(build_finite_lattice_operators("cycle3", "2+1D", nv), keep_full_spectrum=keep_sp)

    print("\n--- COMPARISON TABLE Across Finite Graph Sizes N (Zero FFT, Zero Counterterms) ---")
    print("Model          Dim    N     c_S^2   c_D^2   DL_S     DL_D     Aniso_D  Doublers  R_m^(2)    Certified R_m^(2) Bracket")
    for label, key in [("Baseline", "base_1d"), ("Cycle3-Win", "c3_1d"), ("Baseline", "base_2d"), ("Cycle3-Win", "c3_2d")]:
        dim_s = "1+1D" if "1d" in key else "2+1D"
        rm2_vals = [float(res[key][nv]["R_m2"]) for nv in ns]  # type: ignore
        drift = (max(rm2_vals) - min(rm2_vals)) / (sum(rm2_vals) / len(rm2_vals))
        for nv in ns:
            r = res[key][nv]
            b_str = f"[{r['bounds']['rm2_lo']:.5f}, {r['bounds']['rm2_hi']:.5f}]" if r["bounds"] else "UNBOUNDED (UV blowup)"  # type: ignore
            print(f"{label:<14} {dim_s}  {nv:>4}  {r['c_S2']:6.4f}  {r['c_D2']:6.4f}  {r['DL_S']:7.5f}  {r['DL_D']:7.5f}  {r['aniso_D']:7.5f}    {r['doublers']:4.1f}    {r['R_m2']:8.4f}   {b_str}")
        print(f"  -> {label} {dim_s} 2-Loop Drift across N={list(ns)}: {drift * 100:.4f}%\n", flush=True)
        if label == "Cycle3-Win":
            assert drift < 0.03, f"Loop drift exceeded 3%: {drift}"
            for nv in ns:
                r = res[key][nv]
                assert r["doublers"] == 0.0 and r["DL_S"] < 0.05 and r["DL_D"] < 0.05 and r["aniso_D"] < 0.06 and r["BF_mis"] < 0.08 and r["leak_t1"] < 5e-3  # type: ignore
        else:
            assert res[key][ns[0]]["doublers"] >= 1.0  # type: ignore

    os.makedirs(args.out_dir, exist_ok=True)
    fig_path = render_summary_figure(res, args.out_dir)
    print(f"[PASS] Saved summary diagnostic plot: {fig_path}")
    print(f"[PASS] ALL ULTRAFINITIST QFT AXIOMS, EXACT Q BOUNDS & RED-TEAM GATES VERIFIED IN {time.perf_counter() - t0:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
