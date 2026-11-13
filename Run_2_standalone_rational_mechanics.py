#!/usr/bin/env python3
"""Standalone Zero-Dependency Ultrafinitist Rational Mechanics Engine & Benchmark Suite.

Self-contained Python 3.9+ stdlib implementation of exact rational (`fractions.Fraction`)
discrete Hamiltonian & manifold-projected integrators across all three benchmark battles:
  - Battle 1: 2D Coupled Chaotic Quartic Oscillator (Pullen-Edmonds / Yang-Mills)
  - Battle 2: 3D Interacting 3-Body SO(3) Molecular Quadrance Cluster (7 exact invariants)
  - Battle 3: 2D Eccentric Newtonian Kepler 1/r Problem (Levi-Civita Conic Feistel)

Architecture:
  - Part I:   Pure-Fraction Ultrafinitist Q-Engine (strictly `+`, `-`, `*`, `/`; zero `float`,
              zero `math`, zero `isinstance`/`type`, zero Heron loops, zero sub-stepping).
  - Part II:  Zero-Handicap IEEE-754 `float64` Continuum Baselines (`Euler`, `Velocity Verlet`,
              `RK4`, and `RK4(dt_ref=1/1024)`).
  - Part III: End-to-End 1,000-Step Verification & Strict `0/1` Conservation / Reversibility
              Assertions across `Delta_t in {1/2, 1/4, 1/8}` + Planar 3-Body Submanifolds.
  - Part IV:  Zero-Dependency Vector SVG + Pure-Stdlib `zlib` Truecolor PNG Plot Generator.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
import zlib
from fractions import Fraction
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple


# ============================================================================
# PART I: PURE-FRACTION ULTRAFINITIST Q-ENGINE (STRICTLY ISOLATED OVER Q)
# Zero float, zero math, zero isinstance/type, zero loops/Heron/sub-stepping.
# ============================================================================

# Canonical macroscopic time steps required across all cycles
CANONICAL_DTS: Tuple[Fraction, Fraction, Fraction] = (
    Fraction(1, 2),
    Fraction(1, 4),
    Fraction(1, 8),
)


# ============================================================================
# Battle 1: 2D Coupled Chaotic Quartic Oscillator (Pullen-Edmonds)
# H(q, p) = 1/2 (p1^2 + p2^2) + 1/2 (q1^2 + q2^2) + 1/2 * alpha * q1^2 * q2^2
# ============================================================================

BATTLE1_ALPHA: Fraction = Fraction(4, 1)
BATTLE1_MASSES: Tuple[Fraction, Fraction] = (Fraction(1, 1), Fraction(1, 1))
BATTLE1_Q0: Tuple[Fraction, Fraction] = (Fraction(1, 1), Fraction(3, 5))
BATTLE1_P0: Tuple[Fraction, Fraction] = (Fraction(2, 5), Fraction(6, 5))


def battle1_potential_q(q: Sequence[Any]) -> Any:
    """Exact rational potential V(q1, q2) = 1/2 (q1^2 + q2^2) + 1/2 * alpha * q1^2 * q2^2."""
    q1, q2 = q[0], q[1]
    q1_sq = q1 * q1
    q2_sq = q2 * q2
    half = Fraction(1, 2)
    return half * (q1_sq + q2_sq) + (half * BATTLE1_ALPHA) * (q1_sq * q2_sq)


def battle1_hamiltonian_q(q: Sequence[Any], p: Sequence[Any]) -> Any:
    """Exact rational Hamiltonian H_1(q, p) = T(p) + V(q) over Q."""
    p1, p2 = p[0], p[1]
    ke = Fraction(1, 2) * (p1 * p1 + p2 * p2)
    return ke + battle1_potential_q(q)


def battle1_force_q(q: Sequence[Any]) -> Tuple[Any, Any]:
    """Exact rational analytical force F(q) = -grad V(q) over Q."""
    q1, q2 = q[0], q[1]
    one = Fraction(1, 1)
    f1 = -q1 * (one + BATTLE1_ALPHA * (q2 * q2))
    f2 = -q2 * (one + BATTLE1_ALPHA * (q1 * q1))
    return (f1, f2)


def battle1_rhs_q(q: Sequence[Any], p: Sequence[Any]) -> Tuple[Tuple[Any, Any], Tuple[Any, Any]]:
    """Exact rational Hamiltonian vector field X_H(q, p) = (dq/dt, dp/dt)."""
    dq = (p[0], p[1])
    dp = battle1_force_q(q)
    return (dq, dp)


def battle1_invariants_q(q: Sequence[Any], p: Sequence[Any]) -> Dict[str, Any]:
    """Conserved invariants dictionary for Battle 1 over Q."""
    return {"H": battle1_hamiltonian_q(q, p)}



# ============================================================================
# Battle 2: 3D Interacting 3-Body SO(3) Molecular Quadrance Cluster
# Exact S_3-symmetric pairwise-quadrance Wildberger-Archimedes potential:
#   V_2(q) = 1/2 * k_s * (Q_12 + Q_23 + Q_31)
#            + 1/4 * g_c * (2*(Q_12*Q_23 + Q_23*Q_31 + Q_31*Q_12)
#                           - (Q_12^2 + Q_23^2 + Q_31^2))
# where Q_ij = ||q_i - q_j||^2 and Archimedes' triangle-area identity gives
#   2*(Q_12*Q_23 + Q_23*Q_31 + Q_31*Q_12) - (Q_12^2 + Q_23^2 + Q_31^2) == 4 * ||u x v||^2.
# Exact 7 conserved invariants over Q: H (1), P (3), L = sum_i q_i x p_i (3).
# ============================================================================

BATTLE2_K_SPRING: Fraction = Fraction(1, 1)
BATTLE2_G_COUPLE: Fraction = Fraction(1, 1)
BATTLE2_MASSES: Tuple[Fraction, ...] = (
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
    Fraction(1, 1),
)

# Generic 3D non-coplanar initial cluster with non-zero P_0 and 3D L_0
BATTLE2_Q0: Tuple[Fraction, ...] = (
    Fraction(2, 3),
    Fraction(-1, 3),
    Fraction(3, 10),
    Fraction(-1, 3),
    Fraction(-8, 15),
    Fraction(-1, 5),
    Fraction(-1, 30),
    Fraction(17, 30),
    Fraction(1, 20),
)

BATTLE2_P0: Tuple[Fraction, ...] = (
    Fraction(69, 200),
    Fraction(307, 300),
    Fraction(67, 120),
    Fraction(179, 200),
    Fraction(-269, 300),
    Fraction(-43, 600),
    Fraction(-47, 50),
    Fraction(-49, 150),
    Fraction(-29, 75),
)


def battle2_quadrances_q(q: Sequence[Any]) -> Tuple[Any, Any, Any]:
    """Compute exact pairwise quadrances (Q_12, Q_23, Q_31) where Q_ij = ||q_i - q_j||^2."""
    d12x = q[0] - q[3]
    d12y = q[1] - q[4]
    d12z = q[2] - q[5]
    q12 = d12x * d12x + d12y * d12y + d12z * d12z

    d23x = q[3] - q[6]
    d23y = q[4] - q[7]
    d23z = q[5] - q[8]
    q23 = d23x * d23x + d23y * d23y + d23z * d23z

    d31x = q[6] - q[0]
    d31y = q[7] - q[1]
    d31z = q[8] - q[2]
    q31 = d31x * d31x + d31y * d31y + d31z * d31z

    return (q12, q23, q31)


def battle2_potential_q(q: Sequence[Any]) -> Any:
    """Exact S_3-symmetric rational pure-quadrance molecular potential V_2(Q_12, Q_23, Q_31)."""
    q12, q23, q31 = battle2_quadrances_q(q)
    half = Fraction(1, 2)
    quarter = Fraction(1, 4)
    two = Fraction(2, 1)
    v_harm = (half * BATTLE2_K_SPRING) * (q12 + q23 + q31)
    arch = two * (q12 * q23 + q23 * q31 + q31 * q12) - (q12 * q12 + q23 * q23 + q31 * q31)
    v_coup = (quarter * BATTLE2_G_COUPLE) * arch
    return v_harm + v_coup


def battle2_hamiltonian_q(q: Sequence[Any], p: Sequence[Any]) -> Any:
    """Exact rational 3-body Hamiltonian H_2(q, p) = 1/2 sum ||p_i||^2 + V_2(q)."""
    ke_sq = (
        p[0] * p[0]
        + p[1] * p[1]
        + p[2] * p[2]
        + p[3] * p[3]
        + p[4] * p[4]
        + p[5] * p[5]
        + p[6] * p[6]
        + p[7] * p[7]
        + p[8] * p[8]
    )
    return Fraction(1, 2) * ke_sq + battle2_potential_q(q)


def battle2_force_q(q: Sequence[Any]) -> Tuple[Any, ...]:
    """Exact S_3-symmetric rational analytical 3-body forces F_i(q) = -grad_{q_i} V_2(q)."""
    q12, q23, q31 = battle2_quadrances_q(q)
    c12 = BATTLE2_K_SPRING + BATTLE2_G_COUPLE * (q23 + q31 - q12)
    c23 = BATTLE2_K_SPRING + BATTLE2_G_COUPLE * (q31 + q12 - q23)
    c31 = BATTLE2_K_SPRING + BATTLE2_G_COUPLE * (q12 + q23 - q31)

    r12x, r12y, r12z = q[0] - q[3], q[1] - q[4], q[2] - q[5]
    r23x, r23y, r23z = q[3] - q[6], q[4] - q[7], q[5] - q[8]
    r31x, r31y, r31z = q[6] - q[0], q[7] - q[1], q[8] - q[2]

    f1x = -c12 * r12x + c31 * r31x
    f1y = -c12 * r12y + c31 * r31y
    f1z = -c12 * r12z + c31 * r31z

    f2x = c12 * r12x - c23 * r23x
    f2y = c12 * r12y - c23 * r23y
    f2z = c12 * r12z - c23 * r23z

    f3x = -c31 * r31x + c23 * r23x
    f3y = -c31 * r31y + c23 * r23y
    f3z = -c31 * r31z + c23 * r23z

    return (f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z)


def battle2_rhs_q(q: Sequence[Any], p: Sequence[Any]) -> Tuple[Tuple[Any, ...], Tuple[Any, ...]]:
    """Exact rational Hamiltonian vector field X_H(q, p) for Battle 2."""
    dq = (p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8])
    return (dq, battle2_force_q(q))


def battle2_linear_momentum_q(q: Sequence[Any], p: Sequence[Any]) -> Tuple[Any, Any, Any]:
    """Exact total linear momentum P = p_1 + p_2 + p_3 over Q^3."""
    return (
        p[0] + p[3] + p[6],
        p[1] + p[4] + p[7],
        p[2] + p[5] + p[8],
    )


def battle2_angular_momentum_q(q: Sequence[Any], p: Sequence[Any]) -> Tuple[Any, Any, Any]:
    """Exact total angular momentum L = sum_{i=1}^3 q_i x p_i over Q^3."""
    lx = (
        (q[1] * p[2] - q[2] * p[1])
        + (q[4] * p[5] - q[5] * p[4])
        + (q[7] * p[8] - q[8] * p[7])
    )
    ly = (
        (q[2] * p[0] - q[0] * p[2])
        + (q[5] * p[3] - q[3] * p[5])
        + (q[8] * p[6] - q[6] * p[8])
    )
    lz = (
        (q[0] * p[1] - q[1] * p[0])
        + (q[3] * p[4] - q[4] * p[3])
        + (q[6] * p[7] - q[7] * p[6])
    )
    return (lx, ly, lz)


def battle2_invariants_q(q: Sequence[Any], p: Sequence[Any]) -> Dict[str, Any]:
    """All 7 exact rational conserved invariants {H, P (3), L (3)} for Battle 2."""
    return {
        "H": battle2_hamiltonian_q(q, p),
        "P": battle2_linear_momentum_q(q, p),
        "L": battle2_angular_momentum_q(q, p),
    }



# ============================================================================
# Battle 3: 2D Newtonian Kepler 1/r Problem (Levi-Civita Spinor Regularization)
# Physical Hamiltonian: H(x, y, px, py) = 1/2 (px^2 + py^2) - mu / sqrt(x^2 + y^2)
# Rational Levi-Civita lift:
#   x = u^2 - v^2,  y = 2*u*v   =>  r = u^2 + v^2 in Q  (x^2 + y^2 == r^2 identically!)
#   px = (u*pu - v*pv) / (2*r),     py = (v*pu + u*pv) / (2*r)
# Canonical eccentric planetary orbit:
#   mu = 1, m = 1, semi-major axis a = 1, energy E_0 = -1/2, eccentricity e = 336/625 = 0.5376
#   Aphelion initial state: (u0, v0) = (31/25, 0), (pu0, pv0) = (0, 34/25)
#   Cartesian initial state: (x0, y0) = (961/625, 0), (px0, py0) = (0, 17/31)
# ============================================================================

BATTLE3_MU: Fraction = Fraction(1, 1)
BATTLE3_MASS: Fraction = Fraction(1, 1)
BATTLE3_MASSES: Tuple[Fraction, Fraction] = (Fraction(1, 1), Fraction(1, 1))

BATTLE3_SPINOR_U0: Tuple[Fraction, Fraction] = (Fraction(31, 25), Fraction(0, 1))
BATTLE3_SPINOR_PU0: Tuple[Fraction, Fraction] = (Fraction(0, 1), Fraction(34, 25))

BATTLE3_Q0: Tuple[Fraction, Fraction] = (Fraction(961, 625), Fraction(0, 1))
BATTLE3_P0: Tuple[Fraction, Fraction] = (Fraction(0, 1), Fraction(17, 31))


def battle3_spinor_to_cartesian_q(
    u: Any, v: Any, pu: Any, pv: Any
) -> Tuple[Tuple[Any, Any], Tuple[Any, Any], Any]:
    """Map Levi-Civita spinor state (u, v, pu, pv) -> ((x, y), (px, py), r) over Q without sqrt."""
    u_sq = u * u
    v_sq = v * v
    r = u_sq + v_sq
    two = Fraction(2, 1)
    x = u_sq - v_sq
    y = two * u * v
    denom = two * r
    px = (u * pu - v * pv) / denom
    py = (v * pu + u * pv) / denom
    return ((x, y), (px, py), r)


def battle3_hamiltonian_spinor_q(u: Any, v: Any, pu: Any, pv: Any) -> Any:
    """Exact rational Kepler energy H = (pu^2 + pv^2) / (8*(u^2 + v^2)) - mu / (u^2 + v^2)."""
    r = u * u + v * v
    return (pu * pu + pv * pv) / (Fraction(8, 1) * r) - BATTLE3_MU / r


def battle3_lz_spinor_q(u: Any, v: Any, pu: Any, pv: Any) -> Any:
    """Exact rational planar angular momentum Lz = 1/2 (u*pv - v*pu) == x*py - y*px."""
    return Fraction(1, 2) * (u * pv - v * pu)


def battle3_lrl_spinor_q(u: Any, v: Any, pu: Any, pv: Any) -> Tuple[Any, Any]:
    """Exact rational Laplace-Runge-Lenz vector A = (py*Lz - mu*x/r, -px*Lz - mu*y/r)."""
    q_cart, p_cart, r = battle3_spinor_to_cartesian_q(u, v, pu, pv)
    x, y = q_cart
    px, py = p_cart
    lz = battle3_lz_spinor_q(u, v, pu, pv)
    ax = py * lz - BATTLE3_MU * x / r
    ay = -px * lz - BATTLE3_MU * y / r
    return (ax, ay)


def battle3_invariants_spinor_q(u: Any, v: Any, pu: Any, pv: Any) -> Dict[str, Any]:
    """Exact rational conserved invariants {H, Lz, LRL} for Battle 3."""
    return {
        "H": battle3_hamiltonian_spinor_q(u, v, pu, pv),
        "Lz": battle3_lz_spinor_q(u, v, pu, pv),
        "LRL": battle3_lrl_spinor_q(u, v, pu, pv),
    }



Vec3 = Tuple[Any, Any, Any]
State4 = Tuple[Vec3, Vec3, Vec3, Vec3]

_ZERO = Fraction(0, 1)
_ONE = Fraction(1, 1)
_NEG_ONE = Fraction(-1, 1)
_HALF = Fraction(1, 2)
_NEG_HALF = Fraction(-1, 2)
_QUARTER = Fraction(1, 4)
_NEG_QUARTER = Fraction(-1, 4)
_TWO = Fraction(2, 1)
_NEG_TWO = Fraction(-2, 1)
_THREE_QUARTERS = Fraction(3, 4)
_THREE_HALVES = Fraction(3, 2)


# ============================================================================
# 1. Scale-Covariant Floating-Dyadic Quantizers & Closed-Form [2/2] Pade Kernels
# ============================================================================


def _q6_bits(mag: Fraction, w: Fraction) -> Tuple[Fraction, Fraction]:
    """Extract 6 binary significand bits starting at dyadic weight `w`."""
    acc = _ZERO
    if mag >= w:
        acc = acc + w
        mag = mag - w
    w = w * _HALF
    if mag >= w:
        acc = acc + w
        mag = mag - w
    w = w * _HALF
    if mag >= w:
        acc = acc + w
        mag = mag - w
    w = w * _HALF
    if mag >= w:
        acc = acc + w
        mag = mag - w
    w = w * _HALF
    if mag >= w:
        acc = acc + w
        mag = mag - w
    w = w * _HALF
    if mag >= w:
        acc = acc + w
        mag = mag - w
    return (acc, mag)


def _normalize_dyadic(mag: Fraction) -> Tuple[Fraction, Fraction]:
    """Scale `mag > 0` by an exact power of 2 into `[1, 2)` across `(2^-64, 2^64)`."""
    scale = _ONE
    if mag >= _ONE:
        if mag >= Fraction(4294967296, 1):
            mag = mag * Fraction(1, 4294967296)
            scale = scale * Fraction(4294967296, 1)
        if mag >= Fraction(65536, 1):
            mag = mag * Fraction(1, 65536)
            scale = scale * Fraction(65536, 1)
        if mag >= Fraction(256, 1):
            mag = mag * Fraction(1, 256)
            scale = scale * Fraction(256, 1)
        if mag >= Fraction(16, 1):
            mag = mag * Fraction(1, 16)
            scale = scale * Fraction(16, 1)
        if mag >= Fraction(4, 1):
            mag = mag * Fraction(1, 4)
            scale = scale * Fraction(4, 1)
        if mag >= _TWO:
            mag = mag * _HALF
            scale = scale * _TWO
    else:
        if mag < Fraction(1, 4294967296):
            mag = mag * Fraction(4294967296, 1)
            scale = scale * Fraction(1, 4294967296)
        if mag < Fraction(1, 65536):
            mag = mag * Fraction(65536, 1)
            scale = scale * Fraction(1, 65536)
        if mag < Fraction(1, 256):
            mag = mag * Fraction(256, 1)
            scale = scale * Fraction(1, 256)
        if mag < Fraction(1, 16):
            mag = mag * Fraction(16, 1)
            scale = scale * Fraction(1, 16)
        if mag < Fraction(1, 4):
            mag = mag * Fraction(4, 1)
            scale = scale * Fraction(1, 4)
        if mag < _HALF:
            mag = mag * _TWO
            scale = scale * _HALF
        if mag < _ONE:
            mag = mag * _TWO
            scale = scale * _HALF
    return (mag, scale)


def quantize_sig29(val: Fraction) -> Fraction:
    """Odd-symmetric 30-significant-bit floating-dyadic quantizer (`Q(-x) == -Q(x)`)."""
    if val == _ZERO:
        return _ZERO
    sgn = _ONE
    mag = val
    if val < _ZERO:
        sgn = _NEG_ONE
        mag = -val
    norm, scale = _normalize_dyadic(mag)
    a1, norm = _q6_bits(norm, _ONE)
    a2, norm = _q6_bits(norm, Fraction(1, 64))
    a3, norm = _q6_bits(norm, Fraction(1, 4096))
    a4, norm = _q6_bits(norm, Fraction(1, 262144))
    a5, norm = _q6_bits(norm, Fraction(1, 16777216))
    rem = Fraction(1, 536870912) if norm >= Fraction(1, 1073741824) else _ZERO
    return sgn * scale * (a1 + a2 + a3 + a4 + a5 + rem)


def quantize_sig17(val: Fraction) -> Fraction:
    """Odd-symmetric 18-significant-bit floating-dyadic quantizer."""
    if val == _ZERO:
        return _ZERO
    sgn = _ONE
    mag = val
    if val < _ZERO:
        sgn = _NEG_ONE
        mag = -val
    norm, scale = _normalize_dyadic(mag)
    a1, norm = _q6_bits(norm, _ONE)
    a2, norm = _q6_bits(norm, Fraction(1, 64))
    a3, norm = _q6_bits(norm, Fraction(1, 4096))
    rem = Fraction(1, 131072) if norm >= Fraction(1, 262144) else _ZERO
    return sgn * scale * (a1 + a2 + a3 + rem)


def quantize_sig14(val: Fraction) -> Fraction:
    """Odd-symmetric 15-significant-bit floating-dyadic quantizer."""
    if val == _ZERO:
        return _ZERO
    sgn = _ONE
    mag = val
    if val < _ZERO:
        sgn = _NEG_ONE
        mag = -val
    norm, scale = _normalize_dyadic(mag)
    a1, norm = _q6_bits(norm, _ONE)
    a2, norm = _q6_bits(norm, Fraction(1, 64))
    b1 = Fraction(1, 4096) if norm >= Fraction(1, 4096) else _ZERO
    if b1 != _ZERO:
        norm = norm - b1
    b2 = Fraction(1, 8192) if norm >= Fraction(1, 8192) else _ZERO
    if b2 != _ZERO:
        norm = norm - b2
    rem = Fraction(1, 16384) if norm >= Fraction(1, 32768) else _ZERO
    return sgn * scale * (a1 + a2 + b1 + b2 + rem)


def quantize_sig11(val: Fraction) -> Fraction:
    """Odd-symmetric 12-significant-bit floating-dyadic quantizer."""
    if val == _ZERO:
        return _ZERO
    sgn = _ONE
    mag = val
    if val < _ZERO:
        sgn = _NEG_ONE
        mag = -val
    norm, scale = _normalize_dyadic(mag)
    a1, norm = _q6_bits(norm, _ONE)
    a2, norm = _q6_bits(norm, Fraction(1, 64))
    rem = Fraction(1, 2048) if norm >= Fraction(1, 4096) else _ZERO
    return sgn * scale * (a1 + a2 + rem)


def quantize_sig8(val: Fraction) -> Fraction:
    """Odd-symmetric 9-significant-bit floating-dyadic quantizer."""
    if val == _ZERO:
        return _ZERO
    sgn = _ONE
    mag = val
    if val < _ZERO:
        sgn = _NEG_ONE
        mag = -val
    norm, scale = _normalize_dyadic(mag)
    a1, norm = _q6_bits(norm, _ONE)
    b1 = Fraction(1, 64) if norm >= Fraction(1, 64) else _ZERO
    if b1 != _ZERO:
        norm = norm - b1
    b2 = Fraction(1, 128) if norm >= Fraction(1, 128) else _ZERO
    if b2 != _ZERO:
        norm = norm - b2
    rem = Fraction(1, 256) if norm >= Fraction(1, 512) else _ZERO
    return sgn * scale * (a1 + b1 + b2 + rem)


def pade22_sqrt_from_one(r_val: Any) -> Any:
    """Closed-form [2/2] Pade rational approximant `(5*r^2 + 10*r + 1)/(r^2 + 10*r + 5)` near `r = 1`."""
    r2 = r_val * r_val
    five = Fraction(5, 1)
    ten = Fraction(10, 1)
    return (five * r2 + ten * r_val + _ONE) / (r2 + ten * r_val + five)


def pade_sqrt_positive(r_val: Fraction) -> Fraction:
    """Range-reduced closed-form [2/2] Pade square-root surrogate on `(0, infty)` (~2.2e-4 rel accuracy)."""
    if r_val <= _ZERO:
        return _ZERO
    s0 = _ONE
    r = r_val
    if r >= _ONE:
        if r >= Fraction(65536, 1):
            r = r * Fraction(1, 65536)
            s0 = s0 * Fraction(256, 1)
        if r >= Fraction(256, 1):
            r = r * Fraction(1, 256)
            s0 = s0 * Fraction(16, 1)
        if r >= Fraction(16, 1):
            r = r * Fraction(1, 16)
            s0 = s0 * Fraction(4, 1)
        if r >= Fraction(4, 1):
            r = r * Fraction(1, 4)
            s0 = s0 * _TWO
    else:
        if r < Fraction(1, 65536):
            r = r * Fraction(65536, 1)
            s0 = s0 * Fraction(1, 256)
        if r < Fraction(1, 256):
            r = r * Fraction(256, 1)
            s0 = s0 * Fraction(1, 16)
        if r < Fraction(1, 16):
            r = r * Fraction(16, 1)
            s0 = s0 * Fraction(1, 4)
        if r < Fraction(1, 4):
            r = r * Fraction(4, 1)
            s0 = s0 * _HALF
        if r < _HALF:
            r = r * Fraction(4, 1)
            s0 = s0 * _HALF
    return s0 * pade22_sqrt_from_one(r)


def pade_sqrt_fine(r_val: Fraction) -> Fraction:
    """Two-stage closed-form [2/2] Pade rational surrogate (~1.5e-19 rel accuracy, zero Heron `x + r/x`)."""
    if r_val <= _ZERO:
        return _ZERO
    s1 = pade_sqrt_positive(r_val)
    return s1 * pade22_sqrt_from_one(r_val / (s1 * s1))


# ============================================================================
# 2. Battle 1: Pullen-Edmonds Chaotic Oscillator
# Separable Zero-Curl Shadow Verlet + Two-Stage Ellipsoid Vieta Projection
# ============================================================================


def b1_shadow_force_q(q1: Any, q2: Any, dt: Any) -> Tuple[Any, Any]:
    """Exact conservative `-grad W_{1,dt}` for separable `W_{1,dt}` (`dF1/dq2 == dF2/dq1` identically)."""
    cdt2 = _HALF * dt * dt
    d1 = _ONE + cdt2 * q1 * q1
    d2 = _ONE + cdt2 * q2 * q2
    f1 = -q1 - (BATTLE1_ALPHA * q1 * q2 * q2) / (d1 * d1 * d2)
    f2 = -q2 - (BATTLE1_ALPHA * q2 * q1 * q1) / (d2 * d2 * d1)
    return (f1, f2)


def vieta_ellipsoid3_kernel(
    c1: Any, c2: Any, c3: Any, pb: Tuple[Any, Any, Any], x_hat: Tuple[Any, Any, Any]
) -> Tuple[Any, Any, Any]:
    """Pure-field 3D ellipsoid Vieta involution preserving `c1*x^2 + c2*y^2 + c3*z^2` identically."""
    s0 = x_hat[0] + pb[0]
    s1 = x_hat[1] + pb[1]
    s2 = x_hat[2] + pb[2]
    num = _TWO * (c1 * pb[0] * s0 + c2 * pb[1] * s1 + c3 * pb[2] * s2)
    den = c1 * s0 * s0 + c2 * s1 * s1 + c3 * s2 * s2
    if den == _ZERO:
        return pb
    lam = num / den
    return (-pb[0] + lam * s0, -pb[1] + lam * s1, -pb[2] + lam * s2)


def _vieta_ellipsoid3_signed(
    c1: Fraction,
    c2: Fraction,
    c3: Fraction,
    p0: Tuple[Fraction, Fraction, Fraction],
    x_hat: Tuple[Fraction, Fraction, Fraction],
) -> Tuple[Fraction, Fraction, Fraction]:
    ip = c1 * p0[0] * x_hat[0] + c2 * p0[1] * x_hat[1] + c3 * p0[2] * x_hat[2]
    sgn = _NEG_ONE if ip < _ZERO else _ONE
    pb = (sgn * p0[0], sgn * p0[1], sgn * p0[2])
    return vieta_ellipsoid3_kernel(c1, c2, c3, pb, x_hat)


def b1_pre_rotate_kernel(
    q1_0: Any, q2_0: Any, p1_0: Any, p2_0: Any, eta1: Any, eta2: Any
) -> Tuple[Any, Any, Any, Any]:
    """Pure-field sequential Cayley pre-rotation conserving `H_1(q, p)` identically for all `(eta1, eta2)`."""
    w1 = _ONE + BATTLE1_ALPHA * q2_0 * q2_0
    d1 = _ONE + w1 * eta1 * eta1
    q1_s = ((_ONE - w1 * eta1 * eta1) * q1_0 + _TWO * eta1 * p1_0) / d1
    p1_s = (-_TWO * w1 * eta1 * q1_0 + (_ONE - w1 * eta1 * eta1) * p1_0) / d1

    w2 = _ONE + BATTLE1_ALPHA * q1_s * q1_s
    d2 = _ONE + w2 * eta2 * eta2
    q2_s = ((_ONE - w2 * eta2 * eta2) * q2_0 + _TWO * eta2 * p2_0) / d2
    p2_s = (-_TWO * w2 * eta2 * q2_0 + (_ONE - w2 * eta2 * eta2) * p2_0) / d2
    return (q1_s, q2_s, p1_s, p2_s)


def _pre_rotate_b1(
    q1_0: Fraction, q2_0: Fraction, p1_0: Fraction, p2_0: Fraction
) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    w1 = _ONE + BATTLE1_ALPHA * q2_0 * q2_0
    e1_twice = p1_0 * p1_0 + w1 * q1_0 * q1_0
    r1_est = pade_sqrt_positive(e1_twice)
    if p1_0 < _ZERO:
        r1_est = -r1_est
    denom1 = p1_0 + r1_est
    eta1 = quantize_sig17(-q1_0 / denom1) if denom1 != _ZERO else _ZERO

    d1 = _ONE + w1 * eta1 * eta1
    q1_s = ((_ONE - w1 * eta1 * eta1) * q1_0 + _TWO * eta1 * p1_0) / d1
    w2 = _ONE + BATTLE1_ALPHA * q1_s * q1_s
    e2_twice = p2_0 * p2_0 + w2 * q2_0 * q2_0
    r2_est = pade_sqrt_positive(e2_twice)
    if p2_0 < _ZERO:
        r2_est = -r2_est
    denom2 = p2_0 + r2_est
    eta2 = quantize_sig17(-q2_0 / denom2) if denom2 != _ZERO else _ZERO
    return b1_pre_rotate_kernel(q1_0, q2_0, p1_0, p2_0, eta1, eta2)


def b1_direct_chart_kernel(
    z_ref: Tuple[Any, Any, Any, Any],
    q1_h: Any,
    q2_h: Any,
    p1_h: Any,
    p2_h: Any,
    rho1: Any,
    rho2: Any,
    sgn1: Any = _ONE,
    sgn2: Any = _ONE,
) -> Tuple[Any, Any, Any, Any]:
    """Pure-field composed 2-stage ellipsoid Vieta projection preserving `H_1 == H_1(z_ref)` identically."""
    q1_0, q2_0, p1_0, p2_0 = z_ref
    c1_a = _HALF * (_ONE + BATTLE1_ALPHA * q2_0 * q2_0)
    pb1 = (sgn1 * q1_0, sgn1 * p1_0, sgn1 * p2_0)
    q1_st, p1_aux, p2_aux = vieta_ellipsoid3_kernel(
        c1_a, _HALF, _HALF, pb1, (q1_h, rho1 * p1_h, rho1 * p2_h)
    )
    c2_b = _HALF * (_ONE + BATTLE1_ALPHA * q1_st * q1_st)
    pb2 = (sgn2 * q2_0, sgn2 * p1_aux, sgn2 * p2_aux)
    q2_st, p1_st, p2_st = vieta_ellipsoid3_kernel(
        c2_b, _HALF, _HALF, pb2, (q2_h, rho2 * p1_h, rho2 * p2_h)
    )
    return (q1_st, q2_st, p1_st, p2_st)


def project_battle1_manifold(
    z_ref: Tuple[Fraction, Fraction, Fraction, Fraction],
    z_hat: Tuple[Fraction, Fraction, Fraction, Fraction],
) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    """Project shadow state `z_hat` onto `M_{H0}` with `0/1` energy conservation and `O(1)` bit-length."""
    q1_0, q2_0, p1_0, p2_0 = z_ref
    if z_hat == z_ref:
        return z_ref
    q1_h = q1_0 + quantize_sig29(z_hat[0] - q1_0)
    q2_h = q2_0 + quantize_sig29(z_hat[1] - q2_0)
    p1_h = p1_0 + quantize_sig29(z_hat[2] - p1_0)
    p2_h = p2_0 + quantize_sig29(z_hat[3] - p2_0)

    h0 = battle1_hamiltonian_q((q1_0, q2_0), (p1_0, p2_0))
    c1_a = _HALF * (_ONE + BATTLE1_ALPHA * q2_0 * q2_0)
    r1_direct = h0 - _HALF * q2_0 * q2_0 - c1_a * q1_h * q1_h
    p_sq_h = p1_h * p1_h + p2_h * p2_h
    p_sq_0 = p1_0 * p1_0 + p2_0 * p2_0
    if p_sq_0 > _ZERO and p_sq_h >= Fraction(1, 8) * h0 and r1_direct >= Fraction(1, 8) * h0:
        rho1_raw = pade22_sqrt_from_one((_TWO * r1_direct) / p_sq_h)
        rho1 = _ONE + quantize_sig29(rho1_raw - _ONE)
        ip1 = c1_a * q1_0 * q1_h + _HALF * rho1 * (p1_0 * p1_h + p2_0 * p2_h)
        sgn1 = _NEG_ONE if ip1 < _ZERO else _ONE
        pb1 = (sgn1 * q1_0, sgn1 * p1_0, sgn1 * p2_0)
        q1_st, p1_aux, p2_aux = vieta_ellipsoid3_kernel(
            c1_a, _HALF, _HALF, pb1, (q1_h, rho1 * p1_h, rho1 * p2_h)
        )
        c2_b = _HALF * (_ONE + BATTLE1_ALPHA * q1_st * q1_st)
        r2_rem = h0 - _HALF * q1_st * q1_st - c2_b * q2_h * q2_h
        if r2_rem > _ZERO:
            rho2 = _ONE + quantize_sig29(pade_sqrt_positive((_TWO * r2_rem) / p_sq_h) - _ONE)
        else:
            rho2 = _ONE
        ip2 = c2_b * q2_0 * q2_h + _HALF * rho2 * (p1_aux * p1_h + p2_aux * p2_h)
        sgn2 = _NEG_ONE if ip2 < _ZERO else _ONE
        return b1_direct_chart_kernel(z_ref, q1_h, q2_h, p1_h, p2_h, rho1, rho2, sgn1, sgn2)

    q1_s, q2_s, p1_s, p2_s = _pre_rotate_b1(q1_0, q2_0, p1_0, p2_0)
    h0_s = battle1_hamiltonian_q((q1_s, q2_s), (p1_s, p2_s))
    p_sq_s = p1_s * p1_s + p2_s * p2_s
    if abs(q1_h) <= abs(q2_h):
        c1_s = _HALF * (_ONE + BATTLE1_ALPHA * q2_s * q2_s)
        k_need_twice = _TWO * (h0_s - _HALF * q2_s * q2_s - c1_s * q1_h * q1_h)
        rho = quantize_sig29(pade_sqrt_positive(k_need_twice / p_sq_s))
        q1_st, p1_aux, p2_aux = _vieta_ellipsoid3_signed(
            c1_s, _HALF, _HALF, (q1_s, p1_s, p2_s), (q1_h, rho * p1_s, rho * p2_s)
        )
        c2_s = _HALF * (_ONE + BATTLE1_ALPHA * q1_st * q1_st)
        r2_rem = h0_s - _HALF * q1_st * q1_st - c2_s * q2_h * q2_h
        rho2 = (
            quantize_sig29(pade_sqrt_positive((_TWO * r2_rem) / p_sq_h))
            if (r2_rem > _ZERO and p_sq_h > _ZERO)
            else _ONE
        )
        q2_st, p1_st, p2_st = _vieta_ellipsoid3_signed(
            c2_s, _HALF, _HALF, (q2_s, p1_aux, p2_aux), (q2_h, rho2 * p1_h, rho2 * p2_h)
        )
        return (q1_st, q2_st, p1_st, p2_st)

    c2_s = _HALF * (_ONE + BATTLE1_ALPHA * q1_s * q1_s)
    k_need_twice = _TWO * (h0_s - _HALF * q1_s * q1_s - c2_s * q2_h * q2_h)
    rho = quantize_sig29(pade_sqrt_positive(k_need_twice / p_sq_s))
    q2_st, p1_aux, p2_aux = _vieta_ellipsoid3_signed(
        c2_s, _HALF, _HALF, (q2_s, p1_s, p2_s), (q2_h, rho * p1_s, rho * p2_s)
    )
    c1_b = _HALF * (_ONE + BATTLE1_ALPHA * q2_st * q2_st)
    r1_rem = h0_s - _HALF * q2_st * q2_st - c1_b * q1_h * q1_h
    rho2 = (
        quantize_sig29(pade_sqrt_positive((_TWO * r1_rem) / p_sq_h))
        if (r1_rem > _ZERO and p_sq_h > _ZERO)
        else _ONE
    )
    q1_st, p1_st, p2_st = _vieta_ellipsoid3_signed(
        c1_b, _HALF, _HALF, (q1_s, p1_aux, p2_aux), (q1_h, rho2 * p1_h, rho2 * p2_h)
    )
    return (q1_st, q2_st, p1_st, p2_st)


def build_state_from_qp_battle1(q: Sequence[Any], p: Sequence[Any]) -> Any:
    """Construct minimal globally-injective Battle 1 state `((qh1, qh2), (ph1, ph2), z_ref)`."""
    return ((q[0], q[1]), (p[0], p[1]), (q[0], q[1], p[0], p[1]))


def initial_state_battle1(dt: Optional[Fraction] = None) -> Any:
    """Return canonical Battle 1 initial state."""
    return build_state_from_qp_battle1(BATTLE1_Q0, BATTLE1_P0)


def step_battle1(
    state: Any, dt: Any, q_fn: Callable[[Any], Any] = quantize_sig29
) -> Any:
    """Single-step symmetric conservative shadow Velocity Verlet with `dt` factored outside `q_fn`."""
    qh_pair, ph_pair, z_ref = state
    qh1, qh2 = qh_pair
    ph1, ph2 = ph_pair
    f_ref1, f_ref2 = b1_shadow_force_q(z_ref[0], z_ref[1], dt)
    fa1, fa2 = b1_shadow_force_q(qh1, qh2, dt)
    half_dt = dt * _HALF
    ph1_m = ph1 + half_dt * (f_ref1 + q_fn(fa1 - f_ref1))
    ph2_m = ph2 + half_dt * (f_ref2 + q_fn(fa2 - f_ref2))
    qh1_n = qh1 + dt * ph1_m
    qh2_n = qh2 + dt * ph2_m
    fb1, fb2 = b1_shadow_force_q(qh1_n, qh2_n, dt)
    ph1_n = ph1_m + half_dt * (f_ref1 + q_fn(fb1 - f_ref1))
    ph2_n = ph2_m + half_dt * (f_ref2 + q_fn(fb2 - f_ref2))
    return ((qh1_n, qh2_n), (ph1_n, ph2_n), z_ref)


def extract_qp_battle1(state: Any) -> Tuple[List[Any], List[Any]]:
    """Extract physical `([q1, q2], [p1, p2])` on `M_{H0}`."""
    qh_pair, ph_pair, z_ref = state
    q1_st, q2_st, p1_st, p2_st = project_battle1_manifold(
        z_ref, (qh_pair[0], qh_pair[1], ph_pair[0], ph_pair[1])
    )
    return ([q1_st, q2_st], [p1_st, p2_st])


def invariants_battle1(state: Any) -> Dict[str, Any]:
    """Evaluate exact physical Hamiltonian `H` on `extract_qp_battle1(state)`."""
    q, p = extract_qp_battle1(state)
    return {"H": battle1_hamiltonian_q(q, p)}


# ============================================================================
# 3. Battle 2: 3D 3-Body SO(3) Molecular Quadrance Cluster
# Zero-Curl Shadow Verlet + Planar & 3D 7-Invariant Vieta Fibration
# ============================================================================


def _dot3(a: Vec3, b: Vec3) -> Any:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross3(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _add3(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub3(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _scale3(s: Any, a: Vec3) -> Vec3:
    return (s * a[0], s * a[1], s * a[2])


def _qvec(v: Vec3, q_fn: Callable[[Any], Any]) -> Vec3:
    return (q_fn(v[0]), q_fn(v[1]), q_fn(v[2]))


def _q_state8(za: State4) -> State4:
    return (
        _qvec(za[0], quantize_sig8),
        _qvec(za[1], quantize_sig8),
        _qvec(za[2], quantize_sig8),
        _qvec(za[3], quantize_sig8),
    )


def _v2_arch_uv(u: Vec3, v: Vec3) -> Any:
    uxv = _cross3(u, v)
    return (
        _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(u, u)
        + BATTLE2_K_SPRING * _dot3(v, v)
        + BATTLE2_G_COUPLE * _dot3(uxv, uxv)
    )


def _h_int_arch(u: Vec3, v: Vec3, pu: Vec3, pv: Vec3) -> Any:
    return _dot3(pu, pu) + _THREE_QUARTERS * _dot3(pv, pv) + _v2_arch_uv(u, v)


def _l_int_q(u: Vec3, v: Vec3, pu: Vec3, pv: Vec3) -> Vec3:
    return _add3(_cross3(u, pu), _cross3(v, pv))


def _pick_ek(vec: Vec3) -> Tuple[Vec3, Fraction]:
    a0, a1, a2 = abs(vec[0]), abs(vec[1]), abs(vec[2])
    if a0 >= a1 and a0 >= a2:
        return ((_ONE, _ZERO, _ZERO), vec[0])
    if a1 >= a2:
        return ((_ZERO, _ONE, _ZERO), vec[1])
    return ((_ZERO, _ZERO, _ONE), vec[2])


def _cartesian_to_jacobi(q: Sequence[Any], p: Sequence[Any]) -> Tuple[Vec3, Vec3, Vec3, Vec3, Vec3, Vec3]:
    third = Fraction(1, 3)
    two_thirds = Fraction(2, 3)
    rc = (
        third * (q[0] + q[3] + q[6]),
        third * (q[1] + q[4] + q[7]),
        third * (q[2] + q[5] + q[8]),
    )
    pt = (p[0] + p[3] + p[6], p[1] + p[4] + p[7], p[2] + p[5] + p[8])
    u = (q[0] - q[3], q[1] - q[4], q[2] - q[5])
    v = (
        q[6] - _HALF * (q[0] + q[3]),
        q[7] - _HALF * (q[1] + q[4]),
        q[8] - _HALF * (q[2] + q[5]),
    )
    pu = (_HALF * (p[0] - p[3]), _HALF * (p[1] - p[4]), _HALF * (p[2] - p[5]))
    pv = (
        two_thirds * p[6] - third * (p[0] + p[3]),
        two_thirds * p[7] - third * (p[1] + p[4]),
        two_thirds * p[8] - third * (p[2] + p[5]),
    )
    return (rc, pt, u, v, pu, pv)


def _jacobi_to_cartesian(
    rc: Vec3, pt: Vec3, u: Vec3, v: Vec3, pu: Vec3, pv: Vec3
) -> Tuple[Tuple[Any, ...], Tuple[Any, ...]]:
    third = Fraction(1, 3)
    two_thirds = Fraction(2, 3)
    q1x = rc[0] + _HALF * u[0] - third * v[0]
    q1y = rc[1] + _HALF * u[1] - third * v[1]
    q1z = rc[2] + _HALF * u[2] - third * v[2]

    q2x = rc[0] - _HALF * u[0] - third * v[0]
    q2y = rc[1] - _HALF * u[1] - third * v[1]
    q2z = rc[2] - _HALF * u[2] - third * v[2]

    q3x = rc[0] + two_thirds * v[0]
    q3y = rc[1] + two_thirds * v[1]
    q3z = rc[2] + two_thirds * v[2]

    p1x = third * pt[0] + pu[0] - _HALF * pv[0]
    p1y = third * pt[1] + pu[1] - _HALF * pv[1]
    p1z = third * pt[2] + pu[2] - _HALF * pv[2]

    p2x = third * pt[0] - pu[0] - _HALF * pv[0]
    p2y = third * pt[1] - pu[1] - _HALF * pv[1]
    p2z = third * pt[2] - pu[2] - _HALF * pv[2]

    p3x = third * pt[0] + pv[0]
    p3y = third * pt[1] + pv[1]
    p3z = third * pt[2] + pv[2]

    return (
        (q1x, q1y, q1z, q2x, q2y, q2z, q3x, q3y, q3z),
        (p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z),
    )


def b2_shadow_force_uv(uu: Vec3, vv: Vec3, dt: Any) -> Tuple[Vec3, Vec3]:
    """Exact conservative `-grad W_{2,dt}` for `W_{2,dt}(S, A) = k*S + g*A / (1 + 1/4*dt^2*S)^2`."""
    uu_sq = _dot3(uu, uu)
    vv_sq = _dot3(vv, vv)
    uv = _dot3(uu, vv)
    s_val = _THREE_QUARTERS * uu_sq + vv_sq
    uxv = _cross3(uu, vv)
    a_val = _dot3(uxv, uxv)
    cdt2 = _QUARTER * dt * dt
    d = _ONE + cdt2 * s_val
    d2 = d * d
    d3 = d2 * d
    coef_s = BATTLE2_K_SPRING - (_TWO * BATTLE2_G_COUPLE * cdt2 * a_val) / d3
    coef_a = BATTLE2_G_COUPLE / d2
    dvu = _add3(
        _scale3(_THREE_HALVES * coef_s + _TWO * coef_a * vv_sq, uu),
        _scale3(-_TWO * coef_a * uv, vv),
    )
    dvv = _add3(
        _scale3(_TWO * coef_s + _TWO * coef_a * uu_sq, vv),
        _scale3(-_TWO * coef_a * uv, uu),
    )
    return (_scale3(_NEG_ONE, dvu), _scale3(_NEG_ONE, dvv))


def apply_dem_gauge(za: State4, eps_u: Any, eps_v: Any, eta_dem: Any) -> State4:
    """Exact rational `O(2)_dem` kinematic-democracy transformation preserving `H_int` and `L_int`."""
    ua, va, pua, pva = za
    u1 = _scale3(eps_u, ua)
    pu1 = _scale3(eps_u, pua)
    v1 = _scale3(eps_v, va)
    pv1 = _scale3(eps_v, pva)
    if eta_dem == _ZERO:
        return (u1, v1, pu1, pv1)
    t2 = _THREE_QUARTERS * eta_dem * eta_dem
    den = _ONE + t2
    b11 = (_ONE - t2) / den
    b12 = (_TWO * eta_dem) / den
    b21 = (Fraction(-3, 2) * eta_dem) / den
    return (
        _add3(_scale3(b11, u1), _scale3(b12, v1)),
        _add3(_scale3(b21, u1), _scale3(b11, v1)),
        _sub3(_scale3(b11, pu1), _scale3(b21, pv1)),
        _add3(_scale3(-b12, pu1), _scale3(b11, pv1)),
    )


def rot_vec_cayley(x: Vec3, a: Vec3, c_rot: Any) -> Vec3:
    """Pure-field `SO(3)` Cayley vector rotation around axis `a` with `c_rot = 2 / (1 + a.a)`."""
    ax = _cross3(a, x)
    return _add3(x, _scale3(c_rot, _add3(ax, _cross3(a, ax))))


def apply_l_rot_q(za: State4, l_vec: Vec3, sig_l: Any) -> State4:
    """Exact rational Cayley rotation around `L_int` preserving `H_int` and `L_int` identically."""
    if sig_l == _ZERO:
        return za
    a = _scale3(sig_l, l_vec)
    aa = _dot3(a, a)
    if aa == _ZERO:
        return za
    c_rot = _TWO / (_ONE + aa)
    return (
        rot_vec_cayley(za[0], a, c_rot),
        rot_vec_cayley(za[1], a, c_rot),
        rot_vec_cayley(za[2], a, c_rot),
        rot_vec_cayley(za[3], a, c_rot),
    )


def apply_l_half_turn_q(za: State4, l_vec: Vec3) -> State4:
    """Exact rational 180-degree spatial rotation around `L_int` preserving `H_int` and `L_int`."""
    ll = _dot3(l_vec, l_vec)
    if ll == _ZERO:
        return za
    two_inv_ll = _TWO / ll
    return (
        _sub3(_scale3(two_inv_ll * _dot3(l_vec, za[0]), l_vec), za[0]),
        _sub3(_scale3(two_inv_ll * _dot3(l_vec, za[1]), l_vec), za[1]),
        _sub3(_scale3(two_inv_ll * _dot3(l_vec, za[2]), l_vec), za[2]),
        _sub3(_scale3(two_inv_ll * _dot3(l_vec, za[3]), l_vec), za[3]),
    )


def transvection_vieta_kernel(za: State4, w: Tuple[Any, Any, Any, Any]) -> State4:
    """Pure-field `SO(3)`-equivariant symplectic Vieta chord preserving `H_int` and `L_int`."""
    u, v, pu, pv = za
    w1, w2, w3, w4 = w
    d = _add3(_add3(_scale3(w1, u), _scale3(w2, v)), _add3(_scale3(w3, pu), _scale3(w4, pv)))
    du = _scale3(w3, d)
    dv = _scale3(w4, d)
    dpu = _scale3(-w1, d)
    dpv = _scale3(-w2, d)
    uxv0 = _cross3(u, v)
    duxv = _add3(_cross3(du, v), _cross3(u, dv))
    c1 = _TWO * (
        _dot3(pu, dpu)
        + _THREE_QUARTERS * _dot3(pv, dpv)
        + _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(u, du)
        + BATTLE2_K_SPRING * _dot3(v, dv)
        + BATTLE2_G_COUPLE * _dot3(uxv0, duxv)
    )
    c2 = (
        _dot3(dpu, dpu)
        + _THREE_QUARTERS * _dot3(dpv, dpv)
        + _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(du, du)
        + BATTLE2_K_SPRING * _dot3(dv, dv)
        + BATTLE2_G_COUPLE * _dot3(duxv, duxv)
    )
    if c2 == _ZERO:
        return za
    lam = -c1 / c2
    return (
        _add3(u, _scale3(lam, du)),
        _add3(v, _scale3(lam, dv)),
        _add3(pu, _scale3(lam, dpu)),
        _add3(pv, _scale3(lam, dpv)),
    )


def kernel_radial_flip(z2: State4) -> State4:
    """Pure-field radial momentum reflection inside `ker(L_int) cap M_{H0}`."""
    u2, v2, pu2, pv2 = z2
    c1 = _TWO * (_dot3(pu2, u2) + _THREE_QUARTERS * _dot3(pv2, v2))
    c2 = _dot3(u2, u2) + _THREE_QUARTERS * _dot3(v2, v2)
    if c2 == _ZERO:
        return z2
    lam = -c1 / c2
    return (u2, v2, _add3(pu2, _scale3(lam, u2)), _add3(pv2, _scale3(lam, v2)))


def move_u_3d_kernel(za: State4, du: Vec3, ek: Vec3, v_ek: Any, beta1: Any) -> State4:
    """Pure-field 3D Stage-`u` Vieta chord conserving `H_int` and `L_int` for arbitrary `(du, beta1)`."""
    ua, va, pua, pva = za
    d1 = _dot3(va, _cross3(ua, du))
    alpha1 = -_dot3(va, _cross3(du, pua)) / d1
    dpu_b = _scale3(alpha1, du)
    r1 = _scale3(_NEG_ONE, _add3(_cross3(du, pua), _scale3(alpha1, _cross3(ua, du))))
    dpv_part = _scale3(_ONE / v_ek, _cross3(r1, ek))
    dpv = _add3(dpv_part, _scale3(beta1, va))
    uxv0 = _cross3(ua, va)
    duxv = _cross3(du, va)
    c1 = _TWO * (
        _dot3(pua, dpu_b)
        + _THREE_QUARTERS * _dot3(pva, dpv)
        + _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(ua, du)
        + BATTLE2_G_COUPLE * _dot3(uxv0, duxv)
    )
    c2 = (
        _dot3(dpu_b, dpu_b)
        + _THREE_QUARTERS * _dot3(dpv, dpv)
        + _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(du, du)
        + BATTLE2_G_COUPLE * _dot3(duxv, duxv)
    )
    if c2 == _ZERO:
        return za
    lam1 = -c1 / c2
    return (
        _add3(ua, _scale3(lam1, du)),
        va,
        _add3(pua, _scale3(lam1, dpu_b)),
        _add3(pva, _scale3(lam1, dpv)),
    )


def move_v_3d_kernel(za: State4, dv: Vec3, ek: Vec3, u_ek: Any, beta2: Any) -> State4:
    """Pure-field 3D Stage-`v` Vieta chord conserving `H_int` and `L_int` for arbitrary `(dv, beta2)`."""
    ua, va, pua, pva = za
    d2 = _dot3(ua, _cross3(va, dv))
    alpha2 = -_dot3(ua, _cross3(dv, pva)) / d2
    dpv_b = _scale3(alpha2, dv)
    r2 = _scale3(_NEG_ONE, _add3(_cross3(dv, pva), _scale3(alpha2, _cross3(va, dv))))
    dpu_part = _scale3(_ONE / u_ek, _cross3(r2, ek))
    dpu = _add3(dpu_part, _scale3(beta2, ua))
    uxv1 = _cross3(ua, va)
    uxdv = _cross3(ua, dv)
    c1 = _TWO * (
        _dot3(pua, dpu)
        + _THREE_QUARTERS * _dot3(pva, dpv_b)
        + BATTLE2_K_SPRING * _dot3(va, dv)
        + BATTLE2_G_COUPLE * _dot3(uxv1, uxdv)
    )
    c2 = (
        _dot3(dpu, dpu)
        + _THREE_QUARTERS * _dot3(dpv_b, dpv_b)
        + BATTLE2_K_SPRING * _dot3(dv, dv)
        + BATTLE2_G_COUPLE * _dot3(uxdv, uxdv)
    )
    if c2 == _ZERO:
        return za
    lam2 = -c1 / c2
    return (
        ua,
        _add3(va, _scale3(lam2, dv)),
        _add3(pua, _scale3(lam2, dpu)),
        _add3(pva, _scale3(lam2, dpv_b)),
    )


def move_u_planar_kernel(za: State4, n0: Vec3, du: Vec3, alpha1: Any, beta1: Any) -> State4:
    """Pure-field planar Stage-`u` Vieta chord conserving `H_int` and `L_int` on `Pi_0 = n0^perp`."""
    ua, va, pua, pva = za
    nn0 = _dot3(n0, n0)
    va_sq = _dot3(va, va)
    dpu = _scale3(alpha1, du)
    wv = _cross3(n0, va)
    l_err = _add3(_scale3(alpha1, _cross3(ua, du)), _cross3(du, pua))
    kap1 = -_dot3(n0, l_err) / (va_sq * nn0)
    dpv = _add3(_scale3(kap1, wv), _scale3(beta1, va))
    uxv0 = _cross3(ua, va)
    duxv = _cross3(du, va)
    c1 = _TWO * (
        _dot3(pua, dpu)
        + _THREE_QUARTERS * _dot3(pva, dpv)
        + _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(ua, du)
        + BATTLE2_G_COUPLE * _dot3(uxv0, duxv)
    )
    c2 = (
        _dot3(dpu, dpu)
        + _THREE_QUARTERS * _dot3(dpv, dpv)
        + _THREE_QUARTERS * BATTLE2_K_SPRING * _dot3(du, du)
        + BATTLE2_G_COUPLE * _dot3(duxv, duxv)
    )
    if c2 == _ZERO:
        return za
    lam = -c1 / c2
    return (
        _add3(ua, _scale3(lam, du)),
        va,
        _add3(pua, _scale3(lam, dpu)),
        _add3(pva, _scale3(lam, dpv)),
    )


def move_v_planar_kernel(za: State4, n0: Vec3, dv: Vec3, alpha2: Any, beta2: Any) -> State4:
    """Pure-field planar Stage-`v` Vieta chord conserving `H_int` and `L_int` on `Pi_0 = n0^perp`."""
    ua, va, pua, pva = za
    nn0 = _dot3(n0, n0)
    ua_sq = _dot3(ua, ua)
    dpv = _scale3(alpha2, dv)
    wu = _cross3(n0, ua)
    l_err = _add3(_scale3(alpha2, _cross3(va, dv)), _cross3(dv, pva))
    kap2 = -_dot3(n0, l_err) / (ua_sq * nn0)
    dpu = _add3(_scale3(kap2, wu), _scale3(beta2, ua))
    uxv0 = _cross3(ua, va)
    uxdv = _cross3(ua, dv)
    c1 = _TWO * (
        _dot3(pua, dpu)
        + _THREE_QUARTERS * _dot3(pva, dpv)
        + BATTLE2_K_SPRING * _dot3(va, dv)
        + BATTLE2_G_COUPLE * _dot3(uxv0, uxdv)
    )
    c2 = (
        _dot3(dpu, dpu)
        + _THREE_QUARTERS * _dot3(dpv, dpv)
        + BATTLE2_K_SPRING * _dot3(dv, dv)
        + BATTLE2_G_COUPLE * _dot3(uxdv, uxdv)
    )
    if c2 == _ZERO:
        return za
    lam = -c1 / c2
    return (
        ua,
        _add3(va, _scale3(lam, dv)),
        _add3(pua, _scale3(lam, dpu)),
        _add3(pva, _scale3(lam, dpv)),
    )


def stage3_kernel_chord_exact(z2: State4, xi1: Any, xi2: Any, xi3: Any) -> State4:
    """Pure-field Stage-3 Vieta chord inside `K(u2, v2)` conserving `H_int` and `L_int` identically."""
    u2, v2, pu2, pv2 = z2
    yu = _add3(_scale3(xi1, u2), _scale3(xi2, v2))
    yv = _add3(_scale3(xi2, u2), _scale3(xi3, v2))
    c1 = _TWO * (_dot3(pu2, yu) + _THREE_QUARTERS * _dot3(pv2, yv))
    c2 = _dot3(yu, yu) + _THREE_QUARTERS * _dot3(yv, yv)
    if c2 == _ZERO:
        return z2
    lam = -c1 / c2
    return (u2, v2, _add3(pu2, _scale3(lam, yu)), _add3(pv2, _scale3(lam, yv)))


def _solve_k_subspace_coords(
    uu: Any, ww: Any, vv: Any, pu: Vec3, pv: Vec3, u_vec: Vec3, v_vec: Vec3
) -> Vec3:
    """Exact `T`-orthogonal projection coordinates `xi = M^-1 b` onto `K(u, v)`."""
    det_g = uu * vv - ww * ww
    s_val = _THREE_QUARTERS * uu + vv
    denom = _THREE_QUARTERS * s_val * det_g
    b1 = _dot3(pu, u_vec)
    b2 = _dot3(pu, v_vec) + _THREE_QUARTERS * _dot3(pv, u_vec)
    b3 = _THREE_QUARTERS * _dot3(pv, v_vec)
    c11 = _THREE_QUARTERS * vv * (vv + _THREE_QUARTERS * uu) - Fraction(9, 16) * ww * ww
    c12 = -_THREE_QUARTERS * vv * ww
    c13 = _THREE_QUARTERS * ww * ww
    c22 = _THREE_QUARTERS * uu * vv
    c23 = -_THREE_QUARTERS * uu * ww
    c33 = det_g + _THREE_QUARTERS * uu * uu
    xi1 = (c11 * b1 + c12 * b2 + c13 * b3) / denom
    xi2 = (c12 * b1 + c22 * b2 + c23 * b3) / denom
    xi3 = (c13 * b1 + c23 * b2 + c33 * b3) / denom
    return (xi1, xi2, xi3)


def _m_quad_norm(uu: Any, ww: Any, vv: Any, xi: Vec3) -> Any:
    x1, x2, x3 = xi
    yu_sq = uu * x1 * x1 + _TWO * ww * x1 * x2 + vv * x2 * x2
    yv_sq = uu * x2 * x2 + _TWO * ww * x2 * x3 + vv * x3 * x3
    return yu_sq + _THREE_QUARTERS * yv_sq


def stage3_match_momenta(
    z2: State4, pu_hat: Vec3, pv_hat: Vec3, q_fn: Callable[[Fraction], Fraction]
) -> State4:
    """Radially-normalized `T`-orthogonal Stage-3 Vieta matching on `K(u2, v2)`."""
    u2, v2, pu2, pv2 = z2
    u2_g = _qvec(u2, quantize_sig29)
    v2_g = _qvec(v2, quantize_sig29)
    pu2_g = _qvec(pu2, quantize_sig29)
    pv2_g = _qvec(pv2, quantize_sig29)
    uu, ww, vv = _dot3(u2_g, u2_g), _dot3(u2_g, v2_g), _dot3(v2_g, v2_g)
    det_g = uu * vv - ww * ww
    if det_g <= _ZERO:
        return z2
    xi2 = _solve_k_subspace_coords(uu, ww, vv, pu2_g, pv2_g, u2_g, v2_g)
    xi_h = _solve_k_subspace_coords(uu, ww, vv, pu_hat, pv_hat, u2_g, v2_g)
    r2_sq = _m_quad_norm(uu, ww, vv, xi2)
    rh_sq = _m_quad_norm(uu, ww, vv, xi_h)
    if r2_sq <= _ZERO or rh_sq <= _ZERO:
        return z2
    rho = pade_sqrt_fine(r2_sq / rh_sq)
    dxi = (rho * xi_h[0] - xi2[0], rho * xi_h[1] - xi2[1], rho * xi_h[2] - xi2[2])
    d_sq = _m_quad_norm(uu, ww, vv, dxi)
    if d_sq <= Fraction(1, 4194304) * r2_sq:
        return z2
    return stage3_kernel_chord_exact(z2, q_fn(dxi[0]), q_fn(dxi[1]), q_fn(dxi[2]))


def _planar_basis_quantize(
    w_target: Vec3, ref_vec: Vec3, u0: Vec3, v0: Vec3, nn0: Fraction, q_fn: Callable[[Fraction], Fraction]
) -> Vec3:
    diff = _sub3(w_target, ref_vec)
    uu = _dot3(u0, u0)
    uv = _dot3(u0, v0)
    vv = _dot3(v0, v0)
    du = _dot3(diff, u0)
    dv = _dot3(diff, v0)
    ca = q_fn((vv * du - uv * dv) / nn0)
    cb = q_fn((uu * dv - uv * du) / nn0)
    return _add3(ref_vec, _add3(_scale3(ca, u0), _scale3(cb, v0)))


def project_planar_b2(z_ref: State4, z_hat: State4, n0: Vec3) -> State4:
    """Exact in-plane `(u0, v0)`-basis Vieta chart for planar 3-body submanifolds (`L_0 || n0`)."""
    u0, v0, _, _ = z_ref
    nn0 = _dot3(n0, n0)
    u_g = _planar_basis_quantize(z_hat[0], u0, u0, v0, nn0, quantize_sig14)
    v_g = _planar_basis_quantize(z_hat[1], v0, u0, v0, nn0, quantize_sig14)
    pu_t = _sub3(z_hat[2], _scale3(_dot3(z_hat[2], n0) / nn0, n0))
    pv_t = _sub3(z_hat[3], _scale3(_dot3(z_hat[3], n0) / nn0, n0))

    c_a = transvection_vieta_kernel(z_ref, (_ONE, _ZERO, _HALF, _ZERO))
    c_b = apply_dem_gauge(c_a, _NEG_ONE, _NEG_ONE, _HALF)
    za0 = c_b if (_dot3(c_b[0], c_b[0]) > Fraction(1, 16) and _dot3(c_b[1], c_b[1]) > Fraction(1, 16)) else c_a
    ua, va, pua, pva = za0
    h0 = _h_int_arch(ua, va, pua, pva)

    du = _sub3(u_g, ua)
    du_sq = _dot3(du, du)
    if du_sq > _ZERO:
        alpha1 = quantize_sig11(_dot3(_sub3(pu_t, pua), du) / du_sq)
        va_sq = _dot3(va, va)
        wv = _cross3(n0, va)
        l_err1 = _add3(_scale3(alpha1, _cross3(ua, du)), _cross3(du, pua))
        kap1 = -_dot3(n0, l_err1) / (va_sq * nn0)
        pu_at1 = _add3(pua, _scale3(alpha1, du))
        pv_base1 = _add3(pva, _scale3(kap1, wv))
        g0 = _dot3(pv_base1, va) / va_sq
        pv_orth1 = _sub3(pv_base1, _scale3(g0, va))
        rem1 = h0 - _v2_arch_uv(u_g, va) - _dot3(pu_at1, pu_at1) - _THREE_QUARTERS * _dot3(pv_orth1, pv_orth1)
        g_pos1 = pade_sqrt_positive(rem1 / (_THREE_QUARTERS * va_sq))
        g_est1 = g_pos1 if _dot3(pv_t, va) >= _ZERO else -g_pos1
        beta1 = quantize_sig11(g_est1 - g0)
        s1 = move_u_planar_kernel(za0, n0, du, alpha1, beta1)
    else:
        s1 = za0

    u1, v1, pu1, pv1 = s1
    dv = _sub3(v_g, v1)
    dv_sq = _dot3(dv, dv)
    if dv_sq > _ZERO:
        alpha2 = quantize_sig14(_dot3(_sub3(pv_t, pv1), dv) / dv_sq)
        u1_sq = _dot3(u1, u1)
        wu = _cross3(n0, u1)
        l_err2 = _add3(_scale3(alpha2, _cross3(v1, dv)), _cross3(dv, pv1))
        kap2 = -_dot3(n0, l_err2) / (u1_sq * nn0)
        pv_at2 = _add3(pv1, _scale3(alpha2, dv))
        pu_base2 = _add3(pu1, _scale3(kap2, wu))
        g0_2 = _dot3(pu_base2, u1) / u1_sq
        pu_orth2 = _sub3(pu_base2, _scale3(g0_2, u1))
        rem2 = h0 - _v2_arch_uv(u1, v_g) - _THREE_QUARTERS * _dot3(pv_at2, pv_at2) - _dot3(pu_orth2, pu_orth2)
        g_pos2 = pade_sqrt_positive(rem2 / u1_sq)
        g_est2 = g_pos2 if _dot3(pu_t, u1) >= _ZERO else -g_pos2
        beta2 = quantize_sig14(g_est2 - g0_2)
        s2 = move_v_planar_kernel(s1, n0, dv, alpha2, beta2)
    else:
        s2 = s1

    return stage3_match_momenta(s2, pu_t, pv_t, quantize_sig14)


def b2_local_chart_kernel(
    z_ref: State4,
    z_hat: State4,
    ek1: Vec3,
    ek2: Vec3,
    q_fn: Callable[[Any], Any] = quantize_sig17,
) -> State4:
    """Pure-field local reverse-transvection chart near `z_ref` (nonsingular at `dt = 0` over `QQJet`)."""
    u0, v0, _, _ = z_ref
    u_t, v_t, pu_t, pv_t = z_hat
    u_g = _add3(u0, _qvec(_sub3(u_t, u0), q_fn))
    v_g = _add3(v0, _qvec(_sub3(v_t, v0), q_fn))

    za_u = transvection_vieta_kernel(z_ref, (_ONE, _ZERO, _NEG_ONE, _ZERO))
    ua, va, pua, pva = za_u
    h0 = _h_int_arch(ua, va, pua, pva)
    du = _sub3(u_g, ua)
    d1 = _dot3(va, _cross3(ua, du))
    alpha1 = -_dot3(va, _cross3(du, pua)) / d1
    dpu_b = _scale3(alpha1, du)
    r1 = _scale3(_NEG_ONE, _add3(_cross3(du, pua), _scale3(alpha1, _cross3(ua, du))))
    v_ek1 = _dot3(va, ek1)
    dpv_part = _scale3(_ONE / v_ek1, _cross3(r1, ek1))
    pu_at1 = _add3(pua, dpu_b)
    pv_base1 = _add3(pva, dpv_part)
    va_sq = _dot3(va, va)
    g0_1 = _dot3(pv_base1, va) / va_sq
    pv_orth1 = _sub3(pv_base1, _scale3(g0_1, va))
    rem1 = h0 - _v2_arch_uv(u_g, va) - _dot3(pu_at1, pu_at1) - _THREE_QUARTERS * _dot3(pv_orth1, pv_orth1)
    g_sq1 = rem1 / (_THREE_QUARTERS * va_sq)
    g_tgt1 = _dot3(pv_t, va) / va_sq
    s1_mid = g_tgt1 * pade22_sqrt_from_one(g_sq1 / (g_tgt1 * g_tgt1))
    g_est1 = s1_mid * pade22_sqrt_from_one(g_sq1 / (s1_mid * s1_mid))
    beta1 = q_fn(g_est1 - g0_1)
    s1 = move_u_3d_kernel(za_u, du, ek1, v_ek1, beta1)

    za_v = transvection_vieta_kernel(s1, (_ZERO, _ONE, _ZERO, _NEG_ONE))
    ua2, va2, pua2, pva2 = za_v
    dv = _sub3(v_g, va2)
    d2 = _dot3(ua2, _cross3(va2, dv))
    alpha2 = -_dot3(ua2, _cross3(dv, pva2)) / d2
    dpv_b = _scale3(alpha2, dv)
    r2 = _scale3(_NEG_ONE, _add3(_cross3(dv, pva2), _scale3(alpha2, _cross3(va2, dv))))
    u_ek2 = _dot3(ua2, ek2)
    dpu_part = _scale3(_ONE / u_ek2, _cross3(r2, ek2))
    pv_at2 = _add3(pva2, dpv_b)
    pu_base2 = _add3(pua2, dpu_part)
    ua2_sq = _dot3(ua2, ua2)
    g0_2 = _dot3(pu_base2, ua2) / ua2_sq
    pu_orth2 = _sub3(pu_base2, _scale3(g0_2, ua2))
    rem2 = h0 - _v2_arch_uv(ua2, v_g) - _THREE_QUARTERS * _dot3(pv_at2, pv_at2) - _dot3(pu_orth2, pu_orth2)
    g_sq2 = rem2 / ua2_sq
    g_tgt2 = _dot3(pu_t, ua2) / ua2_sq
    s2_mid = g_tgt2 * pade22_sqrt_from_one(g_sq2 / (g_tgt2 * g_tgt2))
    g_est2 = s2_mid * pade22_sqrt_from_one(g_sq2 / (s2_mid * s2_mid))
    beta2 = q_fn(g_est2 - g0_2)
    s2 = move_v_3d_kernel(za_v, dv, ek2, u_ek2, beta2)

    s2_f = kernel_radial_flip(s2)
    u2, v2, pu2, pv2 = s2_f
    uu, ww, vv = _dot3(u2, u2), _dot3(u2, v2), _dot3(v2, v2)
    xi_f = _solve_k_subspace_coords(uu, ww, vv, pu2, pv2, u2, v2)
    xi_s2 = _solve_k_subspace_coords(uu, ww, vv, s2[2], s2[3], u2, v2)
    xi_h = _solve_k_subspace_coords(uu, ww, vv, pu_t, pv_t, u2, v2)
    r2_sq = _m_quad_norm(uu, ww, vv, xi_f)
    rh_sq = _m_quad_norm(uu, ww, vv, xi_h)
    rho_1 = pade22_sqrt_from_one(r2_sq / rh_sq)
    rho = rho_1 * pade22_sqrt_from_one((r2_sq / rh_sq) / (rho_1 * rho_1))
    dxi_0 = (xi_s2[0] - xi_f[0], xi_s2[1] - xi_f[1], xi_s2[2] - xi_f[2])
    dxi_raw = (rho * xi_h[0] - xi_f[0], rho * xi_h[1] - xi_f[1], rho * xi_h[2] - xi_f[2])
    scale_inv = -_ONE / dxi_0[0]
    y1 = -_ONE + q_fn(scale_inv * dxi_raw[0] - (-_ONE))
    y2 = q_fn(scale_inv * dxi_raw[1])
    y3 = -_ONE + q_fn(scale_inv * dxi_raw[2] - (-_ONE))
    return stage3_kernel_chord_exact(s2_f, y1, y2, y3)


def _near_zref_runtime(z_ref: State4, z_hat: State4) -> State4:
    ek1, _ = _pick_ek(z_ref[1])
    ek2, _ = _pick_ek(z_ref[0])
    return b2_local_chart_kernel(z_ref, z_hat, ek1, ek2, quantize_sig17)


def _eval_uv_fast8(za_q: State4, u_q: Vec3, v_q: Vec3, sgn1: Fraction) -> Fraction:
    ua, va, pua, pva = za_q
    du = _sub3(u_q, ua)
    dv = _sub3(v_q, va)
    du_sq = _dot3(du, du)
    dv_sq = _dot3(dv, dv)
    if du_sq == _ZERO or dv_sq == _ZERO:
        return _NEG_ONE
    d1 = _dot3(va, _cross3(ua, du))
    d2 = _dot3(u_q, _cross3(va, dv))
    if d1 * d1 < Fraction(1, 20000) * du_sq or d2 * d2 < Fraction(1, 20000) * dv_sq:
        return _NEG_ONE
    h0 = _h_int_arch(ua, va, pua, pva)
    alpha1 = -_dot3(va, _cross3(du, pua)) / d1
    dpu_b = _scale3(alpha1, du)
    r1 = _scale3(_NEG_ONE, _add3(_cross3(du, pua), _scale3(alpha1, _cross3(ua, du))))
    ek, v_ek = _pick_ek(va)
    if v_ek == _ZERO:
        return _NEG_ONE
    dpv_part = _scale3(_ONE / v_ek, _cross3(r1, ek))
    pu_at1 = _add3(pua, dpu_b)
    pv_base1 = _add3(pva, dpv_part)
    va_sq = _dot3(va, va)
    g0 = _dot3(pv_base1, va) / va_sq
    pv_orth1 = _sub3(pv_base1, _scale3(g0, va))
    rem1 = h0 - _v2_arch_uv(u_q, va) - _dot3(pu_at1, pu_at1) - _THREE_QUARTERS * _dot3(pv_orth1, pv_orth1)
    if rem1 <= Fraction(1, 2000):
        return _NEG_ONE
    g_est = sgn1 * pade_sqrt_positive(rem1 / (_THREE_QUARTERS * va_sq))
    pv1 = _add3(pv_orth1, _scale3(g_est, va))
    alpha2 = -_dot3(u_q, _cross3(dv, pv1)) / d2
    dpv_b2 = _scale3(alpha2, dv)
    r2 = _scale3(_NEG_ONE, _add3(_cross3(dv, pv1), _scale3(alpha2, _cross3(va, dv))))
    ek2, u_ek2 = _pick_ek(u_q)
    if u_ek2 == _ZERO:
        return _NEG_ONE
    dpu_part2 = _scale3(_ONE / u_ek2, _cross3(r2, ek2))
    pv_at2 = _add3(pv1, dpv_b2)
    pu_base2 = _add3(pu_at1, dpu_part2)
    ut_sq = _dot3(u_q, u_q)
    g0_2 = _dot3(pu_base2, u_q) / ut_sq
    pu_orth2 = _sub3(pu_base2, _scale3(g0_2, u_q))
    rem2 = h0 - _v2_arch_uv(u_q, v_q) - _THREE_QUARTERS * _dot3(pv_at2, pv_at2) - _dot3(pu_orth2, pu_orth2)
    if rem2 <= Fraction(1, 2000):
        return _NEG_ONE
    return rem1 if rem1 < rem2 else rem2


def _eval_vu_fast8(za_q: State4, u_q: Vec3, v_q: Vec3, sgn1: Fraction) -> Fraction:
    ua, va, pua, pva = za_q
    du = _sub3(u_q, ua)
    dv = _sub3(v_q, va)
    du_sq = _dot3(du, du)
    dv_sq = _dot3(dv, dv)
    if du_sq == _ZERO or dv_sq == _ZERO:
        return _NEG_ONE
    d1 = _dot3(ua, _cross3(va, dv))
    d2 = _dot3(v_q, _cross3(ua, du))
    if d1 * d1 < Fraction(1, 20000) * dv_sq or d2 * d2 < Fraction(1, 20000) * du_sq:
        return _NEG_ONE
    h0 = _h_int_arch(ua, va, pua, pva)
    alpha1 = -_dot3(ua, _cross3(dv, pva)) / d1
    dpv_b = _scale3(alpha1, dv)
    r1 = _scale3(_NEG_ONE, _add3(_cross3(dv, pva), _scale3(alpha1, _cross3(va, dv))))
    ek, u_ek = _pick_ek(ua)
    if u_ek == _ZERO:
        return _NEG_ONE
    dpu_part = _scale3(_ONE / u_ek, _cross3(r1, ek))
    pv_at1 = _add3(pva, dpv_b)
    pu_base1 = _add3(pua, dpu_part)
    ua_sq = _dot3(ua, ua)
    g0 = _dot3(pu_base1, ua) / ua_sq
    pu_orth1 = _sub3(pu_base1, _scale3(g0, ua))
    rem1 = h0 - _v2_arch_uv(ua, v_q) - _THREE_QUARTERS * _dot3(pv_at1, pv_at1) - _dot3(pu_orth1, pu_orth1)
    if rem1 <= Fraction(1, 2000):
        return _NEG_ONE
    g_est = sgn1 * pade_sqrt_positive(rem1 / ua_sq)
    pu1 = _add3(pu_orth1, _scale3(g_est, ua))
    alpha2 = -_dot3(v_q, _cross3(du, pu1)) / d2
    dpu_b2 = _scale3(alpha2, du)
    r2 = _scale3(_NEG_ONE, _add3(_cross3(du, pu1), _scale3(alpha2, _cross3(ua, du))))
    ek2, v_ek2 = _pick_ek(v_q)
    if v_ek2 == _ZERO:
        return _NEG_ONE
    dpv_part2 = _scale3(_ONE / v_ek2, _cross3(r2, ek2))
    pu_at2 = _add3(pu1, dpu_b2)
    pv_base2 = _add3(pv_at1, dpv_part2)
    vt_sq = _dot3(v_q, v_q)
    g0_2 = _dot3(pv_base2, v_q) / vt_sq
    pv_orth2 = _sub3(pv_base2, _scale3(g0_2, v_q))
    rem2 = h0 - _v2_arch_uv(u_q, v_q) - _dot3(pu_at2, pu_at2) - _THREE_QUARTERS * _dot3(pv_orth2, pv_orth2)
    if rem2 <= Fraction(1, 2000):
        return _NEG_ONE
    return rem1 if rem1 < rem2 else rem2


BestCfg = Tuple[State4, bool, Fraction, Fraction]


def _chk_anchor1(
    z_b: State4,
    is_rot: bool,
    eu: Fraction,
    ev: Fraction,
    et: Fraction,
    u_q: Vec3,
    v_q: Vec3,
    best_m: Fraction,
    best_cfg: BestCfg,
) -> Tuple[Fraction, BestCfg]:
    if best_m > Fraction(1, 20):
        return (best_m, best_cfg)
    za = apply_dem_gauge(z_b, eu, ev, et)
    za_q = _q_state8(za)
    m1 = _eval_uv_fast8(za_q, u_q, v_q, _ONE)
    if m1 > best_m:
        best_m, best_cfg = m1, (za, is_rot, _ZERO, _ONE)
        if m1 > Fraction(1, 20):
            return (best_m, best_cfg)
    m2 = _eval_uv_fast8(za_q, u_q, v_q, _NEG_ONE)
    if m2 > best_m:
        best_m, best_cfg = m2, (za, is_rot, _ZERO, _NEG_ONE)
        if m2 > Fraction(1, 20):
            return (best_m, best_cfg)
    m3 = _eval_vu_fast8(za_q, u_q, v_q, _ONE)
    if m3 > best_m:
        best_m, best_cfg = m3, (za, is_rot, _ONE, _ONE)
        if m3 > Fraction(1, 20):
            return (best_m, best_cfg)
    m4 = _eval_vu_fast8(za_q, u_q, v_q, _NEG_ONE)
    if m4 > best_m:
        best_m, best_cfg = m4, (za, is_rot, _ONE, _NEG_ONE)
    return (best_m, best_cfg)


def _chk_dem5(
    z_b: State4,
    is_rot: bool,
    eu: Fraction,
    ev: Fraction,
    u_q: Vec3,
    v_q: Vec3,
    best_m: Fraction,
    best_cfg: BestCfg,
) -> Tuple[Fraction, BestCfg]:
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _ZERO, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _HALF, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _NEG_HALF, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _NEG_ONE, u_q, v_q, best_m, best_cfg)
    return (best_m, best_cfg)


def _chk_dem5_extra(
    z_b: State4,
    is_rot: bool,
    eu: Fraction,
    ev: Fraction,
    u_q: Vec3,
    v_q: Vec3,
    best_m: Fraction,
    best_cfg: BestCfg,
) -> Tuple[Fraction, BestCfg]:
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _QUARTER, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _NEG_QUARTER, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _TWO, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _NEG_TWO, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_anchor1(z_b, is_rot, eu, ev, _THREE_HALVES, u_q, v_q, best_m, best_cfg)
    return (best_m, best_cfg)


def _chk_rot_base(
    r_b: State4, u_q: Vec3, v_q: Vec3, best_m: Fraction, best_cfg: BestCfg
) -> Tuple[Fraction, BestCfg]:
    if best_m > Fraction(1, 100):
        return (best_m, best_cfg)
    best_m, best_cfg = _chk_dem5(r_b, True, _ONE, _ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5(r_b, True, _ONE, _NEG_ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5(r_b, True, _NEG_ONE, _NEG_ONE, u_q, v_q, best_m, best_cfg)
    return (best_m, best_cfg)


def _search_rot_only(z_ref: State4, u_q: Vec3, v_q: Vec3) -> Tuple[Fraction, BestCfg]:
    best_m = _NEG_ONE
    l_vec = _l_int_q(z_ref[0], z_ref[1], z_ref[2], z_ref[3])
    r_ht = apply_l_half_turn_q(z_ref, l_vec)
    best_cfg: BestCfg = (r_ht, True, _ZERO, _ONE)
    best_m, best_cfg = _chk_rot_base(r_ht, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_rot_base(apply_l_rot_q(z_ref, l_vec, _QUARTER), u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_rot_base(apply_l_rot_q(z_ref, l_vec, _NEG_QUARTER), u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_rot_base(apply_l_rot_q(z_ref, l_vec, _HALF), u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_rot_base(apply_l_rot_q(z_ref, l_vec, _NEG_HALF), u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_rot_base(apply_l_rot_q(z_ref, l_vec, _ONE), u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_rot_base(apply_l_rot_q(z_ref, l_vec, _NEG_ONE), u_q, v_q, best_m, best_cfg)
    return (best_m, best_cfg)


def _search_full_atlas(z_ref: State4, u_q: Vec3, v_q: Vec3) -> Tuple[Fraction, BestCfg]:
    best_m = _NEG_ONE
    best_cfg: BestCfg = (z_ref, False, _ZERO, _ONE)
    best_m, best_cfg = _chk_dem5(z_ref, False, _ONE, _ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5(z_ref, False, _NEG_ONE, _NEG_ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5(z_ref, False, _ONE, _NEG_ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5(z_ref, False, _NEG_ONE, _ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5_extra(z_ref, False, _ONE, _ONE, u_q, v_q, best_m, best_cfg)
    best_m, best_cfg = _chk_dem5_extra(z_ref, False, _ONE, _NEG_ONE, u_q, v_q, best_m, best_cfg)
    if best_m > _ZERO:
        return (best_m, best_cfg)
    return _search_rot_only(z_ref, u_q, v_q)


def _exec_3d_chart(
    za: State4,
    is_rot: bool,
    fine_q: bool,
    u_g: Vec3,
    v_g: Vec3,
    pu_t: Vec3,
    pv_t: Vec3,
    mode: Fraction,
    sgn1: Fraction,
) -> State4:
    q1_fn = quantize_sig29 if fine_q else (quantize_sig8 if is_rot else quantize_sig11)
    q2_fn = quantize_sig29 if fine_q else (quantize_sig8 if is_rot else quantize_sig11)
    rad_fn = pade_sqrt_fine if fine_q else pade_sqrt_positive
    qg0_fn = quantize_sig29 if fine_q else (quantize_sig8 if is_rot else quantize_sig14)

    if mode == _ZERO:
        ua, va, pua, pva = za
        h0 = _h_int_arch(ua, va, pua, pva)
        du = _sub3(u_g, ua)
        d1 = _dot3(va, _cross3(ua, du))
        alpha1 = -_dot3(va, _cross3(du, pua)) / d1
        dpu_b = _scale3(alpha1, du)
        r1 = _scale3(_NEG_ONE, _add3(_cross3(du, pua), _scale3(alpha1, _cross3(ua, du))))
        ek1, v_ek1 = _pick_ek(va)
        dpv_part = _scale3(_ONE / v_ek1, _cross3(r1, ek1))
        pu_at1 = _add3(pua, dpu_b)
        pv_base1 = _add3(pva, dpv_part)
        va_sq = _dot3(va, va)
        g0_1 = _dot3(pv_base1, va) / va_sq
        pv_orth1 = _sub3(pv_base1, _scale3(g0_1, va))
        rem1 = h0 - _v2_arch_uv(u_g, va) - _dot3(pu_at1, pu_at1) - _THREE_QUARTERS * _dot3(pv_orth1, pv_orth1)
        g_est1 = sgn1 * rad_fn(rem1 / (_THREE_QUARTERS * va_sq))
        beta1 = q1_fn(g_est1 - qg0_fn(g0_1))
        s1 = move_u_3d_kernel(za, du, ek1, v_ek1, beta1)

        ua2, va2, pua2, pva2 = s1
        dv = _sub3(v_g, va2)
        d2 = _dot3(ua2, _cross3(va2, dv))
        alpha2 = -_dot3(ua2, _cross3(dv, pva2)) / d2
        dpv_b = _scale3(alpha2, dv)
        r2 = _scale3(_NEG_ONE, _add3(_cross3(dv, pva2), _scale3(alpha2, _cross3(va2, dv))))
        ek2, u_ek2 = _pick_ek(ua2)
        dpu_part = _scale3(_ONE / u_ek2, _cross3(r2, ek2))
        pv_at2 = _add3(pva2, dpv_b)
        pu_base2 = _add3(pua2, dpu_part)
        ua2_sq = _dot3(ua2, ua2)
        g0_2 = _dot3(pu_base2, ua2) / ua2_sq
        pu_orth2 = _sub3(pu_base2, _scale3(g0_2, ua2))
        rem2 = h0 - _v2_arch_uv(ua2, v_g) - _THREE_QUARTERS * _dot3(pv_at2, pv_at2) - _dot3(pu_orth2, pu_orth2)
        g_pos2 = rad_fn(rem2 / ua2_sq)
        g_est2 = g_pos2 if _dot3(pu_t, ua2) >= _ZERO else -g_pos2
        beta2 = q2_fn(g_est2 - qg0_fn(g0_2))
        s2 = move_v_3d_kernel(s1, dv, ek2, u_ek2, beta2)
    else:
        ua, va, pua, pva = za
        h0 = _h_int_arch(ua, va, pua, pva)
        dv = _sub3(v_g, va)
        d1 = _dot3(ua, _cross3(va, dv))
        alpha1 = -_dot3(ua, _cross3(dv, pva)) / d1
        dpv_b = _scale3(alpha1, dv)
        r1 = _scale3(_NEG_ONE, _add3(_cross3(dv, pva), _scale3(alpha1, _cross3(va, dv))))
        ek1, u_ek1 = _pick_ek(ua)
        dpu_part = _scale3(_ONE / u_ek1, _cross3(r1, ek1))
        pv_at1 = _add3(pva, dpv_b)
        pu_base1 = _add3(pua, dpu_part)
        ua_sq = _dot3(ua, ua)
        g0_1 = _dot3(pu_base1, ua) / ua_sq
        pu_orth1 = _sub3(pu_base1, _scale3(g0_1, ua))
        rem1 = h0 - _v2_arch_uv(ua, v_g) - _THREE_QUARTERS * _dot3(pv_at1, pv_at1) - _dot3(pu_orth1, pu_orth1)
        g_est1 = sgn1 * rad_fn(rem1 / ua_sq)
        beta1 = q1_fn(g_est1 - qg0_fn(g0_1))
        s1 = move_v_3d_kernel(za, dv, ek1, u_ek1, beta1)

        ua2, va2, pua2, pva2 = s1
        du = _sub3(u_g, ua2)
        d2 = _dot3(va2, _cross3(ua2, du))
        alpha2 = -_dot3(va2, _cross3(du, pua2)) / d2
        dpu_b = _scale3(alpha2, du)
        r2 = _scale3(_NEG_ONE, _add3(_cross3(du, pua2), _scale3(alpha2, _cross3(ua2, du))))
        ek2, v_ek2 = _pick_ek(va2)
        dpv_part = _scale3(_ONE / v_ek2, _cross3(r2, ek2))
        pu_at2 = _add3(pua2, dpu_b)
        pv_base2 = _add3(pva2, dpv_part)
        va2_sq = _dot3(va2, va2)
        g0_2 = _dot3(pv_base2, va2) / va2_sq
        pv_orth2 = _sub3(pv_base2, _scale3(g0_2, va2))
        rem2 = h0 - _v2_arch_uv(u_g, va2) - _dot3(pu_at2, pu_at2) - _THREE_QUARTERS * _dot3(pv_orth2, pv_orth2)
        g_pos2 = rad_fn(rem2 / (_THREE_QUARTERS * va2_sq))
        g_est2 = g_pos2 if _dot3(pv_t, va2) >= _ZERO else -g_pos2
        beta2 = q2_fn(g_est2 - qg0_fn(g0_2))
        s2 = move_u_3d_kernel(s1, du, ek2, v_ek2, beta2)

    if is_rot and not fine_q:
        return s2
    q3_fn = quantize_sig17 if fine_q else quantize_sig11
    return stage3_match_momenta(s2, pu_t, pv_t, q3_fn)


def project_battle2_manifold(z_ref: State4, z_hat: State4) -> State4:
    """Project Jacobi internal state `z_hat` onto `M_{H0, L0}` over `Q` across all planar & 3D regimes."""
    if z_hat == z_ref:
        return z_ref
    u0, v0, pu0, pv0 = z_ref
    n0 = _cross3(u0, v0)
    nn0 = _dot3(n0, n0)
    if nn0 > _ZERO and _dot3(pu0, n0) == _ZERO and _dot3(pv0, n0) == _ZERO:
        return project_planar_b2(z_ref, z_hat, n0)

    u_t, v_t, pu_t, pv_t = z_hat
    du0 = _sub3(u_t, u0)
    dv0 = _sub3(v_t, v0)
    dist_sq = _dot3(du0, du0) + _dot3(dv0, dv0)

    if dist_sq < Fraction(1, 1048576):
        return _near_zref_runtime(z_ref, z_hat)

    if dist_sq < Fraction(1, 64):
        u_g = _add3(u0, _qvec(du0, quantize_sig29))
        v_g = _add3(v0, _qvec(dv0, quantize_sig29))
        u_q = _qvec(u_g, quantize_sig8)
        v_q = _qvec(v_g, quantize_sig8)
        _, best_cfg = _search_rot_only(z_ref, u_q, v_q)
        za_win, is_rot, mode_win, sgn_win = best_cfg
        return _exec_3d_chart(za_win, is_rot, True, u_g, v_g, pu_t, pv_t, mode_win, sgn_win)

    h0 = _h_int_arch(u0, v0, pu0, pv0)
    v_cur = _v2_arch_uv(u_t, v_t)
    k_hat = _dot3(pu_t, pu_t) + _THREE_QUARTERS * _dot3(pv_t, pv_t)
    v_cap = h0 - (k_hat if k_hat >= Fraction(1, 40) else Fraction(1, 40))
    if v_cur > v_cap and v_cap > _ZERO:
        s_t = _THREE_QUARTERS * _dot3(u_t, u_t) + _dot3(v_t, v_t)
        uxv_t = _cross3(u_t, v_t)
        a_t = _dot3(uxv_t, uxv_t)
        disc = pade_sqrt_positive(s_t * s_t + Fraction(4, 1) * a_t * v_cap)
        ps_sq = (_TWO * v_cap) / (s_t + disc)
        pull_scale = quantize_sig14(pade_sqrt_positive(ps_sq))
        u_eff = _scale3(pull_scale, u_t)
        v_eff = _scale3(pull_scale, v_t)
    else:
        u_eff, v_eff = u_t, v_t

    u_g14 = _qvec(u_eff, quantize_sig14)
    v_g14 = _qvec(v_eff, quantize_sig14)
    u_q = _qvec(u_g14, quantize_sig8)
    v_q = _qvec(v_g14, quantize_sig8)

    best_m, best_cfg = _search_full_atlas(z_ref, u_q, v_q)
    if best_m <= _ZERO:
        u_in = _scale3(Fraction(15, 16), u_g14)
        v_in = _scale3(Fraction(15, 16), v_g14)
        best_m, best_cfg = _search_full_atlas(
            z_ref, _qvec(u_in, quantize_sig8), _qvec(v_in, quantize_sig8)
        )
        u_eff, v_eff = u_in, v_in

    za_win, is_rot, mode_win, sgn_win = best_cfg
    u_g = _qvec(u_eff, quantize_sig8 if is_rot else quantize_sig14)
    v_g = _qvec(v_eff, quantize_sig8 if is_rot else quantize_sig14)
    return _exec_3d_chart(za_win, is_rot, False, u_g, v_g, pu_t, pv_t, mode_win, sgn_win)


def build_state_from_qp_battle2(q: Sequence[Any], p: Sequence[Any]) -> Any:
    """Construct minimal globally-injective Battle 2 Jacobi state `(rc, pt, uh, vh, puh, pvh, z_ref)`."""
    rc, pt, u, v, pu, pv = _cartesian_to_jacobi(q, p)
    return (rc, pt, u, v, pu, pv, (u, v, pu, pv))


def initial_state_battle2(dt: Optional[Fraction] = None) -> Any:
    """Return canonical Battle 2 initial state."""
    return build_state_from_qp_battle2(BATTLE2_Q0, BATTLE2_P0)


def step_battle2(
    state: Any, dt: Any, q_fn: Callable[[Any], Any] = quantize_sig29
) -> Any:
    """Single-step symmetric conservative shadow Jacobi Verlet with exact `0/1` time-reversibility."""
    rc, pt, uh, vh, puh, pvh, z_ref = state
    rc_n = _add3(rc, _scale3(Fraction(1, 3) * dt, pt))
    f_ru, f_rv = b2_shadow_force_uv(z_ref[0], z_ref[1], dt)
    fu0, fv0 = b2_shadow_force_uv(uh, vh, dt)
    half_dt = dt * _HALF
    puh_m = (
        puh[0] + half_dt * (f_ru[0] + q_fn(fu0[0] - f_ru[0])),
        puh[1] + half_dt * (f_ru[1] + q_fn(fu0[1] - f_ru[1])),
        puh[2] + half_dt * (f_ru[2] + q_fn(fu0[2] - f_ru[2])),
    )
    pvh_m = (
        pvh[0] + half_dt * (f_rv[0] + q_fn(fv0[0] - f_rv[0])),
        pvh[1] + half_dt * (f_rv[1] + q_fn(fv0[1] - f_rv[1])),
        pvh[2] + half_dt * (f_rv[2] + q_fn(fv0[2] - f_rv[2])),
    )
    two_dt = _TWO * dt
    th_dt = _THREE_HALVES * dt
    uh_n = (uh[0] + two_dt * puh_m[0], uh[1] + two_dt * puh_m[1], uh[2] + two_dt * puh_m[2])
    vh_n = (vh[0] + th_dt * pvh_m[0], vh[1] + th_dt * pvh_m[1], vh[2] + th_dt * pvh_m[2])
    fu1, fv1 = b2_shadow_force_uv(uh_n, vh_n, dt)
    puh_n = (
        puh_m[0] + half_dt * (f_ru[0] + q_fn(fu1[0] - f_ru[0])),
        puh_m[1] + half_dt * (f_ru[1] + q_fn(fu1[1] - f_ru[1])),
        puh_m[2] + half_dt * (f_ru[2] + q_fn(fu1[2] - f_ru[2])),
    )
    pvh_n = (
        pvh_m[0] + half_dt * (f_rv[0] + q_fn(fv1[0] - f_rv[0])),
        pvh_m[1] + half_dt * (f_rv[1] + q_fn(fv1[1] - f_rv[1])),
        pvh_m[2] + half_dt * (f_rv[2] + q_fn(fv1[2] - f_rv[2])),
    )
    return (rc_n, pt, uh_n, vh_n, puh_n, pvh_n, z_ref)


def extract_qp_battle2(state: Any) -> Tuple[List[Any], List[Any]]:
    """Extract physical Cartesian `([q1x..q3z], [p1x..p3z])` conserving all 7 invariants `0/1`."""
    rc, pt, uh, vh, puh, pvh, z_ref = state
    u_p, v_p, pu_p, pv_p = project_battle2_manifold(z_ref, (uh, vh, puh, pvh))
    qq, pp = _jacobi_to_cartesian(rc, pt, u_p, v_p, pu_p, pv_p)
    return (list(qq), list(pp))


def invariants_battle2(state: Any) -> Dict[str, Any]:
    """Evaluate all 7 exact physical invariants `(H, P, L)` on `extract_qp_battle2(state)`."""
    q, p = extract_qp_battle2(state)
    return battle2_invariants_q(q, p)


# ============================================================================
# 4. Battle 3: 2D Eccentric Newtonian Kepler 1/r (Levi-Civita Conic Feistel)
# Single Uniform-Time Step (Zero Sub-Stepping), Global Injectivity, Zero Precession
# ============================================================================


def spinor_from_cs(
    z_ref: Tuple[Any, Any, Any, Any], c_val: Any, s_val: Any, kappa: Any
) -> Tuple[Any, Any, Any, Any]:
    """Pure-field linear symplectic spinor flow on `c^2 + 4*kappa*s^2 == 1`."""
    u0, v0, pu0, pv0 = z_ref
    eight = Fraction(8, 1)
    u = c_val * u0 + _HALF * s_val * pu0
    pu = -eight * kappa * s_val * u0 + c_val * pu0
    v = c_val * v0 + _HALF * s_val * pv0
    pv = -eight * kappa * s_val * v0 + c_val * pv0
    return (u, v, pu, pv)


def project_battle3_conic_kernel(
    c_hat: Any, s_hat: Any, kappa: Any, sgn: Any = _ONE
) -> Tuple[Any, Any]:
    """Pure-field antipodal Vieta chord onto `c^2 + 4*kappa*s^2 == 1`."""
    four_k = Fraction(4, 1) * kappa
    dc = c_hat + sgn
    ds = s_hat
    den = dc * dc + four_k * ds * ds
    if den == _ZERO:
        return (sgn, _ZERO)
    lam = (_TWO * sgn * dc) / den
    return (-sgn + lam * dc, lam * ds)


def project_battle3_conic(c_hat: Fraction, s_hat: Fraction, kappa: Fraction) -> Tuple[Fraction, Fraction]:
    """Sign-switched antipodal Vieta projection onto `c^2 + 4*kappa*s^2 == 1`."""
    sgn = _NEG_ONE if c_hat < _ZERO else _ONE
    return project_battle3_conic_kernel(c_hat, s_hat, kappa, sgn)


def _sign_nonneg(val: Fraction) -> Fraction:
    return _NEG_ONE if val < _ZERO else _ONE


def cs_rhs_kernel(
    z_ref: Tuple[Any, Any, Any, Any],
    c_val: Any,
    s_val: Any,
    kappa: Any,
    dt: Any,
    sgn: Any = _ONE,
) -> Tuple[Any, Any]:
    """Pure-field physical-time Kepler `(dc/dt, ds/dt)` with even-in-`dt` perihelion regularization."""
    c_on, s_on = project_battle3_conic_kernel(c_val, s_val, kappa, sgn)
    u, v, _, _ = spinor_from_cs(z_ref, c_on, s_on, kappa)
    r_val = u * u + v * v
    reg_den = Fraction(64, 1) * r_val * r_val + dt * dt
    dc_dt = -(Fraction(128, 1) * kappa * r_val * s_on) / reg_den
    ds_dt = (Fraction(32, 1) * r_val * c_on) / reg_den
    return (dc_dt, ds_dt)


def build_state_from_spinor_battle3(
    u: Any, v: Any, pu: Any, pv: Any, t_phys: Any = _ZERO
) -> Any:
    """Construct minimal globally-injective Battle 3 state `(t_phys, ca, sa, cb, sb, z_ref)`."""
    return (t_phys, _ONE, _ZERO, _ONE, _ZERO, (u, v, pu, pv))


def initial_state_battle3(dt: Optional[Fraction] = None) -> Any:
    """Return canonical Battle 3 aphelion initial state."""
    return build_state_from_spinor_battle3(
        BATTLE3_SPINOR_U0[0],
        BATTLE3_SPINOR_U0[1],
        BATTLE3_SPINOR_PU0[0],
        BATTLE3_SPINOR_PU0[1],
        _ZERO,
    )


def step_battle3(
    state: Any,
    dt: Any,
    q_fn: Callable[[Any], Any] = quantize_sig29,
    sgn_fn: Callable[[Any], Any] = _sign_nonneg,
) -> Any:
    """Single-step symmetric uniform-physical-time Feistel midpoint map (zero sub-stepping)."""
    t_phys, ca, sa, cb, sb, z_ref = state
    energy = battle3_hamiltonian_spinor_q(z_ref[0], z_ref[1], z_ref[2], z_ref[3])
    kappa = -_HALF * energy
    dc_ref, ds_ref = cs_rhs_kernel(z_ref, _ONE, _ZERO, kappa, dt, _ONE)

    half_dt = dt * _HALF
    dca1, dsa1 = cs_rhs_kernel(z_ref, cb, sb, kappa, dt, sgn_fn(cb))
    ca_m = ca + half_dt * (dc_ref + q_fn(dca1 - dc_ref))
    sa_m = sa + half_dt * (ds_ref + q_fn(dsa1 - ds_ref))

    dcb, dsb = cs_rhs_kernel(z_ref, ca_m, sa_m, kappa, dt, sgn_fn(ca_m))
    cb_n = cb + dt * (dc_ref + q_fn(dcb - dc_ref))
    sb_n = sb + dt * (ds_ref + q_fn(dsb - ds_ref))

    dca2, dsa2 = cs_rhs_kernel(z_ref, cb_n, sb_n, kappa, dt, sgn_fn(cb_n))
    ca_n = ca_m + half_dt * (dc_ref + q_fn(dca2 - dc_ref))
    sa_n = sa_m + half_dt * (ds_ref + q_fn(dsa2 - ds_ref))
    return (t_phys + dt, ca_n, sa_n, cb_n, sb_n, z_ref)


def extract_phys_time_battle3(state: Any) -> Any:
    """Extract exact uniform physical time `t_phys == n * dt`."""
    return state[0]


def extract_spinor_battle3_kernel(
    state: Any, sgn: Any = _ONE, q_fn: Callable[[Any], Any] = quantize_sig29
) -> Tuple[Any, Any, Any, Any]:
    """Pure-field spinor extraction kernel on `c^2 + 4*kappa*s^2 == 1`."""
    _, ca, sa, cb, sb, z_ref = state
    energy = battle3_hamiltonian_spinor_q(z_ref[0], z_ref[1], z_ref[2], z_ref[3])
    kappa = -_HALF * energy
    c_mid = _ONE + q_fn(_HALF * (ca + cb) - _ONE)
    s_mid = q_fn(_HALF * (sa + sb))
    c_st, s_st = project_battle3_conic_kernel(c_mid, s_mid, kappa, sgn)
    return spinor_from_cs(z_ref, c_st, s_st, kappa)


def extract_spinor_battle3(state: Any) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    """Extract exact Levi-Civita spinor `(u, v, pu, pv)` on the invariant Kepler torus."""
    _, ca, sa, cb, sb, _ = state
    sgn = _NEG_ONE if (ca + cb) < _ZERO else _ONE
    return extract_spinor_battle3_kernel(state, sgn, quantize_sig29)


def extract_qp_battle3(state: Any) -> Tuple[List[Any], List[Any]]:
    """Project Battle 3 state -> physical Cartesian `([x, y], [px, py])` over `Q` without `sqrt`."""
    u, v, pu, pv = extract_spinor_battle3(state)
    q_cart, p_cart, _ = battle3_spinor_to_cartesian_q(u, v, pu, pv)
    return ([q_cart[0], q_cart[1]], [p_cart[0], p_cart[1]])


def invariants_battle3(state: Any) -> Dict[str, Any]:
    """Evaluate exact Kepler invariants `{'H', 'Lz', 'LRL'}` on `extract_spinor_battle3(state)`."""
    u, v, pu, pv = extract_spinor_battle3(state)
    return battle3_invariants_spinor_q(u, v, pu, pv)


# ============================================================================
# PART II: ZERO-HANDICAP IEEE-754 FLOAT64 CONTINUUM BASELINES
# Strictly isolated from Part I's pure-Fraction Q-Engine call graph.
# Implements Forward Euler O(dt), Symplectic Velocity Verlet O(dt^2),
# Classical RK4 O(dt^4), and fine-grid reference RK4(dt_ref = 1/1024).
# ============================================================================


def battle1_hamiltonian_float(q_f: Sequence[float], p_f: Sequence[float]) -> float:
    """Continuum float64 Hamiltonian H_1(q, p) for Pullen-Edmonds oscillator."""
    q1, q2 = q_f[0], q_f[1]
    p1, p2 = p_f[0], p_f[1]
    alpha_f = float(BATTLE1_ALPHA)
    q1_sq = q1 * q1
    q2_sq = q2 * q2
    return 0.5 * (p1 * p1 + p2 * p2) + 0.5 * (q1_sq + q2_sq) + 0.5 * alpha_f * q1_sq * q2_sq


def battle1_force_float(q_f: Sequence[float]) -> List[float]:
    """Continuum float64 analytical force F_1(q) = -grad V_1(q)."""
    q1, q2 = q_f[0], q_f[1]
    alpha_f = float(BATTLE1_ALPHA)
    return [
        -q1 * (1.0 + alpha_f * q2 * q2),
        -q2 * (1.0 + alpha_f * q1 * q1),
    ]


def battle1_rhs_float(q_f: Sequence[float], p_f: Sequence[float]) -> Tuple[List[float], List[float]]:
    """Continuum float64 ODE RHS (dq/dt, dp/dt) for Battle 1."""
    return ([p_f[0], p_f[1]], battle1_force_float(q_f))


def battle2_hamiltonian_float(q_f: Sequence[float], p_f: Sequence[float]) -> float:
    """Continuum float64 Hamiltonian H_2(q, p) for 3D 3-Body SO(3) Quadrance Cluster."""
    ks = float(BATTLE2_K_SPRING)
    gc = float(BATTLE2_G_COUPLE)
    d12x, d12y, d12z = q_f[0] - q_f[3], q_f[1] - q_f[4], q_f[2] - q_f[5]
    d23x, d23y, d23z = q_f[3] - q_f[6], q_f[4] - q_f[7], q_f[5] - q_f[8]
    d31x, d31y, d31z = q_f[6] - q_f[0], q_f[7] - q_f[1], q_f[8] - q_f[2]
    q12 = d12x * d12x + d12y * d12y + d12z * d12z
    q23 = d23x * d23x + d23y * d23y + d23z * d23z
    q31 = d31x * d31x + d31y * d31y + d31z * d31z
    ke = 0.5 * sum(pc * pc for pc in p_f)
    arch = 2.0 * (q12 * q23 + q23 * q31 + q31 * q12) - (q12 * q12 + q23 * q23 + q31 * q31)
    pot = 0.5 * ks * (q12 + q23 + q31) + 0.25 * gc * arch
    return ke + pot


def battle2_force_float(q_f: Sequence[float]) -> List[float]:
    """Continuum float64 analytical forces F_i(q) = -grad_{q_i} V_2(q)."""
    ks = float(BATTLE2_K_SPRING)
    gc = float(BATTLE2_G_COUPLE)
    r12x, r12y, r12z = q_f[0] - q_f[3], q_f[1] - q_f[4], q_f[2] - q_f[5]
    r23x, r23y, r23z = q_f[3] - q_f[6], q_f[4] - q_f[7], q_f[5] - q_f[8]
    r31x, r31y, r31z = q_f[6] - q_f[0], q_f[7] - q_f[1], q_f[8] - q_f[2]
    q12 = r12x * r12x + r12y * r12y + r12z * r12z
    q23 = r23x * r23x + r23y * r23y + r23z * r23z
    q31 = r31x * r31x + r31y * r31y + r31z * r31z
    c12 = ks + gc * (q23 + q31 - q12)
    c23 = ks + gc * (q31 + q12 - q23)
    c31 = ks + gc * (q12 + q23 - q31)
    return [
        -c12 * r12x + c31 * r31x,
        -c12 * r12y + c31 * r31y,
        -c12 * r12z + c31 * r31z,
        c12 * r12x - c23 * r23x,
        c12 * r12y - c23 * r23y,
        c12 * r12z - c23 * r23z,
        -c31 * r31x + c23 * r23x,
        -c31 * r31y + c23 * r23y,
        -c31 * r31z + c23 * r23z,
    ]


def battle2_rhs_float(q_f: Sequence[float], p_f: Sequence[float]) -> Tuple[List[float], List[float]]:
    """Continuum float64 ODE RHS (dq/dt, dp/dt) for Battle 2."""
    return (list(p_f), battle2_force_float(q_f))


def battle2_linear_momentum_float(q_f: Sequence[float], p_f: Sequence[float]) -> Tuple[float, float, float]:
    """Continuum float64 total linear momentum P."""
    return (p_f[0] + p_f[3] + p_f[6], p_f[1] + p_f[4] + p_f[7], p_f[2] + p_f[5] + p_f[8])


def battle2_angular_momentum_float(q_f: Sequence[float], p_f: Sequence[float]) -> Tuple[float, float, float]:
    """Continuum float64 total angular momentum L = sum q_i x p_i."""
    lx = (
        (q_f[1] * p_f[2] - q_f[2] * p_f[1])
        + (q_f[4] * p_f[5] - q_f[5] * p_f[4])
        + (q_f[7] * p_f[8] - q_f[8] * p_f[7])
    )
    ly = (
        (q_f[2] * p_f[0] - q_f[0] * p_f[2])
        + (q_f[5] * p_f[3] - q_f[3] * p_f[5])
        + (q_f[8] * p_f[6] - q_f[6] * p_f[8])
    )
    lz = (
        (q_f[0] * p_f[1] - q_f[1] * p_f[0])
        + (q_f[3] * p_f[4] - q_f[4] * p_f[3])
        + (q_f[6] * p_f[7] - q_f[7] * p_f[6])
    )
    return (lx, ly, lz)


def battle3_hamiltonian_float(q_f: Sequence[float], p_f: Sequence[float]) -> float:
    """Continuum float64 Kepler Hamiltonian H_3(q, p) (math.sqrt ONLY in float baselines)."""
    x, y = q_f[0], q_f[1]
    px, py = p_f[0], p_f[1]
    r = math.sqrt(x * x + y * y)
    return 0.5 * (px * px + py * py) - float(BATTLE3_MU) / r


def battle3_force_float(q_f: Sequence[float]) -> List[float]:
    """Continuum float64 Newtonian gravitational force F_3(q) = -mu * q / r^3."""
    x, y = q_f[0], q_f[1]
    r2 = x * x + y * y
    r = math.sqrt(r2)
    r3 = r2 * r
    mu_f = float(BATTLE3_MU)
    return [-mu_f * x / r3, -mu_f * y / r3]


def battle3_rhs_float(q_f: Sequence[float], p_f: Sequence[float]) -> Tuple[List[float], List[float]]:
    """Continuum float64 ODE RHS (dq/dt, dp/dt) for Battle 3."""
    return ([p_f[0], p_f[1]], battle3_force_float(q_f))


def battle3_angular_momentum_float(q_f: Sequence[float], p_f: Sequence[float]) -> Tuple[float]:
    """Continuum float64 scalar planar angular momentum (Lz,)."""
    return (q_f[0] * p_f[1] - q_f[1] * p_f[0],)


def baseline_forward_euler_step(
    q_f: Sequence[float],
    p_f: Sequence[float],
    dt_f: float,
    rhs_fn: Callable[[Sequence[float], Sequence[float]], Tuple[List[float], List[float]]],
) -> Tuple[List[float], List[float]]:
    """Standard 1st-order explicit Forward Euler step."""
    dq, dp = rhs_fn(q_f, p_f)
    return (
        [qi + dt_f * dqi for qi, dqi in zip(q_f, dq)],
        [pi + dt_f * dpi for pi, dpi in zip(p_f, dp)],
    )


def baseline_velocity_verlet_step(
    q_f: Sequence[float],
    p_f: Sequence[float],
    dt_f: float,
    force_fn: Callable[[Sequence[float]], List[float]],
    inv_masses: Sequence[float],
) -> Tuple[List[float], List[float]]:
    """Standard 2nd-order symplectic Velocity Verlet (kick-drift-kick) step."""
    half_dt = 0.5 * dt_f
    f0 = force_fn(q_f)
    p_half = [pi + half_dt * fi for pi, fi in zip(p_f, f0)]
    q_next = [qi + dt_f * inv_m * phi for qi, inv_m, phi in zip(q_f, inv_masses, p_half)]
    f1 = force_fn(q_next)
    p_next = [phi + half_dt * fi for phi, fi in zip(p_half, f1)]
    return (q_next, p_next)


def baseline_rk4_step(
    q_f: Sequence[float],
    p_f: Sequence[float],
    dt_f: float,
    rhs_fn: Callable[[Sequence[float], Sequence[float]], Tuple[List[float], List[float]]],
) -> Tuple[List[float], List[float]]:
    """Classical 4-stage 4th-order Runge-Kutta (RK4) step."""
    half_dt = 0.5 * dt_f
    sixth_dt = dt_f / 6.0
    k1_q, k1_p = rhs_fn(q_f, p_f)
    q2 = [qi + half_dt * kq for qi, kq in zip(q_f, k1_q)]
    p2 = [pi + half_dt * kp for pi, kp in zip(p_f, k1_p)]
    k2_q, k2_p = rhs_fn(q2, p2)
    q3 = [qi + half_dt * kq for qi, kq in zip(q_f, k2_q)]
    p3 = [pi + half_dt * kp for pi, kp in zip(p_f, k2_p)]
    k3_q, k3_p = rhs_fn(q3, p3)
    q4 = [qi + dt_f * kq for qi, kq in zip(q_f, k3_q)]
    p4 = [pi + dt_f * kp for pi, kp in zip(p_f, k3_p)]
    k4_q, k4_p = rhs_fn(q4, p4)
    q_next = [
        qi + sixth_dt * (a + 2.0 * b + 2.0 * c + d)
        for qi, a, b, c, d in zip(q_f, k1_q, k2_q, k3_q, k4_q)
    ]
    p_next = [
        pi + sixth_dt * (a + 2.0 * b + 2.0 * c + d)
        for pi, a, b, c, d in zip(p_f, k1_p, k2_p, k3_p, k4_p)
    ]
    return (q_next, p_next)


def get_battle_baseline_steppers() -> Dict[
    str, Callable[[int, Sequence[float], Sequence[float], float], Tuple[List[float], List[float]]]
]:
    """Expose 1-step baseline hooks for independent Gate 4 auditing."""
    rhs_map = {1: battle1_rhs_float, 2: battle2_rhs_float, 3: battle3_rhs_float}
    force_map = {1: battle1_force_float, 2: battle2_force_float, 3: battle3_force_float}
    invm_map = {1: [1.0, 1.0], 2: [1.0] * 9, 3: [1.0, 1.0]}
    return {
        "euler": lambda bid, q_f, p_f, dt_f: baseline_forward_euler_step(q_f, p_f, dt_f, rhs_map[bid]),
        "verlet": lambda bid, q_f, p_f, dt_f: baseline_velocity_verlet_step(
            q_f, p_f, dt_f, force_map[bid], invm_map[bid]
        ),
        "rk4": lambda bid, q_f, p_f, dt_f: baseline_rk4_step(q_f, p_f, dt_f, rhs_map[bid]),
    }


def _sim_wrap_pi(ang: float) -> float:
    twopi = 2.0 * math.pi
    w = (ang + math.pi) % twopi - math.pi
    return w + twopi if w <= -math.pi else w


def _sim_lrl_drift_series(qs: Sequence[Sequence[float]], ps: Sequence[Sequence[float]], mu_f: float = 1.0) -> List[float]:
    """Compute cumulative unwrapped osculating LRL perihelion angle drift (degrees)."""
    drifts_deg: List[float] = []
    th0 = 0.0
    prev_raw = 0.0
    unwrapped = 0.0
    for idx, (q_f, p_f) in enumerate(zip(qs, ps)):
        x, y = q_f[0], q_f[1]
        px, py = p_f[0], p_f[1]
        r = math.hypot(x, y)
        if r <= 1.0e-15 or not math.isfinite(r):
            break
        lz = x * py - y * px
        ax = py * lz - mu_f * (x / r)
        ay = -px * lz - mu_f * (y / r)
        ecc = math.hypot(ax, ay) / mu_f
        raw = math.atan2(ay, ax) if ecc > 1.0e-10 else 0.0
        if idx == 0:
            th0 = raw
            prev_raw = raw
            unwrapped = raw
        else:
            unwrapped += _sim_wrap_pi(raw - prev_raw)
            prev_raw = raw
        drifts_deg.append(abs(math.degrees(unwrapped - th0)))
    return drifts_deg


def _sim_count_perihelia(qs: Sequence[Sequence[float]], ps: Sequence[Sequence[float]]) -> int:
    s_vals = [q_f[0] * p_f[0] + q_f[1] * p_f[1] for q_f, p_f in zip(qs, ps)]
    count = 0
    for idx in range(len(s_vals) - 1):
        if s_vals[idx] < 0.0 and s_vals[idx + 1] >= 0.0:
            count += 1
    return count


def _sim_rk4_reference_at_horizons(
    q0: Sequence[float],
    p0: Sequence[float],
    rhs_fn: Callable[[Sequence[float], Sequence[float]], Tuple[List[float], List[float]]],
    horizons: Sequence[float],
    dt_ref: float = 1.0 / 1024.0,
) -> Dict[float, Tuple[List[float], List[float]]]:
    q_c = list(q0)
    p_c = list(p0)
    t_c = 0.0
    out: Dict[float, Tuple[List[float], List[float]]] = {}
    for h_t in sorted(horizons):
        n_sub = max(1, int(round((h_t - t_c) / dt_ref)))
        dt_sub = (h_t - t_c) / float(n_sub)
        for _ in range(n_sub):
            q_c, p_c = baseline_rk4_step(q_c, p_c, dt_sub, rhs_fn)
        t_c = h_t
        out[h_t] = (list(q_c), list(p_c))
    return out


def _sim_run_float_baseline(
    battle_id: int,
    method: str,
    dt_f: float,
    n_steps: int,
    ref_horizons: Dict[float, Tuple[List[float], List[float]]],
) -> Dict[str, Any]:
    q0_map = {
        1: [float(x) for x in BATTLE1_Q0],
        2: [float(x) for x in BATTLE2_Q0],
        3: [float(x) for x in BATTLE3_Q0],
    }
    p0_map = {
        1: [float(x) for x in BATTLE1_P0],
        2: [float(x) for x in BATTLE2_P0],
        3: [float(x) for x in BATTLE3_P0],
    }
    ham_map = {1: battle1_hamiltonian_float, 2: battle2_hamiltonian_float, 3: battle3_hamiltonian_float}
    rhs_map = {1: battle1_rhs_float, 2: battle2_rhs_float, 3: battle3_rhs_float}
    force_map = {1: battle1_force_float, 2: battle2_force_float, 3: battle3_force_float}
    invm_map = {1: [1.0, 1.0], 2: [1.0] * 9, 3: [1.0, 1.0]}

    q_c = list(q0_map[battle_id])
    p_c = list(p0_map[battle_id])
    ham_fn = ham_map[battle_id]
    rhs_fn = rhs_map[battle_id]
    force_fn = force_map[battle_id]
    invm = invm_map[battle_id]

    h0 = ham_fn(q_c, p_c)
    p_mom0 = battle2_linear_momentum_float(q_c, p_c) if battle_id == 2 else None
    l_mom0 = (
        battle2_angular_momentum_float(q_c, p_c)
        if battle_id == 2
        else (battle3_angular_momentum_float(q_c, p_c) if battle_id == 3 else None)
    )

    qs: List[List[float]] = [list(q_c)]
    ps: List[List[float]] = [list(p_c)]
    times: List[float] = [0.0]
    dh_series: List[float] = [0.0]
    dp_series: List[float] = [0.0]
    dl_series: List[float] = [0.0]

    blew_up = False
    blowup_step: Optional[int] = None
    escaped = False
    escape_step: Optional[int] = None
    steps_done = 0

    for step_idx in range(1, n_steps + 1):
        try:
            if method == "euler":
                q_n, p_n = baseline_forward_euler_step(q_c, p_c, dt_f, rhs_fn)
            elif method == "verlet":
                q_n, p_n = baseline_velocity_verlet_step(q_c, p_c, dt_f, force_fn, invm)
            else:
                q_n, p_n = baseline_rk4_step(q_c, p_c, dt_f, rhs_fn)
        except (OverflowError, ZeroDivisionError, ValueError):
            blew_up = True
            blowup_step = step_idx
            break

        if not all(math.isfinite(v) and abs(v) <= 1.0e100 for v in q_n + p_n):
            blew_up = True
            blowup_step = step_idx
            break

        try:
            hn = ham_fn(q_n, p_n)
        except (OverflowError, ZeroDivisionError, ValueError):
            blew_up = True
            blowup_step = step_idx
            break
        if not math.isfinite(hn):
            blew_up = True
            blowup_step = step_idx
            break

        q_c, p_c = q_n, p_n
        steps_done = step_idx
        t_cur = step_idx * dt_f
        qs.append(list(q_c))
        ps.append(list(p_c))
        times.append(t_cur)
        dh_series.append(abs(hn - h0))

        if p_mom0 is not None:
            pm = battle2_linear_momentum_float(q_c, p_c)
            dp_series.append(math.sqrt(sum((a - b) ** 2 for a, b in zip(pm, p_mom0))))
        if l_mom0 is not None:
            lm = (
                battle2_angular_momentum_float(q_c, p_c)
                if battle_id == 2
                else battle3_angular_momentum_float(q_c, p_c)
            )
            dl_series.append(math.sqrt(sum((a - b) ** 2 for a, b in zip(lm, l_mom0))))

        if battle_id == 3 and not escaped:
            r_cur = math.hypot(q_c[0], q_c[1])
            rdot = q_c[0] * p_c[0] + q_c[1] * p_c[1]
            if hn >= 0.0 and r_cur > 10.0 and rdot > 0.0:
                escaped = True
                escape_step = step_idx

    horizon_errs: Dict[str, Optional[float]] = {}
    for h_t, (q_ref, p_ref) in ref_horizons.items():
        idx_target = int(round(h_t / dt_f))
        if idx_target < len(qs) and abs(times[idx_target] - h_t) <= 1.0e-12:
            err = math.sqrt(
                sum((a - b) ** 2 for a, b in zip(qs[idx_target], q_ref))
                + sum((a - b) ** 2 for a, b in zip(ps[idx_target], p_ref))
            )
            horizon_errs[f"T={h_t:g}"] = err
        else:
            horizon_errs[f"T={h_t:g}"] = None

    lrl_drifts = _sim_lrl_drift_series(qs, ps) if battle_id == 3 else []
    peri_count = _sim_count_perihelia(qs, ps) if battle_id == 3 else 0

    return {
        "method": method,
        "dt_float": dt_f,
        "steps_completed": steps_done,
        "blew_up": blew_up,
        "blowup_step": blowup_step,
        "escaped": escaped,
        "escape_step": escape_step,
        "max_dh": max(dh_series, default=0.0),
        "max_dp": max(dp_series, default=0.0),
        "max_dl": max(dl_series, default=0.0),
        "horizon_errors": horizon_errs,
        "max_lrl_drift_deg": max(lrl_drifts, default=0.0),
        "perihelia": peri_count,
        "qs": qs,
        "ps": ps,
        "dh_series": dh_series,
        "dp_series": dp_series,
        "dl_series": dl_series,
        "lrl_drifts_deg": lrl_drifts,
    }


# ============================================================================
# PART III: END-TO-END 1,000-STEP EXACT RATIONAL SUITE & ASSERTION HARNESS
# Asserts strict Fraction closure, <= 10,000-bit Ultrafinitist ceiling,
# exact 0/1 conservation of ALL physical invariants, and exact 0/1
# discrete time-reversibility at every single step n = 1..1000.
# ============================================================================


def _sim_inspect_fraction_bits(obj: Any) -> Tuple[int, int, bool]:
    """Recursively verify strict Fraction leaf types and return (max_bits, leaf_count, strict_ok)."""
    max_b = 0
    count = 0
    strict_ok = True

    def _walk(node: Any) -> None:
        nonlocal max_b, count, strict_ok
        if node is None or isinstance(node, (bool, str)):
            return
        if type(node) is Fraction:
            count += 1
            nb = node.numerator.bit_length()
            db = node.denominator.bit_length()
            mb = nb if nb >= db else db
            if mb > max_b:
                max_b = mb
            return
        if isinstance(node, (tuple, list)):
            for item in node:
                _walk(item)
            return
        if isinstance(node, dict):
            for val in node.values():
                _walk(val)
            return
        strict_ok = False

    _walk(obj)
    return (max_b, count, strict_ok)


def _sim_run_exact_q_battle(
    battle_id: int,
    dt_q: Fraction,
    n_steps: int,
    ref_horizons: Dict[float, Tuple[List[float], List[float]]],
) -> Dict[str, Any]:
    t0_wall = time.perf_counter()
    if battle_id == 1:
        st = initial_state_battle1(dt_q)
        q0_q, p0_q = extract_qp_battle1(st)
        aux0_q: Any = ()
        inv0 = invariants_battle1(st)
    elif battle_id == 2:
        st = initial_state_battle2(dt_q)
        q0_q, p0_q = extract_qp_battle2(st)
        aux0_q = ()
        inv0 = invariants_battle2(st)
    else:
        st = initial_state_battle3(dt_q)
        q0_q, p0_q = extract_qp_battle3(st)
        aux0_q = extract_spinor_battle3(st)
        inv0 = invariants_battle3(st)

    b0_comb, _, ok0 = _sim_inspect_fraction_bits((st, q0_q, p0_q, aux0_q))
    b0_state, _, _ = _sim_inspect_fraction_bits(st)
    assert ok0, f"Initial state non-Fraction leaf in Battle {battle_id}"

    bit_series: List[int] = [b0_comb]
    state_bit_series: List[int] = [b0_state]
    qs_flt: List[List[float]] = [[float(x) for x in q0_q]]
    ps_flt: List[List[float]] = [[float(x) for x in p0_q]]
    times_flt: List[float] = [0.0]

    prev_q_q, prev_p_q = q0_q, p0_q

    for step_idx in range(1, n_steps + 1):
        prev_st = st
        if battle_id == 1:
            st = step_battle1(prev_st, dt_q)
            qn_q, pn_q = extract_qp_battle1(st)
            auxn_q: Any = ()
            hn_q = battle1_hamiltonian_q(qn_q, pn_q)
            assert hn_q - inv0["H"] == Fraction(0, 1), (
                f"Battle 1 H drift at step {step_idx}: {hn_q - inv0['H']}"
            )
            if step_idx in (1, n_steps):
                assert invariants_battle1(st)["H"] - inv0["H"] == Fraction(0, 1)
            back_st = step_battle1(st, -dt_q)
            assert back_st == prev_st, f"Battle 1 reversibility failed at step {step_idx}"
            if step_idx in (1, 10, 100, 500, n_steps):
                assert extract_qp_battle1(back_st) == (prev_q_q, prev_p_q)
            tn_flt = float(step_idx * dt_q)

        elif battle_id == 2:
            st = step_battle2(prev_st, dt_q)
            qn_q, pn_q = extract_qp_battle2(st)
            auxn_q = ()
            inv_n = battle2_invariants_q(qn_q, pn_q)
            assert inv_n["H"] - inv0["H"] == Fraction(0, 1), (
                f"Battle 2 H drift at step {step_idx}"
            )
            assert all(
                a - b == Fraction(0, 1) for a, b in zip(inv_n["P"], inv0["P"])
            ), f"Battle 2 P drift at step {step_idx}"
            assert all(
                a - b == Fraction(0, 1) for a, b in zip(inv_n["L"], inv0["L"])
            ), f"Battle 2 L drift at step {step_idx}"
            if step_idx in (1, n_steps):
                assert invariants_battle2(st) == inv0
            back_st = step_battle2(st, -dt_q)
            assert back_st == prev_st, f"Battle 2 reversibility failed at step {step_idx}"
            if step_idx in (1, 10, 100, 500, n_steps):
                assert extract_qp_battle2(back_st) == (prev_q_q, prev_p_q)
            tn_flt = float(step_idx * dt_q)

        else:
            st = step_battle3(prev_st, dt_q)
            qn_q, pn_q = extract_qp_battle3(st)
            un_q, vn_q, pun_q, pvn_q = extract_spinor_battle3(st)
            auxn_q = (un_q, vn_q, pun_q, pvn_q)
            tn_q = extract_phys_time_battle3(st)
            assert tn_q - step_idx * dt_q == Fraction(0, 1), (
                f"Battle 3 non-uniform time at step {step_idx}"
            )
            rn_q = un_q * un_q + vn_q * vn_q
            assert qn_q[0] * qn_q[0] + qn_q[1] * qn_q[1] == rn_q * rn_q, (
                f"Battle 3 radius not in Q at step {step_idx}"
            )
            inv_n = battle3_invariants_spinor_q(un_q, vn_q, pun_q, pvn_q)
            assert inv_n["H"] - inv0["H"] == Fraction(0, 1), (
                f"Battle 3 H drift at step {step_idx}"
            )
            assert inv_n["Lz"] - inv0["Lz"] == Fraction(0, 1), (
                f"Battle 3 Lz drift at step {step_idx}"
            )
            assert (
                inv_n["LRL"][0] - inv0["LRL"][0] == Fraction(0, 1)
                and inv_n["LRL"][1] - inv0["LRL"][1] == Fraction(0, 1)
            ), f"Battle 3 LRL drift at step {step_idx}"
            if step_idx in (1, n_steps):
                assert invariants_battle3(st) == inv0
            back_st = step_battle3(st, -dt_q)
            assert back_st == prev_st, f"Battle 3 reversibility failed at step {step_idx}"
            if step_idx in (1, 10, 100, 500, n_steps):
                assert extract_qp_battle3(back_st) == (prev_q_q, prev_p_q)
            tn_flt = float(tn_q)

        b_comb, _, strict_ok = _sim_inspect_fraction_bits((st, qn_q, pn_q, auxn_q))
        b_st, _, _ = _sim_inspect_fraction_bits(st)
        assert strict_ok, f"Non-Fraction leaf at Battle {battle_id} step {step_idx}"
        assert b_comb <= 10000, (
            f"Gate 2 bit-length ceiling exceeded ({b_comb} > 10000) at Battle {battle_id} step {step_idx}"
        )

        bit_series.append(b_comb)
        state_bit_series.append(b_st)
        qs_flt.append([float(x) for x in qn_q])
        ps_flt.append([float(x) for x in pn_q])
        times_flt.append(tn_flt)
        prev_q_q, prev_p_q = qn_q, pn_q

    wall_sec = time.perf_counter() - t0_wall
    dt_f = float(dt_q)
    horizon_errs: Dict[str, Optional[float]] = {}
    for h_t, (q_ref, p_ref) in ref_horizons.items():
        idx_target = int(round(h_t / dt_f))
        if idx_target < len(qs_flt) and abs(times_flt[idx_target] - h_t) <= 1.0e-12:
            err = math.sqrt(
                sum((a - b) ** 2 for a, b in zip(qs_flt[idx_target], q_ref))
                + sum((a - b) ** 2 for a, b in zip(ps_flt[idx_target], p_ref))
            )
            horizon_errs[f"T={h_t:g}"] = err
        else:
            horizon_errs[f"T={h_t:g}"] = None

    lrl_drifts = _sim_lrl_drift_series(qs_flt, ps_flt) if battle_id == 3 else []
    peri_count = _sim_count_perihelia(qs_flt, ps_flt) if battle_id == 3 else 0

    return {
        "dt": f"{dt_q.numerator}/{dt_q.denominator}",
        "dt_float": dt_f,
        "steps_completed": n_steps,
        "peak_combined_bits": max(bit_series),
        "peak_state_bits": max(state_bit_series),
        "final_combined_bits": bit_series[-1],
        "exact_invariants_0_over_1": True,
        "exact_reversibility_0_over_1": True,
        "horizon_errors": horizon_errs,
        "max_lrl_drift_deg": max(lrl_drifts, default=0.0),
        "perihelia": peri_count,
        "wall_sec": wall_sec,
        "bit_series": bit_series,
        "qs": qs_flt,
        "ps": ps_flt,
        "lrl_drifts_deg": lrl_drifts,
    }


def _sim_verify_planar_submanifolds() -> Dict[str, Any]:
    """Assert exact 0/1 7-invariant conservation & reversibility on planar 3-body submanifolds."""
    cases = {
        "equilateral_breathing_from_rest": (
            [
                Fraction(1, 1), Fraction(0, 1), Fraction(0, 1),
                Fraction(0, 1), Fraction(1, 1), Fraction(0, 1),
                Fraction(0, 1), Fraction(0, 1), Fraction(1, 1),
            ],
            [Fraction(0, 1)] * 9,
        ),
        "planar_z0_motion": (
            [
                Fraction(1, 1), Fraction(0, 1), Fraction(0, 1),
                Fraction(-1, 2), Fraction(4, 5), Fraction(0, 1),
                Fraction(-1, 2), Fraction(-3, 4), Fraction(0, 1),
            ],
            [
                Fraction(0, 1), Fraction(1, 2), Fraction(0, 1),
                Fraction(-2, 5), Fraction(-1, 4), Fraction(0, 1),
                Fraction(2, 5), Fraction(-1, 4), Fraction(0, 1),
            ],
        ),
    }
    dt_p = Fraction(1, 4)
    out: Dict[str, Any] = {}
    for cname, (qe, pe) in cases.items():
        st = build_state_from_qp_battle2(qe, pe)
        inv0 = invariants_battle2(st)
        peak_b = 0
        for _ in range(100):
            prev = st
            st = step_battle2(st, dt_p)
            qp = extract_qp_battle2(st)
            b_cur, _, ok = _sim_inspect_fraction_bits((st, qp))
            assert ok and b_cur <= 10000
            if b_cur > peak_b:
                peak_b = b_cur
            inv_c = battle2_invariants_q(qp[0], qp[1])
            assert inv_c == inv0, f"Planar case {cname} invariant drift!"
            assert step_battle2(st, -dt_p) == prev, f"Planar case {cname} reversibility failed!"
        out[cname] = {"steps_verified": 100, "peak_bits": peak_b, "exact_0_over_1": True}
    return out


# ============================================================================
# PART IV: ZERO-DEPENDENCY PURE-STDLIB SVG + ZLIB PNG DIAGNOSTIC PLOT ENGINE
# Uses ONLY Python built-in `.to_bytes(4, "big")` + `zlib` (zero `struct`,
# zero `matplotlib`, zero `PIL`) to emit both crisp vector `.svg` and
# truecolor raster `.png` diagnostic plots into `plots/`.
# ============================================================================

_PLOT_FONT_5X7: Dict[str, Tuple[int, ...]] = {
    " ": (0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00),
    "0": (0x0E, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0E),
    "1": (0x04, 0x0C, 0x04, 0x04, 0x04, 0x04, 0x0E),
    "2": (0x0E, 0x11, 0x01, 0x06, 0x08, 0x10, 0x1F),
    "3": (0x1F, 0x02, 0x04, 0x06, 0x01, 0x11, 0x0E),
    "4": (0x02, 0x06, 0x0A, 0x12, 0x1F, 0x02, 0x02),
    "5": (0x1F, 0x10, 0x1E, 0x01, 0x01, 0x11, 0x0E),
    "6": (0x06, 0x08, 0x10, 0x1E, 0x11, 0x11, 0x0E),
    "7": (0x1F, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08),
    "8": (0x0E, 0x11, 0x11, 0x0E, 0x11, 0x11, 0x0E),
    "9": (0x0E, 0x11, 0x11, 0x0F, 0x01, 0x02, 0x0C),
    "A": (0x0E, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11),
    "B": (0x1E, 0x11, 0x11, 0x1E, 0x11, 0x11, 0x1E),
    "C": (0x0E, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0E),
    "D": (0x1C, 0x12, 0x11, 0x11, 0x11, 0x12, 0x1C),
    "E": (0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x1F),
    "F": (0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x10),
    "G": (0x0E, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0F),
    "H": (0x11, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11),
    "I": (0x0E, 0x04, 0x04, 0x04, 0x04, 0x04, 0x0E),
    "J": (0x07, 0x02, 0x02, 0x02, 0x02, 0x12, 0x0C),
    "K": (0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11),
    "L": (0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1F),
    "M": (0x11, 0x1B, 0x15, 0x15, 0x11, 0x11, 0x11),
    "N": (0x11, 0x19, 0x15, 0x13, 0x11, 0x11, 0x11),
    "O": (0x0E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E),
    "P": (0x1E, 0x11, 0x11, 0x1E, 0x10, 0x10, 0x10),
    "Q": (0x0E, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0D),
    "R": (0x1E, 0x11, 0x11, 0x1E, 0x14, 0x12, 0x11),
    "S": (0x0F, 0x10, 0x10, 0x0E, 0x01, 0x01, 0x1E),
    "T": (0x1F, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04),
    "U": (0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E),
    "V": (0x11, 0x11, 0x11, 0x11, 0x11, 0x0A, 0x04),
    "W": (0x11, 0x11, 0x11, 0x15, 0x15, 0x1B, 0x11),
    "X": (0x11, 0x11, 0x0A, 0x04, 0x0A, 0x11, 0x11),
    "Y": (0x11, 0x11, 0x0A, 0x04, 0x04, 0x04, 0x04),
    "Z": (0x1F, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1F),
    ".": (0x00, 0x00, 0x00, 0x00, 0x00, 0x0C, 0x0C),
    ",": (0x00, 0x00, 0x00, 0x00, 0x0C, 0x04, 0x08),
    ":": (0x00, 0x0C, 0x0C, 0x00, 0x0C, 0x0C, 0x00),
    "-": (0x00, 0x00, 0x00, 0x1F, 0x00, 0x00, 0x00),
    "+": (0x00, 0x04, 0x04, 0x1F, 0x04, 0x04, 0x00),
    "=": (0x00, 0x00, 0x1F, 0x00, 0x1F, 0x00, 0x00),
    "/": (0x01, 0x02, 0x02, 0x04, 0x08, 0x08, 0x10),
    "(": (0x02, 0x04, 0x08, 0x08, 0x08, 0x04, 0x02),
    ")": (0x08, 0x04, 0x02, 0x02, 0x02, 0x04, 0x08),
    "[": (0x0E, 0x08, 0x08, 0x08, 0x08, 0x08, 0x0E),
    "]": (0x0E, 0x02, 0x02, 0x02, 0x02, 0x02, 0x0E),
    "<": (0x02, 0x04, 0x08, 0x10, 0x08, 0x04, 0x02),
    ">": (0x08, 0x04, 0x02, 0x01, 0x02, 0x04, 0x08),
    "|": (0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04),
    "^": (0x04, 0x0A, 0x11, 0x00, 0x00, 0x00, 0x00),
    "_": (0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1F),
    "*": (0x00, 0x04, 0x15, 0x0E, 0x15, 0x04, 0x00),
    "~": (0x00, 0x00, 0x08, 0x15, 0x02, 0x00, 0x00),
    "%": (0x19, 0x19, 0x02, 0x04, 0x08, 0x13, 0x13),
    "@": (0x0E, 0x11, 0x17, 0x15, 0x17, 0x10, 0x0E),
    "&": (0x0C, 0x12, 0x14, 0x08, 0x15, 0x12, 0x0D),
    "#": (0x0A, 0x1F, 0x0A, 0x0A, 0x1F, 0x0A, 0x00),
}


def _plot_png_chunk(chunk_type: bytes, data: bytes) -> bytes:
    crc = zlib.crc32(chunk_type + data) & 0xFFFFFFFF
    return len(data).to_bytes(4, "big") + chunk_type + data + crc.to_bytes(4, "big")


class _PlotCanvas:
    """Dual raster RGB PNG + vector SVG scientific plotting canvas (zero external deps)."""

    def __init__(self, width: int, height: int, bg: Tuple[int, int, int] = (250, 251, 253)) -> None:
        self.width = width
        self.height = height
        self.buf = bytearray(bytes(bg) * (width * height))
        self.svg_items: List[str] = [
            f'<rect x="0" y="0" width="{width}" height="{height}" fill="rgb({bg[0]},{bg[1]},{bg[2]})"/>'
        ]

    def _set_px(self, x: int, y: int, color: Tuple[int, int, int]) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            idx = (y * self.width + x) * 3
            self.buf[idx] = color[0]
            self.buf[idx + 1] = color[1]
            self.buf[idx + 2] = color[2]

    def fill_rect(
        self,
        x0: int,
        y0: int,
        w: int,
        h: int,
        fill: Tuple[int, int, int],
        border: Optional[Tuple[int, int, int]] = None,
    ) -> None:
        row_bytes = bytes(fill) * max(0, min(self.width, x0 + w) - max(0, x0))
        xa = max(0, x0)
        xb = min(self.width, x0 + w)
        if xb > xa:
            for yy in range(max(0, y0), min(self.height, y0 + h)):
                start = (yy * self.width + xa) * 3
                self.buf[start : start + (xb - xa) * 3] = row_bytes
        stroke_attr = (
            f' stroke="rgb({border[0]},{border[1]},{border[2]})" stroke-width="1"'
            if border is not None
            else ""
        )
        self.svg_items.append(
            f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" '
            f'fill="rgb({fill[0]},{fill[1]},{fill[2]})"{stroke_attr}/>'
        )
        if border is not None:
            for xx in range(x0, x0 + w):
                self._set_px(xx, y0, border)
                self._set_px(xx, y0 + h - 1, border)
            for yy in range(y0, y0 + h):
                self._set_px(x0, yy, border)
                self._set_px(x0 + w - 1, yy, border)

    def _draw_segment_px(
        self, x0: int, y0: int, x1: int, y1: int, color: Tuple[int, int, int], thickness: int = 2
    ) -> None:
        dx = abs(x1 - x0)
        dy = -abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        rad = max(0, thickness - 1)
        cx, cy = x0, y0
        while True:
            if rad == 0:
                self._set_px(cx, cy, color)
            else:
                for ox in range(-rad, rad + 1):
                    for oy in range(-rad, rad + 1):
                        if ox * ox + oy * oy <= rad * rad + 1:
                            self._set_px(cx + ox, cy + oy, color)
            if cx == x1 and cy == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                cx += sx
            if e2 <= dx:
                err += dx
                cy += sy

    def draw_line(
        self,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        color: Tuple[int, int, int],
        thickness: int = 2,
        dashed: bool = False,
    ) -> None:
        ix0, iy0, ix1, iy1 = int(round(x0)), int(round(y0)), int(round(x1)), int(round(y1))
        if dashed:
            dist = max(1.0, math.hypot(x1 - x0, y1 - y0))
            steps = max(1, int(dist // 8))
            for s_idx in range(0, steps, 2):
                ta = s_idx / float(steps)
                tb = min(1.0, (s_idx + 1) / float(steps))
                self._draw_segment_px(
                    int(round(x0 + ta * (x1 - x0))),
                    int(round(y0 + ta * (y1 - y0))),
                    int(round(x0 + tb * (x1 - x0))),
                    int(round(y0 + tb * (y1 - y0))),
                    color,
                    max(1, thickness - 1),
                )
            dash_attr = ' stroke-dasharray="6,5"'
        else:
            self._draw_segment_px(ix0, iy0, ix1, iy1, color, thickness)
            dash_attr = ""
        self.svg_items.append(
            f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
            f'stroke="rgb({color[0]},{color[1]},{color[2]})" stroke-width="{thickness}"{dash_attr}/>'
        )

    def draw_polyline(
        self,
        pts: Sequence[Tuple[float, float]],
        color: Tuple[int, int, int],
        thickness: int = 2,
        clip_box: Optional[Tuple[int, int, int, int]] = None,
    ) -> None:
        if len(pts) < 2:
            return
        filtered: List[Tuple[float, float]] = []
        for px, py in pts:
            if not (math.isfinite(px) and math.isfinite(py)):
                continue
            if clip_box is not None:
                bx, by, bw, bh = clip_box
                px = max(float(bx), min(float(bx + bw), px))
                py = max(float(by), min(float(by + bh), py))
            filtered.append((px, py))
        if len(filtered) < 2:
            return
        for idx in range(len(filtered) - 1):
            self._draw_segment_px(
                int(round(filtered[idx][0])),
                int(round(filtered[idx][1])),
                int(round(filtered[idx + 1][0])),
                int(round(filtered[idx + 1][1])),
                color,
                thickness,
            )
        pts_str = " ".join(f"{px:.1f},{py:.1f}" for px, py in filtered)
        self.svg_items.append(
            f'<polyline fill="none" stroke="rgb({color[0]},{color[1]},{color[2]})" '
            f'stroke-width="{thickness}" stroke-linejoin="round" stroke-linecap="round" points="{pts_str}"/>'
        )

    def draw_circle(self, cx: float, cy: float, radius: int, color: Tuple[int, int, int]) -> None:
        icx, icy = int(round(cx)), int(round(cy))
        for ox in range(-radius, radius + 1):
            for oy in range(-radius, radius + 1):
                if ox * ox + oy * oy <= radius * radius:
                    self._set_px(icx + ox, icy + oy, color)
        self.svg_items.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{radius}" fill="rgb({color[0]},{color[1]},{color[2]})"/>'
        )

    def draw_text(
        self,
        x: int,
        y: int,
        text: str,
        color: Tuple[int, int, int] = (25, 32, 45),
        scale: int = 2,
        anchor: str = "start",
    ) -> None:
        clean = text.upper()
        char_w = 6 * scale
        total_w = max(0, len(clean) * char_w - scale)
        start_x = x
        if anchor == "middle":
            start_x = x - total_w // 2
        elif anchor == "end":
            start_x = x - total_w

        for c_idx, ch in enumerate(clean):
            glyph = _PLOT_FONT_5X7.get(ch, _PLOT_FONT_5X7[" "])
            gx = start_x + c_idx * char_w
            for row_i, row_bits in enumerate(glyph):
                for col_i in range(5):
                    if (row_bits >> (4 - col_i)) & 1:
                        for sy in range(scale):
                            for sx in range(scale):
                                self._set_px(gx + col_i * scale + sx, y + row_i * scale + sy, color)

        svg_anchor = {"start": "start", "middle": "middle", "end": "end"}.get(anchor, "start")
        font_px = 7 * scale + 2
        esc = (
            text.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
        )
        self.svg_items.append(
            f'<text x="{x}" y="{y + 6 * scale}" fill="rgb({color[0]},{color[1]},{color[2]})" '
            f'font-family="monospace, sans-serif" font-size="{font_px}" font-weight="600" '
            f'text-anchor="{svg_anchor}">{esc}</text>'
        )

    def save_svg_and_png(self, svg_path: str, png_path: str) -> None:
        svg_doc = (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" '
            f'viewBox="0 0 {self.width} {self.height}">\n  '
            + "\n  ".join(self.svg_items)
            + "\n</svg>\n"
        )
        with open(svg_path, "w", encoding="utf-8") as sf:
            sf.write(svg_doc)

        raw_rows = bytearray()
        stride = self.width * 3
        for yy in range(self.height):
            raw_rows.append(0)  # Filter type 0 (None)
            raw_rows.extend(self.buf[yy * stride : (yy + 1) * stride])

        ihdr = (
            self.width.to_bytes(4, "big")
            + self.height.to_bytes(4, "big")
            + bytes([8, 2, 0, 0, 0])
        )
        compressed = zlib.compress(bytes(raw_rows), level=6)
        png_bytes = (
            b"\x89PNG\r\n\x1a\n"
            + _plot_png_chunk(b"IHDR", ihdr)
            + _plot_png_chunk(b"IDAT", compressed)
            + _plot_png_chunk(b"IEND", b"")
        )
        with open(png_path, "wb") as pf:
            pf.write(png_bytes)


def _plot_draw_axes_panel(
    canvas: _PlotCanvas,
    box: Tuple[int, int, int, int],
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    x_ticks: Sequence[Tuple[float, str]],
    y_ticks: Sequence[Tuple[float, str]],
    panel_title: str,
    x_label: str,
    y_label: str,
) -> Callable[[float, float], Tuple[float, float]]:
    bx, by, bw, bh = box
    canvas.fill_rect(bx, by, bw, bh, fill=(255, 255, 255), border=(185, 194, 208))
    xmin, xmax = x_range
    ymin, ymax = y_range

    def _map_pt(xv: float, yv: float) -> Tuple[float, float]:
        px = bx + ((xv - xmin) / max(1e-30, xmax - xmin)) * bw
        py = by + bh - ((yv - ymin) / max(1e-30, ymax - ymin)) * bh
        return (px, py)

    for xv, xlbl in x_ticks:
        px, _ = _map_pt(xv, ymin)
        canvas.draw_line(px, by, px, by + bh, (233, 237, 243), thickness=1)
        canvas.draw_line(px, by + bh, px, by + bh + 5, (110, 120, 135), thickness=1)
        canvas.draw_text(int(round(px)), by + bh + 9, xlbl, color=(65, 75, 90), scale=1, anchor="middle")

    for yv, ylbl in y_ticks:
        _, py = _map_pt(xmin, yv)
        canvas.draw_line(bx, py, bx + bw, py, (233, 237, 243), thickness=1)
        canvas.draw_line(bx - 5, py, bx, py, (110, 120, 135), thickness=1)
        canvas.draw_text(bx - 9, int(round(py)) - 4, ylbl, color=(65, 75, 90), scale=1, anchor="end")

    canvas.draw_text(bx + bw // 2, by - 22, panel_title, color=(20, 28, 42), scale=2, anchor="middle")
    canvas.draw_text(bx + bw // 2, by + bh + 26, x_label, color=(40, 50, 68), scale=2, anchor="middle")
    canvas.draw_text(bx + 8, by + 8, y_label, color=(70, 82, 102), scale=1, anchor="start")
    return _map_pt


def _plot_draw_legend(
    canvas: _PlotCanvas,
    x0: int,
    y0: int,
    items: Sequence[Tuple[str, Tuple[int, int, int]]],
) -> None:
    box_w = max(210, max(len(lbl) for lbl, _ in items) * 7 + 34)
    box_h = len(items) * 18 + 10
    canvas.fill_rect(x0, y0, box_w, box_h, fill=(248, 250, 253), border=(205, 212, 224))
    for idx, (lbl, col) in enumerate(items):
        yy = y0 + 11 + idx * 18
        canvas.draw_line(x0 + 8, yy + 3, x0 + 26, yy + 3, col, thickness=2)
        canvas.draw_text(x0 + 32, yy, lbl, color=(30, 38, 52), scale=1, anchor="start")


def _plot_generate_all(results: Dict[str, Any], plots_dir: str) -> List[str]:
    os.makedirs(plots_dir, exist_ok=True)
    written: List[str] = []

    # 1. Bit-Length Trajectory Across All 3 Battles (n = 0..1000)
    c1 = _PlotCanvas(1240, 680)
    c1.draw_text(
        620, 18, "ULTRAFINITIST ALGEBRAIC BIT-LENGTH PLATEAU (N = 0..1000 STEPS)", scale=2, anchor="middle"
    )
    m_bits = _plot_draw_axes_panel(
        c1,
        box=(95, 78, 1095, 535),
        x_range=(0.0, 1000.0),
        y_range=(0.0, 11000.0),
        x_ticks=[(0, "0"), (200, "200"), (400, "400"), (600, "600"), (800, "800"), (1000, "1000")],
        y_ticks=[
            (0, "0"),
            (2000, "2000"),
            (4000, "4000"),
            (6000, "6000"),
            (8000, "8000"),
            (10000, "10000"),
        ],
        panel_title="MAX BIT-LENGTH ACROSS STATE + EXTRACTED PHYSICAL (Q, P) VS GATE 2 10,000-BIT CAP",
        x_label="STEP N",
        y_label="BITS",
    )
    g2_a = m_bits(0.0, 10000.0)
    g2_b = m_bits(1000.0, 10000.0)
    c1.draw_line(g2_a[0], g2_a[1], g2_b[0], g2_b[1], (205, 45, 45), thickness=2, dashed=True)
    c1.draw_text(int(g2_a[0]) + 14, int(g2_a[1]) - 13, "GATE 2 CEILING: 10,000 BITS", color=(195, 35, 35), scale=1)

    # Naive Rational Verlet O(9.9^n) reference explosion from Cycle 1
    naive_pts = [m_bits(0, 4), m_bits(1, 58), m_bits(2, 559), m_bits(3, 5551), m_bits(4, 11000)]
    c1.draw_polyline(naive_pts, (150, 50, 175), thickness=2, clip_box=(95, 78, 1095, 535))

    b1_bits = results["battles"]["battle1"]["q_runs"]["1/2"]["bit_series"]
    b2_bits_half = results["battles"]["battle2"]["q_runs"]["1/2"]["bit_series"]
    b2_bits_quar = results["battles"]["battle2"]["q_runs"]["1/4"]["bit_series"]
    b3_bits = results["battles"]["battle3"]["q_runs"]["1/2"]["bit_series"]

    c1.draw_polyline(
        [m_bits(float(i), float(b)) for i, b in enumerate(b2_bits_quar)],
        (230, 130, 20),
        thickness=1,
        clip_box=(95, 78, 1095, 535),
    )
    c1.draw_polyline(
        [m_bits(float(i), float(b)) for i, b in enumerate(b2_bits_half)],
        (25, 115, 210),
        thickness=2,
        clip_box=(95, 78, 1095, 535),
    )
    c1.draw_polyline(
        [m_bits(float(i), float(b)) for i, b in enumerate(b1_bits)],
        (20, 155, 95),
        thickness=2,
        clip_box=(95, 78, 1095, 535),
    )
    c1.draw_polyline(
        [m_bits(float(i), float(b)) for i, b in enumerate(b3_bits)],
        (210, 45, 115),
        thickness=2,
        clip_box=(95, 78, 1095, 535),
    )
    _plot_draw_legend(
        c1,
        835,
        118,
        [
            ("NAIVE Q-VERLET O(9.9^N) STALL@4", (150, 50, 175)),
            ("B2 3-BODY Q-STEPPER (DT=1/4, <=7033B)", (230, 130, 20)),
            ("B2 3-BODY Q-STEPPER (DT=1/2, <=3963B)", (25, 115, 210)),
            ("B1 CHAOTIC Q-STEPPER (DT=1/2, <=1603B)", (20, 155, 95)),
            ("B3 KEPLER Q-STEPPER (DT=1/2, <=174B)", (210, 45, 115)),
        ],
    )
    svg1 = os.path.join(plots_dir, "bit_length_vs_step.svg")
    png1 = os.path.join(plots_dir, "bit_length_vs_step.png")
    c1.save_svg_and_png(svg1, png1)
    written.extend([svg1, png1])

    # 2. Battle 3 Kepler Orbit Closure & LRL Precession (dt = 1/4, 40 Orbits)
    c2 = _PlotCanvas(1280, 640)
    c2.draw_text(
        640, 16, "BATTLE 3: ECCENTRIC KEPLER 1/R ORBIT CLOSURE & APSIDAL PRECESSION (DT=1/4, N=1000)", scale=2, anchor="middle"
    )
    m_orb = _plot_draw_axes_panel(
        c2,
        box=(75, 78, 520, 495),
        x_range=(-1.9, 1.9),
        y_range=(-1.9, 1.9),
        x_ticks=[(-1.5, "-1.5"), (-0.75, "-0.75"), (0.0, "0"), (0.75, "0.75"), (1.5, "1.5")],
        y_ticks=[(-1.5, "-1.5"), (-0.75, "-0.75"), (0.0, "0"), (0.75, "0.75"), (1.5, "1.5")],
        panel_title="ORBITAL PLANE (X, Y) OVER 40 ORBITS",
        x_label="X",
        y_label="Y",
    )
    b3_q_qs = results["battles"]["battle3"]["q_runs"]["1/4"]["qs"]
    b3_vv_qs = results["battles"]["battle3"]["baselines"]["verlet"]["1/4"]["qs"]
    b3_rk_qs = results["battles"]["battle3"]["baselines"]["rk4"]["1/4"]["qs"][:350]

    c2.draw_polyline(
        [m_orb(q_f[0], q_f[1]) for q_f in b3_vv_qs],
        (235, 140, 35),
        thickness=1,
        clip_box=(75, 78, 520, 495),
    )
    c2.draw_polyline(
        [m_orb(q_f[0], q_f[1]) for q_f in b3_rk_qs],
        (195, 55, 55),
        thickness=1,
        clip_box=(75, 78, 520, 495),
    )
    c2.draw_polyline(
        [m_orb(q_f[0], q_f[1]) for q_f in b3_q_qs],
        (18, 105, 210),
        thickness=2,
        clip_box=(75, 78, 520, 495),
    )
    fx, fy = m_orb(0.0, 0.0)
    c2.draw_circle(fx, fy, 5, (20, 25, 35))
    _plot_draw_legend(
        c2,
        88,
        492,
        [
            ("Q-STEPPER: CLOSED ELLIPSE (0/1 LRL)", (18, 105, 210)),
            ("VELOCITY VERLET: 565 DEG PRECESSION", (235, 140, 35)),
            ("CLASSICAL RK4: PLUNGES & ESCAPES@334", (195, 55, 55)),
        ],
    )

    m_prec = _plot_draw_axes_panel(
        c2,
        box=(685, 78, 545, 495),
        x_range=(0.0, 1000.0),
        y_range=(-16.0, 4.0),
        x_ticks=[(0, "0"), (250, "250"), (500, "500"), (750, "750"), (1000, "1000")],
        y_ticks=[(-15, "1E-15"), (-10, "1E-10"), (-5, "1E-5"), (0, "1 DEG"), (3, "1000")],
        panel_title="LOG10 LRL PERIHELION ANGLE DRIFT (DEG)",
        x_label="STEP N",
        y_label="LOG10(DEG)",
    )
    for m_key, col in (("euler", (145, 65, 175)), ("verlet", (235, 140, 35)), ("rk4", (195, 55, 55))):
        dr = results["battles"]["battle3"]["baselines"][m_key]["1/4"]["lrl_drifts_deg"]
        c2.draw_polyline(
            [m_prec(float(i), math.log10(max(1e-16, d))) for i, d in enumerate(dr)],
            col,
            thickness=2,
            clip_box=(685, 78, 545, 495),
        )
    # Q-Stepper exact 0/1 LRL vector conservation -> plot at floor -15.5
    c2.draw_polyline(
        [m_prec(0.0, -15.5), m_prec(1000.0, -15.5)],
        (18, 105, 210),
        thickness=3,
        clip_box=(685, 78, 545, 495),
    )
    _plot_draw_legend(
        c2,
        925,
        240,
        [
            ("FORWARD EULER (ESCAPES@788)", (145, 65, 175)),
            ("VELOCITY VERLET (565.1 DEG)", (235, 140, 35)),
            ("CLASSICAL RK4 (ESCAPES@334)", (195, 55, 55)),
            ("Q-STEPPER: EXACT 0/1 LRL DRIFT", (18, 105, 210)),
        ],
    )
    svg2 = os.path.join(plots_dir, "battle3_kepler_orbits_and_precession.svg")
    png2 = os.path.join(plots_dir, "battle3_kepler_orbits_and_precession.png")
    c2.save_svg_and_png(svg2, png2)
    written.extend([svg2, png2])

    # 3. Battle 1 Pullen-Edmonds Chaotic Phase Space & Energy Conservation
    c3 = _PlotCanvas(1280, 640)
    c3.draw_text(
        640, 16, "BATTLE 1: PULLEN-EDMONDS CHAOTIC QUARTIC OSCILLATOR (DT=1/2, N=1000)", scale=2, anchor="middle"
    )
    m_b1p = _plot_draw_axes_panel(
        c3,
        box=(75, 78, 520, 495),
        x_range=(-2.3, 2.3),
        y_range=(-2.3, 2.3),
        x_ticks=[(-2.0, "-2"), (-1.0, "-1"), (0.0, "0"), (1.0, "1"), (2.0, "2")],
        y_ticks=[(-2.0, "-2"), (-1.0, "-1"), (0.0, "0"), (1.0, "1"), (2.0, "2")],
        panel_title="CONFINED CHAOTIC (Q1, Q2) ORBIT ON M_H0",
        x_label="Q1",
        y_label="Q2",
    )
    b1_qs = results["battles"]["battle1"]["q_runs"]["1/2"]["qs"]
    c3.draw_polyline(
        [m_b1p(q_f[0], q_f[1]) for q_f in b1_qs],
        (20, 145, 90),
        thickness=1,
        clip_box=(75, 78, 520, 495),
    )
    m_b1e = _plot_draw_axes_panel(
        c3,
        box=(685, 78, 545, 495),
        x_range=(0.0, 1000.0),
        y_range=(-17.0, 6.0),
        x_ticks=[(0, "0"), (250, "250"), (500, "500"), (750, "750"), (1000, "1000")],
        y_ticks=[(-16, "0/1"), (-10, "1E-10"), (-4, "1E-4"), (0, "1E0"), (4, "1E4")],
        panel_title="LOG10 ENERGY ERROR |H_N - H_0| (DT=1/2 & 1/4)",
        x_label="STEP N",
        y_label="LOG10(|DH|)",
    )
    for m_key, dt_k, col in (
        ("euler", "1/2", (145, 65, 175)),
        ("verlet", "1/2", (215, 50, 50)),
        ("verlet", "1/4", (235, 140, 35)),
        ("rk4", "1/2", (90, 105, 130)),
    ):
        dhs = results["battles"]["battle1"]["baselines"][m_key][dt_k]["dh_series"]
        c3.draw_polyline(
            [m_b1e(float(i), min(6.0, math.log10(max(1e-17, v)))) for i, v in enumerate(dhs)],
            col,
            thickness=2,
            clip_box=(685, 78, 545, 495),
        )
    c3.draw_polyline(
        [m_b1e(0.0, -16.2), m_b1e(1000.0, -16.2)],
        (20, 145, 90),
        thickness=3,
        clip_box=(685, 78, 545, 495),
    )
    _plot_draw_legend(
        c3,
        890,
        245,
        [
            ("EULER DT=1/2 (BLEWUP@13)", (145, 65, 175)),
            ("VERLET DT=1/2 (BLEWUP@27)", (215, 50, 50)),
            ("VERLET DT=1/4 (|DH|=0.129)", (235, 140, 35)),
            ("RK4 DT=1/2 (|DH|=1.591)", (90, 105, 130)),
            ("Q-STEPPER: EXACT 0/1 |DH|=0", (20, 145, 90)),
        ],
    )
    svg3 = os.path.join(plots_dir, "battle1_chaotic_phase_and_energy.svg")
    png3 = os.path.join(plots_dir, "battle1_chaotic_phase_and_energy.png")
    c3.save_svg_and_png(svg3, png3)
    written.extend([svg3, png3])

    # 4. Battle 2 3D 3-Body SO(3) Cluster COM Orbits & 7-Invariant Conservation
    c4 = _PlotCanvas(1280, 640)
    c4.draw_text(
        640, 16, "BATTLE 2: 3D 3-BODY SO(3) MOLECULAR CLUSTER & 7-INVARIANT CONSERVATION (DT=1/4)", scale=2, anchor="middle"
    )
    m_b2o = _plot_draw_axes_panel(
        c4,
        box=(75, 78, 520, 495),
        x_range=(-1.6, 1.6),
        y_range=(-1.6, 1.6),
        x_ticks=[(-1.5, "-1.5"), (-0.75, "-0.75"), (0.0, "0"), (0.75, "0.75"), (1.5, "1.5")],
        y_ticks=[(-1.5, "-1.5"), (-0.75, "-0.75"), (0.0, "0"), (0.75, "0.75"), (1.5, "1.5")],
        panel_title="COM-FRAME (X, Y) ORBITS (BODIES 1-3, N<=250)",
        x_label="X - X_COM",
        y_label="Y - Y_COM",
    )
    b2_qs = results["battles"]["battle2"]["q_runs"]["1/4"]["qs"][:251]
    body_cols = [(25, 110, 210), (220, 60, 60), (25, 155, 95)]
    for b_idx, b_col in enumerate(body_cols):
        pts_b = []
        for q_f in b2_qs:
            cx = (q_f[0] + q_f[3] + q_f[6]) / 3.0
            cy = (q_f[1] + q_f[4] + q_f[7]) / 3.0
            pts_b.append(m_b2o(q_f[3 * b_idx] - cx, q_f[3 * b_idx + 1] - cy))
        c4.draw_polyline(pts_b, b_col, thickness=2, clip_box=(75, 78, 520, 495))
    _plot_draw_legend(
        c4,
        90,
        505,
        [("BODY 1", body_cols[0]), ("BODY 2", body_cols[1]), ("BODY 3", body_cols[2])],
    )

    m_b2e = _plot_draw_axes_panel(
        c4,
        box=(685, 78, 545, 495),
        x_range=(0.0, 1000.0),
        y_range=(-17.0, 6.0),
        x_ticks=[(0, "0"), (250, "250"), (500, "500"), (750, "750"), (1000, "1000")],
        y_ticks=[(-16, "0/1"), (-12, "1E-12"), (-6, "1E-6"), (0, "1E0"), (4, "1E4")],
        panel_title="LOG10 MAX(|DH|, ||DP||, ||DL||) DRIFT (DT=1/4)",
        x_label="STEP N",
        y_label="LOG10(ERR)",
    )
    for m_key, col in (("euler", (145, 65, 175)), ("verlet", (235, 140, 35)), ("rk4", (195, 55, 55))):
        br = results["battles"]["battle2"]["baselines"][m_key]["1/4"]
        comb_err = [
            max(dh, dp, dl)
            for dh, dp, dl in zip(br["dh_series"], br["dp_series"], br["dl_series"])
        ]
        c4.draw_polyline(
            [m_b2e(float(i), min(6.0, math.log10(max(1e-17, v)))) for i, v in enumerate(comb_err)],
            col,
            thickness=2,
            clip_box=(685, 78, 545, 495),
        )
    c4.draw_polyline(
        [m_b2e(0.0, -16.2), m_b2e(1000.0, -16.2)],
        (25, 110, 210),
        thickness=3,
        clip_box=(685, 78, 545, 495),
    )
    _plot_draw_legend(
        c4,
        875,
        245,
        [
            ("EULER DT=1/4 (BLEWUP@17)", (145, 65, 175)),
            ("VERLET DT=1/4 (|DH|=0.176)", (235, 140, 35)),
            ("RK4 DT=1/4 (|DH|=2.62, ||DL||=0.95)", (195, 55, 55)),
            ("Q-STEPPER: ALL 7 INVS EXACT 0/1", (25, 110, 210)),
        ],
    )
    svg4 = os.path.join(plots_dir, "battle2_3body_cluster_and_invariants.svg")
    png4 = os.path.join(plots_dir, "battle2_3body_cluster_and_invariants.png")
    c4.save_svg_and_png(svg4, png4)
    written.extend([svg4, png4])

    return written


def _sim_fmt_sci(val: Optional[float]) -> str:
    if val is None:
        return "N/A"
    if val == 0.0:
        return "0.0"
    if not math.isfinite(val):
        return str(val)
    return f"{val:.3e}"


def run_standalone_suite(n_steps: int = 1000, plots_dir: Optional[str] = None) -> Dict[str, Any]:
    """Run all 3 battles over n_steps, assert exact 0/1 laws & reversibility, print tables, and plot."""
    t_suite0 = time.perf_counter()
    horizons = (1.0, 2.0, 4.0, 8.0)
    battle_titles = {
        1: "Battle 1: 2D Coupled Chaotic Quartic Oscillator (Pullen-Edmonds)",
        2: "Battle 2: 3D Interacting 3-Body SO(3) Molecular Quadrance Cluster",
        3: "Battle 3: 2D Eccentric Newtonian Kepler 1/r (Levi-Civita Conic Feistel)",
    }
    q0_map = {
        1: [float(x) for x in BATTLE1_Q0],
        2: [float(x) for x in BATTLE2_Q0],
        3: [float(x) for x in BATTLE3_Q0],
    }
    p0_map = {
        1: [float(x) for x in BATTLE1_P0],
        2: [float(x) for x in BATTLE2_P0],
        3: [float(x) for x in BATTLE3_P0],
    }
    rhs_map = {1: battle1_rhs_float, 2: battle2_rhs_float, 3: battle3_rhs_float}

    print("=" * 108)
    print("ULTRAFINITIST RATIONAL MECHANICS — STANDALONE ZERO-DEPENDENCY VERIFICATION SUITE")
    print(f"Python stdlib fractions.Fraction engine | Steps per run: {n_steps} | Step sizes: 1/2, 1/4, 1/8")
    print("=" * 108)

    battles_out: Dict[str, Any] = {}
    for bid in (1, 2, 3):
        ref_h = _sim_rk4_reference_at_horizons(
            q0_map[bid], p0_map[bid], rhs_map[bid], horizons=horizons, dt_ref=1.0 / 1024.0
        )
        q_runs: Dict[str, Any] = {}
        b_runs: Dict[str, Dict[str, Any]] = {"euler": {}, "verlet": {}, "rk4": {}}

        print(f"\n{battle_titles[bid]}")
        print("-" * 108)
        print(
            f"{'Method':<18} | {'dt':^5} | {'Steps':^9} | {'Status / Peak Bits':<22} | "
            f"{'Max |dH|':<14} | {'Max ||dP||':<12} | {'Max ||dL||':<14} | {'Err @ T=2':<11}"
        )
        print("-" * 108)

        for dt_q in CANONICAL_DTS:
            lbl = f"{dt_q.numerator}/{dt_q.denominator}"
            qr = _sim_run_exact_q_battle(bid, dt_q, n_steps, ref_h)
            q_runs[lbl] = qr
            q_stat = f"Q<=10k ({qr['peak_combined_bits']}b/{qr['peak_state_bits']}s)"
            q_dh = "0/1 (EXACT)"
            q_dp = "0/1 (EXACT)" if bid == 2 else "—"
            q_dl = "0/1 (EXACT)" if bid in (2, 3) else "—"
            print(
                f"{'Exact Q-Stepper':<18} | {lbl:^5} | {qr['steps_completed']:>4}/{n_steps:<4} | "
                f"{q_stat:<22} | {q_dh:<14} | {q_dp:<12} | {q_dl:<14} | "
                f"{_sim_fmt_sci(qr['horizon_errors'].get('T=2')):<11}"
            )

            for m_key, m_pretty in (
                ("euler", "Forward Euler"),
                ("verlet", "Velocity Verlet"),
                ("rk4", "Classical RK4"),
            ):
                br = _sim_run_float_baseline(bid, m_key, float(dt_q), n_steps, ref_h)
                b_runs[m_key][lbl] = br
                if br["blew_up"]:
                    b_stat = f"BlewUp@{br['blowup_step']}"
                elif br["escaped"]:
                    b_stat = f"Escaped@{br['escape_step']}"
                else:
                    b_stat = "Stable (float64)"
                b_dp = _sim_fmt_sci(br["max_dp"]) if bid == 2 else "—"
                b_dl = _sim_fmt_sci(br["max_dl"]) if bid in (2, 3) else "—"
                print(
                    f"{m_pretty:<18} | {lbl:^5} | {br['steps_completed']:>4}/{n_steps:<4} | "
                    f"{b_stat:<22} | {_sim_fmt_sci(br['max_dh']):<14} | {b_dp:<12} | "
                    f"{b_dl:<14} | {_sim_fmt_sci(br['horizon_errors'].get('T=2')):<11}"
                )

        e12 = q_runs["1/2"]["horizon_errors"].get("T=1")
        e14 = q_runs["1/4"]["horizon_errors"].get("T=1")
        e18 = q_runs["1/8"]["horizon_errors"].get("T=1")
        ord_a = math.log2(e12 / e14) if (e12 and e14 and e14 > 0.0) else 0.0
        ord_b = math.log2(e14 / e18) if (e14 and e18 and e18 > 0.0) else 0.0
        print(
            f"  -> Q-Stepper T=1 empirical convergence orders: "
            f"log2(err(1/2)/err(1/4)) = {ord_a:.3f}, log2(err(1/4)/err(1/8)) = {ord_b:.3f}"
        )
        if bid == 3:
            for lbl in ("1/2", "1/4", "1/8"):
                qr3 = q_runs[lbl]
                vv3 = b_runs["verlet"][lbl]
                print(
                    f"  -> Kepler dt={lbl} ({qr3['perihelia']} orbits): "
                    f"Q-Stepper max LRL drift = {qr3['max_lrl_drift_deg']:.2e} deg (exact 0/1 in Q) "
                    f"vs Velocity Verlet = {vv3['max_lrl_drift_deg']:.2f} deg"
                )

        battles_out[f"battle{bid}"] = {"q_runs": q_runs, "baselines": b_runs}

    planar_summary = _sim_verify_planar_submanifolds()
    print("\nBattle 2 Planar & Equilateral Breathing Submanifold Assertions:")
    for cname, cinfo in planar_summary.items():
        print(
            f"  -> {cname}: {cinfo['steps_verified']} steps verified, "
            f"peak_bits = {cinfo['peak_bits']}, exact 0/1 7-invariants & reversibility = PASS"
        )

    results_bundle: Dict[str, Any] = {
        "n_steps": n_steps,
        "battles": battles_out,
        "planar_submanifolds": planar_summary,
        "total_wall_sec": time.perf_counter() - t_suite0,
    }

    plot_files: List[str] = []
    if plots_dir is not None:
        plot_files = _plot_generate_all(results_bundle, plots_dir)
        summary_json_path = os.path.join(plots_dir, "standalone_summary.json")
        compact_summary = {
            "n_steps": n_steps,
            "total_wall_sec": results_bundle["total_wall_sec"],
            "all_assertions_passed": True,
            "planar_submanifolds": planar_summary,
            "plots_generated": [os.path.basename(p) for p in plot_files],
            "battles": {
                bkey: {
                    dt_k: {
                        "steps_completed": qr["steps_completed"],
                        "peak_combined_bits": qr["peak_combined_bits"],
                        "peak_state_bits": qr["peak_state_bits"],
                        "exact_invariants_0_over_1": qr["exact_invariants_0_over_1"],
                        "exact_reversibility_0_over_1": qr["exact_reversibility_0_over_1"],
                        "horizon_errors": qr["horizon_errors"],
                        "max_lrl_drift_deg": qr["max_lrl_drift_deg"],
                        "perihelia": qr["perihelia"],
                        "wall_sec": qr["wall_sec"],
                    }
                    for dt_k, qr in bval["q_runs"].items()
                }
                for bkey, bval in battles_out.items()
            },
        }
        with open(summary_json_path, "w", encoding="utf-8") as jf:
            json.dump(compact_summary, jf, indent=2)
        print(f"\nGenerated {len(plot_files)} SVG+PNG diagnostic plots in: {plots_dir}")
        for pf in plot_files:
            print(f"  - {pf}")

    print("\n" + "=" * 108)
    print(
        f"ALL EXACT 0/1 CONSERVATION, 0/1 REVERSIBILITY, AND <=10,000-BIT ASSERTIONS PASSED "
        f"IN {results_bundle['total_wall_sec']:.2f}s."
    )
    print("=" * 108)
    return results_bundle


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Standalone Zero-Dependency Ultrafinitist Rational Mechanics Suite"
    )
    parser.add_argument("--steps", type=int, default=1000, help="Steps per run (default: 1000)")
    parser.add_argument(
        "--plots-dir",
        type=str,
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots"),
        help="Directory for diagnostic SVG and PNG plots",
    )
    parser.add_argument("--no-plots", action="store_true", help="Skip plot generation")
    args = parser.parse_args()
    pdir = None if args.no_plots else args.plots_dir
    run_standalone_suite(n_steps=args.steps, plots_dir=pdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
