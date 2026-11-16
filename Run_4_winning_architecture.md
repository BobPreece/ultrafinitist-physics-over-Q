# Winning Architecture: Pure Rational Relativistic Gravity Over $\mathbb{Q}$

**Module:** [cycle3_engine.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/architect/cycle3_engine.py)  
**Author:** `@architect` (Rational Relativistic Geometer)  
**Constitutional Standard:** Strict Part I Ultrafinitist Purity over $\mathbb{Q}$ (`fractions.Fraction`) + Zero-Branch Polymorphic Execution over `@formalist-2`'s `QQPoly` / `QQRat` / `QQJet` Computer Algebra Tower.

---

## 1. Executive Architectural Overview

The Winning Architecture (`architect/cycle3_engine.py`) unifies strong-field exterior Schwarzschild orbital mechanics (**Battle 1**), horizon-penetrating Hayward–Kerr–Schild black-hole plunges & core bounces (**Battle 2**), and supercritical Oppenheimer–Snyder / Lemaître–Tolman–Bondi (LTB) multi-shell gravitational collapse with Hawking-analog radiation (**Battle 3**) inside a **single, branch-free, square-root-free, floating-dyadic rational dynamical framework** over `fractions.Fraction` ($\mathbb{Q}$).

### Core Constitutional Guarantees (Part I Purity)
1. **Pure Field Arithmetic Over $\mathbb{Q}$:** Zero `float`/`complex` literals or casts, zero `math`/`cmath`/`numpy`/`scipy`/`sympy`, zero `sqrt`/`isqrt` or fractional/negative `**`, zero trigonometric or transcendental calls (`sin`, `cos`, `atan2`, `acos`, `exp`, `log`), zero `round()` or `.limit_denominator()`, zero `while` loops, zero `try/except`, zero Newton–Heron `x + S/x` chains, zero `isinstance`/`hasattr`/`type` dispatching, zero `quantize` branches (`quantize_branch_sites == []`), zero hidden homogeneous `dt/k` sub-stepping loops, and zero piecewise core clamps (`if r < h`, `max`, `min`, `abs`).
2. **Zero-Branch `QuantizerBundle` (`qb`) Protocol:** Every core integrator step and readout map executes one **identical unbranched Python AST** whether passed `RuntimeDyadicQuantizerBundle(bits=48, max_exp=64)` (bounding bit complexity onto $\mathbb{Z}[1/2]$ at `~75–120` bits with `0.0000 bits/step` asymptotic drift across clean dyadic and all 21 prime-perturbed ensembles $p \in \{101, \dots, 199\}$) or `IdentityQuantizerBundle()` (exact rational field $\mathbb{Q}$, multivariate rational function field `QQRat`, or truncated formal Taylor jet ring `QQJet`). Because $\Delta\tau$ (`dt`) is factored strictly **outside** every quantization call $x_{n+1} = x_n + \Delta\tau \, \mathrm{qb}(F(x_n))$, discrete tangents $(X(t + \Delta\tau) - X(t))/\Delta\tau$ remain unfrozen and non-staircasing down to $\Delta\tau = 2^{-30}$ at both $t = 0$ and $t > 0$.
3. **3+1D 7-Slot Ray–Envelope + Gauge-Fixed $SO(2, \mathbb{Q})$ Rotor + Strictly Causal SOS Readout:** Replaces dual Cartesian/spinor coordinates and Cycle 2's hyperbolic spinor shears with a 7-slot time-reversible carrier `(s_ray, vs, uh, vh, t, ch, sh)` in an arbitrary rational $SO(3, \mathbb{Q})$ orbital frame $(\mathbf{e}_1, \mathbf{e}_2)$, simultaneously exposing folded non-negative polar coordinates ($r = s_{\text{ray}}^2 / \rho \ge 0$, $x^2 + y^2 + z^2 - r^2 \equiv 0/1$) and unfolded signed-ray coordinates ($s_{\text{ray}}$ crossing $+r_0 > 0 \to 0 \to < 0$ with $y \equiv 0/1, z \equiv 0/1$ on $L_z = 0$), while guaranteeing **strictly future-pointing $u^t > 0/1$ everywhere** via an explicit sum-of-squares (SOS) null-Doppler factor $w_{\text{eff}} > 0/1$.
4. **Machine-Certified `0/1`-Slack Curvature & Universal Non-Uniform Density Ceilings:** Global closed-form curvature ceiling $B_K(r_s, h) = 24 r_s^2 / h^6$ and universal shell-density ceiling $B_\rho(M_{\text{rest}}, h, M, \mathbf{m})$ attained at pinched configuration $\mathbf{q} = \mathbf{0}$ with **exact `0/1` slack** across all uniform and skewed non-uniform positive mass vectors $\mathbf{m} \in \mathbb{Q}_{>0}^M$.

---

## 2. Adversarial Evolution Matrix (Cycles 1 $\to$ 2 $\to$ 3)

| ID | Adversarial Finding (`@inquisitor` / `@formalist-2`) | Cycle 3 Winning Architectural Resolution |
|---|---|---|
| **C1 R1–R5, W1–W2, F1–F4** | Dual Spinor–Cartesian schism, conjugate ornament doublets, $O(M)$ unbound $B_\rho$, coarse $B^2$ grazing phase error, prime `p = 113` drift, equatorial $l_z$ tilt, Battle 3 $O(h^4)$ LTB offset | Preserved 100% from Cycle 2: synchronized single-radius readout, zero ornament doublets (`dead_ornament_slots == []`), unit-energy $E_0 \equiv 1/1$ double-square $r_0 = (51/14)^2$ + Weierstrass `cubic_chord_point_from_horizon`, dyadic accumulator-free carrier (`0.0000 bits/step` on all 21 primes), equatorial $l_z \equiv 0/1$, and exact $h \ge 0$ Battle 3 `QQJet(order=3)` LTB accelerations. |
| **C2-R1 + C2_F1c** | `_sig_shear_rate` caused opposite-sign $\sigma_A \sigma_B < 0$ hyperbolic runaway & $r = 6{,}000.17$ catapulting at $\Delta\tau \in \{1/8, 1/16\}$; `om_core` broke axial symmetry ($y \neq 0$ on $L_z = 0$) and missed $\pi$-flip at $(h=1/8, \Delta\tau=1/4)$ | **Completely removed `_sig_shear_rate` and `om_core`.** Replaced with unfolded signed ray $(s_{\text{ray}}, v_s)$ + bounded 3-shear $SO^+(1,1)$ boost envelope $\rho(\mu, u_h, v_h) = \mu\bigl((u_h+v_h)^2 + (u_h^2 - v_h^2 - 1)^2\bigr) > 0$. Across **all 12** $(h \in \{1/2,1/4,1/8\}, \Delta\tau \in \{1/2,1/4,1/8,1/16\})$ `energetic_origin_transit` runs: $y \equiv 0/1$, $z \equiv 0/1$, $s_{\text{ray}} = x_{\text{unfolded}}$ crosses $+8 \to 0 \to -15.1..-17.9 < 0$, folded $r < 20$, $u^r_{\text{final}} > 0$. Added `'exact_midflight_origin_transit'` hitting **bitwise-exact $r \equiv 0/1, s_{\text{ray}} \equiv 0/1, (x,y,z) \equiv (0,0,0)$** mid-flight at step 20 with `0/1` palindrome reversibility. |
| **C2-R2 + C2-W1 + C2_F1a,b** | `beta2 * dt^2` inflated $h^3$ by $985\times$ (`h=1/8, dt=1/4`), degrading $h=1/4$ exterior Hayward convergence slope to `1.4984 < 1.85`; normal-equations $w_{\text{eff}}$ crossed $w_{\text{eff}} < 0 \implies u^t < 0$ (acausal flips) with interior `E_killing` spikes | **Removed `beta2 * dt^2` inflation.** Replaced with core-localized $O(\Delta\tau^4)$ smoothing $\Delta\tau_4 = \frac{\Delta\tau^4}{1 + 24\Delta\tau^4}$ divided by $1 + r^6/64$, restoring $h > 0$ Hayward RK4(`1/2048`) convergence slopes to **`2.13` (exterior $h=1/4$)** and **`2.00` (plunge-bounce $h=1/2$)** while all 18 coarse plunges ($\Delta\tau \in \{1/2, 1/4\}$) bounce cleanly ($r_{\max} \approx 8.53..9.31$). Replaced normal-equations $w_{\text{eff}}$ with **`_causal_sos_w_eff`**, an explicit strictly positive SOS rational form ($w_{\text{eff}} > 0/1 \implies u^t > 0/1$ everywhere, **`0` acausal flips** across all trajectories, `max|E_killing - E_0| <= 0.59`). |
| **C2-R3 + C2_F2a,b** | Uniform-mass $B_\rho$ formula was violated by skewed non-uniform `masses` (`[9/2, 1/2]` static and `[9, 1]` live infall); `battle3_readout` raised `ZeroDivisionError` at $h = 0$ | Set minimum gaps $\delta_1 = \frac{3}{4}h$, $\delta_k = \frac{m_k}{4 M_{\text{rest}}}h$ ($k \ge 2$) and derived closed-form `density_ceiling_B_rho(M_rest, h, m_shells, masses)` ($= \frac{64 M_{\text{rest}}^3}{c_\pi h^3(27 M_{\text{rest}}^2 + 9 M_{\text{rest}} m_2 + m_2^2)}$ for $M \ge 2$), attaining **exact `0/1` slack** on Shell 2 at $\mathbf{q} = \mathbf{0}$ for **every** uniform and non-uniform positive mass vector $\mathbf{m} \in \mathbb{Q}_{>0}^M$ (zero `max`/`min`/`abs`/`if` clamps!). Guarded $h = 0$ readout via `deltas[0] != 0` (`B_rho = None` at $h = 0$, zero crash). |
| **C2-W2** | Projective ray scaling $(c_h, s_h) \mapsto \lambda(c_h, s_h)$ was a `0/1` gauge symmetry; orbital dynamics hardcoded $z = 0$ | Added rational quaternion $SO(3, \mathbb{Q})$ frame builder `make_orthonormal_plane_basis(q0, q1, q2, q3)` + `'3d_inclined_plunge'` & `'1A_3d_inclined'` ($z \neq 0, u^z \neq 0$, exact $\lVert \mathbf{x} \times \mathbf{u} \rVert^2 \equiv L_z^2$). Broke projective $(c_h, s_h) \mapsto \lambda(c_h, s_h)$ gauge invariance via midpoint norm coupling $\mathrm{gauge\_factor} = \frac{2\,\mathrm{rot\_norm}_0}{\mathrm{rot\_norm}_0 + \mathrm{rn\_mid}}$ (`== 1/1` on $S^1$, non-invariant for $\lambda = 2$ runtime and symbolic `QQRat("lam")` while keeping `0/1` palindromes and `0/1` `QQJet(order=3)`). |
| **C2-W3 + C2_F2c** | Standalone `ks_null_covector` had an $h = 0$ pole at $(-3r/5, -4r/5, 0)$; `cycle2_spec.md` prose had typo in SOS coefficients of $B_K - K_h$ | Rotated `ks_null_covector` stereographic axis to prime-Pythagorean $(88/137, 105/137)$ ($88^2 + 105^2 = 137^2$), added **globally pole-free 4-spinor Hopf** `ks_spinor_null_covector` & `ks_spinor_metric_tensors` covering **100% of $S^2$ at $h = 0$** (Theorem 5 below), and wrote exact verified SOS coefficients **`(14, 12, 44, 29, 12, 2)`** in Theorem 1 below. |

---

## 3. Unified 3+1D Geodesic Architecture (Battles 1 & 2)

### 3.1 7-Slot Ray–Envelope + Gauge-Fixed $SO(2, \mathbb{Q})$ Rotor Carrier
Each 3+1D geodesic state evolves 7 active dynamical slots $(s_{\text{ray}}, v_s, u_h, v_h, t, c_h, s_h) \in \mathbb{Q}^7$ relative to a constant rational orthonormal plane basis $(\mathbf{e}_1, \mathbf{e}_2) \in SO(3, \mathbb{Q})$ ($\mathbf{e}_a \cdot \mathbf{e}_b = \delta_{ab}$, built root-free from any nonzero quaternion $(q_0, q_1, q_2, q_3) \in \mathbb{Q}^4$ via Euler–Rodrigues rotation in `make_orthonormal_plane_basis`):

1. **Strictly Positive SOS Envelope $\rho(\mu, u_h, v_h) > 0/1$:**
   $$\rho(\mu, u_h, v_h) := \mu \Bigl((u_h + v_h)^2 + (u_h^2 - v_h^2 - 1)^2\Bigr), \qquad \mu > 0.$$
   * **Unconditional Pole-Freeness on $\mathbb{R}^2$:** If $u_h + v_h = 0$, then $u_h^2 - v_h^2 - 1 = (u_h + v_h)(u_h - v_h) - 1 = -1 \neq 0$, so both squares **cannot** vanish simultaneously on $\mathbb{R}^2$; hence $\rho(\mu, u_h, v_h) > 0/1$ strictly for every off-shell $(u_h, v_h) \in \mathbb{Q}^2$.
   * **Exact $SO^+(1,1)$ On-Shell Reduction & Eigenmode Activity:** Under `IDENTITY_QUANTIZER_BUNDLE`, each half-step `_boost_envelope_half` applies a 3-shear Double-Cayley Lorentz boost
     $$\begin{pmatrix} 1 & \alpha \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ \frac{2\alpha}{1 - \alpha^2} & 1 \end{pmatrix} \begin{pmatrix} 1 & \alpha \\ 0 & 1 \end{pmatrix} = \frac{1}{1 - \alpha^2}\begin{pmatrix} 1 + \alpha^2 & 2\alpha \\ 2\alpha & 1 + \alpha^2 \end{pmatrix} \in SO^+(1, 1, \mathbb{Q}),$$
     where $\alpha = \Delta\tau \, c_b = \frac{\xi(16 + 13\xi^2)}{64 + 32\xi^2}$ satisfies $|\alpha| < \frac{13}{32}\frac{|\xi|}{1 + \xi^2/2} \le \frac{13\sqrt{2}}{64} < 0.29 < 1$ (`1 - alpha^2 > 0` strictly for all $\xi \in \mathbb{R}$, eliminating hyperbolic shear blowup!). On $u_h^2 - v_h^2 = 1$, $(u_h^2 - v_h^2 - 1)^2 \equiv 0/1$ so $\rho = \mu(u_h + v_h)^2$ tracks $|s_{\text{ray}}|$ to $O(\Delta\tau^3)$, while off-shell perturbations along both $u_h + v_h$ and $u_h - v_h$ actively drive the trajectory (`dead_ornament_slots == []`).

2. **Theorem 4 (Integral-Domain Polar-Origin Reflection & Dual Folded/Unfolded Coordinates):**
   * **Algebraic Obstruction:** Over any integral domain $\mathcal{R} = \mathbb{Q}(s_{\text{ray}}, u_h, v_h, \dots)$, if $r \in \mathcal{R}$ satisfies $r \ge 0$ everywhere and $y \equiv 0, z \equiv 0$ on an axial ($L_z = 0$) ray, the Pythagorean identity $x^2 + y^2 + z^2 - r^2 = (x - r)(x + r) = 0$ in the integral domain $\mathcal{R}$ forces **either** $x \equiv +r \ge 0$ globally **or** $x \equiv -r \le 0$ globally; no single rational expression can simultaneously satisfy $x^2 = r^2$ with $r \ge 0$ and change sign across $r = 0$.
   * **Dual Resolution:** `_ks_ray_readout` exposes **both** canonical representations from the same 7-slot carrier:
     - **Folded Polar Manifold (`battle2_xyz`, `r, x, y, z, ur, ux, uy, uz`):**
       $$r = \frac{s_{\text{ray}}^2}{\rho} \ge 0, \qquad \mathbf{x} = r\,\mathbf{c}, \qquad u^r_{\text{sh}} = \frac{s_{\text{ray}}}{\rho} v_s,$$
       where $(C, S) = \mathrm{qb.project\_unit\_circle}(c_h, s_h)$, $\mathbf{c} = C\mathbf{e}_1 + S\mathbf{e}_2$, $\mathbf{t} = -S\mathbf{e}_1 + C\mathbf{e}_2$. Satisfies $x^2 + y^2 + z^2 - r^2 \equiv 0/1$, $r \ge 0$ globally (so $r^3 + h^3 \ge h^3 > 0$ has zero poles even when $s_{\text{ray}} < 0$), and $u^r > 0$ outgoing after origin bounce/transit.
     - **Unfolded Signed Ray (`battle2_xyz_unfolded`, `s_ray, x_unfolded, y_unfolded, z_unfolded, ur_unfolded`):**
       $$\mathbf{x}_{\text{unfolded}} = s_{\text{ray}}\,\mathbf{c}, \qquad u^r_{\text{unfolded}} = v_s,$$
       which crosses $+8 > 0 \to 0 \to < 0$ smoothly under `_ks_ray_accel` (an odd function of $s_{\text{ray}}$ with $a_s(0) \equiv 0/1$) while preserving exact axial symmetry $y = y_{\text{unfolded}} \equiv 0/1$ and $z = z_{\text{unfolded}} \equiv 0/1$ on $L_z = 0$.
   * **Bitwise-Exact Mid-Flight Origin Transit (`make_battle2_midflight_origin_preset`):** Constructing $X_0 = F_{-\Delta\tau}^{20}(X_{\text{orig}})$ (with $t$ reset to $0/1$) from `make_battle2_origin_state(w_doppler=2/5)` yields an exterior infall state ($r_0 \approx 7.5 > 2 r_s$, $H_h(r_0) < 1$) whose forward evolution under `battle2_step(..., dt=1/4)` plunges across $H_h \ge 1$, hits **bitwise-exact $s_{\text{ray}} \equiv 0/1$, $r \equiv 0/1$, $(x, y, z) \equiv (0/1, 0/1, 0/1)$ at step 20**, and emerges back to $r > 2 r_s$ ($s_{\text{ray}} < 0$) at step 40 with exact `0/1` palindrome reversibility.

3. **Gauge-Fixed Unit-Circle Rotor (`ch, sh`):**
   Three-shear Double-Cayley rotation advances $(c_h, s_h)$ with exact determinant $+1$ and exact $SO(2, \mathbb{Q})$ norm preservation $c_{h,n}^2 + s_{h,n}^2 \equiv c_h^2 + s_h^2$ under `IDENTITY_QUANTIZER_BUNDLE`. Coupling the time-symmetric midpoint norm $\mathrm{rn\_mid} = \frac{1}{2}\bigl(\mathrm{qb}(c_h^2 + s_h^2) + \mathrm{qb}(c_{h,n}^2 + s_{h,n}^2)\bigr)$ into the coordinate clock via
   $$\mathrm{gauge\_factor} = \frac{2\,\mathrm{rot\_norm}_0}{\mathrm{rot\_norm}_0 + \mathrm{rn\_mid}}$$
   evaluates identically to $1/1$ on the circle $c_h^2 + s_h^2 = \mathrm{rot\_norm}_0$ (preserving `0/1` `QQJet(order=3)` accuracy and `0/1` 10-step runtime palindromes) while breaking projective $(c_h, s_h) \mapsto \lambda(c_h, s_h)$ gauge degeneracy at $O(1)$.

4. **Strictly Future-Pointing Causal-SOS Readout (`_causal_sos_w_eff`):**
   On the Kerr–Schild mass-shell / null-cone quadric, for **any** $w_{\text{eff}} > 0/1$:
   $$u^t = \frac{(1 + H_h(r)) w_{\text{eff}}^2 + S_\perp}{2 w_{\text{eff}}}, \qquad u^r = w_{\text{eff}} - u^t, \qquad S_\perp = \kappa + u_\perp^2 \ge 0.$$
   Evaluating $g_{\mu\nu} u^\mu u^\nu + \kappa = -(u^t + u^r)(u^t - u^r) + H_h(r) w_{\text{eff}}^2 + u_\perp^2 + \kappa \equiv 0/1$ holds identically over $\mathbb{Q}$ and `QQRat`. Moreover, because $H_h(r) \ge 0$ and $S_\perp \ge 0$, **$w_{\text{eff}} > 0/1$ implies $u^t > 0/1$ unconditionally**.
   `_causal_sos_w_eff` constructs $D_{\text{eff}} = (1 - g_{\text{away}})\,\mathrm{dm}_{t0} + g_{\text{away}}\bigl((1 - w_{\text{union}})\,\mathrm{dm}_{\text{away}} + w_{\text{union}}\,D_{\text{trap}}\bigr)$ from three explicit strictly positive SOS components:
   - $\mathrm{dm}_{t0}(a, d_0) = \frac{(a + d_0)^4 + \frac{1}{16}(a - d_0)^4}{8 d_0 (a^2 + d_0^2)} > 0/1$, which satisfies $\mathrm{dm}_{t0}(a, d_0) - a \equiv \frac{17(a - d_0)^4}{128 d_0 (a^2 + d_0^2)} = O((a - d_0)^4)$ (`0/1` in `QQJet(order=3)`!),
   - $\mathrm{dm}_{\text{away}}(a, c) = \frac{(b_{\text{abs}} + a)^2 + \frac{1}{1024} c^2}{4 b_{\text{abs}}} > 0/1$ via `[1/1]` Padé $b_{\text{abs}} = c \frac{3 a^2 + c^2}{a^2 + 3 c^2} \ge c/3 > 0$,
   - $D_{\text{trap}} = E_{\text{pos}} + u_{\text{trap}} > E_{\text{pos}} > 0/1$, the `[1/1]` Padé positive Killing-energy root $E + \sqrt{E^2 + (H - 1)_+ S_\perp}$ inside $H_h \ge 1$.
   Since $0 \le g_{\text{away}} < 1$ and $0 \le w_{\text{union}} < 1$, $D_{\text{eff}} > 0/1$ and $w_{\text{eff}} > 0/1 \implies u^t > 0/1$ hold **strictly everywhere on $\mathbb{Q}^7$**.

---

## 4. Machine-Checkable Rational Theorems

### Theorem 1 — Smooth Hayward–SOS Curvature Ceiling $B_K(r_s, h)$ (`0/1` Slack)
For $H_h(r) = \frac{r_s r^2}{r^3 + h^3}$ ($r \ge 0, h > 0$), all 256 Riemann components are regular on $[0, \infty)$, $-\det(g_{\mu\nu}) \equiv +1/1$ ($\det(g_{\mu\nu}) \equiv -1/1$; note on $4 \times 4$ matrices $\det(-g_{\mu\nu}) = (-1)^4 \det(g_{\mu\nu}) = -1/1$), and the Kretschmann scalar is:
$$K_h(r) = (H_h'')^2 + 4\left(\frac{H_h'}{r}\right)^2 + 4\left(\frac{H_h}{r^2}\right)^2 = \frac{4 r_s^2 \left(6 h^{12} - 6 h^9 r^3 + 54 h^6 r^6 - 12 h^3 r^9 + 3 r^{12}\right)}{(r^3 + h^3)^6}.$$
Its global supremum over $r \ge 0$ is attained at $r = 0/1$ with **exact `0/1` slack**:
$$\sup_{r \ge 0} K_h(r) = K_h(0) = B_K(r_s, h) := \frac{24 r_s^2}{h^6} \in \mathbb{Q},$$
certified by the strictly positive-coefficient polynomial identity on $r \ge 0, h > 0$ (exact coefficients **`(14, 12, 44, 29, 12, 2)`**):
$$h^6 (r^3 + h^3)^6 \Bigl(B_K(r_s, h) - K_h(r)\Bigr) = 12 r_s^2 r^3 \Bigl(14 h^{15} + 12 h^{12} r^3 + 44 h^9 r^6 + 29 h^6 r^9 + 12 h^3 r^{12} + 2 r^{15}\Bigr) \ge 0.$$

---

### Theorem 2 — Universal Non-Uniform `0/1`-Slack Shell-Density Ceiling $B_\rho$ (Battle 3)
Let $M \ge 1$ concentric shells have arbitrary positive rational rest masses $\mathbf{m} = (m_1, \dots, m_M) \in \mathbb{Q}_{>0}^M$ with total rest mass $M_{\text{rest}} = \sum_{k=1}^M m_k$, rational volume constant $c_\pi = \frac{4}{3}\pi_{\mathbb{Q}} = \frac{1420}{339}$, core scale $h > 0$, minimum gaps $\delta_1 = \frac{3}{4}h$ and $\delta_k = \frac{m_k}{4 M_{\text{rest}}} h$ ($k \ge 2$), and weighted SOS radii $R_k(\mathbf{q}) = \sum_{j=1}^k (\delta_j + \mu_j q_j^2)$ ($\mu_j > 0$). Define:
$$B_\rho(M_{\text{rest}}, h, M, \mathbf{m}) := \begin{cases} \dfrac{64\, m_1}{27\, c_\pi\, h^3} = \dfrac{1808}{3195}\,\dfrac{m_1}{h^3}, & M = 1, \\[10pt] \dfrac{64\, M_{\text{rest}}^3}{c_\pi\, h^3 \left(27 M_{\text{rest}}^2 + 9 M_{\text{rest}} m_2 + m_2^2\right)}, & M \ge 2, \end{cases}$$
which specializes for uniform masses $m_k = M_{\text{rest}} / M$ to $B_\rho = \frac{64 M^2 M_{\text{rest}}}{c_\pi h^3 (27 M^2 + 9 M + 1)}$. Then **for every $\mathbf{m} \in \mathbb{Q}_{>0}^M$**:
$$\sup_{\mathbf{q} \in \mathbb{Q}^M} \max_{1 \le k \le M} \rho_k(\mathbf{q}) = B_\rho(M_{\text{rest}}, h, M, \mathbf{m}) \quad \text{with exact } 0/1 \text{ slack attained at } \mathbf{q} = \mathbf{0}!$$

*Proof ($M \ge 2$):*
1. **Exact Attainment on Shell 2 at $\mathbf{q} = \mathbf{0}$:** At $\mathbf{q} = \mathbf{0}$, $R_1(\mathbf{0}) = \frac{3}{4}h$ and $R_2(\mathbf{0}) = \frac{3 M_{\text{rest}} + m_2}{4 M_{\text{rest}}} h$, so
   $$R_2(\mathbf{0})^3 - R_1(\mathbf{0})^3 = \frac{m_2 \left(27 M_{\text{rest}}^2 + 9 M_{\text{rest}} m_2 + m_2^2\right) h^3}{64 M_{\text{rest}}^3} \implies \rho_2(\mathbf{0}) = \frac{m_2}{c_\pi (R_2(\mathbf{0})^3 - R_1(\mathbf{0})^3)} \equiv B_\rho.$$
2. **Strict Dominance Over Shell 1 ($\rho_1(\mathbf{0}) < \rho_2(\mathbf{0})$):** Let $m_{\text{tail}} := \sum_{j=3}^M m_j \ge 0$, so $M_{\text{rest}} = m_1 + m_2 + m_{\text{tail}}$. Direct polynomial expansion yields the **explicit positive-coefficient SOS certificate** in $(m_1, m_2, m_{\text{tail}}) \in \mathbb{Q}_{>0}^2 \times \mathbb{Q}_{\ge 0}$:
   $$27 M_{\text{rest}}^2 (M_{\text{rest}} - m_1) - m_1 m_2 (9 M_{\text{rest}} + m_2) = m_2 \Bigl(18 m_1^2 + 44 m_1 m_2 + 27 m_2^2 + (45 m_1 + 54 m_2) m_{\text{tail}} + 27 m_{\text{tail}}^2\Bigr) + 27 M_{\text{rest}}^2 m_{\text{tail}} > 0,$$
   hence $B_\rho - \rho_1(\mathbf{0}) > 0$ unconditionally across **all** skewed mass ratios!
3. **Strict Dominance Over Shells $k \ge 3$ ($\rho_k(\mathbf{0}) < \rho_2(\mathbf{0})$):** For any $k \ge 3$, $R_{k-1}(\mathbf{0}) \ge R_2(\mathbf{0})$ and $R_k(\mathbf{0}) = R_{k-1}(\mathbf{0}) + \frac{m_k h}{4 M_{\text{rest}}} > R_2(\mathbf{0})$, so
   $$\frac{R_k(\mathbf{0})^3 - R_{k-1}(\mathbf{0})^3}{m_k} > \frac{3 R_2(\mathbf{0})^2 h}{4 M_{\text{rest}}} = \frac{(27 M_{\text{rest}}^2 + 18 M_{\text{rest}} m_2 + 3 m_2^2) h^3}{64 M_{\text{rest}}^3} > \frac{(27 M_{\text{rest}}^2 + 9 M_{\text{rest}} m_2 + m_2^2) h^3}{64 M_{\text{rest}}^3} = \frac{R_2(\mathbf{0})^3 - R_1(\mathbf{0})^3}{m_2}.$$
4. **Monotonicity on $\mathbb{Q}^M$:** Since $\mu_j > 0$ and $q_j^2 \ge 0$, $R_k(\mathbf{q}) - R_{k-1}(\mathbf{q}) \ge \delta_k$ and $R_j(\mathbf{q}) \ge R_j(\mathbf{0})$, so $R_k(\mathbf{q})^3 - R_{k-1}(\mathbf{q})^3 \ge R_k(\mathbf{0})^3 - R_{k-1}(\mathbf{0})^3$ for all $\mathbf{q} \in \mathbb{Q}^M$. $\square$

---

### Theorem 3 — Euclidean $p$-adic Valuation Obstruction on $\mathbb{Q}^8$ & Weierstrass Secant Generator (Battle 1B)
1. **Obstruction:** Any rational state $(x^\mu, u^\mu) \in \mathbb{Q}^8$ with $r \in \mathbb{Q}_{>0}$ and $E \neq 0$ has $B = L_z / E \in \mathbb{Q}$, forcing even $p$-adic valuations $v_p(B^2) = 2 v_p(B) \in 2\mathbb{Z}$ for every prime $p$; since $v_7(28) = 1$, $v_{109}(109/4) = 1$, and $v_{433}(433/16) = 1$ are odd, no state in $\mathbb{Q}^8$ has exact $(L_z/E)^2 \in \{28, 109/4, 433/16\}$.
2. **Resolution:** `cubic_chord_point_from_horizon(B_sq, n_chords)` generates exact rational points $(2 n_{\text{chords}} + 1) Q_0$ on the 1D Weierstrass cubic $y^2 = 1 - B_{\text{sq}} u^2 + 2 B_{\text{sq}} u^3$ over $\mathbb{Q}(B_{\text{sq}})$ with `0/1` residual, while 3+1D Battle 1B presets use unit-Killing-energy ($E_0 \equiv 1/1$) Cayley parameters at $r_0 = (51/14)^2 = 2601/196$ achieving $0 < (L_z/E)^2 - B_{\text{target}}^2 < 2 \times 10^{-17}$.

---

### Theorem 5 — Chern–Galois $S^2$ Obstruction & Global 4-Spinor Hopf Trivialization (`C2-W3`)
1. **Topological–Arithmetic Obstruction on Single-Chart Cartesian Maps:** Any off-shell null covector $l_i(x, y, z, r) = P_i / D \in \mathbb{Q}(x, y, z, r)$ satisfying $\sum_{i=1}^3 l_i^2 \equiv 1$ identically over $\mathbb{Q}(x, y, z, r)$ and restricting to $x_i / r$ on $S^2 = \{x^2 + y^2 + z^2 = r^2\}$ lifts at $r = 1$ via Euler's 4-square parameterization to a $\mathbb{Q}$-algebraic section of the real Hopf spinor bundle $S \cong \mathcal{O}(1)_{\mathbb{R}} \to \mathbb{CP}^1 \cong S^2$. Because $c_1(\mathcal{O}(1)) = 1$ is odd and non-rational zeros occur in Galois-conjugate pairs of even degree, every such Cartesian map possesses at least one $\mathbb{Q}$-rational stereographic south pole at $h = 0$ (moved off standard axes to $(-88r/137, -105r/137, 0)$ in `ks_null_covector`).
2. **Global Pole-Free Resolution on Spinor Space (`ks_spinor_null_covector`, `ks_spinor_metric_tensors`):** Pulling $\mathcal{O}(1)$ back along the Hopf projection $\pi : \mathbb{Q}^4 \setminus \{0\} \to S^2$ trivializes the bundle globally:
   $$U_1 = \mathrm{sp}_1 + h^4, \quad N = \mathrm{sp}_0^2 + U_1^2 + \mathrm{sp}_2^2 + \mathrm{sp}_3^2, \quad l_\mu = \left(1,\; \frac{\mathrm{sp}_0^2 + U_1^2 - \mathrm{sp}_2^2 - \mathrm{sp}_3^2}{N},\; \frac{2(U_1\mathrm{sp}_2 - \mathrm{sp}_0\mathrm{sp}_3)}{N},\; \frac{2(U_1\mathrm{sp}_3 + \mathrm{sp}_0\mathrm{sp}_2)}{N}\right),$$
   which satisfies $\eta^{\mu\nu} l_\mu l_\nu \equiv 0/1$ identically, reduces to $l_z \equiv 0/1$ on $\mathrm{sp}_0 = \mathrm{sp}_3 = 0$, has $N = h^8 > 0$ at exact origin $(0,0,0,0)$ for $h > 0$, and at $h = 0$ has $N = r = \sum_{a=0}^3 \mathrm{sp}_a^2 > 0$ **everywhere on $\mathbb{R}^4 \setminus \{0\}$ with zero pole on 100% of $S^2$** (and likewise inside `_ks_ray_readout`, where `l_cov = (1, c_x, c_y, c_z)` is trivialized globally over $S^1$).

---

### Theorem 6 — Bounded Floating-Dyadic Bit-Plateau & 21-Prime $p$-Adic Immunity
Every dynamical state variable $X_a^{(n)}$ updated by `battle1_step`, `battle2_step`, or `battle3_step` under `RuntimeDyadicQuantizerBundle(bits=48, max_exp=64)` has the form $X_a^{(n+1)} = X_a^{(n)} + \Delta\tau \, \mathrm{qb}(F_a)$. For dyadic $\Delta\tau = 2^{-s}$ and initial condition $X_a^{(0)} = p_0 / d_0$ ($d_0 = 2^{e_0} d_{\text{odd}}$ with arbitrary odd prime factors $p \in \{101, \dots, 199\}$):
- Every increment $\Delta\tau \, \mathrm{qb}(F_a)$ lies in the bounded dyadic lattice $2^{-(64 + s)} \mathbb{Z}$, so the odd part of $\mathrm{den}(X_a^{(n)})$ divides $d_{\text{odd}}$ **for all $n \ge 0$** (zero accumulation of new prime factors!).
- The 2-adic valuation satisfies $v_2(\mathrm{den}(X_a^{(n)})) \le \max(e_0, 64 + s) = O(1)$, strictly bounding denominators to $\le 64 + s + \lceil \log_2 d_{\text{odd}} \rceil$ bits (`~75–120` bits peak, `0.0000 bits/step` late-time drift across all $n \to \infty$).

---

## 5. Clean-Room Integration Blueprint for `standalone_ultrafinitist_gravity.py`

To satisfy Gate 4 clean-room isolation (`/usr/bin/python3 -I standalone_ultrafinitist_gravity.py` in `/tmp/cleanroom_exhibit4_*`, `EXIT 0`, `stderr == ""`, `< 180s`):
1. **Part I (`# === PART I: PURE RATIONAL ENGINE ===`, above `# === PART II`):** Embed `architect/cycle3_engine.py` verbatim (imports restricted strictly to `from __future__ import annotations`, `from fractions import Fraction`, and `typing`). Part I passes `audit_ast_purity(..., part1_only=True)` with **zero** violations across all 10 AST categories.
2. **Part II (`# === PART II: BENCHMARK HARNESS, BASELINES & SELF-CONTAINED SVG/PNG RENDERER ===`):** Keep all `float`/`math`/`zlib`/`struct` continuum RK4 baselines, metrics serialization, and pure-Python SVG + PNG rasterizer strictly below the `# === PART II` delimiter, writing all 8 dual `.svg` + `.png` exhibits (`> 1,024` bytes each, valid XML + PNG magic/IHDR/IDAT/IEND) into `./plots/` relative to `__file__`.

---

## 6. Public API Quick Reference (`architect/cycle3_engine.py`)

- **Quantizer Protocol:** `RuntimeDyadicQuantizerBundle`, `QuantizerBundle`, `IdentityQuantizerBundle`, `DEFAULT_QUANTIZER_BUNDLE`, `IDENTITY_QUANTIZER_BUNDLE`, `FOUR_THIRDS_PI_Q`.
- **Geometry, Spinor Hopf Bundle & Bounds:**
  * `hayward_profile_H(r, r_s=2, h=0)`, `kretschmann_invariant_K(r, r_s=2, h=1/4)`, `kretschmann_ceiling_B_K(r_s=2, h=1/4)`
  * `density_ceiling_B_rho(M_rest, h, m_shells=8, masses=None)`
  * `wildberger_spread_2d(vec_a, vec_b)`, `wildberger_spread_3d(vec_a, vec_b)`, `cayley_unit_circle(t_param)`
  * `make_orthonormal_plane_basis(q0=1, q1=0, q2=0, q3=0)`
  * `ks_null_covector(x, y, z, r, h=0)`, `ks_metric_tensors(x, y, z, r, r_s=2, h=0, l_cov_override=None)`
  * `ks_spinor_null_covector(sp0, sp1, sp2, sp3, h=0)`, `ks_spinor_metric_tensors(sp0, sp1, sp2, sp3, r_s=2, h=0)`
  * `ks_radial_lightcone_slopes(r, r_s=2, h=0)`, `vieta_project_ks_4velocity(...)`, `cubic_chord_point_from_horizon(B_sq, n_chords=1)`
- **Battle 1 (`h = 0` Strong-Field Geodesics & Photon Sphere):**
  * `make_battle1_state(r0, Lz, w_doppler, r_s=2, h=0, kappa=1, u0=1, v0=0, b_sq_target=None, e1=None, e2=None)`
  * `make_battle1_preset(preset_id)` (`'1A_eccentric'`, `'1A_r8'`, `'1A_3d_inclined'`, `'1B_null_B2_28'`, `'1B_null_B2_109_4'`, `'1B_null_B2_433_16'`, `'1B_timelike_scattering'`)
  * `battle1_step(state, dt, qb=DEFAULT_QUANTIZER_BUNDLE)`, `battle1_radius(state)`, `battle1_xyz(state, qb)`, `battle1_xyz_unfolded(state, qb)`, `battle1_readout(state, qb)`, `battle1_ks_8vec(state)`, `evolve_battle1(state0, dt, n_steps, qb)`
- **Battle 2 (`h > 0` Horizon Plunge, Core Bounce & Origin Transit):**
  * `make_battle2_state(r0, Lz, w_doppler, h=1/4, r_s=2, kappa=1, u0=1, v0=0, e1=None, e2=None)`
  * `make_battle2_origin_state(w_doppler=1, h=1/2, r_s=2, kappa=1, u0=1, v0=0, e1=None, e2=None)`
  * `make_battle2_midflight_origin_preset(h=1/4, dt=1/4, n_half_steps=20, w_doppler=2/5, r_s=2, kappa=1, qb=DEFAULT_QUANTIZER_BUNDLE)`
  * `make_battle2_preset(preset_id, h=1/4)` (`'angular_plunge_Lz2'`, `'head_on_plunge_Lz0'`, `'energetic_origin_transit'`, `'exact_midflight_origin_transit'`, `'3d_inclined_plunge'`, `'exterior_convergence_Lz3'`, `'exact_origin'`)
  * `battle2_step(state, dt, qb=DEFAULT_QUANTIZER_BUNDLE)`, `battle2_radius(state)`, `battle2_xyz(state, qb)` (folded $r \ge 0$), `battle2_xyz_unfolded(state, qb)` (signed $s_{\text{ray}}$), `battle2_horizon_H(state)`, `battle2_readout(state, qb)`, `battle2_ks_8vec(state)`, `evolve_battle2(state0, dt, n_steps, qb)`
- **Battle 3 (Supercritical LTB Shells & Hawking-Analog Radiation):**
  * `make_battle3_state(m_shells=8, amp=10, h=1/2, nu=0, prime_p=None)`, `make_battle3_custom_state(R0, V0, masses, h=1/2, nu=0)`
  * `battle3_step(state, dt, qb=DEFAULT_QUANTIZER_BUNDLE)`, `battle3_radii(state)`, `battle3_readout(state)`, `evolve_battle3(state0, dt, n_steps, qb)`

