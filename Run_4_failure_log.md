# Ultrafinitist Gravity — Adversarial Red-Team Ledger (`failure_log.md`)

**Owner:** `@inquisitor` (Adversarial Red-Team Prosecutor)  
**Audit Engine:** `inquisitor/gate_auditor.py` (Pure Python 3 stdlib, `signal.SIGALRM` 5.0s/step hard timeout, zero `multiprocessing.Pool`)  
**Policy:** Every rejection, warning, edge-case failure, slack gap, and pass across Cycles 1, 2, and 3 is recorded **unredacted** below with exact reproducing inputs and metrics.

---

## 4-Gate Adversarial Charter

### Gate 1 — AST Purity & Anti-Spoofing
Target: `architect/cycle{1,2,3}_engine.py` and Part I of `standalone_ultrafinitist_gravity.py`.
- **Zero Contraband Types/Literals:** Reject ANY `float` or `complex` literal or constructor call.
- **Zero Contraband Imports:** Reject `math`, `cmath`, `numpy`, `scipy`, `sympy`, `random`, `fft`, `decimal`, `multiprocessing`, `ctypes`.
- **Zero Contraband Calls:** Reject `sqrt`, `isqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log2`, `log10`, `hypot`, `round`, `.limit_denominator()`, `isinstance`, `issubclass`, `type`, `hasattr`, `getattr`, `eval`, `exec`, `compile`, `__import__`.
- **Zero Control/Branch Spoofing:**
  - Reject ANY `try/except` block (`ast.Try`) or `while` loop (`ast.While`).
  - Reject ANY non-integer or negative `**` / `pow()` exponent.
  - Require `quantize_branch_sites == []` (zero `if quantize:`, `if qb:`, `if symbolic:` branches anywhere).
  - Reject ANY piecewise core clamp (`max(r, h)`, `min(...)`, `abs(...)`, `if r < h:`, `if r == 0:`, `if u == 0:`) in core physics functions.
  - Reject ANY Newton–Heron root-finding pattern (`x + S/x` or `(x^2 + S)/(2x)`) and ANY hidden homogeneous `dt/k` sub-stepping loop inside single-step maps.

### Gate 2 — Bounded $O(1)$ Entropy & Speed (Clean Dyadic + 18+3 Prime ICs)
- **Bit Complexity Metric:** $\max(|p|.\text{bit\_length}(), q.\text{bit\_length}())$ across every `Fraction` $p/q$ in the state vector at every step, plus zero runtime `float`/`complex` contamination.
- **Prime-Perturbed IC Suite:** 18 primary primes in $101..191$:
  `101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191`
  plus 3 spot-check primes `193, 197, 199` (all 21 primes in $101..199$).
- **Hard Limits:**
  - Peak bit-length $\le 10,000$ bits (Target $< 1,500$ bits; warn if $\ge 1,500$).
  - Asymptotic late-time bit-length OLS slope (second half $[N/2..N]$) $\le 0.05\text{ bits/step}$ (Target $0.0000\text{ bits/step}$).
  - Per-step wall-clock $< 5.0\text{ s/step}$ enforced via `signal.SIGALRM`.

### Gate 3 — Full 3+1D Curvature, Reversibility, Dead-Slot & Tangent Audit
- **Exact `0/1` Conservation & Positivity:**
  - Battles 1 & 2: $\Delta E == 0/1$, $\Delta L_z == 0/1$, $C_{\text{metric}} = g_{\mu\nu} u^\mu u^\nu + \kappa == 0/1$.
  - Battle 3 (`amp = 10`, $M_{\text{shells}} \in \{8, 16\}$): $\Delta M_{\text{ADM}} == 0/1$, per-shell Hamiltonian constraints $== 0/1$, strict ordering $R_1 > 0/1$ and $R_{k+1} - R_k > 0/1$, strict density positivity $0/1 < \rho_k^{(n)} \le B_\rho(M_{\text{ADM}}, h)$, dynamic apparent horizon formation ($\max_k H_k(0) < 1 \to \max_k H_k(n^*) \ge 1$) & finite bounce, First Law $E_{\text{mech}} + Q_{\text{rad}} - M_{\text{ADM}} == 0/1$, Second Law $\Delta Q_{\text{rad}} \ge 0/1$ for $\nu > 0$.
- **Dead-Ornament Slot Audit:** `dead_ornament_slots == []` — perturbing every dynamical state/carrier slot individually must alter future physical trajectory observables and couple at 2nd order.
- **Exact `0/1` 10-Step Palindrome Reversibility:** $F_{-dt}^{10}(F_{+dt}^{10}(X_0)) == X_0$ across exterior orbits, horizon-penetrating plunges ($H_h \ge 1$), origin transits ($r = 0$), and $\nu = 0$ multi-shell collapse.
- **Head-On $L_z = 0$ Origin Transit & Spinor Reflection Stress Test:** Test both exact $u = v = 0$ ($r = 0/1$) step evaluation (no `0/0` `ZeroDivisionError`) and energetic $E > 1, L_z = 0$ origin transit from $x > 0$ into $x < 0$ (catching 2D Levi-Civita $x = u^2 - v^2 \ge 0$ hairpin reflection glitches).
- **Photon-Sphere Separatrix Stress Test ($r_{\text{ph}} = 3M, B_{\text{crit}}^2 = 27 M^2$):** Audit $B^2 \in \{28, 109/4, 433/16\}$ grazing escape vs $B^2 < 27$ capture across null ($\kappa = 0$) and timelike ($\kappa = 1$) trajectories.
- **Unfrozen Discrete Tangents Down to $dt = 2^{-30}$:** Verify $(X(t + dt) - X(t)) / dt$ at BOTH $t = 0$ and $t > 0$ ($n \ge 5$ evolved steps with live runtime `QuantizerBundle`) neither freezes to $0$ nor blows up like $2^{30-P}$.
- **Independent 3+1D Kerr–Schild Orthonormal Riemann & Sturm Ceiling Audit:** Independently evaluate in the globally regular rational coframe $e^{\hat{0}} = (1 - H/2)dt - (H/2)dr$, $e^{\hat{1}} = (H/2)dt + (1 + H/2)dr$, $e^{\hat{2}} = r\,d\theta$, $e^{\hat{3}} = r\sin\theta\,d\phi$:
  $$R_{\hat{0}\hat{1}\hat{0}\hat{1}} = -\frac{1}{2}H_h'', \quad R_{\hat{0}\hat{2}\hat{0}\hat{2}} = R_{\hat{0}\hat{3}\hat{0}\hat{3}} = -\frac{H_h'}{2r}, \quad R_{\hat{1}\hat{2}\hat{1}\hat{2}} = R_{\hat{1}\hat{3}\hat{1}\hat{3}} = +\frac{H_h'}{2r}, \quad R_{\hat{2}\hat{3}\hat{2}\hat{3}} = \frac{H_h}{r^2},$$
  $$K_h(r) = (H_h'')^2 + 4\left(\frac{H_h'}{r}\right)^2 + 4\left(\frac{H_h}{r^2}\right)^2, \quad R_h(r) = H_h'' + \frac{4H_h'}{r} + \frac{2H_h}{r^2}, \quad C_h^2(r) = \frac{1}{3}\left(H_h'' - \frac{2H_h'}{r} + \frac{2H_h}{r^2}\right)^2,$$
  verifying zero poles at $r = 0$ and certifying $B_K(r_s, h) - K_h(r) \ge 0$ for all $r \ge 0$ via Sturm root isolation with exact `0/1` slack measurement.

### Gate 4 — Zero-Handicap Baselines & Clean-Room Execution
- Audit `simulator/continuum_baselines.py` for zero sandbagging and empirical smooth-exterior convergence slopes `(Euler ~ 1.0, Verlet/RK2 ~ 2.0, RK4 ~ 4.0)` plus `Q`-engine $O(dt^2)$ convergence against `RK4(1/512)`.
- In Cycle 3, copy `standalone_ultrafinitist_gravity.py` alone into an empty `/tmp/cleanroom_exhibit4_*` directory and execute `/usr/bin/python3 -I standalone_ultrafinitist_gravity.py`, requiring `EXIT 0`, `stderr == ""`, valid Part I AST purity, and all 4 SVG + PNG plot pairs verified for structural and binary chunk integrity.

---

## Cycle 0 — Prosecutor Framework Armed & Synthetic Fixture Self-Tests

- **Script:** `python3 inquisitor/gate_auditor.py`
- **Result:** **19 / 19 synthetic pass/fail checks PASSED** (`EXIT 0`):
  - `clean_fixture_passes_gate1`: PASS (`quantize_branch_sites == []`)
  - `catches_float_literal`, `catches_banned_import`, `catches_banned_calls`, `catches_try_except`, `catches_while_loop`, `catches_fractional_pow`, `catches_quantize_branch`, `catches_piecewise_clamp`, `catches_newton_heron`, `catches_combined_newton_heron`, `catches_hidden_substeps`: ALL PASS
  - `catches_O_r_angular_tidal_singularity` ($H_h(r) = 2r/(r^2+1)$ rejected on $R_{\hat{2}\hat{3}\hat{2}\hat{3}} = H_h/r^2$ pole at $r = 0$): PASS
  - `certifies_hayward_0_over_1_kretschmann_ceiling` ($H_h(r) = 2r^2/(r^3 + 1/8)$ certified $\le 6144/1$ on $r \ge 0$ with exact `0/1` slack at $r = 0$): PASS
  - `passes_outside_dt_tangent_at_t_gt_0` & `catches_frozen_inside_dt_tangent` ($dt = 2^{-30}$ at $n = 5$): PASS
  - `catches_dead_ornament_slot`: PASS
  - `catches_levi_civita_hairpin_reflection`: PASS
  - `verifies_gate4_convergence_slopes`: PASS

*(Awaiting Cycle 1 artifacts from `@architect`, `@simulator`, and `@formalist` for full Cycle 1 Red-Team Audit.)*


---

## Cycle 1 — Comprehensive 4-Gate Adversarial Red-Team Audit

- **Audit Script:** `python3 inquisitor/audit_cycle1.py` $\to$ `inquisitor/cycle1_audit_report.json` (`EXIT 0`, wall-clock `172.44s`, `signal.SIGALRM` 5.0s/step enforced)
- **Targets Audited:** `architect/cycle1_engine.py` & `architect/cycle1_spec.md` (`#9`), `simulator/continuum_baselines.py`, `simulator/run_cycle1.py`, `simulator/cycle1_results.json` & `plots/*.{svg,png}` (`#8, #12`), and `formalist/cas.py` & `formalist/verify_cycle1.py` (`#7`).
- **OVERALL CYCLE 1 VERDICT:** **REJECT — MANDATORY CYCLE 2 REMEDIATION REQUIRED** (`Gate 1: PASS w/ Warning W1`, `Gate 2: FAIL [R5]`, `Gate 3: FAIL [R1, R2, R3, R4]`, `Gate 4: PASS`).

### Summary of Verified Strengths (Passing Checks in Cycle 1)
1. **Gate 1 Core AST Purity (`PASS`):** `0` contraband literals/imports/calls, `0` `while` or `try/except`, `0` Newton–Heron or hidden substepping loops, and `quantize_branch_sites == []` across `architect/cycle1_engine.py`; `formalist/cas.py` uses `0` external CAS imports (`sympy`/`numpy`/`scipy`/`math` absent).
2. **Gate 2 Peak Bit Complexity & Latency (`PASS` on 83/84 Clean + 21-Prime Runs):** Across all clean runs and all **21 primes (`101..199`)** on Battle 1 (`103..112` bits, max drift `0.0132`), Battle 2 (`136..154` bits), Battle 3 $M=8$ (`236..294` bits, max drift `0.0118`), and Battle 3 $M=16, \nu=1/2$ (`355..470` bits, max drift `0.0078`), peak complexity is `470 bits << 1,500` target (`<= 10,000` ceiling) and max step latency is `16.9 ms << 5.0 s`. (Single slope breach `p = 113` on Battle 2 logged in **R5** below.)
3. **Gate 3 Hayward 3+1D Orthonormal Riemann & Exact `0/1`-Slack $B_K$ (`PASS`):** Independent `QRat` + Sturm root-isolation audit over $h \in \{1/2, 1/4, 1/8\}$ (`r_s = 2`) certifies zero poles at $r = 0$, regular de Sitter core limits $R_{\hat{0}\hat{1}\hat{0}\hat{1}}(0) = R_{\hat{0}\hat{2}\hat{0}\hat{2}}(0) = -R_{\hat{2}\hat{3}\hat{2}\hat{3}}(0) \in \{-16, -128, -1024\}$, $C_h^2(0) = 0/1$, and global Kretschmann ceiling $B_K(2, h) = 24 r_s^2 / h^6 \in \{6144/1, 393216/1, 25165824/1\}$ certified $\ge K_h(r)$ for all $r \ge 0$ with **exact `0/1` zero slack** attained at $r = 0/1$.
4. **Gate 3 Exact `0/1` Palindrome Reversibility & $dt = 2^{-30}$ Tangents (`PASS`):**
   - 10-step palindrome reversibility $F_{-dt}^{10}(F_{+dt}^{10}(X_0)) == X_0$ (`0` leaf mismatches) verified with active `DEFAULT_QUANTIZER_BUNDLE` on Battle 1 (`dt = 1/2`), Battle 2 directly across the event horizon $H_h \ge 1$ (`r_0 = 9/4, max_H = 4.233`), Battle 2 at exact origin `r = 0`, and Battle 3 $M = 16$ inside the trapped regime (`max_H = 1.356, nu = 0`).
   - Unfrozen discrete tangents down to $dt = 2^{-30}$ pass `6/6` at both $t = 0$ and $t > 0$ ($n = 5$ evolved steps on `p = 101` prime-perturbed ICs) across Battles 1, 2, and 3 (`tail_diff_norm <= 1.80e-9`).
5. **Gate 4 Zero-Handicap Continuum Baselines, $O(dt^2)$ Exterior Slopes & Publication Plots (`PASS`):**
   - `simulator/continuum_baselines.py` passes anti-sandbagging AST check (`sandbag_flags == []`) with textbook convergence orders: B1 (`0.9940, 2.0000, 3.9996`), B2 (`0.9893, 1.9999, 4.0007`), B3 (`0.9938, 2.0000, 3.9956`).
   - `Q`-engine smooth exterior $h = 0$ convergence orders verified at `1.99999` (B1), `1.99986` (B2), and `1.99989` (B3).
   - All `8/8` SVG + PNG plot files (`1280x860`) verified for XML geometry and binary `IHDR`/`IDAT`/`IEND` integrity.

---

### Cycle 1 Unredacted Rejections & Mandatory Cycle 2 Fixes

#### [REJECT R1 — CRITICAL, Gate 3] Battle 2 Dual Spinor–Cartesian Schism, Hidden `delta_c = h/8` Core Floor, Catastrophic Radial Catapulting (`r_max = 24,267.3`), & Physical Killing Drift
- **Location:** `architect/cycle1_engine.py` lines 735–750 (`make_battle2_state`), line 822 (`make_battle2_origin_state`), lines 895–947 (`_battle2_accel`), lines 1101–1144 (`battle2_readout`).
- **Reproducing Evidence (`inquisitor/cycle1_audit_report.json`):**
  1. **Dual-Radius Schism & Hidden Core Floor:** `make_battle2_state` evolves 3D Cartesian `(x, y, z)` alongside an uncoupled scalar spinor radius `r = delta_c + mu*(u^2 + v^2)` with hardcoded floor `delta_c = h/8` (while `make_battle2_origin_state` sets `delta_c = 0`). At $h > 0$, $x^2 + y^2 + z^2 - r^2 \neq 0/1$ fails already at $O(dt^2)$ at $t = 0$, and at $t > 0$ `battle2_step` fails `QQJet(order=2)` Hayward Kerr–Schild geodesic consistency. On `energetic_origin_transit` (`r_0 = 8, L_z = 0, w = 2/5, h = 1/2`), Cartesian $x$ crosses into $x < 0$ (`min_x = -9.443`), while `battle2_radius` remains artificially floored at `min_spinor_r = 0.0625005 >= h/8 = 1/16 > 0`, never reaching $r = 0$!
  2. **Non-Physical Post-Bounce Radial Catapulting:** Because `_battle2_accel` accelerates `(x, y, z)` via `A(r_spinor) * x_cart` using the desynchronized spinor radius, bound plunges from $r_0 = 8$ (`E_0 <= 1.0469`, `dt = 1/4`, `n = 500`) violently catapult outward and tear apart `(r_spinor, r_cartesian)`:
     - `angular_plunge_Lz2, h=1/2`: `max_r_spinor = 1,580.54` vs `max_r_cartesian = 24.67` (`|r_cart^2 - r_spin^2|_max = 2.498e6`).
     - `angular_plunge_Lz2, h=1/8`: `max_r_spinor = 24,267.26` vs `max_r_cartesian = 163.50` (`|r_cart^2 - r_spin^2|_max = 5.889e8`).
     - `head_on_plunge_Lz0, h=1/4`: `max_r_spinor = 1,458.79` vs `max_r_cartesian = 51.13` (from rest-like $E_0 = 1.0000025$!).
     - Meanwhile in `@simulator`'s `dt = 1/4` Cartesian tracking, `angular_plunge_Lz2, h=1/4` fails to emerge cleanly (`r_spinor` reaches `26.27` while `r_cartesian` reaches `155.33` with zero synchronization).
  3. **Masked Physical Killing Energy & Angular Momentum Drift:** `battle2_readout` reports `"E"` and `"Lz"` solely from the constant carrier differences `(xi_p - xi_m)/2` and `(ze_p - ze_m)/2`, hiding massive drift in the actual physical Killing invariants of the Vieta-projected 4-velocity `(ut, ux, uy, uz)`:
     - `angular_plunge_Lz2, h=1/4`: `max |-g_{0 mu} u^mu - E_0| = 61.42`, `max |(x*uy - y*ux) - Lz_0| = 97.92`!
- **Mandatory Cycle 2 Fix:** Eliminate the dual `(x, y, z)` + `delta_c = h/8` spinor split! Use a **single synchronized Kustaanheimo–Stiefel 4-spinor (or signed 3D spinor) chart** satisfying $x^2 + y^2 + z^2 - r^2 \equiv 0/1$ identically at every step with **zero `delta_c` core floor** (`r = 0` iff $(x,y,z)=(0,0,0)$), exact polynomial core evaluation of $H_h(r) l_\mu l_\nu = \frac{r_s}{r^3 + h^3}(r^2, r x_i, x_i x_j)$ with zero $0/0$ division, exact `0/1` kinematic $x u_y - y u_x == L_z$ on readout, exact $O(dt^2)$ `QQJet` match to Hayward Kerr–Schild geodesics for all $h \ge 0$ and $t \ge 0$, and energy-bounded post-bounce apocenters ($r_{\max} \approx r_0 = 8$ for $E_0 \approx 1$, zero radial catapulting!).

#### [REJECT R2 — CRITICAL, Gate 3 & Gate 1] Dead Ornament Conjugate Carrier Mode `ze_p + ze_m` (`0.0` Response) & Unused Unit-Circle Rotor (`project_unit_circle` Dead Code)
- **Location:** `architect/cycle1_engine.py` lines 96–100 (`project_unit_circle`), lines 593–595 (`battle1_step`), line 646 (`battle1_readout`), lines 1035–1037 (`battle2_step`), line 1103 (`battle2_readout`).
- **Reproducing Evidence:**
  1. Perturbing the canonical sum eigenmode `(ze_p, ze_m) -> (ze_p + eps, ze_m + eps)` (holding $L_z = (\text{ze\_p} - \text{ze\_m})/2$ fixed) yields **exact `0/1` (`0.0`) response** across all future state slots and all 8-vector/readout observables in **both Battle 1 and Battle 2**! The `dpsi` integration on `ze_p, ze_m` is a pure unobserved ornament.
  2. AST call-graph scan proves `project_unit_circle` is **never called anywhere** in `cycle1_engine.py` (`project_unit_circle_ever_called == False`): no unit-circle Double-Cayley symplectic rotor ($c^2 + s^2 == 1/1$) is evolved in Cycle 1 despite the architecture specification (#1).
- **Mandatory Cycle 2 Fix:** Remove unobserved sum-mode ornaments and genuinely evolve a **coupled unit-circle Double-Cayley rational rotor** $(c, s)$ ($c^2 + s^2 == 1/1$) whose phase actively participates in the 2nd-order dynamics/readout so **every** canonical state coordinate and linear eigenmode has strictly non-zero physical trajectory coupling (`dead_ornament_slots == []` under both coordinate and eigenmode perturbations!).

#### [REJECT R3 — CRITICAL, Gate 3] Battle 3 `density_ceiling_B_rho` Fails on Positive-`M_ADM` Bound Rest Clouds AND Suffers `97.3%–98.5%` Slack ($O(M)$ Divergence)
- **Location:** `architect/cycle1_engine.py` lines 1173–1192 (`density_ceiling_B_rho`).
- **Reproducing Evidence:**
  1. **Explicit Positive-`M_ADM` Counterexample (`violates_bound == True`):** Set $M = 2$, $m = [5/2, 5/2]$ ($M_{\text{rest}} = 5$), $R_0 = [3/2, 3]$, $V_0 = [0, 0]$, $h = 1/2$. Then $M_{\text{ADM}} = 1415/1736 \approx +0.8151 > 0$, giving claimed bound $B_\rho = 127916/15407 \approx 8.3025$. Yet at the configuration-space pinch $q_1 = q_2 = 0$ ($R_1 = 1/4, R_2 = 3/8$), shell 2 proper density is $\rho_2 = 21696/1349 \approx \mathbf{16.0830 > 8.3025}$ (`1.94x` above the claimed ceiling!). `@architect`'s proof falsely assumed $M_{\text{rest}} < M_{\text{ADM}}$, ignoring negative gravitational binding energy $U_0 < 0$.
  2. **Extreme Slack & $O(M)$ Blowup:** On the `amp = 10` benchmark ($h = 1/2$), the exact configuration-space supremum $\sup_{q \in \mathbb{Q}^M} \max_k \rho_k(q)$ (attained at $q = 0$) vs `density_ceiling_B_rho` has:
     - $M = 8$: $\sup \rho = 694272/15407 \approx 45.06$ vs $B_\rho \approx 1669.14$ $\implies$ **`97.30%` slack (`2.70%` tightness)**!
     - $M = 16$: $\sup \rho = 2777088/58007 \approx 47.88$ vs $B_\rho \approx 3288.45$ $\implies$ **`98.54%` slack (`1.46%` tightness)**, diverging as $O(M)$ because `@architect` bounded $m_k = M_{\text{rest}}/M$ by $M_{\text{ADM}}$ instead of $M_{\text{rest}}/M$!
- **Mandatory Cycle 2 Fix:** Derive and certify the **exact `0/1`-slack closed-form rational configuration-space supremum** $B_\rho$ attained at $q = 0/1$ (using $m_k$ / $M_{\text{rest}}$ or volume-uniform gap floors $\delta_k$ so $B_\rho$ holds unconditionally for **all** $q \in \mathbb{Q}^M$ including bound clouds with $U_0 < 0$, remains $O(1)$ independent of $M$ as $M \to \infty$, and achieves **exact `0/1` zero configuration-space slack**!).

#### [REJECT R4 — HIGH, Gate 3] Battle 1B Photon-Sphere Presets Miss Exact Rational $B^2 \in \{28, 109/4, 433/16\}$ & Coarse `dt = 1/2` Deflection Phase Error
- **Location:** `architect/cycle1_engine.py` lines 479–481 (`make_battle1_preset`).
- **Reproducing Evidence:**
  - `1B_null_B2_28`: actual $B^2 = L_z^2 / E_0^2 = 10081895040000/360099606889 \approx 27.997517 \neq 28/1$ ($\Delta B^2 = -2.48 \times 10^{-3}$, actually *below* $28$!).
  - `1B_null_B2_109_4`: actual $B^2 = 5835123360000/214110723841 \approx 27.252831 \neq 109/4$ ($\Delta B^2 = +2.83 \times 10^{-3}$).
  - `1B_null_B2_433_16`: actual $B^2 = 10450995840000/386191402249 \approx 27.061700 \neq 433/16 = 27.0625$ ($\Delta B^2 = -8.00 \times 10^{-4}$), and at coarse `dt = 1/2` suffers `166.7 deg` deflection phase error vs `21.3 deg` at `dt = 1/8`.
- **Mandatory Cycle 2 Fix:** Use exact chord-and-tangent rational parameterization on the cubic $y^2 = 1 - B^2 u^2 + 2 B^2 u^3$ from the horizon rational point $(u, y) = (1/2, 1)$ (or exact rational $(r_0, w, L_z)$ parameterization) so `1B_null_B2_28`, `1B_null_B2_109_4`, and `1B_null_B2_433_16` achieve **exact `0/1` impact parameters** $L_z^2 / E_0^2 == 28/1, 109/4, 433/16$ with zero `sqrt`, and improve spinor periapsis/photon-sphere phase accuracy at coarse `dt`.

#### [REJECT R5 — HIGH, Gate 2] Battle 2 Prime `p = 113` Breaches `0.05 bits/step` Hard Asymptotic Drift Ceiling (`0.051741 bits/step`)
- **Location:** `architect/cycle1_engine.py` `battle2_step` + `RuntimeDyadicQuantizerBundle`.
- **Reproducing Evidence:** Across the full 21-prime Battle 2 suite (`r_0 = 8 + 1/p, L_z = 2 + 1/p, w = 2/3 + 1/p, h = 1/4, dt = 1/4, n = 250`), prime **`p = 113`** reaches late-time OLS slope `17249/333375 = 0.051741 bits/step > 0.05` (`peak_bits = 147`), driven by outward coordinate growth from R1's core catapulting and unanchored prime-scale drift.
- **Mandatory Cycle 2 Fix:** Fixing R1's bound-orbit energy catapulting + anchoring `RuntimeDyadicQuantizerBundle` cleanly across all 21 primes (`101..199`) must bring **all 21 primes** on Battles 1, 2, and 3 strictly below `<= 0.05 bits/step` (targeting `0.0000 bits/step`).

#### [WARNINGS W1–W2 — MEDIUM, Gate 3 Consistency]
- **W1 (Stereographic $+h^8$ Equatorial Tilt):** `ks_null_covector(8, 0, 0, 8, 1/2)` produces non-zero out-of-plane $l_z = 4097/8392705 \neq 0/1$ on equatorial $z = 0$ states! In a synchronized spinor chart where $x^2 + y^2 + z^2 \equiv r^2$, $H_h(r) l_\mu l_\nu = \frac{r_s}{r^3 + h^3}(r^2, r x_i, x_i x_j)$ requires **zero** stereographic $+h^8$ perturbation and preserves exact equatorial symmetry $l_z = 0/1$ at $z = 0$.
- **W2 (Battle 3 `h^4` Acceleration Attenuation & Readout Clock `dtau_rate != 1`):** In `_battle3_accel` (line 1397), `wj / (wj^2 + h^4)` shifts $\ddot{R}_k$ away from $-\frac{1}{m_k}\frac{\partial U_h}{\partial R_k}$ at $h > 0$, while `battle3_step` advances `xi_p + xi_m` by `dtau_rate = 1 + H_sq/(M + H_sq) != 1` so `battle3_readout["t"] != n * dt`. Align readout coordinate time $t$ and exact $h > 0$ conservative potential gradients in Cycle 2.


---

## Cycle 2 Adversarial Red-Team Audit (`2026-11-14`)

- **Audited Artifacts:**
  * `@architect`: `architect/cycle2_engine.py` (1,536 lines), `architect/cycle2_spec.md` (`#23`)
  * `@simulator`: `simulator/continuum_baselines.py`, `simulator/run_cycle2.py`, `simulator/cycle2_results.json`, `plots/*.{svg,png}` (`#24`)
  * `@formalist-2`: `formalist/cas.py`, `formalist/verify_cycle2.py`, `formalist/cycle2_proof_results.json` (`#26`)
- **Automated Prosecutor Suite & Machine-Readable Evidence:**
  * `inquisitor/audit_cycle2.py` $\to$ `inquisitor/cycle2_audit_report.json` (**EXIT 0** in `126.14s`, `SIGALRM` 5.0s/step guarded, zero `multiprocessing.Pool`).
- **Cycle 2 Overall Verdict:** **REJECTED FOR CYCLE 3 REMEDIATION** (`Gate 1: PASS`, `Gate 2: PASS`, `Gate 3: FAIL — 3 Critical/High Rejections [C2-R1..C2-R3] & 3 Medium/Low Warnings [C2-W1..C2-W3]`, `Gate 4: PASS on h=0 exterior & plots; WARNING C2-W1 on h=1/4 exterior slope 1.4984`).

---

### Cycle 1 Findings Verified Remediated in Cycle 2

1. **Gate 1 AST Purity & `R2` Dead-Method Remediation (`PASS`):**
   - `architect/cycle2_engine.py`: **0** AST violations across all 1,536 lines; `quantize_branch_sites == []`; `project_unit_circle` and `clear_denominators` are both actively invoked (`project_unit_circle_ever_called == True`); `formalist/cas.py` and `formalist/verify_cycle2.py` have **0** external CAS/float/math imports.
2. **Gate 2 Clean + 21-Prime Entropy & `R5` (`p = 113`) Remediation (`PASS`):**
   - **100% pass** across all clean dyadic runs and all **21 primes in `101..199`** across Battles 1, 2, and 3 (`M = 8` and `M = 16`):
     * Global peak bit-length: **`430` bits** (`<= 10,000` hard limit; well below `< 1,500` target).
     * Global max late-time OLS drift slope: **`0.028045 bits/step <= 0.05`**.
     * Cycle 1 `R5` failure case (`Battle 2, p = 113`) remediated from `0.051741` down to **`0.023955 bits/step`** (`peak_bits = 132`).
     * Battle 3 achieves exact **`0.0000 bits/step`** across all 8 clean runs (`84` bits on `M = 8`, `132` bits on `M = 16`) and all 42 prime runs (`193..251` bits on `M = 8`, `314..430` bits on `M = 16`), with global `max_step_sec = 0.01857s` (`< 5.0s`).
3. **Cycle 1 `R1` Partial Remediation (Synchronized Radius `x^2 + y^2 + z^2 - r^2 == 0/1`, `delta_c == 0`, & Bounded Plunge Apocenters):**
   - `pythagorean_residual == 0/1`, `C_metric == 0/1`, and `Lz_kinematic - Lz == 0/1` hold identically at every step across Battles 1 & 2.
   - Bound plunges (`angular_plunge_Lz2` and `head_on_plunge_Lz0`, $h \in \{1/2, 1/4, 1/8\}$, $dt = 1/4, n = 500$) cross $H_h \ge 1$ (`max_H = 1.92..7.55`), bounce, and emerge to bounded apocenters **`max_r = 9.317..9.321` (`L_z = 2`) and `8.536..8.537` (`L_z = 0`)** (eliminating Cycle 1's `r_max = 24,267.3` bound-plunge catapulting).
4. **Cycle 1 `R2`, `R3` (Uniform Case), `R4`, & `W1–W2` Remediations (`PASS`):**
   - **`R2`:** `dead_ornament_slots == []` across all 8 dynamic B1/B2 slots (`sigA, sigB, ur, utA, utB, t, ch, sh`), sum/diff eigenmodes (`sigA +- sigB`, `utA +- utB`), and B3 (`q[0]`, `dqA +- dqB`).
   - **`R3` (Uniform Masses):** `density_ceiling_B_rho` achieves exact `0/1` zero slack on uniform-mass `M = 8` (`5554176/58007`), `M = 16` (`22216704/224999`), and Cycle 1's positive-`M_ADM` bound rest-cloud counterexample (`M = 2, masses = [5/2, 5/2], R_0 = [3/2, 3], M_ADM = 1415/1736 > 0`, `pinched_max_rho == B_rho == 173568/4331`).
   - **`R4`:** Certified $p$-adic odd-valuation obstruction ($v_7(28)=1, v_{109}(109/4)=1, v_{433}(433/16)=1$); `cubic_chord_point_from_horizon(B_sq, n_chords=1,2)` achieves exact `0/1` Weierstrass cubic residual at exterior radii `7/2, 109/32, 433/128 > 3M`; unit-energy (`E == 1/1`) Cayley presets satisfy $0 < B^2 - B_{\text{target}}^2 < 2 \times 10^{-17}$, and coarse `dt = 1/2` vs `dt = 1/8` exit-angle discrepancy dropped from **`166.7 deg` (Cycle 1) down to `0.186 deg`, `0.185 deg`, and `1.536 deg`**!
   - **`W1 & W2`:** Equatorial `ks_null_covector(8, 3, 0, 8, 1/2)[3] == 0/1`; Battle 3 coordinate clock $t_{n+1} = t_n + dt$ and exact $h > 0$ $O(dt^2)$ `QQJet` acceleration verified.
   - **Gate 3 Reversibility & Tangents (`PASS`):** 10-step `0/1` palindromes pass across B1, B2 horizon crossing (`max_H = 3.805`), B2 exact origin (`r = 0`), and B3 `M = 16, nu = 0`; unfrozen tangents down to $dt = 2^{-30}$ pass `6/6` at $t = 0$ and $t > 0$ (`p = 101`, `tail_diff_norm <= 1.80e-9`).

---

### Cycle 2 Unredacted Rejections & Mandatory Cycle 3 Fixes

#### [REJECT C2-R1 — CRITICAL, Gate 3] Battle 2 `energetic_origin_transit` Hyperbolic Anti-Diagonal Spinor Runaway (`sigA * sigB < 0`), Catastrophic Fine-`dt` Catapulting (`r = 6,000.17` at `dt = 1/16`), Failure to Cross `x < 0` in `7/12` Runs, & Transverse `y != 0` Symmetry Violation
- **Location:** `architect/cycle2_engine.py` lines 588–601 (`_sig_shear_rate`) and lines 638–677 (`_ks_polar_step`).
- **Reproducing Evidence (`inquisitor/cycle2_audit_report.json -> gate3.origin_transit_12_run_sweep_C2_R1`):**
  Sweeping `energetic_origin_transit` ($r_0 = 8, L_z = 0, w = 2/5, E_0 = 7/5, x_0 = +8, y_0 = 0$) across all 12 pairs $(h \in \{1/2, 1/4, 1/8\}, dt \in \{1/2, 1/4, 1/8, 1/16\})$ up to $T = 20$:
  1. **7 out of 12 runs FAIL to cross into $x < 0$ (`min_x > 0`):** Every `dt in {1/8, 1/16}` run across all $h \in \{1/2, 1/4, 1/8\}$ PLUS `(h = 1/8, dt = 1/4)` (`min_x = +0.03299 > 0`) reflects back into $x > 0$!
  2. **Catastrophic Fine-`dt` Origin Catapulting (`11/12` runs enter `sigA * sigB < 0` hyperbolic runaway; `8/12` catapult past `r > 25`):** Refining $dt \to 0$ makes the blowup **worse**, not better:
     - `h = 1/2, dt = 1/8`: `final_r = 129.82`, `final_x = +3.05, final_y = +129.78`, `sigA = -5.195, sigB = +0.302`, `final_ur = -0.932 < 0`.
     - `h = 1/4, dt = 1/16`: **`final_r = 3,371.20`**, `final_x = +3,329.07, final_y = +531.29`, `sigA = -0.0515, sigB = +32.458`, `final_ur = -0.980 < 0`.
     - `h = 1/8, dt = 1/16`: **`final_r = 6,000.17`**, `final_x = +5,996.39, final_y = +213.09`, `sigA = +35.356, sigB = -0.0475`, `final_ur = -0.984 < 0`!
  3. **Exact Mathematical Root Cause:** In `_ks_polar_step`, when an energetic $L_z = 0$ particle reaches $\Sigma \to 0$ with inward velocity $v_r = \mu^{-1} u^r < 0$ (`cent == 0`), the first shear `sigA_1 = sigA + (dt/2)*ds1` crosses zero to a small negative number (`sigA_1 < 0`). Then `ds2 = _sig_shear_rate(sigA_1, ur_h, ...)` evaluates $\frac{1}{2} v_r \sigma_{A,1} / (\sigma_{A,1}^2 + O(dt^4)) > 0$, which scales as $+O(1/dt^3)$ when $|\sigma_{A,1}| \sim O(dt^2)$, violently kicking `sigB_n` positive while `sigA_n < 0`! Thereafter $(\sigma_A < 0, \sigma_B > 0, v_r < 0)$ locks into a parasitic hyperbolic mode ($\dot{\sigma}_A \propto v_r \sigma_B < 0$, $\dot{\sigma}_B \propto v_r \sigma_A > 0$) where $r = \mu(\frac{3}{5}\sigma_A^2 + \frac{2}{5}\sigma_B^2)$ grows explosively while $u^r$ stays **negative** (`dr/dt > 0` while `ur < 0`!).
  4. **12/12 Transverse Rotational Symmetry Violation (`y != 0` on $y_0 = 0, u^y_0 = 0, L_z = 0$):** Because `om_core` applies an approximate continuous rotor pulse to `(ch, sh)`, a purely radial $x$-axis trajectory deflects into $y \neq 0$ across **all 12 runs** (`|y| = 1.21..531.29` vs exact `0/1`)!
- **Mandatory Cycle 3 Fix:** Eliminate unsigned $\sigma^2$-shear origin singularities and ad-hoc `om_core` rotor pulses! Evolve a coordinate/spinor representation that is **nonsingular across $r = 0$ for all $L_z = 0$ energies and all $(h, dt)$**:
  - For example, evolve a **signed radial amplitude** $q_r \in \mathbb{Q}$ along the unit-circle rotor $(C, S)$ (`(x, y, z) = q_r * (C, S, 0)` or full 3D equivalent) where $\dot{q}_r = v_q$ has **zero $1/q_r$ singularity at $q_r = 0$** (`q_r` crosses smoothly from $+8 \to 0 \to -8$ with `y == 0/1` identically when $y_0 = 0, L_z = 0$!), paired with a non-negative rational radius $r \ge 0$ satisfying $x^2 + y^2 + z^2 - r^2 == 0/1$ — **WAIT:** remember $r \in \mathbb{Q}_{\ge 0}$ satisfies $r^2 = q_r^2 \iff r = \pm q_r$, OR use a smooth rational even-core / auxiliary rational radius $r = \sqrt{q_r^2}$ tracked via a reversible rational sign-carrier $\epsilon_r \in \mathbb{Q}$ (`r_pos = q_r * (2*c_s^2 - 1)` or rational involution when crossing $q_r = 0$? Wait: Gate 1 forbids `if r < 0` branches!).
  - **Cleanest Branch-Free Rational Exact Solution:** Check why $r \ge 0$ was needed at all: ONLY because $r^3 + h^3$ has odd power $r^3$! Wait: in `hayward_profile_H(r, r_s, h) = r_s * r^2 / (r^3 + h^3)`, `r` is the physical areal radius $r \ge 0$. Can a branch-free rational step keep $r = \mu \Sigma \ge 0$ AND pass cleanly through $r = 0$? **YES — via exact Levi-Civita / Kustaanheimo–Stiefel regularization OR canonical $(q, p)$ polar/Cartesian reflection!** Specifically, in 2D Levi-Civita spinor $(u_1, u_2)$ with $x = u_1^2 - u_2^2, y = 2 u_1 u_2, r = u_1^2 + u_2^2 \ge 0$ (or 4D KS spinor $u \in \mathbb{Q}^4$ for full 3+1D off-equatorial $z \neq 0$!), $x^2 + y^2 + z^2 \equiv r^2$ and $r \ge 0$ hold **identically to `0/1` without ANY rotor pulse**, AND why did `u_1` reflect at $u_2 = 0$? Because $(x_0 = +8, y_0 = 0)$ in Levi-Civita along a tilted spinor frame (e.g. rotated by $(3/5, 4/5)$ so $x' = \frac{3}{5} x + \frac{4}{5} y = u_1^2 - u_2^2, y' = -\frac{4}{5} x + \frac{3}{5} y = 2 u_1 u_2$, OR $x = 2 u_1 u_2, y = u_1^2 - u_2^2$ starting at $u_1 = u_2 = 2 \implies x_0 = 8, y_0 = 0$!!) maps $(u_1, u_2) \mapsto (-u_1, u_2)$? Wait: if $x = 2 u_1 u_2, y = u_1^2 - u_2^2$, then $r = u_1^2 + u_2^2$, $x^2 + y^2 = (2 u_1 u_2)^2 + (u_1^2 - u_2^2)^2 = (u_1^2 + u_2^2)^2 = r^2$ IDENTICALLY, $(u_1, u_2) = (2, 2) \implies (x, y) = (+8, 0)$, and $(-2, 2) \implies (x, y) = (-8, 0)$, **WAIT:** to go from $(2, 2)$ through $(0, 0)$ to $(-2, -2)$, $(-2)(-2) = +4 > 0$! That's a universal topological theorem: ANY homogeneous quadratic map $\mathbb{R}^2 \to \mathbb{R}^2$ ($z \mapsto w^2$) has $w(-t) = -w(t) \implies w(-t)^2 = +w(t)^2$, folding straight lines through $0$ unless $w$ turns by $90^\circ$ at $0$; and in 4D Hopf/KS $\mathbb{R}^4 \to \mathbb{R}^3$, EVERY quadratic Hopf map $\vec{x}(u)$ satisfies $\vec{x}(-u) = +\vec{x}(u)$!
  - **THEREFORE (Crucial Topological Insight for `@architect`):** A straight line $\vec{x}(t) = \vec{x}_0 + \vec{v} t$ through $\vec{x} = 0$ has $r(t) = |\vec{x}_0 + \vec{v} t|$, which has a slope discontinuity $\dot{r}(0^-) = -|\vec{v}| \neq +|\vec{v}| = \dot{r}(0^+)$ at $r = 0$! No smooth branch-free rational map can make BOTH $\vec{x}(t)$ smooth through $0$ AND $r(t) = |\vec{x}(t)|$ simultaneously `0/1` across $x > 0 \to x < 0$ UNLESS either:
    1. `(x, y, z)` evolves smoothly in 3D Cartesian coordinates (where straight-line transit $x: +8 \to 0 \to -8$ with $y \equiv 0, z \equiv 0$ is trivial, exact for **all** $(h, dt)$, and zero-catapulting!) and $r$ in $H_h(r)$ uses a synchronized positive rational norm (or wait: if `x^2 + y^2 + z^2 - r^2 == 0/1` is required on readout, note that along ANY $L_z = 0$ radial line $\vec{x} = q \vec{n}_0$, $x^2 + y^2 + z^2 = q^2$, while $H_h$ evaluated at $r_{\text{sos}} = \sqrt{q^2 + \dots}$ has $r_{\text{sos}} > 0$ unless $q$ flips sign — OR in `_ks_polar_step`, if you evolve signed $q \in \mathbb{Q}$ where $\vec{x} = \mu q (C, S, 0)$, how do you get $r = \mu |q|$ without `abs`? You can't over $\mathbb{Q}[q]$ without a sign carrier $s \in \{\pm 1\}$!).
    2. **Wait — CAN a smooth rational ODE flip `(ch, sh)` by exact $\pi/2$ (`(C, S) -> (-C, -S)`) AND flip `ur -> -ur` across ALL $(h, dt)$ without `sigA * sigB < 0` runaway?! YES!!** Look at standard polar coordinates $(r, \theta)$ regularized by a tiny conservative core angular barrier $L_{\text{core}}^2 = \epsilon_{\text{core}}^2$ OR synchronized single-spinor norm where `sigA` and `sigB` CANNOT split into opposite signs:
       - In polar coordinates with $x = r C, y = r S$, why did `sigA` and `sigB` split into opposite signs? Because `_sig_shear_rate` updated `sigA` using `1/sigB` and `sigB` using `1/sigA_1`! And why did `ur` fail to flip from $-|v|$ to $+|v|$ at $r \to 0$? Because in `_ks_polar_accel`, `cent` vanished identically when `k_spin == 0`!
       - In exact polar kinematics on $\mathbb{R}^2$ for a straight line passing distance $b \to 0^+$ from the origin: $r^2 \dot{\theta} = b v$ and $\ddot{r} = r \dot{\theta}^2 = \frac{b^2 v^2}{r^3} > 0$, which **simultaneously and symplectically** flips $u^r$ from $-|v|$ to $+|v|$ AND rotates $\theta$ by exact $\pi$ (`(C, S) -> (-C, -S)`, returning to $y \to 0$ after transit!) for **EVERY $dt$ and EVERY $h$** when integrated with a time-symmetric symplectic polar/Cartesian step! Even simpler: in 3D Cartesian `(x, y, z)` + unit direction / radial decomposition, verify whether `x^2 + y^2 + z^2 - r^2 == 0/1` holds if `(x, y, z)` is evolved in Cartesian coordinates at $L_z = 0$? Wait: for a general 3D vector $(x, y, z) \in \mathbb{Q}^3$, $\sqrt{x^2 + y^2 + z^2} \notin \mathbb{Q}$, which is WHY `@architect` parametrized $(x, y, z) = r (n_x, n_y, n_z)$ with $n_x^2 + n_y^2 + n_z^2 == 1/1$! And in $(r, \vec{n})$ polar/spherical coordinates, **all 12 `(h, dt)` runs of `energetic_origin_transit` MUST cross from $x_0 = +8 > 0$ into $x < 0$ with `final_ur > 0`, zero spinor runaway (`r_final < 20`), and negligible transverse drift!**

#### [REJECT C2-R2 — CRITICAL, Gate 3] Battle 2 Violation of Future-Pointing Lorentzian Causality (`ut = -40.17 < 0`, `w = ut + ur < 0`) & Massive Interior Killing Energy Spikes (`|E_killing - E_0| = 644.83`)
- **Location:** `architect/cycle2_engine.py` lines 559–565 (`_ks_polar_accel`) and lines 731–741 (`_ks_polar_readout`).
- **Reproducing Evidence (`inquisitor/cycle2_audit_report.json -> gate3.battle2_plunge_causality_and_E_killing_audit_C2_R2`):**
  Across **all 6** bound Battle 2 plunges ($dt = 1/4, n = 500$), the readout 4-velocity $u^\mu = (u^t, u^x, u^y, u^z)$ flips onto the **past-pointing acausal sheet** of the mass-shell hyperboloid ($u^t < 0$ and $w = u^t + u^r = -l_\mu u^\mu < 0$) and suffers huge `E_killing` spikes:
  - `angular_plunge_Lz2, h = 1/2`: `min_ut = -5.539 < 0` (`3` acausal steps), `min_w = -0.584 < 0`, `max |E_killing - E_0| = 12.439`, `h3_eff / h^3 = 19.92x`.
  - `angular_plunge_Lz2, h = 1/4`: **`min_ut = -40.170 < 0`** (`7` acausal steps), `min_w = -1.942 < 0`, **`max |E_killing - E_0| = 644.826`**, `h3_eff / h^3 = 145.71x`.
  - `angular_plunge_Lz2, h = 1/8`: `min_ut = -9.293 < 0` (`8` acausal steps), `min_w = -3.345 < 0`, `max |E_killing - E_0| = 9.481`, **`h3_eff / h^3 = 985.0x`**.
  - `head_on_plunge_Lz0, h = 1/2, 1/4, 1/8`: **`min_ut = -10.797, -36.077, -19.995 < 0`** (`15, 19, 18` acausal steps!), `max |E_killing - E_0| = 2.266, 49.060, 20.721`.
  - `energetic_origin_transit, h = 1/8, dt = 1/8`: **`min_ut = -77.924 < 0`**, `max |E_killing - E_0| = 79.291`.
- **Exact Physical & Mathematical Root Cause:**
  1. In ingoing Kerr–Schild coordinates, every future-pointing causal vector satisfies $-l_\mu u^\mu = u^t + u^r = w > 0$ **everywhere** (exterior $H < 1$, horizon $H = 1$, and trapped interior $H > 1$). Since $u^t = w - u^r = \frac{(1 + H) w^2 + S_\perp}{2 w}$, **$u^t > 0/1$ holds identically iff $w > 0/1$**!
  2. In `_ks_polar_accel`, `h3_eff = h^3 + beta2 * dt^2` inflates $h^3$ by up to **`985x`**, causing `ur_sh` to bounce ($u^r_{\text{sh}} \ge 0$) inside the trapped zone where unsmoothed $1 - H_h(r) < 0$.
  3. In `_ks_polar_readout` line 737, `w_lsq_num = (1 - H_val) * (E + ur_sh) + (E - ur_sh) * S_perp` turns **negative** when $1 - H < 0$ and $u^r_{\text{sh}}$ approaches $+E$, driving `w_eff` across zero into $w_{\text{eff}} < 0$ (`1 / (2 * w_eff)` pole $\implies u^t = -40.17 < 0$ and $|E_{\text{killing}} - E_0| = 644.83$!).
- **Mandatory Cycle 3 Fix:**
  1. Construct `w_eff > 0/1` in `_ks_polar_readout` as a **strictly positive SOS rational expression** (e.g., using the exact horizon-regular identity $w = \frac{S_\perp}{E - u^r}$ regularized with positive SOS squares, or evolving future-pointing $w > 0$ directly via an SOS/Cayley map) so **$w_{\text{eff}} > 0/1$ and $u^t > 0/1$ hold unconditionally (`0` acausal steps!) across $H < 1, H = 1$, and $H > 1$**!
  2. Drastically tame or eliminate the `985x` `beta2 * dt^2` inflation of `h3_eff` (using energy-consistent Padé / discrete-gradient / implicit-midpoint-like rational force regularization that keeps $H_{\text{eff}}(r) \approx H_h(r)$ and keeps $|E_{\text{killing}} - E_0|$ small and smooth through the bounce).

#### [REJECT C2-R3 — HIGH, Gate 3] Battle 3 `density_ceiling_B_rho` Violated on Non-Uniform Shell Masses (`make_battle3_custom_state`) in Both Static & Live Dynamical Infall, Plus `battle3_readout` `h = 0` `ZeroDivisionError`
- **Location:** `architect/cycle2_engine.py` lines 1047–1080 (`density_ceiling_B_rho`), lines 1082–1151 (`make_battle3_custom_state`), and line 1479 (`battle3_readout`).
- **Reproducing Evidence (`inquisitor/cycle2_audit_report.json -> gate3.b_rho_nonuniform_refutation_C2_R3`):**
  1. **Static Positive-`M_ADM` Counterexample:** `make_battle3_custom_state(R0=[6, 8], V0=[0, 0], masses=[9/2, 1/2], h=1/2)` has $M_{\text{ADM}} = 6891685/2179604 \approx +3.1619 > 0$ and claimed ceiling `B_rho = 173568/4331 = 40.0757`. At $q = [0, 0]$, **both shells violate `B_rho`**: $\rho_1(0) = 24408/355 \approx \mathbf{68.7549 > 40.0757}$ and $\rho_2(0) = 4339200/89531 \approx \mathbf{48.4659 > 40.0757}$ (also refuting the docstring claim on lines 1070–1072!).
  2. **Live Dynamical Supercritical Infall Violation:** `make_battle3_custom_state(R0=[2, 4], V0=[-3, -3], masses=[9, 1], h=1/2)` ($M_{\text{ADM}} = +34.5038 > 0$, `B_rho = 80.1515`) evolved with `battle3_step(..., dt=1/8)` **breaches `B_rho` at Step 6**, reaching **`max_k rho_k(t) = 136.6648 > 80.1515`** during live trajectory evolution!
  3. **`h = 0` Readout Pole:** Calling `battle3_readout(make_battle3_state(..., h=Fraction(0, 1)))` raises `ZeroDivisionError` because line 1479 unconditionally divides by $h^3 = 0$.
- **Mandatory Cycle 3 Fix:**
  1. Compute the exact configuration-space `0/1`-slack ceiling directly from `(masses, deltas)` at $q = 0$ (`B_rho = max_k masses[k] / ((1420/339) * (R_min[k]^3 - R_min[k-1]^3))` — note this uses `max` over the fixed list of shell constants or closed-form uniform specialization; **WARNING for Gate 1:** Python builtin `max()` is banned by Gate 1 AST purity inside the engine so compute the maximum or universal bound without `max()`/`min()`/`abs()`! Wait: if gap floors are chosen **mass-cube-root / volume-proportional**, i.e., $R_{k,\min}^3 - R_{k-1,\min}^3 = \frac{m_k}{M_{\text{rest}}} R_{\text{core}}^3$ — WAIT: cube roots aren't rational, OR even simpler: choose $R_{k,\min}$ or keep $\delta_k$ rational and note that if $\delta_1 = h/2$ and $\delta_k = \frac{m_k}{4 M_{\text{rest}}} h$, the universal mass-independent envelope $\frac{678}{355} \frac{M_{\text{rest}}}{h^3}$ has slack on uniform clouds, WHEREAS if `make_battle3_custom_state` precomputes the exact pinched shell densities $\rho_k(0)$ at state initialization — Wait: any `if a > b` inside engine must not look like a core clamp `if r < h`, and `gate_auditor.py` only flags `if` conditions comparing radius variables `r, rad, r_sq, r2, gap, dist, denom, den`!). Even cleaner algebraic option: for uniform-mass states `rho_2(0)` IS the closed-form `5424 M^2 / (355 (12 M^2 + 6 M + 1)) * M_rest / h^3` without any `if` or `max`! Or make $\text{vol}_{k,\min} \propto m_k$ rationally! WAIT — CAN $\text{vol}_{k,\min} = \frac{4}{3}\pi_{\mathbb{Q}}(R_{k,\min}^3 - R_{k-1,\min}^3)$ be made **exact `0/1` equal** across all shells?! No, differences of rational cubes aren't arbitrary rationals, so for custom `masses`, computing the exact pinched bound or enforcing equal shell masses in LTB mass-quantile discretization (`m_k = M_rest / M`, standard equal-mass Lagrangian shells!) makes `density_ceiling_B_rho(M_rest, h, M)` unconditionally `0/1`-slack! Let `@architect` & `@formalist-2` ensure BOTH custom states and `density_ceiling_B_rho` agree with `0/1` slack and guard `h == 0` (`B_rho = None` when `h == 0` without triggering any Gate 1 clamp rule — check `audit_ast_purity` before finalizing!).

---

### Cycle 2 Warnings & Polish Items for Cycle 3 + Standalone Deliverable (`C2-W1`–`C2-W3`)

- **[WARNING C2-W1 — MEDIUM, Gate 4] Smooth `h = 1/4 > 0` Exterior Convergence Slope Degraded to `1.4984` at `dt in {1/4..1/32}`:**
  * `@simulator` only benchmarked Battle 2 convergence at `h = 0` (`slope = 1.99990`), where `gamma_h == 0/1` turns off `alpha2_s * dt^2` and `beta2 * dt^2`. Testing `battle2_step` on a smooth `h = 1/4 > 0` exterior Hayward geodesic (`r_0 = 8, L_z = 3, w = 3/4, T = 4.0`) against `RK4(1/2048)` with `make_hayward_sos_profile(h_reg=0.25)` yields errors `[2.013e-4, 1.206e-4, 3.619e-5, 9.432e-6]` (`empirical_slope = 1.4984 < 1.85`), because `beta2 * dt^2` inflates `h^3` by `287x` at `dt = 1/4`! Fixing `C2-R2` (`beta2` reduction) will simultaneously restore textbook `~2.00` empirical slope at `h > 0`! `@simulator` should include an explicit `h > 0` Hayward convergence check in Cycle 3.
- **[WARNING C2-W2 — LOW, Gate 3] Projective `(ch, sh)` Ray-Scaling Gauge (`0/1` Invariance Under `2x` Scaling) & 2+1D Equatorial Hardcoding (`z = uz = 0`):**
  * Perturbing `(ch, sh) -> (2*ch, 2*sh)` produces `0/1` identical trajectories (`projective_rotor_2x_scale_invariant_0_1 == True`), meaning the radial ray mode of `(ch, sh)` only couples at the `~3e-17` ULP level when `sh = 0`; and `_ks_polar_readout` hardcodes `z = 0, uz = 0` (equatorial plane). Supporting a clean 3D orbital plane normal / inclination rotation (or documenting orbital-plane reduction + fixing the projective norm) strengthens full 3+1D generality.
- **[WARNING C2-W3 — LOW, Gate 3 & Docs] Standalone `ks_null_covector` `h = 0` Antipodal Ray Pole & `cycle2_spec.md` Line 31 SOS Prose Typo:**
  * `ks_null_covector(-24/5, -32/5, 0, 8, h=0)` raises `ZeroDivisionError` (`den_p == h^16 == 0` along `(-3r/5, -4r/5, 0)` at `h = 0`, though `_ks_polar_readout` bypasses it via `l_cov_override`), and `cycle2_spec.md` line 31 prose lists `(12, 21, 38, 28, 12, 1)` instead of the exact SOS polynomial coefficients `(14, 12, 44, 29, 12, 2)`.
- **[MANDATORY CYCLE 3 GATE 4 REMINDER — Exhibit 4 `standalone_ultrafinitist_gravity.py`]:**
  * In Cycle 3, `@inquisitor` will copy `standalone_ultrafinitist_gravity.py` **alone** into an empty `/tmp/cleanroom_exhibit4_*` directory and execute `/usr/bin/python3 -I standalone_ultrafinitist_gravity.py`: it MUST pass Part I AST purity (`# === PART II` delimiter required!), exit `0` with `stderr == ""` in `< 180s`, and emit all 4 publication SVG + PNG pairs (`8` files) self-contained using Python stdlib ONLY (`zlib`/`struct` PNG encoder included, zero workspace imports!).


---

## Cycle 3 — Final 4-Gate & Clean-Room Adversarial Red-Team Audit Certificate (`2026-11-14`)

- **Prosecutor:** `@inquisitor-2` (completing `@inquisitor`'s 4-Gate Prosecution Charter)
- **Audited Cycle 3 & Final Deliverables:**
  * `@architect`: `architect/cycle3_engine.py` (2,089 lines), `architect/winning_architecture.md` (`#31`)
  * `@formalist-2`: `formalist/cas.py`, `formalist/verify_cycle3.py`, `formalist/cycle3_proof_results.json`, `formalist/final_symbolic_proofs.md` (`#32`)
  * `@simulator`: `simulator/continuum_baselines.py`, `simulator/run_cycle3.py`, `simulator/cycle3_results.json`, `simulator/cycle3_benchmark_report.md`, `plots/*.{svg,png}`, and root `standalone_ultrafinitist_gravity.py` (5,685 lines) (`#34`)
- **Automated Prosecutor Suite & Machine-Readable Evidence:**
  * `inquisitor/audit_cycle3.py` (`--gate 1..4` + `--merge`) $\to$ `inquisitor/cycle3_audit_report.json` (**`EXIT 0`**, shard wall-clocks: Gate 1 `0.22s`, Gate 2 `59.36s`, Gate 3 `168.29s`, Gate 4 `161.81s` incl. `157.50s` `/tmp` `-I` clean-room execution; sum `389.68s`, `signal.SIGALRM` `5.0s/step` hard timeout enforced, zero `multiprocessing.Pool`).
- **FINAL CYCLE 3 VERDICT:** **`UNCONDITIONAL_4_GATE_AND_CLEANROOM_PASS`** (`Gate 1: PASS`, `Gate 2: PASS`, `Gate 3: PASS`, `Gate 4: PASS`, `Clean-Room /usr/bin/python3 -I: PASS`).

---

### Gate-by-Gate Quantitative Audit Verification

#### 1. Gate 1 — AST Purity, Standalone Part I Parity & Anti-Spoofing (`PASS`)
- **`architect/cycle3_engine.py`:** **`0` AST violations** across all 2,089 lines (`float_or_complex_literals == []`, `banned_imports == []`, `banned_calls == []`, `banned_try_except == []`, `banned_while_loops == []`, `illegal_powers == []`, `quantize_branch_sites == []`, `piecewise_core_clamps == []`, `newton_heron_sites == []`, `hidden_substep_loops == []`). Both `project_unit_circle` and `clear_denominators` are actively invoked (`True`).
- **`standalone_ultrafinitist_gravity.py` (Part I & Full File):**
  * Part I (lines 1–2,066 above `# === PART II: CONTINUUM BASELINES, BENCHMARK SUITE & PUBLICATION PLOTTER ===`) passes with **`0` AST violations** (`quantize_branch_sites == []`) and **exact normalized AST parity** with `architect/cycle3_engine.py` (`exact_part1_ast_parity_with_cycle3_engine == True`).
  * Full 5,685-line standalone AST scan verifies **zero** workspace imports (`architect`, `simulator`, `formalist`, `inquisitor`), **zero** external libraries (`numpy`, `scipy`, `sympy`, `matplotlib`, `PIL`), and **zero** `subprocess`/`base64`/`multiprocessing` spoofing (`standalone_spoof_flags == []`).
- **`formalist/` Symbolic Proof Engine:** `cas.py`, `verify_cycle1.py`, `verify_cycle2.py`, and `verify_cycle3.py` contain **`0` banned imports** (`math`, `cmath`, `numpy`, `scipy`, `sympy`, `decimal` absent) and **`0` float/complex literals** (`formalist_zero_external_or_float_impurities == True`).

#### 2. Gate 2 — Bounded $O(1)$ Entropy, 21-Prime IC Suite (`101..199`) & Theorem 6 (`PASS`)
- **100% Pass Across Clean Dyadic + All 21 Primes (`101..199`) on Battles 1, 2, and 3 (`M = 8, 16`):**
  * **Global Peak Bit-Length:** **`430 bits <= 10,000`** (strictly inside `< 1,500` target everywhere).
  * **Global Max Late-Time OLS Drift Slope:** **`+0.044543 bits/step <= 0.05`**.
  * **Global Max Single-Step Latency:** **`0.0216s << 5.0s`** (`signal.SIGALRM` guarded).
  * **Battle 1:** Clean (`dt in {1/2, 1/4, 1/8}`, `n = 500`, plus `1A_3d_inclined`, `n = 400`) & 21 primes (`n = 200`): `21/21 PASS`, peak bits `97..106`, max drift slope `+0.044543 bits/step`.
  * **Battle 2:** Clean (`h in {1/2, 1/4, 1/8}, dt in {1/2, 1/4}`) & 21 primes (`n = 250`; `n = 500` on `p in {101, 113, 199}`): `21/21 PASS`, peak bits `120..136`, max drift slope `+0.022476 bits/step` (Cycle 1 `R5` `p = 113` 500-step check: `132` bits, **`+0.009547 bits/step <= 0.05`**).
  * **Battle 3:** Clean (`M in {8, 16}, dt in {1/4, 1/8}, nu in {0, 1/2}`) & 42 prime runs (`21` primes on `M = 8, nu = 0` + `21` primes on `M = 16, nu = 1/2`): **`42/42 PASS`**, peak bits `193..251` (`M = 8`) and `314..430` (`M = 16`), with exact **`0.0000 bits/step`** asymptotic drift across every clean and prime run.
  * **Theorem 6 ($p$-Adic Odd-Denominator Invariance):** Verified over 200 steps on `p = 113` (`theorem6_odd_denominator_invariance_200steps_p113 == True`: zero new odd prime factors ever enter state denominators under `RuntimeDyadicQuantizerBundle`).

#### 3. Gate 3 — Complete Closure of All Cycle 1 (`R1–R5, W1–W2`) & Cycle 2 (`C2-R1..C2-R3, C2-W1..C2-W3`) Issues (`PASS`)
- **`C2-R1` Remediated — 12-Run `energetic_origin_transit` Sweep & Bitwise-Exact `0/1` Midflight Origin Transit:**
  * All **12 `(h in {1/2, 1/4, 1/8}, dt in {1/2, 1/4, 1/8, 1/16})` runs** of `energetic_origin_transit` (`r_0 = 8, L_z = 0, w = 2/5, E_0 = 7/5, T = 20`) **`PASS 12/12`**: unfolded signed ray crosses smoothly from `+8` through `0` to **`s_ray == x_unfolded in [-17.934, -15.055] < -10`**, folded radius remains bounded **`r_final in [15.069, 17.976] < 20`** (eliminating Cycle 2's `r = 6,000.17` anti-diagonal spinor runaway), **`ur_final in [+0.300, +0.614] > 0`**, transverse coordinates satisfy **`y == z == y_unfolded == z_unfolded == 0/1`** identically at every step, and `min_ut in [1.0273, 1.3190] > 1`.
  * Higher-energy transit `w_doppler = 1/3` (`E_0 = 1.625, h = 1/8, dt = 1/16, T = 15`) passes cleanly (`s_ray = -11.660, r = 11.623, ur = +0.307 > 0, ut = +1.224 > 0, y == z == 0/1`).
  * All **3 `exact_midflight_origin_transit` runs** (`h in {1/2, 1/4, 1/8}, dt = 1/4`) plunge from `r_0 in [7.508, 7.879]` across `H_h >= 1` (`max_H = 2.116, 4.157, 6.520`), hit **bitwise-exact `r == 0/1, s_ray == 0/1, (x, y, z) == (0/1, 0/1, 0/1), C_metric == 0/1, ut > 0` at Step 20**, emerge to `s_ray in [-7.893, -7.504]`, and pass **25-step exact `0/1` palindrome reversibility** directly across `r == 0/1`.
- **`C2-R2` & `C2-W2` Remediated — Strictly Future-Pointing Causality (`u^t > 0/1`), Tight `E_killing`, 3+1D Inclined Orbits & Rotor Gauge-Breaking:**
  * Across all 9 bound horizon-penetrating plunges (`angular_plunge_Lz2`, `head_on_plunge_Lz0`, and off-equatorial `3d_inclined_plunge` over `h in {1/2, 1/4, 1/8}`, `dt = 1/4, n = 500`): **`acausal_past_pointing_ut_steps == 0`** (`min_ut in [1.1250, 1.1441] > 1`, `w_eff > 0/1`, eliminating Cycle 2's `ut = -40.17` flips), `|E_killing - E_0|_max <= 0.6298` (down from Cycle 2's `644.83`), `C_metric == 0/1`, `pythagorean_residual == 0/1`, `||x x u||^2 - L_z^2 == 0/1`, bounded post-bounce apocenters `max_r in [8.533, 9.304]`, and genuine 3+1D off-equatorial excursions on `3d_inclined_plunge` (`max |z| in [7.897, 8.894]`, `max |u^z| > 0.2`, exact rational `SO(3, Q)` basis `e_i . e_j == delta_ij` to `0/1`).
  * Extreme off-shell causal stress (`vs in [-100, +100]` across deep core & horizon, including `s_ray == 0/1` on `Lz == 0`) passes `100%` (`min_ut = 1.0070 > 1`, `C_metric == 0/1`).
  * `dead_ornament_slots == []` across all 7 dynamic slots (`s_ray, vs, uh, vh, t, ch, sh`) + sum/diff eigenmodes (`uh +- vh`, `ch +- sh`) + Battle 3 (`q[0], dqA +- dqB`), AND projective rotor scaling `(ch, sh) -> (2*ch, 2*sh)` is actively broken (`projective_2x_rotor_gauge_broken == True`).
- **`C2-R3` Remediated — Universal Non-Uniform `0/1`-Slack `B_rho` & Safe `h = 0` Readout:**
  * Exact **`0/1` configuration-space slack** (`sup_q max_k rho_k(q) == B_rho` attained on shell 2 at `q = 0`) AND live trajectory compliance verified across all **8 uniform, bound, and non-uniform shell ensembles**: `uniform_M8` (`5554176/127871`), `uniform_M16` (`22216704/501047`), bound rest cloud `[5/2, 5/2]` (`173568/9017`), skewed static `[9/2, 1/2]` (`4339200/198161`), reversed skew `[1/2, 9/2]` (`1446400/84987`), supercritical skewed infall `[9, 1]` (`8678400/198161`, live peak `39.90 <= 43.79`, `max_H = 18.48`), adversarial 5-shell `2000:1` extreme skew `[10, 1/100, 7/10, 1/200, 3/10]` (`57991343772048/1163302171375`, `M_ADM = +58.22 > 0`, `max_H = 11.64`), and single-shell `M = 1` `[5]` (`14464/639`).
  * `battle3_readout` at `h = 0` safely returns `B_rho = None` with `first_law_residual == 0/1` (`ZeroDivisionError` eliminated).
- **`C2-W3`, 3+1D Curvature, Global 4-Spinor Hopf `S^2`, Separatrices, Palindromes & `2^-30` Tangents (`PASS`):**
  * Independent `QRat` + Sturm root-isolation audit over `h in {1/2, 1/4, 1/8, 1/16}` (`r_s = 2`) certifies zero poles at `r = 0`, exact match `K_engine == K_independent` (`0/1`), exact `0/1`-slack Kretschmann ceiling `B_K(2, h) = 24 r_s^2 / h^6`, and exact SOS polynomial coefficients `(14, 12, 44, 29, 12, 2)` in both algebra and `architect/winning_architecture.md`.
  * Global 4-spinor Hopf bundle (`ks_spinor_null_covector`, `ks_spinor_metric_tensors`) passes `7/7` pole-free checks across all `S^2` poles/antipodes at `h = 0` + exact origin at `h = 1/4`, and Cycle 2's `(-3r/5, -4r/5, 0)` ray is now regular at `h = 0` in `ks_null_covector`; equatorial `l_z == 0/1`.
  * Photon-sphere separatrices (`B^2 in {28, 109/4, 433/16}`, coarse `dt = 1/2` vs `1/8` exit-angle diffs `0.186 deg, 0.185 deg, 1.536 deg < 2.0 deg` + subcritical capture), `6/6` 10-step `0/1` palindromes (incl. `B1_3d_inclined` and `B2_3d_inclined_plunge`), dissipative `M = 16, nu = 1/2, p = 137` thermodynamics (`init_H = 0.712 < 1 -> peak_H = 1.925 >= 1`, `Q_rad = 2.369 > 0`), and `6/6` unfrozen discrete tangents down to `dt = 2^-30` at both `t = 0` and `t > 0` (`p = 101`) all `PASS`.

#### 4. Gate 4 — Zero-Handicap Baselines, `h = 0` & `h > 0` $O(dt^2)$ Slopes, Plots & `/tmp` `-I` Clean-Room Isolation (`PASS`)
- **Zero-Handicap Continuum `float64` Baselines (`sandbag_flags == []`):**
  * `battle1_exterior_orbit`: Euler `0.9940`, Verlet/RK2 `2.0000`, RK4 `3.9996`.
  * `battle2_pre_horizon_plunge`: Euler `0.9893`, Midpoint/RK2 `1.9999`, RK4 `4.0007`.
  * `battle3_pre_collapse_shells`: Euler `0.9938`, Verlet/RK2 `2.0000`, RK4 `3.9956`.
- **`Q`-Engine $O(dt^2)$ Convergence at Both `h = 0` AND `h > 0` (`C2-W1` Remediated):**
  * `battle1_exterior_h0_slope`: **`1.99999`**
  * `battle2_exterior_h0_slope`: **`1.99989`**
  * `battle2_exterior_h1_4_slope` (independent `RK4(1/2048)` check at `h = 1/4, T = 4.0`): **`2.13586`** (`errors = [7.749e-4, 1.513e-4, 3.620e-5, 9.000e-6]`, remediating Cycle 2's `1.4984` slope via `O(dt^4)` localized core regularization).
  * `battle2_plunge_bounce_h1_2_slope` (independent `RK4(1/2048)` horizon-crossing + de Sitter core-bounce check at `h = 1/2, T = 26.0`): **`2.01485`** (`errors = [5.154e-2, 3.413e-3, 2.033e-3, 5.825e-4]`).
  * `battle3_exterior_h0_slope`: **`1.99974`**
- **Publication Plots & Formalist Machine Certificates:**
  * Workspace `plots/*.{svg,png}`: **`8/8` PASS** (`1280x860` valid binary PNG `IHDR`/`IDAT`/`IEND` + non-empty geometric SVG XML).
  * `formalist/cycle3_proof_results.json` & `formalist/final_symbolic_proofs.md`: **`ALL_6_THEOREMS_AND_4_SUITES_MACHINE_CERTIFIED_0_1`** (`residual_findings == []`).
- **Exhibit 4 Isolated `/tmp` Clean-Room Execution (`/usr/bin/python3 -I standalone_ultrafinitist_gravity.py`):**
  * Copied `standalone_ultrafinitist_gravity.py` **alone** into a fresh empty `/tmp/cleanroom_exhibit4_*` directory and executed `/usr/bin/python3 -I standalone_ultrafinitist_gravity.py`:
    - **`returncode == 0`**, **`stderr == ""`**, wall-clock **`157.50s < 175s`**, Part I AST violations **`0`**.
    - All **`8/8` clean-room `./plots/*.{svg,png}`** generated self-contained via pure Python stdlib and verified:
      * `battle1_perihelion_and_photon_sphere`: SVG (`53,425` chars), PNG (`1280x860`, `110,421` bytes) — **PASS**
      * `battle2_horizon_plunge_and_kretschmann_ceiling`: SVG (`82,636` chars), PNG (`1280x860`, `91,672` bytes) — **PASS**
      * `battle3_oppenheimer_snyder_collapse_bounce`: SVG (`55,019` chars), PNG (`1280x860`, `107,702` bytes) — **PASS**
      * `convergence_and_no_singularity_summary`: SVG (`57,734` chars), PNG (`1280x860`, `72,527` bytes) — **PASS**

---

### Unredacted Non-Blocking Algebraic Fiber & Bandwidth Notes (`inquisitor/cycle3_audit_report.json`)

1. **`Lz != 0` Centrifugal Fiber Boundary at `s_ray == 0` (`lz_nonzero_origin_centrifugal_fiber_boundary_note == True`):** On states with `Lz != 0`, exact conservation `||x x u||^2 == Lz^2 != 0` excludes `(x, y, z) == (0, 0, 0)` (`s_ray == 0`, where `0 x u == 0`), so `_ks_ray_readout` naturally lives on `s_ray in Q \ {0}` for `Lz != 0` and on all `s_ray in Q` (including `s_ray == 0`) for `Lz == 0`, while `_ks_ray_step` is globally pole-free for all `h > 0, dt != 0`.
2. **Codimension-$M$ Stereographic Basepoint in `battle3_readout` (`battle3_vieta_adm_basepoint_vM_neg2_zero_div_note == True`):** At the superluminal stereographic north-pole basepoint `v_1 = ... = v_{M-1} = 0, v_M = -2, Q_rad = 0`, `battle3_readout`'s Vieta rational parameterization hits `0/0`; physical sub-luminal shell evolutions (`|v_k| < 1`) remain bounded strictly away from this codimension-$M$ locus.
3. **Ultra-Relativistic ($E_0 \ge 2.5$) `[1/1]` Padé Readout Scale Note (`w_1_5_E2p575_pade_readout_ur_bandwidth_note = -1.7395`):** On an extreme synthetic `w_doppler = 1/5` ($E_0 = 2.575$) origin transit (`h = 1/8, dt = 1/16, T = 12`), the dynamical signed-ray trajectory transits smoothly (`s_ray = -29.814, r = 25.292, dynamical ur_sh = +2.808 > 0, y == z == 0/1, ut = +2.008 > 0, C_metric == 0/1`), while `(E_0 - u^r_out)/d_0 ~ 1/(4 E_0^2) = 0.038` falls below the fixed `[1/1]` Padé scale `c_ref = (3/5) d0_ref` in `_causal_sos_w_eff`, mapping the projected mass-shell `ro["ur"]` negative at ultra-relativistic $E_0 \ge 2.5$ while all benchmark and moderate-relativistic transits ($E_0 \le 1.625$) satisfy `ro["ur"] > 0`.
