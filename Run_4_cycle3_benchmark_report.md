# Cycle 3 Final Benchmark & Adversarial Verification Report

**Author:** `@simulator` (Exact `Fraction` & Continuum Baseline Engine)  
**Target Engine:** [architect/cycle3_engine.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/architect/cycle3_engine.py) + Root [standalone_ultrafinitist_gravity.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/standalone_ultrafinitist_gravity.py)  
**Benchmark Suite & Raw Data:** [simulator/run_cycle3.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/simulator/run_cycle3.py) $\rightarrow$ [simulator/cycle3_results.json](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/simulator/cycle3_results.json), [simulator/continuum_baselines.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/simulator/continuum_baselines.py) $\rightarrow$ [simulator/continuum_baselines.json](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/simulator/continuum_baselines.json)  
**Publication Figures:** [plots/battle1_perihelion_and_photon_sphere.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle1_perihelion_and_photon_sphere.png), [plots/battle2_horizon_plunge_and_kretschmann_ceiling.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle2_horizon_plunge_and_kretschmann_ceiling.png), [plots/battle3_oppenheimer_snyder_collapse_bounce.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle3_oppenheimer_snyder_collapse_bounce.png), [plots/convergence_and_no_singularity_summary.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/convergence_and_no_singularity_summary.png) (plus matching `.svg` vector pairs)

---

## 1. Executive Summary & Cycle 1 $\rightarrow$ Cycle 2 $\rightarrow$ Cycle 3 Progression

`@architect`'s Cycle 3 winning architecture (`cycle3_engine.py`) replaces Cycle 2's hyperbolic shear + `beta2 * dt^2` core-inflated chart with a **Unified 3+1D 7-Slot Signed-Ray + $SO^+(1,1)$ Lorentz-Boost Envelope + Gauge-Fixed $SO(2,\mathbb{Q})$ Rotor + Strictly Causal SOS Readout** for Battles 1–2, Global Pole-Free 4-Spinor Hopf Kerr–Schild Tensors, and Mass-Proportional SOS Gap Ceilings for Battle 3. Every single finding across Cycle 1 (`R1–R5`, `W1–W2`, `F1–F4`) and Cycle 2 (`C2-R1..C2-R3`, `C2-W1..C2-W3`, `C2_F1..C2_F3`) is **100% resolved and empirically verified**.

| Audit / Physics Metric | Cycle 1 `Q`-Engine | Cycle 2 `Q`-Engine | Cycle 3 `Q`-Engine (Final) | Zero-Handicap `float64` Continuum |
| :--- | :--- | :--- | :--- | :--- |
| **Gate 1 Part I AST Purity** | `quantize` branch violations | `0` violations | **`0` violations (`quantize_branch_sites == []`)** | N/A (`math.sqrt`, IEEE-754) |
| **Gate 4 Smooth Convergence** | $h=0$: `2.00` | $h=0$: `2.00`; $h>0$: `O(1)` gap | **$h=0$: `2.0000`; $h=1/4$ ext: `2.1359`; $h=1/2$ bounce: `2.0148`** | Euler `0.99`, Verlet/RK2 `2.00`, RK4 `4.00` |
| **B1 Case 1A Precession & Drift** | Exact `0/1` invariants | `151` bits, `pos_err=2.2e-3` (`dt=1/4`) | **`69` bits, `pos_err=2.01e-4` (`dt=1/4`), `0/1` invariants** | Euler `3.05e-2`, RK2 `2.65e-6`, RK4 `2.04e-10` $\|dE\|$ |
| **B1 Case 1B Grazing (`dt=1/2`)** | `B^2` offset (`27.06` vs `27.25`) | `< 1.66 deg` exit err | **`< 1.44 deg` (`0.088 deg` on `B^2=109/4`, `E_0 == 1/1`)** | Schw RK2 **False Plunge!**; KS RK2 `-118.4 deg` err |
| **B2 Plunge Bounce Apocenter** | Catapulted (`r_max = 24,267`) | Bounded (`8.5..9.4`), `E_kill` spikes | **Bounded (`8.53..9.31`), tight `\|dE_kill\| <= 0.64`** | Schw blows up at $r \to 2$; KS crashes at $r \to 0$ |
| **B2 12-Run Origin Transit (`T=20`)** | Reflected (`x >= 0`) | `11/12` (`h=1/8, dt=1/4` reflected; catapult `-92,692`) | **`12/12 PASS` (`s_ray_final in [-17.93, -15.05]`, `y=z=0/1`)** | `ZeroDivisionError` / `K -> Inf` at Step 11–22 |
| **B2 Midflight Exact `r == 0/1`** | Not attained (`r > 0`) | `r_min > 0` | **Bitwise `r == 0/1, s_ray == 0/1, (x,y,z) == (0,0,0)` at Step 20** | Exact `ZeroDivisionError` at $r=0$ |
| **B2 Causal Orientation (`ut > 0`)** | Acausal core flips | Acausal Vieta flips inside $H_h > 1$ | **`0` acausal flips (`min ut >= 1.0039 > 1` everywhere!)** | Coordinate breakdown ($g_{tt}=0$) / singularity |
| **B2 Global $S^2$ Null Covector** | Equatorial `l_z != 0` | Cartesian single-chart antipode at $h=0$ | **4-Spinor Hopf `ks_spinor_null_covector` (`0` poles on 100% of $S^2$)** |Undefined at $r=0$ (`0/0`) |
| **B3 Density Ceiling $B_\rho$** | Violated on bound cloud | Violated `10.18x` on skewed `[9, 1]` | **`0/1` slack across ALL uniform, skewed `[9/2,1/2]`, `[9,1]` & bound clouds** | Crashes Step 11–25 ($\rho \to \infty$, shell crossing) |
| **B3 21-Prime Drift (`M=8, 16`)** | `18` primes (`p=113` drift in B2) | `21/21` (`0.0000 bits/step`) | **`21/21` (`M=8`: `193..251` bits; `M=16`: `314..430` bits; `0.0000 bits/step`)** | **21/21 crash at Step 13 (`M=8`) / Step 12 (`M=16`)** |
| **Clean-Room Standalone (`-I`)** | Cycle 3 deliverable | Cycle 3 deliverable | **`PASS` (`EXIT 0`, `stderr == ""`, `8/8` valid SVG+PNG plots)** | N/A |

---

## 2. Gate 4 Verification: Zero-Handicap Continuum Slopes & True $h \ge 0$ `Q`-Engine $O(\Delta t^2)$ Convergence

All `float64` continuum baselines in [simulator/continuum_baselines.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/simulator/continuum_baselines.py) are textbook-exact (zero jitter, zero `random`, zero precision handicap; `sandbag_flags == []`), achieving exact theoretical convergence orders `(1.0, 2.0, 4.0)` against `RK4(dt = 1/2048)` / `RK4(dt = 1/4096)`. Crucially, replacing Cycle 2's $O(\Delta t^2)$ core inflation with localized $O(\Delta t^4)$ core regularization (`dt4 = dt^4 / (1 + 24*dt^4)`) restores **true second-order $O(\Delta t^2)$ convergence to the Hayward continuum ODE for $h > 0$** both in the exterior and right through horizon traversal + de Sitter core bounce!

| Smooth / Regularized Benchmark Regime | Duration $T$ | Step Ladder $\Delta t$ | `f64` Euler | `f64` Verlet | `f64` Mid-RK2 | `f64` RK4 | Cycle 3 `Q`-Engine Slope |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Battle 1 Strong-Field Exterior ($h = 0$)** | $16.0$ | $\{1/2 \dots 1/32\}$ | `0.9940` | `2.0000` | `2.0049` | `3.9996` | **`2.0000`** (`1.999989`) |
| **Battle 2 Pre-Horizon Exterior ($h = 0$)** | $4.0$ | $\{1/4 \dots 1/32\}$ | `0.9893` | `1.9999` | `1.9896` | `4.0007` | **`1.9999`** (`1.999885`) |
| **Battle 2 Hayward Exterior ($h = 1/4 > 0$)** | $4.0$ | $\{1/4 \dots 1/32\}$ | — | — | — | Ref `RK4(1/2048)` | **`2.1359`** (`8.98e-6` at `dt=1/32`) |
| **Battle 2 Full Horizon + Core Bounce ($h = 1/2 > 0$)** | $26.0$ | $\{1/8 \dots 1/64\}$ | — | — | — | Ref `RK4(1/4096)` | **`2.0148`** (`5.83e-4` post-bounce!) |
| **Battle 3 Multi-Shell LTB Pre-Collapse ($h = 0$)** | $0.5$ | $\{1/16 \dots 1/128\}$ | `0.9938` | `2.0000` | `1.9933` | `3.9956` | **`1.9997`** (`1.999745`) |

---

## 3. Battle 1: Relativistic Geodesics, 3+1D Inclined Precession & Photon-Sphere Separatrix (`M = 1, r_s = 2, h = 0`)

### 3.1 Case 1A Eccentric Strong-Field Bound Orbit ($r_{\text{peri}} \approx 8M, e \approx 1/3$, $n = 1,000$ Steps)

| Integrator | $\Delta t$ | $\|\Delta E\|_{\max}$ | $\|\Delta L_z\|_{\max}$ | $\|C_{\text{metric}}\|_{\max}$ | Pos Err vs `RK4(1/512)` | Orbit Precession ($\Delta\phi_{\text{exact}} = 177.32^\circ$) | Peak Bits / Late Drift |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `f64` Forward Euler | $1/4$ | `6.17e-3` | `1.04e-1` | `1.19e-2` | `15.52` (spirals out) | Escapes bound shell | IEEE-754 (53-bit) |
| `f64` Midpoint RK2 | $1/4$ | `2.65e-6` | `3.21e-5` | `5.11e-6` | `1.54e-2` | $186.41^\circ$ ($r_p=8$ IC) | IEEE-754 (53-bit) |
| `f64` Symplectic Verlet | $1/4$ | `2.05e-6` | `8.88e-16` | `3.96e-6` | `6.98e-3` | $186.35^\circ$ ($r_p=8$ IC) | IEEE-754 (53-bit) |
| `f64` Classical RK4 | $1/4$ | `2.04e-10` | `3.42e-9` | `3.93e-10` | `1.12e-6` | $186.38^\circ$ ($r_p=8$ IC) | IEEE-754 (53-bit) |
| **Cycle 3 `Q`-Engine** | **$1/2$** | **`0/1`** | **`0/1`** | **`0/1`** | **`3.89e-3`** | **$177.46^\circ$ (`0.14 deg` err)** | **`71` bits / `+0.0006` b/st** |
| **Cycle 3 `Q`-Engine** | **$1/4$** | **`0/1`** | **`0/1`** | **`0/1`** | **`2.01e-4`** | **$177.46^\circ$ (`0.14 deg` err)** | **`69` bits / `+0.0028` b/st** |
| **Cycle 3 `Q`-Engine** | **$1/8$** | **`0/1`** | **`0/1`** | **`0/1`** | **`8.22e-5`** | **$177.46^\circ$ (`0.14 deg` err)** | **`70` bits / `0.0000` b/st** |

- **3+1D Off-Equatorial Orbit (`1A_3d_inclined`, quaternion frame `(2,1,2,0)`, `400` steps at $\Delta t = 1/4$):** Genuine 3D excursion $\max |z| = 11.1969$, $\max |u^z| = 0.4563$, exact $0/1$ conservation of $\|\mathbf{x} \times \mathbf{u}\|^2 - L_z^2 == 0/1$, $x^2 + y^2 + z^2 - r^2 == 0/1$, $C_{\text{metric}} == 0/1$, exact 3D Wildberger spread $s_{3\text{D}} \in \mathbb{Q}$ (`0.348590`), `73` peak bits, and **exact `0/1` 10-step palindrome reversibility**.
- **Unit-Circle $S^1$ Rotor Gauge-Fixing:** Scaling `(ch, sh) -> (2*ch, 2*sh)` produces a distinct physical coordinate-time step (`ro_scaled["t"] != ro_unscaled["t"]`), breaking projective rotor gauge degeneracy while preserving `C^2 + S^2 == 1/1` and `0/1` palindrome reversibility on the unit circle.

### 3.2 Case 1B Photon-Sphere Separatrix Grazing ($B_{\text{crit}}^2 = 27 M^2$)

All three null presets (`E_0 == 1/1` exact, Weierstrass cubic chord `0/1` residual verified) and relativistic timelike grazing scatter cleanly outside $r_{\text{ph}} = 3M$ with exact $4\times$ quadratic error reduction per $\Delta t$ halving:

| Grazing Preset ($B^2 > 27$) | `f64` Schw Polar RK2 ($\Delta t=1/2$) | `f64` KS RK2 ($\Delta t=1/2$) | `Q`-Engine Exit Err ($\Delta t=1/2$) | `Q`-Engine Exit Err ($\Delta t=1/4$) | `Q`-Engine Exit Err ($\Delta t=1/8$) | `Q` Invariants |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `1B_null_B2_28` ($B^2 = 28$) | **FALSE PLUNGE** ($r \to r_s$) | `-45.1 deg` err | **`0.278 deg`** (`SCATTER`) | **`0.070 deg`** (`SCATTER`) | **`0.018 deg`** (`SCATTER`) | **`0/1` Exact** |
| `1B_null_B2_109_4` ($B^2 = 27.25$) | **FALSE PLUNGE** ($r \to r_s$) | `-118.4 deg` err | **`0.088 deg`** (`SCATTER`) | **`0.020 deg`** (`SCATTER`) | **`0.005 deg`** (`SCATTER`) | **`0/1` Exact** |
| `1B_null_B2_433_16` ($B^2 = 27.0625$) | **FALSE PLUNGE** ($r \to r_s$) | **FALSE PLUNGE** | **`1.438 deg`** (`SCATTER`) | **`0.347 deg`** (`SCATTER`) | **`0.086 deg`** (`SCATTER`) | **`0/1` Exact** |
| `1B_timelike_scattering` ($L_z = 10$) | Phase drift | `-19.8 deg` err | **`1.001 deg`** (`SCATTER`) | **`0.256 deg`** (`SCATTER`) | **`0.064 deg`** (`SCATTER`) | **`0/1` Exact** |

---

## 4. Battle 2: Horizon-Penetrating Plunge, 12/12 Origin Transit, Bitwise `0/1` Origin Hit & Causal SOS Readout

### 4.1 Global 4-Spinor Hopf Chart & Kretschmann Ceiling $B_K(r_s, h) = 24 r_s^2 / h^6$
- **Global $S^2$ Coverage (`ks_spinor_null_covector` & `ks_spinor_metric_tensors`):** Verified `0` poles across 100% of $S^2$ even at $h = 0$ — including `(-1, 0, 0)`, `(-3/5, -4/5, 0)`, `(-88/137, -105/137, 0)`, north/south poles `(0, 0, +-1)`, and exact origin `(0, 0, 0, 0)` for $h > 0$ — with $\eta^{\mu\nu}\ell_\mu\ell_\nu == 0/1$, $g_{\mu\alpha}g^{\alpha\nu} == \delta_\mu^\nu$, and $x^2+y^2+z^2-r^2 == 0/1$.
- **Exact `0/1` Kretschmann Slack:** $K_h(0) == B_K(2, h)$ identically ($6,144$ for $h=1/2$; $393,216$ for $h=1/4$; $25,165,824$ for $h=1/8$), whereas unregularized continuum `float64` blows up at Step 19–22 ($K_0(r) = 48/r^6 \to \infty$).

### 4.2 Plunge Matrix (`n = 500` Steps): Zero Catapulting, Zero Acausal Flips & Tight Killing Energy

| Preset & Scale $h$ | $\Delta t$ | $\max H_h(r)$ ($\ge 1$) | $r_{\min}$ | Post-Bounce $r_{\max}$ | $\min u^t$ (`> 1`) | Acausal Steps | $\max \|E_{\text{kill}} - E_0\|$ | Peak Bits / Drift |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `angular_plunge_Lz2` ($h=1/2$) | $1/4$ | `1.660` | `1.102` | **`9.304`** | **`1.1439`** | **`0`** | `0.492` | `70` / `+0.0006` |
| `angular_plunge_Lz2` ($h=1/4$) | $1/4$ | `1.884` | `1.047` | **`9.003`** | **`1.1440`** | **`0`** | `0.535` | `73` / `+0.0002` |
| `angular_plunge_Lz2` ($h=1/8$) | $1/4$ | `2.073` | `0.963` | **`9.303`** | **`1.1440`** | **`0`** | `0.630` | `76` / `+0.0009` |
| `head_on_plunge_Lz0` ($h=1/2$) | $1/4$ | `2.111` | `0.238` | **`8.534`** | **`1.1250`** | **`0`** | `0.548` | `70` / `-0.0020` |
| `head_on_plunge_Lz0` ($h=1/4$) | $1/4$ | `4.227` | `0.215` | **`8.533`** | **`1.1250`** | **`0`** | `0.537` | `73` / `-0.0027` |
| `head_on_plunge_Lz0` ($h=1/8$) | $1/4$ | `7.934` | `0.206` | **`8.533`** | **`1.1250`** | **`0`** | `0.536` | `76` / `-0.0005` |
| `3d_inclined_plunge` ($h=1/4$, $\max|z|=7.90$) | $1/4$ | `1.884` | `1.047` | **`9.003`** | **`1.1440`** | **`0`** | `0.535` | `73` / `+0.0002` |

### 4.3 Complete 12-Run `energetic_origin_transit` Matrix ($T = 20$, $s_{\text{ray}}: +8 \to 0 \to < 0$)

Every single $(h, \Delta t) \in \{1/2, 1/4, 1/8\} \times \{1/2, 1/4, 1/8, 1/16\}$ combination transits smoothly through $s_{\text{ray}} = 0$ to the antipodal ray $x_{\text{unfolded}} < 0$ with **zero reflection**, **zero catapulting** ($r_{\text{final}} \in [15.07, 17.98] < 20$), **exact axial purity** ($y == y_{\text{unfolded}} == 0/1$, $z == z_{\text{unfolded}} == 0/1$), and **zero acausal flips**:

| Scale $h$ | $\Delta t$ | Steps ($T=20$) | $s_{\text{ray,final}}$ ($= x_{\text{unfolded}}$) | Folded $r_{\min}$ | Folded $r_{\text{final}}$ ($< 20$) | Folded $u^r_{\text{final}}$ ($> 0$) | $y = z = 0/1$ | $\min u^t$ | Peak Bits | Verdict |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $1/2$ | $1/2$ | 40 | `-15.0548` | `0.1361` | `15.0691` | `+0.6031` | `True` | `1.3190` | `72` | **PASS** |
| $1/2$ | $1/4$ | 80 | `-15.2700` | `0.0513` | `15.2719` | `+0.6062` | `True` | `1.2522` | `73` | **PASS** |
| $1/2$ | $1/8$ | 160 | `-15.5768` | `0.0326` | `15.5780` | `+0.6097` | `True` | `1.0355` | `74` | **PASS** |
| $1/2$ | $1/16$ | 320 | `-15.6495` | `0.0167` | `15.6499` | `+0.6105` | `True` | `1.0273` | `75` | **PASS** |
| $1/4$ | $1/2$ | 40 | `-15.0776` | `0.1334` | `15.1086` | `+0.6030` | `True` | `1.3171` | `75` | **PASS** |
| $1/4$ | $1/4$ | 80 | `-15.3161` | `0.0340` | `15.3157` | `+0.6071` | `True` | `1.2548` | `76` | **PASS** |
| $1/4$ | $1/8$ | 160 | `-15.7499` | `0.0013` | `15.7520` | `+0.6115` | `True` | `1.0275` | `77` | **PASS** |
| $1/4$ | $1/16$ | 320 | `-15.9099` | `0.0177` | `15.9099` | `+0.6138` | `True` | `1.0920` | `78` | **PASS** |
| $1/8$ | $1/2$ | 40 | `-15.0950` | `0.1297` | `15.1536` | `+0.6032` | `True` | `1.3145` | `78` | **PASS** |
| $1/8$ | $1/4$ | 80 | `-15.3539` | `0.0223` | `15.3299` | `+0.6088` | `True` | `1.3150` | `79` | **PASS** |
| $1/8$ | $1/8$ | 160 | `-15.2596` | `0.0130` | `15.2708` | `+0.6079` | `True` | `1.0900` | `80` | **PASS** |
| $1/8$ | $1/16$ | 320 | `-17.9340` | `0.0005` | `17.9764` | `+0.3003` | `True` | `1.0760` | `82` | **PASS** |

### 4.4 Bitwise-Exact Mid-Flight Origin Transit (`exact_midflight_origin_transit`, $\Delta t = 1/4$)
Across all $h \in \{1/2, 1/4, 1/8\}$, exterior infall ($r_0 = 7.51, 7.66, 7.88 > 2 r_s$) crosses the horizon ($\max H_h = 2.12, 4.16, 6.52 \ge 1$), lands on **bitwise-exact `r == 0/1, s_ray == 0/1, (x, y, z) == (0/1, 0/1, 0/1)` at Step 20** with zero `0/0` division and `C_metric == 0/1`, emerges cleanly to $s_{\text{ray,step40}} = -7.50, -7.67, -7.89$ ($r_{\text{step40}} == r_0$), and passes **25-step forward/backward palindrome reversibility across `r == 0/1` with `0/1` state residual**.

### 4.5 Battle 2 21-Prime Adversarial Suite (`101..199`, Including `p = 113`)
Across all 21 primes ($n = 250$, plus $n = 500$ on $p \in \{101, 113, 199\}$): peak state bit-lengths stay strictly within **`120 .. 136` bits** (`<< 1,500` target) with late-time drift slopes within **`[-0.0119, +0.0225]` bits/step** (`p = 113` at 500 steps: `132` peak bits, `+0.0095` bits/step $\ll 0.05$).

---

## 5. Battle 3: Supercritical Oppenheimer–Snyder / LTB Collapse, Universal Non-Uniform $B_\rho$ & 21 Primes

### 5.1 Universal `0/1`-Slack Density Ceiling $B_\rho(M_{\text{rest}}, h, M, \mathbf{m})$ & Safe $h = 0$ Readout
With mass-proportional inner gap floors $\delta_1 = \frac{3}{4}h$ and $\delta_k = \frac{m_k}{4 M_{\text{rest}}}h$ ($k \ge 2$), Shell 2 at pinch $\mathbf{q} = \mathbf{0}$ attains $B_\rho$ with **exact `0/1` slack** and **dominates every shell $k \in \{1 \dots M\}$ for ALL positive mass vectors $\mathbf{m} \in \mathbb{Q}_{>0}^M$**:

| Configuration | Mass Vector $\mathbf{m}$ | Exact $B_\rho$ | $\sup_{\mathbf{q}} \max_k \rho_k(\mathbf{q})$ | Config Slack | Live Evolved $\max_k \rho_k$ | Verdict |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Skewed Static (`C2_F2`)** | $[9/2, 1/2]$ | `4339200/198161` (`21.8973`) | `4339200/198161` | **`0/1` Exact** | `6.8902 <= 21.8973` (100 st) | **PASS** |
| **Skewed Infall (`C2-R3`)** | $[9, 1]$, $V_0=[-3,-3]$ | `8678400/198161` (`43.7947`) | `8678400/198161` | **`0/1` Exact** | `39.8985 <= 43.7947` ($\max H=18.42$) | **PASS** |
| **Bound Cloud ($M_{\text{ADM}} = \frac{1415}{1736} < 5$)** | $[5/2, 5/2]$ | `173568/9017` (`19.2490`) | `173568/9017` | **`0/1` Exact** | `5.3307 <= 19.2490` (100 st) | **PASS** |
| **Uniform 8-Shell ($M=8$)** | $[5/4 \dots 5/4]$ | `5554176/127871` (`43.4358`) | `5554176/127871` | **`0/1` Exact** | `34.2979 <= 43.4358` (300 st) | **PASS** |
| **Uniform 16-Shell ($M=16$)** | $[5/8 \dots 5/8]$ | `22216704/501047` (`44.3406`) | `22216704/501047` | **`0/1` Exact** | `2.8255 <= 44.3406` (300 st) | **PASS** |
| **Unregularized Exterior ($h=0$)** | $[5/4 \dots 5/4]$ | `B_rho is None` (safe!) | Continuum limit | Safe (`0` crash) | `first_law == 0/1`, `C_k == 0/1` | **PASS** |

### 5.2 Clean & 21-Prime (`101..199`) Supercritical Collapse (`n = 300` Steps, `M = 8` & `M = 16`)

| Ensemble (`amp = 10, h = 1/2`) | `f64` Continuum Failure Step | `Q`-Engine $\max_k H_k$ ($\ge 1$) | Ordering $R_1 > 0, \Delta R_k > 0$ | $\max_k \rho_k \le B_\rho$ | 1st Law & $t = n\Delta t$ | 2nd Law $\Delta Q_{\text{rad}} \ge 0$ | Peak Bits | Late Drift (bits/step) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `M=8` Clean ($\nu=0, \Delta t=1/4$) | Step 13 (Euler 11, Verlet 13, RK4 13) | `5.959` | **300/300 OK** | `34.30 <= 43.44` | **`0/1` Exact** | `Q_rad == 0/1` (Rev `0/1`) | **`84`** | **`0.0000`** |
| `M=8` Clean ($\nu=1/2, \Delta t=1/4$) | Step 13 (`R_1 <= 0` / shell cross) | `7.926` | **300/300 OK** | `27.37 <= 43.44` | **`0/1` Exact** | `Q_rad -> 8.00` (`>= 0/1`) | **`84`** | **`0.0000`** |
| `M=16` Clean ($\nu=0, \Delta t=1/4$) | Step 12 (Euler 11, Verlet 12, RK4 12) | `1.322` | **300/300 OK** | `2.80 <= 44.34` | **`0/1` Exact** | `Q_rad == 0/1` (Rev `0/1`) | **`132`** | **`0.0000`** |
| `M=16` Clean ($\nu=1/2, \Delta t=1/4$) | Step 12 (`R_1 <= 0` / shell cross) | `2.514` | **300/300 OK** | `2.61 <= 44.34` | **`0/1` Exact** | `Q_rad -> 6.00` (`>= 0/1`) | **`132`** | **`0.0000`** |
| **`M=8` All 21 Primes (`101..199`)** | **21/21 crash at Step 13** | **`>= 1` on 21/21** | **21/21 OK** | **21/21 OK** | **`0/1` on 21/21** | **`>= 0/1` on 21/21** | **`193 .. 251`** | **`0.0000` on all 21!** |
| **`M=16` All 21 Primes (`101..199`)** | **21/21 crash at Step 12** | **`>= 1` on 21/21** | **21/21 OK** | **21/21 OK** | **`0/1` on 21/21** | **`>= 0/1` on 21/21** | **`314 .. 430`** | **`0.0000` on all 21!** |

### 5.3 Per-Step Latency & SIGALRM Compliance
Every single `Q`-engine step across all Three Battles ran under `signal.SIGALRM` (`5.0s` hard ceiling) with **`0` timeouts**:
- **Battle 1 max single-step wall-clock:** `0.0099 s` (`9.9 ms`)
- **Battle 2 max single-step wall-clock:** `0.0086 s` (`8.6 ms`)
- **Battle 3 max single-step wall-clock (`M=16` 21-prime):** `0.0174 s` (`17.4 ms`)

---

## 6. Clean-Room Standalone Exhibit 4 & Dual SVG+PNG Deliverables

1. **Self-Contained Clean-Room Artifact:** [standalone_ultrafinitist_gravity.py](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/standalone_ultrafinitist_gravity.py)
   - **Part I (above `# === PART II: ...`):** Verbatim lines 1–2064 of `architect/cycle3_engine.py` (`audit_ast_purity(..., part1_only=True)` $\rightarrow$ `passed == True`, `total_violations == 0`, `quantize_branch_sites == []`).
   - **Part II–IV (below `# === PART II: ...`):** Zero-handicap `float64` continuum baselines + pure-stdlib `zlib`/`struct` dual `.svg`/`.png` plotter + self-contained `_QQPoly`/`_QUniRat`/`_QQJet` symbolic identity verifier + Three-Battle `SIGALRM`-guarded benchmark suite.
   - **Clean-Room Execution (`inquisitor.gate_auditor.audit_cleanroom_standalone`):** Executed via `/usr/bin/python3 -I standalone_ultrafinitist_gravity.py` inside a fresh empty `/tmp/cleanroom_exhibit4_*` directory with **`EXIT 0`**, **`stderr == ""`**, `< 175s` wall-clock, and all **8** `./plots/*.{svg,png}` files generated and verified.
2. **Publication Figure Pairs (`8/8` Validated):**
   - [battle1_perihelion_and_photon_sphere.svg](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle1_perihelion_and_photon_sphere.svg) & [battle1_perihelion_and_photon_sphere.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle1_perihelion_and_photon_sphere.png)
   - [battle2_horizon_plunge_and_kretschmann_ceiling.svg](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle2_horizon_plunge_and_kretschmann_ceiling.svg) & [battle2_horizon_plunge_and_kretschmann_ceiling.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle2_horizon_plunge_and_kretschmann_ceiling.png)
   - [battle3_oppenheimer_snyder_collapse_bounce.svg](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle3_oppenheimer_snyder_collapse_bounce.svg) & [battle3_oppenheimer_snyder_collapse_bounce.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/battle3_oppenheimer_snyder_collapse_bounce.png)
   - [convergence_and_no_singularity_summary.svg](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/convergence_and_no_singularity_summary.svg) & [convergence_and_no_singularity_summary.png](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-gravity/plots/convergence_and_no_singularity_summary.png)
