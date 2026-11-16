# Cycle 3 Empirical Benchmark, Zero-Handicap Continuum Comparison & Clean-Room Verification Report (`@simulator`)

**Author:** `@simulator` (Exact `Fraction` & Zero-Handicap Continuum Engine)  
**Deliverables Covered:**
- `simulator/continuum_baselines.py` (Honest Full-Precision IEEE-754 `float64` Continuum Baselines; `sandbag_flags == []`)
- `simulator/plotter.py` & `plots/*.{svg,png}` (Pure-Stdlib Dual Vector `.svg` + `zlib`/`struct` 24-bit RGB `.png` Scientific Renderer; `8/8` Plots)
- `standalone_ultrafinitist_statmech.py` (Self-Contained Clean-Room Artifact: `5,164`-line Pure-$\mathbb{Q}$ Part I + `2,769`-line Benchmark/Plotter Part II; `7,933` lines total)
- `checkpoints/simulator_cycle{1,2,3}.json` (3-Cycle Empirical Telemetry Ledger)

---

## 1. Executive Summary & Adversarial Zero-Handicap Protocol

All empirical benchmarks in `simulator/` adhere strictly to two non-negotiable constitutional principles:

1. **Zero Sandbagging of Continuum Baselines (`sandbag_flags == []`):**  
   Every continuum comparator in `simulator/continuum_baselines.py` uses full-precision 64-bit IEEE-754 floating-point arithmetic (`float`, `math.exp`, `math.log`) and textbook numerical integration schemes (`Forward Euler`, symplectic `Velocity Verlet` / Störmer-Verlet, `Midpoint RK2`, and `Classical RK4`) with zero artificial precision truncation (`float32`/`round`), zero injected noise, and exact independent agreement (`< 10^{-15}`) with `@inquisitor`'s Gate 4 reference integrators.
2. **Strict Separation of Pure-$\mathbb{Q}$ Engine and Continuum Diagnostics:**  
   In `standalone_ultrafinitist_statmech.py`, **Part I** (`5,164` lines) embeds `@architect`'s combinatorial, canonical, and birational Hamiltonian modules alongside `@formalist`'s computer-algebra and proof suites using **exclusively** `int` and `fractions.Fraction` (`0` AST violations under Gate 1's static purity scanner). **Part II** executes all continuum comparisons, quantitative assertions, and pure-stdlib `.svg`/`.png` plot generation inside an isolated `/tmp` clean-room via `/usr/bin/python3 -I` (`EXIT 0`, `stderr == ""`, `38.595s < 180s`).

### 1.1 Cross-Cycle Adversarial Benchmark Progression (`Cycle 1 -> Cycle 2 -> Cycle 3`)

| Benchmark Metric | Cycle 1 (`simulator_cycle1.json`) | Cycle 2 (`simulator_cycle2.json`) | Cycle 3 (`simulator_cycle3.json`) | Gate Requirement / Target |
| :--- | :---: | :---: | :---: | :---: |
| **Part I Pure-$\mathbb{Q}$ AST Violations** | `0` (`4,496` lines) | `0` (`5,173` lines) | **`0`** (`5,164` lines) | Hard `== 0` |
| **Continuum Baseline Slopes (`Euler, Verlet, RK2, RK4`)** | `1.03, 2.00, 2.08, 3.98` | `0.95, 2.00, 2.06, 4.03` | **`0.95, 2.00, 2.06, 4.03`** | `(1.0, 2.0, 2.0, 4.0) +- 0.15` |
| **Pure-$\mathbb{Q}$ Engine Smooth $O(\Delta t^2)$ Slope** | `0.0142` (`FAIL` `[C1-G4-001]`) | `1.9450` (Gate 4: `1.9613`) | **`1.9450`** (Gate 4: `1.9619`) | `[1.85, 2.15]`, `err(1/128) < 5e-3` |
| **Carrier Peak Bit-Length (`n = 300`, Dyadic + 21 Primes)** | `34` (`44,170` at $\Delta t = 1/4$) | `62` (`63` global across 47 runs) | **`62`** (`63` global across 47 runs) | Target `< 1,500`, Hard `<= 10,000` |
| **21-Prime Worst Late Entropy Drift (`n = 100..300`)** | `+0.0020 bits/step` | `0.000000 bits/step` | **`0.000000 bits/step`** | Hard `<= 0.05 bits/step` |
| **$N = 64$ Cooling Quench (`240` steps across $E_c = 160$)** | `div/0` @ step `210` (`[C1-G3-002]`) | `240/240` (`E: 228.20 -> 57.06`) | **`240/240`** (`E: 228.20 -> 57.06`) | `0` `div/0`, $E_{\text{end}} < 160$, $0/1$ 1st Law |
| **Discrete Tangent Error at $\Delta t = 2^{-30}$ (`t = 0` / `t > 0`)** | `4.66e-7` / `0.2197` (`[C1-G3-003]`) | `4.49e-9` / `4.37e-9` | **`4.49e-9` / `4.37e-9`** (`slope = 1.01`) | Unfrozen $O(\Delta t)$ down to $2^{-30}$ |
| **`/tmp` `python3 -I` Clean-Room Standalone Audit** | `EXIT 0` (`29.66s`, `8/8` plots) | `EXIT 0` (`38.55s`, `8/8` plots) | **`EXIT 0` (`38.60s`, `8/8` plots)** | `EXIT 0`, `stderr == ""`, `< 180s` |
| **Official `@inquisitor` Red-Team Verdict** | `FAIL` (`15` defects; `0` on `@simulator`) | `PASS` (`0` defects across Gates 1–4) | **Ready (`0` preflight defects)** | `0` defects across Gates 1–4 |

---

## 2. Battle 1: Exact Combinatorial Microcanonical vs. `float64` Canonical Ensemble

### 2.1 Exact Transfer-Frontier State Counting & $2^{16}$ Brute-Force Cross-Check
Across periodic $L \times L$ square tori for $L \in \{4, 6, 8\}$ ($N \in \{16, 36, 64\}$), both Case 1A (2nd-order ferromagnetic Ising torus, $J = 1$) and Case 1B ($C_4$-symmetric $2 \times 2$ plaquette-bond Ising torus, $J_1 = 2, J_2 = 1, J_4 = -1$) execute in **`0.9957s` combined wall-clock time** (`< 5.0s` `SIGALRM` guard):

| Lattice Size $L \times L$ ($N$) | Total Microstates $\sum_k \Omega_N(E_k) = 2^N$ | Case 1A Exact Inflection $\beta_c(N) \in \mathbb{Q}$ | Case 1B Backbending Pairs ($\Delta\beta_{\mathbb{Q}} > 0$) | Case 1B $\min_k (C_{\text{micro}}(E_k) / N)$ | Canonical $\min_\beta (C_{\text{can}}(\beta) / N)$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **$4 \times 4$ ($N = 16$)** | `65,536` (`max |dOmega| == 0` vs $2^{16}$ brute-force) | `3217 / 14224` (`0.226167`) | `3` shell pairs (`5` $\kappa_S > 0$ shells) | **`-2.17031`** | `> 0.0` (`C_can >= 0` strictly) |
| **$6 \times 6$ ($N = 36$)** | `68,719,476,736` ($2^{36}$) | `39034882895 / 170394708754` (`0.229085`) | `8` shell pairs (`7` $\kappa_S > 0$ shells) | **`-10.58211`** | `> 0.0` (`C_can >= 0` strictly) |
| **$8 \times 8$ ($N = 64$)** | `18,446,744,073,709,551,616` ($2^{64}$) | `17652068955825 / 75904121287972` (`0.232557`) | `11` shell pairs (`11` $\kappa_S > 0$ shells) | **`-0.09393`** | **`+0.000071`** ($\ge 0$ everywhere) |

### 2.2 Ensemble Inequivalence & Genuine IEEE-754 `OverflowError` Cliff
1. **Erasure of First-Order Phase Coexistence in the Canonical Ensemble:**  
   In Case 1B, interfacial domain-wall tension creates a certified finite-$N$ **convex entropy intruder** ($\kappa_S(E_k) > 0/1$), producing caloric S-curve backbending ($\beta_{\mathbb{Q}}(E_{k+1}) > \beta_{\mathbb{Q}}(E_k)$) and strictly negative microcanonical specific heat $C_{\text{micro}}(E_k) < 0/1$ inside the phase-coexistence band. By contrast, even with overflow-safe log-sum-exp shifting (`canonical_thermodynamics_shifted_float64`), the canonical specific heat $C_{\text{can}}(\beta) = \beta^2\,\mathrm{Var}_\beta(E) \ge 0$ is algebraically non-negative everywhere, replacing the physical S-curve across coexistence with a smooth monotonic Maxwell bridge that smears out latent-heat coexistence.
2. **Unshifted IEEE-754 `float64` Partition Sum Overflow (`-beta * E_min > 709.7827`):**  
   On $N = 64$ Case 1A ($E_{\min} = -128$), evaluating the textbook unshifted canonical partition sum $Z_{\text{can}}(\beta) = \sum_k \Omega_{64}(E_k)\,\exp(-\beta E_k)$ (`canonical_sum_unshifted_strict`) hits genuine IEEE-754 `OverflowError: math range error` for every inverse temperature above
   $$\beta_{\text{overflow}}(64) = \frac{\ln(\text{DBL\_MAX})}{|E_{\min}|} = \frac{709.782712893384}{128} = 5.54518$$
   (verified via `OverflowError` trap at $\beta = 6.0$), whereas the pure-$\mathbb{Q}$ Padé generating polynomial $Z_{64}(u)$ evaluates cleanly at arbitrary rational $\beta \in \mathbb{Q}_{>0}$ without overflow or underflow.

---

## 3. Battle 2: SOS `[2/2]` Padé Ising Convergence, Critical Ceilings & Lee-Yang Quadrance

### 3.1 Asymptotic Convergence of SOS `[2/2]` Padé $R_{2,2}(x)$ vs. `math.exp(-x)`
Using exact $N = 16$ Ising degeneracies ($\Delta E = 4$), `pade_boltzmann_convergence_table` benchmarks $R_{2,2}(x) = \frac{(x - 3)^2 + 3}{(x + 3)^2 + 3}$ against IEEE-754 `math.exp(-x)` over $x = \beta\Delta E \in \{2^{-2}, \dots, 2^{-7}\}$:

| $x = \beta\Delta E$ | $R_{2,2}(x)$ (`float64` view) | $\exp(-x)$ (`math.exp`) | Pointwise Error $|R_{2,2}(x) - e^{-x}|$ | Derivative Error $|R'_{2,2}(x) + e^{-x}|$ | Observable Error $|\langle E\rangle_{\text{Padé}} - \langle E\rangle_{\exp}|$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| `1/4` (`0.2500000`) | `0.778801843318` | `0.778800783071` | `1.060247e-06` | `2.017614e-05` | `8.062105e-04` |
| `1/8` (`0.1250000`) | `0.882496940024` | `0.882496902585` | `3.743988e-08` | `1.460712e-06` | `5.227996e-05` |
| `1/16` (`0.0625000`) | `0.939413064058` | `0.939413062813` | `1.244587e-09` | `9.833161e-08` | `3.328292e-06` |
| `1/32` (`0.0312500`) | `0.969233234516` | `0.969233234476` | `4.012091e-11` | `6.379380e-09` | `2.099575e-07` |
| `1/64` (`0.0156250`) | `0.984496437007` | `0.984496437005` | `1.273426e-12` | `4.062387e-10` | `1.318365e-08` |
| `1/128` (`0.0078125`) | `0.992217938260` | `0.992217938260` | `4.007905e-14` | `2.562872e-11` | `8.259059e-10` |
| **OLS Log-Log Slope** | — | — | **`4.93631`** ($O(x^5)$, $c_5 = \frac{1}{720}$) | **`3.92327`** ($O(x^4)$, $c_4 = \frac{1}{144}$) | **`3.98086`** ($O(x^4)$ observable) |

### 3.2 Finite-$N$ Critical Intervals, Rational Ceilings vs. Landau Pole, and Lee-Yang Positivity

1. **Certified Rational Critical Intervals (Converging Toward Onsager $u_c(\infty) = (\sqrt{2}-1)^2 \approx 0.1715729$):**
   - **Binder Cumulant Crossings $U_{4,L_1}(u^*) = U_{4,L_2}(u^*)$:**
     - $(4, 6)$: $u^* \in \left[\frac{28413177787}{171798691840},\, \frac{11365271115}{68719476736}\right] \approx [0.16538646,\, 0.16538646]$
     - $(4, 8)$: $u^* \in \left[\frac{57418078807}{343597383680},\, \frac{7177259851}{42949672960}\right] \approx [0.16710860,\, 0.16710860]$
     - $(6, 8)$: $u^* \in \left[\frac{29004656913}{171798691840},\, \frac{58009313827}{343597383680}\right] \approx [0.16882932,\, 0.16882932]$
   - **Specific-Heat Peaks $\frac{dc_{V,N}}{du}(u_{\text{peak}}) = 0$:**  
     $u_{\text{peak}}(16) \approx 0.1949230$, $u_{\text{peak}}(36) \approx 0.1880567$, $u_{\text{peak}}(64) \approx 0.1842213$ (monotonically approaching $u_c(\infty)$).
2. **Tight Finite-$N$ Ceilings vs. Continuum Landau Mean-Field Divergence:**
   - **Susceptibility Ceiling $B_\chi(N) = N$:** Attained with **$0/1$ exact slack** at $u = 0$ ($\chi_{\text{red},N}(0) = N$) across $N \in \{16, 36, 64\}$, whereas continuum Landau mean-field susceptibility $\chi_{\text{MF}}(\beta) = \frac{1}{zJ(\beta_c^{\text{MF}} - \beta)}$ (`landau_mean_field_susceptibility`) suffers an unphysical `ZeroDivisionError` pole at $\beta = \beta_c^{\text{MF}} = 0.25$.
   - **Specific-Heat Ceiling $B_c(L) = \frac{2}{3}H_{L/2} - \frac{3}{16}$:** Bounds peak $c_{V,N}(u)$ tightly across all sizes (`max_grid / B_c`: `0.807354 <= 13/16 = 0.812500` on $N=16$; `1.027814 <= 149/144 = 1.034722` on $N=36$; `1.175595 <= 173/144 = 1.201389` on $N=64$; $> 97.8\%$–$99.4\%$ empirical tightness!).
3. **Lee-Yang Unit-Circle Cayley Quadrance (`lee_yang_quadrance`):**  
   Across all 12 rational unit-circle probes $t \in [0, 5/4]$ at near-critical activity $u = 17/100$, $|Z_N(u, z(t))|^2 / Z_N(u, 1)^2 \in \mathbb{Q}_{>0}$ remains strictly positive (`min log10 ratio = -2.2450, -2.7585, -5.0835` for $N = 16, 36, 64$), confirming zero real phase-transition poles at finite $N$ while tracking the deepening Lee-Yang pinch toward $t \to 0^+$ as $N$ grows.

---

## 4. Battle 3: Coupled 2D $\phi^4$ Lattice Dynamics ($\lambda = 10, \kappa = 20$)

### 4.1 Smooth-Regime Trajectory Convergence vs. `RK4(dt = 1/512)` (`T_final = 1/4, L = 4`)
Evaluated on the smooth $C_4$-asymmetric subcritical lattice wavepacket (`default_smooth_initial_state_q(4)`) up to $T_{\text{final}} = 1/4$ against high-resolution `Classical RK4(dt = 1/512)`:

| Step Size $\Delta t$ | `Forward Euler` $L^2$ Error | `Velocity Verlet` $L^2$ Error | `Midpoint RK2` $L^2$ Error | `Classical RK4` $L^2$ Error | **Pure-$\mathbb{Q}$ Engine** $L^2$ Error |
| :---: | :---: | :---: | :---: | :---: | :---: |
| `1/16` (`0.0625000`) | `2.717722e-02` | `1.562664e-03` | `3.482395e-03` | `1.567881e-05` | **`2.542322e-02`** |
| `1/32` (`0.0312500`) | `1.468541e-02` | `3.879375e-04` | `8.117154e-04` | `9.350128e-07` | **`6.890377e-03`** (`3.69x`) |
| `1/64` (`0.0156250`) | `7.517506e-03` | `9.681671e-05` | `1.951847e-04` | `5.805008e-08` | **`1.772960e-03`** (`3.89x`) |
| `1/128` (`0.0078125`) | `3.790823e-03` | `2.419374e-05` | `4.780560e-05` | `3.618909e-09` | **`4.467196e-04`** (`3.97x`) |
| **OLS Log-Log Slope** | **`0.94915`** (`~1.0`) | **`2.00422`** (`~2.0`) | **`2.06164`** (`~2.0`) | **`4.02525`** (`~4.0`) | **`1.94503`** (`~2.0`, $O(\Delta t^2)$!) |

### 4.2 Stiff Strong-Quartic Coarse-$\Delta t$ Blowup vs. Exact $0/1$ Energy Conservation
On the stiff checkerboard initial state (`default_stiff_quartic_state_q(4)`, $\omega_{\max} = \sqrt{180} \approx 13.42$, $\omega_{\max}\Delta t = 3.35 > 2\sqrt{2}$ at $\Delta t = 1/4$) integrated for $n = 300$ coarse steps:

| Solver (`L = 4, lam = 10, kap = 20`) | $\Delta t = 1/4$ Status (`n = 300`) | $\Delta t = 1/8$ Status (`n = 300`) | Max Energy Drift $|H_n - E_0|$ ($\Delta t = 1/4$) | Quartic Bound $(\phi_i^2 - 1)^2 \le B_\phi$ |
| :--- | :---: | :---: | :---: | :---: |
| `float64` Forward Euler | **BLOWUP (`Inf`/`NaN`) @ step `12`** | **BLOWUP (`Inf`/`NaN`) @ step `16`** | $\infty$ | Violated |
| `float64` Velocity Verlet | **BLOWUP (`Inf`/`NaN`) @ step `6`** | Stable (oscillatory energy drift) | $\infty$ | Violated |
| `float64` Midpoint RK2 | **BLOWUP (`Inf`/`NaN`) @ step `7`** | **BLOWUP (`Inf`/`NaN`) @ step `12`** | $\infty$ | Violated |
| `float64` Classical RK4 | **BLOWUP (`Inf`/`NaN`) @ step `6`** | Stable (dissipative energy drift) | $\infty$ | Violated |
| **Pure-$\mathbb{Q}$ Hurwitz Engine** | **`300 / 300` Steps Completed** | **`300 / 300` Steps Completed** | **`0/1` (Exact Zero!)** | **`100%` Respected (`True`)** |

### 4.3 Topological $\mathbb{Z}_2$ Sector Locking, Palindrome Reversibility & 21-Prime Bit-Entropy Plateau
1. **Finite-$N$ Spontaneous Symmetry Breaking & Topological Sector Locking ($L = 4, N = 16, E_c(16) = 40$):**
   - **Subcritical Trajectory ($E_0 = 5.8959 < 40$):** Order parameter $M_n = \frac{1}{16}\sum_i \phi_{n,i}$ remains strictly locked inside the positive well $M_n \ge +M_{\min}(E_0) = 0.784904 > 0$ across all $n = 0 \dots 300$ steps (`subcritical_locked_all_300 == True`).
   - **Supercritical Trajectory ($E_0 = 114.3840 > 40$):** Surmounts the $M = 0$ saddle barrier and oscillates symmetrically between $+M$ and $-M$ wells (`supercritical_crossed_zero == True`).
   - **Exact Palindrome Reversibility ($\nu = 0$):** Advancing $10$ steps at $+\Delta t = +1/8$ followed by $10$ steps at $-\Delta t = -1/8$ recovers $(\phi_0, \pi_0)$ with **`0/1` max error**.
2. **Stationary $O(1)$ Bit-Length Entropy Across Clean Dyadic & All 21 Primes (`101..199`):**  
   Over $n = 300$ steps at $\Delta t = 1/8$ on $22$ distinct ensembles (clean dyadic + prime perturbations $\pm 1/p$ for every odd prime $p \in \{101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199\}$):
   - **Clean Dyadic Carrier:** Peak state bit-length = **`62` bits**; late-horizon OLS slope ($n \in [100, 300]$) = **`0.000000 bits/step`**.
   - **All 21 Odd Primes (`101..199`):** Peak state bit-length = **`62` bits** (`min = max = 62` across steps $1 \dots 300$); worst late-horizon OLS slope across all 21 primes = **`0.000000 bits/step`**.

### 4.4 $N = 64$ ($8 \times 8$) Deterministic SOS Cooling Quench & Unfrozen $2^{-30}$ Tangents
- **$N = 64$ Cooling Quench (`make_quench_config()`, $240$ steps, $\Delta t = 1/8, \nu = 1/4, E_c(64) = 160$):**
  - Completes all **`240 / 240` steps** without momentum zero-crossing (`quench_zero_div_step == null`).
  - Cools disordered supercritical initial state ($E_{\text{sys},0} = 228.2004 > 160$, $|M_0| = 0.02979$) smoothly across the critical barrier $E_c(64) = 160$ down to **$E_{\text{sys},240} = 57.0584 < 160$** ($Q_{\text{bath},240} = 171.1420$), spontaneously ordering into **$|M_{240}| = 1.09825 > M_{\min}(E_{\text{sys},240}) = 0.63469$**.
  - **Thermodynamic Laws:** Exact First Law $|(E_{\text{sys},n} + Q_{\text{bath},n}) - E_{\text{tot},0}| = \mathbf{0/1}$ across all checkpoints and strict Second Law monotonicity $\Delta Q_{\text{bath},n} \ge 0/1$ at every step (`second_law_monotonic == True`).
- **Unfrozen Discrete Tangents Down to $\Delta t = 2^{-30} \approx 9.31 \times 10^{-10}$:**  
  Because $\Delta t$ is factored strictly outside `qb(...)`, discrete finite-difference quotients $(R(F_{\Delta t}(X)) - R(X))/\Delta t$ converge cleanly to the Hamiltonian vector field $X_H(R(X))$ across $\Delta t \in \{2^{-4}, 2^{-8}, \dots, 2^{-30}\}$ with **`err(t=0, 2^-30) = 4.49e-9` (`slope = 1.0050`)** and **`err(t>0, 2^-30) = 4.37e-9` (`slope = 1.0064`)**.

---

## 5. Pure-Stdlib Scientific Plot Suite & `/tmp` `-I` Clean-Room Verification

### 5.1 Scientific Visualization Inventory (`plots/*.{svg,png}`)
Generated in `< 1.0s` by `simulator/plotter.py` (`DualCanvas` + `SubplotPanel`) using **zero** `matplotlib`/`PIL` dependencies (pure XML `.svg` + `struct`/`zlib` 24-bit RGB `.png` with embedded `7x10` bitmap typography, clipped viewports, dashed lines, shaded coexistence/envelope bands, and legends):

| Plot Stem (`plots/<stem>.{svg,png}`) | Resolution | `.svg` Size | `.png` Size | Ink Ratio | Scientific Panels (`2 x 2` Multi-Panel Layout) |
| :--- | :---: | :---: | :---: | :---: | :--- |
| `battle1_microcanonical_backbending_and_specific_heat` | `1440x980` | `46,940 B` | `52,563 B` | `0.55686` | (a) Case 1A monotonic $\beta_{\mathbb{Q}}(E/N)$; (b) Case 1B S-curve backbending vs canonical smearing; (c) $C_{\text{micro}} < 0$ vs $C_{\text{can}} \ge 0$; (d) IEEE-754 `OverflowError` cliff vs pure-$\mathbb{Q}$ $Z_{64}$. |
| `battle2_rational_pade_ising_and_critical_ceilings` | `1440x980` | `61,916 B` | `58,986 B` | `0.55828` | (a) Binder $U_{4,N}(u)$ & Sturm crossings; (b) $B_\chi(N)=N$ ceilings vs Landau `div/0` pole; (c) $c_{V,N}(u)$ peaks & $B_c(L)$; (d) Lee-Yang unit-circle quadrance $|Z_N(z(t))|^2 > 0/1$. |
| `battle3_hamiltonian_symmetry_breaking_and_quench` | `1440x980` | `48,611 B` | `54,312 B` | `0.55757` | (a) Subcritical $\mathbb{Z}_2$ sector locking vs supercritical crossing; (b) Stiff $\Delta t = 1/4$ `float64` blowup vs $0/1$ $\mathbb{Q}$ conservation; (c) 62-bit plateau across dyadic + 21 primes; (d) $N=64$ cooling quench across $E_c(64)=160$. |
| `convergence_and_finite_N_phase_summary` | `1440x980` | `33,040 B` | `49,623 B` | `0.55778` | (a) Smooth ODE slopes (`0.95, 2.00, 2.06, 4.03` & $\mathbb{Q}$-engine `1.95`); (b) Padé $R_{2,2}$ slopes (`4.94, 3.92, 3.98`); (c) Finite-$N$ critical scaling vs $1/N$; (d) Unfrozen discrete tangents down to $\Delta t = 2^{-30}$. |

### 5.2 Clean-Room `/usr/bin/python3 -I` Audit Summary (`standalone_ultrafinitist_statmech.py`)
Copied in complete isolation into a fresh empty `/tmp/uf_statmech_clean_*` directory and executed via `/usr/bin/python3 -I standalone_ultrafinitist_statmech.py`:

- **Total Artifact Size:** `295,400 bytes` (`7,933` lines total; `5,164` lines in Part I with `part1_ast_violations == []`).
- **Return Code & Diagnostics:** `returncode == 0`, `stderr == ""` (`0` warnings, `0` errors).
- **Wall-Clock Execution Time:** **`38.595 seconds`** (`4.66x` faster than the `180.0s` clean-room budget).
- **Clean-Room Plot Validation:** **`8 / 8` files** (`4` XML-validated `.svg` + `4` CRC32/IHDR/IDAT-decompressed non-blank `.png` files) verified inside `/tmp/.../plots/`.
