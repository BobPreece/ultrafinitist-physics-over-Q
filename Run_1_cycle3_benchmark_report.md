# Quantitative Benchmark Report — Exhibit I (`ultrafinitist-qft` v2.0)

**Benchmark Engine:** `standalone_ultrafinitist_qft.py` & `simulator/v2_benchmarks_and_plotter.py`  
**Audit Output:** `results/cycle4_v2_audit_report.json` (`elapsed_sec = 51.043s`, `status = ALL_4_GATES_PASSED`)  
**Publication Plots:** `plots/battle1_gaussian_torus_and_no_doublers.{svg,png}`, `plots/battle2_varga_positivity_and_fair_uv_loops.{svg,png}`, `plots/battle3_cayley_wavepacket_and_block_spin_rg.{svg,png}`, `plots/convergence_and_ultrafinitist_qft_summary.{svg,png}` (+ `ultrafinitist_qft_summary.{svg,png}`)

---

## 1. Battle 1: Dispersion Isotropy, Doubler Elimination & Cubic IR Chiral Scaling

Across $N = 2^m \in \{64, 128, 256, 512, 1024\}$ ($m \in \{6, 7, 8, 9, 10\}$) with bare mass $m_0 = 1/2$, Cycle 3 couples 4th-order Symanzik difference stencils ($K_{\text{Sym}}, D_{\text{Sym}}$) to the Equal-Corner Sum-of-Squares Biharmonic Wilson operator $W_{\text{flat}} = (L_x - \frac{1}{2}L_y)^2 + \frac{3}{4}L_y^2$ on the $C_4$-automorphic Gaussian-integer quotient torus $V_m = \mathbb{Z}[i]/\langle(1+i)^m\rangle$.

### 1.1 $(1+1)\text{D}$ Dispersion, Doublers & Cubic $8\times$ IR Chiral Suppression
| $N$ | Naive Doublers | Cycle 3 Doublers | Naive $\Delta_{\text{Lorentz, Dirac}}$ | Cycle 3 $\Delta_{\text{Lorentz, scalar}}$ | Cycle 3 $\Delta_{\text{Lorentz, Dirac}}$ | Cycle 2 $c_{\text{Dirac}}^2$ | Cycle 3 $c_{\text{Dirac}}^2$ | Cycle 3 IR Chiral $|W(p_1)|/|p_1|$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | 1 | **0** | 0.068593 | **0.000754** | **0.007455** | 1.5000 | **1.0000** | `5.904e-05` |
| **128** | 1 | **0** | 0.068593 | **0.000754** | **0.007455** | 1.5000 | **1.0000** | `7.390e-06` ($8.0\times\downarrow$) |
| **256** | 1 | **0** | 0.077276 | **0.000916** | **0.007908** | 1.5000 | **1.0000** | `9.240e-07` ($8.0\times\downarrow$) |
| **512** | 1 | **0** | 0.081845 | **0.001006** | **0.008111** | 1.5000 | **1.0000** | `1.155e-07` ($8.0\times\downarrow$) |
| **1024** | 1 | **0** | 0.081845 | **0.001006** | **0.008111** | 1.5000 | **1.0000** | `1.444e-08` ($8.0\times\downarrow$) |

### 1.2 $(2+1)\text{D}$ Gaussian Quotient Torus $V_m$ vs Naive Cubic ($N = 1024$ Summary)
- **BZ Corner Fermion Doublers:** Naive Cubic leaves **3 exact BZ-corner doublers** ($(\pi,0), (0,\pi), (\pi,\pi)$ at $E^2 = 1/4$), whereas $W_{\text{flat}} = (0, 16, 16, 16)$ lifts all 3 corners equally to $E_{\text{Dirac}}^2 = 4/1$ (**`0` doublers**).
- **Boson–Fermion Lightcone Speed Parity:** Cycle 2's 1st-order Wilson term inflated $c_{\text{Dirac}}^2 = 1.4952 \ne c_{\text{scalar}}^2 = 0.9968$; Cycle 3's $O(|p|^4)$ biharmonic $W_{\text{flat}}$ restores **$c_{\text{scalar}}^2 = 1.0000$ and $c_{\text{Dirac}}^2 = 1.0035$** (`1/1` in `QQJet(order=4)`).
- **Isotropy Gain:** In the Cycle 3 intrinsic spectral lower-20% window, Dirac $\Delta_{\text{Lorentz}}$ drops from `0.3336` (Naive Cubic) to `0.0156` (**`21.4x` Dirac isotropy gain**), and on the full analytic BZ lower-20% dual ball it achieves **`16.2x` gain** while eliminating odd-$m$ rectangular spikes at $N = 128, 512$.

---

## 2. Battle 2: Varga Positivity, Pure-`Q` Quadrature Brackets & Fair-Fight $1/N$ Loop Tables

### 2.1 Propagator Entrywise Positivity ($N = 64$ over Pure `Q`)
- **Cycle 3 Symanzik ($P = \frac{1}{12}(8I - A_{\text{nn}})^2$ Varga Splitting):** $R^{(1+1)} = \frac{11}{4} I \ge 0$ and $R^{(2+1)} = \frac{5}{12} I + \frac{1}{6} A_x A_y \ge 0$ at $m_0^2 = \frac{1}{4}$, certifying **$G_{xy} > 0/1$ strictly entrywise** ($G_{\min}^{(2+1)} = +0.002976 > 0$).
- **Unfactored Biharmonic Higher-Derivative Regulator ($\frac{1}{4}\sum_\mu L_\mu^2$):** Violates Varga splitting, suffering oscillatory negative spatial ringing **$G_{\min}^{(2+1)} = -0.035569 < 0$** (`-0.0870 < 0` in $(1+1)\text{D}$).

### 2.2 Remediated Fair-Fight $1/N$ Comparison vs Continuum Box-Shrinking $a \to 0$
Both **Cycle 3** and **Naive Cubic** use the **exact same $1/N$ per-site normalization** ($\tau_1 = \frac{1}{N}\operatorname{Tr}(G)$, $s_2 = \frac{1}{N}\sum_{x,y} G_{xy}^3$, $\lambda_{\text{site}} = \lambda_0 s_1 / \kappa_W(4)$ — **zero `* N` multiplier on Naive Cubic!**). They are contrasted against the **Continuum Box-Shrinking Regime** ($L_{\text{box}} = 8.0$, $a(N) = L_{\text{box}}/N^{1/d} \to 0$, $\Lambda(N) = \pi/a(N) \to \infty$):

#### $(1+1)\text{D}$ Fair-Fight $1/N$ Loop Table ($m_0 = 1/2, \lambda_0 = 1/2, d = 4$)
| $N$ | Cycle 3 $\tau_1(N)$ | Fair Cubic $\tau_1(N)$ | Pure-`Q` Mid-54 Tightness | Cycle 3 $R_m^{(1)}(N)$ | Cycle 3 $R_m^{(2)}(N)$ | Fair Cubic $R_m^{(2)}(N)$ | Continuum $a = 8/N \to 0$ $R_m^{(1)}(\Lambda)$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | 0.924989 | 0.970143 | **99.83%** | 2.156236 | **1.803138** | 1.728815 | `8.3379` |
| **128** | 0.924989 | 0.970143 | **99.83%** | 2.156236 | **1.803138** | 1.728815 | `9.4410` ($+\ln 2$ div) |
| **256** | 0.924989 | 0.970143 | **99.83%** | 2.156236 | **1.803138** | 1.728815 | `10.5442` ($+\ln 2$ div) |
| **512** | 0.924989 | 0.970143 | **99.83%** | 2.156236 | **1.803138** | 1.728815 | `11.6473` ($+\ln 2$ div) |
| **1024** | 0.924989 | 0.970143 | **99.83%** | 2.156236 | **1.803138** | 1.728815 | `12.7505` (Log UV blowup!) |

#### $(2+1)\text{D}$ Fair-Fight $1/N$ Loop Table ($m_0 = 1/2, \lambda_0 = 1/2, d = 4$)
| $N$ | Cycle 3 $\tau_1(N)$ | Fair Cubic $\tau_1(N)$ | Pure-`Q` Mid-54 Tightness | Cycle 3 $R_m^{(1)}(N)$ | Cycle 3 $R_m^{(2)}(N)$ | Fair Cubic $R_m^{(2)}(N)$ | Continuum $a = 8/\sqrt{N} \to 0$ $R_m^{(1)}(\Lambda)$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **64** | 0.346501 | 0.385501 | **98.53%** | 1.866252 | **1.792999** | 1.713913 | `2.0668` |
| **128** | 0.339323 | 0.380662 | **99.16%** | 1.848308 | **1.784648** | 1.708746 | `2.5800` ($\propto \sqrt{N}$ div) |
| **256** | 0.338173 | 0.376958 | **99.35%** | 1.845432 | **1.783191** | 1.704386 | `3.3090` ($\propto \sqrt{N}$ div) |
| **512** | 0.338080 | 0.376905 | **99.37%** | 1.845199 | **1.783072** | 1.704323 | `4.3422` ($\propto \sqrt{N}$ div) |
| **1024** | 0.338077 | 0.376855 | **99.37%** | 1.845192 | **1.783068** | 1.704264 | `5.8050` (Linear UV blowup!) |

- **2-Loop Drift Summary ($N = 64 \to 1024$):** Cycle 3 achieves **`0.0000%` (`1.13e-12%`)** in $(1+1)\text{D}$ and **`0.5569%`** in $(2+1)\text{D}$ (strictly below the $3\%$ gate), with smooth monotone convergence on $V_m$ across both even and odd $m$ (whereas Fair Cubic exhibits rectangular aspect-ratio steps between $N = 128, 256, 512$).

---

## 3. Battle 3: Coarse-$\Delta t$ `0/1` Conservation, 21 Primes (`101..199`), Slopes & Block-Spin RG

### 3.1 Coarse-$\Delta t$ Norm & Energy Drift Over 300 Steps ($N = 16$)
| Step Size $\Delta t$ | Steps $n$ | Pure-`Q` Double-Cayley + Vieta $\|\Delta\text{Norm}\|, \|\Delta E\|$ | `float64` `Euler` Max Drift | `float64` `Midpoint RK2` Max Drift | `float64` `RK4` Max Drift |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **$\Delta t = 1/4$** | 300 | **`0/1`, `0/1` (`0.00e+00`)** | `1.11e+12` (Blowup!) | `2.89e-01` (Instability) | `8.42e-04` (Dissipative) |
| **$\Delta t = 1/8$** | 300 | **`0/1`, `0/1` (`0.00e+00`)** | `6.34e+03` (Blowup!) | `7.76e-03` (Drift) | `1.34e-05` (Drift) |

### 3.2 Smooth Convergence Slopes ($T = 1/2$ vs `RK4(dt = 1/512)`) & Unfrozen `2^-30` Tangents
- **Global $L_2$ Convergence Slopes ($\Delta t \in \{1/4, 1/8, 1/16, 1/32, 1/64\}$):**
  - `float64` `Euler`: **`1.0091`** (textbook 1st-order)
  - `float64` `Midpoint RK2`: **`1.9996`** (textbook 2nd-order)
  - `float64` `RK4`: **`3.9996`** (textbook 4th-order)
  - **Pure-`Q` `step_wavepacket` + `vieta_readout`:** **`2.0755`** (textbook 2nd-order Strang/Cayley)
- **Unfrozen Discrete Tangents Down to $\Delta t = 2^{-30}$:**
  - At $t = 0$: $\|(\Psi(\Delta t) - \Psi(0))/\Delta t + i H_D \Psi(0)\|_2 = \mathbf{2.9653 \times 10^{-10}}$ (exact $O(\Delta t)$ slope `-1.00`, zero freeze).
  - At $t > 0$ ($n = 5, t = 5/8$): $\|(\Psi(t + \Delta t) - \Psi(t))/\Delta t + i H_D \Psi(t)\|_2 = \mathbf{2.9651 \times 10^{-10}}$ (zero freeze).
- **10-Step Palindrome, Causality & `(1+i)`-Adic Block-Spin RG:**
  - $F_{-\Delta t}^{10}(F_{+\Delta t}^{10}(\Psi_0)) - \Psi_0 \equiv \mathbf{0/1}$; single-step impulse outside $|x| \le 12$ on $\mathbb{Z}_{64}$ $\equiv \mathbf{0/1}$; `dead_ornament_slots == []`.
  - $R_{\text{block}} : V_{m+2} \to V_m$ across **`1024 -> 256`, `512 -> 128`, `256 -> 64`**: `0/1` charge conservation, `0/1` unit quantum norm, `0/1` $C_4$ covariance, and `0/1` IR lightcone speed $c_{\text{RG}}^2 = 1/1$ in `QQJet(order=4)`.

### 3.3 Bounded Bit Complexity & Theorem 6 Across Clean Dyadic + All 21 Primes (`101..199`)
Over $n = 300$ coarse steps at $\Delta t = 1/8$ under `RuntimeDyadicQuantizerBundle(bits=48, max_exp=64)`:

| Initial Wavepacket Case | Prime $p$ | Peak Bit-Length (`< 150` Cap) | Final Bit-Length ($n = 300$) | Asymptotic Slope (`steps 150..300`) | Theorem 6 `odd_part(den) in {1, p}` |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`clean_dyadic`** | — | **70** | 70 | **`0.000000 bits/step`** | **PASS (`{1}`)** |
| **Primes `101..127` (6 primes)** | `101, 103, 107, 109, 113, 127` | **76–77** | 76–77 | **`0.000000 bits/step`** | **PASS (`{1, p}`)** |
| **Primes `131..199` (15 primes)** | `131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199` | **77** | 77 | **`0.000000 bits/step`** | **PASS (`{1, p}`)** |
| **`2+1D` Gaussian Torus Check** | `101` | **77** | 77 | **`0.000000 bits/step`** | **PASS (`{1, 101}`)** |
