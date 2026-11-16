# No-Go Inquisitor Red-Team Failure & Audit Log — Exhibit I (`ultrafinitist-qft` v1.0–v2.0)

**Audit Artifacts:** `inquisitor/cycle1_audit_report.json`, `inquisitor/cycle2_audit_report.json`, `inquisitor/cycle3_audit_report.json`, `results/cycle4_v2_audit_report.json`  
**Constitutional Verifier:** `standalone_ultrafinitist_qft.py` & `formalist/verify_proofs.py`

---

## Audit Gate Specifications (Cycles 1–3 Red-Team Charter)

1. **Gate 1 — Continuum Contraband Audit:** Zero continuous Poisson sprinkling or pseudo-random number generators (`random`, `np.random`), zero continuous numerical quadrature (`scipy.integrate`), zero continuous FFT assumptions (`fft`) replacing discrete graph eigenvectors, zero ad-hoc renormalization counterterms (`delta_m`, `delta_lambda`, `sigma_ref`), and strictly rational combinatorial coupling scalings over $\mathbb{Q}$.
2. **Gate 2 — Nielsen–Ninomiya Fermion Doubler & Chirality Audit:** Zero spurious UV zero-energy or low-energy doubler modes (`doubler_count == 0` at normalized Dirichlet roughness $\tilde{R}_n \ge 0.25$ and $|E| < 2 m_0 = 1.0$), exact discrete Dirac Hermiticity ($\|H_D - H_D^\dagger\| < 10^{-10}$), and explicit quantification of chiral/axial symmetry breaking across IR modes.
3. **Gate 3 — Lorentz Dispersion, Isotropy & Finite Lightcone Causality Audit:** Proper geometric graph diameter scaling $\operatorname{diam}(G) = \Theta(N^{1/d_s})$, lower-20% relativistic dispersion violation $\Delta_{\text{Lorentz}} = \max |E^2(p) - (p^2 + m_0^2)| / E^2(p) < 0.05$, $(2+1)\text{D}$ directional shell anisotropy $< 0.06$, Boson–Fermion lightcone mismatch $|c_{\text{Dirac}}^2 - c_{\text{scalar}}^2| / c_{\text{scalar}}^2 < 0.08$, and zero superluminal step leakage ($W_{\text{step}}(T=1)$ leak outside causal stencil $= 0.0$).
4. **Gate 4 — Finite Weyl–Schwinger Algebra & Zero-Subtraction Loop Stability Audit:** Exact finite-$d$ Weyl–Schwinger clock-and-shift relations ($U^d = V^d = I_d$, $U V = \omega_d V U$ for $d \in \{3, 4, 5\}$ on $D = d^N < \infty$) and $< 3\%$ relative drift in 1-loop and 2-loop effective mass ratios $R_m^{(1)}(N), R_m^{(2)}(N)$ across $N \in \{64, 128, 256, 512, 1024\}$ without counterterm subtraction or trivial coupling decoupling.

---

## Cycle 0 Baseline Audit: Naive Periodic Cubic Lattice (`architect/baseline_cubic.py`) — **REJECTED**

- **Candidates Audited:** `baseline_cubic_1p1D` ($\mathbb{Z}_N$) and `baseline_cubic_2p1D` ($\mathbb{Z}_{L_x} \times \mathbb{Z}_{L_y}$, rectangular $8\times 16$ at $N = 128$ and $16\times 32$ at $N = 512$) with 1-hop symmetric derivative $D_{1,\mu}$ ($\sin p_\mu$) and standard graph Laplacian $L_\mu$ ($2 - 2\cos p_\mu$).
- **Red-Team Rejections:**
  - **Gate 2 FAIL (Nielsen–Ninomiya Fermion Doublers):** $\sin p_\mu = 0$ vanishes at every Brillouin-zone boundary $p_\mu = \pi$, producing **`1.0` exact doubler** in $(1+1)\text{D}$ ($p = \pi$, $|E| = m_0 = 0.500$) and **`3.0` exact doublers** in $(2+1)\text{D}$ at $(\pi, 0), (0, \pi), (\pi, \pi)$ ($|E| = m_0 = 0.500$).
  - **Gate 3 FAIL (Grid Anisotropy & Odd-$m$ Rectangular Bias):** At $N = 1024$, $(2+1)\text{D}$ scalar $\Delta_{\text{Lorentz}} = 0.0705 > 0.05$, Dirac $\Delta_{\text{Lorentz}} = 0.3336 \gg 0.05$, and Dirac shell anisotropy $= 0.2884 \gg 0.06$, exacerbated at odd $m$ ($N = 128, 512$) where $L_x \ne L_y$ breaks $90^\circ$ $C_4$ symmetry.
  - **Gate 4 Note (v1.0 Unscaled Mode-Sum Convention vs v2.0 Fair-Fight Remediation):** In v1.0 (`unscaled_mode_sum=True`), summing unnormalized modes without the $1/N$ volume factor caused naive $R_m^{(2)}$ to blow up $O(N)$ (`47.64 -> 747.31` in $(1+1)\text{D}$, `24.76 -> 374.24` in $(2+1)\text{D}$, `~241%` drift). *Constitutional Remediation in Cycle 4 (v2.0):* Removing the artificial `* N` multiplier and applying the **identical $1/N$ per-site normalization** to both Naive Cubic and Cycle 3 shows that per-site loop traces at fixed $a = 1$ are bounded on *both* finite graphs (while Naive Cubic still fails Gate 2 doublers, Gate 3 $O(p^2)$ anisotropy, and odd-$m$ $C_4$ symmetry, and continuum box-shrinking $a = L/N^{1/d} \to 0$ drives genuine logarithmic/linear UV divergence).

---

## Cycle 1 Audit: Deterministic Number-Theoretic Chord Graphs (`architect/cycle1_candidates.py`) — **HARD REJECTED**

**Report:** `inquisitor/cycle1_audit_report.json` (`total_elapsed_s = 210.722`)  
**Candidates Audited:**
- `cycle1_modular_qr_chord_1p1D`: Deterministic odd-square modular chords $S = \{\pm 1, \pm 9, \pm 25, \pm 49\} \bmod N$ on $\mathbb{Z}_N$ ($k = 8$).
- `cycle1_gaussian_prime_chord_2p1D`: Norm-1, norm-9, and norm-13 Gaussian-integer chords ($\{\pm 1, \pm i, \pm 3, \pm 3i, \pm 3 \pm 2i, \pm 2 \pm 3i\}$, $k = 16$) on $V_m = \mathbb{Z}[i]/\langle(1+i)^m\rangle$.

### Cycle 1 Verdict Summary
- **Gate 1 (Contraband):** **PASS** (`10` files scanned, `0` violations).
- **Gate 2 (Doublers):** **HARD FAIL**
  - Any odd-parity antisymmetric shift sum $\sum_{s \in S} w_s \sin(s \cdot p)$ vanishes identically at every half-period $p \in \{0, \pi\}^d$, leaving **`1.0` exact BZ-boundary doubler** (`1+1D` local) / **`19.0–35.0` chord-oscillation doublers** (`1+1D` full chord spectrum) and **`3.0` exact BZ-corner doublers** in $(2+1)\text{D}$ ($|E| = m_0 = 0.500$).
- **Gate 3 (Lorentz Dispersion, Boson–Fermion Speed & Lightcone Causality):** **CATASTROPHIC FAIL**
  - Unweighted long-range scalar chord Laplacians inflated the scalar lightcone speed squared to $c_{\text{scalar}}^2 = \sum_{s > 0} s^2 = 1^2 + 9^2 + 25^2 + 49^2 = 3108$ (`3088.7` measured in `1+1D`) and $c_{\text{scalar}}^2 = 35.12$ in `2+1D` while $c_{\text{Dirac}}^2 \approx 1.00$, yielding $\Delta_{\text{Lorentz, scalar}} = 0.9255$ (`1+1D`) / $0.8981$ (`2+1D`), $\Delta_{\text{Lorentz, Dirac}} = 0.1761$ (`1+1D`) / $2.1322$ (`2+1D`), Boson–Fermion speed mismatch $= 0.9361$ / $0.9616$, and `2+1D` shell anisotropy $= 0.2037$ / $0.3166$.
  - Long-range number-theoretic chords collapsed the graph hop diameter (`diam = 4..14` in `1+1D`, `4..8` in `2+1D`) and caused superluminal single-step jumps across `84` sites (`1+1D`) and `16` sites (`2+1D`), leaking **`100.0%` (`1+1D`) and `97.16%` (`2+1D`)** of $W_{\text{step}}(T=1)$ probability outside the local causal ball $r \le 2$.
- **Gate 4 (Loop Stability):** **PARTIAL FAIL** — Modular wrapping of length-49 and norm-13 chords on small tori ($N = 64, 128$) produced **`10.44% / 7.65%`** (`1+1D` 1L/2L) and **`24.58% / 14.97%`** (`2+1D` 1L/2L) mass-ratio drift ($> 3\%$ cap).

---

## Cycle 2 Audit: Strictly Local Isotropic Moore–Wilson Operators (`architect/cycle2_candidates.py`) — **REJECTED**

**Report:** `inquisitor/cycle2_audit_report.json` (`total_elapsed_s = 67.493`)  
**Candidates Audited:**
- `cycle2_local_wilson_ring_1p1D`: Nearest-neighbor ring $\mathbb{Z}_N$ with 1st-order Wilson mass $M_W = m_0 I + \frac{r_W}{2} L_1$ ($r_W = 1$).
- `cycle2_isotropic_moore_wilson_2p1D`: 9-point isotropic scalar stencil $K_{\text{iso}} = \frac{2}{3}(L_x + L_y) + \frac{1}{6}L_{\text{diag}}$ on $V_m$ with 1st-order Wilson mass $M_W = m_0 I + \frac{r_W}{2} K_{\text{iso}}$ ($r_W = 1$) and 1-hop Dirac derivative $\sin p_\mu$.

### Cycle 2 Verdict Summary
- **Gate 1 (Contraband):** **PASS** (`12` files scanned, `0` violations).
- **Gate 2 (Doublers):** **PASS** — `doubler_count = 0.0` in both `1+1D` ($\min |E|_{\text{UV}} = 1.325$) and `2+1D` ($\min |E|_{\text{UV}} = 1.803$).
- **Gate 4 (Weyl & Loop Stability):** **PASS** — `0.0000%` 1L/2L drift in `1+1D` and `0.8435% / 0.5353%` 1L/2L drift in `2+1D` ($< 3\%$); single-step $W_{\text{step}}(T=1)$ causal leak $= 0.0$; `2+1D` scalar shell anisotropy reduced to `0.00562`.
- **Gate 3 (Boson–Fermion Lightcone Speed & Dirac Isotropy):** **REJECTED**
  1. **1st-Order Wilson Cross-Term Speed Inflation ($c_{\text{Dirac}}^2 = 1.50 \ne c_{\text{scalar}}^2 = 1.00$):** Because $M_W(p) = m_0 + \frac{r_W}{2}|p|^2 + O(|p|^4)$, squaring $M_W(p)^2 = m_0^2 + m_0 r_W |p|^2 + O(|p|^4)$ injected an unwanted $+m_0 r_W |p|^2 = +0.50 |p|^2$ cross-term into $E_{\text{Dirac}}^2(p) = m_0^2 + (1 + m_0 r_W)|p|^2 + O(|p|^4)$, inflating $c_{\text{Dirac}}^2 = 1.5000$ (`1+1D`) and $1.4952$ (`2+1D`) against $c_{\text{scalar}}^2 = 1.0000 / 0.9968$ (`DL_Dirac = 0.1132 / 0.2396`, `BF_mismatch = 0.1412 / 0.4462`, $O(|p|)$ IR chiral breaking ratio `0.3430`).
  2. **Unimproved 1-Hop Dirac Anisotropy & Scalar Radial Droop:** Unimproved $\sin^2 p_x + \sin^2 p_y = |p|^2 - \frac{1}{3}(p_x^4 + p_y^4)$ left `2+1D` Dirac shell anisotropy at `0.1370 > 0.06`, while the $-\frac{1}{12}|p|^4$ radial droop of $K_{\text{iso}}$ left `2+1D` scalar $\Delta_{\text{Lorentz}} = 0.0941 > 0.05$.

---

## Cycle 3 Audit: 4th-Order Symanzik + Equal-Corner SOS Biharmonic Wilson on $V_m$ — **CERTIFIED PASS**

**Report:** `inquisitor/cycle3_audit_report.json` (`all_gates_passed = true`, `total_elapsed_s = 119.295`)  
**Candidates Certified:** `cycle3_symanzik_biharmonic_1p1D` & `cycle3_gaussian_symanzik_biharmonic_2p1D`.

- **Gate 1:** **PASS** (`16` files scanned, `0` violations).
- **Gate 2:** **PASS** — `doubler_count = 0.0` across all $N \in \{64..1024\}$; $W_{\text{flat}} = (L_x - \frac{1}{2}L_y)^2 + \frac{3}{4}L_y^2$ lifts all 3 nonzero BZ corners equally (`16 : 16 : 16`, $|E| = 2.00 = 4 m_0$); IR chiral breaking $|W_{\text{op}}(p_1)|/|p_1|$ drops by an exact factor of $8 = 2^3$ per doubling (`5.90e-5 -> 7.39e-6 -> 9.24e-7 -> 1.16e-7 -> 1.44e-8` in `1+1D`).
- **Gate 3:** **PASS** — $c_{\text{scalar}}^2 = 1.0000$, $c_{\text{Dirac}}^2 = 1.0000$ (`1+1D`) / $1.0035$ (`2+1D`); $\Delta_{\text{Lorentz, scalar}} = 0.00017$ (`1+1D`) / $0.01094$ (`2+1D`); $\Delta_{\text{Lorentz, Dirac}} = 0.00257$ (`1+1D`) / $0.01559$ (`2+1D`, **`21.4x` isotropy gain** over Naive Cubic's `0.3336`); Boson–Fermion mismatch $= 0.00276$ / $0.01192 < 0.08$; `2+1D` shell anisotropy $= 0.02547$ / $0.01978 < 0.06$; $W_{\text{step}}(T=1)$ leak $= 0.0$.
- **Gate 4:** **PASS** — Varga squared-$M$-matrix regular splitting $P = \frac{1}{12}(8I - A_{\text{nn}})^2$ certifies $G_{xy} > 0$ entrywise; 1L/2L mass-ratio drift $= 0.0000\% / 0.0000\%$ (`1+1D`, $R_m^{(2)} = 1.803138$) and $1.1383\% / 0.5562\%$ (`2+1D`, $R_m^{(2)}(1024) = 1.783068$) inside certified rational quadrature brackets.

---

## Cycle 4 (v2.0) Constitutional 4-Gate & `/tmp` `-I` Clean-Room Audit Certificate — **CERTIFIED PASS**

**Scope:** Surgical upgrade of Exhibit I to **Strict Part I Pure-`Q` Constitutional Parity**, adding Battle 3 float-free Strang-split 3-shear Double-Cayley + birational Vieta wavepacket dynamics over $\mathbb{Q}^2$, exact `(1+i)`-adic block-spin decimation $R_{\text{block}} : V_{m+2} \to V_m$, pure-stdlib `cas.py` + `verify_proofs.py` (`EXIT 0`), `1440x980` `DualCanvas` `.svg`+`.png` publication plotting, and **Fair-Fight $1/N$ Remediation** (removing v1.0's `* N` cubic baseline handicap).

| Constitutional Gate | Audited Criteria | v2.0 Result | Status |
| :--- | :--- | :--- | :--- |
| **Gate 1: Strict Part I & `formalist/` AST Purity** | Only `__future__`, `fractions`, `typing` in Part I; `0` `float`/`complex` literals or calls; `0` `math`/`cmath`/`ctypes`/`numpy`/`scipy`/`sympy`; `0` `sqrt`/`isqrt`/variable `**`/`round`/`limit_denominator`; `0` `while`, `0` `try/except`, `0` Newton–Heron, `0` `isinstance`/`type`/`hasattr`, `0` `max`/`min`/`abs`; `quantize_branch_sites == []` | `total_violations == 0` across Part I of `standalone_ultrafinitist_qft.py`, `architect/v2_pure_q_engine.py`, `formalist/cas.py`, `formalist/symbolic_cas_and_bounds.py`, and `formalist/verify_proofs.py` | **PASS (`0` violations)** |
| **Gate 2: Bounded Bit-Growth & 21-Prime Invariance** | 300 steps ($\Delta t = 1/8$) across clean dyadic + all **21 primes (`101..199`)** under `RuntimeDyadicQuantizerBundle(bits=48, max_exp=64)` (`5.0s/step` `SIGALRM`); peak `< 150` bits; asymptotic slope `<= 0.005 b/st`; Theorem 6 `odd_part(den) in {1, p}` | Peak `70` bits (clean) / `76–77` bits (`101..199`) `< 150`; asymptotic slope **`0.000000 bits/step`** across all 22 runs; `21/21` primes satisfy `odd_part(den) in {1, p}` (`v_p >= -1`) | **PASS (`0.0000 b/st`)** |
| **Gate 3: Symbolic `0/1` Proofs, Ornaments, `2^-30` Tangents & Block-Spin RG** | `/usr/bin/python3 -I formalist/verify_proofs.py` (`EXIT 0`); `0/1` Vieta norm & energy conservation; `QQJet(order=3)` Taylor PDE match; `0/1` 10-step palindrome & lightcone causality; `dead_ornament_slots == []`; unfrozen tangents down to $\Delta t = 2^{-30}$ at $t=0$ & $n=5$; `(1+i)`-adic $R_{\text{block}} : V_{m+2} \to V_m$ (`1024->256`, `512->128`, `256->64`) | All symbolic residuals `== 0/1`; `dead_ornament_slots == []`; `dt = 2^-30` tangent errors `~10^-10` at $t=0$ and $n=5$ (slope `-1.00`); `R_block` preserves `0/1` charge, `0/1` unit norm, `0/1` $C_4$ covariance, and `0/1` $c_{\text{RG}}^2 = 1/1$ in `QQJet(4)` | **PASS (`0/1` exact)** |
| **Gate 4: Fair-Fight $1/N$ Baselines, Coarse-$\Delta t$ Drift & Smooth Slopes** | Identical $1/N$ normalization on Naive Cubic & Cycle 3 (zero `* N` multiplier!) + continuum $a = L/N^{1/d} \to 0$ UV blowup; coarse $\Delta t \in \{1/4, 1/8\}$ 300-step `0/1` conservation vs `float64` `Euler`/`RK2`/`RK4` drift/blowup; smooth $T=1/2$ slopes vs `RK4(1/512)` | Fair $1/N$ parity verified across $N \in \{64..1024\}$ (`99.82%` / `98.53%` pure-`Q` brackets); `Euler`/`RK2` blow up $> 10^{12}$ and `RK4` drifts at $\Delta t = 1/4$ vs `0/1` Pure-`Q`; smooth slopes `Euler = 1.01`, `RK2 = 2.00`, `RK4 = 4.00`, `Pure-Q = 2.07` | **PASS** |
| **Clean-Room `/tmp` `-I` Isolation Audit** | Copy `standalone_ultrafinitist_qft.py` alone into fresh `/tmp/cleanroom_exhibit1_v2` and run `/usr/bin/python3 -I standalone_ultrafinitist_qft.py` | `returncode == 0`, `stderr == ""`, `53.76s` (`< 90s` cap), all 4 `1440x980` `.svg`+`.png` pairs + `ultrafinitist_qft_summary.{svg,png}` valid | **PASS** |
