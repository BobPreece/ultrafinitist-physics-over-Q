# Methodological Transparency: Adversarial Swarm Constitutions & Prompts (`Exhibits I–V`)

To allow mathematicians, physicists, philosophers of science, and AI researchers to audit **what was human-specified as constitutional law** versus **what the AI agents invented autonomously**, this document records the exact multi-agent operating protocol, the 4-Gate red-team enforcement rubric, and the Exhibit I–V mission specifications executed on **Gemini 4 Carbon (High)** in Google Antigravity.

---

## 1. Why a 4-Agent Adversarial Separation of Powers?

Asking a single conversational LLM to "formulate physics without real numbers or infinity" reliably produces hand-waving prose, hidden `float64` casts, `.limit_denominator()` truncation hacks, or sandbagged baselines.

Instead, every exhibit was executed across **three mandatory adversarial cycles** (`Cycle 1: Propose -> Simulate & Prove -> Red-Team Audit` $\to$ `Cycle 2: Remediate & Unify` $\to$ `Cycle 3: Hardening & Clean-Room Synthesis`) with strict directory walls between four competing agents plus an orchestrator/stall-monitor:

| Agent Role | Exclusive Directory | Mandate & Incentive Structure |
| :--- | :--- | :--- |
| **`@architect`** | `architect/` | **Constructive Rational Geometer.** Designs the discrete topologies, rational operators, dyadic shadow carriers, and birational Vieta/Hurwitz/Hopf manifold readouts strictly over `fractions.Fraction` ($\mathbb{Q}$). Forbidden from editing test harnesses or symbolic provers. |
| **`@formalist`** | `formalist/` | **Zero-Dependency Symbolic CAS & Proof Engine.** Builds pure-Python computer algebra (`QQPoly`, `QQRat`, `QQJet`, Subresultant PRS GCD, primitive `sturm_chain`, `LDL^T` Cholesky & Pólya–Lukacs Sum-of-Squares) to certify universal `0/1` identities over indeterminate variables—wiring directly into `@architect`'s live Python AST via `IdentityQuantizerBundle`. |
| **`@simulator`** | `simulator/`, `plots/` | **Zero-Handicap IEEE-754 Baselines & Clean-Room Packager.** Implements full-precision textbook `float64` continuum comparators (`Forward Euler`, `Velocity Verlet`, `Midpoint RK2`, `Classical RK4`, canonical `math.exp`/`math.log` sums; `sandbag_flags == []`), benchmarks all 21 prime-perturbed ensembles, renders pure-stdlib `1440x980` `.svg`/`.png` plots, and assembles the single-file `standalone_*.py` artifact. |
| **`@inquisitor`** | `inquisitor/`, `failure_log.md` | **Hostile Red-Team Prosecutor.** Armed with unconditional **VETO** power across Gates 1–4. Scans ASTs for contraband or dead ornaments, injects off-axis/planar/lopsided/prime counterexamples, enforces `5.0s/step` `signal.SIGALRM` timers, logs every defect unredacted in `failure_log.md`, and runs `/usr/bin/python3 -I` in an empty `/tmp` sandbox. |

---

## 2. The Universal 4-Gate Constitutional Rubric (Evolved Across Exhibits I–V)

Each time `@architect` invented a loophole in an early cycle, that loophole was permanently codified into `@inquisitor`'s automated AST and physics gates for all subsequent cycles and exhibits:

1. **Gate 1 — Scope-Aware AST Contraband & Zero-Branch Protocol:**
   - **Part I Pure-$\mathbb{Q}$ Isolation:** Above `# === PART II`, standalone verifiers and `architect/` + `formalist/` modules may import *only* `__future__.annotations`, `fractions.Fraction`, and `typing`.
   - **Banned Operations:** Zero `float` or `complex` literals/casts, zero `math`/`cmath`/`numpy`/`scipy`/`sympy`/`ctypes`, zero `sqrt`/`isqrt` or fractional/negative `**`, zero `exp`/`log`/`sin`/`cos`/`atan2`/`fft`, zero `while` loops or `try/except`, zero `.limit_denominator()`/`round()`/`//` coordinate surgery, zero unrolled Newton–Heron `x + S/x` iterations (`heron_sites == []`), zero hidden `dt/k` sub-stepping (`substep_sites == []`), zero piecewise clamps (`max`, `min`, `abs`; `clamp_sites == []`), and zero `isinstance`/`type`/`hasattr` runtime-vs-symbolic branching (`quantize_branch_sites == []`).
2. **Gate 2 — Feasible $O(1)$ Bit-Entropy, `5.0s` `SIGALRM` Speed & 21-Prime Torture Suite:**
   - Every macro-step must complete in `< 5.0s` under a kernel `signal.SIGALRM` timer (no `multiprocessing.Pool`).
   - Across `250–1,000` coarse steps ($\Delta t \in \{1/2, 1/4, 1/8\}$) on clean dyadic initial states AND **all 21 odd prime-perturbed ensembles** ($p \in \{101, 103, \dots, 199\}$), state bit-length must stay strictly below hard ceilings with late-time OLS drift slope **`<= 0.005 bits/step`** (`0.0000 bits/step` plateau via `Theorem 6` odd-denominator invariance).
3. **Gate 3 — Independent Physics, Ornaments, Palindromes, `2^-30` Tangents & Symbolic Proofs:**
   - **No Hood Ornaments (`dead_ornament_slots == []`):** Partial-derivative sensitivity probes across both individual slots and sum/difference eigenmodes verify every state component actively couples to physical observables.
   - **Exact Reversibility & Unfrozen Tangents:** Exact `0/1` 10-step palindrome $F_{-\Delta t}^{10}(F_{+\Delta t}^{10}(X_0)) == X_0$, plus discrete tangent convergence $(F_{\Delta t}(X) - X)/\Delta t \to \text{RHS}(X)$ down to $\Delta t = 2^{-30}$ at **both** $t = 0$ and $t > 0$ ($n = 5$)—requiring $\Delta t$ to sit strictly *outside* `QuantizerBundle` calls.
4. **Gate 4 — Zero-Handicap Baselines, Smooth $O(\Delta t^2)$ Convergence & `/tmp` `-I` Isolation:**
   - Continuum `float64` baselines must achieve textbook convergence slopes (`Euler ~ 1.0`, `Verlet/RK2 ~ 2.0`, `RK4 ~ 4.0`) against ultra-fine `RK4(dt=1/512..1/2048)` in smooth regimes, while the pure-$\mathbb{Q}$ engine achieves `slope in [1.85, 2.15]`.
   - Copying `Run_*_standalone_*.py` alone into a fresh empty `/tmp` folder and executing `/usr/bin/python3 -I` must exit `0` with `stderr == ""` and emit all 4 valid `.svg` + `.png` plot pairs.

---

## 3. What Was Human-Seeded vs. What the Swarms Discovered Autonomously

| Exhibit | Human Prompt Constraints & Seed Hints | Autonomous Mathematical Discoveries by the Carbon Swarms |
| :--- | :--- | :--- |
| **Exhibit I (QFT & v2.0 Parity)** ([`Run_1_winning_architecture.md`](Run_1_winning_architecture.md)) | Finite Weyl–Schwinger algebra ($U^d = V^d = I$), discrete rings/tori, eliminate Nielsen–Ninomiya doublers & UV loop blow-ups without counterterms; v2.0 mandated pure-$\mathbb{Q}^2$ spinors, fair $1/N$ cubic normalization, Double-Cayley wavepackets & block-spin RG. | **Gaussian-Integer Quotient Torus** $V_m = \mathbb{Z}[i]/\langle(1+i)^m\rangle$ (exact $C_4$ symmetry at odd $m$); **Equal-Corner SOS Biharmonic Wilson Mass** $W_{\text{flat}} = (L_x - \frac{1}{2}L_y)^2 + \frac{3}{4}L_y^2$ (`16:16:16` corner lift, $c_{\text{Dirac}}^2 = c_{\text{scalar}}^2 = 1/1$, `21.4x` isotropy); **Varga Regular Splitting** $P = \frac{1}{12}(8I - A_{\text{nn}})^2$ proving $G_{xy} > 0/1$ for $0 < m_0^2 \le 2/3$; 54th-order rational Gauss–Lobatto loop brackets (`99.83%`); Pythagorean $U(1)$ Vieta readout (`70–77` bits) & `(1+i)`-adic block-spin RG. |
| **Exhibit II (Classical Mechanics)** ([`Run_2_winning_architecture.md`](Run_2_winning_architecture.md)) | Chaotic oscillator, 3D $S_3$ 3-body problem, and $1/r$ Kepler orbits at coarse $\Delta t = 1/2$ over $\mathbb{Q}$ without `sqrt` or `.limit_denominator()`; seeded Wildberger quadrances, Levi-Civita spinors & Vieta jumping as candidate directions. | **Conservative Shadow-Guided Birational Gauge-Fibration** (`15` zero-curl & `3` zero-torque polynomial identities) locking all 7 $E(3)$ invariants $(H, \mathbf{P}, \mathbf{L})$ to `0/1` without planar equilateral `0/0` singularities; synchronized physical-time Feistel conic kernel + antipodal Vieta readout locking $(H_3, L_z, A_x, A_y)$ to `0/1` (`76` orbits in `174` bits, $0^\circ$ precession); symbolic Ge–Marsden $O(\Delta t^3)$ defect proof. |
| **Exhibit III (3D Turbulence)** ([`Run_3_winning_architecture.md`](Run_3_winning_architecture.md)) | 2D enstrophy cascade, 3D Taylor–Green vortex stretching ($N \in \{216, 512\}$), and `amp=10` compressible 3D vortex-tube collisions; prove a closed-form finite rational vorticity ceiling. | **Sharp `0/1`-Slack 3D No-Blowup Theorem** ($\max \lVert u \rVert^2 \le \frac{2(2N+5)}{3N}E_0$, $\max \lVert \omega \rVert^2 \le \frac{65}{18h^2}E_0$ via discrete Parseval shell orthogonality); pole-free **Double-Cayley rotors**; tenth-power `[2/2]` Padé SOS acoustic majorants ($\rho_{\min} > 0/1$); proof that $\mu_{\text{conf}} = 0$ eliminates the symplectic eigenvalue boomerang (`0.0000 bits/step`, `260x` common-denominator speedup). |
| **Exhibit IV (Black Hole GR)** ([`Run_4_winning_architecture.md`](Run_4_winning_architecture.md)) | Relativistic perihelion shift + photon sphere, horizon-penetrating $r \to 0$ plunge with closed-form Kretschmann ceiling, and `amp=10` Oppenheimer–Snyder shell collapse over $\mathbb{Q}$. | **4-Spinor Hopf $S^2$ Fibration** + cubic Weierstrass chord points ($B_{\text{crit}}^2 = 27$); **7-Slot Ray–Envelope** ($\rho > 0/1$) achieving **bitwise-exact Step-20 `(0,0,0)` black hole origin transit** with 25-step `0/1` palindrome; **Causal-SOS Readout** ($u^t > 1$, `0` acausal flips); `0/1`-slack Kretschmann SOS ceiling $24 r_s^2 / h^6$ (`14, 12, 44, 29, 12, 2`) & universal `2000:1` mass-hierarchy density bound $B_\rho$. |
| **Exhibit V (StatMech & SSB)** ([`Run_5_winning_architecture.md`](Run_5_winning_architecture.md)) | Exact finite-$N$ Ising & 1st-order coexistence combinatorics ($L \in \{4,6,8\}$), rational Boltzmann factor + Lee–Yang/susceptibility ceilings, and 2D $\phi^4$ Hamiltonian symmetry breaking & quench without $N \to \infty$. | **Dihedral Meet-in-the-Middle Karatsuba Transfer Engine** counting $2^{64}$ microstates in `0.996s` & proving microcanonical S-curve backbending ($C_{\text{micro}} < 0$); Joukowsky half-degree Sturm isolation of Binder crossings (`< 2^-40`) & `0/1`-slack $B_\chi(N) = N$; `LDL^T` Cholesky SOS Finite-$N$ SSB barrier $E_c(N) = N\lambda/4$; **Degree-$(2,4)$ Hurwitz Quaternion Vieta Identity** $1 + 8t^7(t-1)$ collapsing an 8th-degree $\phi^4$ root to $t = -c_7/c_8 \in \mathbb{Q}$ (`63` peak bits). |

---

## 4. Verbatim Swarm Prompt Template (Reproducible Across Multi-Agent Harnesses)

To reproduce or adapt any exhibit on a new physical system (e.g., Yang–Mills gauge lattices, Maxwell electrodynamics, or magnetohydrodynamics), initialize a 4-agent workspace with the following constitutional template:

```markdown
# 3-CYCLE ADVERSARIAL SWARM MISSION: [Physical Domain] Over Pure Field `Q`

Initialize a 4-agent adversarial research swarm (`@architect`, `@simulator`, `@formalist`, `@inquisitor`) plus stall monitor (`@monitor`) running strictly on **Gemini 4 Carbon (High)**.

## I. NON-NEGOTIABLE ULTRAFINITIST CONSTITUTION
1. **Strict Part I Pure-`Q` Isolation:** Above `# === PART II`, import ONLY `__future__.annotations`, `fractions.Fraction`, and `typing`. Zero `float`/`complex`, zero `math`/`numpy`/`scipy`/`sympy`/`ctypes`, zero `sqrt`/`isqrt`, zero `exp`/`log`/trig/FFT, zero `while` or `try/except`, zero Newton–Heron `x + S/x`, zero `.limit_denominator()`/`round()`, zero `isinstance` dispatch (`quantize_branch_sites == []`), and zero piecewise clamps (`max`, `min`, `abs`).
2. **Zero-Branch `QuantizerBundle` & 21-Prime Immunity:** Factor `dt` strictly OUTSIDE `qb(...)` ($Z_{n+1} = Z_n + \Delta t \cdot \mathrm{qb}(\text{rate})$) with bitwise `d & (-d)` denominator clearing. Across 300–1,000 coarse steps on clean dyadic + all 21 primes `101..199`, enforce `< 5.0s/step` `SIGALRM` latency and `<= 0.005 bits/step` late-time drift (`Theorem 6`).
3. **Three Battles & Closed-Form `0/1`-Slack Ceilings:**
   - **Battle 1:** [Core geometric/kinematic invariant test]
   - **Battle 2:** [Singular/critical regime + closed-form rational ceiling `B(h, N) in Q` proved with `0/1` slack in `formalist/cas.py`, diverging iff `h -> 0` or `N -> inf`]
   - **Battle 3:** [Strongly nonlinear coupled dynamics + thermodynamic/RG quench]
4. **Zero-Handicap Baselines & `/tmp` `-I` Clean-Room Deliverable:** Compare honestly against full-precision `float64` `Euler`, `Verlet`, `Midpoint RK2`, and `RK4` (`sandbag_flags == []`). Package `standalone_*.py` so copying it alone to `/tmp` and running `/usr/bin/python3 -I` exits `0` (`stderr == ""`) and renders 4 `.svg` + `.png` `1440x980` `DualCanvas` plot pairs.
```

---

## 5. Interactive Red-Teaming Quickstart

Every claim in this repository is backed by executable assertions rather than prose:

```sh
# 1. Run all 5 standalone verifiers in isolated /tmp clean-room sandboxes:
bash verify_all_cleanroom.sh

# 2. Stress-test Exhibit II over 5,000 coarse steps (dt = 1/2):
/usr/bin/python3 -I Run_2_standalone_rational_mechanics.py --steps 5000

# 3. Quick-check Exhibit I v2.0 QFT into a custom plot directory:
/usr/bin/python3 -I Run_1_standalone_ultrafinitist_qft.py --quick --out-dir /tmp/qft_test_plots
```

For complete line-by-line histories of every rejected shortcut across Cycles 1–3, read [`Run_1_failure_log.md`](Run_1_failure_log.md), [`Run_2_failure_log.md`](Run_2_failure_log.md), [`Run_3_failure_log.md`](Run_3_failure_log.md), [`Run_4_failure_log.md`](Run_4_failure_log.md), and [`Run_5_failure_log.md`](Run_5_failure_log.md).
