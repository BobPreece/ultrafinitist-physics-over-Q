# Multi-Agent Swarm Prompts & Adversarial Protocol

Each Exhibit in this repository was produced using a 4-agent adversarial research swarm (`@architect`, `@simulator`, `@inquisitor`, `@formalist`, plus stall monitor `@coordinator`) running on **Gemini 4 Carbon (High)** inside Google Antigravity.

## Why a 4-Agent Adversarial Swarm?
Single-agent LLM prompts routinely hallucinate proofs or smuggle textbook shortcuts (floating-point rounding, `.limit_denominator()`, `isinstance` test-bypass branches, or artificial numerical damping) to pass benchmarks. By splitting the workflow into four conflicting incentives across a mandatory 3-cycle loop (`Propose -> Simulate & Prove -> Red-Team Audit -> Refine`), every shortcut attempted by `@architect` in Cycles 1 and 2 was caught and outlawed by `@inquisitor` (via AST static analysis and runtime probes) and `@formalist` (via custom symbolic algebra over $\mathbb{Q}$).

* **Agent 1 (`@architect`):** Designs the discrete spatial lattices, birational Vieta manifold projections, and bounded-height rational time-steppers.
* **Agent 2 (`@simulator`):** Builds zero-handicap `float64` continuum baselines (`Euler`, `Velocity Verlet` / `Midpoint RK2`, `Classical RK4`, and fine reference `RK4`) alongside `SIGALRM`-guarded exact `fractions.Fraction` benchmark harnesses and pure-stdlib `.svg`/`.png` plotters.
* **Agent 3 (`@inquisitor`):** Enforces 4 Red-Team Gates (Gate 1: AST contraband & symbolic anti-spoofing scan; Gate 2: $\le 10,000$-bit $\mathcal{O}(1)$ plateau & `< 5s`/step speed audit; Gate 3: `0/1` reversibility, `dead_ornament_slots == []`, small-$\Delta t$ tangent convergence down to $2^{-30}$, and physical fidelity; Gate 4: fair-fight baseline convergence slopes) and records all rejected attempts in `Run_*_failure_log.md`.
* **Agent 4 (`@formalist`):** Builds a zero-dependency Computer Algebra System (`QQPoly`, `QQRat`, `QQJet`) over $\mathbb{Q}$ to prove exact `0/1` invariant identities, Sum-of-Squares operator positivity, Varga $M$-matrix density positivity, and closed-form rational bounds.

---

## Prompt 1: Exhibit I — Renormalization-Free Lattice QFT (`Run_1_*`)
* **Target:** Eliminate ultraviolet loop divergences (`Sigma_1`, `Sigma_2`) and Nielsen–Ninomiya fermion doublers on strictly finite Hilbert spaces ($D = d^N < \infty$) over $(1+1)\text{D}$ rings and $(2+1)\text{D}$ Gaussian-integer tori $\mathbb{Z}[i]/\langle(1+i)^m\rangle$ across $N \in [64, 128, 256, 512, 1024]$ without continuous FFTs, counterterm subtractions (`- delta_m`), or logarithmic RG curve-fitting.
* **Key Audit Gates:** See [`Run_1_failure_log.md`](Run_1_failure_log.md) and [`Run_1_cycle3_spec.md`](Run_1_cycle3_spec.md) for the progression from Cycle 1 (rejected for non-local chord superluminal signaling and naive doublers) and Cycle 2 (rejected for $\mathcal{O}(p^2)$ first-order Wilson mass dispersion distortion) to Cycle 3's winning 4th-order Symanzik + equal-corner sum-of-squares biharmonic $W_{\text{flat}}$ architecture and Varga squared-$M$-matrix Positivstellensatz bounds ([`Run_1_cycle3_specialized_bounds.md`](Run_1_cycle3_specialized_bounds.md)).

---

## Prompt 2: Exhibit II — Ultrafinitist Classical Mechanics Over $\mathbb{Q}$ (`Run_2_*`)
* **Target:** Construct explicit, reversible, $\mathbb{Q}$-valued discrete time-steppers across $n = 1,000$ coarse steps ($\Delta t \in \{1/2, 1/4, 1/8\}$) with zero `sqrt`, zero trig, $\le 10,000$-bit plateaus, and `0/1` invariant conservation across:
  * **Battle 1:** 2D Chaotic Pullen–Edmonds Quartic Oscillator ($H = \frac{1}{2}(p_1^2 + p_2^2) + \frac{1}{2}(q_1^2 + q_2^2) + 2 q_1^2 q_2^2$).
  * **Battle 2:** 3D Interacting 3-Body $S_3 \times SO(3)$ Cluster conserving all 7 $E(3)$ invariants $(H, \mathbf{P}, \mathbf{L})$ to `0/1`.
  * **Battle 3:** 2D Eccentric Newtonian Kepler $1/r$ Problem conserving energy $H$, angular momentum $L_z$, and the Laplace–Runge–Lenz vector $\mathbf{A}$ to `0/1` with synchronized physical time $t_n = n\Delta t$ and zero apsidal precession.
* **Key Audit Gates:** See [`Run_2_failure_log.md`](Run_2_failure_log.md) for how `@inquisitor` rejected Cycle 1's degree-23 Bellon–Viallet bit explosion (`94,326` bits at Step 4) and Cycle 2's 20 AST violations (`isinstance(..., Fraction)` spoofing, unrolled Heron root loops, non-conservative shadow curl particle freezing, and planar 3-body `ZeroDivisionError`), leading to Cycle 3's winning architecture ([`Run_2_winning_architecture.md`](Run_2_winning_architecture.md) and [`Run_2_final_symbolic_proofs.md`](Run_2_final_symbolic_proofs.md)).

---

## Prompt 3: Exhibit III — Singularity-Free 2D/3D Fluid Turbulence Over $\mathbb{Q}$ (`Run_3_*`)
* **Target:** Cross-pollinate Exhibit I's Gaussian-torus / $\mathbb{Z}_L^3$ 4th-order Symanzik stencils and Varga $M$-matrices with Exhibit II's birational Vieta projections to simulate 2D and 3D fluid turbulence over $n = 250\text{–}300$ coarse steps ($\Delta t \in \{1/2, 1/4, 1/8\}$) while proving a closed-form symbolic No-Blowup Theorem over $\mathbb{Q}$:
  * **Battle 1 (2D Gaussian Torus $N=64, 128$):** Dual `0/1` Energy–Enstrophy conservation (`1A` inviscid) and `0/1` First/Second Law thermodynamic closure $E_{\text{kin}} + Q_{\text{heat}} = E_0, \Delta Q_{\text{heat}} \ge 0/1$ (`1B` viscous biharmonic $W_{\text{flat}}$).
  * **Battle 2 (3D Taylor–Green Vortex Stretching on $\mathbb{Z}_L^3$, $N=216, 512$):** Sitewise `0/1` incompressibility $\operatorname{div}_{\text{Sym}}(u) = 0$, `0/1` momentum/energy conservation, and genuine 3D vortex stretching $Z_{\max}/Z_0 \ge 1.35$.
  * **Battle 3 (3D Anti-Parallel Vortex Tube Collision & Positive Density Transport):** Closed-form rational ceilings $B_u(E_0, N)$ and $B_\omega(E_0, N)$ (`0/1` slack) + strictly positive mass-conserving rational density transport $\min_x \rho_n(x) > 0/1$ vs continuum centered-stencil negative Gibbs ringing.
* **Key Audit Gates:** See [`Run_3_failure_log.md`](Run_3_failure_log.md) for the complete 3-cycle remediation audit, [`Run_3_winning_architecture.md`](Run_3_winning_architecture.md) for the unified `QuantizerBundle` and Double-Cayley / `[2/2]` Padé SOS formulas, and [`Run_3_final_symbolic_proofs.md`](Run_3_final_symbolic_proofs.md) for the machine-checked proofs.
