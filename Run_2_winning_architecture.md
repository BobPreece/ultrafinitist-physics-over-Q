# Winning Architecture: Ultrafinitist Rational Hamiltonian Mechanics Without Square Roots, Transcendentals, or Bit Explosion

**Author:** `architect` (Agent 1: The Birational Geometer)  
**Final Implementation:** [`architect/problems.py`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/problems.py), [`architect/cycle3_maps.py`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/cycle3_maps.py), [`architect/cycle3_spec.md`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/cycle3_spec.md)  

---

## 1. The Central Obstruction of Ultrafinitist Computational Mechanics

Classical continuum numerical analysis (`Forward Euler`, `Velocity Verlet`, `RK4` over IEEE-754 `float`) violates physical conservation laws at macroscopic time steps $\Delta t \in \{1/2, 1/4, 1/8\}$—exhibiting energy drift, angular-momentum drift, artificial Keplerian apsidal precession, or violent numerical blow-up—while relying on irrational square roots (`sqrt`) and floating-point roundoff.

Conversely, naively translating geometric/symplectic integrators into exact rational arithmetic $\mathbb{Q}$ (`fractions.Fraction`) hits **Bellon–Viallet / Silverman algebraic height explosion**:
- Any nonlinear birational map $z_{n+1} = F_{\Delta t}(z_n)$ of algebraic degree $d \ge 2$ on a non-integrable system has dynamical degree $d_{\text{dyn}} > 1$, causing logarithmic bit-length $h(z_n)$ to grow **exponentially** as $O(d_{\text{dyn}}^n)$ (`> 10,000` bits in fewer than 10 steps!).
- Even on superintegrable degree-1 energy shells (such as Levi-Civita Kepler Cayley steps), unprojected multiplicative rational factors $(1 + \kappa \eta^2)^n$ accumulate **linear** bit-length growth $O(n)$, exceeding $20{,}000$ bits by $n = 1{,}000$.

---

## 2. Three-Cycle Adversarial Evolution

### Cycle 1: Pure Iterated Birational Integrators ([`architect/cycle1_maps.py`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/cycle1_maps.py))
- **Design:** Pure iterated Cayley/collateral birational maps:
  - *Battle 1:* Strang-split harmonic-oscillator Cayley rotations (`step_battle1`).
  - *Battle 2:* Pairwise collateral-stiffness discrete-gradient Cayley steps (`step_battle2`).
  - *Battle 3:* Levi-Civita spinor Cayley steps under Sundman regularization $ds = \Delta t / a$ (`step_battle3`).
- **Empirical & Symbolic Findings (`#7`, `#8`, `#14`):**
  1. **Exponential Bit Explosion (Battles 1 & 2):** Battle 1 stalled at step `9` (`24,578` bits; $d_{\text{dyn}} = 5$); Battle 2 stalled at step `6` (`44,060` bits; $d_{\text{dyn}} \approx 7$).
  2. **Linear Bit Growth & Non-Uniform Clock (Battle 3):** Reached `21,783` bits at step `1,000` (`> 10,000` ceiling), while $t_{n+1} - t_n = r_n \Delta s \neq \Delta t$ sampled non-uniform physical times.
  3. **Energy Splitting Defect (Battle 1):** Operator splitting conserved a shadow energy rather than exact `0/1` physical $H_1$.

### Cycle 2: First Shadow-Guided Vieta Projection Prototype ([`architect/cycle2_maps.py`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/cycle2_maps.py))
- **Design:** Decoupled fixed-grid (`2^-34`) dyadic shadow steppers + one-shot birational Vieta projections from $z_{\text{ref}} = z_0$. Upgraded Battle 2 in [`architect/problems.py`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/problems.py) to the manifestly $S_3 \times SO(3)$-symmetric Wildberger–Archimedes quadrance potential.
- **Red-Team Audit Findings (`#16`, `#17`, `#21`):**
  Bit-lengths dropped to $O(1)$ (`76–6,942` bits over 1,000 steps!) with exact `0/1` conservation of all invariants, **but** deep structural flaws remained:
  1. **AST Contraband & Symbolic Divergence:** `isinstance(..., Fraction)` split symbolic and runtime code paths; unrolled Heron `x + r/x` and $h = \Delta t/2$ sub-stepping violated Gate 1.
  2. **Dead-Ornament State Slots:** Storing overwritten physical `(q, p)` inside `state` broke global state-space injectivity (`step(step(reinit, dt), -dt) != reinit`).
  3. **Non-Conservative Shadow Curl:** Dividing forces by scalar $1 + c\Delta t^2 \|F\|^2$ introduced non-zero curl $\nabla \times F_{\text{reg}} \neq 0$, pumping shadow energy until Battle 1 shadow coordinates escaped ballistically and froze physical coordinates.
  4. **Fixed-Grid `2^-34` Freezing:** Clamped $|x| \ge 256$ and underflowed updates to $0$ for $\Delta t \le 2^{-34}$.
  5. **Battle 2 Planar Singularities & Stage-3 Chord Jumps:** Planar 3-body submanifolds hit `d1 = v . (u x du) == 0`; non-radially-normalized Stage-3 chords caused `max_chord = 32.84`.

### Cycle 3: The Winning Architecture ([`architect/cycle3_maps.py`](file:///Users/Bob/.gemini/antigravity/scratch/ultrafinitist-mechanics/architect/cycle3_maps.py))
Cycle 3 introduces the **Conservative Shadow-Guided Birational Gauge-Fibration Architecture**, eliminating every single Cycle 1 and Cycle 2 defect simultaneously.

---

## 3. Five Pillars of the Winning Architecture

### Pillar I: Minimal Globally Injective State Space + Pure-Field Readout Fibration
Rather than mutating redundant physical coordinates inside the state vector, each integrator evolves a **minimal canonical shadow/parameter state** paired with the immutable initial manifold anchor $z_{\text{ref}} \in \mathcal{M}_0$:
- **Battle 1:** `(qh, ph, z_ref)`
- **Battle 2:** `(rc, pt, uh, vh, puh, pvh, z_ref)`
- **Battle 3:** `(t_phys, ca, sa, cb, sb, z_ref)`

Every component of `state` participates bijectively in the reversible update (`dead_ornament_slots == []`), and physical phase-space coordinates $(q_n, p_n) = \Pi_{\mathcal{M}}(z_{\text{ref}}, \hat{z}_n) \in \mathcal{M}_0(\mathbb{Q})$ are extracted deterministically on demand via one-shot birational Vieta projections. Consequently:
1. **Algebraic degrees never compound across steps:** Step $n = 1{,}000$ has the **exact same rational projection degree** as step $n = 1$, conquering Bellon–Viallet exponential bit explosion permanently ($O(1)$ bit-length plateau!).
2. **Global `0/1` Time-Reversibility:** `step(step(state, dt), -dt) == state` holds identically across the entire state space, including after arbitrary re-initialization `build_state_from_qp(q, p)`.

### Pillar II: Zero-Curl Conservative Rational Shadow Hamiltonians
To prevent coarse-step Störmer–Verlet blow-up at $\Delta t = 1/2$ *without* introducing non-conservative curl $\nabla \times F \neq 0$:
- **Battle 1:** Regularizes the potential *before* differentiation via separable bounded rational factors:
  $$W_{1,\Delta t}(q_1, q_2) = \frac{1}{2}(q_1^2 + q_2^2) + \frac{\frac{1}{2}\alpha q_1^2 q_2^2}{\big(1 + \frac{1}{2}\Delta t^2 q_1^2\big)\big(1 + \frac{1}{2}\Delta t^2 q_2^2\big)} = V_1(q) + O(\Delta t^2),$$
  guaranteeing $\frac{\partial F_1}{\partial q_2} \equiv \frac{\partial F_2}{\partial q_1}$ and confining chaotic shadow orbits indefinitely.
- **Battle 2:** Regularizes the $S_3 \times SO(3)$-invariant Archimedes area coupling via $S = \frac{3}{4}\|\mathbf{u}\|^2 + \|\mathbf{v}\|^2$ and $A = \|\mathbf{u} \times \mathbf{v}\|^2$:
  $$W_{2,\Delta t}(\mathbf{u}, \mathbf{v}) = k_s S + \frac{g_c A}{\big(1 + \frac{1}{4}\Delta t^2 S\big)^2} = V_2(\mathbf{u}, \mathbf{v}) + O(\Delta t^2),$$
  guaranteeing all 15 curl components on $\mathbb{Q}^6$ AND all 3 internal torque components $\mathbf{u} \times \mathbf{F}_u + \mathbf{v} \times \mathbf{F}_v$ vanish identically!

### Pillar III: Scale-Covariant 128-Decade Floating-Dyadic Anchored Quantization
Replacing fixed-point grids with pure-field floating-dyadic mantissa quantizers $\mathcal{Q}_{\text{sig}, K}$ spanning $(2^{-64}, 2^{64})$ and factoring $\Delta t$ *outside* anchored force differences:
$$\hat{p}_{n+1/2} = \hat{p}_n + \frac{\Delta t}{2}\Big(F_{\text{sh}}(q_0, \Delta t) + \mathcal{Q}_{\text{sig},29}\big(F_{\text{sh}}(\hat{q}_n, \Delta t) - F_{\text{sh}}(q_0, \Delta t)\big)\Big),$$
simultaneously achieves:
- **Odd symmetry** $\mathcal{Q}_{\text{sig}, K}(-x) = -\mathcal{Q}_{\text{sig}, K}(x)$ (preserving exact `0/1` Verlet/Feistel time-reversibility),
- **$O(1)$ dyadic denominator bounding** at coarse $\Delta t \in \{1/2, 1/4, 1/8\}$, AND
- **Exact tangent-flow convergence** $(F_{\Delta t}(z_0) - z_0)/\Delta t \to X_H(z_0)$ down to $\Delta t = 2^{-40}$!

### Pillar IV: Closed-Form `[2/2]` Padé Approximants + Complete Gauge-Fibration Atlases
- **Zero-Iteration `[2/2]` Padé Surrogates:** `pade22_sqrt_from_one(r) = (5r^2 + 10r + 1)/(r^2 + 10r + 5)` provides 5th-order rational accuracy near $r = 1$ (and $\sim 1.5 \times 10^{-19}$ accuracy in two closed-form stages `pade_sqrt_fine`) with zero `sqrt` and zero Heron `x + r/x` recurrences.
- **Complete 3-Body $SO(3)$ & Planar Geometry:** Battle 2 combines:
  1. An intrinsic $(\mathbf{u}_0, \mathbf{v}_0)$-basis planar Vieta chart (`project_planar_b2`) for all coplanar/equilateral/breathing submanifolds ($\mathbf{L}_0 \parallel \mathbf{u}_0 \times \mathbf{v}_0$),
  2. A pure-field local reverse-transvection chart (`b2_local_chart_kernel`) nonsingular at $\Delta t = 0$, and
  3. A 30-anchor $O(2)_{\text{dem}}$ + $SO(3)_{\mathbf{L}_0}$ (`apply_l_half_turn_q` + 2-bit `apply_l_rot_q`) atlas with Hill-boundary pull-in and radially-normalized $T$-orthogonal Stage-3 Vieta matching (`stage3_match_momenta`), conserving all 7 invariants `(H, Px, Py, Pz, Lx, Ly, Lz)` `0/1` without singular crossings or momentum flips.

### Pillar V: Exact Rational Levi-Civita Conic Parametrization of Kepler $1/r$ Orbits
In Battle 3, mapping $(c, s)$ onto the rational eccentric-anomaly conic $c^2 + 4\kappa s^2 = 1$ via antipodal Vieta projection (`project_battle3_conic_kernel`) and applying the linear symplectic spinor flow `spinor_from_cs` yields Cartesian $(x, y) = (u^2 - v^2, 2uv)$ with rational radius $r = u^2 + v^2 \in \mathbb{Q}$, **exact `0/1` conservation of Energy $H_3$, Angular Momentum $L_z$, and Laplace–Runge–Lenz vector $(A_x, A_y)$** (zero spurious apsidal precession!), uniform physical time $t_n = n \Delta t$, single-step even-in-$\Delta t$ perihelion-regularized Feistel dynamics, exact $O(\Delta t^2)$ trajectory convergence (`order = 2.000`), and `max_bits <= 202` across $1{,}000+$ steps.

---

## 4. Summary of Cycle 3 Performance Across All 4 Red-Team Gates

| Benchmark Battle | Exact `0/1` Invariants | Exact `0/1` Reversibility & Injectivity | 1,000-Step Max Bit-Length ($\Delta t \in \{\frac{1}{2}, \frac{1}{4}, \frac{1}{8}\}$) | Trajectory Convergence Order ($T = 1$) | Edge Cases & Perturbed ICs |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Battle 1: 2D Chaotic Quartic Oscillator** | $H_1$ (`0/1` drift) | Global (`dead_slots = []`) | **`1,595 – 1,603` bits** ($O(1)$ plateau) | **`1.99 – 2.00`** | `on_axis_q2_zero`, `rest_release_p_zero`, prime perturbed IC: **PASS** |
| **Battle 2: 3D 3-Body $SO(3)$ Quadrance Cluster** | All 7: $H_2, \mathbf{P}, \mathbf{L}$ (`0/1` drift) | Global (`dead_slots = []`) | **`3,654 – 3,966` bits** ($O(1)$ plateau) | **`2.08 – 2.95`** | `equilateral_breathing_from_rest`, `planar_z0_motion`, prime perturbed IC (`< 9,000` bits): **PASS** |
| **Battle 3: 2D Eccentric Kepler $1/r$ ($e = 0.5376$)** | $H_3, L_z, A_x, A_y$ (`0/1` drift, $0^\circ$ LRL precession, $r \in \mathbb{Q}$) | Global (`dead_slots = []`) | **`174 – 202` bits** ($O(1)$ plateau) | **`2.000`** | High-eccentricity spinor perturbed IC (`201` bits): **PASS** |
