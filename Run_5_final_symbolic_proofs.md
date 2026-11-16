# Machine-Checked Pure-$\mathbb{Q}$ Symbolic Proofs for Ultrafinitist Statistical Mechanics

**Author:** `@formalist` (Zero-Dependency Symbolic CAS & Proof Verifier)  
**CAS & Proof Suite:** `formalist/cas.py`, `formalist/proof_engines.py`, `formalist/verify_proofs.py`  
**Constitutional Compliance (Gate 1 AST Verified):** Pure $\mathbb{Q}$ (`int` and `fractions.Fraction` only; zero `float`/`complex`, zero `math.*`/`cmath.*`/`random`/`numpy`/`scipy`/`sympy`, zero `while` loops, zero `try/except`, zero `abs`/`max`/`min`/`round`/`limit_denominator`/`as_integer_ratio`, zero `isinstance`/`type`/`hasattr`/`getattr`, zero branch predicates on quantizer or critical thresholds).

---

## 0. Symbolic CAS Architecture (`formalist/cas.py`)

All algebraic identities, root isolations, and positivity certificates are machine-checked from first principles over $\mathbb{Q}$ using zero external packages:

1. **Polymorphic Rational-Tag Tower (`int`/`Fraction` $\subset$ `QQPoly` $\subset$ `QQRat` $\subset$ `QQJet` and `QQMVPoly` $\subset$ `QQMVRat`):**
   - Exact univariate polynomials `QQPoly` and reduced rational functions `QQRat` over $\mathbb{Q}[x]$ and $\mathbb{Q}(x)$ using Collins signed subresultant polynomial remainder sequences (PRS) accelerated by bounded Hensel $2$-adic modular inversion (`_subresultant_step_2adic`) and repeated-factor derivative reduction (`gcd(D, D')`).
   - Sparse exact $n$-variate polynomials `QQMVPoly` and rational functions `QQMVRat` over $\mathbb{Q}[x_0, \dots, x_{n-1}]$ (certifying full $16$-variable quartic SOS identities on $\mathbb{Z}_4^2$ and $63$-variable hyperplane identities on $\mathbb{Z}_8^2$ to exact $0/1$ monomial-by-monomial remainder).
   - Arbitrary-order truncated Taylor jets `QQJet(coeffs, order)` over any base ring over $\mathbb{Q}$, supporting exact jet addition, Cauchy convolution multiplication, valuation-shifted formal power-series division, and bounded binary exponentiation.
   - Zero-introspection polymorphic dispatch via `.numerator` / `.denominator` rational-tag invariants (`denominator >= 1` for `int`/`Fraction`; `-1..-5` for CAS types), requiring zero `isinstance`, `type`, or `hasattr` calls.
2. **Exact Subresultant Sturm Chains & Bisection Root Isolation (`sturm_chain`, `isolate_roots_sturm`):**
   - Single-pass primitive subresultant Sturm sequence construction with bounded `for` loops, certifying exact root counts on open intervals $(a, b)$, the positive half-line $(0, +\infty)$, and $(-\infty, +\infty)$ (`sturm_count_interval`, `sturm_count_halfline_pos`, `sturm_count_real`).
3. **Exact Rational Sum-of-Squares (SOS) & Positivity Engines:**
   - `sos_quad_univariate`: Exact Lagrange completion of squares over $\mathbb{Q}[x]$.
   - `cholesky_ldlt_sos`: Exact square-root-free symmetric pivoting $P^\top Q P = L D L^\top$ over $\mathbb{Q}$ (`0/1` matrix reconstruction residual).
   - `lukacs_positivity_halfline` & `lukacs_positivity_interval`: Exact Pólya–Lukacs SOS certificates $(1 + u)^r P(u) = \sum_j w_j g_j(u)^2 + u \sum_k v_k h_k(u)^2$ on $[0, +\infty)$ and Möbius-mapped compact intervals $[a, b]$.
4. **Unbranched AST Symbolic Execution (`IdentityQuantizerBundle`):**
   - Passes `QQPoly`, `QQRat`, and `QQJet` objects straight through `@architect`'s unbranched Python functions (`R_22`, `canonical_observables`, `init_state`, `step_phi4`, `readout`, `hamiltonian`) with `quantize_branch_sites == []`.

---

## 1. Theorem 1 — SOS $[2/2]$ Padé Boltzmann Algebra

Executed symbolically on `architect.battle2.R_22` over `QQPoly("x")`, `QQRat("x")`, and `QQJet(order=6)`:

$$R_{2,2}(x) = \frac{x^2 - 6x + 12}{x^2 + 6x + 12} = \frac{(x - 3)^2 + 3}{(x + 3)^2 + 3} \in \mathbb{Q}(x).$$

1. **Global Strict Positivity via SOS ($\forall x \in \mathbb{Q}$):**
   - `sos_quad_univariate` certifies with `0/1` remainder:
     $$N(x) = 1 \cdot (x - 3)^2 + 3 \ge 3 > \frac{0}{1}, \qquad D(x) = 1 \cdot (x + 3)^2 + 3 \ge 3 > \frac{0}{1}.$$
   - Hence $D(x)$ never vanishes on $\mathbb{R}$ and $R_{2,2}(x) > 0/1$ strictly for every $x \in \mathbb{Q}$.
2. **Exact Inversion Reciprocity in $\mathbb{Q}(x)$:**
   - Symbolic rational simplification proves:
     $$R_{2,2}(x) \cdot R_{2,2}(-x) - 1 = \frac{0}{1} \quad \text{identically in } \mathbb{Q}(x).$$
3. **Order-$4$ Exponential Consistency & Exact $5$th-Order Peano Coefficient:**
   - Evaluating `architect.battle2.R_22(QQJet.var(0, order=6))` yields exact coefficients:
     $$R_{2,2}(x) = 1 - x + \frac{1}{2}x^2 - \frac{1}{6}x^3 + \frac{1}{24}x^4 - \frac{1}{144}x^5 + O(x^6).$$
   - Comparing against $\sum_{m=0}^5 \frac{(-x)^m}{m!} = 1 - x + \frac{1}{2}x^2 - \frac{1}{6}x^3 + \frac{1}{24}x^4 - \frac{1}{120}x^5$ certifies `0/1` match through $O(x^4)$ and exact $5$th-order coefficient difference:
     $$[x^5]\bigl(R_{2,2}(x) - e^{-x}\bigr) = -\frac{1}{144} - \left(-\frac{1}{120}\right) = +\frac{1}{720}, \qquad [x^5]\Bigl(R_{2,2}\bigl(\tfrac{x}{2}\bigr)^2 - e^{-x}\Bigr) = +\frac{1}{2880}.$$
4. **Strict Sub-Unity & Monotonicity Certificates:**
   - $D(x) - N(x) = 12x > 0$ on $(0, +\infty)$ (`lukacs_positivity_halfline` remainder `0/1`), proving $0 < R_{2,2}(x) < 1$ for all $x \in \mathbb{Q}_{>0}$.
   - $-R_{2,2}'(x) = \frac{12(12 - x^2)}{(x^2 + 6x + 12)^2} > 0$ on $[0, 3]$ (`lukacs_positivity_interval` remainder `0/1`).

---

## 2. Theorem 2 — Finite-$N$ Criticality, Sturm Root Isolation, Lee-Yang Circle Theorem & Exact Ceilings

Executed symbolically on `architect.battle1.solve_case_1a(L)` and `architect.battle2.canonical_observables(L, u_sym)` for $L \in \{4, 6, 8\}$ ($N \in \{16, 36, 64\}$, $2^{16}$, $2^{36}$, $2^{64}$ exact microstates, $E_k = -2N + 4k$, $u = R_{2,2}(4\beta J) \in (0, 1)$):

### 2.1 Sturm-Isolated Rational Binder Cumulant Crossings $U_{4,L_1}(u) = U_{4,L_2}(u)$
For each pair $(L_1, L_2) \in \{(4, 6), (6, 8), (4, 8)\}$, `isolate_roots_sturm` on the numerator of $U_{4,L_1}(u) - U_{4,L_2}(u)$ (up to degree $194$) isolates the unique ferromagnetic crossing root in $(1/200, 95/100)$ to width $< 2^{-40}$ with strict opposite signs at $[u_{\text{lo}}, u_{\text{hi}}]$:

| Lattice Pair $(L_1, L_2)$ | Certified Lower Bound $u_{\text{lo}} \in \mathbb{Q}$ | Certified Upper Bound $u_{\text{hi}} \in \mathbb{Q}$ | Interval Width $u_{\text{hi}} - u_{\text{lo}}$ | Sturm Root Count |
| :--- | :--- | :--- | :--- | :--- |
| $(4, 6)$ | `36368867567419/219902325555200` | `4546108445951/27487790694400` | `189/219902325555200` ($< 2^{-40}$) | `1` (unique) |
| $(6, 8)$ | `9281490212173/54975581388800` | `37125960848881/219902325555200` | `189/219902325555200` ($< 2^{-40}$) | `1` (unique) |
| $(4, 8)$ | `7349514087299/43980465111040` | `9186892609171/54975581388800` | `189/219902325555200` ($< 2^{-40}$) | `1` (unique) |

Notice the finite-size drift toward the Onsager fixed point ($\sqrt{2} - 1$ squared activity $( \sqrt{2}-1 )^2 \approx 0.171573$): $u_\times(4,6) \approx 0.165386 < u_\times(4,8) \approx 0.167109 < u_\times(6,8) \approx 0.168829$, captured purely in finite-$N$ rational intervals without $N \to \infty$.

### 2.2 Joukowsky Half-Degree Reduction & Specific-Heat Peak Roots $\frac{d c_{V,N}}{du} = 0$
By bipartite energy symmetry $\Omega_N(-E) = \Omega_N(E)$ (`u^N Z_N(1/u) == Z_N(u)`) and anti-reciprocity $L_{2,2}(1/u) = -L_{2,2}(u)$, the canonical specific heat $c_{V,N}(u)$ satisfies $c_{V,N}(1/u) = c_{V,N}(u)$ identically in $\mathbb{Q}(u)$, so $c_{V,N}(u) = \widetilde{c}_{V,N}(w)$ is an exact rational function of half degree in the Joukowsky variable:
$$w(u) = \frac{(u - 1)^2}{u} \in [0, +\infty) \quad \forall u \in (0, +\infty),$$
verified to `0/1` remainder in `QQRat("u")` for all $L \in \{4, 6, 8\}$. Isolating $\frac{d \widetilde{c}_{V,N}}{dw} = 0$ via `isolate_roots_sturm` and cross-checking `QQJet(order=1)` directly on `architect.battle2.canonical_observables(L, QQJet.var(u, order=1))` certifies strict maximum sign changes $\frac{d c_{V,N}}{du}(u_{\text{lo}}) > 0/1 > \frac{d c_{V,N}}{du}(u_{\text{hi}})$:

| Side $L$ ($N = L^2$) | Certified Peak Interval $[u_{\text{lo}}, u_{\text{hi}}]$ | Certified Lower Peak Witness $\underline{c}_V \in \mathbb{Q}$ |
| :--- | :--- | :--- |
| $L = 4$ ($N = 16$) | `[247094293497643787680373736053187/1267650600228229401496703205376000, 3953508695980332826001882687426107/20282409603651670423947251286016000]` | `6786277/8388608` ($\approx 0.808987$) |
| $L = 6$ ($N = 36$) | `[3814243667971012347914642663874381/20282409603651670423947251286016000, 7628487335975501599997111478870371/40564819207303340847894502572032000]` | `17244717/16777216` ($\approx 1.027865$) |
| $L = 8$ ($N = 64$) | `[186822579766479361689671151290409/1014120480182583521197362564300800, 14945806381382504544203846510291081/81129638414606681695789005144064000]` | `19778637/16777216` ($\approx 1.178899$) |

### 2.3 Machine-Checked Lee-Yang Circle Theorem ($L = 4$) & Unit-Circle Quadrance ($L \in \{4, 6, 8\}$)
Parameterizing the unit circle $|z| = 1$ rationally via $z(t) = \frac{1 - t^2}{1 + t^2} + i \frac{2t}{1 + t^2}$ ($t \in \mathbb{Q}$):
- On $L = 4$ ($N = 16$, $u = 1/2$), `cayley_circle_quadrance_poly` expands the exact $65{,}536$-state fugacity polynomial $P_{16}(z)$ into $|P_{16}(z(t))|^2 = \frac{\text{Re}(t)^2 + \text{Im}(t)^2}{(1 + t^2)^{16}}$ (`0/1` match with `architect.battle2.lee_yang_quadrance(4, 1/2, t)` at $t \in \{0, 1/32, 1/16\}$).
- Computing $G_{\text{LY}}(t) = \gcd(\text{Re}(t), \text{Im}(t))$ in $\mathbb{Q}[t]$ yields $\deg(G_{\text{LY}}) = 16$, and `sturm_count_real(sturm_chain(G_LY)) == 16` proves that **all $16$ complex zeros of $Z_{16}(1/2, z)$ lie strictly on $|z| = 1$** with finite-$N$ edge spectral gap:
  $$t_1(16) \in \left[\frac{20432765}{134217728},\, \frac{10216385}{67108864}\right] > \frac{1}{8} > \frac{0}{1}.$$
- Across all $L \in \{4, 6, 8\}$ and $t \in \{0, 1/32, 1/16\}$, $|Z_N(1/2, z(t))|^2 \in \mathbb{Q}_{>0}$ is strictly positive, forbidding zero-field phase singularities at finite $N$.

### 2.4 Global No-Divergence Ceilings $B_\chi(N)$ (`0/1` Slack) & $B_c(N, J)$ ($> 98\%$ Sturm-Certified)
1. **Susceptibility Ceiling $B_\chi(N) = N$ (`0/1` Zero Slack):**
   - $\chi_{\text{red},N}(0) = \frac{2 N^2}{N \cdot 2} = N$ (`0/1` attainment slack at $u = 0$).
   - For all $u \in (0, +\infty)$, $N - \chi_{\text{red},N}(u) = \frac{\sum_k (N^2 \Omega_N(E_k) - \sum_{\mathcal{S}(E_k)} M^2) u^k}{N \sum_k \Omega_N(E_k) u^k} \ge \frac{0}{1}$ is certified by `lukacs_positivity_halfline` at Pólya degree `0` (`0/1` remainder).
2. **Onsager-Harmonic Specific-Heat Ceiling $B_c(L) = \frac{2}{3} H_{L/2} - \frac{3}{16}$:**
   - Using the Joukowsky bijection $u \in (0, +\infty) \mapsto w = \frac{(u-1)^2}{u} \in [0, +\infty)$, `sturm_count_halfline_pos_poly` certifies **zero roots** of $B_c(L) - \widetilde{c}_{V,N}(w) = 0$ on $[0, +\infty)$ for every $L \in \{4, 6, 8\}$, proving $\sup_{u > 0} c_{V,N}(u) < B_c(L)$ globally over $\mathbb{R}_{>0}$ with $> 98\%$ certified tightness:
     - $L = 4$: $B_c(4) = \frac{13}{16} = 0.8125$, certified tightness $\ge \frac{6786277}{6815744} \approx 99.57\%$ (`B_c_tight = 1696739/2097152`, $> 99.98\%$).
     - $L = 6$: $B_c(6) = \frac{149}{144} \approx 1.03472$, certified tightness $\ge \frac{155202453}{156237824} \approx 99.34\%$ (`B_c_tight = 17246443/16777216`, $> 99.98\%$).
     - $L = 8$: $B_c(8) = \frac{173}{144} \approx 1.20139$, certified tightness $\ge \frac{178007733}{181403648} \approx 98.13\%$ (`B_c_tight = 2472577/2097152`, $> 99.98\%$).

---

## 3. Theorem 3 — Finite-$N$ Spontaneous Symmetry Breaking Theorem (Zero $N \to \infty$)

For the periodic 2D lattice $\phi^4$ potential on $\mathbb{Z}_L^2$ ($N = L^2 \in \{16, 64\}$, $\lambda = 10$, $\kappa = 20$, undirected nearest-neighbor bonds):
$$V(\phi) = \sum_{i=1}^N \frac{\lambda}{4}(\phi_i^2 - 1)^2 + \frac{\kappa}{2}\sum_{\langle i, j\rangle} (\phi_i - \phi_j)^2, \qquad M(\phi) = \frac{1}{N}\sum_{i=1}^N \phi_i, \qquad \psi_i = \phi_i - M(\phi).$$

1. **Exact Rational Laplacian $LDL^\top$ Cholesky SOS on $M(\phi) = 0$:**
   - **$L = 4$ ($N = 16$):** Restricting $\Delta_4 - c_4 I_{16}$ to the $15$-dimensional hyperplane $\sum_{i=1}^{16} \psi_i = 0$ at $c_4 = 2$, `cholesky_ldlt_sos` computes an exact rational $P^\top B_4 P = L_4 D_4 L_4^\top$ of rank `11` (reflecting the $4$-fold degeneracy of the minimal nonzero Fourier modes $(1, 0), (3, 0), (0, 1), (0, 3)$ with eigenvalue $4 - 2\cos(\pi/2) - 2 = 2$) with `0/1` matrix reconstruction residual, while $c_4 + \frac{1}{1000}$ turns strictly indefinite. Thus $c_4 = \frac{2}{1}$ with `0/1` slack, and $\kappa c_4 = 40 \ge \lambda = 10$.
   - **$L = 8$ ($N = 64$):** Restricting $\Delta_8 - c_8 I_{64}$ to the $63$-dimensional hyperplane $\sum_{i=1}^{64} \psi_i = 0$ at rational Poincaré constant $c_8 = \frac{239}{408} < 2 - \sqrt{2} < \frac{577}{985}$ (width $\frac{1}{401880}$), `cholesky_ldlt_sos` computes a full-rank (`rank = 63`, strictly positive definite, max pivot bit-length $\le 264$) rational $LDL^\top$ decomposition with `0/1` reconstruction error, and expands all $2{,}016$ quadratic monomials in `QQMVPoly` to verify $\sum_{\langle i,j\rangle}(\psi_i - \psi_j)^2 - \frac{239}{408}\sum_i \psi_i^2 = \sum_{k=1}^{63} d_k \ell_k(\psi)^2$ with `0/1` polynomial remainder. Moreover $\kappa c_8 = 20 \cdot \frac{239}{408} = \frac{1195}{102} > 10 = \lambda$.
2. **Exact Master Quartic SOS Identity & Critical Barrier $E_c(N) = \frac{N\lambda}{4}$:**
   - Letting $S_2 = \frac{1}{N}\sum_{i=1}^N \phi_i^2 = M^2 + R$ where $R = \frac{1}{N}\sum_{i=1}^N \psi_i^2 \ge 0$:
     $$V(\phi) - \frac{N\lambda}{4}(1 - M^2)^2 = \frac{\lambda}{4}\sum_{i=1}^N (\phi_i^2 - S_2)^2 + \frac{\kappa}{2}\sum_{k=1}^r d_k \ell_k(\psi)^2 + \frac{N\lambda}{4} R^2 + \frac{N(\kappa c_L - \lambda)}{2} R + \frac{N\lambda}{2} M^2 R \ge \frac{0}{1}.$$
   - On $\mathbb{Z}_4^2$, `prove_finite_n_ssb_sos` expands all **$4{,}012$ degree-$4$ monomials in $16$ variables** in `QQMVPoly` and verifies `lhs - rhs == 0` (`0/1` remainder!).
   - Setting $M = 0$ gives $V(\phi)\big|_{M=0} - \frac{N\lambda}{4} \ge \frac{0}{1}$, attained with `0/1` slack at $\phi = 0$. Hence the exact finite-$N$ symmetry-breaking energy threshold is $E_c(N) = \frac{N\lambda}{4}$ ($E_c(16) = 40$, $E_c(64) = 160$).
3. **Strict Topological Disjointness of $\mathbb{Z}_2$ Sectors for Every $E_0 < E_c(N)$:**
   - Whenever $H(\phi, \pi) \le E_0 < E_c(N)$, $V(\phi) \le E_0$ implies $(1 - M(\phi)^2)^2 \le \frac{E_0}{E_c(N)} < 1$, yielding a strictly positive certified rational lower bound (`ssb_m_min_sq_bound`):
     $$M(\phi)^2 \ge M_{\min}^2(E_0) > \frac{0}{1} \quad \text{(e.g. at } E_0 = \tfrac{3}{4}E_c(N),\ M_{\min}^2 > \tfrac{287250861}{2147483648} > \tfrac{1}{8}\text{)},$$
     splitting the rational energy shell $\mathcal{M}_{E_0}(\mathbb{Q})$ into two disjoint $\mathbb{Z}_2$-conjugate sectors $\mathcal{M}_{E_0}^+(\mathbb{Q}) = \{M > 0\}$ and $\mathcal{M}_{E_0}^-(\mathbb{Q}) = \{M < 0\}$ separated by the forbidden hyperplane $M = 0$ at strictly finite $N \in \{16, 64\}$.
4. **Exact $0/1$-Slack Single-Site Quartic Ceiling $B_\phi(E_0, \lambda, N) = \frac{4 E_0}{\lambda}$:**
   - $\frac{4}{\lambda} V(\phi) - (\phi_i^2 - 1)^2 = \sum_{j \ne i}(\phi_j^2 - 1)^2 + \frac{2\kappa}{\lambda}\sum_{\langle u, v\rangle}(\phi_u - \phi_v)^2 \ge \frac{0}{1}$ (`0/1` remainder in `QQMVPoly`), with `0/1` slack attained at the ordered ground states $\phi = \pm \mathbf{1}$ ($E_0 = 0, B_\phi = 0$) and at $E_0 = E_c(16)$ for homogeneous zero-mode excitations.

---

## 4. Theorems 4, 5 & 6 — Symbolic Battle 3 Stepper Verification & $21$-Prime $p$-adic Invariance

Executed directly on `@architect`'s Cycle 3 unbranched `architect/battle3.py` AST (`init_state`, `step_phi4`, `readout`, `hamiltonian`) — incorporating the $C_4$-covariant orbit-paired Hurwitz quaternion frame $(I, J, K)$, global SOS cubic-tamed Verlet carrier, single-stage degree-$(2, 4)$ Hurwitz birational Vieta energy-shell readout $t = -c_7 / c_8$ with 2-stage $[4/3]$ binomial Padé gauge scalars $\frac{1 + 28w + 70w^2 + 28w^3 + w^4}{8 + 56w + 56w^2 + 8w^3}$ (`_b3_pade_sqrt_pos`), and `FractionQuantizerBundle(shift_bits=56)` (`DEFAULT_QB`) — via `IdentityQuantizerBundle`, `JetTangentBundle(DEFAULT_QB)`, and `DEFAULT_QB`:

1. **Theorem 4 (Order-$3$ Taylor Jet Continuum $\phi^4$ ODE Consistency at $t = 0$ and $t > 0$):**
   - **Exact $t = 0$ Continuum ODE Match (`0/1`):** Passing `dt = QQJet.var(0, order=3)` at $\nu = 0$ through `step_phi4` and `readout` via `IdentityQuantizerBundle` on an asymmetric rational state on $\mathcal{M}_{E_0}(\mathbb{Q})$ certifies for every lattice site $i \in \{0, \dots, 15\}$:
     $$\phi_i(\Delta t) = \phi_i(0) + \pi_i(0)\,\Delta t - \frac{1}{2}\nabla_{\phi_i} V(\phi(0))\,\Delta t^2 + O(\Delta t^3),$$
     $$\pi_i(\Delta t) = \pi_i(0) - \nabla_{\phi_i} V(\phi(0))\,\Delta t - \frac{1}{2}\bigl(D^2 V(\phi(0))[\pi(0)]\bigr)_i\,\Delta t^2 + O(\Delta t^3),$$
     with `0/1` jet-coefficient residuals (`c0`, `c1`, `c2`) and `0/1` energy conservation across all four jet orders (`c0..c3`), proving $O(\Delta t^3)$ local truncation error and $O(\Delta t^2)$ global consistency.
   - **Evolved $t > 0$ Unfrozen Jet & $O(\Delta t_{\text{warm}}^2)$ Tangent Convergence:** Lifting dyadic trajectory states at $t = 1/4$ (after coarse warmup $\Delta t_{\text{warm}} = 1/8$ vs fine warmup $\Delta t_{\text{warm,fine}} = 1/32$, plus $\Delta t_{\text{warm}} = 1/256$) into `QQJet(order=3)` via `JetTangentBundle(DEFAULT_QB)` certifies unbranched execution, `0/1` base-point (`c0`) agreement, strictly non-zero unfrozen velocity tangent (`||c1||^2 > 0`), exact `0/1` symplectic energy conservation $H(R(F_{\varepsilon}(\text{st}_t))) - E_0 = (0/1, 0/1, 0/1, 0/1)$ across all four jet orders `c0, c1, c2, c3` (for both `hamiltonian` and independent continuum $H_{\text{ref}}$), and quadratic $O(\Delta t_{\text{warm}}^2)$ convergence of the exact jet tangent `c1` to $X_H(R(\text{st}_t))$ ($\sim 160\times$ squared-error reduction from $\Delta t_{\text{warm}} = 1/8 \to 1/32$, and $< 10^{-4}$ squared L2 error at $\Delta t_{\text{warm}} = 1/256$).
2. **Theorem 5 (Non-Perturbative Rational Identities in $\mathbb{Q}(\Delta t)$, $\mathbb{Q}(\nu)$, and Formal Power Series $\mathbb{Q}[[\varepsilon]]/(\varepsilon^5)$):**
   - **Exact $\mathbb{Z}_2$ Odd Involution ($\forall \Delta t \in \mathbb{Q}$):** Evaluating over `QQRat.var("dt")` on `step_phi4`, over `QQJet(order=4)` via `IdentityQuantizerBundle` on the full composite $F_{\Delta t} = R \circ S_{\Delta t}$, and around generic finite $\Delta t = \frac{1}{8} + \varepsilon$ certifies $F_{\Delta t}(-X) + F_{\Delta t}(X) = \mathbf{0}$ (`0/1` remainder on all $32$ phase-space coordinates, both at $\nu = 0$ and $\nu > 0$).
   - **Exact $C_4$ Lattice Rotation Equivariance ($\forall \Delta t \in \mathbb{Q}$):** Certifies $S_{\Delta t}(R_{90} X) - R_{90} S_{\Delta t}(X) = \mathbf{0}$ identically in $\mathbb{Q}(\Delta t)$, $F_{\varepsilon}(R_{90} X) - R_{90} F_{\varepsilon}(X) = \mathbf{0}$ in `QQJet(order=4)` (`0/1` across `c0..c4`), and $F_{\frac{1}{8}+\varepsilon}(R_{90} X) - R_{90} F_{\frac{1}{8}+\varepsilon}(X) = \mathbf{0}$ at finite step.
   - **Exact Palindrome Time-Reversibility at $\nu = 0$ ($\forall \Delta t \in \mathbb{Q}$):** Certifies $F_{-\Delta t}(F_{+\Delta t}(X)) - X = \mathbf{0}$ directly in $\mathbb{Q}(\Delta t)$ (and in `QQJet(order=4)` via `IdentityQuantizerBundle`).
   - **Exact Symplectic Energy Conservation at $\nu = 0$:** Certifies anchor conservation in $\mathbb{Q}(\Delta t)$ and exact `0/1` energy conservation $H(F_{\Delta t}(X)) - H(X) = 0/1$ (for both `hamiltonian` and independent continuum $H_{\text{ref}}$) across all jet orders `c0..c4` in `QQJet(order=4)` and around finite $\Delta t = \frac{1}{8} + \varepsilon$.
   - **Exact First Law & SOS Second Law at $\nu > 0$ ($\forall \nu \in \mathbb{Q}_{\ge 0}$):** Evaluating over `QQRat.var("nu")` certifies exact anchor First Law conservation in $\mathbb{Q}(\nu)$ (plus composite readout First Law in `QQJet(order=4)` and around finite $\nu = \frac{1}{4} + \varepsilon$), while $\Delta Q_{\text{bath}}(\nu) = Q_{\text{bath}}' - Q_{\text{bath}} \ge 0/1$ for all $\nu \ge 0$ is certified algebraically in $\mathbb{Q}(\nu)$ at Pólya degree $0$ (`polya_lukacs_sos_deg_0`) via `lukacs_positivity_halfline` on both numerator and denominator (`0/1` polynomial remainder).
3. **Theorem 6 ($p$-adic Odd-Denominator Valuation Invariance Across All $21$ Primes $p \in [101, 199]$):**
   - For every odd prime $p \in \{101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199\}$, evolving prime-perturbed initial conditions $X_0^{(p)} \in \mathbb{Z}_{(2)} \oplus p^{-1}\mathbb{Z}$ for $36$ steps under `DEFAULT_QB = FractionQuantizerBundle(56)` verifies non-negative asymptotic $p$-adic valuation drift ($\min_{t \in \{19,27,35,36\}} v_p(X_t) = -1 \ge \min_{t \in \{1,2,9,18\}} v_p(X_t) = -1 \ge -2$) and stationary $O(1)$ coordinate bit-lengths across all $21$ primes (e.g. $p = 101$: early peak `4532` bits $\to$ late peak `4530` bits; $p = 199$: early peak `4617` bits $\to$ late peak `4625` bits).

---

## 5. Cycle 3 Verification Ledger & Exact Reproduction

All machine-checked symbolic proofs execute end-to-end with zero external packages via:

```bash
/usr/bin/python3 formalist/verify_proofs.py --cycle 3
```

| Verification Check | Target Module / Theorem | Certified Residual / Result | Status |
| :--- | :--- | :--- | :--- |
| **Gate 1 Pure-$\mathbb{Q}$ AST Audit** | `formalist/{__init__,cas,proof_engines,verify_proofs}.py` | `0` banned calls, `0` `while`/`try`, `0` vacuous asserts | **PASS (`0` defects)** |
| **Theorem 1: Padé Boltzmann Algebra** | `architect.battle2.R_22` | SOS $> 0/1$, reciprocity `0/1`, `QQJet(4)` error `0/1` | **PASS (`0/1`)** |
| **Theorem 2: Criticality, Lee-Yang & Ceilings** | `architect.battle1`, `architect.battle2` | $3$ Binder + $3$ $c_V$ Sturm intervals ($< 2^{-40}$), $16/16$ Lee-Yang roots on $|z|=1$, $B_\chi = N$ (`0/1` slack), $B_c$ ($> 98\%$ tight) | **PASS (`0/1`)** |
| **Theorem 3: Finite-$N$ SSB Cholesky SOS** | `architect.battle3` ($\mathbb{Z}_4^2, \mathbb{Z}_8^2$) | $c_4 = 2$, $c_8 = 239/408$, $4{,}012$-monomial `QQMVPoly` remainder `0/1`, $E_c(N) = N\lambda/4$, $M_{\min}^2(E_0) > 0/1$, $B_\phi = 4E_0/\lambda$ | **PASS (`0/1`)** |
| **Theorems 4–6: Symbolic Stepper & $21$ Primes** | `architect.battle3` (`DEFAULT_QB(56)`, 2-stage `[4/3]` Padé) | $t=0$ & $t>0$ `QQJet(3)` `0/1` $H == E_0$, $O(\Delta t_{\text{warm}}^2)$ tangent, $\mathbb{Z}_2$/$C_4$/palindrome `0/1` in $\mathbb{Q}(\Delta t)$, 1st Law `0/1` & 2nd Law `polya_lukacs_sos_deg_0` in $\mathbb{Q}(\nu)$, $21/21$ primes $v_p = -1$ | **PASS (`0/1`, `34s`)** |

