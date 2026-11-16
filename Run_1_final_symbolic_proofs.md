# Machine-Certified Pure-`Q` Symbolic Proofs — Exhibit I (`ultrafinitist-qft` v2.0)

**CAS Engine:** `formalist/cas.py` (`0` Gate 1 AST violations, `0` external dependencies)  
**Executable Proof Verifier:** `formalist/verify_proofs.py` (`/usr/bin/python3 -I formalist/verify_proofs.py` $\implies$ `EXIT 0`, all residuals `0/1`)

---

## 1. Zero-Dependency Pure-`Q` Computer Algebra System (`formalist/cas.py`)

`formalist/cas.py` implements an exact rational computer algebra system strictly over Python `int` and `fractions.Fraction` (`Q`) using `.denominator` negative-integer ring tags (`-1..-5`) instead of `isinstance` dispatch, binary `pow_int`/`pow_ring` instead of variable `**`, and fixed-bound primitive subresultant loops instead of `while` or `try/except`:

- **`QQPoly` & `QQRat`:** Exact univariate polynomials and reduced monic-denominator rational functions over $\mathbb{Q}$ via subresultant primitive polynomial GCD (`gcd_poly`).
- **`QQMVPoly` & `QQMVRat`:** Sparse exact multivariate polynomials and reduced rational functions (`sym_vars`) with exact factor division (`exact_div_or_none`).
- **`QQJet` & `jet_cos_sin`:** Arbitrary-order truncated Taylor series over $\mathbb{Q}$ or `QQMVPoly` (`QQJet(order=k)`), supporting exact valuation-shifted jet division and factorial trigonometric Taylor series without `math.cos`/`math.sin`.
- **`sturm_chain` & `sturm_count_interval`:** Exact signed primitive subresultant Sturm sequences and rational interval root isolation.
- **`cholesky_ldlt_sos` & `verify_sos_cone_certificate`:** Square-root-free symmetric pivoting $P^\top Q P = L D L^\top$ Sum-of-Squares (SOS) certificates and Krivine–Stengle Positivstellensatz cone identity verifiers (`lhs - rhs == cone_sum / positive_den` to `0/1`).
- **`dyadic_valuation`, `odd_part`, `padic_valuation`:** Exact arithmetic valuation verifiers across `PRIMES_101_199` (all 21 primes in $[101, 199]$).

---

## 2. Theorem 1: Gaussian Quotient Torus $V_m$, Symanzik $O(p^4)$ Jets & $W_{\text{flat}}$ `16:16:16` Doubler Elimination

**Machine-Certified in `prove_theorem1_topology_and_no_doublers()` (`0/1` residual):**

1. **Finite Weyl–Schwinger $d = 4$ over $\mathbb{Q}^2$:** Direct $4\times 4$ matrix multiplication over $\mathbb{Q}^2$ certifies $U^4 - I_4 = 0$, $V^4 - I_4 = 0$, $U V - i V U = 0$, $\frac{1}{4}\operatorname{Tr}(\Phi^2) = \frac{1}{2}$, $\frac{1}{4}\operatorname{Tr}(\Phi^4) = \frac{1}{2}$, and $\kappa_W(4) = 2$.
2. **Exact $C_4$ Automorphism on $V_m = \mathbb{Z}[i]/\langle(1+i)^m\rangle$ ($m \in \{6, 7, 8, 9, 10\}$):** For all 5 sizes $N \in \{64, 128, 256, 512, 1024\}$ (both even and odd $m$), `c4_map` ($z \mapsto i z$) satisfies `c4_map[0] == 0`, `c4_map^4 == id`, `c4_map^2 != id`, $c_4(z + 1) = c_4(z) + i$, and $c_4(z + i) = c_4(z) - 1$ identically.
3. **4th-Order Symanzik & Biharmonic Wilson `QQJet` Identities:**
   - Univariate `QQJet(order=6)` in $(1+1)\text{D}$:
     $$K_{\text{Sym}}(\epsilon) = \epsilon^2 - \frac{1}{90}\epsilon^6, \qquad P_{\text{Sym}}(\epsilon) = \epsilon - \frac{1}{30}\epsilon^5, \qquad W_{\text{1D}}(\epsilon) = \epsilon^4 + O(\epsilon^6),$$
     yielding $c_{\text{scalar}}^2 = c_{\text{Dirac}}^2 = \frac{1}{1}$ (whereas Cycle 2 yields $c_{\text{Dirac, Cycle 2}}^2 = \frac{3}{2}$).
   - Bivariate `QQJet(order=4)` over `QQMVPoly(("a", "b"))` at $(p_x, p_y) = (a\epsilon, b\epsilon)$ in $(2+1)\text{D}$:
     $$W_{\text{flat}}(a\epsilon, b\epsilon) = (a^4 - a^2 b^2 + b^4)\epsilon^4, \qquad E_{\text{scalar}}^2(a\epsilon, b\epsilon) = \frac{1}{4} + (a^2 + b^2)\epsilon^2 + 0\epsilon^4, \qquad E_{\text{Dirac}}^2(a\epsilon, b\epsilon) = \frac{1}{4} + (a^2 + b^2)\epsilon^2 + O(\epsilon^4).$$
4. **$W_{\text{flat}}$ `LDL^T` SOS, Sturm BZ Isolation & `16:16:16` Doubler Elimination:**
   - `cholesky_ldlt_sos([[1, -1/2], [-1/2, 1]])` returns positive pivots $(1, \frac{3}{4})$ with reconstruction error `0/1`, proving $W_{\text{flat}} = (L_x - \frac{1}{2}L_y)^2 + \frac{3}{4}L_y^2 \succeq 0$.
   - Writing $P_{\text{Sym}}(p_\mu) = \sin p_\mu (1 + \frac{1}{3}u_\mu)$ with $u_\mu = 1 - \cos p_\mu \in [0, 2]$, `sturm_count_interval(sturm_chain(1 + u/3), 0, 2) == 0`, proving `P_Sym` vanishes **only** at $p_\mu \in \{0, \pi\}$.
   - Evaluating at $(0,0), (\pi,0), (0,\pi), (\pi,\pi)$ gives $W_{\text{flat}} = (0, 16, 16, 16)$ and $E_{\text{Dirac}}^2 = (\frac{1}{4}, 4, 4, 4)$, eliminating all 3 BZ-corner fermion doublers.
5. **Cubic $8 = 2^3$ IR Chiral Suppression:** In `QQJet(order=6)`, $\lim_{\epsilon \to 0} \frac{W_{\text{op}}(\epsilon)/\epsilon}{W_{\text{op}}(\epsilon/2)/(\epsilon/2)} = \frac{8}{1}$.

---

## 3. Theorem 2: Varga Positivity, Positivstellensatz Quadrature Certificates & Global Loop Ceilings

**Machine-Certified in `prove_theorem2_varga_and_quadrature_ceilings()` (`0/1` residual):**

1. **Algebraic Varga Splitting Identities over `QQMVPoly`:**
   Using the abelian shift-operator identities $L_\mu = 2I - A_\mu$ and $L_{2\mu} = 2I - A_{2\mu} = 4I - A_\mu^2$:
   $$\frac{1}{12}(8I - A_1)^2 - (K_{\text{Sym}}^{(1+1)} + m_0^2 I) = (3 - m_0^2) I,$$
   $$\frac{1}{12}(8I - A_x - A_y)^2 - (K_{\text{Sym}}^{(2+1)} + m_0^2 I) = \left(\frac{2}{3} - m_0^2\right) I + \frac{1}{6} A_x A_y.$$
   Verified entrywise on $64 \times 64$ matrices at $m_0^2 \in \{\frac{1}{4}, \frac{2}{3}\}$: $R \ge 0$ and $G_{xy} > 0/1$ strictly, whereas unfactored biharmonic $K_{\text{hd}} = \frac{1}{4}\sum_\mu L_\mu^2$ exhibits negative ringing $G_{\min} < 0/1$.
2. **Six Symbolic Positivstellensatz SOS Cone Identities (`verify_sos_cone_certificate` over `QQMVRat`):**
   On $[0, M]$ with $u \ge 0$, $M - u \ge 0$, $m^2 > 0$:
   - **Global Ceiling (`0/1` slack at $u = 0$):** $\frac{1}{m^2} - \frac{1}{u + m^2} = \frac{u}{m^2(u + m^2)} \ge 0$, giving $\tau_1 \le \frac{1}{m_0^2}$ and $s_2 = \sum_x G_{0x}^3 \le G_{00}\sum_x G_{0x}^2 = \tau_1 \tau_2 \le \frac{1}{m_0^6}$ (saturated with `0/1` slack at zero hopping $K = 0$).
   - **Secant Upper Bound:** $\left(\frac{1}{m^2} - \frac{u}{m^2(M + m^2)}\right) - \frac{1}{u + m^2} = \frac{u(M - u)}{m^2(M + m^2)(u + m^2)} \ge 0$.
   - **Jensen Lower Bound:** $\frac{1}{u + m^2} - \left(\frac{1}{s_1 + m^2} - \frac{u - s_1}{(s_1 + m^2)^2}\right) = \frac{(u - s_1)^2}{(s_1 + m^2)^2(u + m^2)} \ge 0$.
   - **Degree-2 Left Gauss–Radau Upper Bound:** $P_{\text{LR},2}(u) - \frac{1}{u + m^2} = \frac{u(u - x_1)^2}{m^2(x_1 + m^2)^2(u + m^2)} \ge 0$.
   - **Degree-2 Right Gauss–Radau Lower Bound:** $\frac{1}{u + m^2} - P_{\text{RR},2}(u) = \frac{(M - u)(u - y_1)^2}{(M + m^2)(y_1 + m^2)^2(u + m^2)} \ge 0$.
   - **Degree-3 Gauss–Lobatto Upper Bound:** $P_{\text{GL},3}(u) - \frac{1}{u + m^2} = \frac{u(M - u)(u - z_1)^2}{m^2(M + m^2)(z_1 + m^2)^2(u + m^2)} \ge 0$.
3. **54th-Order Midpoint Brackets & Fair-Fight $1/N$ Parity:**
   Writing $K = \frac{M}{2}I - B_{\text{mid}}$ with $\|B_{\text{mid}}\|_{\text{op}} \le \frac{M}{2}$, the zero-mode-projected Neumann series over $\mathbb{Q}$ ($\rho = \frac{M/2}{M/2 + m_0^2} = \frac{32}{35}$ in 1D, $\frac{64}{67}$ in 2D) brackets $\tau_1$ to **`99.82%` tightness in $(1+1)\text{D}$** and **`98.53%` tightness in $(2+1)\text{D}$**, containing both `cycle3` and `baseline` under identical $1/N$ normalization.
4. **Dyadic-Shell Proof of Continuum $a \to 0$ UV Divergence over $\mathbb{Q}$:**
   Partitioning $[m_{\text{phys}}, 2^J m_{\text{phys}}]$ into $J$ dyadic shells $[2^j m_{\text{phys}}, 2^{j+1} m_{\text{phys}}]$ ($\Lambda = 2^J m_{\text{phys}} = \pi/a \to \infty$), exact rational lower sums certify $\int_{m_{\text{phys}}}^\Lambda \frac{k\,dk}{k^2 + m_{\text{phys}}^2} \ge \sum_{j=0}^{J-1} \frac{4^j}{4^{j+1} + 1} \ge \frac{J}{5} \to +\infty$ and $\int_{m_{\text{phys}}}^\Lambda \frac{k^2\,dk}{k^2 + m_{\text{phys}}^2} \ge m_{\text{phys}}\sum_{j=0}^{J-1} \frac{8^j}{4^{j+1} + 1} \ge \frac{2^J - 1}{5} m_{\text{phys}} \to +\infty$.

---

## 4. Theorem 3 & Theorem 6: Double-Cayley + Vieta `0/1` Conservation, Block-Spin RG & 21-Prime Invariance

**Machine-Certified in `prove_theorem3_and_theorem6_dynamics_and_rg()` (`0/1` residual):**

1. **`QQJet(order=3)` Taylor PDE Consistency & `0/1` Manifold Invariants:**
   Evaluating `step_wavepacket` and `vieta_readout` over `QQJet(order=3)` under `IDENTITY_QUANTIZER_BUNDLE` in both $(1+1)\text{D}$ and $(2+1)\text{D}$ certifies componentwise to `0/1`:
   $$\Psi_{\text{step}}(\epsilon) = \Psi_0 - i\epsilon H_D \Psi_0 - \frac{1}{2}\epsilon^2 H_D^2 \Psi_0 + O(\epsilon^3), \qquad \Psi_{\text{Vieta}}(\epsilon) = \Psi_0 - i\epsilon H_D \Psi_0 - \frac{1}{2}\epsilon^2 H_D^2 \Psi_0 + O(\epsilon^3),$$
   together with $\|\Psi_{\text{Vieta}}(\epsilon)\|_2^2 - \|\Psi_0\|_2^2 \equiv 0/1$ and $\langle \Psi_{\text{Vieta}}(\epsilon), H_D \Psi_{\text{Vieta}}(\epsilon)\rangle - E_0 \equiv 0/1$ in `QQJet(order=3)` AND over `Q` at coarse $\Delta t \in \{\frac{1}{4}, \frac{1}{8}\}$.
2. **`0/1` Palindrome Time-Reversal & Finite Lightcone:**
   $F_{-\Delta t}^{10}(F_{+\Delta t}^{10}(\Psi_0)) - \Psi_0 = 0/1$ identically under `RuntimeDyadicQuantizerBundle`, and a single step on $\mathbb{Z}_{64}$ from a point impulse at $x = 0$ produces `0/1` amplitude at all 39 sites outside $|x| \le 12$.
3. **Exact `(1+i)`-Adic Block-Spin Decimation ($1024 \to 256$, $512 \to 128$, $256 \to 64$):**
   Certifies to `0/1` across $m_{\text{coarse}} \in \{6, 7, 8\}$ (both even and odd $m$): column sums $= 1/1$ (exact spinor charge conservation), unit norm $\|\Psi_{\text{coarse}}\|_2^2 = \|\Psi_{\text{fine}}\|_2^2 = 1/1$, exact $C_4$ equivariance $R_{\text{block}} C_4^{(m+2)} = C_4^{(m)} R_{\text{block}}$, and `QQJet(order=4)` coarse symbol $K_{\text{RG}}^{(m)}(P_x\epsilon, P_y\epsilon) = (P_x^2 + P_y^2)\epsilon^2 + O(\epsilon^4)$ ($c_{\text{RG}}^2 = 1/1$).
4. **Theorem 6 ($p$-Adic Odd-Denominator Invariance Across All 21 Primes `101..199`):**
   Because `RuntimeDyadicQuantizerBundle` emits values $m / 2^e \in \mathbb{Z}[1/2]$ and $\Delta t \in \{\frac{1}{4}, \frac{1}{8}\} \subset \mathbb{Z}[1/2]$, every step increment lies in $\mathbb{Z}[1/2]$. Hence for any initial wavepacket perturbed by $\frac{1}{4p}$ with $p \in \{101, 103, \dots, 199\}$, $v_p(\Psi_n) \ge -1$ and $\operatorname{odd\_part}(\operatorname{den}(\Psi_{n,k})) \in \{1, p\}$ hold for all $n \ge 0$, machine-verified across all 21 primes with `0` exceptions.
