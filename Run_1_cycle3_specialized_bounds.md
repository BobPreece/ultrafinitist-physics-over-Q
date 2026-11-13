# Constructive Formalist Cycle 3 Specialization: Exact Positivstellensatz Certificates, Squared-$M$-Matrix Positivity, and Rational Loop Envelopes

**Author:** Agent 4 (`formalist`, Bounds & Proof Verifier)  
**Target Winning Architectures (`architect/cycle3_candidates.py`):**
- **`(1+1)D`:** `cycle3_symanzik_biharmonic_1p1D` on $\mathbb{Z}_N$ ($N \in \{64, 128, 256, 512, 1024\}$)
- **`(2+1)D`:** `cycle3_gaussian_symanzik_biharmonic_2p1D` on the $D_4$-symmetric Gaussian-integer quotient ring $\mathbb{Z}[i]/\langle(1+i)^m\rangle$ ($N = 2^m \in \{64, 128, 256, 512, 1024\}$)

All 20 formal multivariate polynomial and Positivstellensatz certificates below are algebraically proved with **zero floating-point error** in `formalist/symbolic_cas_and_bounds.py` via both pure-Python `SymPoly`/`SymRat` over $\mathbb{Q}$ (`fractions.Fraction`) and native `sympy` (`prove_all_identities_sympy()`), with `0` AST continuum/counterterm contraband violations.

---

## 1. Theorem C3.1 — Exact Weyl-Schwinger ($d = 4$) Operator Norms & Many-Body Boundedness

For local Hilbert dimension $d = 4$ ($d \equiv 0 \pmod 4$) on $\mathcal{H} = (\mathbb{C}^4)^{\otimes N}$ ($\dim \mathcal{H} = 4^N$):
1. **Local Field Spectrum & Moments:** $\Phi_x = \frac{U_x - U_x^\dagger}{2i}$ has exact spectrum $\{\sin(\frac{\pi q}{2})\}_{q=0}^3 = \{0, +1, 0, -1\}$, exact operator norm
   $$\|\Phi_x\|_{\mathrm{op}} = B(4) = 1 \in \mathbb{Q}, \qquad B(4)^4 = 1,$$
   and exact rational normalized moments $\frac{1}{4}\operatorname{Tr}(\Phi^2) = \frac{1}{2}$, $\frac{1}{4}\operatorname{Tr}(\Phi^4) = \frac{1}{2}$, giving exact Weyl kurtosis $\kappa_W(4) = \frac{1/2}{(1/2)^2} = 2 \in \mathbb{Q}$.
2. **Exact Rational Couplings (`lambda_0 = 1/2`, `m_0^2 = 1/4`):**
   $$\lambda_{\mathrm{site}} = \lambda_0 \frac{k_{\mathrm{eff}}}{\kappa_W(4)} = \begin{cases} \frac{1}{2}\cdot\frac{5/2}{2} = \frac{5}{8}, & (1+1)\mathrm{D}\ (k_{\mathrm{eff}} = \frac{5}{2}), \\[4pt] \frac{1}{2}\cdot\frac{5}{2} = \frac{5}{4}, & (2+1)\mathrm{D}\ (k_{\mathrm{eff}} = 5), \end{cases} \qquad \lambda_{\mathrm{mode}}(N) = \frac{\lambda_{\mathrm{site}}}{N}.$$
3. **Exact Tensor-Product Saturation & Total Hamiltonian Norm:** On the product clock eigenstate $|1\rangle^{\otimes N} \in (\mathbb{C}^4)^{\otimes N}$,
   $$\|H_{\mathrm{int}}\|_{\mathrm{op}} = \left\|\lambda_{\mathrm{site}} \sum_{x=0}^{N-1} \Phi_x^4\right\|_{\mathrm{op}} = N \lambda_{\mathrm{site}} B(4)^4 = \begin{cases} \frac{5}{8} N < \infty, & (1+1)\mathrm{D}, \\ \frac{5}{4} N < \infty, & (2+1)\mathrm{D}, \end{cases}$$
   and $\|H_0 + H_{\mathrm{int}}\|_{\mathrm{op}} \le N\big(2 + \frac{1}{2}(M + m_0^2) + \lambda_{\mathrm{site}}\big) = \frac{65}{12}N$ in `(1+1)D` and $\frac{209}{24}N$ in `(2+1)D`.

---

## 2. Theorem C3.2 — Positivstellensatz Operator Certificates for Symanzik & Equal-Corner Biharmonic Wilson

Let $T_\mu$ be the commuting shift unitaries ($T_\mu^\dagger = T_\mu^{-1}$, $[T_x, T_y] = 0$), $A_\mu = T_\mu + T_\mu^\dagger$ ($-2 I \preceq A_\mu \preceq 2 I$), $L_\mu = 2 I - A_\mu = (I - T_\mu)^\dagger(I - T_\mu)$ ($0 \preceq L_\mu \preceq 4 I$), and $L_{2\mu} = 2 I - (T_\mu^2 + (T_\mu^\dagger)^2) = 4 I - A_\mu^2 = 4 L_\mu - L_\mu^2$.

1. **Basis-Independent Symanzik Positivity (`Certificate 15`):**
   $$K_{\mathrm{Sym},\mu} = \frac{4}{3} L_\mu - \frac{1}{12} L_{2\mu} \equiv L_\mu + \frac{1}{12} L_\mu^2 = L_\mu^{1/2}\!\left(I + \frac{1}{12} L_\mu\right)\!L_\mu^{1/2} \succeq 0.$$
   Hence $K_{\mathrm{Sym}}^{(1+1)} = K_{\mathrm{Sym},x} \succeq 0$ and $K_{\mathrm{Sym}}^{(2+1)} = K_{\mathrm{Sym},x} + K_{\mathrm{Sym},y} \succeq 0$ WITHOUT invoking any Fourier transform.
2. **Exact Sharp Spectral Caps $M = \frac{16}{3}$ and $M = \frac{32}{3}$ (`Certificate 16`):**
   For any eigenvalue $\ell \in [0, 4]$ of $L_\mu$:
   $$\frac{16}{3} - \left(\ell + \frac{1}{12}\ell^2\right) = \frac{1}{12}(4 - \ell)(16 + \ell) = \frac{4}{3}(4 - \ell) + \frac{1}{12}\ell(4 - \ell) \ge 0.$$
   Equality holds at the staggered UV mode $\ell = 4$ ($p_\mu = \pi$), proving $\|K_{\mathrm{Sym}}^{(1+1)}\|_2 = \frac{16}{3}$ and $\|K_{\mathrm{Sym}}^{(2+1)}\|_2 = \frac{32}{3}$ (matching the Gershgorin row-sum bound $\frac{5}{2} + \frac{8}{3} + \frac{1}{6} = \frac{16}{3}$ and $5 + \frac{16}{3} + \frac{1}{3} = \frac{32}{3}$ with **zero slack**).
3. **Equal-Corner SOS Biharmonic Wilson Operator (`Certificate 17`):**
   In `(2+1)D`, $W_{\mathrm{flat}} = L_x^2 - L_x L_y + L_y^2$ admits the exact commuting operator Sum-of-Squares decomposition and $[0, 4]^2$ global cap certificate:
   $$W_{\mathrm{flat}} = \left(L_x - \frac{1}{2}L_y\right)^2 + \frac{3}{4}L_y^2 \succeq 0, \qquad 16 I - W_{\mathrm{flat}} = (4 I - L_x)(4 I - L_y) + L_x(4 I - L_x) + L_y(4 I - L_y) \succeq 0.$$
   Every cone summand in $16 - W_{\mathrm{flat}}(\ell_x, \ell_y)$ is non-negative on $[0, 4]^2$ and vanishes **if and only if** $(\ell_x, \ell_y) \in \{(4, 0), (0, 4), (4, 4)\}$, proving simultaneously that $0 \preceq W_{\mathrm{flat}} \preceq 16 I$ and that all three UV doubler corners $(p_x, p_y) \in \{(\pi, 0), (0, \pi), (\pi, \pi)\}$ receive the exact isotropic lift $m_0 + \frac{3}{32}\cdot 16 = \frac{1}{2} + \frac{3}{2} = 2$ (`16 : 16 : 16`), while $W_{\mathrm{flat}} = O(|p|^4)$ preserves $c_{\mathrm{Dirac}}^2 = c_{\mathrm{scalar}}^2 = 1$ identically in the IR.

---

## 3. Theorem C3.3 — Varga Squared-$M$-Matrix Regular Splitting ($G_{xy} > 0$ Entrywise Despite Repulsive 2-Hop Stencil)

Because $K_{\mathrm{Sym}}$ has positive off-diagonal entries $+\frac{1}{12}$ at axial distance $2$, $K_{\mathrm{Sym}} + m_0^2 I$ is **not** an $M$-matrix. Nevertheless, let $A_{\mathrm{nn}}$ be the nearest-neighbor adjacency matrix ($A_x$ in `(1+1)D`; $A_x + A_y$ in `(2+1)D`) and define $B = 8 I - A_{\mathrm{nn}}$.

1. Since $A_{\mathrm{nn}} \ge 0$ has row sum $2 < 8$ (`1+1D`) or $4 < 8$ (`2+1D`), $B = 8\big(I - \frac{1}{8}A_{\mathrm{nn}}\big)$ is a strictly diagonally dominant symmetric $M$-matrix with strictly positive inverse $B^{-1} = \frac{1}{8}\sum_{s=0}^\infty (\frac{1}{8}A_{\mathrm{nn}})^s > 0$ (entrywise on connected graphs). Therefore $P := \frac{1}{12} B^2$ is invertible with strictly positive inverse $P^{-1} = 12 (B^{-1})^2 > 0$.
2. Using $K_{\mathrm{Sym},x} = \frac{7}{3}I - \frac{4}{3}A_x + \frac{1}{12}A_x^2$ (since $L_x = 2I - A_x$ and $L_{2x} = 4I - A_x^2$), we obtain the exact polynomial operator identities (`Certificates 18 & 19`):
   $$R^{(1+1)} := \frac{1}{12}(8 I - A_x)^2 - \big(K_{\mathrm{Sym}}^{(1+1)} + m_0^2 I\big) \equiv (3 - m_0^2) I \ge 0 \qquad (\forall\, 0 < m_0^2 \le 3),$$
   $$R^{(2+1)} := \frac{1}{12}(8 I - A_x - A_y)^2 - \big(K_{\mathrm{Sym}}^{(2+1)} + m_0^2 I\big) \equiv \left(\frac{2}{3} - m_0^2\right) I + \frac{1}{6} A_x A_y \ge 0 \qquad \left(\forall\, 0 < m_0^2 \le \frac{2}{3}\right).$$
3. At Cycle 3's bare mass $m_0^2 = \frac{1}{4}$, $R^{(1+1)} = \frac{11}{4} I \ge 0$ and $R^{(2+1)} = \frac{5}{12} I + \frac{1}{6} A_x A_y \ge 0$ entrywise! Since $K_{\mathrm{Sym}} + m_0^2 I = P - R$ is symmetric positive definite with $P^{-1} > 0$ and $R \ge 0$, Varga's Regular Splitting Theorem guarantees $\rho(P^{-1} R) < 1$ and convergent entrywise-positive Neumann expansion:
   $$G = (K_{\mathrm{Sym}} + m_0^2 I)^{-1} = \sum_{\ell=0}^\infty (P^{-1} R)^\ell P^{-1} > 0 \quad \text{entrywise on all } N \times N \text{ pairs } (x, y)!$$

---

## 4. Theorem C3.4 — Exact Rational Moments, Wrapping Identity & High-Order Loop Bounds

### 4.1 Closed-Form Normalized Spectral Moments $s_r = \frac{1}{N}\operatorname{Tr}(K_{\mathrm{Sym}}^r)$
- **`(1+1)D` (`N in {64, 128, 256, 512, 1024}`):** Since a range-2 stencil walk of $\le 4$ steps spans at most $\pm 8$ sites and $N \ge 64 > 16$, the normalized moments $s_1..s_{14}$ are **strictly independent of $N$**:
  $$s_1 = \frac{5}{2}, \qquad s_2 = \frac{707}{72}, \qquad s_3 = \frac{6233}{144}, \qquad s_4 = \frac{77209}{384}.$$
- **`(2+1)D` (`N in {64, 128, 256, 512, 1024}`):**
  $$s_1 = 5, \qquad s_2 = \frac{1157}{36}, \qquad s_3 = \frac{8419}{36}, \qquad s_4(N) = \begin{cases} \dfrac{4785749}{2592}, & N = 64, \\[8pt] \dfrac{3190499}{1728}, & N \in \{128, 256, 512, 1024\}. \end{cases}$$
  **Exact Toroidal Wrapping Proof for $N = 64$:** On $\mathbb{Z}[i]/\langle(1+i)^6\rangle$, $(1+i)^6 = -8i$ generates the $8 \times 8$ square lattice $8\mathbb{Z}[i]$. Four identical axial 2-hop steps $(\pm 2, 0)^4$ or $(0, \pm 2)^4$ reach $(\pm 8, 0) \equiv (0, 0)$ and $(0, \pm 8) \equiv (0, 0)$, contributing exactly $4 \times (+\frac{1}{12})^4 = \frac{4}{20736} = \frac{1}{5184} = \frac{4785749}{2592} - \frac{3190499}{1728}$! For $N = 2^m \ge 128$ ($m \ge 7$), $(1+i)^7 = 8(1 - i)$ has shortest non-zero $L^1$-displacement $16$, requiring $8$ range-2 steps, so $s_1, \dots, s_7$ are **identically constant** across all $N \in \{128, 256, 512, 1024\}$.

### 4.2 High-Order Positivstellensatz Quadrature & Sunset Bounds (`s_0 .. s_{14}`)
Using machine-verified shift transitivity (`verify_shift_transitivity_and_moments`), `formalist` extracts exact `Fraction` moments $s_0, \dots, s_{14}$ and evaluates:
1. **Order-$r$ Gauss-Lobatto & Left-Gauss-Radau Upper Bounds (`Lobatto_{2r+1}`, $r = 6 \implies \deg 13$; `Left_Radau_{2r}`, $r = 7 \implies \deg 14$):**
   $$P_{2r+1}^{\mathrm{Lob}}(\mu; c) - \frac{1}{\mu + m_0^2} = \frac{\mu(M - \mu)\,Q_r(\mu + m_0^2)^2}{m_0^2 (M + m_0^2)(\mu + m_0^2)} \ge 0 \quad \forall \mu \in [0, M], \quad Q_r(0) = 1.$$
2. **Zero-Mode-Peeled Right-Gauss-Radau Lower Bound (`ZeroMode_Right_Radau_{2r}`, $r = 7 \implies \deg 14$):**
   Separates the exact zero mode $\mu_0 = 0$ (weight $\frac{1}{N}$, value $\frac{1}{m_0^2}$) and applies $\frac{1}{\mu + m_0^2} - P_{2r}^{\mathrm{R}}(\mu; d) = \frac{(M - \mu) Q_r(\mu + m_0^2)^2}{(M + m_0^2)(\mu + m_0^2)} \ge 0$ on the remaining $N - 1$ modes.
3. **Squared Lobatto Bubble Bound (`Lobatto_7^2`, $\deg 14$) + Spectral-Convexity Off-Diagonal Gap + Dyson-Jensen Sunset Lower Bound:**
   - Since $P_7^{\mathrm{Lob}}(\mu; c) \ge (\mu + m_0^2)^{-1} > 0$ on $[0, M]$, $\tau_2 = \frac{1}{N}\operatorname{Tr}(G^2) \le \mathbb{E}\big[(P_7^{\mathrm{Lob}}(\mu; c))^2\big]$.
   - For every $y \neq 0$, let $v = \frac{1}{\sqrt{2}}(e_0 - e_y)$ ($\|v\|_2 = 1$). Convexity of $(\mu + m_0^2)^{-1}$ yields $G_{00} - G_{0y} = v^\top G v \ge (K_{00} - K_{0y} + m_0^2)^{-1} \ge \Delta_{\mathrm{gap}}$ ($\Delta_{\mathrm{gap}} = \frac{12}{49}$ in `(1+1)D`, $\frac{12}{79}$ in `(2+1)D`), tightening $\max_{y \neq 0} G_{0y} \le \tau_1^{\mathrm{up}} - \Delta_{\mathrm{gap}}$.
   - Since $G_{xy} > 0$ (Theorem C3.3) and $K_{0w} > 0$ on 2-hop sites $w \in \mathcal{N}_2(0)$, Dyson's identity $(K_{00} + m_0^2) G_{00} - 1 = \frac{4}{3}\sum_{v \in \mathcal{N}_1(0)} G_{0v} - \frac{1}{12}\sum_{w \in \mathcal{N}_2(0)} G_{0w} < \frac{4}{3}\sum_{v \in \mathcal{N}_1(0)} G_{0v}$ combined with Jensen's inequality on $|\mathcal{N}_1(0)| = 2 D_s$ gives the strictly positive rational sunset lower bound:
     $$\mathcal{S}_2^{(p=0)}(N) \ge (\tau_1^{\mathrm{low}})^3 + \frac{1}{(2 D_s)^2}\left[\frac{3}{4}\big((k_{\mathrm{eff}} + m_0^2)\tau_1^{\mathrm{low}} - 1\big)\right]^3 > 0.$$
4. **Midpoint-Shifted Exponential Stabilization Without Counterterms:**
   Writing $K_{\mathrm{Sym}} + m_0^2 I = \alpha I - B_{\mathrm{mid}}$ with $\alpha = \frac{M}{2} + m_0^2$, $\|B_{\mathrm{mid}}\|_2 = \frac{M}{2}$, and contraction factor $\rho_{\mathrm{mid}} = \frac{M}{M + 2m_0^2} < 1$ ($\rho_{\mathrm{mid}} = \frac{32}{35}$ in `(1+1)D`, $\frac{64}{67}$ in `(2+1)D`), all local walk and 3-banana diagrams up to $L_{\mathrm{match}} = \lceil L_{\mathrm{wrap}}/2 \rceil - 1$ steps match identically across lattice sizes $N_1 < N_2$, proving exponential convergence $|R_m(N_1) - R_m(N_2)| = O(\rho_{\mathrm{mid}}^{L_{\mathrm{match}}+1})$ with **zero counterterm subtraction**, in stark contrast to the unbounded continuum UV divergences proved in Certificates 13–14.

---

## 5. Certified Tightness Summary Across `N in [64, 128, 256, 512, 1024]` (`m_0^2 = 1/4`, `d = 4`, `lambda_0 = 1/2`)

| Dimension | $N$ | $\tau_1^{\mathrm{low}}$ | $\tau_1^{\mathrm{actual}}$ | $\tau_1^{\mathrm{up}}$ | $\tau_1$ Tightness | $\mathcal{S}_2^{\mathrm{low}}$ | $\mathcal{S}_2^{\mathrm{actual}}$ | $\mathcal{S}_2^{\mathrm{up}}$ | $R_m^{(2)\,\mathrm{low}}$ | $R_m^{(2)\,\mathrm{actual}}$ | $R_m^{(2)\,\mathrm{up}}$ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `(1+1)D` | 64 | 0.922499 | 0.924989 | 0.927922 | **99.684%** | 1.16791 | 1.35590 | 1.74707 | 1.69816 | 1.80314 | 1.85576 |
| `(1+1)D` | 128 | 0.922265 | 0.924989 | 0.927922 | **99.684%** | 1.16683 | 1.35590 | 1.74707 | 1.69786 | 1.80314 | 1.85604 |
| `(1+1)D` | 256 | 0.922161 | 0.924989 | 0.927922 | **99.684%** | 1.16635 | 1.35590 | 1.74707 | 1.69773 | 1.80314 | 1.85617 |
| `(1+1)D` | 512 | 0.922112 | 0.924989 | 0.927922 | **99.684%** | 1.16612 | 1.35590 | 1.74707 | 1.69767 | 1.80314 | 1.85622 |
| `(1+1)D` | 1024 | 0.922088 | 0.924989 | 0.927922 | **99.684%** | 1.16601 | 1.35590 | 1.74707 | 1.69764 | 1.80314 | 1.85625 |
| `(2+1)D` | 64 | 0.346443 | 0.346501 | 0.348050 | **99.555%** | 0.05606 | 0.07032 | 0.11911 | 1.74204 | 1.79300 | 1.81173 |
| `(2+1)D` | 128 | 0.338854 | 0.339323 | 0.343869 | **98.678%** | 0.05137 | 0.06111 | 0.11578 | 1.72653 | 1.78465 | 1.80616 |
| `(2+1)D` | 256 | 0.336426 | 0.338173 | 0.343529 | **98.441%** | 0.04994 | 0.05975 | 0.11555 | 1.72070 | 1.78319 | 1.80680 |
| `(2+1)D` | 512 | 0.335453 | 0.338080 | 0.343528 | **98.414%** | 0.04937 | 0.05964 | 0.11555 | 1.71827 | 1.78307 | 1.80739 |
| `(2+1)D` | 1024 | 0.335006 | 0.338077 | 0.343528 | **98.413%** | 0.04912 | 0.05964 | 0.11555 | 1.71715 | 1.78307 | 1.80766 |
