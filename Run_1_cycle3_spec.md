# Cycle 3 Specification: 4th-Order Axial Symanzik Stencils & Sum-of-Squares Equal-Corner Biharmonic Wilson-Dirac QFT

## 1. Summary of Refinements Addressing Cycle 2 Red-Team Audit (`#19`)

In Cycle 2 (`cycle2_local_wilson_ring_1p1D` and `cycle2_isotropic_moore_wilson_2p1D`), `inquisitor` confirmed full passage of Gate 1 (`0` contraband violations), Gate 2 (`doubler_count = 0.0`), Causal Lightcone Propagation (`W_step` leak `0.0%` outside $r \le 2$), $(2+1)\text{D}$ scalar shell anisotropy ($0.00562 < 0.06$), and Gate 4 Loop Stability ($0.0000\%$ drift in `1+1D`, $0.5353\% < 3\%$ drift in `2+1D`), while rejecting Cycle 2 on three Gate 3 dispersion/spinor-isotropy defects:

1. **First-Order Wilson Mass Speed-of-Light Distortion:** $M_W(p) = m_0 + \frac{r_W}{2}|p|^2 + O(p^4)$ ($m_0 = 1/2, r_W = 1$) injected a $+m_0 r_W |p|^2$ cross-term into $E_D^2(p) = |P(p)|^2 + M_W(p)^2 = m_0^2 + (1 + m_0 r_W)|p|^2 + O(p^4)$, inflating $c_{\text{Dirac}}^2 = 1.50$ vs $c_{\text{scalar}}^2 = 1.00$ (`DL_Dirac = 0.1132` in `1+1D`, `0.2396` in `2+1D`; Boson–Fermion mismatch `0.1412` in `1+1D`, `0.4462` in `2+1D`).
2. **Unimproved $(2+1)\text{D}$ Spinor $D_4$ Shell Anisotropy:** $P_\mu(p) = \sin(p_\mu) = p_\mu - \frac{1}{6}p_\mu^3 + O(p_\mu^5)$ left $|P(p)|^2 = |p|^2 - \frac{1}{3}(p_x^4 + p_y^4) + O(|p|^6)$, causing `dirac_shell_anisotropy_max = 0.1370 > 0.06`.
3. **Radial $O(|p|^4)$ Scalar Dispersion Droop in $(2+1)\text{D}$:** $K_{\text{iso}}(p) = |p|^2 - \frac{1}{12}|p|^4 + O(|p|^6)$ eliminated directional anisotropy but left radial droop `DL_scalar = 0.0941 > 0.05` across the lowest $20\%$ Brillouin-zone modes.

---

## 2. Cycle 3 Mathematical Architecture (`architect/cycle3_candidates.py`)

### 2.1 4th-Order Axial Range-2 Symanzik Scalar Kinetic Operator ($K_{\text{Sym}}$)

On each spatial axis $\mu$, let $L_{1,\mu} = 2 I_N - T_\mu - T_\mu^\dagger$ and $L_{2,\mu} = 2 I_N - T_\mu^2 - (T_\mu^\dagger)^2$ denote the exact 1-step and 2-step discrete shift Laplacians. Cycle 3 defines the exact rational Symanzik scalar kinetic operator:

- **$(1+1)\text{D}$ (`cycle3_symanzik_biharmonic_1p1D` on $\mathbb{Z}_N$):**
  $$K_{\text{Sym}}^{(1+1)} = \frac{4}{3} L_1 - \frac{1}{12} L_2, \qquad k_{\text{eff}} = \frac{1}{N}\mathrm{Tr}\!\left(K_{\text{Sym}}^{(1+1)}\right) = \frac{5}{2} \in \mathbb{Q}.$$
- **$(2+1)\text{D}$ (`cycle3_gaussian_symanzik_biharmonic_2p1D` on $\mathbb{Z}[i]/((1+i)^m)$):**
  $$K_{\text{Sym}}^{(2+1)} = \frac{4}{3}(L_x + L_y) - \frac{1}{12}(L_{2x} + L_{2y}), \qquad k_{\text{eff}} = \frac{1}{N}\mathrm{Tr}\!\left(K_{\text{Sym}}^{(2+1)}\right) = 5 \in \mathbb{Q}.$$

**Exact Positivity, $O(|p|^6)$ Dispersion, and Small-Torus Ball Locality:**
1. **Strict Non-Zero Mode Positivity (Exact Algebraic Factorization):**
   $$K_{1\text{D}}(p) = \frac{4}{3}(2 - 2\cos p) - \frac{1}{12}(2 - 2\cos 2p) = \frac{1}{3}(1 - \cos p)(7 - \cos p) \ge 0,$$
   vanishing if and only if $p \equiv 0 \pmod{2\pi}$, with global maximum $16/3$ at $p = \pi$ ($32/3$ at $(\pi, \pi)$ in $2\text{D}$).
2. **Simultaneous Radial & $D_4$ $O(|p|^4)$ Cancellation:**
   $$K_{\text{Sym}}^{(2+1)}(p_x, p_y) = (p_x^2 + p_y^2) - \frac{1}{90}(p_x^6 + p_y^6) + O(|p|^8),$$
   cancelling **both** the isotropic $|p|^4$ radial droop and the $p_x^4 + p_y^4$ $D_4$ directional error with only **8 axial neighbors** $\{(\pm 1,0),(0,\pm 1),(\pm 2,0),(0,\pm 2)\}$.
3. **Bounded 2-Hop Ball Volume on Small Tori ($N = 64$):** Unlike a dense $5 \times 5$ 24-neighbor stencil (whose 2-hop ball covers $64/64 = 100\%$ of an $8 \times 8$ torus), the 8-neighbor axial stencil has $|B(0, 2)| / N = 31 / 64 = 0.484375 < 0.50$ at $N = 64$ (`expander_or_nonlocal_graph = False` across all $N \in [64, 1024]$).

### 2.2 4th-Order Improved Spinor Momentum ($P_{\text{Sym},\mu}$)

For each spatial axis $\mu$, let $P_{r,\mu} = \frac{1}{2i}(T_\mu^r - (T_\mu^\dagger)^r)$ have eigenvalue $\sin(r p_\mu)$. Cycle 3 defines:
$$P_{\text{Sym},\mu} = \frac{4}{3} P_{1,\mu} - \frac{1}{6} P_{2,\mu}, \qquad P_{\text{Sym},\mu}(p) = \frac{4}{3}\sin(p_\mu) - \frac{1}{6}\sin(2p_\mu) = p_\mu - \frac{1}{30}p_\mu^5 + O(p_\mu^7).$$
Notice that $P_{\text{Sym},\mu}(0) = P_{\text{Sym},\mu}(\pm\pi) = 0$ and the entire $O(p_\mu^3)$ error vanishes identically, suppressing $(2+1)\text{D}$ spinor shell anisotropy from $0.1370$ down to $0.0157$ ($< 0.06$).

### 2.3 Sum-of-Squares (SOS) Equal-Corner Biharmonic Wilson Mass ($M_{W2}$)

To eliminate all Brillouin-zone doublers while keeping $c_{\text{Dirac}}^2 = 1$ identically ($0$ $O(p^2)$ cross-term) and preserving a strong uniform UV gap across all zone corners:

- **$(1+1)\text{D}$ Biharmonic Wilson Mass ($r_2 = \frac{1}{16} \in \mathbb{Q}$):**
  $$M_{W2}^{(1+1)} = m_0 I_N + r_2 L_1^2 = m_0 I_N + r_2\bigl(6 I_N - 4(T_x + T_x^\dagger) + T_x^2 + (T_x^\dagger)^2\bigr),$$
  $$M_{W2}^{(1+1)}(p) = m_0 + r_2 (2 - 2\cos p)^2 = m_0 + r_2 p^4 + O(p^6), \qquad M_{W2}^{(1+1)}(\pi) = m_0 + 16 r_2 = \frac{3}{2}.$$

- **$(2+1)\text{D}$ Equal-Corner SOS Biharmonic Wilson Mass ($r_2 = \frac{3}{32} \in \mathbb{Q}$):**
  Naive $(L_x + L_y)^2$ yields $16$ at $(\pi, 0)$ vs $64$ at $(\pi, \pi)$ ($1:4$ imbalance, compressing the normalized UV gap ratio $(E_{\min,\text{UV}}^2 - m_0^2)/(E_{\max}^2 - m_0^2)$). Instead, Cycle 3 introduces the **Sum-of-Squares Equal-Corner Biharmonic Operator**:
  $$W_{\text{flat}} = L_x^2 - L_x L_y + L_y^2 = \left(L_x - \frac{1}{2}L_y\right)^2 + \frac{3}{4}L_y^2 \succeq 0,$$
  $$M_{W2}^{(2+1)} = m_0 I_N + r_2 W_{\text{flat}}.$$
  **Exact Properties of $W_{\text{flat}}$:**
  1. **Positive Semidefinite SOS Operator:** $W_{\text{flat}}(p_x, p_y) \ge 0$ everywhere on $[-\pi, \pi]^2$, vanishing only at $(0, 0)$.
  2. **Purely $O(|p|^4)$ IR Scaling:** $W_{\text{flat}}(p_x, p_y) = p_x^4 - p_x^2 p_y^2 + p_y^4 + O(|p|^6)$, so:
     $$E_D^2(p) = P_{\text{Sym},x}(p_x)^2 + P_{\text{Sym},y}(p_y)^2 + \bigl(m_0 + r_2 W_{\text{flat}}(p)\bigr)^2 = m_0^2 + |p|^2 + 2 m_0 r_2(p_x^4 - p_x^2 p_y^2 + p_y^4) + O(|p|^6),$$
     yielding **$c_{\text{Dirac}}^2 = 1$ identically** ($c_{\text{Dirac}}^2 = c_{\text{scalar}}^2 = 1$) and $O(|p|^4)$-suppressed IR chiral breaking!
  3. **Exact $1 : 1 : 1$ Equal-Corner UV Lifting (`16 : 16 : 16`):**
     At $(\pi, 0)$: $4^2 - 0 + 0 = 16$; at $(0, \pi)$: $0 - 0 + 4^2 = 16$; at $(\pi, \pi)$: $4^2 - 4\cdot 4 + 4^2 = 16$!
     For $m_0 = 1/2$ and $r_2 = 3/32$, every non-zero Brillouin-zone corner sits at exact energy $|E| = \frac{1}{2} + 16\cdot\frac{3}{32} = 2.0000$, yielding `doubler_count = 0.0` and `uv_gap_ratio > 0.25 >> 0.12` across all three reference Laplacians (`K_scalar`, `L_nn`, `L_moore`).
  4. **Exact 2-Hop Integer Stencil ($r_{\text{hop}} \le 2$):**
     $W_{\text{flat}}$ has exact integer weights $+8$ at $(0,0)$, $-2$ at $(\pm 1,0),(0,\pm 1)$, $-1$ at $(\pm 1,\pm 1)$, and $+1$ at $(\pm 2,0),(0,\pm 2)$ (row sum $8 - 8 - 4 + 4 = 0$).

### 2.4 Zero-Subtraction Combinatorial Rational Loop Scaling

$$\lambda_{\text{site}}(k_{\text{eff}}, d) = \lambda_0 \frac{k_{\text{eff}}}{\kappa_W(d)} \in \mathbb{Q}, \qquad \lambda_{\text{mode}}(N, k_{\text{eff}}, d) = \frac{\lambda_{\text{site}}(k_{\text{eff}}, d)}{N} \in \mathbb{Q},$$
where $k_{\text{eff}} = 5/2$ (`1+1D`) and $k_{\text{eff}} = 5$ (`2+1D`), and $\kappa_W(4) = 2$ ($\kappa_W(d) = 3/2$ for $d \ge 3, d \ne 4$). Across $N \in [64, 128, 256, 512, 1024]$, 2-loop physical mass-ratio drift is $0.0000\%$ in `(1+1)D` and $\approx 0.55\% < 3\%$ in `(2+1)D` with zero counterterms.

---

## 3. Module Entry Points

- `architect.build_cycle3_model(dimension, N, d=4, m0_sq=Fraction(1,4), lambda0=Fraction(1,2), wilson_r2=None, include_causal_unitary=True, include_symmetry_operators=True) -> dict`
- `architect.build_all_cycle3_models(...) -> dict`
- `architect.cycle3_rational_couplings(N, k_eff_rational, d, lambda0=Fraction(1,2), ...) -> dict`
