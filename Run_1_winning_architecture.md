# Winning Ultrafinitist QFT Architecture — Exhibit I (`ultrafinitist-qft` v2.0)

**Primary Implementation:** `architect/v2_pure_q_engine.py` & Part I of `standalone_ultrafinitist_qft.py`  
**Symbolic Verifier:** `formalist/cas.py` & `formalist/verify_proofs.py` (`/usr/bin/python3 -I`, `EXIT 0`)  
**Constitutional Standard:** Strict Part I Pure-`Q` Parity (`0` Gate 1 AST violations, `quantize_branch_sites == []`, `dead_ornament_slots == []`, Fair-Fight $1/N$ Loop Normalization)

---

## 1. Strict Part I Pure-`Q` Constitutional Foundation

Above `# === PART II: CONTINUUM BASELINES, BENCHMARK SUITE & PUBLICATION PLOTTER ===` in `standalone_ultrafinitist_qft.py` (and throughout `architect/v2_pure_q_engine.py` and `formalist/`), every operator, state vector, propagator solve, quadrature bound, time-evolution step, and block-spin coarse-graining projection executes **exclusively over $\mathbb{Q}$** (`fractions.Fraction` and Python `int`):

1. **Allowed Part I Imports Only:** `from __future__ import annotations`, `from fractions import Fraction`, and `typing`. Zero `float` or `complex` literals/constructors, zero `math`/`cmath`/`ctypes`/`Accelerate.framework`/`numpy`/`scipy`/`sympy`, zero `sqrt`/`isqrt`, zero fractional/negative or variable `**`, zero `round()`/`.limit_denominator()`, zero `while` loops, zero `try/except`, zero Newton–Heron `x + S/x`, zero `isinstance`/`type`/`hasattr`, and zero `max`/`min`/`abs` clamps.
2. **Exact Rational Spinor Field Algebra over $\mathbb{Q}^2$:** Complex quantum/Dirac amplitudes live in $\mathbb{Q}^2$ as pairs $(a, b) \equiv a + i b$ with exact bilinear multiplication $(a, b) \cdot (c, d) = (ac - bd, ad + bc)$ and Hermitian inner product $\langle \Psi, \Phi \rangle_{\text{Re}} = \sum_k (\psi_{\text{re},k}\phi_{\text{re},k} + \psi_{\text{im},k}\phi_{\text{im},k})$.
3. **Zero-Branch `QuantizerBundle` Protocol:**
   - `RuntimeDyadicQuantizerBundle(bits=48, max_exp=64)` maps nonzero rationals into the dyadic subring $\mathbb{Z}[1/2]$ ($m / 2^e$ with $m < 2^{48}, |e| \le 64$) while `dt` is factored **strictly outside** `bundle(...)` at every integrator stage (`up1 = up0 + qdt * bundle(gam * k1)`). Consequently, discrete tangents $(\Psi(t + \Delta t) - \Psi(t))/\Delta t$ remain unfrozen down to $\Delta t = 2^{-30}$ at both $t = 0$ and $t > 0$ ($n = 5$).
   - `clear_denominators(xs)` clears rational denominators into fast Python `int` stencil arithmetic using loop-free bitwise dyadic valuation `lowbit = d & (-d)` and odd-factor pooling.
   - `IdentityQuantizerBundle` executes the **identical unbranched AST** as an exact algebraic pass-through over `Fraction`, `QQPoly`, `QQRat`, `QQMVPoly`, `QQMVRat`, and `QQJet`.
4. **Pole-Free `[4/3]` Binomial Padé Square Root (`pade_sqrt_pos`):** Evaluates $\sqrt{w}$ around $w = 1$ via the strictly positive rational Padé map
   $$R_{4,3}(w) = \frac{1 + 28w + 70w^2 + 28w^3 + w^4}{8 + 56w + 56w^2 + 8w^3}, \qquad R_{4,3}(w)^2 - w = \frac{(w - 1)^8}{64(1 + 7w + 7w^2 + w^3)^2}$$
   composed in two fixed multiplicative stages ($O((w-1)^{16})$ contact at $w = 1$, strictly positive denominator on $[0, \infty)$, zero Newton–Heron `x + S/x`).

---

## 2. Battle 1: Gaussian-Integer Quotient Torus $V_m$, 4th-Order Symanzik Stencils & $W_{\text{flat}}$ Doubler Elimination

### 2.1 Finite Weyl–Schwinger Clock-and-Shift Field Algebra ($d = 4$ over $\mathbb{Q}^2$)
At each site $x \in V$ ($|V| = N$, total Hilbert dimension $D = d^N < \infty$), canonical scalar conjugates are finite unitary $d \times d$ matrices $U, V$. At $d = 4$, $\omega_4 = i = (0, 1) \in \mathbb{Q}^2$, so every matrix entry lies strictly in $\mathbb{Q}^2$:
$$U = \operatorname{diag}(1, i, -1, -i), \qquad V_{r,c} = \delta_{(r-c) \bmod 4, 1}, \qquad \Phi = \frac{U - U^\dagger}{2i} = \operatorname{diag}(0, 1, 0, -1), \qquad \Pi = \frac{V - V^\dagger}{2i}.$$
Machine-verified to `0/1`: $U^4 = V^4 = I_4$, $U V = i V U$, $\frac{1}{4}\operatorname{Tr}(\Phi^2) = \frac{1}{2}$, $\frac{1}{4}\operatorname{Tr}(\Phi^4) = \frac{1}{2}$, and Weyl kurtosis $\kappa_W(4) = \frac{\operatorname{Tr}(\Phi^4)/4}{(\operatorname{Tr}(\Phi^2)/4)^2} = 2$.

### 2.2 Deterministic $(1+1)\text{D}$ Ring $\mathbb{Z}_N$ and $(2+1)\text{D}$ Gaussian Quotient Torus $V_m = \mathbb{Z}[i]/\langle(1+i)^m\rangle$
- **$(1+1)\text{D}$:** Cyclic ring $\mathbb{Z}_N$ across $N \in \{64, 128, 256, 512, 1024\}$ with exact spatial inversion $x \mapsto -x \bmod N$.
- **$(2+1)\text{D}$:** Quotient ring $V_m = \mathbb{Z}[i]/\langle(1+i)^m\rangle$ of exact norm $N = |(1+i)^m|^2 = 2^m \in \{64, 128, 256, 512, 1024\}$ ($m \in \{6, 7, 8, 9, 10\}$):
  - Even $m = 2k$ ($L = 2^k$): $(1+i)^{2k} = (2i)^k$ generates the square ideal $L\mathbb{Z}[i]$.
  - Odd $m = 2k + 1$ ($L = 2^k, N = 2L^2 \in \{128, 512\}$): $(1+i)^{2k+1}$ generates the twisted $45^\circ$ index-2 sublattice $\{(X, Y) \in (L\mathbb{Z})^2 : X \equiv Y \pmod{2L}\}$.
  - Because $i \in \mathbb{Z}[i]^\times$ is a ring unit, multiplication $z \mapsto i z$ ($(x, y) \mapsto (-y, x) \bmod \langle(1+i)^m\rangle$) is an **exact order-4 $90^\circ$ $C_4$ ring automorphism** at **every** integer $m$ (eliminating the $2:1$ rectangular grid anisotropy that afflicts naive $\mathbb{Z}_{L_x} \times \mathbb{Z}_{L_y}$ lattices at $N = 128, 512$).

### 2.3 4th-Order Symanzik Operators & Equal-Corner SOS Biharmonic Wilson Mass $W_{\text{flat}}$
Let $(L_\mu f)(z) = 2f(z) - f(z + e_\mu) - f(z - e_\mu)$ ($L_\mu(p) = 2 - 2\cos p_\mu$) and $(D_{s,\mu} f)(z) = \frac{f(z + s e_\mu) - f(z - s e_\mu)}{2s}$:
- **4th-Order Symanzik Scalar Kinetic Operator:**
  $$K_{\text{Sym}} = \sum_\mu \left(L_\mu + \frac{1}{12} L_\mu^2\right) = \sum_\mu \left(\frac{4}{3} L_\mu - \frac{1}{12} L_{2,\mu}\right), \qquad K_{\text{Sym}}(p) = |p|^2 - \frac{1}{90}\sum_\mu p_\mu^6 + O(|p|^8).$$
- **4th-Order Symanzik Antisymmetric Spinor Derivative:**
  $$D_{\text{Sym},\mu} = \frac{4}{3} D_{1,\mu} - \frac{1}{3} D_{2,\mu}, \qquad P_{\text{Sym},\mu}(p) = \frac{4}{3}\sin p_\mu - \frac{1}{6}\sin(2p_\mu) = p_\mu - \frac{1}{30}p_\mu^5 + O(p_\mu^7).$$
- **Equal-Corner Sum-of-Squares Biharmonic Wilson Operator (`16 : 16 : 16`):**
  - In $(1+1)\text{D}$: $W_{\text{1D}} = L_1^2 \succeq 0$, $M_W^{(1+1)} = m_0 I + r_2 L_1^2$ with $m_0 = \frac{1}{2}, r_2 = \frac{1}{16}$ ($M_W(\pi) = \frac{1}{2} + \frac{16}{16} = \frac{3}{2}$).
  - In $(2+1)\text{D}$:
    $$W_{\text{flat}} = L_x^2 - L_x L_y + L_y^2 = \left(L_x - \frac{1}{2}L_y\right)^2 + \frac{3}{4}L_y^2 \succeq 0, \qquad M_W^{(2+1)} = m_0 I + r_2 W_{\text{flat}}, \quad r_2 = \frac{3}{32}.$$
  - **Exact Corner & IR Properties (Certified in `QQJet(order=6)`):**
    1. At the 4 Brillouin-zone corners $(0,0), (\pi,0), (0,\pi), (\pi,\pi)$, $W_{\text{flat}} = (0, 16, 16, 16)$ (`16:16:16` equal lift!), raising all 3 nonzero corners to $E_{\text{Dirac}}^2 = (\frac{1}{2} + \frac{3}{32}\cdot 16)^2 = 4$ (`0` Nielsen–Ninomiya doublers, vs `3` spurious doublers on the naive cubic lattice).
    2. Because $W_{\text{flat}}(p) = p_x^4 - p_x^2 p_y^2 + p_y^4 + O(|p|^6)$ vanishes at $O(|p|^4)$, $M_W(p)^2 = m_0^2 + 2 m_0 r_2 W_{\text{flat}}(p) + O(|p|^8)$ contributes **zero $O(|p|^2)$ cross-term**, locking $c_{\text{Dirac}}^2 = c_{\text{scalar}}^2 = \frac{1}{1}$ to `0/1` (fixing Cycle 2's $c_{\text{Dirac}}^2 = 1 + m_0 = \frac{3}{2}$ lightcone mismatch), delivering a **`21.4x` Dirac dispersion isotropy gain** at $N = 1024$ (`0.3336 -> 0.0156`), and suppressing IR chiral breaking cubically by **$8 = 2^3$ per doubling** of linear scale.

---

## 3. Battle 2: Varga Positivity, Rational Quadrature Brackets & Remediated Fair-Fight $1/N$ UV Loops

### 3.1 Varga Squared-$M$-Matrix Regular Splitting ($G_{xy} > 0/1$ Entrywise)
Although $K_{\text{Sym}}$ has positive 2-hop stencil weights $+\frac{1}{12}$ (so $K_{\text{Sym}} + m_0^2 I$ is not a naive Stieltjes $M$-matrix), it admits the exact **Varga Regular Splitting** $K_{\text{Sym}} + m_0^2 I = P - R$ with factored squared-$M$-matrix $P = \frac{1}{12}(8 I - A_{\text{nn}})^2$:
$$R^{(1+1)} = (3 - m_0^2) I_N \ge 0 \quad (0 < m_0^2 \le 3), \qquad R^{(2+1)} = \left(\frac{2}{3} - m_0^2\right) I_N + \frac{1}{6} A_x A_y \ge 0 \quad \left(0 < m_0^2 \le \frac{2}{3}\right).$$
Since $B = 8I - A_{\text{nn}}$ is a strictly diagonally dominant irreducible $M$-matrix ($B^{-1} > 0 \implies P^{-1} = 12 B^{-2} > 0$) and $R \ge 0$, the Neumann series $G = (I - P^{-1}R)^{-1}P^{-1} = \sum_{k=0}^\infty (P^{-1}R)^k P^{-1} > 0$ proves **$G_{xy} > 0/1$ strictly entrywise** for all $0 < m_0^2 \le \frac{2}{3}$ (including $m_0^2 = \frac{1}{4}$), avoiding the oscillatory negative spatial ringing ($G_{\min} = -0.0356 < 0$) of unfactored biharmonic higher-derivative regulators.

### 3.2 Fair-Fight $1/N$ Loop Normalization & Continuum $a \to 0$ Contrast
In Exhibit I v2.0, **both** Naive Cubic and Cycle 3 are evaluated with the **exact same $1/N$ per-site normalization** (zero `* N` multiplier on the cubic baseline) and identical rational combinatoric coupling rule $\lambda_{\text{site}} = \lambda_0 k_{\text{eff}} / \kappa_W(d)$ ($k_{\text{eff}} = s_1 = \frac{1}{N}\operatorname{Tr}(K)$):
$$\tau_1(N) = \frac{1}{N}\operatorname{Tr}(G) = G_{00}, \qquad s_2(N) = \frac{1}{N}\sum_{x,y \in V} G_{xy}^3 = \sum_{x \in V} G_{0x}^3,$$
$$\Sigma^{(1)} = \frac{1}{2}\lambda_{\text{site}} \tau_1(N), \qquad \Sigma^{(2)} = -\frac{1}{6}\lambda_{\text{site}}^2 s_2(N), \qquad R_m^{(L)}(N) = \frac{m_0^2 + \sum_{\ell=1}^L \Sigma^{(\ell)}}{m_0^2}.$$
- **Global `0/1`-Slack Ceilings:** $\tau_1(N) \le \frac{1}{m_0^2} = 4$ and $s_2(N) \le \frac{1}{m_0^6} = 64$ for all finite $N \ge 1$ (`0/1` slack at zero hopping $K = 0$).
- **54th-Order Pure-`Q` Midpoint + Gauss–Lobatto Brackets:** Certifies $\tau_1(N)$ over $\mathbb{Q}$ to **`99.82%` tightness in $(1+1)\text{D}$** and **`98.53%` tightness in $(2+1)\text{D}$**.
- **Where UV Divergences Truly Come From:** At fixed ultrafinitist fundamental scale $a = 1, m_0 = \frac{1}{2}$, per-site loop traces remain strictly finite for *both* finite lattices (with $V_m$ removing odd-$m$ rectangular aspect-ratio kinks, doublers, and $O(p^2)$ anisotropy). Ultraviolet infinities appear **if and only if** one shrinks the lattice spacing $a = L_{\text{box}}/N^{1/d} \to 0$ ($\Lambda = \pi/a \to \infty$), causing continuum momentum integrals to diverge logarithmically in $(1+1)\text{D}$ and linearly in $(2+1)\text{D}$.

---

## 4. Battle 3: Unitary Double-Cayley Wavepacket Dynamics, Vieta Readout & `(1+i)`-Adic Block-Spin RG

### 4.1 Pole-Free Strang-Split 3-Shear Double-Cayley Step (`step_wavepacket`)
Decomposing the Hermitian Cycle 3 Dirac Hamiltonian $H_D = A + i B$ into real symmetric $A = A^\top$ ($\sigma_z \otimes M_W - \sigma_x \otimes i D_{\text{Sym},y}$) and real skew-symmetric $B = -B^\top = \sigma_x \otimes (-S_x)$ ($S_x = D_{\text{Sym},x}$), Schrödinger–Dirac evolution $\dot{\Psi} = -i H_D \Psi$ on $\Psi = u + i v \in \mathbb{Q}^{2N} \oplus i\mathbb{Q}^{2N}$ splits into $\dot{\Psi} = (B - i A)\Psi$. With unconditionally bounded rational damping factor $\gamma(\Delta t) = \frac{1}{1 + \frac{1}{16}\omega_{\text{cap}}^4 \Delta t^4} = 1 + O(\Delta t^4)$, `step_wavepacket` composes:
1. Half-step $\Delta t / 2$ spinor-bipartite 3-shear leapfrog on $B$ (`up -> dn -> up`),
2. Full-step $\Delta t$ real/imag-bipartite 3-shear leapfrog on $-i A$ (`u -> v -> u`),
3. Endpoint half-step $\Delta t / 2$ spinor-bipartite 3-shear leapfrog on $B$.

Because every shear is triangular and symmetric under $\Delta t \mapsto -\Delta t$, `step_wavepacket` requires **zero dense matrix inversions**, strictly respects a **finite stencil lightcone** (`0/1` leak outside $|x| \le 12$ per step), satisfies **exact `0/1` palindrome time-reversal** $F_{-\Delta t}^{10}(F_{+\Delta t}^{10}(\Psi_0)) = \Psi_0$, matches $e^{-i\Delta t H_D}\Psi_0$ to `0/1` through $O(\Delta t^2)$ in `QQJet(order=3)`, and obeys **Theorem 6 ($p$-adic odd-denominator invariance)**: since `RuntimeDyadicQuantizerBundle` emits increments in $\mathbb{Z}[1/2]$ and $\Delta t \in \{1/4, 1/8\} \subset \mathbb{Z}[1/2]$, $\Psi_n - \Psi_0 \in \mathbb{Z}[1/2]^{4N}$ for all $n \ge 0$, locking `odd_part(den) in {1, p}` (`76 < 150` peak bits, **`0.0000 bits/step`** asymptotic drift across all 21 primes `101..199`).

### 4.2 Birational Vieta Norm/Energy Manifold Readout (`vieta_readout`)
Given initial invariants $N_0 = \|\Psi_0\|_2^2$ and $E_0 = \langle \Psi_0, H_D \Psi_0 \rangle$ ($\widehat{H} = H_D - \frac{E_0}{N_0} I$, $\Xi_0 = \widehat{H}\Psi_0$), `vieta_readout` constructs a fixed reference state $\Phi_{\text{ref}}$ on the intersection manifold $\mathcal{M}_{N_0, E_0} = \{\|\Psi\|_2^2 = N_0\} \cap \{\langle \Psi, \widehat{H}\Psi \rangle = 0\}$, rotates its $U(1)$ phase $\Phi_{\text{mid}} = C\Phi_{\text{ref}} + S(i\Phi_{\text{ref}})$ via a singularity-free Pythagorean map centered at $b_0 = \operatorname{Re}\langle i\Phi_{\text{ref}}, \Xi_0\rangle < 0$ so $\operatorname{Re}\langle \Phi_{\text{mid}}, \widehat{H}\Psi_1\rangle \approx 0$, and performs an exact linear Vieta cone intersection + stereographic sphere lift:
$$\|\Psi_{\text{phys}}\|_2^2 - N_0 \equiv \frac{0}{1}, \qquad \langle \Psi_{\text{phys}}, H_D \Psi_{\text{phys}} \rangle - E_0 \equiv \frac{0}{1},$$
while preserving exact `0/1` $O(\Delta t^2)$ Taylor agreement in `QQJet(order=3)` and coupling all spinor slots (`dead_ornament_slots == []`).

### 4.3 Exact Rational `(1+i)`-Adic Block-Spin Decimation ($R_{\text{block}} : V_{m+2} \to V_m$)
Since $(1+i)^2 = 2i$ scales linear lattice coordinates by $2$ with $0^\circ$ rotation modulo units, the Kadanoff–bilinear restriction operator $R_{\text{block}} : V_{m+2} \to V_m$ ($1024 \to 256$, $512 \to 128$, $256 \to 64$) over $\mathbb{Q}$ with **zero counterterms** is:
$$(R_{\text{block}} f)(Z) = f(2Z) + \frac{1}{2}\sum_{s \in \{\pm 1, \pm i\}} f(2Z + s) + \frac{1}{4}\sum_{s \in \{\pm 1 \pm i\}} f(2Z + s).$$
Every fine site belongs to a unique $(1+i)^2$ coset with column weight sum $1/1$ (`0/1` partition of unity $\implies$ exact charge conservation $\sum_Z (R_{\text{block}} f)(Z) = \sum_z f(z)$), $R_{\text{block}} C_4^{(m+2)} = C_4^{(m)} R_{\text{block}}$ holds to `0/1` at both even and odd $m$, Householder projection through the $C_4$-invariant origin pole $e_{Z=0,\uparrow}$ preserves unit quantum norm $\|\Psi_{\text{coarse}}\|_2^2 = \|\Psi_{\text{fine}}\|_2^2 = 1/1$ to `0/1`, and the coarse kinetic operator $K_{\text{RG}}^{(m)} = R_{\text{block}} K_{\text{Sym}}^{(m+2)} R_{\text{block}}^\top$ satisfies in `QQJet(order=4)`:
$$K_{\text{RG}}^{(m)}(P_x, P_y) = \frac{1}{1}(P_x^2 + P_y^2) - \frac{1}{24}(P_x^4 + P_y^4) - \frac{1}{4}P_x^2 P_y^2 + O(|P|^6), \qquad c_{\text{RG}}^2 = \frac{1}{1}.$$
