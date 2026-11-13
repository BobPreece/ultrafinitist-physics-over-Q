# Cycle 3 Final Quantitative Simulation & Zero-Handicap Baseline Report

- **Total Wall-Clock Runtime:** `1622.43 s`
- **Total Q-Engine Trajectories Executed:** `38`
- **SIGALRM (> 5.0s) Stalls:** `0`
- **Peak Bit-Length Across All Runs:** `2954 <= 10,000`

## 1. Head-to-Head Comparison Table (`Q`-Engine vs Zero-Handicap `float64` Continuum Integrators)

| Battle & Regime (`dt`) | Integrator | Stability / Steps Completed | Invariant Drift (`E_kin`, `Z`, `P`, `div`, `mass`) | Density `min rho` / Bounds | Convergence Order | Reversibility |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Battle 1A** (`N=128`, `dt=1/4`, `n=300`) | **Exact `Q`-Engine** | `300/300` (`1938b` plateau) | **`0/1` exact** (`Gamma, E_kin, Z, Pars`) | Bounded `E_0` | **`O(dt^2)`** | **`0/1` exact** |
| **Battle 1A** (`N=128`, `dt=1/4`, `n=300`) | `float64` Forward Euler | **BLOWUP @ step 148** | `NaN / Inf` | Unbounded | `O(dt^1)` (`1.002`) | Irreversible |
| **Battle 1A** (`N=128`, `dt=1/4`, `n=300`) | `float64` Midpoint/RK2 | `300/300` | `dE=1.97e-01, dZ=6.68e-01` | Unbounded | `O(dt^2)` (`2.000`) | Irreversible |
| **Battle 1A** (`N=128`, `dt=1/4`, `n=300`) | `float64` Classical RK4 | `300/300` | `dE=3.84e-04, dZ=1.33e-03` | Unbounded | `O(dt^4)` (`4.000`) | Irreversible |
| **Battle 2A** (`L=6`, `dt=1/2`, `n=250`) | **Exact `Q`-Engine** | `250/250` (`168b` plateau) | **`0/1` exact** (`div, P, E_kin, Pars`) | `Z_max/Z0=1.722` | **`O(dt^2)`** | **`0/1` exact** |
| **Battle 2A** (`L=6`, `dt=1/2`, `n=250`) | `float64` Forward Euler | **BLOWUP @ step 7** | `NaN / Inf` | Blow-up | `O(dt^1)` (`1.071`) | Irreversible |
| **Battle 2A** (`L=6`, `dt=1/2`, `n=250`) | `float64` Midpoint/RK2 | **BLOWUP @ step 5** | `NaN / Inf` | Blow-up | `O(dt^2)` (`2.013`) | Irreversible |
| **Battle 2A** (`L=6`, `dt=1/2`, `n=250`) | `float64` Classical RK4 | **BLOWUP @ step 9** | `NaN / Inf` | Blow-up | `O(dt^4)` (`4.056`) | Irreversible |
| **Battle 2A** (`L=6`, `dt=1/512` ref) | `float64` `RK4(dt=1/512)` | Fine reference | `~1e-13` float roundoff | `Z_max/Z0=1.6059` | `O(dt^4)` ref | `~1e-13` |
| **Battle 3** (`L=6`, `dt=1/4`, `n=250`) | **Exact `Q`-Engine** | `250/250` (`188b` plateau) | **`0/1` exact** (`div, P, E_kin, mass`) | **`min_rho=0.0013 > 0`**, `||u||^2 <= B_u` | **`O(dt^2)`** (`u`) | **`0/1`** (vel) |
| **Battle 3** (`L=6`, `dt=1/4`, `n=250`) | `float64` Forward Euler (centered) | **BLOWUP @ step 6** | `NaN / Inf` | **`min_rho < 0` (Gibbs)** | `O(dt^1)` (`1.061`) | Irreversible |
| **Battle 3** (`L=6`, `dt=1/4`, `n=250`) | `float64` Midpoint/RK2 (centered) | **BLOWUP @ step 4** | `NaN / Inf` | **`min_rho < 0` (Gibbs)** | `O(dt^2)` (`2.001`) | Irreversible |
| **Battle 3** (`L=6`, `dt=1/4`, `n=250`) | `float64` Classical RK4 (centered) | **BLOWUP @ step 2** | `NaN / Inf` | **`min_rho < 0` (Gibbs)** | `O(dt^4)` (`3.991`) | Irreversible |
| **Battle 3** (`L=6`, `dt=1/512` ref) | `float64` `RK4(dt=1/512)` (centered) | Smooth `t in [0, 1]` | `~1e-12` | **`min_rho = -269.175 < 0`** | `O(dt^4)` ref | Irreversible |

## 2. Gate 4 Zero-Handicap Continuum `float64` Slope Verification (`RK4(dt=1/512)` Ref)

| Config | `max |rhs_f64 - rhs_Q| / ||rhs_Q||` | Forward Euler Order | Midpoint/RK2 Order | Classical RK4 Order | Passed |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `B1A_m6` | `1.53e-16` | `0.9995` | `2.0001` | `3.9915` | **PASS** |
| `B1A_m7` | `6.90e-16` | `1.0015` | `1.9999` | `4.0002` | **PASS** |
| `B1B_m6` | `1.53e-16` | `1.0005` | `2.0006` | `3.9956` | **PASS** |
| `B2A_L6` | `1.01e-15` | `1.0708` | `2.0127` | `4.0557` | **PASS** |
| `B2A_L8` | `1.87e-15` | `1.0348` | `1.9949` | `3.9975` | **PASS** |
| `B2B_L6` | `8.64e-16` | `1.0731` | `2.0123` | `4.0516` | **PASS** |
| `B3_L6_centered` | `1.04e-15` | `1.0608` | `2.0009` | `3.9907` | **PASS** |
| `B3_L6_varga_gen` | `1.04e-15` | `1.0611` | `2.0013` | `3.9950` | **PASS** |
| `B3_L8_centered` | `9.20e-16` | `1.0395` | `2.0004` | `3.9975` | **PASS** |
| `B3_L8_varga_gen` | `3.67e-16` | `1.0393` | `2.0035` | `4.0028` | **PASS** |

## 3. Cycle 3 `Q`-Engine Trajectory Matrix (`0/1` Invariants, `O(1)` Bit Plateaus, & `< 5s` Speed)

| Run Tag | Steps | Peak Bits (`state` / `phys`) | Late Bit Slope | Max Step (s) | Exact `0/1` Invariants | Extra Gate Metrics |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| `B1A_m6_dt1_2` | `300/300` | `1244 (66 / 1244)` | `-0.0144` | `0.021s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=0.00, IR: 0.179->0.652, rev10=True |
| `B1A_m6_dt1_4` | `300/300` | `1254 (75 / 1254)` | `+0.0241` | `0.020s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=0.00, IR: 0.179->0.360, rev10=True |
| `B1A_m6_dt1_8` | `300/300` | `1268 (70 / 1268)` | `+0.0447` | `0.020s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=0.00, IR: 0.179->0.466, rev10=True |
| `B1A_m6_dt1_4_prime` | `300/300` | `2199 (74 / 2199)` | `+0.0829` | `0.026s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=0.00, IR: 0.179->0.360, rev10=True |
| `B1A_m7_dt1_2` | `300/300` | `1952 (68 / 1952)` | `+0.0596` | `0.046s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=19.18, IR: 0.655->0.810, rev10=True |
| `B1A_m7_dt1_4` | `300/300` | `1938 (72 / 1938)` | `+0.0093` | `0.052s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=19.18, IR: 0.655->0.809, rev10=True |
| `B1A_m7_dt1_8` | `300/300` | `1933 (77 / 1933)` | `+0.0181` | `0.057s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=19.18, IR: 0.655->0.750, rev10=True |
| `B1A_m7_dt1_4_prime` | `300/300` | `2930 (72 / 2930)` | `-0.0115` | `0.069s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=19.18, IR: 0.655->0.809, rev10=True |
| `B1A_m7_dt1_2_prime` | `300/300` | `2954 (68 / 2954)` | `+0.0283` | `0.064s` | `0/1 (Gamma,E,Z,Pars)` | R1_0=19.18, IR: 0.655->0.810, rev10=True |
| `B1B_m6_dt1_2` | `300/300` | `284 (150 / 284)` | `-0.0172` | `0.233s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=1.032e-03 |
| `B1B_m6_dt1_4` | `300/300` | `287 (152 / 287)` | `+0.0265` | `0.234s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=1.120e-03 |
| `B1B_m6_dt1_8` | `300/300` | `291 (154 / 291)` | `+0.0006` | `0.190s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=3.219e-03 |
| `B1B_m6_dt1_4_prime` | `300/300` | `427 (206 / 427)` | `+0.0070` | `0.279s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=1.120e-03 |
| `B1B_m7_dt1_2` | `300/300` | `417 (187 / 417)` | `+0.0274` | `0.919s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=2.346e-03 |
| `B1B_m7_dt1_4` | `300/300` | `413 (188 / 413)` | `+0.0259` | `0.937s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=3.920e-03 |
| `B1B_m7_dt1_8` | `300/300` | `415 (190 / 415)` | `+0.0058` | `0.735s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=1.319e-02 |
| `B1B_m7_dt1_4_prime` | `300/300` | `548 (248 / 548)` | `+0.0049` | `1.036s` | `0/1 (Gamma,E+Q), dQ>=0` | min dQ=3.920e-03 |
| `B2A_L6_dt1_2` | `250/250` | `168 (67 / 168)` | `+0.0274` | `0.043s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=1.7215, UV: 0.000->0.999, rev=True |
| `B2A_L6_dt1_4` | `250/250` | `166 (69 / 166)` | `+0.0707` | `0.048s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=1.6120, UV: 0.000->0.996, rev=True |
| `B2A_L6_dt1_8` | `250/250` | `167 (67 / 167)` | `+0.0523` | `0.042s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=1.6906, UV: 0.000->1.000, rev=True |
| `B2A_L6_dt1_4_prime` | `250/250` | `315 (67 / 315)` | `+0.0403` | `0.052s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=1.6120, UV: 0.000->0.996, rev=True |
| `B2A_L6_dt1_4_galilean` | `250/250` | `169 (67 / 169)` | `+0.0170` | `0.054s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=1.4175, UV: 0.000->0.877, rev=True |
| `B2B_L6_dt1_2` | `250/250` | `208 (119 / 208)` | `+0.0001` | `0.188s` | `0/1 (div,P,E+Q), dQ>=0` | Z_max/Z0=1.5859, min dQ=6.943e-06 |
| `B2B_L6_dt1_4` | `250/250` | `211 (121 / 211)` | `-0.0222` | `0.214s` | `0/1 (div,P,E+Q), dQ>=0` | Z_max/Z0=1.2473, min dQ=7.300e-05 |
| `B2B_L6_dt1_8` | `250/250` | `210 (123 / 210)` | `+0.0312` | `0.212s` | `0/1 (div,P,E+Q), dQ>=0` | Z_max/Z0=1.1825, min dQ=4.972e-02 |
| `B2B_L6_dt1_4_prime` | `250/250` | `365 (167 / 365)` | `-0.0252` | `0.582s` | `0/1 (div,P,E+Q), dQ>=0` | Z_max/Z0=1.2473, min dQ=2.691e-04 |
| `B2A_L8_dt1_2` | `250/250` | `205 (67 / 205)` | `+0.0456` | `0.155s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=2.1804, UV: 0.000->0.700, rev=True |
| `B2A_L8_dt1_4` | `250/250` | `213 (69 / 213)` | `+0.0021` | `0.153s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=2.1580, UV: 0.000->0.314, rev=True |
| `B2A_L8_dt1_8` | `250/250` | `209 (67 / 209)` | `+0.0782` | `0.168s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=2.7400, UV: 0.000->0.542, rev=True |
| `B2A_L8_dt1_4_prime` | `250/250` | `469 (92 / 469)` | `+0.0680` | `0.227s` | `0/1 (div,P,E,Pars)` | Z_max/Z0=2.1580, UV: 0.000->0.314, rev=True |
| `B2B_L8_dt1_4` | `250/250` | `256 (138 / 256)` | `+0.0497` | `0.826s` | `0/1 (div,P,E+Q), dQ>=0` | Z_max/Z0=1.5917, min dQ=5.413e-02 |
| `B3_L6_dt1_2` | `250/250` | `186 (67 / 186)` | `-0.0063` | `0.246s` | `0/1 (div,P,E,mass)` | min_rho=0.00400>0, rho_max(t=1)=9.784, ||u||^2/B_u=0.140, ||om||^2/B_om=0.073, rev4(vel)=True |
| `B3_L6_dt1_4` | `250/250` | `188 (68 / 188)` | `+0.0240` | `0.240s` | `0/1 (div,P,E,mass)` | min_rho=0.00126>0, rho_max(t=1)=6.516, ||u||^2/B_u=0.088, ||om||^2/B_om=0.050, rev4(vel)=True |
| `B3_L6_dt1_4_prime` | `250/250` | `316 (77 / 316)` | `-0.0084` | `0.277s` | `0/1 (div,P,E,mass)` | min_rho=0.00126>0, rho_max(t=1)=6.516, ||u||^2/B_u=0.088, ||om||^2/B_om=0.050, rev4(vel)=True |
| `B3_L8_dt1_2` | `250/250` | `300 (68 / 300)` | `+0.0205` | `0.644s` | `0/1 (div,P,E,mass)` | min_rho=0.00201>0, rho_max(t=1)=9.896, ||u||^2/B_u=0.036, ||om||^2/B_om=0.029, rev4(vel)=True |
| `B3_L8_dt1_4` | `250/250` | `299 (67 / 299)` | `+0.0418` | `0.671s` | `0/1 (div,P,E,mass)` | min_rho=0.00117>0, rho_max(t=1)=11.195, ||u||^2/B_u=0.041, ||om||^2/B_om=0.050, rev4(vel)=True |
| `B3_L8_dt1_4_prime` | `250/250` | `475 (95 / 475)` | `-0.0066` | `0.562s` | `0/1 (div,P,E,mass)` | min_rho=0.00117>0, rho_max(t=1)=11.195, ||u||^2/B_u=0.041, ||om||^2/B_om=0.050, rev4(vel)=True |
| `B3_L8_dt1_4_boost` | `250/250` | `300 (68 / 300)` | `+0.0413` | `0.692s` | `0/1 (div,P,E,mass)` | min_rho=0.00139>0, rho_max(t=1)=11.195, ||u||^2/B_u=0.040, ||om||^2/B_om=0.025, rev4(vel)=True |

## 4. Exact `Q`-Engine `O(dt^2)` Convergence Against `RK4(dt=1/512)` Reference

| Config | Observed `L2` Orders (`dt` halving) | Smallest `dt` `L2` RMS Error | Extra Notes |
| :--- | :---: | :---: | :--- |
| `B2A_L6_T1` | `[2.6671, 2.0817, 1.9681]` | `6.773e-03` | T1, Z/Z0_Q=1.2189 vs ref 1.2189 |
| `B2A_L6_T2` | `[2.1461, 1.7948, 1.9455]` | `4.956e-02` | T2, Z/Z0_Q=1.5319 vs ref 1.5368 |
| `B2A_L8_T1` | `[2.6156, 2.0929]` | `2.284e-02` | T1, Z/Z0_Q=1.3423 vs ref 1.3446 |
| `B1A_m6` | `[1.9448, 1.986, 1.9965]` | `5.136e-05` | T=1.0 |
| `B1A_m7` | `[2.9565, 2.1911, 2.0133]` | `1.545e-04` | T=1.0 |
| `B1B_m6` | `[2.0152, 2.0073, 2.0036]` | `8.774e-07` | T=1.0 |
| `B3_L6` | `[0.0945, 0.4563, 1.6189]` | `1.238e+00` | T=0.5 |

## 5. Cycle 2 Defect Remediation Checks (`C2-2` Shift Equivariance & `C2-3` Off-Manifold Tangents)

- **`C2-2` (`L=8` Shift-Equivariant Spectrum):** `max_abs_shift_diff = 0` (`exact_zero = True`)
- **`C2-3` Off-Manifold `B1A_m6_t2` (`t=2.0`):** `carrier_split_linf = 3.627e-04`, `rel_tangent_err(2^-30) = 4.070e-07`
- **`C2-3` Off-Manifold `B2A_L6_t2` (`t=2.0`):** `carrier_split_linf = 6.186e-01`, `rel_tangent_err(2^-30) = 6.812e-04`
- **`C2-3` Off-Manifold `B3_L6_t2` (`t=2.0`):** `carrier_split_linf = 8.544e+00`, `rel_tangent_err(2^-30) = 1.299e-02`
