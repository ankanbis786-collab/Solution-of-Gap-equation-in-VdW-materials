# Figure 3 Reproduction Code  
## Interlayer Charge-Transfer Ferroelectric Fluctuations as a Pairing Mechanism in vdW Superconductors

This repository contains the numerical code used to generate **Figure 3** of the paper:

> *Interlayer Charge-Transfer Ferroelectric Fluctuations as a Pairing Mechanism in van der Waals Superconductors*  
> Ankan Biswas et al. (2026)

The paper develops a theory of superconductivity mediated by interlayer ferroelectric (FE) fluctuations in van der Waals materials.

---

# check_splitting.m

## Description

`check_splitting.m` numerically solves the **linearized Eliashberg equation in the ordered ferroelectric (FE) phase**, implementing **Eq. (B23)** of the paper.

The script computes the layer-resolved superconducting transition temperatures  
\( T_{c,t} \) and \( T_{c,b} \), and extracts the FE-induced splitting

\[
\Delta T_c = T_{c,t} - T_{c,b}.
\]

It reproduces the ordered-phase results shown in Figure 3 and verifies the quadratic scaling

\[
\Delta T_c / T_c \sim r^2
\]

near the FE quantum critical point.

---

## Theory Implemented

The script implements the ordered-state linearized gap equation:

\[
\Delta_l(k_0) =
\frac{\lambda \nu_F}{k_{F,l} a}
T \sum_{p_0 \neq k_0}
\frac{
[d(p_0-k_0) + d(p_0+k_0)]
}{
p_0
}
\frac{\Delta_l(p_0)}
{
1 + \frac{\lambda \nu_F}{k_{F,l} a}
\frac{1}{k_0}
\sum_{p_0' \neq k_0}
[d(p_0'-k_0) - d(p_0'+k_0)]
}
\]

where:

- \( l = t, b \) labels the layers
- \( d(q_0) \) is the momentum-integrated bosonic propagator
- The ordered phase enters via:
  - \( \mu_l = \mu \pm \delta\mu \)
  - \( r \rightarrow 2|r| \) in the bosonic propagator
  - Layer-dependent \( k_{F,l}, v_{F,l} \)

The zero-frequency divergence is removed using the standard transformation:

\[
\Delta_l(k_0) = \frac{k_0 \Phi_l(k_0)}{\tilde{\Sigma}_l(k_0)}.
\]

---

## Numerical Procedure

1. **Discretize Matsubara frequencies**
   \[
   k_0 = (2n+1)\pi T
   \]

2. **Construct bosonic kernel**
   \[
   d(q_0) = \int_0^\Lambda
   \frac{dq}{
   2|r| + q^2 a^2 +
   \frac{4\lambda \nu_F |q_0|}{v_{\mathrm{eff}} q}
   }
   \]

3. **Compute layer-dependent self-energy**
   \[
   \Sigma_l(k_0)
   \]

4. **Assemble pairing matrix** from Eq. (B23)

5. **Solve eigenvalue problem**
   - Largest eigenvalue \( \lambda_{\max}(T) \)
   - Identify \( T_c \) from condition:
     \[
     \lambda_{\max}(T_c) = 1
     \]

6. Repeat separately for \( l = t, b \)

---

## Physical Content

- FE order generates a layer chemical potential shift
  \[
  \delta\mu = \lambda u_0
  \]
- This induces small differences in:
  - \( k_{F,l} \)
  - \( v_{F,l} \)
  - pairing kernel prefactor
- The resulting Tc splitting is parametrically small:
  \[
  \Delta T_c \propto r^2
  \]
  consistent with Appendix C.

The script numerically confirms this scaling.

---

## Outputs

- \( T_{c,t}(r) \)
- \( T_{c,b}(r) \)
- \( \Delta T_c(r) \)
- Ordered vs disordered phase comparison
- Data used for Figure 3

---

## Parameters

Typical inputs include:

- \( k_F \)
- \( m^* \)
- \( \lambda \)
- \( r \) (distance from QCP)
- UV frequency cutoff \( \sim E_F \)

Convergence with respect to Matsubara and momentum cutoffs must be checked.

---

## Usage

Run in MATLAB:

```matlab
check_splitting

