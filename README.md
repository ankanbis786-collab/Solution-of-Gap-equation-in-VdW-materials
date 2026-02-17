# Figure 3 Reproduction Code  
## Interlayer Charge-Transfer Ferroelectric Fluctuations as a Pairing Mechanism in vdW Superconductors

This repository contains the numerical code used to generate **Figure 3** of the paper:

> *Interlayer Charge-Transfer Ferroelectric Fluctuations as a Pairing Mechanism in van der Waals Superconductors*  
> Ankan Biswas et al. (2026)

The paper develops a theory of superconductivity mediated by interlayer ferroelectric (FE) fluctuations in van der Waals materials.

---

# check_splitting.m

## Description

``check_splitting.m" numerically solves the linearized Eliashberg equation in the ordered ferroelectric (FE) phase ( Eq. (B23) in Appendix B)

The script computes the **layer-resolved superconducting transition temperatures** $T_{c,t}$ and $T_{c,b}$, and extracts the FE-induced splitting

$$
\Delta T_c = T_{c,t} - T_{c,b}.
$$

It reproduces the ordered-phase results shown in **Figure 3** and verifies the quadratic scaling

$$
\frac{\Delta T_c}{T_c} \sim r^2
$$

near the FE quantum critical point.

---

## Theory Implemented

The script implements the ordered-state linearized gap equation (Eq. B23):

$$
\Delta_l(k_0) =
\frac{\lambda \nu_F}{k_{F,l} a}
T \sum_{p_0 \neq k_0}
\frac{
d(p_0-k_0) + d(p_0+k_0)
}{
p_0
}
\frac{\Delta_l(p_0)}
{
1 + \frac{\lambda \nu_F}{k_{F,l} a}
\frac{1}{k_0}
\sum_{p_0' \neq k_0}
\left[d(p_0'-k_0) - d(p_0'+k_0)\right]
}
$$

where:

- $l \in \{t,b\}$ labels the two layers (top/bottom)
- $k_0, p_0$ are fermionic Matsubara frequencies
- $\nu_F$ is the density of states at the Fermi level
- $k_{F,l}$ is the layer-dependent Fermi momentum (shifted in the FE state)
- $a$ is the microscopic length scale used in the bosonic sector
- $\lambda$ is the effective electron–boson coupling
- $d(q_0)$ is the **momentum-integrated bosonic propagator** (kernel)

The ordered FE phase enters through:

1. **Layer chemical potential splitting**
   $$
   \mu_l = \mu \pm \delta\mu,
   $$
   which induces layer-dependent $k_{F,l}$ (and $v_{F,l}$).

2. **Ordered-phase boson mass replacement**
   $$
   r \rightarrow 2|r|
   $$
   in the bosonic propagator.

---

## Bosonic Kernel

The pairing interaction is controlled by the frequency-dependent kernel $d(q_0)$ obtained by integrating out momentum in the bosonic propagator. Schematically,

$$
d(q_0)=\int_0^\Lambda dq \; D(q,q_0),
$$

where $D(q,q_0)$ contains the FE mass term and Landau damping (see Appendix B for the precise form used in the paper). This kernel enters Eq. (B23) through the combinations
$d(p_0-k_0)\pm d(p_0+k_0)$.

---

## Regularization / Gap Transformation

To remove the apparent divergence associated with the static ($q_0=0$) contribution, the code uses the standard transformation (Appendix B):

$$
\Delta_l(k_0) = \frac{k_0 \, \Phi_l(k_0)}{\tilde{\Sigma}_l(k_0)} ,
$$

where $\Phi_l$ is the pairing vertex and $\tilde{\Sigma}_l$ is the (renormalized) normal self-energy factor entering the linearized equation.

---

## Numerical Method

At a given control parameter $r$ (distance to the FE QCP), the script:

1. Sets up a Matsubara grid:
   $$
   k_0=(2n+1)\pi T
   $$
   with UV cutoff of order $E_F$ (see code parameters).

2. Builds the kernel $d(q_0)$ on the Matsubara grid.

3. Constructs the linearized pairing operator from Eq. (B23) separately for $l=t$ and $l=b$.

4. Solves the resulting eigenvalue problem. The transition temperature is determined by the condition:
   $$
   \lambda_{\max}(T_c)=1,
   $$
   where $\lambda_{\max}$ is the largest eigenvalue of the linearized gap kernel.

5. Extracts $T_{c,t}$ and $T_{c,b}$ and then computes $\Delta T_c$.

---

## Outputs

Depending on the settings inside the script, typical outputs include:

- $T_{c,t}(r)$ and $T_{c,b}(r)$ (layer-resolved critical temperatures)
- $\Delta T_c(r)$ and/or $\Delta T_c/T_c$
- Plots/data used to generate the ordered-phase splitting features in **Figure 3**

---

## Usage

Run directly in MATLAB from the folder containing the script:

```matlab
check_splitting

