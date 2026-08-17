
This code calculates the anharmonic part of the Helmholtz free energy, internal energy, entropy and heat capacity. It includes the contributions from baseline shifts, renormalized phonons and higher order terms. This program must use force constants from sTDEP/SSCHA (self-consistent phonons). Force cosntants fit to MD/PIMD data are incompatible with this theory.

Only the first and second-cumulant corrections are computed by this code. The constant corrections (i.e., $U_0$ and its derivatives) must be implemented elsewhere as it requires interfacing with a force-calculator like LAMMPS. Please refer to [^Meitz2026] and a reference implementation [here](https://github.com/ejmeitz/CrystalCumulants.jl).

### Longer summary

Per the (quasi)harmonic approximation the free energy is determined as

$$
F(T,V) = U(V)+F_{ph}(T,V)
$$

where $U$ is the energy of the static lattice, and $F_{ph}$ is the free energy of the phonons. The only temperature dependence enters as occupation numbers in the phonon free energy expression. In the TDEP formalism it is only slightly more involved.

As established in [extract_forceconstants.md](extract forceconstants) the TDEP free energy, to lowest order is given by

$$
F(T,V) = U_0(T,V)+F_{ph}(T,V)
$$

where $F_{ph}(T,V)$ is the free energy of the (effective) phonons, and $U_0$ is a renormalized baseline energy, given by

$$
U_0 = \left\langle
	U - \frac{1}{2}\sum_{ij} u_i\Phi_{ij}u_j -
	\frac{1}{2}\sum_{ij} u_i \Phi^{\text{polar}}_{ij}u_j
	\right\rangle
$$

that is, the potential energy from the molecular dynamics/stochastic sampling minus the potential energy of the force constant model. It's important to remember to also subtract the electrostatic energy since it is included in $F_{ph}$.

#### Higher order corrections

The free energy can be improved by explicitly including the contribution from the higher order terms in the Hamiltonian. The free energy then reads

$$
\begin{equation}
    F = U_0 + F_{\textrm{ph}} + \Delta F^{3\textrm{ph}} + \Delta F^{4\textrm{ph}}
\end{equation}
$$

where $F_{\textrm{ph}}$ is the phonon free energy and $U_0$ the renormalized baseline free energy. However, since the derivations of $\Delta F^{3\mathrm{ph}}$ and $\Delta F^{4\mathrm{ph}}$ assume self-consistent phonons the average is with respect to harmonic dynamics (i.e., $\langle \cdot \rangle_\mathrm{H}$). The expression is given by

$$
\begin{equation}
    U_0 = \left\langle U -
    \frac{1}{2!}\sum_{\substack{ ij\\ \alpha\beta } }\overset{\textrm{lr}}{\Phi}_{ij}^{\alpha\beta}
u_i^\alpha u_j^\beta -
    \frac{1}{2!}\sum_{\substack{ ij\\ \alpha\beta } }\Phi_{ij}^{\alpha\beta}
u_i^\alpha u_j^\beta -
 \frac{1}{3!}
\sum_{\substack{ijk\\ \alpha\beta\gamma}}\Phi_{ijk}^{\alpha\beta\gamma}
u_i^\alpha u_j^\beta u_k^\gamma -
\frac{1}{4!}
	\sum_{\substack{
	ijkl\\
	\alpha\beta\gamma\delta
	}}
\Phi_{ijkl}^{\alpha\beta\gamma\delta}
u_i^\alpha u_j^\beta u_k^\gamma u_l^\delta
    \right\rangle_\mathrm{H}
\end{equation}
$$

Here it is important to note that we have to subtract the long-ranged polar interactions to avoid double-counting them. The explicit anharmonic contributions are given via[^Leibfried1961][^Cowley1963][^wallace1998thermodynamics]

$$
\begin{equation}%\label{eq:deltaF3}
	\Delta F^{3\textrm{ph}} = -\frac{\hbar^2}{48N}\sum_{\substack{\bm k\bm k^\prime \bm k^{\prime\prime} \\ \lambda\lambda^\prime\lambda^{\prime\prime}}}\frac{\left|\Phi_{\lambda \lambda^\prime \lambda^{\prime\prime}} ^ {\bm k \bm k^\prime \bm k^{\prime\prime}}\right|^2}{\omega\omega^\prime\omega^{\prime\prime}}  \left[\frac{(n+1)(n^\prime + n^{\prime\prime} + 1) + n^\prime n^{\prime\prime}}{\omega + \omega^\prime + \omega^{\prime\prime}}
    + 3\frac{n^{\prime\prime}(n+n^\prime+1)-nn^\prime}{\omega + \omega^\prime - \omega^{\prime\prime}}\right]
\end{equation}
$$

and

$$
\begin{equation}%\label{eq:deltaF4}
	\Delta F^{4\textrm{ph}} =
	\frac{\hbar^2}{8N}\sum_{\bm k\ \bm k^\prime \lambda\lambda^\prime}\frac{\Phi_{\lambda\lambda\lambda^\prime\lambda^\prime}^{\bm k,-\bm k,\bm k^\prime,-\bm k^\prime}}{\omega\omega^{\prime}}\left(n+\frac{1}{2}\right)\left(n^\prime + \frac{1}{2}\right)
\end{equation}
$$

[^Meitz2026]: [Meitz, E., Castellano, A., Wang, G. J. & McGaughey, A. J. H. (2026) Quantum Anharmonic Phonon Thermodynamics from the Free Energy Cumulant Expansion. npj Computational Materials (in review)](https://www.researchsquare.com/article/rs-10541555/v1)

[^Leibfried1961]: Leibfried, G. & Ludwig, W. (1961) Theory of Anharmonic Effects in Crystals. Solid State Phys 12, 275–444.

[^Cowley1963]: Cowley, R. A. (1963) The lattice dynamics of an anharmonic crystal. Adv Phys 12, 421–480.

[^wallace1998thermodynamics]: Wallace, D.C. Thermodynamics of Crystals. (Dover Publications).

