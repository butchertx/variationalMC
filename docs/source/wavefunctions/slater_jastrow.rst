

Slater-Jastrow Wavefunctions
#################################################

The Slater-Jastrow wavefunction is a product of a Slater determinant and a Jastrow factor. The Slater determinant is constructed from single-particle orbitals, which are determined by the chosen mean-field ansatz and the associated variational parameters. The Jastrow factor is a correlation factor that accounts for additional long-range correlations and onsite quantum fluctuations.

The wavefunction is given by:

.. math::

    |\Psi\rangle = \hat{P}^G \mathcal{J} |\psi_{MF}\rangle

The Gutzwiller projection :math:`\hat{P}^G`
============================================

The Gutzwiller projector :math:`\hat{P}^G` simply projects out the unphysical states of the wavefunction, which are those with double occupancy or onsite vacancies. This projection is handled exactly by the Monte Carlo updates, which only sample over states in the computational basis:

.. math::

    |x\rangle = |(r_1,m_1),(r_2,m_2),\ldots,(r_N,m_N)\rangle

where :math:`r_i` is the position of the :math:`i^{\text{th}}` lattice site and :math:`m_i` is the magnetic quantum number of the particle on that site. The Gutzwiller projector is then given by:

.. math::

    \hat{P}^G = \sum_{m_1 = -S}^{S}\cdots \sum_{m_N = -S}^{S} |m_1, m_2, \ldots, m_N\rangle \langle m_1, m_2, \ldots, m_N|

which is a sum over all possible single-occupancy spin-:math:`S` states. The notation of the basis states is the same as that of the computational basis for qudits with :math:`d = 2S + 1`, where the lattice indices are dropped because they are implicit in the Gutzwiller projector. This projection is applied exactly by the Monte Carlo updates, and does not have to be explicitly included in the wavefunction. 

The Jastrow factor :math:`\mathcal{J}`
========================================

The Slater Determinant :math:`|\psi_{MF}\rangle`
========================================

The Slater determinant is constructed from single-particle orbitals, which are determined by the chosen mean-field ansatz and the associated variational parameters. The Slater determinant is given by: