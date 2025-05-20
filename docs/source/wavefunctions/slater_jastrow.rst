

Slater-Jastrow Wavefunctions
#################################################

The Slater-Jastrow wavefunction is a product of a Slater determinant and a Jastrow factor. The Slater determinant is constructed from single-particle orbitals, which are determined by the chosen mean-field ansatz and the associated variational parameters. The Jastrow factor is a correlation factor that accounts for additional long-range correlations and onsite quantum fluctuations.

The wavefunction is given by:

.. math::

    |\Psi\rangle = \hat{P}^G \mathcal{J} |\psi_{MF}\rangle

The Gutzwiller projection :math:`\hat{P}^G`
=============================================================

The Gutzwiller projector :math:`\hat{P}^G` simply projects out the unphysical states of the wavefunction, which are those with double occupancy or onsite vacancies. This projection is handled exactly by the Monte Carlo updates, which only sample over states in the computational basis:

.. math::

    |x\rangle = |(r_1,m_1),(r_2,m_2),\ldots,(r_N,m_N)\rangle

where :math:`r_i` is the position of the :math:`i^{\text{th}}` lattice site and :math:`m_i` is the magnetic quantum number of the particle on that site. The Gutzwiller projector is then given by:

.. math::

    \hat{P}^G = \sum_{m_1 = -S}^{S}\cdots \sum_{m_N = -S}^{S} |m_1, m_2, \ldots, m_N\rangle \langle m_1, m_2, \ldots, m_N|

which is a sum over all possible single-occupancy spin-:math:`S` states. The notation of the basis states is the same as that of the computational basis for qudits with :math:`d = 2S + 1`, where the lattice indices are dropped because they are implicit in the Gutzwiller projector. This projection is applied exactly by the Monte Carlo updates, and does not have to be explicitly included in the wavefunction. 

The Jastrow factor :math:`\mathcal{J}`
=============================================================

The Slater Determinant :math:`|\psi_{MF}\rangle`
=============================================================

The Slater determinant is constructed from single-particle orbitals, which are determined by the chosen mean-field ansatz and the associated variational parameters. The Slater determinant is given by:

.. math::

    \psi(x) = \langle x | \psi_{MF}\rangle = \text{det}
    \begin{pmatrix}
        \phi_1(x_1) & \phi_2(x_1) & \cdots & \phi_N(x_1) \\
        \phi_1(x_2) & \phi_2(x_2) & \cdots & \phi_N(x_2) \\
        \vdots & \vdots & \ddots & \vdots \\
        \phi_1(x_N) & \phi_2(x_N) & \cdots & \phi_N(x_N)
    \end{pmatrix},

where the :math:`\phi_i(x_j)` are the occupied single-particle states of the mean-field Hamiltonian :math:`H_{MF}`. :math:`H_{MF}` encodes the chosen ansatz and the associated variational parameters. It takes the form

.. math::

    H_{MF} = \sum_{i,j,\sigma} t_{ij}^\sigma c_{i\sigma}^\dagger c_{j\sigma} + h\sum_{i\sigma\sigma'} d_i^{\sigma\sigma'}c_{i\sigma}^\dagger c_{i\sigma'} + \sum_{ij\sigma\sigma'} \Delta_{ij}^{\sigma\sigma'} c^\dagger_{i\sigma} c^\dagger_{j\sigma'} + \text{h.c.}.

For more details on the mean-field Hamiltonian, see :doc:`mean_field`.


Constructing the Slater determinant
=============================================================

Diagonalizing :math:`H_{MF}` gives the available single-particle orbitals:

.. math::

    \Phi = \begin{pmatrix}
        \phi_1(r_1, m_1) & \phi_2(r_1, m_1) & \cdots & \phi_{dN}(r_1, m_1) \\
        \phi_1(r_2, m_2) & \phi_2(r_2, m_2) & \cdots & \phi_{dN}(r_2, m_2) \\
        \vdots & \vdots & \ddots & \vdots \\
        \phi_1(r_{dN}, m_{dN}) & \phi_2(r_{dN}, m_{dN}) & \cdots & \phi_{dN}(r_{dN}, m_{dN})
    \end{pmatrix}

Note that this matrix has dimensions :math:`dN \times dN`, where :math:`d` is the local dimension of the spins and :math:`N` is the number of particles/sites. This list of orbitals is assumed to be sorted in order of increasing mean-field energy, so we will only use elements from the first :math:`N` orbitals to construct the Slater determinant.

To combine the available orbitals with a basis state and initialize the Slater determinant, we first choose a basis state :math:`|x\rangle` from the computational basis. In general this is done by randomly assigning a spin :math:`m_i` to each lattice site, while ensuring the total-:math:`S^z = \sum_i m_i` is zero. Given a configuration :math:`|x\rangle`, the Slater determinant is then the :math:`N` rows of :math:`\Phi` that correspond to the values of :math:`m_i` in the basis state. For example, suppose we have a spin-1/2 system with 2 lattice sites. Suppose the chosen initial configuration is :math:`|x\rangle = |\uparrow, \downarrow\rangle`. The available single-particle orbitals are given by:

.. math::

    \Phi = \begin{pmatrix}
        \phi_1(r_1, \uparrow) & \phi_2(r_1, \uparrow) & \phi_3(r_1, \uparrow) & \phi_4(r_1, \uparrow)\\
        \phi_1(r_2, \uparrow) & \phi_2(r_2, \uparrow) & \phi_3(r_2, \uparrow) & \phi_4(r_2, \uparrow)\\
        \phi_1(r_1, \downarrow) & \phi_2(r_1, \downarrow) & \phi_3(r_1, \downarrow) & \phi_4(r_1, \downarrow)\\
        \phi_1(r_2, \downarrow) & \phi_2(r_2, \downarrow) & \phi_3(r_2, \downarrow) & \phi_4(r_2, \downarrow)
    \end{pmatrix}

By convention, the lattice sites are grouped so that the full spatial wavefunction for a given value of :math:`m_i` takes up :math:`N` consecutive rows in each column. From the available orbitals above, the Slater determinant corresponding to the configuration :math:`|x\rangle = |\uparrow, \downarrow\rangle` is given by:

.. math::

    \psi(x) = \text{det}
    \begin{pmatrix}
        \phi_1(r_1, \uparrow) & \phi_2(r_1, \uparrow)\\
        \phi_1(r_2, \downarrow) & \phi_2(r_2, \downarrow)
    \end{pmatrix}.


Computing a proposed Monte Carlo update probability
=============================================================

Suppose we want to compute the probability of swapping the two particles in the above example. This probability is given by the ratio of the new and old Slater determinants:

.. math::

    P(x\rightarrow x') = \frac{\psi(x')}{\psi(x)} = \frac{\text{det}
    \begin{pmatrix}
        \phi_1(r_1, \downarrow) & \phi_2(r_1, \downarrow)\\
        \phi_1(r_2, \uparrow) & \phi_2(r_2, \uparrow)
    \end{pmatrix}}{\text{det}
    \begin{pmatrix}
        \phi_1(r_1, \uparrow) & \phi_2(r_1, \uparrow)\\
        \phi_1(r_2, \downarrow) & \phi_2(r_2, \downarrow)
    \end{pmatrix}}.

In this small example, the full determinant could be computed at each Monte Carlo step. However, for larger systems, the determinant is too large to compute directly. Instead, the primary object we work with is the matrix :math:`W(x)`, defined as

.. math::

    W(x) = \Phi \cdot I_{dN\times N} \cdot (S(x))^{-1}

where :math:`I_{dN\times N}` is the identity matrix of shape :math:`(dN, N)`, which selects only the first :math:`N` columns of :math:`\Phi`, i.e. the occupied orbitals. :math:`S(x)` is the Slater matrix, which is just the matrix of occupied orbitals that forms the Slater determinant, such that :math:`\psi(x) = \text{det}(S(x))`. 

With this definition in mind, let us now rewrite the ratio of determinants in a more convenient form. Suppose the update :math:`x\rightarrow x'` is a rank-2 update generated by swapping the particles at two sites. Then, the new Slater matrix can be written as

.. math::

    S(x') = S(x) + UV

where :math:`U` is a :math:`(N, 2)` permutation matrix of the form

.. math::

    U = \begin{pmatrix}
        0 & 0 \\
        \vdots & \vdots \\
        1 & 0 \\
        \vdots & \vdots \\
        0 & 1 \\
        \vdots & \vdots \\
        0 & 0
    \end{pmatrix}

The only nonzero entries of :math:`U` are the rows indexed by :math:`i_1` and :math:`i_2`, which are the rows of :math:`S(x)` corresponding to the two particles being swapped. The matrix :math:`V` is a :math:`(2, N)` matrix of the differencing the new and old particle occupations, which are rows of :math:`S(x')` and :math:`S(x)` respectively. For example, if the two particles being swapped are at sites :math:`i_1` and :math:`i_2`, then :math:`V` is given by

.. math::

    V = \begin{pmatrix}
        \phi_1(x_{i_2}) - \phi_1(x_{i_1}) & \phi_2(x_{i_2}) - \phi_2(x_{i_1}) & \cdots & \phi_N(x_{i_2}) - \phi_N(x_{i_1}) \\
        \phi_1(x_{i_1}) - \phi_1(x_{i_2}) & \phi_2(x_{i_1}) - \phi_2(x_{i_2}) & \cdots & \phi_N(x_{i_1}) - \phi_N(x_{i_2})
    \end{pmatrix}.

The notation :math:`x_{i}` refers to the configuration of the :math:`i^{\text{th}}` particle in the Slater determinant, i.e. :math:`x_{i} = (r_i, m_i)`.

With the notation in place, the update probability can be rewritten as

.. math::

    P(x\rightarrow x') = \text{det}(S(x')\cdot S^{-1}(x)) = \text{det}((S(x) + UV)\cdot S^{-1}(x)) = \text{det}(I_{N\times N} + UV\cdot S^{-1}(x)).

in terms of :math:`W(x)`:

.. math::
    
    P(x\rightarrow x') = \frac{\text{det}(S(x'))}{\text{det}(S(x))}
    = \text{det}(S(x')\cdot S^{-1}(x))


.. Instead, we use the following identity:
    .. math::

        \frac{\text{det}(A)}{\text{det}(B)} = \frac{\text{det}(A + UV)}{\text{det}(B + UV)}\cdot\frac{\text{det}(B)}{\text{det}(A)}.

    where :math:`U` and :math:`V` are matrices giving a rank-:math:`k` update to the Slater matrix :math:`S(x)`.

Updating the Slater determinant
=============================================================

We can use the Woodbury matrix identity to compute the ratio of determinants when only changing a small number of rows and columns. The Woodbury matrix identity is:

.. math::

    (A + UV)^{-1} = A^{-1}U(I + VA^{-1}U)^{-1}V A^{-1}.