.. variationalMC documentation master file, created by
   sphinx-quickstart on Wed Apr 16 21:02:47 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

variationalMC documentation
===========================

``variationalMC`` is a C++ library for performing Variational Monte Carlo (VMC)
simulations. It is primarily designed for studying the spin-1/2 and spin-1
Heisenberg models in one, two, and three dimensions. It uses Slater-Jastrow
wavefunctions with an exact Gutzwiller projection, and Stochastic Reconfiguration
for variational parameter optimization. However, the Monte Carlo and optimization
routines are general, and the library can be easily extended to be used with
addiitonal wavefunction types.

- To get started with installation, see :doc:`the Installation Guide </installation/index>`.
- To learn about what the library can do, see :doc:`the Project Overview </overview/index>`.

.. toctree::
   :maxdepth: 2
   :caption: Full Table of Contents:

   installation/index
   overview/index
   examples/index
   lattices/index
   wavefunctions/index
   models/index
   montecarlo/index
   project_spec/index
   testing/index
   contributing/index
   references/index