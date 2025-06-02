# variationalMC

C++ implementation of VMC for spin-1/2, spin-1, and SU(3) lattice models. This software uses the Slater-Jastrow approach as pioneered by Sorella and Becca, which can compute and optimize the wavefunctions as classified by the Projective Symmetry Group (PSG) approach by Xiao-Gang Wen. In addition to Wen's quantum orders, we can compute optimized states with magnetic long-range order using the Jastrow factor and its associated variational parameters. See Chapters 3 and 6 of [my thesis](https://repository.rice.edu/items/b1646003-2c86-44d0-a61e-e2fe3007b9e3) and the references therein for additional details on the theory and implementation, as well as an example of how to use this software to conduct quantum spin liquid research.

### Documentation

- (Work in progress...) [Documentation Homepage](https://butchertx.github.io/variationalMC/index.html)

### Recent Updates

This project is currently undergoing a major upgrade to the build, testing, and documentation. The `main` branch has all the legacy functionality intact. Meanwhile, the `dev` branch is getting CMake integration, Sphinx/Doxygen documentation (see the link above), a test suite using Google Test, and additional implementations to handle all the possible Slater-Jastrow wavefunctions for spin-1/2 and spin-1 systems.

### Quick Links

- [Installation instructions](https://butchertx.github.io/variationalMC/installation/index.html): See this page for step-by-step instructions on installing the software.
- [Git and Github tutorial](https://butchertx.github.io/variationalMC/contributing/git_tutorial.html): See this page if you're not already comfortable using git, branches, and creating pull requests.
