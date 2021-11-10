#!/bin/bash

module purge
ml icc/2019.3.199-GCC-8.3.0 impi/2019.4.243 FFTW/3.3.8 

icpc -c -std=c++20 ../src/vmctype.cpp
icpc -c -std=c++20 ../src/vmc_io.cpp
icpc -c -std=c++20 ../src/Lattice.cpp
icpc -c -std=c++20 ../src/MeanFieldAnsatz.cpp
icpc -c -std=c++20 ../src/Wavefunction.cpp
icpc -c -std=c++20 ../src/ProjectedState.cpp
icpc -c -std=c++20 ../src/VariationalMonteCarlo.cpp
icpc -c -std=c++20 ../src/Lai_CSL.cpp
icpc -L$LD_LIBRARY_PATH -I$CPATH -lgsl -lgslcblas vmctype.o vmc_io.o Lattice.o MeanFieldAnsatz.o Wavefunction.o ProjectedState.o VariationalMonteCarlo.o Lai_CSL.o -lfftw3 -o ../bin/lai_csl.exe
