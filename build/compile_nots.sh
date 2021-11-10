#!/bin/bash

module purge
ml iccifort/2020.1.217 impi/2019.7.217 imkl/2020.1.217 

icpc -c -std=c++20 ../src/vmctype.cpp
icpc -c -std=c++20 ../src/vmc_io.cpp
icpc -c -std=c++20 ../src/Lattice.cpp
icpc -c -std=c++20 ../src/MeanFieldAnsatz.cpp
icpc -c -std=c++20 ../src/Wavefunction.cpp
icpc -c -std=c++20 ../src/ProjectedState.cpp
icpc -c -std=c++20 ../src/VariationalMonteCarlo.cpp
icpc -c -std=c++20 ../src/main.cpp
icpc -L$LD_LIBRARY_PATH -I$CPATH -lgsl -lgslcblas -lmkl -liomp5 -lpthread vmctype.o vmc_io.o Lattice.o MeanFieldAnsatz.o Wavefunction.o ProjectedState.o VariationalMonteCarlo.o main.o -o ../bin/vmc++.exe
