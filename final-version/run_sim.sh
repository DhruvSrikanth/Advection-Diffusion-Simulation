#!/bin/bash

# module load openmpi

g++-11 ./final-version/advection_simulation.cpp -o ./final-version/advection_simulation -I /opt/homebrew/Cellar/open-mpi/4.1.3/include -O3 -ffast-math -mtune=native -march=native -fopenmp -L /opt/homebrew/Cellar/open-mpi/4.1.3/lib -L /opt/homebrew/opt/libevent/lib -lmpi

mpiexec --bind-to none -np 4 ./final-version/advection_simulation 400 20000 1.0 1.0e6 5.0e-7 2.85e-7

rm ./final-version/advection_simulation

python3 ./final-version/generate_plots.py