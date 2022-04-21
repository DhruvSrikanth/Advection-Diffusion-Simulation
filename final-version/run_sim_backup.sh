#!/bin/bash

# module load openmpi

mpi++ ./final-version/advection_simulation.cpp -o ./final-version/advection_simulation -O3 -ffast-math -mtune=native -march=native -fopenmp -lmpi

mpiexec --bind-to none -np 4 ./final-version/advection_simulation 400 20000 1.0 1.0e6 5.0e-7 2.85e-7 2 "LAX"

rm ./final-version/advection_simulation

python3 ./final-version/generate_plots.py