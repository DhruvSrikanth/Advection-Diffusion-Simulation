#!/bin/bash

mkdir ./visualizing/outputs

g++-11 ./visualizing/advection_simulation.cpp -o ./visualizing/advection_simulation -fopenmp -O5

echo Enter the number of timesteps:
read NT
./visualizing/advection_simulation 400 $NT 1.0 1.0e6 5.0e-7 2.85e-7

rm ./visualizing/advection_simulation

python3 ./visualizing/generate_plots.py

rm -rf ./visualizing/outputs/

echo "Done!"