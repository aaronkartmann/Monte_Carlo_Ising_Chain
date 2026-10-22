# Monte Carlo Simulation of the 1D Ising Model

## Overview

This project implements a Monte Carlo simulation of the one-dimensional Ising model using the Metropolis algorithm. The notebook `IsingChain_ExampleSimulation.ipynb` demonstrates the simulation by investigating how temperature affects thermodynamic observables of a spin chain and comparing their Monte Carlo estimates with analytical predictions.

This project is based on the Week 1 tutorial for the Computational Quantum Physics course at Radboud University Nijmegen. It extends the original exercise by providing detailed analytical derivations in `Derivations.ipynb` and implementing the Ising model as a reusable Python class in `IsingChain.py`, alongside additional analysis of the simulation results.

## Physical Model

The system consists of \(N\) interacting spins, where each spin can take the value \(s_i = \pm 1\). The model uses open boundary conditions, nearest-neighbor interactions, and no external magnetic field.

The energy of a spin configuration is
\[
E = -J \sum_{i=1}^{N-1} s_i s_{i+1},
\]
where \(J=1\) and the Boltzmann constant is set to \(k_B=1\).

The Metropolis algorithm generates configurations according to the Boltzmann distribution. Thermodynamic observables are calculated from the sampled configurations and compared with analytical results.

## Features

- Monte Carlo simulation of Ising chains using the Metropolis algorithm.
- Calculation of energy, magnetization, heat capacity, and magnetic susceptibility.
- Statistical uncertainty analysis using binning.
- Comparison of Monte Carlo estimates with analytical predictions.
- Investigation of temperature-dependent behavior and sampling limitations.

## Project Structure

- `IsingChain.py`: Implements the Ising chain and Metropolis algorithm.
- `analysis.py`: Contains analysis functions, including binning and analytical calculations.
- `Derivations.ipynb`: Derives the expressions used in the project.
- `IsingChain_ExampleSimulation.ipynb`: Example simulation: runs the simulation for an \(N=100\) spin chain and analyzes the results.

## Requirements

- Python 3
- NumPy
- Matplotlib
- Jupyter Notebook

## How to Run

1. Clone or download this repository.
2. Install the required packages.
3. Open `IsingChain_ExampleSimulation.ipynb` and run the cells in order.

The notebook imports functions and classes from `IsingChain.py` and `analysis.py`. Keep these files in the same directory as the notebook.
