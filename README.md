# bAIes-IDP
Atomic resolution ensemble predictions of Intrinsically Disordered and Multi-domain Proteins with Alphafold-2

Here you can find scripts and tutorials to perform ensemble prediction of IDPs using Alphafold-2, as introduced in:

V. Schnapka, T. Morozova, S. Sen, M. Bonomi. Atomic resolution ensembles of intrinsically disordered and multi-domain proteins with Alphafold. BiorXiv (2025). doi: [https://doi.org/10.1101/2025.06.18.660298](https://doi.org/10.1101/2025.06.18.660298)

This repository is organized in the following two directories:
* `scripts`: python scripts used for preprocessing and preparations of the bAIes simulations.
* `tutorials`: complete tutorials for IDP ensemble preparation.
* `benchmark`: The input files to reproduce our simulations.

## Software requirements

You will need several tools and a python environment:

### A conda environment containing the intermol library

You can install this environment with conda by using the provided yml file and running the following in a terminal:

`conda env create -f baies.yml`

### GROMACS

GROMACS can be downloaded and installed from [here](https://manual.gromacs.org/current/download.html)

## LAMMPS with PLUMED

LAMMPS version 2 Aug. 2023 source code can be downloaded [here](https://download.lammps.org/tars/index.html)

For bAIes, LAMMPS must be patched with the file `patch_cmap.txt` provided here. After downloading the source code of LAMMPS, go in the source code main directory and run:

     `patch ./src/MOLECULE/fix_cmap.cpp < patch_cmap.txt`

Then, LAMMPS can be compiled using CMake (described [here](https://docs.lammps.org/Build_cmake.html)) or using make (described [here](https://docs.lammps.org/Build_make.html)).

The implementation of PLUMED is described [here](https://docs.lammps.org/Build_extras.html#plumed).

We recommend the use of OpenMP for parallelization.
