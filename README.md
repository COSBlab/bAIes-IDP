# bAIes-IDP
Atomic resolution ensemble predictions of Intrinsically Disordered and Multi-domain Proteins with Alphafold-2

Here you can find scripts and tutorials to perform ensemble prediction of IDPs using Alphafold-2, as introduce in:

(Publication)

This repository is organized in the following two directories:
* `scripts`: python scripts used for preprocessing and preparations of the bAIes simulations
* `tutorials`: complete tutorials for IDP ensemble preparation

## LAMMPS with PLUMED installation

LAMMPS version 2 Aug. 2023 source code can be downloaded [here](https://download.lammps.org/tars/index.html)

For bAIes, LAMMPS must be patched with the file `patch_cmap.txt` provided here. After downloading the source code of LAMMPS, go in the source code main directory and run:

     `patch ./src/MOLECULE/fix_cmap.cpp < patch_cmap.txt`

Then, LAMMPS can be compiled using CMake (described [here](https://docs.lammps.org/Build_cmake.html)) or using make (described [here](https://docs.lammps.org/Build_make.html)).

The implementation of PLUMED is described [here](https://docs.lammps.org/Build_extras.html#plumed).

## GROMACS installation

GROMACS can be downloaded and installed from [here](https://manual.gromacs.org/current/download.html)
