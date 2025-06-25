# bAIes-IDP - Ensembles of disordered proteins and multidomain proteins with Alphafold-2

<p align="center">
  <img src="baies-idp.png" width="1000">
</p>

Atomic resolution ensemble predictions of Intrinsically Disordered and Multi-domain Proteins with Alphafold-2

Here you can find scripts and tutorials to perform ensemble prediction of IDPs using Alphafold-2, as introduced in:

V. Schnapka, T. Morozova, S. Sen, M. Bonomi. Atomic resolution ensembles of intrinsically disordered and multi-domain proteins with Alphafold. BioRxiv (2025). doi: [https://doi.org/10.1101/2025.06.18.660298](https://doi.org/10.1101/2025.06.18.660298)

This repository is organized in the following directories:
* `scripts`: python scripts used for preprocessing and preparations of the bAIes simulations.
* `tutorials`: complete tutorials for IDP ensemble preparation.
* `benchmark`: The input files to reproduce our simulations.
* `installation`: Some useful files to install the necessary software to run bAIes-IDP.

## Software requirements

To run bAIes-IDP, you will need a linux workstation, a python environment and some softwares:

### Alphafold-2

You can get Alphafold-2 [here](https://github.com/google-deepmind/alphafold).

Alternatively, if you cannot run Alphafold-2 locally, you can use a version of [Colabfold](https://github.com/sokrypton/ColabFold) that outputs the distance distributions. An implementation can be found [here](https://github.com/zshengyu14/ColabFold_distmats/blob/main/AlphaFold2.ipynb).

### A conda environment containing the intermol library

You can easily install this environment with conda by using the provided yml file (see `installation`) and running the following in a terminal:

`conda env create -f baies.yml`

Alternatively, the intermol library can be found [here](https://github.com/shirtsgroup/InterMol)

### GROMACS

GROMACS can be downloaded and installed from [here](https://manual.gromacs.org/current/download.html)

### LAMMPS with PLUMED

LAMMPS version 2 Aug. 2023 source code can be downloaded [here](https://download.lammps.org/tars/index.html)

For bAIes, LAMMPS must be patched with the file `patch_cmap.txt` provided in `installation`. After downloading the source code of LAMMPS, go in the source code main directory and run:

     `patch ./src/MOLECULE/fix_cmap.cpp < patch_cmap.txt`

Then, LAMMPS can be compiled using CMake (described [here](https://docs.lammps.org/Build_cmake.html)) or using make (described [here](https://docs.lammps.org/Build_make.html)).

The implementation of PLUMED is described [here](https://docs.lammps.org/Build_extras.html#plumed).

We recommend the use of OpenMP for parallelization.
