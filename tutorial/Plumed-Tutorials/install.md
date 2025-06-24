## Hardware requirements

To complete this tutorial, you need:

* A linux workstation;
* Ideally, a multicore CPU (8 cores is sufficient).

## Software requirements

Make sure you have access to at least one of the following:
* Local Alphafold-2 installation. You can get Alphafold-2 [here](https://github.com/google-deepmind/alphafold);
* Colabfold. You can use Colabfold [here](https://github.com/sokrypton/ColabFold).

Make sure you have installed in your system:
* [Conda](https://www.anaconda.com) to install the environment needed by the conversion scripts
* [GROMACS](https://www.gromacs.org)
* [LAMMPS (2 Aug. 2023)](https://download.lammps.org/tars/index.html) with [PLUMED](https://www.plumed.org): see installation instructions below.

## Setting up the Conda environment for the LAMMPS input preparation

The conversion to LAMMPS requires the use of the python library [intermol](https://github.com/shirtsgroup/InterMol).
You can make this environment by running on the terminal:

`conda env create -f baies.yml`

Where `baies.yml` is provided in the installation directory. The content is:

```name: baies
channels:
  - conda-forge
  - omnia
  - defaults
dependencies:
  - _libgcc_mutex=0.1=main
  - _openmp_mutex=5.1=1_gnu
  - ca-certificates=2024.3.11=h06a4308_0
  - ld_impl_linux-64=2.38=h1181459_1
  - libffi=3.4.4=h6a678d5_0
  - libgcc-ng=11.2.0=h1234567_1
  - libgomp=11.2.0=h1234567_1
  - libstdcxx-ng=11.2.0=h1234567_1
  - ncurses=6.4=h6a678d5_0
  - openssl=3.0.13=h7f8727e_0
  - pip=23.3.1=py38h06a4308_0
  - python=3.8.19=h955ad1f_0
  - readline=8.2=h5eee18b_0
  - setuptools=68.2.2=py38h06a4308_0
  - sqlite=3.41.2=h5eee18b_0
  - tk=8.6.12=h1ccaba5_0
  - wheel=0.41.2=py38h06a4308_0
  - xz=5.4.6=h5eee18b_0
  - zlib=1.2.13=h5eee18b_0
  - pip:
      - intermol==0.1.0.dev0
      - numpy==1.24.4
      - par==1.3.4
      - parm==1.11
      - parmed==3.4.4
      - six==1.16.0
```

### GROMACS

GROMACS can be downloaded and installed from [here](https://manual.gromacs.org/current/download.html)

## Setting up LAMMPS with PLUMED

LAMMPS version 2 Aug. 2023 source code can be downloaded [here](https://download.lammps.org/tars/index.html)

For bAIes, LAMMPS must be patched with the file `patch_cmap.txt` provided in `installation`. After downloading the source code of LAMMPS, go in the source code main directory and run:

     `patch ./src/MOLECULE/fix_cmap.cpp < patch_cmap.txt`

Then, LAMMPS can be compiled using CMake (described [here](https://docs.lammps.org/Build_cmake.html)) or using make (described [here](https://docs.lammps.org/Build_make.html)).

The implementation of PLUMED is described [here](https://docs.lammps.org/Build_extras.html#plumed).

We recommend the use of OpenMP for parallelization.


