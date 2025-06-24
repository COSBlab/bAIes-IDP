#!/usr/bin/bash
# takes GROMACS and proprocessing files and outputs all the input files for the bAIes simulations.
# `conda activate baies` or activate the environment that contains intermol

# inputs
gro=${1}
pdb=${2}
top=${3}

# file containing the residue-specific dihedral correction maps.
cmap=cmap_20240524.cmap


# generic name for intermediary and output files
name=idp

# gromacs to lammps conversion (`conda activate baies`):
python -m intermol.convert --gro_in ${gro} ${top} --lammps >> ${name}_conversion.log


# Force field modification amber99SB-ILDN -> Random coil
./remove_nonbonded_cmap.py -i ${name}_converted.input \
                                  -top ${name}_converted.lmp \
                                  -pdb ${pdb} -cmap ${cmap} \
                                  -oin ${name}_nvt.in \
                                  -otop ${name}_nvt.data \
                                  -cube 200.0 \
                                  -oxtc traj_${name}.xtc

# remove useless intermediary files:
rm ${name}_converted.input ${name}_converted.lmp ${name}_conversion.log
