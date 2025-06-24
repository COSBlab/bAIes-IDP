#!/usr/bin/bash
# takes a pdb model as input and makes gromacs files.
name=idp

echo 6 | gmx pdb2gmx -f ${1} -water none -o ${name}.gro -p ${name}.top -i ${name}.itp -ignh
echo 0 | gmx trjconv -f ${name}.gro -s ${name}.gro -o ${name}.pdb
