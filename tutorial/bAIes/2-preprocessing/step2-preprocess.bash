#!/usr/bin/bash
# takes AF2 and topology information and outputs plumed files for bAIes simulations

# inputs
mdpdb=${1}
pdb=${2}
dist=${3}

# outputs
ndx=atom_list.ndx
dat=baies_params.dat
out=plumed.dat

# noise model for bAIes. JEFFREYS for bAIes, NONE for bAIes-N.
prior=JEFFREYS


# preprocess: generate plumed file
./preprocess_bAIes.py -pdb ${pdb}   \
    -mdpdb ${mdpdb} \
    -pkl ${dist}    \
    -out ${dat}   \
    -model gauss \
    -cutoff matrix \
    -ndxout ${ndx} \
    --verbose
    
    
# Generate the PLUMED file to give to the simulation engine:
touch ${out}
echo "#MOLINFO STRUCTURE=${mdpdb}" >> ${out}
echo "batoms: GROUP NDX_FILE=${ndx} NDX_GROUP=batoms" >> ${out}
echo "baies: BAIES ATOMS=batoms DATA_FILE=${dat} PRIOR=${prior} TEMP=2.478541306" >> ${out}
echo "PRINT ARG=baies.ene FILE=COLVAR STRIDE=500" >> ${out}
echo "bbias: BIASVALUE ARG=baies.ene STRIDE=2" >> ${out}



