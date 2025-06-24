# Scripts for preprocessing and preparation for bAIes ensemble predictions

Here you will find some python and shell scripts useful to prepare bAIes sampling on your system. 
For more info about the python scripts, just type:

`python script.py -h`

For more details, you can have a look at the tutorials.

Brief overview:
* `step1-prepare_gmx.bash`: Two lines script that shows how to prepare the GROMACS inputs;
* `step2-preprocess.bash`: Script that reads Alphafold outputs and calls the `preprocess_bAIes.py` to generate the PLUMED files for the simulations;
* `preprocess_bAIes.py`: Python script that reads the AF2 pdb model and distogram file and generates the necessary plumed files for the bAIes sampling;
* `step3-conversion.bash`: Script that converts the inputs into the LAMMPS files necessary for the bAIes simulations;
* `remove_nonbonded_cmap_plumed.py`: Python script called by `step3-conversion.bash` that performs the force field simplification and generate the final LAMMPS files.
