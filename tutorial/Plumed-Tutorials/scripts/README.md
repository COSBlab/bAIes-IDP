# Scripts for preprocessing and preparation for bAIes ensemble predictions

Here you will find the necessary python scripts to prepare bAIes on your system. 
For more info about the python scripts, just type:

`python script.py -h`

For more details, you can have a look at the tutorials.

Brief overview:
* `preprocess_bAIes.py`: Python script that reads the AF2 pdb model and distogram file and generates the necessary plumed files for the bAIes sampling.
* `remove_nonbonded_cmap_plumed.py`: Python script that performs the force field simplification and generate the final LAMMPS files for bAIes.
* `remove_nonbonded_cmap.py`: Python script that performs the force field simplification and generate the final LAMMPS files for the random coil simulations.

## **Software installation**

To perform step 3, you need a specific conda environment.
To create it on your system, you can use the following file:
* `baies.yml`: file containing all the information about python libraries and versions to reproduce the python environment used for bAIes;

To create this environment, just run:

`conda env create -f baies.yml`

To activate the environment in your terminal:

`conda activate baies`
