#!/usr/bin/env python
import numpy as np
from argparse import ArgumentParser
from pickle import load
import re


def build_parser():
    parser = ArgumentParser(description="Make plumed files for random coil / bAIes simulations")
    parser.add_argument('-i', type=str, help="Input lammps file")
    parser.add_argument('-top', type=str, help="Data topology lammps file")
    parser.add_argument('-pdb', type=str, help="PDB file", default="protein.pdb")
    parser.add_argument('-cmap', type=str, help="CMAP data file", default="ff.cmap")
    parser.add_argument('-oin', type=str, help="Output lammps input file")
    parser.add_argument('-otop', type=str, help="Output topology lammps file", default="conff.data")
    parser.add_argument('-cube', type=float, help="Cube side size (Angstroms) default is 400.0", default=400.0)
    parser.add_argument('-oxtc', type=str, help="xtc file outputed by the lammps simulation", default="traj.xtc")
    parser.add_argument('-plumed', type=str, help="plumed file", default="plumed.dat")
    return parser
    
    
def read_pdb(input_pdb="protein.pdb"):
    """
    Reads a pdb file and outputs a dictionnary of cmaps indexes for lammps data file:
    {cmap_idx: ((restype, restype+1), C-1, N, CA, C, N+1)}
    """
    # Extract the data:
    with open(input_pdb, "r") as f:
        data = [el.rstrip("\n") for el in f.readlines() if el.split()[0] == 'ATOM']
    # Atom numbers in the PDB file:
    atom_idx = [(int(el[5:11])) for el in data]
    # Atom types in the PDB file:
    atom_type = [el[12:17].rstrip(" ").lstrip(" ") for el in data]
    # Residue type in the PBD file:
    res_type = [el[17:20] for el in data]
    # original residue indexes
    res_idx = [int(el[23:26]) for el in data]
    # Storing the atom number of CB (or CA if GLY) and the residue type for each residue:
    read_chain_idx = [el[21] for el in data]

    # create a consistent labelling of the protein chains from the pdb:
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    current_chain = None

    chain_N = 1
    chain_idx = []
    for i, el in enumerate(data):
        if i == 0:
            resN = res_idx[i]
            if read_chain_idx[i] == ' ':
                # create a chain label if no chain label specified
                chain_idx.append('X'+alphabet[chain_N-1])
            else:
                # take the specified chain label
                chain_idx.append(read_chain_idx[i])
            continue
        # if the next residue is smaller, or that the chain column is different from the previous one, it implies that we change chain.
        if res_idx[i] < resN or read_chain_idx[i] != read_chain_idx[i-1]:
            # changement de chaine
            chain_N += 1
        if read_chain_idx[i] == ' ':
            # create a chain label if no chain label specified
            chain_idx.append('X'+alphabet[chain_N-1])
        else:
            # take the specified chain label
            chain_idx.append(read_chain_idx[i])
        resN = res_idx[i]

    # get cmap by atom number:
    atoms_CA = {res_idx[i]: atom_idx[i] for i, el in enumerate(atom_type) if el == "CA"}
    atoms_C = {res_idx[i]: atom_idx[i] for i, el in enumerate(atom_type) if el == "C"}
    atoms_N = {res_idx[i]: atom_idx[i] for i, el in enumerate(atom_type) if el == "N"}
    restypes = {res_idx[i]: res_type[i] for i, el in enumerate(atom_type) if el == "CA"}
    cmaps = dict()
    nb_cmaps = 0
    for i, res in enumerate(res_idx):
        if atom_type[i] == "CA" and res in atoms_CA.keys() and res in atoms_N.keys() and res+1 in atoms_N.keys() and res in atoms_C.keys() and res-1 in atoms_C.keys():
            nb_cmaps += 1
            cmaps[nb_cmaps] = ((res, restypes[res], restypes[res+1]), atoms_C[res-1], atoms_N[res], atoms_CA[res], atoms_C[res], atoms_N[res+1])
    print(cmaps)
    return cmaps
    
    
def write_cmaps(cmaps, output_file):
    # cmap[typenb] = (title, bins, matrix)
    with open(output_file, "w") as f:
        f.write("# Backbone Phi/Psi Dihedrals correction map\n")
        f.write("\n")
        for typenb in cmaps.keys():
            title, bins, matrix = cmaps[typenb][0], cmaps[typenb][1], cmaps[typenb][2]
            f.write("# {}, type {}\n\n".format(title, typenb))
            for i, line in enumerate(matrix):
                f.write("# {} \n".format(bins[i]))
                for j, el in enumerate(line):
                    f.write("{:.6f} ".format(el))
                    if (j+1)%5 == 0 or (j+1) == matrix.shape[0]:
                        f.write("\n")
                f.write("\n")
    return 0
    
    
def read_input_lammps(filename, data_out_file, cmap_file, outxtc, plumedfile):
    two_pow1over6 = 2**(1/6)
    pair_coeffs = dict()
    output_data = []
    with open(filename, "r") as f:
        data = [el for el in f.readlines()]
    for i, line in enumerate(data):
        if re.match("special_bonds", line):
            continue
        if re.match("read_data", line):
            #output_data.append("read_data {}\n".format(data_out_file))
            continue
        elif re.match("pair_style", line):
            output_data.append("pair_style lj/cut 2.0\n")
            # cmap:
            output_data.append("\n")
            output_data.append("fix drycmap all cmap {}\n".format(cmap_file))
            output_data.append("read_data {} fix drycmap crossterm CMAP\n".format(data_out_file))
            output_data.append("fix_modify drycmap energy yes\n\n")
            #output_data.append("read_data {} add append fix cmap crossterm CMAP\n".format(data_out_file))
            continue
        elif re.match("pair_coeff", line):
            line_content = line.split()
            atom_i, atom_j = int(line_content[1]), int(line_content[2]) # j doesn't matter caus atom i is atom j in converted amber 99SB-ILDN
            eps_kcalpermol = float(line_content[3])
            sigma_dist = float(line_content[4])
            pair_coeffs[(atom_i, atom_j)] = (eps_kcalpermol, sigma_dist, two_pow1over6*sigma_dist)
            continue
        if re.match("pair_modify mix", line):
            # write the mixed pair coeffs:
            I = [el[0] for el in pair_coeffs.keys()]
            J = [el[0] for el in pair_coeffs.keys()]
            for i in I:
                for j in J:
                    #print(i, j)
                    if (i, j) in pair_coeffs.keys() or (j, i) in pair_coeffs.keys():
                        continue
                    if line.split()[-1] == "arithmetic":
                        eps = np.sqrt(pair_coeffs[(i, i)][0]*pair_coeffs[(j, j)][0]) # sqrt(eps_i*eps_j)
                        sigma = 0.5*(pair_coeffs[(i, i)][1] + pair_coeffs[(j, j)][1]) # 0.5*(sigma_i + sigma_j)
                        cutoff = two_pow1over6*sigma
                    elif line.split()[-1] == "geometric":
                        print("you need to implement the geometric pair combination yourself.")
                        return None
                    else:
                        print("Pair combination not recognized.")
                        return None
                    pair_coeffs[(i, j)] = (eps, sigma, cutoff)
            for el in pair_coeffs.keys():
                line = "pair_coeff {} {}   {:.7f}   {:.7f}   {:.7f}\n".format(el[0], el[1], pair_coeffs[el][0], pair_coeffs[el][1], pair_coeffs[el][2])
                output_data.append(line)
            output_data.append("pair_modify shift yes\n")
            continue
        elif re.match("pair_modify", line):
            continue
        elif re.match("kspace_style", line):
            continue
        elif re.match("thermo_style", line):
            output_data.append("thermo_style custom step ecoul evdwl ebond eangle edihed f_drycmap eimp\n")
            #output_data.append("thermo_style custom ebond eangle edihed eimp pe\n")
            continue
        elif re.match("run", line):
            continue
        output_data.append(line)
    # write some more stuff for the simulation (can be modified afterwards):
    output_data.append("neighbor 4.0 bin\n") # sets 4.0 extra angstroms cutoff for the list of neighbors (more efficient).
    output_data.append("neigh_modify every 1 delay 1 check yes\n\n")
    output_data.append("timestep 1.0\n") # time step in fs
    output_data.append("thermo 50\n\n") # output thermo every x steps
    output_data.append("minimize 1.0e-4 1.0e-6 10000 100000\n\n") # Energy minimization Etol Ftol maxiter maxeval
    #output_data.append("fix 1 all nvt temp 298.0 298.0 $(100.0*dt)\n")
    output_data.append("fix 1 all nve\n")
    #output_data.append("fix 2 all temp/csvr 298.1 298.1 $(100.0*dt) {}\n\n".format(np.random.randint(1, 10000000)))
    output_data.append("fix 2 all temp/csvr 298.1 298.1 $(100.0*dt) {}\n\n".format(np.random.randint(1, 10000000)))
    output_data.append("fix pl all plumed plumedfile {} outfile p.log\n\n".format(plumedfile))
    output_data.append("dump 1 all xtc 10000 {}\n".format(outxtc)) # save every 10 ps
    output_data.append("run    2000000000\n") # 2 us
    return output_data


def read_data_lammps(filename, pdb_file, cmap_file, cube_side_size=400.0):
    two_pow1over6 = 2**(1/6)
    output_data = []
    with open(filename, "r") as f:
        data = [el for el in f.readlines()]
    # read pdb for cmap topology in data file:
    pdb_cmaps_dict = read_pdb(pdb_file)
    Ncrossterms = len(list(pdb_cmaps_dict.keys()))
    for i, line in enumerate(data):
        if re.search("dihedrals", line):
            output_data.append(line)
            output_data.append("{} crossterms\n".format(Ncrossterms))
            continue 
        cube_half_side = cube_side_size/2
        if re.search("xlo xhi", line):
            output_data.append("{:.7f}  {:.7f} xlo xhi\n".format(-cube_half_side, cube_half_side))
            continue
        elif re.search("ylo yhi", line):
            output_data.append("{:.7f}  {:.7f} ylo yhi\n".format(-cube_half_side, cube_half_side))
            continue
        elif re.search("zlo zhi", line):
            output_data.append("{:.7f}  {:.7f} zlo zhi\n".format(-cube_half_side, cube_half_side))
            continue
        output_data.append(line)
    # CMAP addition:
    output_data.append("\nCMAP\n\n")
    # read cmap file for getting the type associated with each situation:
    with open(cmap_file, "r") as f:
        data = [el.rstrip("\n") for el in f.readlines() if re.search("type", el)]
    cmap_types = dict()
    for el in data:
        pair, cmtype = el.split()[1].rstrip(","), int(el.split()[3])
        pair = (pair.split("-")[0], pair.split("-")[1].rstrip(")").lstrip("("))
        cmap_types[pair] = cmtype
    # loop over the cmaps to write in data file and get pairs and type.
    for i, k in enumerate(pdb_cmaps_dict.keys()):
        # get restype and restype+1
        res_pair = (pdb_cmaps_dict[k][0][1], pdb_cmaps_dict[k][0][2]) if pdb_cmaps_dict[k][0][2] == "PRO" else (pdb_cmaps_dict[k][0][1], "XXX")
        # add info about cmap in lammps data top file
        output_data.append("{} {} {} {} {} {} {}\n".format(i+1, cmap_types[res_pair], pdb_cmaps_dict[k][1], pdb_cmaps_dict[k][2], pdb_cmaps_dict[k][3], pdb_cmaps_dict[k][4], pdb_cmaps_dict[k][5]))
    return output_data


if __name__ == "__main__":
    args = build_parser().parse_args()
    converted_input = read_input_lammps(args.i, args.otop, args.cmap, args.oxtc, args.plumed)
    with open(args.oin, "w") as f:
        for el in converted_input:
            f.write(el)
    # data file and CMAP editing:
    converted_data = read_data_lammps(args.top, args.pdb, args.cmap, args.cube)
    with open(args.otop, "w") as f:
        for el in converted_data:
            f.write(el)
         
