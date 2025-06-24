#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import random as rd
from argparse import ArgumentParser
from argparse import Namespace
from scipy.special import softmax, expit
from scipy.optimize import curve_fit
from lmfit import Minimizer, Parameters
from pickle import load
import re
import sys
import warnings


def build_parser():
    parser = ArgumentParser(description="////////   bAIes - PLUMED Preprocessing script   ////////")
    parser.add_argument('-pdb', type=str, help="Input alphafold PDB or CIF file")
    parser.add_argument('-mdpdb', type=str, help="MD simulation PDB or CIF file if the atom indexes were modified during the simulation set up (compared to the alphafold file).", default='same')
    parser.add_argument('-pkl', type=str, help="Pickle file containing the distograms. If using Colabfold, you can provide the .npy file that ends with '_prob_distributions.npy' in the 'distmat' folder here")
    parser.add_argument('-out', type=str, help="Output file to use as input for bAIes in PLUMED", default="baies_params.dat")
    parser.add_argument('-cutoff', help="cutoff distance in Angstroms (default 8 Angstroms). Use 'matrix' for a residue pair specific cutoff", default=8.0)
    parser.add_argument('-model', type=str, help="distogram fitting model (gauss or lognorm)", default="gauss")
    parser.add_argument('-seqsep', type=int, help="Minimum separation in the sequence. Removes the neighboring residues. default=3 (best to keep helix information.)", default=3)
    parser.add_argument('-chains', type=str, help="For multichain systems, consider only intramolecular: 'intra' or intermolecular: 'inter'. Default is 'all'", default="all")
    parser.add_argument('-plotout', type=str, help="Output file for the plots if --plots is set.", default="distograms.pdf") ###############################################################################
    parser.add_argument('-ndxout', type=str, help="Output index file containing the list of atoms for PLUMED.", default="atom_list.ndx")
    parser.add_argument('--plots', help="Put this flag to plot the fitting results during preprocessing for quality check", action='store_true') #############################################3
    parser.add_argument('--verbose', help="Put this flag to print more information on the terminal during the proprocessing (Recommended)", action='store_true')
    return parser


# The distance cutoff matrix for residue-residue interactions from Baker lab
cutoff_matrix = {
("GLY", "GLY"): (4.467,0.017),
("GLY", "ALA"): (5.201,0.269),
("GLY", "SER"): (5.51,0.153),
("GLY", "VAL"): (5.671,0.107),
("GLY", "CYS"): (5.777,0.129),
("GLY", "THR"): (5.619,0.12),
("GLY", "PRO"): (6.14,0.245),
("GLY", "ASP"): (6.135,0.193),
("GLY", "ASN"): (6.321,0.169),
("GLY", "ILE"): (6.413,0.179),
("GLY", "LEU"): (6.554,0.125),
("GLY", "GLU"): (7.036,0.249),
("GLY", "GLN"): (7.297,0.216),
("GLY", "MET"): (7.383,0.255),
("GLY", "HIS"): (7.472,0.206),
("GLY", "LYS"): (8.216,0.358),
("GLY", "PHE"): (7.966,0.219),
("GLY", "TYR"): (9.098,0.267),
("GLY", "ARG"): (9.166,0.334),
("GLY", "TRP"): (8.966,0.239),
("ALA", "ALA"): (5.381,0.262),
("ALA", "SER"): (5.829,0.291),
("ALA", "VAL"): (5.854,0.312),
("ALA", "CYS"): (6.057,0.394),
("ALA", "THR"): (5.982,0.378),
("ALA", "PRO"): (6.412,0.399),
("ALA", "ASP"): (6.388,0.289),
("ALA", "ASN"): (6.766,0.349),
("ALA", "ILE"): (6.587,0.214),
("ALA", "LEU"): (6.707,0.25),
("ALA", "GLU"): (7.124,0.34),
("ALA", "GLN"): (7.583,0.356),
("ALA", "MET"): (7.605,0.394),
("ALA", "HIS"): (7.591,0.38),
("ALA", "LYS"): (8.327,0.55),
("ALA", "PHE"): (8.162,0.26),
("ALA", "TYR"): (9.121,0.443),
("ALA", "ARG"): (9.365,0.485),
("ALA", "TRP"): (9.252,0.29),
("SER", "SER"): (6.19,0.292),
("SER", "VAL"): (6.567,0.205),
("SER", "CYS"): (6.59,0.24),
("SER", "THR"): (6.45,0.214),
("SER", "PRO"): (6.937,0.321),
("SER", "ASP"): (6.76,0.323),
("SER", "ASN"): (7.081,0.305),
("SER", "ILE"): (7.142,0.342),
("SER", "LEU"): (7.394,0.287),
("SER", "GLU"): (7.483,0.446),
("SER", "GLN"): (7.807,0.408),
("SER", "MET"): (8.01,0.369),
("SER", "HIS"): (8.051,0.435),
("SER", "LYS"): (8.792,0.445),
("SER", "PHE"): (8.694,0.394),
("SER", "TYR"): (9.594,0.467),
("SER", "ARG"): (9.753,0.483),
("SER", "TRP"): (9.77,0.497),
("VAL", "VAL"): (6.759,0.145),
("VAL", "CYS"): (6.941,0.173),
("VAL", "THR"): (6.791,0.138),
("VAL", "PRO"): (7.063,0.298),
("VAL", "ASP"): (6.972,0.287),
("VAL", "ASN"): (7.219,0.232),
("VAL", "ILE"): (7.441,0.242),
("VAL", "LEU"): (7.633,0.179),
("VAL", "GLU"): (7.404,0.51),
("VAL", "GLN"): (8.008,0.359),
("VAL", "MET"): (8.335,0.295),
("VAL", "HIS"): (8.179,0.383),
("VAL", "LYS"): (8.077,0.634),
("VAL", "PHE"): (9.057,0.246),
("VAL", "TYR"): (9.442,0.535),
("VAL", "ARG"): (9.513,0.514),
("VAL", "TRP"): (10.021,0.271),
("CYS", "CYS"): (6.426,0.178),
("CYS", "THR"): (6.801,0.181),
("CYS", "PRO"): (7.157,0.259),
("CYS", "ASP"): (6.985,0.299),
("CYS", "ASN"): (7.205,0.24),
("CYS", "ILE"): (7.476,0.295),
("CYS", "LEU"): (7.685,0.206),
("CYS", "GLU"): (7.449,0.538),
("CYS", "GLN"): (7.962,0.347),
("CYS", "MET"): (8.265,0.439),
("CYS", "HIS"): (8.422,0.203),
("CYS", "LYS"): (8.494,0.521),
("CYS", "PHE"): (9.026,0.286),
("CYS", "TYR"): (9.362,0.585),
("CYS", "ARG"): (9.46,0.491),
("CYS", "TRP"): (9.752,0.417),
("THR", "THR"): (6.676,0.188),
("THR", "PRO"): (7.062,0.32),
("THR", "ASP"): (6.971,0.307),
("THR", "ASN"): (7.159,0.262),
("THR", "ILE"): (7.442,0.259),
("THR", "LEU"): (7.642,0.19),
("THR", "GLU"): (7.628,0.409),
("THR", "GLN"): (8.055,0.378),
("THR", "MET"): (8.397,0.292),
("THR", "HIS"): (8.221,0.417),
("THR", "LYS"): (8.715,0.464),
("THR", "PHE"): (9.03,0.264),
("THR", "TYR"): (9.813,0.43),
("THR", "ARG"): (9.764,0.477),
("THR", "TRP"): (9.98,0.315),
("PRO", "PRO"): (7.288,0.339),
("PRO", "ASP"): (7.321,0.416),
("PRO", "ASN"): (7.497,0.334),
("PRO", "ILE"): (7.554,0.336),
("PRO", "LEU"): (7.751,0.317),
("PRO", "GLU"): (7.938,0.475),
("PRO", "GLN"): (8.308,0.41),
("PRO", "MET"): (8.247,0.388),
("PRO", "HIS"): (8.537,0.457),
("PRO", "LYS"): (9.198,0.55),
("PRO", "PHE"): (8.895,0.425),
("PRO", "TYR"): (9.965,0.506),
("PRO", "ARG"): (10.266,0.506),
("PRO", "TRP"): (9.719,0.462),
("ASP", "ASP"): (8.001,0.392),
("ASP", "ASN"): (7.672,0.337),
("ASP", "ILE"): (7.472,0.341),
("ASP", "LEU"): (7.696,0.348),
("ASP", "GLU"): (8.945,0.354),
("ASP", "GLN"): (8.601,0.357),
("ASP", "MET"): (8.401,0.361),
("ASP", "HIS"): (8.634,0.325),
("ASP", "LYS"): (9.306,0.343),
("ASP", "PHE"): (9.111,0.351),
("ASP", "TYR"): (9.979,0.676),
("ASP", "ARG"): (10.123,0.327),
("ASP", "TRP"): (9.867,0.475),
("ASN", "ASN"): (7.682,0.249),
("ASN", "ILE"): (7.631,0.341),
("ASN", "LEU"): (7.889,0.279),
("ASN", "GLU"): (8.485,0.423),
("ASN", "GLN"): (8.502,0.373),
("ASN", "MET"): (8.55,0.31),
("ASN", "HIS"): (8.672,0.289),
("ASN", "LYS"): (9.319,0.398),
("ASN", "PHE"): (9.168,0.393),
("ASN", "TYR"): (10.039,0.586),
("ASN", "ARG"): (10.135,0.372),
("ASN", "TRP"): (9.976,0.458),
("ILE", "ILE"): (8.096,0.321),
("ILE", "LEU"): (8.342,0.261),
("ILE", "GLU"): (7.949,0.453),
("ILE", "GLN"): (8.302,0.406),
("ILE", "MET"): (8.874,0.327),
("ILE", "HIS"): (8.523,0.379),
("ILE", "LYS"): (8.329,0.582),
("ILE", "PHE"): (9.602,0.347),
("ILE", "TYR"): (9.719,0.589),
("ILE", "ARG"): (9.746,0.557),
("ILE", "TRP"): (10.47,0.397),
("LEU", "LEU"): (8.522,0.198),
("LEU", "GLU"): (8.077,0.475),
("LEU", "GLN"): (8.48,0.411),
("LEU", "MET"): (9.122,0.318),
("LEU", "HIS"): (8.676,0.401),
("LEU", "LYS"): (8.479,0.591),
("LEU", "PHE"): (9.9,0.26),
("LEU", "TYR"): (9.889,0.611),
("LEU", "ARG"): (9.852,0.578),
("LEU", "TRP"): (10.707,0.331),
("GLU", "GLU"): (9.863,0.389),
("GLU", "GLN"): (9.328,0.45),
("GLU", "MET"): (8.87,0.511),
("GLU", "HIS"): (9.454,0.443),
("GLU", "LYS"): (9.842,0.434),
("GLU", "PHE"): (9.403,0.512),
("GLU", "TYR"): (10.544,0.469),
("GLU", "ARG"): (10.713,0.363),
("GLU", "TRP"): (10.303,0.493),
("GLN", "GLN"): (9.074,0.436),
("GLN", "MET"): (9.102,0.498),
("GLN", "HIS"): (9.391,0.401),
("GLN", "LYS"): (9.667,0.521),
("GLN", "PHE"): (9.506,0.451),
("GLN", "TYR"): (10.534,0.547),
("GLN", "ARG"): (10.61,0.535),
("GLN", "TRP"): (10.429,0.49),
("MET", "MET"): (9.53,0.457),
("MET", "HIS"): (9.396,0.342),
("MET", "LYS"): (9.096,0.611),
("MET", "PHE"): (10.253,0.377),
("MET", "TYR"): (10.4,0.661),
("MET", "ARG"): (10.25,0.641),
("MET", "TRP"): (11.11,0.397),
("HIS", "HIS"): (10.606,0.333),
("HIS", "LYS"): (9.582,0.714),
("HIS", "PHE"): (9.602,0.542),
("HIS", "TYR"): (10.843,0.554),
("HIS", "ARG"): (10.879,0.595),
("HIS", "TRP"): (10.661,0.458),
("LYS", "LYS"): (10.662,0.738),
("LYS", "PHE"): (9.344,0.441),
("LYS", "TYR"): (10.627,0.704),
("LYS", "ARG"): (11.322,0.648),
("LYS", "TRP"): (10.136,0.47),
("PHE", "PHE"): (10.903,0.46),
("PHE", "TYR"): (10.999,0.767),
("PHE", "ARG"): (10.577,0.738),
("PHE", "TRP"): (11.758,0.447),
("TYR", "TYR"): (11.536,0.855),
("TYR", "ARG"): (11.615,0.822),
("TYR", "TRP"): (11.807,0.684),
("ARG", "ARG"): (12.05,0.704),
("ARG", "TRP"): (11.355,0.889),
("TRP", "TRP"): (12.806,0.473)}


def read_mmCIF(file_name):
    '''
    Reads an mmCIF file and returns the associated dictionnary.
    '''
    output_dic = {}
    with open(file_name, "r") as f:
        file_data = [el.rstrip("\n") for el in f.readlines()]
    data_block = "None"
    in_loop = False
    loop_status = "None"
    # in case we don't need to read a given line:
    jump_lines = 0
    for l, line in enumerate(file_data):
        # if asked to jump a line, jump a line
        if jump_lines > 0:
            jump_lines -= 1
            continue
        # read data blocks:
        if re.match('data_', line):
            data_block = line.lstrip('data_')
            output_dic[data_block] = {}
        # read loops:
        elif re.match('loop_', line):
            in_loop = True
            loop_status = "Header"
            loop_header = None
        # read the content of a loop
        elif in_loop:
            # if a # comment is read, the loop is ended:
            if loop_status == "Header" and re.match("_", line):
                [key1, key2] = line.lstrip("_").rstrip(" ").split(".")
                if key1 not in output_dic[data_block].keys():
                    output_dic[data_block][key1] = {}
                    loop_header = key1
                output_dic[data_block][key1][key2] = []
            elif re.match("#", line) and in_loop:
                in_loop = False
                loop_status = "None"
                loop_header = None
            elif loop_status == "Rows" or (loop_status == "Header" and not re.match("_", line)):
                loop_status = "Rows"
                keys = list(output_dic[data_block][loop_header].keys())
                line = line.split()
                for i, k in enumerate(keys):
                    if i >= len(line):
                        value = file_data[l+1]
                        jump_lines = 1
                        if re.match(";", value):
                            value = value.lstrip(";")
                            between_semicolons = True
                            idx = 1
                            while between_semicolons:
                                if not re.match(";", file_data[l+1+idx]):
                                    value = value + file_data[l+1+idx]
                                    idx += 1
                                else:
                                    between_semicolons = False
                                jump_lines += 1
                        output_dic[data_block][loop_header][k].append(value)
                    else:
                        output_dic[data_block][loop_header][k].append(line[i])
        # read _ key entries outside a loop
        elif re.match("_", line) and loop_status == "None":
            line = line.split()
            [key1, key2] = line[0].lstrip("_").split(".")
            if key1 not in output_dic[data_block].keys():
                output_dic[data_block][key1] = {}
            if len(line) == 1:
                value = file_data[l+1]
                output_dic[data_block][key1][key2] = value
                jump_lines = 1
            else:
                output_dic[data_block][key1][key2] = line[1]
    return output_dic


def read_cif(input_file):
    '''
    Takes an mmCIF file, reads it and extracts the important information in a dictionnary as:
    {distogram residue index: (atom nb, residue type, chain label, residue nb in the pdb)}
    '''
    cif_dic = read_mmCIF(input_file)
    # for now we consier only the first data block. If there is several the latter data blocks will be ignored.
    block = list(cif_dic.keys())[0]
    # get the atom indices:
    atom_idx = [int(el) for i, el in enumerate(cif_dic[block]["atom_site"]["id"]) if cif_dic[block]["atom_site"]["group_PDB"][i] == "ATOM"]
    # Atom types in the mmCIF file:
    atom_type = [el for i, el in enumerate(cif_dic[block]["atom_site"]["label_atom_id"]) if cif_dic[block]["atom_site"]["group_PDB"][i] == "ATOM"]
    # Residue type in the mmCIF file:
    res_type = [el for i, el in enumerate(cif_dic[block]["atom_site"]["label_comp_id"]) if cif_dic[block]["atom_site"]["group_PDB"][i] == "ATOM"]
    # Residue indexes in the mmCIF file:
    res_idx = [int(el) for i, el in enumerate(cif_dic[block]["atom_site"]["label_seq_id"]) if cif_dic[block]["atom_site"]["group_PDB"][i] == "ATOM"]
    # Chain index in the mmCIF file:
    chain_idx = [el for i, el in enumerate(cif_dic[block]["atom_site"]["label_asym_id"]) if cif_dic[block]["atom_site"]["group_PDB"][i] == "ATOM"]
    res_data = {(k, chain_idx[i]): (atom_idx[i], res_type[i], chain_idx[i]) for i, k in enumerate(res_idx) if atom_type[i] == "CB" or (atom_type[i] == "CA" and res_type[i] == "GLY")}
    dist_res_idx = [i+1 for i in range(len([el[0] for el in res_data.keys()]))]
    res_data = {dist_res_idx[i]: (res_data[k][0], res_data[k][1], k[1], k[0]) for i, k in enumerate(res_data.keys())}
    print(res_data)
    return res_data


def read_pdb(input="protein.pdb"):
    """
    Reads a pdb file and outputs a dictionnary as:
    {distogram residue index: (atom nb, residue type, chain label, residue nb in the pdb)}
    """
    # Extract the data:
    with open(input, "r") as f:
        data = [el.rstrip("\n") for el in f.readlines() if el.split()[0] == 'ATOM']
    # Atom numbers in the PDB file:
    atom_idx = [(int(el[5:11])) for el in data]
    # Atom types in the PDB file:
    atom_type = [el[12:17].rstrip(" ").lstrip(" ") for el in data]
    # Residue type in the PBD file:
    res_type = [el[17:20] for el in data]
    # original residue indexes
    res_idx = [int(el[22:26]) for el in data]
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

    res_data = {(k, chain_idx[i]): (atom_idx[i], res_type[i], chain_idx[i]) for i, k in enumerate(res_idx) if atom_type[i] == "CB" or (atom_type[i] == "CA" and res_type[i] == "GLY")}
    dist_res_idx = [i+1 for i in range(len([el[0] for el in res_data.keys()]))]
    res_data = {dist_res_idx[i]: (res_data[k][0], res_data[k][1], k[1], k[0]) for i, k in enumerate(res_data.keys())}
    # output: {residue idx associated with the distogram[idx][j] or [i][idx]: atom idx, residue type, chain label, residue idx (from the pdb)}
    return res_data


def read_input(filename):
    '''
    Reads the pdb of mmCIF input file and returns the important information as:
    {distogram residue index: (atom nb, residue type, chain label, residue nb in the pdb)}
    '''
    if filename.split(".")[-1] == "pdb":
        return read_pdb(filename)
    elif filename.split(".")[-1] == "pdb":
        return read_cif(filename)
    else:
        print("Input file format not understood. Please enter a .pdb or a .cif file.")
        sys.exit(1)


def read_pkl(pkl_file, convf="softmax"):
    """
    Reads the distograms pkl file outputed from Alphafold and returns the bins (63) and distograms (NxNx64) in the format adapted for the following preprocessing.
    """
    with open(pkl_file, "rb") as f:
        data = load(f)['distogram']
    # Bins -> distance axis: initially 64 bins from 2 to 22 but the last one covers also above 22.
    # 20/64 = 0.3125 = delta_bins. first bin (2.3125) covers 2 to 2.3125. 
    # center of it is therefore 2.3125 - delta_bins/2
    bins = data['bin_edges']
    delta_bins = bins[1]-bins[0]
    # Shift the bins such that the distribution data point are assumed in the center of the bin:
    # (And converting the bin values from Angstrom to nm)
    bins = [(el-(delta_bins/2))*0.1 for el in bins]
    distograms = data['logits']
    # Remove last bin since it is the above 22 Angstrom garbage ([:-1])
    # And apply conversion to probability distribution function to obtain the distribution:
    # distograms = np.array([[softmax(ell[:]) for ell in el] for el in distograms])
    if convf == "softmax":
        distograms = np.array([[softmax(ell[:]) for ell in el] for el in distograms])
    elif convf == "sigmoid":
        distograms = np.array([[expit(ell[:]) for ell in el] for el in distograms])
        distograms = np.array([[ell[:]/np.sum(ell) for ell in el] for el in distograms])
    else:
        print("conversion function not recognized. Please select softmax (default) or sigmoid.")
        sys.exit(1)
    return bins, distograms


def read_npy(pkl_file):
    """
    Reads the distograms npy file outputed from Colabfold and returns the bins (63) and distograms (NxNx64) in the format adapted for the following preprocessing.
    """
    data = np.load(pkl_file) # distogram matrix.
    bins = [0.215625, 0.246875, 0.278125, 0.309375, 0.340625, 0.371875, 0.403125, 0.434375, 0.465625, 0.496875, 0.5281250000000001, 0.5593750000000001, 0.5906250000000001, 0.6218750000000001, 0.6531250000000001, 0.6843750000000001,\
            0.7156250000000001, 0.7468750000000001, 0.7781250000000001, 0.8093750000000001, 0.8406250000000001, 0.8718750000000001, 0.9031250000000001, 0.9343750000000001, 0.9656249046325684, 0.9968750000000001, 1.028125, 1.059375,\
            1.090625, 1.121875, 1.153125, 1.184375, 1.2156250000000002, 1.2468750000000002, 1.2781250000000002, 1.3093750000000002, 1.3406250000000002, 1.3718750000000002, 1.4031250000000002, 1.4343749046325684, 1.4656250000000002,\
            1.4968750000000002, 1.5281250000000002, 1.5593750000000002, 1.5906250000000002, 1.6218750000000002, 1.6531250000000002, 1.6843748092651367, 1.7156248092651367, 1.7468750000000002, 1.7781250000000002, 1.8093750000000002,\
            1.8406250000000002, 1.8718750000000002, 1.9031250000000002, 1.9343750000000002, 1.9656250000000002, 1.9968750000000002, 2.0281248092651367, 2.059375, 2.090625, 2.121875, 2.153125]
    return bins, data


def dmax(bins, distogram):
    """
    Obtain highest probability distance from distogram
    """
    # add last bin corresponding to last bin + >22 Angstroms garbage bin.
    # make a copy of the bins and add 2.2 (nm) corresponding to the garbage:
    dmax_bins = bins[:]
    dmax_bins.append(2.2)
    # Get dmax:
    Ymax = np.max(distogram)
    Dmax = dmax_bins[list(distogram).index(Ymax)]
    return Dmax


### MODEL FUNCTIONS


def gauss_model(x, mu, sigma2, scaling=1, sigma_is_squared=True):
    # Gaussian function with scaling factor
    if sigma_is_squared:
        return scaling*(1/(np.sqrt(2*np.pi*sigma2)))*np.exp(-((x-mu)**2)/(2*sigma2))
    return scaling*(1/(sigma2*np.sqrt(2*np.pi)))*np.exp(-((x-mu)**2)/(2*sigma2**2))


def lognormal_model(x, mu, sigma2, scaling=1, sigma_is_squared=True):
    # Lognormal function with scaling factor (by default we take sigma squared to have the absolute value in the end (otherwise can lead to negative values))
    if sigma_is_squared:
        return scaling*(1/(1*np.sqrt(2*np.pi*sigma2)))*np.exp(-((np.log(x)-mu)**2)/(2*sigma2))
    return scaling*(1/(1*sigma2*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-mu)**2)/(2*sigma2**2))


def gausss_model(x, mus, sigmas2, scalings, sigmas_are_squared=True):
    # Multigaussian function with scaling factors (sigma squared fitted by default)
    if sigmas_are_squared:
        return np.sum([scalings[i]*(1/(np.sqrt(2*np.pi*sigmas2[i])))*np.exp(-((x-mus[i])**2)/(2*sigmas2[i])) for i, el in enumerate(mus)])
    return np.sum([scalings[i]*(1/(sigmas2[i]*np.sqrt(2*np.pi)))*np.exp(-((x-mus[i])**2)/(2*sigmas2[i]**2)) for i, el in enumerate(mus)])


def lognormals_model(x, mus, sigmas2, scalings, sigmas_are_squared=True):
    # Multilognormal function with scaling factors (by default we fit sigma squared to have the absolute value in the end (otherwise can lead to negative values))
    if sigmas_are_squared:
        return np.sum([scalings[i]*(1/(1*np.sqrt(2*np.pi*sigmas2[i])))*np.exp(-((np.log(x)-mus[i])**2)/(2*sigmas2[i])) for i, el in enumerate(mus)])
    return np.sum([scalings[i]*(1/(1*sigmas2[i]*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-mus[i])**2)/(2*sigmas2[i]**2)) for i, el in enumerate(mus)])


### INNER FITTING FUNCTIONS ###

## MULTIMODAL FIT FNCS ##

def make_params(modes=1, model_func=gausss_model, dmax=None):
    """
    Parameter making function for lmfit
    """
    # generate parameters for the multimodal fit taking random initial conditions for the fit.
    params = Parameters()
    for i in range(modes):
        if model_func == gausss_model:
            # if a dmax is inputed, the dmax value is set as initial condition for the first gaussian to ease the fit.
            # random values are used for the rest, with appropriate bounds
            mu_value = dmax if (dmax is not None and i == 0) else rd.uniform(0.2, 2.2)
            params.add('mu_'+str(i+1), value=mu_value, min=0, max=3.0)
            params.add('sigma2_'+str(i+1), value=rd.uniform(0.0009, 0.2), min=0, max=1)
            params.add('scaling_'+str(i+1), value=rd.uniform(0.01, 0.07), min=0)
        elif model_func == lognormals_model:
            # if a dmax is inputed, the dmax value is used in the initial condition for the first lognorm to ease the fit.
            # random values are used for the rest, with appropriate bounds
            sigma2_value = rd.uniform(0.0009, 0.1)
            mu_value = np.log(dmax)+sigma2_value if dmax and i == 0 else rd.uniform(-1.4, 0.9)
            params.add('mu_'+str(i+1), value=mu_value, min=-2, max=1)
            params.add('sigma2_'+str(i+1), value=sigma2_value, min=0, max=1)
            params.add('scaling_'+str(i+1), value=rd.uniform(0.01, 0.07), min=0)
        else:
            print("Fitting model invalid. Please select gauss or lognorm")
            sys.exit(1)
    return params


def minfnc(params, Xbins, model_fnc, distogram, tolerr=0.005):
    """
    Minimizing function for lmfit
    """
    # minimizing function for multimodal fit
    # get the parameters
    mus = [el.value for el in list(params.values()) if el.name[:2] == 'mu']
    sigmas2 = [el.value for el in list(params.values()) if el.name[:3] == 'sig']
    scalings = [el.value for el in list(params.values()) if el.name[:3] == 'sca']
    # Calculate model function with given parameters
    model = np.array([model_fnc(x, mus, sigmas2, scalings, sigmas_are_squared=True) for x in Xbins])
    # get data
    data = np.array(distogram)
    return (data - model) / np.array([tolerr for el in data])


def best_fit(results_dic):
    """
    Get best fit from a dictionnary of lmfit fitting results. Based on reduced chi2.
    Returns the best lmfit result
    """
    best_chi2 = None
    output = {}
    for k, el in results_dic.items():
        if not best_chi2:
            best_chi2 = results_dic[k].redchi
            output = results_dic[k]
        else:
            if results_dic[k].redchi < best_chi2:
                best_chi2 = results_dic[k].redchi
                output = results_dic[k]
    return output


def fit_model(Xbins, distogram, modes=1, model_func=gausss_model, func2min=minfnc, method='leastsq', Niterations=5, dmax=None, tolerr=0.005):
    """
    Fits a distogram to a single specified model.
    Returns the best result (best chi2) for the given model.
    """
    # Fitting function for a multimode model; Several iterations to sample initial conditions and improve the overall robustness of the results.
    results = {}
    # For each iteration, try to fit a model, and register the fit in the result dictionary. At the end, the best fit is kept. 
    # (Sometimes the initial conditions are very important so we try different conditions and keep the one that gives the best fit.)
    for i in range(Niterations):
        try:
            params = make_params(modes=modes, model_func=model_func, dmax=dmax)
            minner = Minimizer(func2min, params, fcn_args=(Xbins, model_func, distogram, tolerr))
            result = minner.minimize(method=method, maxfev=5000)
            results[i] = result
        # if the fit didn't work, we try again several times with one center on Dmax in the initial conditions:
        except ValueError:
            for j in range(10):
                try:
                    params = make_params(modes=modes, model_func=model_func, dmax=dmax)
                    minner = Minimizer(func2min, params, fcn_args=(Xbins, model_func, distogram, tolerr))
                    result = minner.minimize(method=method, maxfev=5000)
                    results[i] = result
                    break
                except ValueError:
                    pass
    return best_fit(results)


def get_parameters_from_result(result):
    """
    Takes an lmfit fitting result object and returns our parameters of interest as:
    [mus], [sigmas2], [scalings], statistics{"stat": value}
    """
    # Extract the parameters from an lmfit fitting result.
    values = result.params.valuesdict()
    mus, sigmas2, scalings = [], [], []
    for k in values.keys():
        if re.match("mu*", k):
            mus.append((int(re.search("_\d+", k).group().lstrip("_")), values[k]))
        elif re.match("sigmas*", k):
            sigmas2.append((int(re.search("_\d+", k).group().lstrip("_")), values[k]))
        elif re.match("scalings*", k):
            scalings.append((int(re.search("_\d+", k).group().lstrip("_")), values[k]))
    mus = [el[1] for el in sorted(mus)]
    sigmas2 = [el[1] for el in sorted(sigmas2)]
    scalings = [el[1] for el in sorted(scalings)]
    stats = {"redchi": result.redchi, "aic": result.aic, "bic": result.bic}
    return mus, sigmas2, scalings, stats


def fit_models(bins, distogram, model="gauss", plots=False, plot_labels=((0,0), (0,0)), nmodes=1, dmax=None, tolerr=0.005, verbose=False):
    """
    Performs the multimodal fittings and the plottings (if asked to).
    Takes the distogram of interest as input and returns the fitted parameters as a tuple (nb of modes, [mu], [sigma2], [scalings]).
    """
    # model functions:
    model_funcs = {"gauss": gausss_model, "lognorm": lognormals_model}
    colors = {1: "darkorange", 2: "red", 3: "purple", 4: "black"}
    if model not in list(model_funcs.keys()):
        print("Model not recognized. Please select between gauss and lognorm")
        sys.exit(0)
    # fitting part stored in a dictionnary {nb of modes (k): result}:
    fitting_results = {}
    for i in range(nmodes):
        if verbose: print("fitting ",i+1," modes")
        fit_result = fit_model(bins, distogram, modes=i+1, model_func=model_funcs[model], dmax=dmax, tolerr=tolerr)
        # if fit_result is not None (not failed), add in the dictionnary
        if fit_result != {}:
            fitting_results[i+1] = fit_result
        else:
            print("Fitting ", i+1, " modes failed.")
    ## Plotting part
    # plot AF distogram if plot asked:
    if plots:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        ax.plot(bins, distogram, color="darkblue", marker='o', label="AF distogram")
        ax.set_xlabel("Pairwise distance (nm)")
    # If the fittings failed, we close the plot and pass this residue pair (Normally should never happen.)
    if fitting_results == {}:
        if plots:
            ax.scatter([], [], label="Fitting failed.", color="white")
            ax.legend(frameon=False)
            ax.set_title("residues: ("+str(plot_labels[0][0])+","+str(plot_labels[0][1])+r"), atoms: ("+str(plot_labels[1][0])+","+str(plot_labels[1][1])+")")
            plots.savefig()
        return None
    # Sort results by the reduced chi2:
    results_by_redchi = []
    # For each fitted model do:
    for k, result in fitting_results.items():
        # get parameters
        mus, sigmas2, scalings, stats = get_parameters_from_result(result)
        # store the model results based on the reduced chi2
        results_by_redchi.append((abs(stats["redchi"]-1), k, mus, sigmas2, scalings))
        # Plot model for each model fitted:
        if plots:
            Ymodel = [model_funcs[model](x, mus, sigmas2, scalings, sigmas_are_squared=True) for x in bins]
            ax.plot(bins, Ymodel, color=colors[k], label=str(k)+" modes: "+r"red$\chi^{2}$="+str(stats["redchi"].round(3))+"  AIC="+str(stats["aic"].round(3))+"  BIC="+str(stats["bic"].round(3)))
    # The models results are sorted based on the reduced chi2 closest to 1
    # And the most suited model (closest to 1) is selected:
    best_model_result = sorted(results_by_redchi)[0]
    # output format: (nb modes, mus, sigmas2, scalings)
    output = (best_model_result[1], best_model_result[2], best_model_result[3], best_model_result[4])
    # Conclude and close the plot file if plot asked:
    if plots:
        ax.scatter([], [], label=str(output[0])+" modes selected", color="white")
        ax.legend(frameon=False)
        ax.set_title("residues: ("+str(plot_labels[0][0])+","+str(plot_labels[0][1])+r"), atoms: ("+str(plot_labels[1][0])+","+str(plot_labels[1][1])+")")
        plots.savefig()
    return output


## SINGLE MODE FIT FUNCS ##


def iterate_fits(model_fnc, bins, distogram, bounds, max_try=100, dmax=0.75):
    """
    Performs fitting of the distogram in an iterative way until a solution is found.
    """
    # fitting more difficult fits
    idx = 0
    popt_out, pcov_out = [], []
    # try fitting with random initial conditions
    for i in range(max_try):
        try:
            if model_fnc == gauss_model:
                #p0 = [rd.uniform(2, 22), rd.uniform(0.1, 5), rd.uniform(1.5, 2.5)]
                p0 = [dmax, rd.uniform(0.0009, 0.2), rd.uniform(0.01, 0.06)]
            elif model_fnc == lognormal_model:
                #p0 = [rd.uniform(0.9, 3.1), rd.uniform(0.001, 0.1), rd.uniform(0.1, 2.0)]
                sigma2 = rd.uniform(0.0009, 0.09)
                p0 = [np.log(dmax)+sigma2, sigma2, rd.uniform(0.01, 0.06)]
            popt, pcov = curve_fit(model_fnc, bins, distogram, p0=p0, maxfev=5000, bounds=bounds)
            popt_out = popt
            pcov_out = pcov
            break
        except RuntimeError:
            pass
    if popt_out == []:
        print("fitting failed due to systematic RuntimeError")
        return None, None
    return popt_out, pcov_out


def fit_single(bins, distogram, model="gauss", plots=False, plot_labels=((0,0), (0,0)), dmax=0.75):
    """
    Fits a distogram to a single distribution model (Gaussian or Lognormal).
    Takes the distogram of interest as input and returns the fitted parameters as a tuple (1, [mu], [sigma2], [1]).
    Also performs the plottings if asked.
    """
    if model == "gauss":
        model_fnc = gauss_model
        bounds = ((0, 0, 0), (3.0, np.inf, np.inf))
        p0 = [dmax, 0.07, 0.05]
    elif model == "lognorm":
        model_fnc = lognormal_model
        bounds = ((-2, 0, 0), (1, np.inf, np.inf))
        p0 = [np.log(dmax)+0.0081, 0.0081, 0.05]
    else:
        print("Model not recognized. Please select between gauss and lognorm")
        sys.exit(0)
    try:
        popt, pcov = curve_fit(model_fnc, bins, distogram, p0=p0, maxfev=5000, bounds=bounds)
    except RuntimeError:
        popt, pcov = iterate_fits(model_fnc, bins, distogram, bounds=bounds, dmax=dmax)
    # plot AF distogram if plot asked:
    if plots:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        ax.set_xlabel("Pairwise distance (nm)")
        ax.plot(bins, distogram, color="darkblue", marker='o', label="AF distogram")
    # If the fittings failed, we close the plot and pass this residue pair (Normally should not happen.)
    if popt is None:
        if plots:
            ax.scatter([], [], ls='', label="Fitting failed.", color="white")
            ax.legend(frameon=False)
            ax.set_title("residues: ("+str(plot_labels[0][0])+","+str(plot_labels[0][1])+r"), atoms: ("+str(plot_labels[1][0])+","+str(plot_labels[1][1])+")")
            plots.savefig()
        return None
    # plot if specified:
    if plots:
        Ymodel = [model_fnc(x, *popt, sigma_is_squared=True) for x in bins]
        ax.plot(bins, Ymodel, color="darkorange", label="fit")
        ax.legend(frameon=False)
        ax.set_title("residues: ("+str(plot_labels[0][0])+","+str(plot_labels[0][1])+r"), atoms: ("+str(plot_labels[1][0])+","+str(plot_labels[1][1])+")")
        plots.savefig()
    mu, sigma2 = popt[0], popt[1]
    # The output is formatted the same way as for the multimodal fit
    return (1, [mu], [sigma2], [1])


### END INNER FITTING FUNCTIONS ###


def fit_distogram(bins, distogram, model="gauss", plots=False, plot_labels=((0,0), (0,0)), nmodes=1, dmax=None, tolerr=0.005, verbose=False):
    """
    The global fit single distogram function. 
    Essentially points towards the right fitting procedure whether we need to perform a model selection ot not.
    """
    if nmodes == 1:
        result = fit_single(bins=bins, distogram=distogram, model=model, plots=plots, plot_labels=plot_labels, dmax=dmax)
    else:
        result = fit_models(bins, distogram, model=model, plots=plots, plot_labels=plot_labels, nmodes=nmodes, dmax=dmax, tolerr=tolerr, verbose=verbose)
    # output: (nmodes, mus, sigmas, scalings)
    return result


def select_and_fit_distograms(bins, distograms, res_data, res_data_md, cutoff, model='gauss', seqsep=3, plots=False, verbose=False, nmodes_max=1, tolerr=0.005, plotout="output.pdf", chains='all'):
    """
    Takes the distogram and bins information along witht the option parameters as input
    Select the relevant atom pairs and launch the distogram fitting for the selected distos.
    Returns the results as a dictionnary with the pairs as keys and the fitted parameters as values.
    """
    # If asked to plot the fitting results, initialization of the pdf containing the plots:
    if plots:
        import matplotlib.backends.backend_pdf
        plots = matplotlib.backends.backend_pdf.PdfPages(plotout)
    # Define a dctionnary where the fit for each selected residue pair will be stored
    results = {}
    # For each distogram do:
    for i in range(len(distograms)):
        for j in range(i, len(distograms[i])):
            ## Ignore the redundant, non existant and unwanted residue pairs:
            # if intramolecular restraints only, removes the intermolecular ones:
            if chains == "intra" and res_data[j+1][2] != res_data[i+1][2]:
                continue
            # if intermolecular restrains only, removes the intramolecular ones:
            elif chains == "inter" and res_data[j+1][2] == res_data[i+1][2]:
                continue
            # If the two residues are neighbors (or the same), and in the same chain, pass:
            if abs(res_data[j+1][3]-res_data[i+1][3]) <= seqsep and res_data[j+1][2] == res_data[i+1][2]:
                continue
            # optional additional verbose printing:
            if verbose: print("residue pair: "+str(i+1)+","+str(j+1))
            # Obtain Dmax:
            Dmax = dmax(bins, distograms[i][j])
            # Cutoff definition. If cutoff defined as 'matrix', do as below:
            if type(cutoff) == str and cutoff == "matrix":
                # residue pair-dependent cutoff distances
                if (res_data[i+1][1], res_data[j+1][1]) not in cutoff_matrix.keys():
                    key = (res_data[j+1][1], res_data[i+1][1])
                else:
                    key = (res_data[i+1][1], res_data[j+1][1])
                defcutoff = float(cutoff_matrix[key][0])
            # else, the cutoff is as defined and non specific to a given type of residue pair.
            else:
                try:
                    defcutoff = float(cutoff)
                except ValueError:
                    print("Program aborted")
                    print("The cutoff must be a number (Angstrom) or 'matrix'.")
                    sys.exit(0)
            # Make sure Dmax is lower than the distance cutoff:
            # The cutoff is converted in nm here (*0.1)
            if Dmax < defcutoff*0.1:
                # optional additional verbose printing:
                if verbose: print("Fitting distribution...")
                # Fit to model and extract the parameters:
                fit_result = fit_distogram(bins, 
                                       distograms[i][j][:-1], 
                                       model=model, 
                                       plots=plots,
                                       plot_labels=((i+1,j+1), (res_data_md[i+1][0], res_data_md[j+1][0])),
                                       nmodes=nmodes_max, 
                                       dmax=Dmax, 
                                       tolerr=tolerr, 
                                       verbose=verbose)
                # fit result = (nb_modes, mus, sigmas, scalings)
                # If the fitting failed (Normally it doesn't happen): 
                if fit_result == None:
                    print("Fitting ", i+1, j+1, " failed.")
                    continue
                # results[(res_i, res_j, atom_i, atom_j)] ==> (Dmax, nb modes, mus, sigmas, weights)
                results[(i+1, j+1, res_data_md[i+1][0], res_data_md[j+1][0])] = (Dmax, *fit_result)
    if plots:
        plots.close()
    return results


def write_output(results, nmodes_max=1, model="gauss", output_file="output.out", ndx_file="atom_list.ndx"):
    """
    Output file writing function
    Takes the fitting results as input and create and write the text file containing:
    the PLUMED-compatible information: the atom pairs and:
    the necessary parameters of the model function of the likelihood.
    The weights are added for multimodal functions. (the sum = 1).
    """
    models = {"gauss":"gaussian", "lognorm": "lognormal"}
    with open(output_file, "w") as f:
        # write header:
        f.write("#! FIELDS Id atom_i atom_j")
        if nmodes_max == 1:
            f.write(" mu sigma")
        else:
            for k in range(nmodes_max):
                f.write(" mu_{} sigma_{} weight_{}".format(k+1, k+1, k+1))
        f.write("\n")
        f.write("#! SET model "+models[model]+"\n")
        # initialize the list of atoms involved in restraints
        atom_list = []
        for i, (k, el) in enumerate(results.items()):
            atom_list.append(k[2])
            atom_list.append(k[3])
            # Extract the parameters from the result:
            Dmax, modes, mus, sigmas2, scalings = el[0], el[1], el[2], el[3], el[4]
            # write corresponding output row with the atom numbers i, j followed by the fitted distogram info Dmax, mu_i, sigmas2_i, weights_i:
            # Id, atoms i and j:
            f.write(str(i+1)+" "+str(k[2])+" "+str(k[3]))
            # info for each mode:
            if nmodes_max == 1:
                f.write(" {:.6f} {:.6f}".format(round(mus[0], 6), round(np.sqrt(sigmas2[0]), 6)))
            else:
                # transform scalings into weights:
                scalings = np.array(scalings)/np.sum(scalings)
                for l in range(nmodes_max):
                    if l < modes:
                        f.write(" {:.6f} {:.6f} {:.6f}".format(round(mus[l], 6), round(np.sqrt(sigmas2[l]), 6), round(scalings[l], 6)))
                    else:
                        f.write(" 0.000000 0.000000 0.000000")
            f.write("\n")
    # Final list of involved atom:
    atom_list = sorted(list(set(atom_list)))
    # write gromacs index file for PLUMED:
    with open(ndx_file, "w") as f:
        f.write("[ batoms ]\n")
        for i, el in enumerate(atom_list):
            f.write(str(el)+" ")
            if (i+1) % 15 == 0 or el == atom_list[-1]:
                f.write("\n")
    return 0


def main(pdb_file, pkl_file, output_file="output.txt", cutoff_A=8, model="gauss", seqsep=3, plots=False, verbose=False, mdpdb='same',
         nmodes=1, tolerr=0.005, convf="softmax", plotout="output.pdf", chains='all', ndxout="atom_list.ndx"):
    """
    Main function
    Reads the input files, process them according to the input parameters and writes the output file.
    """
    # optional additional verbose printing:
    if verbose: print("Reading pdb/cif file...")
    # obtain info on CB (or CA for GLY) atom index per residue in the MD pdb file:
    res_data = read_input(pdb_file)
    # check if an MD pdb file was given and read it if yes:
    if mdpdb != 'same':
        res_data_md = read_input(mdpdb)
    else:
        res_data_md = res_data
    # optional additional verbose printing:
    if verbose: print("Reading pkl file...")
    # obtain the bins and the distograms from the pickle file:
    dist_format = 'pkl' if pkl_file.endswith('.pkl') else 'npy' if pkl_file.endswith('.npy') else 'other'
    if dist_format == 'pkl':
        bins, distograms = read_pkl(pkl_file, convf=convf)
    elif dist_format == 'npy':
        bins, distograms = read_npy(pkl_file)
    else:
        raise ValueError("distograms file format not recognized. Try with a .pkl file or a .npy file.")
    # selecting and fitting the selected distogram pairs:
    results = select_and_fit_distograms(bins, distograms, res_data, res_data_md, cutoff=cutoff_A, model=model, seqsep=seqsep, plots=plots, verbose=verbose, nmodes_max=nmodes, tolerr=tolerr, plotout=plotout, chains=chains)
    # output function:
    write_output(results, nmodes_max=nmodes, model=model, output_file=output_file, ndx_file=ndxout)
    if verbose: print("Done.")
    return 0


if __name__ == "__main__":
    # Get rid of matplotlib multiplot RuntimeWarning:
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    # Read input arguments
    parser = build_parser()
    args = parser.parse_args()
    # display help if no input:
    if args.pdb == None or args.pkl == None:
        parser.print_help()
        sys.exit(0)
    # run main function:
    main(args.pdb, args.pkl, args.out, args.cutoff, args.model, args.seqsep, args.plots, args.verbose, args.mdpdb, 1, 0.005, 'softmax', args.plotout, args.chains, args.ndxout)
