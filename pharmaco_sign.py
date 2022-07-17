from itertools import combinations
import numpy as np


def translate_3aa1(three_letter):
    """Converts 3 letter amino-acid to 1 letter amino-acid"""
    if len(three_letter) > 3:
        three_letter = three_letter[1:4]
    trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
             'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
             'TYR': 'Y', 'VAL': 'V', 'APHE': 'F', 'ATRP':'W', 'ASER':'S', 'AGLU': 'E'}
    return (trans[three_letter])


def atom_coordinates(PDB_block):
    '''Coordinates for polar and hydrophobic Ca atoms from the PDB_Block'''
    #print(PDB_block)
    polar_group = ('R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S', 'T', 'Y')
    hydrophobic_group = ('A', 'C', 'I', 'L', 'M', 'F', 'W', 'V')
    label_coordinates = []
    for line in PDB_block:
        list_elem = line.split()
        if list_elem[2] == 'CA':
            x_cord = 0.0
            y_cord = 0.0
            z_cord = 0.0
            if len(list_elem[6]) > 8:
                xy_comb = list_elem[6].split('-')
                if xy_comb[0] == '':
                    x_cord = float('-' + xy_comb[1])
                    y_cord = float('-' + xy_comb[2])
                else:
                    x_cord = float(xy_comb[0])
                    y_cord = float('-' + xy_comb[1])
            else:
                x_cord = float(list_elem[6])
            if len(list_elem[7]) > 8:
                yz_comb = list_elem[7].split('-')
                if yz_comb[0] == '':
                    y_cord = float('-' + yz_comb[1])
                    z_cord = float('-' + yz_comb[2])
                else:
                    y_cord = float(yz_comb[0])
                    z_cord = float('-' + yz_comb[1])
            else:
                y_cord = float(list_elem[7])
            if z_cord == 0.0:
                z_cord = float(list_elem[8])
            if translate_3aa1(list_elem[3]) in hydrophobic_group:
                label_coordinates.append(('H', (float(x_cord), float(y_cord), float(z_cord))))
            if translate_3aa1(list_elem[3]) in polar_group:
                label_coordinates.append(('D', (float(x_cord), float(y_cord), float(z_cord))))
    return label_coordinates


def distance_euc(a, b):
    """Euclidean distance between two points"""
    return np.linalg.norm(a-b)


def distance_xyz(p1, p2):
    """Euclidean distance between labels with two points."""
    label1, x1y1z1 = p1
    x1, y1, z1 = x1y1z1
    label2, x2y2z2 = p2
    x2, y2, z2 = x2y2z2
    return (label1+'_'+label2, distance_euc(np.array([x1, y1, z1]), np.array([x2, y2, z2])))


def getFrequency(distMatrix, d_step):
    """Frquency of atom classes with cutoff value"""
    atom_classes1 = {'A_A':0, 'A_D':1, 'D_A':1, 'A_H':2, 'H_A':2,'A_N':3, 'N_A':3,'A_a':4, 'a_A':4, 'A_P':5, 'P_A':5,
                    'D_D':6, 'D_H':7, 'H_D':7,'D_N':8, 'N_D':8,'D_a':9, 'a_D':9, 'D_P':10, 'P_D':10,
                    'H_H':11,'H_N':12, 'N_H':12,'H_a':13, 'a_H':13, 'H_P':14, 'P_H':14,
                    'N_N':15,'N_a':16, 'a_N':16, 'N_P':17, 'P_N':17,
                    'a_a':18, 'a_P':19, 'P_a':19,
                    'P_P':20}
    cmbnation = 21
    atom_classes2 = {'A_A':0, 'A_D':0, 'D_A':0, 'A_N':0, 'N_A':0, 'A_P':0, 'P_A':0,'A_H':1, 'H_A':1,
                    'D_D':0, 'D_N':0, 'N_D':0,'D_P':0, 'P_D':0,'D_H':1, 'H_D':1,
                    'N_N':0, 'N_P':0, 'P_N':0, 'N_H':1, 'H_N':1,
                    'P_P':0, 'P_H':1, 'H_P':1,
                    'H_H':2}
    atom_classes3 = {'A_A':0, 'A_D':0, 'D_A':0, 'A_H':1, 'H_A':1,
                    'D_D':0, 'D_H':1, 'H_D':1,
                    'H_H':2}
    cmbnation = 3
    scan_matrix = [0] * cmbnation * int(8/d_step)
    i = 0
    for dist in range(2, 8+1, d_step):
        for rec in distMatrix:
            if rec[1] <= dist and rec[0] in atom_classes3:
                scan_matrix[atom_classes3.get(rec[0])+i*cmbnation] += 1
        i += 1
    return scan_matrix


def get_pharmaco_sign(pdb_block_list):
    pharmaco_coords = atom_coordinates(PDB_block=pdb_block_list)
    distances = []
    for combo in combinations(pharmaco_coords, 2):
        distances.append(distance_xyz(*combo))
    mCSM_sign = getFrequency(distances, 2)
    return mCSM_sign
