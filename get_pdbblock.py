import sys

import numpy as np
from pharmaco_sign import translate_3aa1

def distance_euc(a, b):
    """Euclidean distance between two points"""
    return np.linalg.norm(a-b)


def extract_residueaswhole_within_tresh(path_data, pdb_id, chain, xyz_mean, min_distance, max_distance):
    """This function reads all the atoms of the residue when CA is in the threshold"""
    residue_details = set()
    for line in open(path_data + pdb_id + '.pdb'):
        line_elem = line.split()
        list_elem = line_elem
        if (line_elem[0] == 'ATOM' and line_elem[4] == chain and line_elem[2] == 'CB') or (
                line_elem[0] == 'ATOM' and line_elem[4] == chain and line_elem[3] == 'GLY' and line_elem[2] == 'CA'):
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
                x_cord = float(line_elem[6])
            if len(list_elem[7]) > 8:
                yz_comb = list_elem[7].split('-')
                if yz_comb[0] == '':
                    y_cord = float('-' + yz_comb[1])
                    z_cord = float('-' + yz_comb[2])
                else:
                    y_cord = float(yz_comb[0])
                    z_cord = float('-' + yz_comb[1])
            else:
                y_cord = float(line_elem[7])
            if z_cord==0.0:
                z_cord = float(line_elem[8])
            #print(x_cord, y_cord, z_cord)
            xyz = np.array((x_cord, y_cord, z_cord))
            if min_distance < distance_euc(xyz, xyz_mean) <= max_distance:
                residue_details.add(line_elem[4]+line_elem[5])
        if len(residue_details) > 0 and (line.split()[0]=='TER' or line.split()[0]=='END'):
            break
    #print(residue_details)
    pdb_block_list = []
    pdb_block = ''
    for linee in open(path_data + pdb_id + '.pdb'):
        line_elem = linee.split()
        if line_elem[0] =='ATOM' and line_elem[4]==chain:
            if line_elem[4]+line_elem[5] in residue_details:
                pdb_block = pdb_block+linee
                pdb_block_list.append(linee)
        if len(pdb_block) > 0 and (linee.split()[0]=='TER' or linee.split()[0]=='END'):
            break
    return pdb_block, pdb_block_list


def read_pdbBlock(path_data, pdb_id, chain, residue, residue_position, min_distance, max_distance):
    """Read a block of protein in the given chain. The geometric mean of the wild residue at given position is calculate and all
    atoms are considered from minimum distance to maximum distance.
    Input: pdb_id, residue, residue_position, min_distance, max_distance
    Output: A block of PDB"""
    pdb_block = ''
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0
    d_count = 0
    for line in open(path_data+pdb_id+'.pdb'):
        list_elem = line.split()
        #print(list_elem)
        if list_elem[0]=='ATOM' and list_elem[5]==residue_position and translate_3aa1(list_elem[3])== residue and list_elem[4]==chain:
            x_cord =''
            y_cord=''
            z_cord=''
            if len(list_elem[6]) > 8:
                xy_comb = list_elem[6].split('-')
                if xy_comb[0]=='':
                    x_cord = '-'+xy_comb[1]
                    y_cord = '-' + xy_comb[2]
                else:
                    x_cord = xy_comb[0]
                    y_cord = '-'+xy_comb[1]
            else:
                x_cord = list_elem[6]
            if len(list_elem[7])>8:
                yz_comb = list_elem[7].split('-')
                if yz_comb[0]=='':
                    y_cord = '-'+yz_comb[1]
                    z_cord = '-' + yz_comb[2]
                else:
                    y_cord = yz_comb[0]
                    z_cord = '-'+yz_comb[1]
            else:
                y_cord = list_elem[7]
            if z_cord=='':
                z_cord = list_elem[8]
            #print(x_cord, y_cord, z_cord)
            x_sum += float(x_cord)
            y_sum += float(y_cord)
            z_sum += float(z_cord)
            d_count += 1
        if d_count > 0 and (line.split()[0]=='TER' or line.split()[0]=='END'):
            break
    if d_count == 0:
        print(
            "Incorrect mutation info, check chain: " + chain + " wild-type: " + residue + " and mutation position: " + residue_position)
        sys.exit()
    xyz_mean = np.array((x_sum/d_count, y_sum/d_count, z_sum/d_count))

    #print(x_sum, y_sum, z_sum)

    pdb_block = extract_residueaswhole_within_tresh(path_data,pdb_id, chain, xyz_mean, min_distance, max_distance)
    return pdb_block