from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def rel_ASA(file_in, chain, position):
    """This function uses DSSP program to calculate RSA, phi, and psi"""
    p = PDBParser(QUIET=True)
    structure = p.get_structure('xxx', file_in)
    model = structure[0]
    dssp = DSSP(model, file_in, dssp='mkdssp')
    #print(dssp[(chain, (' ', position, ' '))][3:6])
    return list(dssp[(chain, (' ', position, ' '))][3:6])

'''p = PDBParser(QUIET=True)
structure = p.get_structure("xxx", 'data/protein_pdbs2648/1lve_modif.pdb')
model = structure[0]

dssp = DSSP(model, 'data/protein_pdbs2648/1lve_modif.pdb', dssp='mkdssp')
# DSSP data is accessed by a tuple (chain_id, res_id)
a_key = list(dssp.keys())[2]
# (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
# NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
# NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
print(a_key)
print(dssp[a_key])
print(dssp)
print(dssp[('A', (' ', 27, ' '))])'''