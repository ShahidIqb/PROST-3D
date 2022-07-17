from rdkit import Chem
from rdkit.Chem import AllChem
import mol_fix    # sanifix4 to set N aromaticity
from pmapper.pharmacophore import Pharmacophore as P
import openbabel

def pharmacophore_count(pdb_p_xyz):
    """This function returns the pharmacophore count"""
    pdb_p_count = pdb_p_xyz.get_features_count()
    return pdb_p_count


def pdbBlock2Mol(pdb_block):
    """This function uses openbabel to convert pdbBlock to mol format"""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "mdl")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, pdb_block)
    obConversion.WriteFile(mol, 'mol.mol')
    outMDL=obConversion.WriteString(mol)
    return


def rdkit_atomclasif(pdb_block, wild_res):
    """This function classify the atoms from the pdbBlock using the six atom classes."""
    pdbBlock2Mol(pdb_block)
    mol = Chem.MolFromMolFile('mol.mol', sanitize=True)
    mol_vis = mol
    '''AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol, 'mol_images/'+wild_res+'okk.png')
    Draw.MolToFile(mol, 'mol_images/okk.png')'''
    mol = mol_fix.fix_mol(mol)
    Chem.rdmolops.RemoveStereochemistry(mol)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    #AllChem.EmbedMolecule(mol, useRandomCoords=True, enforceChirality=False)
    AllChem.EmbedMultipleConfs(mol, numConfs=1, useRandomCoords=True, enforceChirality=False)
    vol = AllChem.ComputeMolVolume(mol)
    # create pharmacophore
    p = P()
    p.load_from_mol(mol)
    p_count = pharmacophore_count(p)
    return mol_vis, mol, p_count, vol

