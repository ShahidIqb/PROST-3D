from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole, MolDrawOptions
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)
IPythonConsole.drawOptions.minFontSize=20


def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

def standardize(mol):
    # follows the steps in
    # https://github.com/greglandrum/RSC_OpenScience_Standardization_202104/blob/main/MolStandardize%20pieces.ipynb
    # as described **excellently** (by Greg) in
    # https://www.youtube.com/watch?v=eWTApNX8dJQ
    #mol = Chem.MolFromSmiles(smiles)
    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    clean_mol = rdMolStandardize.Cleanup(mol)

    # if many fragments, get the "parent" (the actual mol we are interested in)
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)

    # try to neutralize molecule
    uncharger = rdMolStandardize.Uncharger()  # annoying, but necessary as no convenience method exists
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)

    # note that no attempt is made at reionization at this step
    # nor at ionization at some pH (rdkit has no pKa caculator)
    # the main aim to to represent all molecules from different sources
    # in a (single) standard way, for use in ML, catalogue, etc.

    te = rdMolStandardize.TautomerEnumerator()  # idem
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)

    return taut_uncharged_parent_clean_mol

def view_difference(mol1, mol2, mutation_id, temp_ph_ddg, dir):
    Chem.SanitizeMol(mol1)
    Chem.SanitizeMol(mol2)
    AllChem.Compute2DCoords(mol1)
    AllChem.Compute2DCoords(mol2)
    Chem.AddHs(mol1)
    Chem.AddHs(mol2)
    #mol1 = standardize(mol1)
    #mol2 = standardize(mol2)
    disable_rdkit_logging()
    mol1 = rdMolStandardize.FragmentParent(mol1)
    mol2 = rdMolStandardize.FragmentParent(mol2)
    #Chem.rdCoordGen.AddCoords(mol1)
    #Chem.rdCoordGen.AddCoords(mol2)
    mcs = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    match1 = mol1.GetSubstructMatch(mcs_mol)
    target_atm1 = []
    for atom in mol1.GetAtoms():
        if atom.GetIdx() not in match1:
            target_atm1.append(atom.GetIdx())
    match2 = mol2.GetSubstructMatch(mcs_mol)
    target_atm2 = []
    for atom in mol2.GetAtoms():
        if atom.GetIdx() not in match2:
            target_atm2.append(atom.GetIdx())
    img = Draw.MolsToGridImage([mol1, mol2], subImgSize=(600, 600), highlightAtomLists=[target_atm1, target_atm2], returnPNG=False, legends=[temp_ph_ddg, ''] * (2000))
    img.save(dir+'/'+mutation_id+'.png')
    img.show()
