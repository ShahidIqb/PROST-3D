from rdkit.Chem import rdFreeSASA
from rdkit import Chem


def compute_sasa(mol):
    """Compute Solvent Accessible Surface Area.
    """
    mol = Chem.AddHs(mol)

    # Get Van der Waals radii (angstrom)
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

    # Compute solvent accessible surface area
    sa = rdFreeSASA.CalcSASA(mol, radii, confIdx=-1)

    return sa


def compute_vdwsa(mol):
    """Compute Van der Waals Surface Area.
    """
    mol = Chem.AddHs(mol)
    radii = rdFreeSASA.classifyAtoms(mol)
    sa = rdFreeSASA.CalcSASA(mol, radii, confIdx=-1)
    return sa
