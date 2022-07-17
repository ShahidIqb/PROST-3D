from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdDepictor, rdMolTransforms
from typing import List


def fix_bond_order(mol: Chem.Mol) -> Chem.Mol:
    """On a Mol where hydrogens are present it guesses bond order."""

    def is_sp2(atom: Chem.Atom) -> bool:
        N_neigh = len(atom.GetBonds())
        symbol = atom.GetSymbol()
        if symbol == 'H':
            return False
        elif symbol == 'N' and N_neigh < 3:
            return True
        elif symbol == 'C' and N_neigh < 4:
            return True
        elif symbol == 'O' and N_neigh < 2:
            return True
        else:
            return False

    def get_other(bond: Chem.Bond, atom: Chem.Atom) -> Chem.Atom:
        """Given an bond and an atom return the other."""
        if bond.GetEndAtomIdx() == atom.GetIdx():  # atom == itself gives false.
            return bond.GetBeginAtom()
        else:
            return bond.GetEndAtom()

    def find_sp2_bonders(atom: Chem.Atom) -> List[Chem.Atom]:
        return [neigh for neigh in find_bonders(atom) if is_sp2(neigh)]

    def find_bonders(atom: Chem.Atom) -> List[Chem.Atom]:
        return atom.GetNeighbours()

    def descr(atom: Chem.Atom) -> str:
        return f'{atom.GetSymbol()}{atom.GetIdx()}'

    ## main body of function
    for atom in mol.GetAtoms():
        # print(atom.GetSymbol(), is_sp2(atom), find_sp2_bonders(atom))
        if is_sp2(atom):
            doubles = find_sp2_bonders(atom)
            if len(doubles) == 1:
                # tobedoubled.append([atom.GetIdx(), doubles[0].GetIdx()])
                b = mol.GetBondBetweenAtoms(atom.GetIdx(), doubles[0].GetIdx())
                if b:
                    b.SetBondType(Chem.rdchem.BondType.DOUBLE)
                else:
                    raise ValueError('Issue with:', descr(atom), descr(doubles[0]))
            elif len(doubles) > 1:
                for d in doubles:
                    b = mol.GetBondBetweenAtoms(atom.GetIdx(), d.GetIdx())
                if b:
                    b.SetBondType(Chem.rdchem.BondType.AROMATIC)
                    b.SetIsAromatic(True)
                else:
                    raise ValueError('Issue with:', descr(atom), descr(d))
            elif len(doubles) == 0:
                print(descr(atom), ' is underbonded!')
        else:
            pass
            # print(descr(atom),' is single', find_bonders(atom))
    return mol