import pdbfixer
import tempfile
from openmm.app import PDBFile
import openmm.app as app


def gen_pdb_mutate(path_data, pdb_id, chain, mut_info, mut_resid):
    """This function mutates the wild residue type to mutant residue and fix the pdb block accordingly. This function
    uses PDBFixer """
    temp_file = path_data + pdb_id.lower() + '.pdb'
    fixer0 = pdbfixer.PDBFixer(filename=temp_file)

    fixer0.findMissingResidues()
    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer0.topology.chains())
    keys = fixer0.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chainA = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chainA.residues()))):
            missingResidues[key] = fixer0.missingResidues[key]
    fixer0.missingResidues = missingResidues

    fixer0.findMissingAtoms()
    fixer0.addMissingAtoms()
    PDBFile.writeFile(fixer0.topology, fixer0.positions, open(temp_file, 'w'), keepIds=True)
    #print(pdb_id.lower(), chain, mut_info)
    temp_mutfile = path_data + pdb_id.lower() + '_' +chain+ '_' + mut_info +'.pdb'
    #print(temp_mutfile)
    '''if os.path.exists(temp_mutfile):
        return'''
    fixer = pdbfixer.PDBFixer(filename=path_data + pdb_id.lower() + '.pdb')
    fixer.applyMutations([mut_info], chain)

    fixer.findMissingResidues()
    # only add missing residues in the middle of the chain, do not add terminal ones
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    missingResidues = dict()
    for key in keys:
        chain = chains[key[0]]
        if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
            missingResidues[key] = fixer.missingResidues[key]
    fixer.missingResidues = missingResidues

    #fixer.findMissingAtoms()
    #fixer.addMissingAtoms()
    with tempfile.NamedTemporaryFile(mode='w+') as temp_pdb:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, temp_pdb)
        temp_pdb.flush()
    PDBFile.writeFile(fixer.topology, fixer.positions, open(temp_mutfile, 'w'), keepIds=True)
