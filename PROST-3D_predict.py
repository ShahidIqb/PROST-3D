import argparse
import sys
from get_pdbblock import read_pdbBlock
import pharmaco_pmapper
import pharmaco_sign
import mutate_pdb
from pdb2fasta_file import pdbTofasta_seqpos
import extract_features
from rdkit import Chem
from subprocess import getstatusoutput
import os
import csv_reader
#Modeller for mutation and building structure
import mutate_model
import pdbmol_prepare
import csv
import mol_difference
global blast_db, hhblits_db, blast_threads

# set these variables accordingly
blast_db = "$HOME/Downloads/uniref50db/uniref50"
hhblits_db = "$HOME/Downloads/uniclust30_2018_08/uniclust30_2018_08"
blast_threads = 12

def read_mutfile(file_name):
    """Reading mutation file"""
    mut_list = []
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            mut_list.append(line.strip().split(' '))
    file.close()
    return mut_list

def get_options():
    """Argument parsing"""
    parser = argparse.ArgumentParser(usage='Invalid arguments.',
                                     description='PROST-3D program for predicting the change in stability (∆∆G(kcal/mol)) upon single point missense mutation.')
    parser.add_argument('-pdb_id','-pdbid','--pdb_id','--pdbid', type=str, dest='pdb_id', help='Input the RCSB PDB ID of protein.')
    parser.add_argument('-file', '--file', action='store', type=str, dest='input_pdb_file', help='Input PDB formatted file of protein structure.')
    parser.add_argument('-mutation','--mutation', nargs='+', type=str, help='Inline mutation (example: A E 32 G 25 7)')
    parser.add_argument('-mutlist','--mutlist','--ml', '--mutation_list', action='store', type=str, dest='ml', help='Mutation-list (as provided in the Input directory).')
    parser.add_argument('-outdir',"--outdir", "--out_dir", action="store",type=str, dest="outdir", help="Output directory")
    parser.add_argument("-out_file", "--out_file", "-outfile" , action="store", type=str, dest="outfile", help="Output result file")
    args = parser.parse_args()

    if args.pdb_id:
        pdb_id = args.pdb_id
    elif not os.path.isfile(args.input_pdb_file):
        print('Error: Incorrect pdb file or no PDB ID provided')
        sys.exit(1)
    else:
        pdb_id = args.input_pdb_file
    print('Protein:', pdb_id)
    mut_file = None
    if args.ml:
        mut_data = read_mutfile(args.ml)
        mol_vis = False
    else:
        mut_data =[]
        mut_data.append(args.mutation)
        mol_vis = True
    outdir = os.getcwd()
    outdir = outdir+'/Result'
    outfile = 'Result'
    if args.outdir: outdir = args.outdir
    if args.outfile: outfile = args.outfile
    return mol_vis, pdb_id, mut_data, outdir, outfile

def translate_1aa3(one_letter):
    """Converts 1 letter amino-acid to 3 letter amino-acid"""
    trans = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS',
             'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
             'Y': 'TYR', 'V': 'VAL'}
    return (trans[one_letter])


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


if __name__ == '__main__':
    print("###############################")
    print("#   PROST-3D is running... #")
    print("###############################")
    mol_vis, pdb_id, mut_data, outdir, outfile = get_options()     # parse arguments
    path_data = "aux_files/"     # storage for supplementary files
    features_list = []
    X_test = []
    smart_ddg = []
    #mut_data = read_csv("test_data_PROST/Myoglobin.csv")
    #mut_data.pop(0)
    a = 0
    last_pdb_id = None
    strtpt = 0
    mut_list = []
    for rec in mut_data:
        #pdb_id = rec[0]
        mut_list.append(rec[0]+' '+rec[1]+' '+rec[2]+' '+rec[3]+' '+rec[4]+' '+rec[5])
        chain = rec[0]
        mut_pos = rec[2]
        mut_seqpos = int(rec[2])+1
        wild_res = rec[1]
        mutant_res = rec[3]
        mut_Temp = rec[4]
        if mut_Temp == '':
            mut_Temp = '25'
        mut_pH = rec[5]
        if mut_pH == '':
            mut_pH = '7'
        mut_info = translate_1aa3(wild_res)+'-'+mut_pos+'-'+translate_1aa3(mutant_res)
        # Fetching pdb structure from RCSB Protein Data Bank using a function from oddt(open drug descovery toolkit)
        if not os.path.exists(path_data + pdb_id.lower() + '.pdb'):
            mol = pdbmol_prepare.FetchStructure(pdb_id)
            Chem.MolToPDBFile(mol, path_data + pdb_id.lower() + '.pdb')

        #save as fasta and return the position of mutation in the sequence
        mut_seqpos = pdbTofasta_seqpos(path_data+pdb_id.lower() + '.pdb', chain, wild_res, mut_pos, path_data)

        wild_PDB_block, wild_PDB_block_list = read_pdbBlock(path_data, pdb_id.lower(), chain, wild_res, mut_pos, 0, 8)
        wild_mol_vis, wild_mol, wild_pdb_p_count, wild_vol = pharmaco_pmapper.rdkit_atomclasif(wild_PDB_block, wild_res)
        wild_pharmaco_sign = pharmaco_sign.get_pharmaco_sign(wild_PDB_block_list)

        # modeller for mutation
        '''sys.stdout = open(os.devnull, 'w')  # blocks output from Modeller
        mutate_model.mutation(pdb_id.lower(), mut_pos, translate_1aa3(mutant_res), chain, mut_info)
        sys.stdout = sys.__stdout__  # allows to print on display'''
        #print(wild_pharm_coord)
        #print(mCSM_sig)
        print(pdb_id, chain, mut_info)
        # PDBFixer for mutation
        if not os.path.exists(path_data + pdb_id.lower() + '_' + chain + '_' + mut_info + '.pdb'):
            mutate_pdb.gen_pdb_mutate(path_data, pdb_id, chain, mut_info, mutant_res)
            # modeller for mutation
            '''sys.stdout = open(os.devnull, 'w')  # blocks output from Modeller
            mutate_model.mutation(pdb_id.lower(), mut_pos, translate_1aa3(mutant_res), chain, mut_info)
            sys.stdout = sys.__stdout__  # allows to print on display'''

        mutant_PDB_block, mutant_PDB_block_list = read_pdbBlock(path_data, pdb_id.lower()+'_'+chain+'_'+mut_info, chain, mutant_res, mut_pos, 0, 8)
        #print(mutant_PDB_block)
        mutant_mol_vis, mutant_mol, mutant_pdb_p_count, mutant_vol = pharmaco_pmapper.rdkit_atomclasif(mutant_PDB_block, mutant_res)
        mutant_pharmaco_sign = pharmaco_sign.get_pharmaco_sign(mutant_PDB_block_list)

        #path_PROSTevo = 'myoglobin'

        #feature extraction for the provided mutation position
        features_list = extract_features.mut_features(pdb_id, chain, mut_pos, mut_pH, mut_Temp, wild_res, mutant_res,
                                                      mut_seqpos, wild_vol, mutant_vol, wild_pharmaco_sign,
                                                      mutant_pharmaco_sign, wild_pdb_p_count, mutant_pdb_p_count,
                                                      wild_mol, mutant_mol, path_data, blast_db, blast_threads, hhblits_db)
        smart_ddg.append([Chem.MolToSmiles(wild_mol), ''])
        X_test.append(features_list)
        '''AllChem.Compute2DCoords(wild_mol)
        Draw.MolToFile(wild_mol, 'mol_images/'+ pdb_id.lower() + '_' + chain + '_' + mut_pos+wild_res+'_wild.png')
        AllChem.Compute2DCoords(mutant_mol)
        Draw.MolToFile(wild_mol, 'mol_images/'+ pdb_id.lower() + '_' + chain + '_' + mut_info +'.png')'''
        if mol_vis:
            wild_PDB_block, wild_PDB_block_list = read_pdbBlock(path_data, pdb_id.lower(), chain, wild_res, mut_pos, 0,
                                                                5)
            wild_mol_vis, wild_mol, wild_pdb_p_count, wild_vol = pharmaco_pmapper.rdkit_atomclasif(wild_PDB_block,
                                                                                                   wild_res)
            mutant_PDB_block, mutant_PDB_block_list = read_pdbBlock(path_data,
                                                                    pdb_id.lower() + '_' + chain + '_' + mut_info,
                                                                    chain, mutant_res, mut_pos, 0, 5)
            # print(mutant_PDB_block)
            mutant_mol_vis, mutant_mol, mutant_pdb_p_count, mutant_vol = pharmaco_pmapper.rdkit_atomclasif(
                mutant_PDB_block, mutant_res)
            mol_difference.view_difference(wild_mol_vis, mutant_mol_vis, pdb_id.lower() + '_' + chain + '_' + mut_info,
                                           "Highlighted Changes ---- Right(Wild-type)    Left (Mutant-Type)   " + pdb_id + '\t' + str(
                                               mut_info), outdir)
            break

    file_mol_ddg = open('temp_mol_ddg.csv', 'w', newline='')
    write_file_mol_ddg = csv.writer(file_mol_ddg)
    write_file_mol_ddg.writerow(['smarts', 'ddG'])
    #for smart_ddg_ent in smart_ddg:
    write_file_mol_ddg.writerows(smart_ddg)
    file_mol_ddg.close()
    file_features = open('temp_features.csv', 'w', newline='')
    write_file_features = csv.writer(file_features)
    write_file_features.writerow('mutation_features')
    write_file_features.writerows(X_test)
    file_features.close()
    out = getstatusoutput('./predict.sh')
    if out[0] != 0:
        print('Error in prediction')
    temp_pred = csv_reader.read_csv('temp_pred.csv')
    temp_pred.pop(0)
    file_pred = open(outdir+'/'+outfile+'.txt', 'w+')
    print('#Chain  Wild-type  Position  Mutant-type  Temp(ºC) pH  predicted ∆∆G(kcal/mol)')
    file_pred.write('#Chain Wild-type Position Mutant-type Temp(ºC) pH\t predicted ∆∆G(kcal/mol)\n')
    for i in range(len(mut_list)):
        print(mut_list[i], temp_pred[i][1])
        file_pred.write("%s\t%s\n" % (mut_list[i], temp_pred[i][1]))
    file_pred.close()
