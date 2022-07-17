from csv_reader import read_csv
from fetch_pdb import fetchstructure
import pdb2fasta_file
import os
import csv

#################################################################################################################
# PPI_dimer_pdb2fasta: The program converts PDB to FASTA and maps pdb's mutation position to fasta's position. #
#################################################################################################################

# Configure these variables accordingly
input_mutfile_csv = 'ppi_mut_data/skempi2_mut_exp_data.csv'  # Input a CSV file
output_mutfile_csv = 'ppi_mut_data/mut_exp_data_seqpos.csv'  # Output csv file, last column with sequence-based mutation position
ppi_pdb_dir = 'ppi_mut_data/ppi_dimer_pdb/'  # Directory to store downloaded pdb files
ppi_fasta_dir = 'ppi_mut_data/ppi_dimer_fasta/'  # Directory to store converted fasta files

if __name__ == '__main__':
    mut_data = read_csv(input_mutfile_csv)
    header = mut_data.pop(0)
    header.append("seq_pos")
    print(header)
    with open(output_mutfile_csv, 'w', encoding='UTF8', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(header)
        for r in mut_data:
            pdbid = r[0]    # pdb id column
            mut_chain = r[1]  # Chain of mutation colum
            wild_res = r[3]   # Wild residue
            mut_pos_pdb = r[4]  # position of mutation in pdb file
            pdbid_partners = r[7]  # dimer information
            partnerA = ''
            partnerB = ''
            # Extracting dimer info
            if len(pdbid_partners)==8:
                partnerA = pdbid_partners.split("_")[1]
                partnerB = pdbid_partners.split("_")[2]
                struct_path = ppi_pdb_dir + pdbid + '.pdb'
                if os.path.exists(struct_path):
                    mut_seq_pos = pdb2fasta_file.pdbTofasta_seqpos(struct_path, mut_chain, wild_res, mut_pos_pdb,
                                                                   partnerA, partnerB, ppi_fasta_dir)
                    r.append(str(mut_seq_pos))
                else:
                    print("Downloading "+pdbid+" from RCSB PDB")
                    pdb_block = fetchstructure(pdbid)
                    with open(struct_path, 'w') as struct_file:
                        struct_file.write(pdb_block)
                    mut_seq_pos = pdb2fasta_file.pdbTofasta_seqpos(struct_path, mut_chain, wild_res, mut_pos_pdb,
                                                                   partnerA, partnerB, ppi_fasta_dir)
                    r.append(str(mut_seq_pos))
                print(r)
                csv_writer.writerow(r)
