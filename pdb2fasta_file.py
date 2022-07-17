import os
def translate_3aa1(three_letter):
    """A dictionary to convert 3 letter amino acids to one letter. """
    if len(three_letter) > 3:
        three_letter = three_letter[1:4]
    trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
             'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
             'TYR': 'Y', 'VAL': 'V'}
    return (trans[three_letter])


#This code is based on the code from https://github.com/dongshuyan/pdb2fasta/blob/master/pdb2fasta_file.py
def pdbTofasta_seqpos(pdbfile, chain, wild_res, mut_pos, dir):
    """Convert PDB to FASTA sequence and map the structure's mutation position to sequence's position."""
    name = (pdbfile.split('.', 1)[0])
    name = name.split('/')[len(name.split('/')) - 1]
    filename = name + '_' + chain + '.fasta'
    if os.path.exists(filename):
        print(filename+" already existed.")
    else:
        Aname = '>' + name + '_' + chain+'\n'
        f = open(dir + filename, "w")
        f.write(Aname)
        prev = '-1'
        input_file = open(pdbfile)
        for line in input_file:
            toks = line.split()
            if len(toks) < 1: continue
            if toks[0] != 'ATOM': continue
            if toks[4] != chain: continue
            if toks[5] != prev:
                f.write('%c' % translate_3aa1(toks[3]))
            prev = toks[5]
        f.write('\n')
        f.close()
    # print name
    # print '>',name[0:len(name)]
    prev = '-1'
    mut_seqpos = -99999
    count_residues = 0
    input_file = open(pdbfile)
    for line in input_file:
        toks = line.split()
        if len(toks) < 1: continue
        if toks[0] != 'ATOM': continue
        if toks[4] != chain: continue
        if toks[5] != prev:
            count_residues += 1
            if mut_pos == toks[5] and wild_res == translate_3aa1(toks[3]):
                mut_seqpos = count_residues
        prev = toks[5]
    input_file.close()
    return mut_seqpos
