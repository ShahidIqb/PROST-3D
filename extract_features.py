from features import *
import dssp_features as dssp
import mol_surf_area as rdk_sa
from ddgun_seq import ddgun_features
from subprocess import getstatusoutput
import os, sys

def mCSM_segDiff(wild_sig, mutant_sig):
    """Difference between wild-type and mutant-type mCSM signatures"""
    return [a_i - b_i for a_i, b_i in zip(wild_sig, mutant_sig)]


def pharmcophore_count_diff(wild_pharm_count, mutant_pharm_count):
    """Calculate the difference between to pharmacophore counts"""
    atom_class = ('A', 'D', 'H', 'N', 'a', 'P')
    pharm_w = [0] * 6
    i = 0
    for item in atom_class:
        if item in wild_pharm_count:
            pharm_w[i] = wild_pharm_count.get(item)
        i += 1
    pharm_m = [0] * 6
    i = 0
    for item in atom_class:
        if item in mutant_pharm_count:
            pharm_m[i] = mutant_pharm_count.get(item)
        i += 1

    return [a_i - b_i for a_i, b_i in zip(pharm_w, pharm_m)]


def job_submit_pssm(file_name, dir, blast_db, blast_threads):
    if os.path.exists(dir+'/'+file_name+'.pssm'):
        return
    print("1) Run PSI-BLAST Search")
    print("2) Build "+file_name+" PSSM...")
    cmd = "psiblast -query "+dir+'/'+file_name+'.fasta'+" -num_threads "+str(blast_threads)+" -db "+blast_db+" -num_iterations 3 -out "+dir+'/'+file_name+".out -out_ascii_pssm "+dir+'/'+file_name+".pssm 2>/dev/null"
    print(cmd)
    out = getstatusoutput(cmd)
    if out[0] != 0:
        print('Error: pssm profile not generated')


def pssm_check(file_name, dir):
    file_name += '.pssm'
    #print(file_name)
    if not os.path.exists(dir+'/'+file_name):
        print("Can't get the PSSM, please check the"+file_name+"sequence in aux_files directory")
        sys.exit()

def mut_features(pdb_id, chain, mut_pos, mut_pH, mut_Temp, wild_res, mutant_res, mut_seqpos,
                        wild_vol, mutant_vol, wild_mCSMFreq_simple, mutant_mCSMFreq_simple, wild_pdb_p_count,
                        mutant_pdb_p_count, wild_mol, mutant_mol, path_data, blast_db, blast_threads, hhblits_db):
    """Extract features for the provided mutation"""
    features = []
    features.append(float(mut_pH))
    #running dssp program to generate RSA, phi, and psi features
    dssp_f = dssp.rel_ASA(path_data + pdb_id.lower() + '.pdb', chain, int(mut_pos))
    features.append(dssp_f[0])
    features.append(dssp_f[1] / 100)
    features.append(dssp_f[2] / 100)
    features.append(float(mut_Temp) / 25)
    features.append((wild_vol-mutant_vol)/100)
    features = features + mCSM_segDiff(wild_mCSMFreq_simple, mutant_mCSMFreq_simple)
    features = features + wild_mCSMFreq_simple
    features = features + pharmcophore_count_diff(wild_pdb_p_count, mutant_pdb_p_count)
    features.append((rdk_sa.compute_sasa(wild_mol) - rdk_sa.compute_sasa(mutant_mol)) / 100)
    features.append((rdk_sa.compute_vdwsa(wild_mol) - rdk_sa.compute_vdwsa(mutant_mol)) / 100)
    features.append(net_volume.net_volume(wild_res, mutant_res))
    features.append(net_hydrophobicity.net_hydrophobicity(wild_res, mutant_res))
    features.append(mutation_hydrophobicity.mutation_hydrophobicity(wild_res, mutant_res))
    features.append(mutation_size.mutation_size(wild_res, mutant_res))
    features.append(mutation_hbonds.mutation_hbonds(wild_res, mutant_res))
    features = features + get_daaph7.get_daaph7(wild_res, mutant_res)
    features = features + z_scale.z_scale(wild_res, mutant_res)
    features.append(modif_gaac.modif_gaac(wild_res, mutant_res))
    job_submit_pssm(pdb_id + '_' + chain, path_data, blast_db, blast_threads)
    pssm_check(pdb_id + '_' + chain, path_data)
    features.append(get_dF.get_dF(pdb_id + '_' + chain, mut_seqpos, wild_res, mutant_res, path_data))
    features = features + mpssmscores.mpssmscores(pdb_id + '_' + chain, mut_seqpos, path_data)
    features = features + psepssm.psepssm(pdb_id + '_' + chain, path_data)
    features = features + ddgun_features(path_data+pdb_id+'_'+chain+'.fasta',wild_res,mut_seqpos,mutant_res,path_data,hhblits_db)

    return features


