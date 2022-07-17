# RPOST-GDNet
A tool for predicting the effects of missense mutation on protein stability changes using 3D structural information of protein. The tool is based on Direct Message Passing Neural Network (a graph convolutional neural network) and Deep Learning Neural Network. PROST-GDNet also visualizes and highlights the change in molecular strucutre for a single missense mutation.

![Screen Shot 2022-07-17 at 11 17 48 pm](https://user-images.githubusercontent.com/48677766/179400249-65018ed2-9f77-4031-8c08-d8031e2abdbe.png)

Requirements: The requirements.yml file is provided.

Installation of anaconda3 is prerffered
1. create PROST-GDNet environment from the provided requirements.yml file.
2. Download the following databases and make them ready for searching
  i) uniref50 (https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50) [make this ready for blast by using the following command]
      makeblastdb -in uniref50.fasta -dbtype prot -out uniref50
  ii) uniclust30_2018_08 (http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz)
  
  2.1) Now, set the paths (in line 21-22 in PROST-GDNet_predict.py) for the installed databases accordingly.
  
3. Activate your virtual environment(PROST-GDNet) and run the python script PROST-GDNet_predict.py

  Command-line arguments for the program:
  
    {-pdbid,--pdbid, -pdb_id, --pdb_id} PDB ID of protein from RCSB Protein Data Bank
    
    {-file,--file}	protein sequence (FASTA format)

    {-mutation, --mutation}	missence mutation (example: A V 8 I 25 7 or A T 77 H 25 7)

    {-mutlist, --mutlist, --ml, --mutation_list}	list of mutations

    {-outdir, --outdir, --out_dir}	directory name for results

    {-out_file, --out_file, -outfile} Name for the result output file

    {-h, --help}	command-line summary
    

Prediction for a single mutation
      
      python PROST-GDNet_predict.py -pdb_id RCSB PDB ID --mutation wild-type position mutant-type temp°C(optional) pH(optional) --outdir (optional)       
      Result --out-file (optional) mutation_result
    python PROST-GDNet_predict.py -file Path_To_PDB_Structure --mutation wild-type position mutant-type temp°C(optional) pH(optional) --outdir (optional) 
    Result --out-file (optional) mutation_result
Prediction for list of mutations
       python PROST-GDNet_predict.py -pdb_id RCSB PDB ID -mutlist Path_To_Mutation_List -outdir(optional) Result -out-file(optional) mut_list_Result
       python PROST-GDNet_predict.py -file Path_To_PDB_Structure -mutlist Path_To_Mutation_List -outdir(optional) Result -out-file(optional) 
       mut_list_Result
       
4. Examples of how to run a single mutation and list of mutations
  1. python PROST-GDNet_predict.py -pdb_id 1aky --mutation A V 8 I 25 7  --outdir Result --out-file 1aky_prediction
  2. python PROST-GDNet_predict.py -pdb_id 1aky --mutation-list Input/1aky_mutlist.txt  --outdir Result --out-file 1aky_prediction

**Auxiliary files will be stored inside aux_files directory.
