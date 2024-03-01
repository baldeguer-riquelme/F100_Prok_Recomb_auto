#!/usr/bin/env python

# Authors: Borja Aldeguer-Riquelme & Roth E. Conrad
# Contact: briquelme3@gatech.edu

'''
This script compiles into a single script the main steps and scripts of the F100 pipeline developed by Roth Conrad.
Ensure all scripts and tools are accessible direclty from the path, the script rely on them to perform its task.
Detailed information about the F100 pipeline can be found at:

https://github.com/rotheconrad/F100_Prok_Recombination/
'''

import os
import sys
import argparse
import glob
import shutil
import subprocess
import time
import pandas as pd
import random
from tqdm import tqdm
from datetime import datetime
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import Process
from zipfile import ZipFile


def call_rename(file):
    rename=subprocess.Popen(f'01a_rename_fasta_v2.py -i {file}', shell=True, stdout=subprocess.DEVNULL)
    rename.wait()

def run_rename(in_dir, out_dir, threads):
    # Create directory to place the input genomes if doesn't exist yet. 
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
        os.mkdir(f"{out_dir}/1_Preprocessing/")
        os.mkdir(f"{out_dir}/1_Preprocessing/Genomes/")
    else:
        # Exit if output directory exists
        print("Output directory already exists, remove or use a different name")
        exit()

    # List input genomes
    genomes_out_dir=glob.glob(f'{in_dir}*')

    # Get only the genome name, remove the remaining path and copy to output directory
    for genome in genomes_out_dir:
        genome_filename = genome.split("/")[-1]
        target = f'{out_dir}/1_Preprocessing/Genomes/{genome_filename}'
        shutil.copy(genome, target)

    # List input genomes in output directory and rename
    list_genomes = glob.glob(f'{out_dir}/1_Preprocessing/Genomes/*')

    ## Multiprocessing
    pool_size = threads
    pool = Pool(pool_size)

    pool.imap_unordered(call_rename, list_genomes)

    pool.close()
    pool.join()

    return(list_genomes)


def call_prodigal(list):
    for file in list:
        file_split = file.split("/")
        genome_filename = file_split[-1].rsplit(".", 1)[0]
        out_dir = '/'.join([file_split[0], file_split[1]])
        pyrodigal=subprocess.call(f'pyrodigal -i {file} -a {out_dir}/Genes_aa/{genome_filename}.faa -d {out_dir}/Genes_nt/{genome_filename}.fnn', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_prodigal(out_dir, threads):
    # Create output directories for genes and proteins
    os.mkdir(f"{out_dir}/1_Preprocessing/Genes_nt/")
    os.mkdir(f"{out_dir}/1_Preprocessing/Genes_aa/")
    
    # Get list of genomes
    list_genomes = glob.glob(f'{out_dir}/1_Preprocessing/Genomes/*')
    
    # Split list genomes
    chunk_size = int(round(len(list_genomes) / threads, 0))
    splitted_list = [list_genomes[i:i + chunk_size] for i in range(0, len(list_genomes), chunk_size)]
    #
    start = time.time()
    # Start processes in parallel
    processes_pyrodigal = [Process(target=call_prodigal, args=(i, )) for i in splitted_list]
    for process in processes_pyrodigal:
        process.start()
    # wait for all processes to complete
    for process in processes_pyrodigal:
        process.join()
    #
    end = time.time()
    print("Elapsed time: %s seconds" % (end - start))


def fastANI(query_list, subject_list, genome_sim_dir, num):
    _ = subprocess.call(f'fastANI --ql {query_list} --rl {subject_list} -o {genome_sim_dir}/fastANI_allV{num}.ani > /dev/null 2>&1', shell=True)


def parallel_fastANI(list_genomes, genome_file_list, genome_sim_dir, threads):
    # Split list genomes
    chunk_size = int(round(len(list_genomes) / threads, 0))
    splitted_list = [list_genomes[i:i + chunk_size] for i in range(0, len(list_genomes), chunk_size)]

    # Save each sublist to a file that will be used as a query list for fastANI
    list_query_files = []
    for x in range(0,len(splitted_list)):
        out_query_file = f'query_list{x}.txt'
        with open(out_query_file,'w') as file:
            for line in splitted_list[x]:
                file.write(line + '\n')
        file.close()
        list_query_files.append(out_query_file)

    # Start processes in parallel
    processes_fastANI = [Process(target=fastANI, args=(sublist, genome_file_list, genome_sim_dir, num, )) for num, sublist in enumerate(list_query_files)]
    for process in processes_fastANI:
        process.start()
    
    # wait for all processes to complete
    for process in processes_fastANI:
        process.join()
    
    # Join all fastANI results in one table
    ani_res = glob.glob(f'{genome_sim_dir}/*.ani')
    final_ani = f'{genome_sim_dir}/fastANI_allV.ani'
    with open(final_ani, 'w') as out:
        for ani in ani_res:
            with open(ani, 'r') as ani_file:
                text = ani_file.read()
                out.write(text)
            ani_file.close()
            os.remove(ani)
    out.close()

    # Remove query files
    for query_file in list_query_files:
        os.remove(query_file)


def run_fastANI(list_genomes, out_dir, threads):
    # Create dir for results
    genome_sim_dir = f'{out_dir}/1_Preprocessing/Genome_similarity'
    os.mkdir(genome_sim_dir)
    
    # Save genome files list
    genome_file_list = f'{genome_sim_dir}/genome_file_list.txt'
    with open(genome_file_list, 'w') as f:
        for file in list_genomes:
            f.write(file + '\n')
    f.close()
    
    # Run fastANI
    if len(list_genomes) < 50:
        num=''
        _ = fastANI(genome_file_list, genome_file_list, genome_sim_dir, num)
    else:
        _ = parallel_fastANI(list_genomes, genome_file_list, genome_sim_dir, threads)

    # Make plots based on ANI values
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 01b_fastANI_scatter_pyGAM.py...\n\n"]))
    _ = subprocess.call(f'01b_fastANI_scatter_pyGAM.py -i {genome_sim_dir}/fastANI_allV.ani -s {out_dir} -o fastANI_allV_sharedfrac_ANI.pdf -m True', shell=True)
    shutil.move('fastANI_allV_sharedfrac_ANI.pdf', genome_sim_dir)
    
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 01c_fastANI_clustermap.py...\n\n"]))
    _ = subprocess.call(f'01c_fastANI_clustermap.py -i {genome_sim_dir}/fastANI_allV.ani -o fastANI_allV_heatmap.pdf', shell=True)
    shutil.move('fastANI_allV_heatmap.pdf', genome_sim_dir)


def makedb(splitted_list):
    for genome in splitted_list:
        subprocess.call(f'makeblastdb -in {genome} -dbtype nucl', shell=True)


def RBM(splitted_list, full_list, num):
    short_list = splitted_list[num]
    if num == 0:
        pos=1
    elif num == 1:
        pos=2+len(splitted_list)
    elif num > 1:
        pos=((2+len(splitted_list))*(num))
    x=1
    for g1 in tqdm(short_list, desc=f'Main {num}', leave=False, position=pos):
        for g2 in tqdm(full_list, desc=f'Secondary {num}-{x}', leave=False, position=pos+x):
            out1 = g1.split('/')[-1].rsplit(".", 1)[0]
            out2 = g2.split('/')[-1].rsplit(".", 1)[0]
            subprocess.call(f'02b_get_RBMs_pipeline_auto.py -g1 {g1} -g2 {g2} -o ./tmp/{out1}-{out2}.rbm > /dev/null 2>&1', shell=True)
        x=x+1


def run_RBM(out_dir, threads):
    # Create dir for results
    RBM_dir = f"{out_dir}/1_Preprocessing/RBM_results"
    os.mkdir(RBM_dir)
    os.mkdir("./tmp")
    
    # Get list of files
    list_genes = glob.glob(f'{out_dir}/1_Preprocessing/Genes_nt/*')
    
    # Calculate chunk size and split list genes
    chunk_size = int(round(len(list_genes) / threads, 0))
    splitted_list = [list_genes[i:i + chunk_size] for i in range(0, len(list_genes), chunk_size)]
    #
    start = time.time()
    #
    print("Building blast databases for input genes..")
    processes_db = [Process(target=makedb, args=(i, )) for i in splitted_list]
    for process in processes_db:
        process.start()
    # wait for all processes to complete
    for process in processes_db:
        process.join()
    #
    end = time.time()
    print("Elapsed time: %s seconds" % (end - start))
    #
    
    print("Running blast...")
    start = time.time()
    processes_blast = [Process(target=RBM, args=(splitted_list, list_genes, i )) for i in range(0, len(splitted_list))]
    for process in processes_blast:
        process.start()
    # wait for all processes to complete
    for process in processes_blast:
        process.join()
    # Remove blast databases
    list_dbs = glob.glob(f'{out_dir}/1_Preprocessing/Genes_nt/.fnn.*')
    for db in list_dbs:
        os.remove(db)
    
    end = time.time()
    print("Elapsed time: %s seconds" % (end - start))


def join_RBM(out_dir):
    rbm_files = glob.glob("./tmp/*.rbm")
    
    with open(f"{out_dir}/1_Preprocessing/RBM_results/RBMs_allV.rbm", "wb") as out_file:
        for f in tqdm(rbm_files):
            with open(f, "rb") as in_file:
                out_file.write(in_file.read())
    
    shutil.rmtree('./tmp/')


def calculate_F100(out_dir, rec_threshold):
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 02c_get_F100.py...\n\n"]))
    subprocess.call(f'02c_get_F100.py -i {out_dir}/1_Preprocessing/RBM_results/RBMs_allV.rbm -o {out_dir} -rec {rec_threshold}', shell=True)
    shutil.move(f'{out_dir}_F100.tsv', f'{out_dir}/1_Preprocessing/RBM_results/')


def run_preprocessing(in_dir, out_dir, threads, rec_threshold):
    # Rename
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Renaming genomes...\n\n"]))
    list_genomes = run_rename(in_dir, out_dir, threads)
    
    # Predict genes and proteins with pyrodigal
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running prodigal...\n\n"]))
    _ = run_prodigal(out_dir, threads)
    
    # Estimate genome similarity with fastANI
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running fastANI...\n\n"]))
    _ = run_fastANI(list_genomes, out_dir, threads)
    
    # Calculate RBM
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Obtaining reciprocal best matches, this might take a while...\n\n"]))
    _ = run_RBM(out_dir, threads)
    
    # Join RBM into a single file
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Joining reciprocal best matches to a single file...\n\n"]))
    _ = join_RBM(out_dir)
    
    # Calculate F100 scores
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Calculating F100 scores...\n\n"]))
    _ = calculate_F100(out_dir, rec_threshold)


def run_models(out_dir, in_file, out_file, plot_title, xaxis_min, xaxis_max, yaxis_min, yaxis_max, model_file, prebuilt_folder, xaxis_step_size, point_size, point_alpha):
    # Create dir for results if doesn't exist
    out_models_dir = f"{out_dir}/2_Models"
    if os.path.isdir(out_models_dir) == False:
        os.mkdir(out_models_dir)
    
    # Run script depending on the variables provided by the user
    if prebuilt_folder is not None:
        # Models in folder
        if os.path.isdir(prebuilt_folder) == True:
            models_analyzed = []
            list_models = glob.glob(f'{prebuilt_folder}/*')
            for model in list_models:
                if os.path.isfile(model) == True and os.path.splitext(model)[0] not in list_models:
                    if os.path.splitext(model)[1] == ".zip":
                        # Extract zip file
                        with ZipFile(model, 'r') as zip_file:
                            list_names = zip_file.namelist()
                            zip_file.extractall(prebuilt_folder)
                        unzipped_model = os.path.splitext(model)[0]
                        models_analyzed.append(unzipped_model)
                        
                        # Rename zip file to match the original zip file
                        os.rename(f'{prebuilt_folder}/{list_names[0]}', unzipped_model)
                        
                        # Set outfile prefix and run script
                        out_prefix = f'{out_file}_{unzipped_model.split("/")[-1].rsplit(".", 1)[0]}'
                        subprocess.call(f'02d_f100_scatter_pyGAM.py -i {unzipped_model} -o {out_prefix} -t {plot_title} -xmin {xaxis_min} -xmax {xaxis_max} -ymin {yaxis_min} -ymax {yaxis_max} -s {xaxis_step_size} -p {point_size} -a {point_alpha} -i2 {in_file}', shell=True)
                    else:
                        # Just analyze the tsv file
                        models_analyzed.append(model)
                        out_prefix = f'{out_file}_{model.split("/")[-1].rsplit(".", 1)[0]}'
                        subprocess.call(f'02d_f100_scatter_pyGAM.py -i {model} -o {out_prefix} -t {plot_title} -xmin {xaxis_min} -xmax {xaxis_max} -ymin {yaxis_min} -ymax {yaxis_max} -s {xaxis_step_size} -p {point_size} -a {point_alpha} -i2 {in_file}', shell=True)
        else:
            print(f'{prebuilt_folder} is not a folder. Please, provide the path to a folder containing the pre-built models in .zip or .tsv')
    else:
        # Model in file
        if os.path.isfile(model_file) == True:
            subprocess.call(f'02d_f100_scatter_pyGAM.py -i {model_file} -o {out_file} -t {plot_title} -xmin {xaxis_min} -xmax {xaxis_max} -ymin {yaxis_min} -ymax {yaxis_max} -s {xaxis_step_size} -p {point_size} -a {point_alpha} -i2 {in_file}', shell=True)
        else:
            print('Exiting, there was no input model. Please, specify one!')
            sys.exit()
    # Move tsv and pdf files to the 2_Models folder
    list_output_files = glob.glob(f"./{out_file}_*")
    for file in list_output_files:
        shutil.move(file,out_models_dir)


def join_genes(in_dir, out_genes_dir):
    gene_files = glob.glob(f"{in_dir}1_Preprocessing/Genes_nt/*.fnn")
    #
    all_genes_file = f"{out_genes_dir}/all_genes.fnn"
    with open(all_genes_file, "wb") as out_file:
        for f in tqdm(gene_files):
            with open(f, "rb") as in_file:
                out_file.write(in_file.read())
    return(all_genes_file)


def join_prots(in_dir, out_genes_dir):
    gene_files = glob.glob(f"{in_dir}1_Preprocessing/Genes_aa/*.faa")
    #
    all_genes_file = f"{out_genes_dir}/all_genes.faa"
    with open(all_genes_file, "wb") as out_file:
        for f in tqdm(gene_files):
            with open(f, "rb") as in_file:
                out_file.write(in_file.read())
    return(all_genes_file)


def mmseqs2(gene_file, out_genes_dir, threads):
    os.mkdir("./mmseqs_tmp")
    subprocess.call(f'mmseqs createdb {gene_file} ./mmseqs_tmp/gene_db -v 1', shell=True)
    subprocess.call(f'mmseqs cluster ./mmseqs_tmp/gene_db ./mmseqs_tmp/DBclustered tempfiles --min-seq-id 0.90 --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign --threads {threads} -v 1', shell=True)
    subprocess.call(f'mmseqs createtsv ./mmseqs_tmp/gene_db ./mmseqs_tmp/gene_db ./mmseqs_tmp/DBclustered {out_genes_dir}/all_genes_CDS_mmseqs_clusters.tsv --threads {threads} -v 1', shell=True)
    subprocess.call(f'mmseqs createsubdb ./mmseqs_tmp/DBclustered ./mmseqs_tmp/gene_db ./mmseqs_tmp/my_rep_seqs -v 1', shell=True)
    rep_gene_file = f'{out_genes_dir}/my_rep_seqs.fnn'
    subprocess.call(f'mmseqs convert2fasta mmseqs_tmp/my_rep_seqs {rep_gene_file} -v 1', shell=True)
    shutil.rmtree('mmseqs_tmp')

    return(rep_gene_file)


def make_gene_plots(out_genes_dir, species):
    # Create binary matrix (genomes vs genes)
    subprocess.call(f'03a_MMSeqsTSV-to-BinaryMatrix.py -i {out_genes_dir}/all_genes_CDS_mmseqs_clusters.tsv -o {out_genes_dir}/pangenome_matrix.tsv', shell=True)
    
    # Plot pangenome plot
    subprocess.call(f'03b_Pangenome_Calculate_Model_Plot.py -b {out_genes_dir}/pangenome_matrix.tsv -o {out_genes_dir}/pangenome_model -n {species}', shell=True)
    
    # Plot heatmap genomes vs genes
    subprocess.call(f'03c_Clustermap_fromBinary.py -b {out_genes_dir}/pangenome_matrix.tsv -o {out_genes_dir}/pangenome_clustermap.pdf', shell=True)


def annotate(in_dir, out_genes_dir, rep_gene_file, tool, threads, species):
    # Join all protein files into a single file
    prot_file = join_prots(in_dir, out_genes_dir)

    # Get representative protein sequences using representative gene sequences file
    subprocess.call(f'03d_get_AA_reps_fasta.py -r {rep_gene_file} -a {prot_file} -o {out_genes_dir}/my_rep_seqs.faa', shell=True)
    
    # Make directory to output results and annotate using the appropiate tool
    if tool == 'eggnog':
        os.mkdir(f'{out_genes_dir}/EggNog')
        subprocess.call(f'emapper.py -i {out_genes_dir}/my_rep_seqs.faa -o {species} --output_dir {out_genes_dir}/EggNog/ --cpu {threads} > /dev/null 2>&1', shell=True)
        shutil.rmtree('./tempfiles')  
    elif tool == 'cogclassifier':
        os.mkdir(f'{out_genes_dir}/COGclassifier')
        subprocess.call(f'COGclassifier -i {out_genes_dir}/my_rep_seqs.faa -o {out_genes_dir}/COGclassifier/ -t {threads} > /dev/null 2>&1', shell=True)
        shutil.rmtree('./tempfiles')
    else:
        print(f'{tool} is not recognised as a valid annotation tool. Please provide one of: eggnog or cogclassifier')
    

def run_genes(in_dir, tool, threads, species):
    # Create directory inside input folder
    out_genes_dir = f"{in_dir}/3_gene_analyses"
    if os.path.isdir(out_genes_dir) == False:
        os.mkdir(out_genes_dir)
    
    # Join genes in one single file
    gene_file = join_genes(in_dir, out_genes_dir)

    # Cluster genes
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Clustering genes with MMseqs2...\n\n"]))
    rep_gene_file = mmseqs2(gene_file, out_genes_dir, threads)

    # Plot pangenome and clustermap
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Plotting pangenome and clustermap...\n\n"]))
    _ = make_gene_plots(out_genes_dir, species)

    # Annotate genes with eggnog or cogclassifier
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Annotating genes with {0}...\n\n".format(tool)]))
    _ = annotate(in_dir, out_genes_dir, rep_gene_file, tool, threads, species)

    # Assign pangenome class to genes
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Assigning pangenome class to genes...\n\n".format(tool)]))
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 03e_Get_Genes_Clusters_PanCat.py...\n\n"]))
    pan_matrix = f'{out_genes_dir}/pangenome_matrix.tsv'
    subprocess.call(f'03e_Get_Genes_Clusters_PanCat.py -b {pan_matrix} -m {out_genes_dir}/all_genes_CDS_mmseqs_clusters.tsv -r {in_dir}1_Preprocessing/RBM_results/RBMs_allV.rbm -o {out_genes_dir}/pancat_file.tsv', shell=True)


def run_one_vs_one(in_dir, rbm, pancat, anot, cds_gA, cds_gB, gA, gB, out_prefix, rec, subdivisions):
    # Create directory inside input folder if that variable is passed
    if in_dir != None:
        out_one_vs_one_dir = f"{in_dir}/4_one-vs-one"
        if os.path.isdir(out_one_vs_one_dir) == False:
            os.mkdir(out_one_vs_one_dir)

    # Run analysis for the pair of genomes selected
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 03f_Recombinant_pair_analysis.py...\n\n"]))
    subprocess.call(f'03f_Recombinant_pair_analysis.py -rbm {rbm} -pc {pancat} -ano {anot} -cA {cds_gA} -cB {cds_gB} -gA {gA} -gB {gB} -o {out_prefix} -rec {rec} -subs {subdivisions}', shell=True)
    
    # Move the output files to the directory if in_dir is provided, otherwise leave it in current folder
    if in_dir != None:
        list_out_files = glob.glob(f"./{out_prefix}*")
        for outfile in list_out_files:
            if os.path.isfile(outfile) == True:
                shutil.copy(outfile, out_one_vs_one_dir)
                os.remove(outfile)


def get_default_clade_list(rbm, out_prefix):
    table = pd.read_table(rbm, usecols=[0,1], header=None)
    merged = pd.concat([table[0], table[1]]).drop_duplicates()
    genomes = merged.str.split('_', expand=True)[0].drop_duplicates()
    df = pd.DataFrame(genomes).assign(clade = 'Clade')
    outfile = f'{out_prefix}_clade_list.txt'
    df.to_csv(outfile, sep = ',', header = False, index = False)
    return(outfile)


def get_default_color_list(clade_list, out_prefix):
    df_clades = pd.read_table(clade_list, header=None, sep=",", usecols=[1])
    df_clades_unique = df_clades.drop_duplicates()
    num_clades = len(df_clades_unique[1])
    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(num_clades)]
    
    df_clade_colors = df_clades_unique.assign(ran_colors = colors)
    
    outfile = f'{out_prefix}_color_clade_list.txt'
    df_clade_colors.to_csv(outfile, sep = ',', header = False, index = False)
    return(outfile)


def run_one_vs_many(in_dir, rbm, pancat, anot, input_list, out_prefix, rec, clade_list, color_list, species):
    # Create directory inside input folder if that variable is passed
    if in_dir != None:
        out_one_vs_many_dir = f"{in_dir}/5_one-vs-many"
        if os.path.isdir(out_one_vs_many_dir) == False:
            os.mkdir(out_one_vs_many_dir)
    else:
        out_one_vs_many_dir = '.'

    # Check if clade_list was provided. If not, create a new one considering all genomes in the same clade
    if clade_list == None:
        clade_list = get_default_clade_list(rbm, out_prefix)

    # Check if color_clade_list was provided. If not create a new one with random colors
    if color_list == None:
        color_list = get_default_color_list(clade_list, out_prefix)

    # Run analysis for one reference genome against many
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 03g_Recombinant_group_analysis.py...\n\n"]))
    subprocess.call(f'03g_Recombinant_group_analysis.py -i {input_list} -o {out_prefix} -rbm {rbm} -pc {pancat} -ano {anot}  -rec {rec}', shell=True)

    # Plot recombinant RBM curve
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 03b_Pangenome_Calculate_Model_Plot.py...\n\n"]))
    subprocess.call(f'03b_Pangenome_Calculate_Model_Plot.py -b {out_prefix}_rbm_matrix.tsv -o {out_one_vs_many_dir}/{out_prefix} -n {species}', shell=True)

    # Plot Recombinant gene clustermap
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 03c_Clustermap_fromBinary.py...\n\n"]))
    subprocess.call(f'03c_Clustermap_fromBinary.py -b {out_prefix}_rbm_matrix.tsv -o {out_one_vs_many_dir}/{out_prefix}_rbmclustermap.pdf', shell=True)

    # Plot recombinant rarefaction plot
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 03i_RBM-Clade_Rarefaction.py...\n\n"]))
    subprocess.call(f'03i_RBM-Clade_Rarefaction.py -r {out_prefix}_rbm_matrix.tsv -l {clade_list} -c {color_list} -o {out_prefix}_rbm_rarefaction', shell=True)

    # Move the output files to the directory if in_dir is provided, otherwise leave it in current folder
    if in_dir != None:
        list_out_files = glob.glob(f"./{out_prefix}*")
        for outfile in list_out_files:
            if os.path.isfile(outfile) == True:
                shutil.copy(outfile, out_one_vs_many_dir)
                os.remove(outfile)


def run_all_vs_all(in_dir, rbm, md, gene_list, out_prefix):
    # Create directory inside input folder if that variable is passed
    if in_dir != None:
        out_all_vs_all_dir = f"{in_dir}/6_all-vs-all"
        if os.path.isdir(out_all_vs_all_dir) == False:
            os.mkdir(out_all_vs_all_dir)
    else:
        out_all_vs_all_dir = '.'

    # Calculate fraction of recombinant genes and make a violin plot summarizing the results
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Calculating rarefaction curves and building violin plots for the fraction of identical genes with script 04a_Allv_RBM_Violinplot.py...\n\n"]))
    subprocess.call(f'04a_Allv_RBM_Violinplot_v2.py -rbm {rbm} -md {md} -gl {gene_list} -o {out_prefix}', shell=True)

    # Calculate PmRm ratio
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Calculating recombination and mutation ratios with script 04b_pmrm_analyses_v2.py...\n\n"]))
    subprocess.call(f'04b_pmrm_analyses_v2.py -md {md} -gl {gene_list} -rbm {rbm} -o {out_prefix}', shell=True)

    # Build PmRm vs ANI plot
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Plotting recombination to mutation ratios (ρ/θ)...\n\n"]))
    table_pmrm = pd.read_table(f'{out_prefix}_pmpr_pairwise_data.tsv', header=0, sep='\t')
    max_pmrm = round(table_pmrm['ρ/θ'].max()) + 1
    y_step_size = max_pmrm / 5
    min_ani =  round(table_pmrm['ani'].min()) - 1 if table_pmrm['ani'].min() < 95 else 95
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 04c_pmrm-ani_pairwise_plot.py...\n\n"]))
    subprocess.call(f'04c_pmrm-ani_pairwise_plot.py -i {out_prefix}_pmpr_pairwise_data.tsv -o {out_prefix} -ymax {max_pmrm} -yt {y_step_size} -xmin {min_ani}', shell=True)

    # Build PmRm vs F100 plot
    print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running script 04d_pmrm-f100_pairwise_plot.py...\n\n"]))
    subprocess.call(f'04d_pmrm-f100_pairwise_plot.py -i {out_prefix}_pmpr_pairwise_data.tsv -o {out_prefix} -ymax {max_pmrm} -yt {y_step_size}', shell=True)

    if in_dir != None:
        list_out_files = glob.glob(f"./{out_prefix}*")
        for outfile in list_out_files:
            if os.path.isfile(outfile) == True:
                shutil.copy(outfile, out_all_vs_all_dir)
                os.remove(outfile)
        if gene_list == 'gene_list.txt':
            shutil.copy(gene_list, out_all_vs_all_dir)
            os.remove(gene_list)


def options(action):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''	''')

    if action == "preprocessing":
        parser.description = "F100 preprocessing"
        parser.add_argument(
            '-i', '--in_dir', 
            default = None, 
            help = 'Directory containing genomes', 
            required = True
            )
        parser.add_argument(
            '-o', '--out_dir', 
            default = None, 
            help = 'Output directory', 
            required = True
            )
        parser.add_argument(
            '-t', '--threads', 
            type = int, 
            default = 1, 
            help = 'Number of threads. Default: 1'
            )
        parser.add_argument(
            '-r', '--rec_threshold', 
            type = float, 
            default = 99.8, 
            help = 'Minimum gene identity for recombining genes. Default: 99.8'
            )

    if action == "models":
        parser.description = "F100 comparison against pre-built models"
        parser.add_argument(
            '-i', '--in_dir',
            help='Please specify the folder with the results from the preprocessing step!',
            #metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-o', '--output_file_prefix',
            help='Please specify the output file prefix!',
            #metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-t', '--plot_title',
            help='OPTIONAL: Please specify the plot title (defalut: F100 vs ANI)!',
            #metavar='',
            type=str,
            nargs='+',
            default=['F100', 'vs', 'ANI'],
            required=False
            )
        parser.add_argument(
            '-xmin', '--xaxis_minimum',
            help='OPTIONAL: Minimum value to plot on x-axis. (Default=95.0)',
            #metavar='',
            type=float,
            default=95.0,
            required=False
            )
        parser.add_argument(
            '-xmax', '--xaxis_maximum',
            help='OPTIONAL: Maximum value to plot on x-axis. (Default=100.0)',
            #metavar='',
            type=float,
            default=100.0,
            required=False
            )
        parser.add_argument(
            '-ymin', '--yaxis_minimum',
            help='OPTIONAL: Minimum value to plot on y-axis. (Default=-8.0)',
            #metavar='',
            type=float,
            default=-8.0,
            required=False
            )
        parser.add_argument(
            '-ymax', '--yaxis_maximum',
            help='OPTIONAL: Maximum value to plot on y-axis. (Default=0.0)',
            #metavar='',
            type=float,
            default=0.0,
            required=False
            )
        parser.add_argument(
            '-mfi', '--model_file',
            help='Path to the pre-built model of interest',
            #metavar='',
            type=str,
            required=False,
            default=None
            )
        parser.add_argument(
            '-mfo', '--prebuilt_model_folder',
            help='Folder containing all pre-built models to compare against in .zip or .tsv format',
            #metavar='',
            type=str,
            required=False,
            default=None
            )
        parser.add_argument(
            '-s', '--xaxis_step_size',
            help='OPTIONAL: X-axis ticks step increment. (Default=1.0)',
            #metavar='',
            type=float,
            default=1.0,
            required=False
            )
        parser.add_argument(
            '-p', '--point_size',
            help='OPTIONAL: Size for i2 plotted points (Default=4.0)',
            #metavar='',
            type=float,
            default=4.0,
            required=False
            )
        parser.add_argument(
            '-a', '--point_alpha',
            help='OPTIONAL: Alpha value for i2 plotted points (Default=0.10)',
            #metavar='',
            type=float,
            default=0.10,
            required=False
            )

    if action == "gene-analysis":
        parser.description = "F100 gene analyses, including clustering, annotation and plotting"
        parser.add_argument(
            '-i', '--in_dir',
            default = None, 
            help = 'Please specify the folder with the results from the preprocessing steps', 
            required = True
            )  
        parser.add_argument(
            '-tool', '--tool', 
            default = 'eggnog', 
            type = str, 
            help = 'One of eggnog or cogclassifier. Default: eggnog'
            )
        parser.add_argument(
            '-sp', '--sp_name', 
            default = None, 
            help = 'Species name'
            )
        parser.add_argument(
            '-t', '--threads', 
            type = int, 
            default = 1, 
            help = 'Number of threads. Default: 1'
            )

    if action == "one-vs-one":
        parser.description = "Genome pair analyses (one vs one). Includes plots for recombinant genes distribution, length of recombined genes, funcional comparison of recombinant and non-recombinant genes, gene identity plots and estimation of recombination to mutation ratios"
        parser.add_argument(
            '-in', '--in_dir',
            help='Please specify the folder with the results from the preprocessing steps',
            metavar='',
            type=str,
            default = None,
            required=False
            )
        parser.add_argument(
            '-rbm', '--RBM_allvall_file',
            help='Please specify the all vs all RBM file! Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-pc', '--pangenome_categories',
            help='Please specify the tsv from 04d_Get_Genes_Clusters_PanCat.py! Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-ano', '--annotation_file',
            help='Specify the representative gene annotation file. Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-cA', '--input_CDS_A',
            help='Please specify the first prodigal CDS in fasta format!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-cB', '--input_CDS_B',
            help='Please specify the second prodigal CDS in fasta format!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-gA', '--input_genome_A',
            help='Please specify the first genome file in fasta format!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-gB', '--input_genome_B',
            help='Please specify the second genome file in fasta format!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-o', '--output_file_prefix',
            help='Please specify a prefix for the output file names!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-rec', '--recombination_cutoff',
            help='(OPTIONAL) Specify recombination cutoff (default = 99.8).',
            metavar='',
            type=float,
            default=99.8,
            required=False
            )
        parser.add_argument(
            '-subs', '--lineplot_subdivisions',
            help='(OPTIONAL)Specify subdivisions for lineplot (default = 10).',
            metavar='',
            type=int,
            default=10,
            required=False
            )

    if action == "one-vs-many":
        parser.description = "Analyses involving one reference genome against many query genomes. Includes plots for recombinant genes distribution, length of recombined genes, funcional comparison of recombinant and non-recombinant genes, gene identity plots and estimation of recombination to mutation ratios"
        parser.add_argument(
            '-in', '--in_dir',
            help='Please specify the folder with the results from the previous steps',
            metavar='',
            type=str,
            default = None,
            required=False
            )
        parser.add_argument(
            '-i', '--input_list',
            help='Please specify the input list (2 column tsv file)!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-o', '--output_file_prefix',
            help='Please specify a prefix for the output file names! Required!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-sp', '--species',
            help='Provide the name of the species. Required!',
            metavar='',
            type=str,
            default=None,
            required=True
            )
        parser.add_argument(
            '-rbm', '--RBM_allvall_file',
            help='Please specify the all vs all RBM file! Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-pc', '--pangenome_categories',
            help='Please specify the tsv from 04d_Get_Genes_Clusters_PanCat.py! Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-ano', '--annotation_file',
            help='Specify the representative gene annotation file. Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-rec', '--recombination_cutoff',
            help='(OPTIONAL) Specify recombination cutoff (default = 99.8).',
            metavar='',
            type=float,
            default=99.8,
            required=False
            )
        parser.add_argument(
            '-cld', '--clade_list',
            help='(OPTIONAL) Specify clade list. Two-columns comma-separated table with genome name and clade. If not provided all genomes will be considered to belong to the same clade.',
            metavar='',
            type=str,
            default=None,
            required=False
            )
        parser.add_argument(
            '-col', '--color_clade_list',
            help='(OPTIONAL) Specify colors for each clade in --clade_list. Two-columns comma-separated table with clade name and hex color code. If not provided random colors will be assigned to each clade.',
            metavar='',
            type=str,
            default=None,
            required=False
            )

    if action == "all-vs-all":
        parser.description = "All vs all analyses. Includes violin plots showing the fraction of recombining genes per genomovar, phylogroup or species as well as calculations and plots of the Purged Mutations (rho) vs Recent Mutations (theta) ratios"
        parser.add_argument(
            '-in', '--in_dir',
            help='Please specify the folder with the results from the previous steps',
            metavar='',
            type=str,
            default = None,
            required=False
            )
        parser.add_argument(
            '-md', '--metadata',
            help='Metadata tsv file with columns Genome, Genomovar, phylogroup, species',
            metavar='',
            type=str,
            default=None,
            required=True
            )
        parser.add_argument(
            '-o', '--output_file_prefix',
            help='Please specify a prefix for the output file names!',
            metavar='',
            type=str,
            required=True
            )
        parser.add_argument(
            '-rbm', '--RBM_allvall_file',
            help='Please specify the all vs all RBM file! Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )
        parser.add_argument(
            '-gl', '--gene_list',
            help='Please specify the file containing all gene ids! Required if input_folder is not provided',
            metavar='',
            type=str,
            required='--in_dir' not in sys.argv and '-in' not in sys.argv
            )

    args, unknown = parser.parse_known_args()
    return(parser, args)


def main():

    list_actions = ["help", "--help", "preprocessing", "models", "gene-analysis", "one-vs-one", "one-vs-many", "all-vs-all"]

    if len(sys.argv) < 2:
        print("Please provide one action of:\n\t", "\n\t".join(list_actions[2:len(list_actions)]))
        sys.exit()
    
    action = sys.argv[1]

    if action not in list_actions:
        print(''.join(["Action '", str(action), "' not recognized.\n"]))
        print("Please provide one action of:\n\t", "\n\t".join(list_actions[2:len(list_actions)]))
        sys.exit()
    
    if action == "help" or action == "--help":
        print("Please provide one action of:\n\t", "\n\t".join(list_actions[2:len(list_actions)]))
        sys.exit()

    parser, args = options(action)

    if action == "preprocessing":
        in_dir = args.in_dir
        out_dir = args.out_dir
        threads = args.threads
        rec_threshold = args.rec_threshold

        print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Starting the preprocessing step...\n\n"]))
        run_preprocessing(in_dir, out_dir, threads, rec_threshold)


    if action == "models":
        in_folder = args.in_dir
        in_file = glob.glob(f"{in_folder}/1_Preprocessing/RBM_results/*_F100.tsv")[0]
        out_file = args.output_file_prefix
        plot_title = args.plot_title
        xaxis_min = args.xaxis_minimum
        xaxis_max = args.xaxis_maximum
        yaxis_min = args.yaxis_minimum
        yaxis_max = args.yaxis_maximum
        model_file = args.model_file
        prebuilt_folder = args.prebuilt_model_folder
        xaxis_step_size = args.xaxis_step_size
        point_size = args.point_size
        point_alpha = args.point_alpha

        print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Running 02d_f100_scatter_pyGAM.py script to build plots...\n\n"]))

        run_models(in_folder, in_file, out_file, plot_title, xaxis_min, xaxis_max, yaxis_min, yaxis_max, model_file, prebuilt_folder, xaxis_step_size, point_size, point_alpha)

    if action == "gene-analysis":
        in_dir = args.in_dir
        tool = args.tool
        threads = args.threads
        if args.sp_name is not None:
            species = args.sp_name
        else:
            species = args.in_dir

        print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Starting gene analysis...\n\n"]))

        run_genes(in_dir, tool, threads, species)

    if action == "one-vs-one":
        if args.in_dir != None:
            in_dir = args.in_dir
            # Find RBM file from preprocessing step (1)
            rbm = f'{in_dir}/1_Preprocessing/RBM_results/RBMs_allV.rbm'
            if os.path.isfile(rbm) == False:
                print(f'RBM file not found in {rbm}, please check your files and provide a valid input folder. Alternatively provide the full path to the RBM file with the argument --RBM_from_AAI')
                sys.exit()
            # Find pancat file from gene analysis step (3)
            pancat = f'{in_dir}/3_gene_analyses/pancat_file.tsv'
            if os.path.isfile(pancat) == False:
                print(f'Pancat file not found in {pancat}, please check your files and provide a valid input folder. Alternatively provide the full path to the pancat file with the argument --pangenome_categories')
                sys.exit()
            # Find annotation file from gene analysis step (3)
            anot_dir_eggnog = f'{in_dir}/3_gene_analyses/EggNog'
            if os.path.isdir(anot_dir_eggnog) == True:
                anot = glob.glob(f'{anot_dir_eggnog}/*emapper.annotations')[0]
                if os.path.isfile(anot) == False:
                    print(f'Annotation file not found in {anot}, please check your files and provide a valid input folder. Alternatively provide the full path to the annotation file with the argument --annotation_file')
                    sys.exit()
            else:
                anot_dir_cog = f'{in_dir}/3_gene_analyses/COGclassifier'
                if os.path.isdir(anot_dir_cog) == True:
                    anot = glob.glob(f'{anot_dir_cog}/classifier_result.tsv')[0] 
                    if os.path.isfile(anot) == False:
                        print(f'Annotation file not found in {anot}, please check your files and provide a valid input folder. Alternatively provide the full path to the annotation file with the argument --annotation_file')
                        sys.exit()
                else:
                    print(f'Annotation folders "EggNog" or "COGclassifier" not found in {in_dir}/3_gene_analyses/. Please, be sure that the gene-analysis action finished successfully before running the one-vs-one action')
                    sys.exit()
        else:
            in_dir = ''
            if args.RBM_allvall_file is not None:
                rbm = args.RBM_allvall_file
            else:
                print("\nArguments --rbm, --pangenome_categories and --annotation_file are REQUIRED when --in_dir is not provided\n")
                sys.exit()
            if args.pangenome_categories is not None:
                pancat = args.pangenome_categories
            else:
                print("\nArguments --rbm, --pangenome_categories and --annotation_file are REQUIRED when --in_dir is not provided\n")
                sys.exit()
            if args.annotation_file is not None:
                anot = args.annotation_file
            else:
                print("\nArguments --rbm, --pangenome_categories and --annotation_file are REQUIRED when --in_dir is not provided\n")
                sys.exit()
        cds_gA = args.input_CDS_A
        cds_gB = args.input_CDS_B
        gA = args.input_genome_A
        gB = args.input_genome_B
        out_prefix = args.output_file_prefix
        rec = args.recombination_cutoff
        subdivisions = args.lineplot_subdivisions

        print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Starting one vs one analysis with script 03f_Recombinant_pair_analysis.py...\n\n"]))
        run_one_vs_one(in_dir, rbm, pancat, anot, cds_gA, cds_gB, gA, gB, out_prefix, rec, subdivisions)

    if action == "one-vs-many":
        if args.in_dir != None:
            in_dir = args.in_dir
            # Find RBM file from preprocessing step (1)
            rbm = f'{in_dir}/1_Preprocessing/RBM_results/RBMs_allV.rbm'
            if os.path.isfile(rbm) == False:
                print(f'RBM file not found in {rbm}, please check your files and provide a valid input folder. Alternatively provide the full path to the RBM file with the argument --RBM_allvall_file')
                sys.exit()
            # Find pancat file from gene analysis step (3)
            pancat = f'{in_dir}/3_gene_analyses/pancat_file.tsv'
            if os.path.isfile(pancat) == False:
                print(f'Pancat file not found in {pancat}, please check your files and provide a valid input folder. Alternatively provide the full path to the pancat file with the argument --pangenome_categories')
                sys.exit()
            # Find annotation file from gene analysis step (3)
            anot_dir_eggnog = f'{in_dir}/3_gene_analyses/EggNog'
            if os.path.isdir(anot_dir_eggnog) == True:
                anot = glob.glob(f'{anot_dir_eggnog}/*emapper.annotations')[0]
                if os.path.isfile(anot) == False:
                    print(f'Annotation file not found in {anot}, please check your files and provide a valid input folder. Alternatively provide the full path to the annotation file with the argument --annotation_file')
                    sys.exit()
            else:
                anot_dir_cog = f'{in_dir}/3_gene_analyses/COGclassifier'
                if os.path.isdir(anot_dir_cog) == True:
                    anot = glob.glob(f'{anot_dir_cog}/classifier_result.tsv')[0]
                    if os.path.isfile(anot) == False:
                        print(f'Annotation file not found in {anot}, please check your files and provide a valid input folder. Alternatively provide the full path to the annotation file with the argument --annotation_file')
                        sys.exit()
                else:
                    print(f'Annotation folders "EggNog" or "COGclassifier" not found in {in_dir}/3_gene_analyses/. Please, be sure that the gene-analysis action finished successfully before running the one-vs-one action')
                    sys.exit()
        else:
            in_dir = ''
            if args.RBM_allvall_file is not None:
                rbm = args.RBM_allvall_file
            else:
                print("\nArguments --RBM_allvall_file, --pangenome_categories and --annotation_file are REQUIRED when --in_dir is not provided\n")
                sys.exit()
            if args.pangenome_categories is not None:
                pancat = args.pangenome_categories
            else:
                print("\nArguments --RBM_allvall_file, --pangenome_categories and --annotation_file are REQUIRED when --in_dir is not provided\n")
                sys.exit()
            if args.annotation_file is not None:
                anot = args.annotation_file
            else:
                print("\nArguments --RBM_allvall_file, --pangenome_categories and --annotation_file are REQUIRED when --in_dir is not provided\n")
                sys.exit()
        input_list = args.input_list
        out_prefix = args.output_file_prefix
        rec = args.recombination_cutoff
        clade_list = args.clade_list
        color_list = args.color_clade_list
        species = args.species

        print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Starting one vs many analysis with script 03g_Recombinant_group_analysis.py...\n\n"]))
        run_one_vs_many(in_dir, rbm, pancat, anot, input_list, out_prefix, rec, clade_list, color_list, species)

    if action == 'all-vs-all':
        if args.in_dir != None:
            in_dir = args.in_dir
            # Find RBM file from preprocessing step (1)
            rbm = f'{in_dir}/1_Preprocessing/RBM_results/RBMs_allV.rbm'
            if os.path.isfile(rbm) == False:
                print(f'RBM file not found in {rbm}, please check your files and provide a valid input folder. Alternatively provide the full path to the RBM file with the argument --RBM_allvall_file')
                sys.exit()
            # Make list of genes
            all_genes = f'{in_dir}/3_gene_analyses/all_genes.fnn'
            if os.path.isfile(all_genes) == True:
                gene_list = 'gene_list.txt'
                subprocess.call(f'grep -i ">" {in_dir}/3_gene_analyses/all_genes.fnn > {gene_list}', shell=True)
            elif os.path.isfile(all_genes) == False:
                list_genes = glob.glob(f'{in_dir}/1_Preprocessing/Genes_nt/*.fnn')
                if len(list_genes) > 0:
                    ids = []
                    for file in list_genes:
                        with open(file, 'r') as f:
                            for line in f:
                                if line.startswith('>'):
                                    ids.append(line.rstrip('\n'))
                    gene_list = 'gene_list.txt'
                    out_gl = open(gene_list,'w')
                    out_gl.writelines(ids)
                    out_gl.close()
                else:
                    print(f'Gene files were not found in {in_dir}/1_Preprocessing/Genes_nt/, please check your files and provide a valid input folder. Alternatively provide the full path to the gene file with the argument --gene_list')
                    sys.exit()
        else:
            in_dir = ''
            if args.RBM_allvall_file is not None:
                rbm = args.RBM_allvall_file
            else:
                print("\nArguments --RBM_allvall_file and --gene_list are REQUIRED when --in_dir is not provided\n")
                sys.exit()
            if args.gene_list is not None:
                gene_list = args.gene_list
            else:
                print("\nArguments --RBM_allvall_file and --gene_list are REQUIRED when --in_dir is not provided\n")
                sys.exit()
        out_prefix = args.output_file_prefix
        md = args.metadata

        print(''.join(["\n\n[", datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "] Starting all vs all analyses...\n\n"]))
        run_all_vs_all(in_dir, rbm, md, gene_list, out_prefix)

if __name__ == "__main__":
    main()
