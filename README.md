### F100_Prok_Recomb_auto
This repository contains the information to run the F100 pipeline to estimate recombination between genomes developed by Roth Conrad (https://github.com/rotheconrad/F100_Prok_Recombination/) calling a single script.

The F100_main.py compile the scripts in the original F100_Prok_Recombination repository into 7 different workflows to make it easier to use. It might be specially useful for those trying the pipeline for the first time or with little experience. However, the users looking for more flexibility might want to run the pipeline script by script, as originally developed.


## Installation

1. Create conda environment and install dependencies

```
conda create -n F100 -c bioconda -c conda-forge prodigal blast diamond mmseqs2 eggnog-mapper cogclassifier python tqdm pandas matplotlib seaborn scipy pygam datashader lmfit
```

```
conda activate F100
```

```
pip install tqdm pyrodigal
```

2. EggNog and/or COGclassifier should be installed (at least one of them):

Install COGclassifier

```
pip install cogclassifier
```

Download EggNog databases

```
download_eggnog_data.py
```

If an error arises be sure that the folder where the script is trying to put the databases exists

3. Download the original F100 repository and add the scripts to the bin folder in the conda environment (important step! Otherwise the F100_main.py script won't work)

```
git clone https://github.com/rotheconrad/F100_Prok_Recombination.git
chmod +x F100_Prok_Recombination/00d_Workflow_Scripts/*.py
py_path=$(which python)
for i in F100_Prok_Recombination/00d_Workflow_Scripts/*.py ; do sed -i "1s?.*?\#\!$py_path?" $i ; done
env_bin_dir=${py_path%%python}
cp F100_Prok_Recombination/00d_Workflow_Scripts/*.py $env_bin_dir
```

4. Clone the repository or download the F100_main.py script

```
git clone https://github.com/baldeguer-riquelme/F100_Prok_Recomb_auto/
```


## Brief description

The pipeline consist of 7 workflows (or actions):

1. preprocessing

2. metadata

3. gene-analyses
              
4. models
              
5. one-vs-one
              
6. one-vs-many
              
7. all-vs-all


The first 3 workflows (preprocessing, metadata and gene-analyses) are key since they generate the outputs needed for the other workflows to run.

The list of scripts and output files produced by each workflow are detailed below.

| Workflow      | Scripts                                         | Output                                         | Explanation                                                                                                                              |
|---------------|-------------------------------------------------|------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------|
| preprocessing | 01a_rename_fasta_v2.py                          | Renamed files                                  | Replace '_' by '.'                                                                                                                       |
|               | Prodigal                                        | Predicted genes and proteins                   | Predicted genes and proteins                                                                                                             |
|               | 02b_get_RBMs_pipeline_auto.py                   | RBMs_allV.rbm                                  | Table with all vs all RBM data                                                                                                           |
|               | 02c_get_F100.py                                 | {out_prefix}_F100.tsv                          | Table with the fraction of shared genes (defined by the identity threshold)                                                              |
| metadata      | fastANI                                         | fastANI_allV.ani                               | All vs all ANI values                                                                                                                    |
|               | 01b_fastANI_scatter_pyGAM.py                    | fastANI_allV_sharedfrac_ANI.pdf                | Dot plot of ANI vs Percentage of shared genome. Each dot is a pair                                                                       |
|               | 01c_fastANI_clustermap_v2.py                    | {out_prefix}_clustermap.pdf                    | Heatmap showing the 100 - ANI values                                                                                                     |
|               |                                                 | {out_prefix}_Legend_clade-X.pdf                | Color legend for each clade in {out_prefix}_clustermap.pdf                                                                               |
|               |                                                 | {out_prefix}_Local_minimums.pdf                | Plot showing the ANI distribution and the minimum and maximum values used to identify clusters                                           |
|               |                                                 | {out_prefix}_predicted_clusters.tsv            | Tab-delimited table indicating the clades each genome belongs to                                                                         |
|               |                                                 | {out_prefix}_meta_colors.tsv                   | Two column tab-delimtied table with the hex color code assigned to each clade                                                            |
|               |                                                 | {out_prefix}_distmat.tsv                       | ANI matrix                                                                                                                               |
| gene-analysis | MMseqs2                                         | all_genes_CDS_mmseqs_clusters.tsv              | Mmseqs2 cluster database                                                                                                                 |
|               |                                                 | my_rep_seqs.fnn                                | Representative genes sequences                                                                                                           |
|               | 03a_MMSeqsTSV-to-BinaryMatrix.py                | pangenome_matrix.tsv                           | Binary matrix file. Columns are genome, rows are genes. 1 = present ; 0=absent                                                           |
|               | 03b_Pangenome_Calculate_Model_Plot.py           | pangenome_model_pangenome_curves.pdf           | Rarefaction curves with pangenome, core and specific genes                                                                               |
|               | 03b_Pangenome_Calculate_Model_Plot.py           | pangenome_model_pangenome_data.tsv             | Data used for the rarefaction plot above                                                                                                 |
|               | 03c_Clustermap_fromBinary.py                    | pangenome_clustermap.pdf                       | Heatmap of presence/absence of genes (rows) in genomes (columns). 1 = present ; 0=absent                                                 |
|               | 03d_get_AA_reps_fasta.py                        | my_rep_seqs.faa                                | Representative protein sequences                                                                                                         |
|               | emapper.py                                      | {out_prefix}.emapper.annotations               | Annotation of representative proteins based on EggNog                                                                                    |
|               | COGclassifier                                   | classifier_result.tsv                          | Annotation of representative proteins based on COGclassifier                                                                             |
|               | 03e_Get_Genes_Clusters_PanCat.py                | pancat_file.tsv                                | Table with information about the cluster, pangenome category, n/N and Distance for each gene                                             |
| models        | 02d_f100_scatter_pyGAM.py                       | {out_prefix}_{model}_GAMplot.pdf               | F100 vs ANI plot for the input data and the pre-built model                                                                              |
|               |                                                 | {out_prefix}_{model}_sig-pairs.tsv             | Pairs outside the 95% confidence interval. Above model: "Recombining" ; Below model: "Non-recombining"                                   |
| one-vs-one    | 03f_Recombinant_pair_analysis.py                | {out_prefix}_A_Adjecent_Rec_Length.pdf         | Histogram of the number of consecutive recombinant and non-recombinant genes in genome A                                                 |
|               |                                                 | {out_prefix}_A_distance.pdf                    | Distance between recombinant genes in genome A                                                                                           |
|               |                                                 | {out_prefix}_A_gene_table.tsv                  | Data (Identity, annotation, etc) for each gene of genome A                                                                               |
|               |                                                 | {out_prefix}_A_posline_out.pdf                 | Sequence identity vs position in the genome A                                                                                            |
|               |                                                 | {out_prefix}_A_Rec_rate_windows.pdf            | Sliding window recombination to mutation rate for genome A                                                                               |
|               |                                                 | {out_prefix}_A_rec_rate_table.tsv              | Data used for the sliding window recombination to mutation rate above                                                                    |
|               |                                                 | {out_prefix}_B_Adjecent_Rec_Length.pdf         | Histogram of the number of consecutive recombinant and non-recombinant genes in genome A                                                 |
|               |                                                 | {out_prefix}_B_distance.pdf                    | Distance between recombinant genes in genome A                                                                                           |
|               |                                                 | {out_prefix}_B_gene_table.tsv                  | Data (Identity, annotation, etc) for each gene of genome A                                                                               |
|               |                                                 | {out_prefix}_B_posline_out.pdf                 | Sequence identity vs position in the genome A                                                                                            |
|               |                                                 | {out_prefix}_B_Rec_rate_windows.pdf            | Sliding window recombination to mutation rate for genome A                                                                               |
|               |                                                 | {out_prefix}_B_rec_rate_table.tsv              | Data used for the sliding window recombination to mutation rate above                                                                    |
|               |                                                 | {out_prefix}_annotations_bar.pdf               | Recombinant vs non-recombinat function annotation                                                                                        |
|               |                                                 | {out_prefix}_genomes.pdf                       | Distribution along the genome of conserved, recombinant and non-recombinant genes                                                        |
| one-vs-many   | 03g_Recombinant_group_analysis.py               | {out_prefix}_Adjecent_Rec_Length.pdf           | Histogram of the number of consecutive recombinant and non-recombinant genes in reference genome                                         |
|               |                                                 | {out_prefix}_Adjecent_Rec_Length.tsv           | Data used for the histogram above                                                                                                        |
|               |                                                 | {out_prefix}_annotations_bar.pdf               | Recombinant vs non-recombinat function annotation                                                                                        |
|               |                                                 | {out_prefix}_gene_distribution.pdf             | Distance between recombinant genes in the reference genome                                                                               |
|               |                                                 | {out_prefix}_group_data.tsv                    | Data (Identity, annotation, etc) for each hit between query and reference genomes                                                        |
|               |                                                 | {out_prefix}_posbar.pdf                        | Recombinant vs non-recombinat function annotation                                                                                        |
|               |                                                 | {out_prefix}_posline.pdf                       | Sequence identity vs position in the reference genome                                                                                    |
|               |                                                 | {out_prefix}_rbm_matrix.tsv                    | Binary matrix file. Columns are genome, rows are genes. 1 = present ; 0=absent                                                           |
|               | 03b_Pangenome_Calculate_Model_Plot.py           | {out_prefix}_pangenome_curves.pdf              | Rarefaction curves with pangenome, core and specific genes                                                                               |
|               |                                                 | {out_prefix}_Pangenome_data.tsv                | Data used for the rarefaction plot above                                                                                                 |
|               | 03c_Clustermap_fromBinary.py                    | {out_prefix}_rbmclustermap.pdf                 | Heatmap of presence/absence of genes (rows) in genomes (columns). 1 = present ; 0=absent                                                 |
|               | workflow (if not provided by user, recommended) | {out_prefix}_clade_list.txt                    | Comma-separated tab with genome name and clade. By default, all genomes are considered in the same clade. Better if provided by the user |
|               |                                                 | {out_prefix}_color_clade_list.txt              | Comma-separated tab with clade and color in hex format. By defautl, random colors are selected. Better if provided by the user           |
|               | 03i_RBM-Clade_Rarefaction.py                    | {out_prefix}_rbm_rarefaction_01_plot.pdf       | Cumulative fraction of identical genes as new genomes are added                                                                          |
|               |                                                 | {out_prefix}_rbm_rarefaction_01_data.tsv       | Data to plot the cumulative fraction of identical genes as new genomes are added                                                         |
| all-vs-all    | 04a_Allv_RBM_Violinplot_v2.py                   | {out_prefix}-figure01-Same_Genomovar.pdf       | Dual violin plot for the cumulative fraction of identical genes within the SAME GENOMOVAR                                                |
|               |                                                 | {out_prefix}-figure02-Different_Genomovar.pdf  | Dual violin plot for the cumulative fraction of identical genes between DIFFERENT GENOMOVAR                                              |
|               |                                                 | {out_prefix}-figure03-Same_Phylogroup.pdf      | Dual violin plot for the cumulative fraction of identical genes within the SAME PHYLOGROUP                                               |
|               |                                                 | {out_prefix}-figure04-Different_Phylogroup.pdf | Dual violin plot for the cumulative fraction of identical genes between DIFFERENT PHYLOGROUP                                             |
|               |                                                 | {out_prefix}-figure05-Species.pdf              | Dual violin plot for the cumulative fraction of identical genes WITHIN and BETWEEN SPECIES                                               |
|               |                                                 | {out_prefix}-figure06-Summary.pdf              | Dual violin plot summarizing the cumulative fraction of identical genes in different groups of genomes                                   |
|               |                                                 | {out_prefix}_genomovar_data.tsv                | Data used for the genomovars plots above                                                                                                 |
|               |                                                 | {out_prefix}_phylogroup_data.tsv               | Data used for the phylogroup plots above                                                                                                 |
|               |                                                 | {out_prefix}_species_data.tsv                  | Data used for the species plots above                                                                                                    |
|               | 04b_pmrm_analyses_v2.py                         | {out_prefix}_pmpr_pairwise_data.tsv            | Table with the purged mutation (Pm) to recent mutation (Rm) ratio for each pairwise comparison                                           |
|               | 04c_pmrm-ani_pairwise_plot.py                   | {out_prefix}_pmrm_ani.pdf                      | Dot plot of ANI vs Pm/Rm ratio                                                                                                           |
|               | 04d_pmrm-f100_pairwise_plot.py                  | {out_prefix}_pmrm_f100.pdf                     | Dot plot of F100 vs Pm/Rm ratio                                                                                                          |


## Run the pipeline

For all workflows but preprocessing, if the --in_dir argument is provided and is valid, the F100_main.py script will search automatically for the needed files to run the workflow. However, the user can also specify a path to each required file. 


This is an example with the minimum requirements to run the different workflows. 

```
conda activate F100

# Run preprocessing
python F100_main.py preprocessing --in_dir genomes/ --out_dir Sruber -t 8


# Run models.
# The Pre-built_models can include the models Subsampled_Genome_Model_Data.tsv.zip and Simulated_Neutral_Model_Data.tsv.zip available in https://github.com/rotheconrad/F100_Prok_Recombination/ and/or those previously built by the user.
# Instead of a folder you can simply specify the path the one specific model with the -mfi argument 
python F100_main.py models --in_dir Sruber --output_file_prefix results -mfo Pre-built_models/


# Run gene-analysis
python F100_main.py gene-analysis --in_dir Sruber -sp Sruber -t 8 --tool eggnog


# Run one vs one
python F100_main.py one-vs-one --in_dir Sruber -cA Sruber/1_Preprocessing/Genes_nt/GCA_905335855.1_Sal._ruber_isolate_CM02_genomic.fnn -cB Sruber/1_Preprocessing/Genes_nt/GCA_905
335875.1_Sal._ruber_isolate_CM05_genomic.fnn -gA Sruber/1_Preprocessing/Genomes/GCA_905335855.1_Sal._ruber_isolate_CM02_genomic.fna -gB Sruber/1_Preprocessing/Genomes/GCA_905335875
.1_Sal._ruber_isolate_CM05_genomic.fna -o results


# Run one vs many
# input_list_one_vs_many is a two column tab separated file with genome fasta paths in first column and gene fasta paths in the second column. The first genome/gene file in the list becomes the main genome the others are compared to. An example can be found in this repository.
python F100_main.py one-vs-many --in_dir Sruber -i input_list_one_vs_many -o results -sp Sruber


# Run all vs all
# metadata.txt is tab separated file with columns Genome, Genomovar, phylogroup, species . An example can be found in this repository
python F100_main.py all-vs-all --in_dir Sruber -md metadata.txt -o results

```
