### F100_Prok_Recomb_auto
This repository contains the information to run the F100 pipeline developed by Roth Conrad (https://github.com/rotheconrad/F100_Prok_Recombination/) to estimate recombination between genomes calling a single script.

The F100_main.py compile the scripts in the original F100_Prok_Recombination repository to make it easier to use. It might be specially useful for those trying the pipeline for the first time or with little experience. However, the users looking for more flexibility might want to run the pipeline script by script, as originally developed.


## Installation

1. Create conda environment and install dependencies

```
conda create -n F100 -c bioconda -c conda-forge prodigal blast diamond mmseqs2 eggnog-mapper cogclassifier python tqdm pandas matplotlib seaborn scipy pygam datashader lmfit
```

```
pip install tqdm pyrodigal cogclassifier
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

3. Download F100 pipeline scripts and add them to the bin folder for the conda environment (critical, otherwise the script won't work)

```
git clone https://github.com/rotheconrad/F100_Prok_Recombination.git
chmod +x F100_Prok_Recombination/00d_Workflow_Scripts/*.py
env_bin_dir=${PATH%%:*}
var=$(echo $env_bin_dir/python)
for i in F100_Prok_Recombination/00d_Workflow_Scripts/*.py ; do sed -i "1s?.*?\#\!$var?" $i ; done
cp F100_Prok_Recombination/00d_Workflow_Scripts/*.py $env_bin_dir
```



## Brief description

The pipeline consist of 6 workflows (or actions):

1. preprocessing
              
2. models

3. gene-analyses
              
4. one-vs-one
              
5. one-vs-many
              
6. all-vs-all


The preprocessing and gene-analyses are key since they produce outputs needed for the other workflows to run.




## Run the pipeline

This is an example with the minimum requirements to run the different workflows.

```
conda activate F100

# Run preprocessing
python F100_main.py preprocessing --in_dir genomes/ --out_dir Sruber -t 8


# Run models.
# The Pre-built_models can include the models Subsampled_Genome_Model_Data.tsv.zip and Simulated_Neutral_Model_Data.tsv.zip available in https://github.com/rotheconrad/F100_Prok_Recombination/ and/or those previously built by the user.
# Instead of a folder you can simply specify the path the one specific model
python F100_main.py models --in_dir Sruber --output_file_prefix results -m Pre-built_models/


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
