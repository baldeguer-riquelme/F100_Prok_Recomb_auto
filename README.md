### F100_Prok_Recomb_auto
This repository contains the information to run the F100 pipeline developed by Roth Conrad (https://github.com/rotheconrad/F100_Prok_Recombination/tree/main) to estimate recombination between genomes calling a single script.

The F100_main.py compile the scripts in the original F100_Prok_Recombination repository to make it easier to use. It might be specially useful for those trying the pipeline for the first time or with little experience. However, the users looking for more flexibility might want to run the pipeline script by script, as originally developed.


## Manual installation

1. Create conda environment and install dependencies

```
conda create -n F100 -c bioconda -c conda-forge prodigal blast diamond=2.1.8 mmseqs2 eggnog-mapper cogclassifier python tqdm pandas matplotlib seaborn scipy pygam datashader lmfit
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
cp F100_Prok_Recombination/00d_Workflow_Scripts/* $env_bin_dir
```

## Automatic installation

Use the install_F100_dependencies.py script to install all dependencies automatically. (TO BE DONE)


## Run the pipeline

The pipeline consist of 6 workflows (or actions):

1. preprocessing
              
2. models

3. gene-analyses
              
4. one-vs-one
              
5. one-vs-many
              
6. all-vs-all


The preprocessing and gene-analyses are key since they produce outputs needed for the other workflows to run.
