#!/bin/bash
# Sergio Alías, 20220928
# Last modified 20231129

# Workflow for the analysis of HPO-related genes COEX values

# What we need (examples using HPA dataset folder):

#   - COEX matrices obtained with COTAN (not automated yet) --> ./HPA/results/tissue/cluster/HPA.tissue-cluster.matrix.COTAN.RDS

#   - Script for converting COEX matrices to TSV (input for the next step) --> ./coex-to-tsv.R

#   - Script for subsetting the matrix --> ./subset-genes.py
#     We need the following annotation files:
#       - HPO-genes relationships --> ./annotations/HPO/hpo_to_genes_20220603.json
#       - Entrez-Ensembl gene names --> ./annotations/HPO/entrez_to_ensembl_20220603.json

#   - Script for performing the coex analysis --> ./plot_coex.R
#     We need the following annotation files:
#       - ./annotations/rna_single_cell_type_tissue.tsv (necessary for next step)
#       - ./annotations/cluster-annotation/HPA/tissue-cluster-annotation (generated in this script if it does not exist yet (HPA))
#       - ./annotations/HPO/hpo-id-to-name

#   - List of HPO terms we wish to analyse. We create it in this script with --> ./annotations/filtered_tis_hpo_ngenes_20220711.tsv

#   - Number of clusters of the tissue --> ./annotations/HPA-nclusters



############################
## Command line arguments ##
############################

TIS=$1                # Tissue name (as in ~/annotations/HPA-nclusters)
HPO_LIST=${2-all}     # HPO names file. Default (all) makes the script work for every HPO related to the tissue
DATASET=${3-HPA}      # Dataset used
ALT_NAME=${4-$TIS}    # Alternative tissue name. Example: blood (as in UBERON), pbmc (as in HPA)



#################
## Main Script ##
#################

### Getting the number of clusters of the tissue

NCLUS_FIL=./annotations/$DATASET-nclusters

NCLUS=$( cat $NCLUS_FIL | grep "$TIS" | cut -f 2 ) # We grep the tissue name with the first letter in mayus


### Listing the COEX matrix files for that tissue

MAT_DIR=./$DATASET/results/$TIS/

MATS=($( find $MAT_DIR | grep ".matrix." | sort -t "-" -k 2 -n )) # We get the COEX matrices paths sorted as an array


### Converting .RDS COEX matrices to TSV files

FINDTSV=($( find $MAT_DIR | grep ".tsv")) # Array of filenames ending with ".tsv"

if [ ${#FINDTSV[@]} -eq 0 ]; then # If TSV files does not exist (therefore $FINDTSV length is zero)
    Rscript ./coex-to-tsv.R -t $TIS -c $NCLUS -d $DATASET # Create them
fi


### Subsetting COEX matrices

if [ "$HPO_LIST" = "all" ]; then
    HPO_FILE=./annotations/filtered_tis_hpo_ngenes_$DATASET.tsv
    HPO_NAMES=($( cat $HPO_FILE | cut -f 1,2 | grep $TIS | cut -f 2 | tr -d ":" ))
    echo "HPO terms list not provided, default: all (${#HPO_NAMES[@]} terms)"
else
    HPO_NAMES=($(cat $HPO_LIST))
    echo "HPO terms list provided (${#HPO_NAMES[@]} terms)"
fi


for HPO in ${HPO_NAMES[@]}; do
    mkdir -p ./coex-analysis/$DATASET/$TIS/$HPO
    FINDSUBSET=($( find ./coex-analysis/$DATASET/$TIS/$HPO | grep "subset_genes"))
    if [ ${#FINDSUBSET[@]} -eq 0 ]; then
        python3.7 ./subset-genes.py $HPO $TIS $NCLUS $DATASET
    fi
done


### Get cluster annotations

mkdir -p ./annotations/cluster-annotation/$DATASET
FINDCLUANN=($( find ./annotations/cluster-annotation/$DATASET | grep $TIS))
if [ ${#FINDCLUANN[@]} -eq 0 ]; then
    echo "Creating $TIS-cluster-annotation..."
    cat ./annotations/rna_single_cell_type_tissue.tsv | grep $ALT_NAME | cut -f 3-5 | head -$NCLUS > ./annotations/cluster-annotation/$DATASET/$TIS-cluster-annotation
fi


### P-val, Wilcoxon test and plots

for HPO in ${HPO_NAMES[@]}; do
    FINDPLOTS=($( find ./coex-analysis/$DATASET/$TIS/$HPO | grep "_plots_"))
    if [ ${#FINDPLOTS[@]} -eq 0 ]; then
        Rscript ./plot_coex.R -t $TIS -c $((NCLUS-1)) -o $HPO -d $DATASET -r 1000
    fi
done
