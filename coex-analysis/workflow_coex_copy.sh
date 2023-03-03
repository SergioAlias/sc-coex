#!/bin/bash
# Sergio Alías, 20220928
# Last modified 20221026

# Workflow for the analysis of HPO-related genes COEX values

# What we need:

#   - COEX matrices obtained with COTAN (not automated yet) --> ~/TFM/COTAN/results/tissue/cluster/HPA.tissue-cluster.matrix.COTAN.RDS

#   - Script for converting COEX matrices to TSV (input for the next step) --> ~/TFM/COTAN/results/coex-to-tsv.R

#   - Script for subsetting the matrix --> ~/TFM/coex-analysis/subset-genes.py
#     We need the following annotation files:
#       - HPO-genes relationships --> ~/TFM/annotations/HPO/hpo_to_genes_20220603.json
#       - Entrez-Ensembl gene names --> ~/TFM/annotations/HPO/entrez_to_ensembl_20220603.json

#   - Script for performing the coex analysis --> ~/TFM/coex-analysis/plot_coex.R
#     We need the following annotation files:
#       - ~/TFM/annotations/rna_single_cell_type_tissue.tsv (necessary for next step)
#       - ~/TFM/annotations/cluster-annotation/tissue-cluster-annotation (generated in this script if it does not exist yet)
#       - ~/TFM/annotations/HPO/hpo-id-to-name

#   - List of HPO terms we wish to analyse. We create it in this script with --> ~/TFM/annotations/filtered_tis_hpo_ngenes_20220711.tsv

#   - Number of clusters of the tissue --> ~/TFM/annotations/HPA-nclusters



############################
## Command line arguments ##
############################

TIS=$1                # Tissue name (as in ~/annotations/HPA-nclusters)
HPO_LIST=${2-all}     # HPO names file. Default (all) makes the script work for every HPO related to the tissue



#################
## Main Script ##
#################

### Getting the number of clusters of the tissue

NCLUS_FIL=~/TFM/annotations/HPA-nclusters

NCLUS=$( cat $NCLUS_FIL | grep "$TIS" | cut -f 2 ) # We grep the tissue name with the first letter in mayus


### Listing the COEX matrix files for that tissue

MAT_DIR=~/TFM/COTAN/results/$TIS/

MATS=($( find $MAT_DIR | grep ".matrix." | sort -t "-" -k 2 -n )) # We get the COEX matrices paths sorted as an array


### Converting .RDS COEX matrices to TSV files

FINDTSV=($( find $MAT_DIR | grep ".tsv")) # Array of filenames ending with ".tsv"

if [ ${#FINDTSV[@]} -eq 0 ]; then # If TSV files does not exist (therefore $FINDTSV length is zero)
    Rscript ~/TFM/COTAN/results/coex-to-tsv.R -t $TIS -c $NCLUS # Create them
fi


### Subsetting COEX matrices

if [ "$HPO_LIST" = "all" ]; then
    HPO_FILE=~/TFM/annotations/filtered_tis_hpo_ngenes_20220711.tsv
    HPO_NAMES=($( cat $HPO_FILE | grep $TIS | cut -f 2 | tr -d ":" ))
    echo "HPO terms list not provided, default: all (${#HPO_NAMES[@]} terms)"
else
    HPO_NAMES=($(cat $HPO_LIST))
    echo "HPO terms list provided (${#HPO_NAMES[@]} terms)"
fi


for HPO in ${HPO_NAMES[@]}; do
    mkdir -p ~/TFM/coex-analysis/$TIS/$HPO
    FINDSUBSET=($( find ~/TFM/coex-analysis/$TIS/$HPO | grep "subset_genes"))
    if [ ${#FINDSUBSET[@]} -eq 0 ]; then
        python3 ~/TFM/coex-analysis/subset-genes.py $HPO $TIS $NCLUS
    fi
done


### Get cluster annotations

mkdir -p ~/TFM/annotations/cluster-annotation
FINDCLUANN=($( find ~/TFM/annotations/cluster-annotation | grep $TIS))
if [ ${#FINDCLUANN[@]} -eq 0 ]; then
    echo "Creating $TIS-cluster-annotation..."
    cat ~/TFM/annotations/rna_single_cell_type_tissue.tsv | grep $TIS | cut -f 3-5 | head -$NCLUS > ~/TFM/annotations/cluster-annotation/$TIS-cluster-annotation
fi


### P-val, Wilcoxon test and plots

for HPO in ${HPO_NAMES[@]}; do
    FINDPLOTS=($( find ~/TFM/coex-analysis/$TIS/$HPO | grep "_plots_"))
    if [ ${#FINDPLOTS[@]} -eq 0 ]; then
        Rscript ~/TFM/coex-analysis/plot_coex.R -t $TIS -c $((NCLUS-1)) -o $HPO -r 1000
    fi
done
