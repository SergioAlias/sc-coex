#!/bin/bash
# Sergio Alías, 20221103
# Last modified 20221103

# Mini-script para correr COTAN-cluster-script.R para todos los clusters de un tissue

TIS=$1          # Tissue name
LAST_CLUS=$2    # Higher cluster number

for ((CLU=0;CLU<=$LAST_CLUS;CLU++)); do
    Rscript ../COTAN-cluster-script.R $TIS-$CLU
done

