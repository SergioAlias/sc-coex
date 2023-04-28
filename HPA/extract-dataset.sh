#!/bin/bash
# Sergio Alías, 20220426
# Last modified 20230428

# Script para extraer los clsters de un tejido del dataset de HPA

# Argumentos

TISSUE=$1       # Tejido (entrecomillar)
LCLUST=$2       # Número del último cluster del tejido (considerando que el primero es 0)

# Ejecución

echo "Executing $0 ..."

TISSFN=$( echo $TISSUE | tr [[:upper:]] [[:lower:]] | tr " " "-")

for ((i=0;i<=LCLUST;i++)); do
	echo "Preparing cluster $i..."
    head -1 rna_single_cell_read_count.tsv > ./HPA/datasets/$TISSFN/$TISSFN-$i-dataset.tsv
    cat rna_single_cell_read_count.tsv | grep -P "$TISSUE\t[0-9]+\t$i\t" >> ./HPA/datasets/$TISSFN/$TISSFN-$i-dataset.tsv
done

echo "Finished! :)"