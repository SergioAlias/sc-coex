#!/bin/bash
# Sergio Alías, 20220426
# Last modified 20230428

# Script para extraer los clsters de un tejido del dataset de HPA

# Argumentos

TISSUE=$1       # Tejido (entrecomillar si es más de una palabra)
LCLUST=$2       # Número del último cluster del tejido

# Ejecución

echo "Executing $0 ..."

TISSFN=$( echo $TISSUE | tr [[:upper:]] [[:lower:]] | tr " " "-")

echo "Filtering tissue $TISSFN..."

TISSLINES=$(cat rna_single_cell_read_count.tsv | grep "$TISSUE")

for ((i=0;i<=LCLUST;i++)); do
    echo "Filtering cluster $i..."
    head -1 rna_single_cell_read_count.tsv > ./HPA/datasets/$TISSFN/$TISSFN-$i-dataset.tsv
    echo $TISSLINES | grep -P "$TISSUE\t[0-9]+\t$i\t" >> ./HPA/datasets/$TISSFN/$TISSFN-$i-dataset.tsv
done

echo "Finished! :)"