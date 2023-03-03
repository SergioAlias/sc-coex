#!/bin/bash
# Sergio Alías, 20220317
# Last modified 20220317

# Script para lanzar children-terms.py en bucle

# Argumentos (PONER LA BARRA AL FINAL DEL PATH)

OBO=$1          # Path del fichero obo
TISSUE_LIST=$2  # Fichero con los tejidos a buscar
OUTDIR=$3       # Directorio de salida

# Ejecución

echo "Executing $0 ..."

TISSUES=$( cat $TISSUE_LIST | tr [[:upper:]] [[:lower:]] )

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

for TIS in $TISSUES; do
    echo "Looking for $TIS children terms..."
    TISNOSPACES=$( echo $TIS | tr ' ' '-')
    python3 ./children-terms.py $OBO $TIS > $OUTDIR$TISNOSPACES-children-terms.txt
    echo "Children terms for $TIS saved --> $OUTDIR$TISNOSPACES-children-terms.txt"
done

IFS=$SAVEIFS

echo "Finished!"

################################################