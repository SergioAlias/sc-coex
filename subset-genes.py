# Sergio Alías, 20220622
# Last modified 20221007

# Script para extraer un subconjunto de genes de una lista de genes
# El output servirá de input para plot_coex.R

###########################################################################################

## Imports

import argparse
import os
import json
# from datetime import datetime


###########################################################################################

## Argumentos de línea de comandos

parser = argparse.ArgumentParser()

parser.add_argument('hpo', type=str,
                    help='HPO term')

parser.add_argument('tissue', type=str,
                    help='Tissue name')

parser.add_argument('n_clusters', type=int,
                    help='Number of clusters of the tissue')

parser.add_argument('dataset', type=str,
                    help='Dataset used')

args = parser.parse_args()


###########################################################################################


## Carga de datos

### Relaciones HPO-genes

hpo_genes_filename = f'annotations/HPO/hpo_to_genes_{args.dataset}.json'

with open(hpo_genes_filename) as f: hpo_genes = json.load(f)

### Relaciones Entrez-Ensembl

ez_embl_filename = f'annotations/HPO/entrez_to_ensembl_{args.dataset}.json'

with open(ez_embl_filename) as f: entrez_ensembl = json.load(f)


###########################################################################################


## Main script

print('Ejecutando subset-genes.py ...\n')

# hpo = 'HP:0006561'
hpo = f'{args.hpo[:2]}:{args.hpo[2:]}'

tissue = args.tissue

n_clusters = args.n_clusters

dataset = args.dataset

dictkey = "GeneName" if dataset == "TabulaSapiens" else "EntrezID"

for c in range(n_clusters):
    cluster = c
    # Genes en la matriz de co-expresión (obtenida en R con COTAN)
    coex_file = f'{dataset}/results/{tissue}/genes-{tissue}-{cluster}.tsv'
    try:
        with open(coex_file, 'r') as f:
            coex_genes = [line.strip('\n').split('\t')[1].strip('"') for line in f.readlines()[1:]]
        print(f'Obteniendo genes asociados a {hpo} en la coex matrix de {tissue}-{cluster}...\n')
        if dataset == "TabulaSapiens":
            genes = [i for i in hpo_genes[hpo][dictkey] if i in coex_genes]
        else:
            genes = [entrez_ensembl[i] for i in hpo_genes[hpo][dictkey] if entrez_ensembl[i] in coex_genes]
        genes = [*dict.fromkeys(genes)] # Usamos un diccionario intermedio para eliminar duplicados

        print(f'Genes asociados a {hpo} presentes en la coex matrix de {tissue}-{cluster}: {len(genes)}\n')

    # date = datetime.now().strftime('%Y%m%d')

        with open(f'coex-analysis/{dataset}/{tissue}/{hpo.replace(":", "")}/{hpo.replace(":", "")}_subset_genes_{tissue}-{cluster}.txt', 'w') as f:
            for g in genes: f.write(f'{g}\n')

        print(f'Genes guardados en coex-analysis/{dataset}/{tissue}/{hpo.replace(":", "")}/{hpo.replace(":", "")}_subset_genes_{tissue}-{cluster}.txt\n')

    except FileNotFoundError:
        print(f'No se encontó {tissue}-{cluster}. Esto podrá (o no) tener consecuencias más adelante. Si las tuviera, revisa pasos previos a este script')

print('Finished! :D\n')
