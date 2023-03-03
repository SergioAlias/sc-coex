# Sergio Alías, 20220602
# Last modified 20221110

# Script para convertir IDs de ENSEMBL a IDs de Entrez

import json
# from datetime import datetime


dataset = "TabulaSapiens"

## Cargamos los IDs de Ensembl de HPA y las relaciones Ensembl-Entrez de BioMart

with open(f'{dataset}-genes-ENSEMBL', 'r') as f: ensembl = [line.strip('\n') for line in f.readlines()]

with open('ENSEMBL_gene_name_Entrez.txt', 'r') as f: ensembl_gname_entrez = [line.strip('\n').split('\t') for line in f.readlines()]


## Creramos diccionarios Entrez-Ensembl y Ensembl-Entrez

temp_dict = {line[0]: line[2] for line in ensembl_gname_entrez}

ensembl_entrez = {e: temp_dict[e] for e in ensembl if e in temp_dict.keys() and temp_dict[e] != ''}

print(f'Ensembl IDs in {dataset}: {len(ensembl)}')
print(f'Ensembl IDs with Entrez ID: {len(ensembl_entrez)}')
print(f'{len(ensembl) - len(ensembl_entrez)} IDs lost')

entrez_ensembl = {v: k for k, v in ensembl_entrez.items()}

## Save the dicts as JSON

# date = datetime.now().strftime('%Y%m%d')

with open(f'ensembl_to_entrez_{dataset}.json', 'w') as f: json.dump(ensembl_entrez, f,  indent = 4)

print(f'Ensembl to Entrez relationships saved --> ensembl_to_entrez_{dataset}.json')

with open(f'entrez_to_ensembl_{dataset}.json', 'w') as f: json.dump(entrez_ensembl, f,  indent = 4)

print(f'Entrez to Ensembl relationships saved --> entrez_to_ensembl_{dataset}.json')

