# Sergio Alías, 20220602
# Last modified 20220603

# Script para obtener las relaciones entre términos HPO y genes (adaptado al dataset de Tabula Sapiens)

import json
# from datetime import datetime


## Carga de datos

### Carga de las relaciones Ensembl ID- Entrez ID

gene_symbol_file = 'TabulaSapiens-genes-GeneSymbol'
ENSEMBL_IDs_file = 'TabulaSapiens-genes-ENSEMBL'

with open(gene_symbol_file) as f: gs = [line.strip('\n') for line in f.readlines()]

with open(ENSEMBL_IDs_file) as f: ens = [line.strip('\n') for line in f.readlines()]

entrez_ensembl = {gs[i]: ens[i] for i in range(0, len(gs))}

### Carga del fichero de anotaciones de HPO- genes (HPO ID \t HPO name \t Entrez ID \t Gene symbol \t ...)

with open('phenotype_to_genes.txt', 'r') as f: phe_to_genes = [line.strip('\n').split('\t') for line in f.readlines()][1:]


## Relaciones HPO- genes

hpo_genes = {line[0]: {'HPOname': line[1], 'EntrezID': [], 'EnsemblID': [], 'GeneName': []} for line in phe_to_genes}

for line in phe_to_genes:
    if line[3] in entrez_ensembl.keys():
        hpo_genes[line[0]]['EntrezID'].append(line[2])
        hpo_genes[line[0]]['EnsemblID'].append(entrez_ensembl[line[3]])
        hpo_genes[line[0]]['GeneName'].append(line[3])

hp_to_remove = []
for hp in hpo_genes:
    if len(hpo_genes[hp]['EntrezID']) == 0: hp_to_remove.append(hp)

for hp in hp_to_remove: del hpo_genes[hp]


## Save the dict as JSON

# date = datetime.now().strftime('%Y%m%d')

with open(f'hpo_to_genes_TabulaSapiens.json', 'w') as f: json.dump(hpo_genes, f,  indent = 4)

print(f'HPO to genes relationships saved --> hpo_to_genes_TabulaSapiens.json')