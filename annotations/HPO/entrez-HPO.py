# Sergio Alías, 20220602
# Last modified 20220603

# Script para obtener las relaciones entre términos HPO y genes

import json
from datetime import datetime


## Carga de datos

### Carga de las relaciones Ensembl ID- Entrez ID

gene_id_file = 'entrez_to_ensembl_20220603.json'

with open(gene_id_file) as f: entrez_ensembl = json.load(f)

### Carga del fichero de anotaciones de HPO- genes (HPO ID \t HPO name \t Entrez ID \t Gene symbol \t ...)

with open('phenotype_to_genes.txt', 'r') as f: phe_to_genes = [line.strip('\n').split('\t') for line in f.readlines()][1:]


## Relaciones HPO- genes

hpo_genes = {line[0]: {'HPOname': line[1], 'EntrezID': [], 'EnsemblID': [], 'GeneName': []} for line in phe_to_genes}

for line in phe_to_genes:
    if line[2] in entrez_ensembl.keys():
        hpo_genes[line[0]]['EntrezID'].append(line[2])
        hpo_genes[line[0]]['EnsemblID'].append(entrez_ensembl[line[2]])
        hpo_genes[line[0]]['GeneName'].append(line[3])

hp_to_remove = []
for hp in hpo_genes:
    if len(hpo_genes[hp]['EntrezID']) == 0: hp_to_remove.append(hp)

for hp in hp_to_remove: del hpo_genes[hp]


## Save the dict as JSON

date = datetime.now().strftime('%Y%m%d')

with open(f'hpo_to_genes_{date}.json', 'w') as f: json.dump(hpo_genes, f,  indent = 4)

print(f'HPO to genes relationships saved --> hpo_to_genes_{date}.json')