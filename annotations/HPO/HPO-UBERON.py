# Sergio Alías, 20220523
# Last modified 20220602

# Script para obtener las relaciones HPO-UBERON

import json
from datetime import datetime


## Carga del fichero

hpowl_filename = 'hp.owl'

with open(hpowl_filename, 'r') as f: f = [line.strip('\n') for line in f.readlines()]

print(f'File {hpowl_filename} loaded\n...')


## Patrones para las líneas con términos HPO, UBERON asociado a HPO y fin de términos HPO

patt_hpo = '<owl:Class rdf:about="http://purl.obolibrary.org/obo/HP_'
patt_uberon = '<owl:someValuesFrom rdf:resource="http://purl.obolibrary.org/obo/UBERON_'
patt_end = '<owl:Class rdf:about="http://purl.obolibrary.org/obo/HsapDv_0000000">'


## División del fichero por términos HPO

line_index  = -1
hpo_begin, hpo_names = [], []


for line in f:
    line_index += 1
    if patt_hpo in line:
        hpo_begin.append(line_index)
        hpo_names.append(line[-12:-2])
    if patt_end in line:
        last_index = line_index


hpo_end = hpo_begin[1:] + [last_index]

hpo_blocks = {k: f[i:j] for i, j, k in zip(hpo_begin, hpo_end, hpo_names)}


## Búsqueda de términos UBERON en cada bloque HPO

hpo_uberon = {key: [] for key in hpo_blocks.keys()}

for key, text in hpo_blocks.items():
    for line in text:
        if patt_uberon in line:
            hpo_uberon[key].append(line[-17:-3])


## Invertimos el diccionario para tener como keys los términos UBERON

uberon_hpo = {} 

for hpo, uberon in hpo_uberon.items():
    for u in uberon:
        uberon_hpo.setdefault(u,[]).append(hpo)


## Save the dicts as JSON files

date = datetime.now().strftime('%Y%m%d')

with open(f'hpo_to_ub_{date}.json', 'w') as f: json.dump(hpo_uberon, f,  indent = 4)

print(f'HPO to UBERON relationships saved --> hpo_to_ub_{date}.json')

with open(f'ub_to_hpo_{date}.json', 'w') as f: json.dump(uberon_hpo, f,  indent = 4)

print(f'UBERON to HPO relationships saved --> ub_to_hpo_{date}.json')