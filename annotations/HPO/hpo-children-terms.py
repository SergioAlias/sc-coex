# Sergio Alías, 20220607
# Last modified 20220607

# Extract children terms of a given term of a ontology (for instance, UBERON)
# JSON output: keys --> HPO id values --> HPO name, children terms IDs, children terms names

import networkx
import obonet
import json
from datetime import datetime

## Carga de la ontología

obo = './hp.obo'

graph = obonet.read_obo(obo) # Graph creation from obo file


## Mapping IDs and names

id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}


## Find edges to children terms

hpo = {id_: {'name': id_to_name[id_], 'ChildrenIDs': [], 'ChildrenNames': []} for id_ in graph.nodes()}

for hp in ['HP:0009875']:
    terms_avaiable = True
    id_list = [hp]
    temp_list = []
    while terms_avaiable:
        print('Start while loop')
        print(f'ID list: {id_list}')
        for h in id_list:
            print(f'Using {h}')
            for parent, child, key in graph.in_edges(h, keys=True):
                hpo[h]['ChildrenIDs'].append(parent)
                hpo[h]['ChildrenNames'].append(id_to_name[parent])
                temp_list.append(parent)
                print(f'------>  found {parent}')
                print(temp_list)
        id_list = temp_list
        temp_list = []
        print('New id list:')
        print(id_list)
        if len(id_list) == 0: terms_avaiable = False


## Save the dict as JSON files

date = datetime.now().strftime('%Y%m%d')

with open(f'tmp_hpo_children_terms_{date}.json', 'w') as f: json.dump(hpo, f,  indent = 4)

print(f'HPO children terms saved --> hpo_children_terms_{date}.json')