# Sergio Alías, 20220923
# Last modified 20220923

# Make a file with all HPO IDs and HPO term names of HPO
# Output used as input for plot_coex.R

import networkx
import obonet


# Main script

graph = obonet.read_obo("hp.obo") # Graph creation from obo file


## Mapping IDs and names

id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}


## Creating the file

with open('hpo-id-to-name', 'w') as f:
    for i in id_to_name:
        f.write(f'{i.replace(":", "")}\t{id_to_name[i]}\n')