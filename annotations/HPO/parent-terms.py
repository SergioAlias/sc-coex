# Sergio Alías, 20220519
# Last modified 20220519

# Extract parent terms of a given term of a ontology
# Tab-separated output: 1 HPO term, 2 term name, 3 HPO parent term, 4 parent term name, 5 relationship

import argparse
import networkx
import obonet


# Arguments
parser = argparse.ArgumentParser()

parser.add_argument('obo', type=str,
                    help='Obo file or link')

parser.add_argument('HPO', type=str,
                    help='HPO term') # 'HP:0010982'

args = parser.parse_args()


# Main script

graph = obonet.read_obo(args.obo) # Graph creation from obo file

## Mapping IDs and names

id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

## Find edges to children terms

HPO_term = args.HPO

def get_file(hpoTerm, exclude = "ignore", *args):
    print('term id\tterm name\tparent term id\tparent term name\tkey')
    id_list = [hpoTerm]
    terms_avaiable = True
    while terms_avaiable:
        for ub_id in id_list:
            subterms = 0
            for child, parent, key in graph.out_edges(ub_id, keys=True):
                if exclude == "ignore": condition = True
                if exclude == True: condition = key not in args
                if exclude == False: condition = key in args
                if condition:
                    print(f'{child}\t{id_to_name[child]}\t{parent}\t{id_to_name[parent]}\t{key}')
                    subterms += 1
                    id_list.append(parent)
            del id_list[0]
            if len(id_list) == 0 or id_list[0] == 'HP:0000001': terms_avaiable = False

get_file(HPO_term)

##########################################################
# Zona de pruebas (ejecutar bajo propia responsabilidad) #
##########################################################

# get_file('HP:0010982')