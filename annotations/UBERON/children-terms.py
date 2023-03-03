# Sergio Alías, 20220316
# Last modified 20220711

# Extract children terms of a given term of a ontology (for instance, UBERON)
# Tab-separated output: 1 UBERON term, 2 term name, 3 UBERON subterm, 4 subterm name, 5 relationship

import argparse
import networkx
import obonet


# Arguments

parser = argparse.ArgumentParser()

parser.add_argument('obo', type=str,
                    help='Obo file or link')

parser.add_argument('tissue', type=str,
                    help='tissue name')

args = parser.parse_args()


# Main script

graph = obonet.read_obo(args.obo) # Graph creation from obo file


## Mapping IDs and names

id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}


## Find edges to children terms

tissue = name_to_id[args.tissue]

def get_file(tissue, keylist = None, exclude = "ignore"):
    print('term id\tterm name\tchild term id\tchild term name\tkey')
    id_list = [tissue]
    used_ids = []
    temp_list = []
    repeated_list = []
    terms_avaiable = True
    while terms_avaiable:
        for ub_id in id_list:
            if ub_id not in used_ids:
                for parent, child, key in graph.in_edges(ub_id, keys=True):
                    if exclude == "ignore": condition = True
                    if exclude == True: condition = key not in keylist
                    if exclude == False: condition = key in keylist
                    if condition:
                        if parent not in repeated_list:
                            print(f'{child}\t{id_to_name[child]}\t{parent}\t{id_to_name[parent]}\t{key}')
                            temp_list.append(parent)
                            repeated_list.append(parent)
                used_ids.append(id_list[0])
        id_list = temp_list
        temp_list = []
        if len(id_list) == 0: terms_avaiable = False


## Llamada a la función

keys = ('part_of', 'is_a')

get_file(tissue, keys, exclude = False)


##########################################################
# Zona de pruebas (ejecutar bajo propia responsabilidad) #
##########################################################

# def mini_get_file(*args):
#     print('term id\tterm name\tsubterm\tsubterm name\tkey')
#     for arg in args:
#         for parent, child, key in graph.in_edges(arg, keys=True):
#             print(f'{child}\t{id_to_name[child]}\t{parent}\t{id_to_name[parent]}\t{key}')
