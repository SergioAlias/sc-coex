# Sergio Alías, 20220602
# Last modified 20230503

# Script para generar dos ficheros de la forma:

# Fichero 1: 
## tissue \t HPO code \t HPO name \t num genes

# Fichero 2: 
# HPO code \t gene code \t gene name

import json
import os
import obonet
import networkx
# from datetime import datetime


dataset = "HPA"


## Carga de datos

### Carga de términos hijos UBERON

ub_subt_dir = f'UBERON/results/{dataset}/'

ub_subt_filenames = [ub_subt_dir + fn for fn in os.listdir(ub_subt_dir)]

ub_subt = {}

for p in ub_subt_filenames:
    with open(p, 'r') as f: f = [line.strip('\n').split('\t') for line in f.readlines()][1:]
    ub_subt[p.split('/')[-1][:-19]] = f

### Carga de las relaciones UBERON-HPO

ub_hpo_filename = 'HPO/ub_to_hpo_20220602.json'

with open(ub_hpo_filename) as f: uberon_hpo = json.load(f)

### Carga de las relaciones HPO-genes

hpo_genes_filename = f'HPO/hpo_to_genes_{dataset}.json'

with open(hpo_genes_filename) as f: hpo_genes = json.load(f)

### Carga de la ontología HPO

hpo_obo = './HPO/hp.obo'

graph = obonet.read_obo(hpo_obo) # Graph creation from obo file

#### Mapping IDs and names

# id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}
# name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}

#### Find edges to children terms (términos hijos HPO)

#### Aquí definimos una función en vez de importar un JSON con todas las relaciones porque
#### para términos muy altos en la jerarquía hay muchísimos términos hijos y creo que es 
#### preferible hacer la búsqueda solo de los HPO necesarios dentro de este mismo script

def child_hpo(hpo):
    """
    Genera una lista con el HPO y todos sus hijos

    Parámetros
    ----------
        hpo : str
            término HPO (ID)

    Returns
    -------
        children : list
            lista con el término HPO orginial y sus hijos, hijos de hijos, etc
    """
    children = [hpo]
    terms_avaiable = True
    id_list = [hpo]
    temp_list = []
    while terms_avaiable:
        for h in id_list:
            for parent, child, key in graph.in_edges(h, keys=True):
                children.append(parent)
                temp_list.append(parent)
        id_list = temp_list
        temp_list = []
        if len(id_list) == 0: terms_avaiable = False
    return children



## Cabeceras de los ficheros

tissue_hpo_ngenes = ['tissue\thpo_code\thpo_name\tnum.genes\n']
hpo_gcode_gname = ['hpo_code\tgene_code\tgene_name\n']



## Diccionario que contendrá UBERON: [hijos de UBERON]

## El input de las relaciones UBERON-términos hijos no es un JSON, sino un TSV por cada tissue
## Aquí los convertimos en un mismo diccionario

t_subt = {t: [] for t in ub_subt.keys()}

for t in ub_subt:
    t_subt[t].append(ub_subt[t][0][0])
    for line in ub_subt[t]:
        if line[2] not in t_subt[t]: t_subt[t].append(line[2])



## Creración de los ficheros

## Usamos un diccionario como paso intermedio para el primer fichero porque antes de construirlo
## me interesa ordenarlo por num_genes y eliminar duplicados dentro de un mismo tissue

file_1 = {t: [] for t in t_subt}

## Para el segundo creamos unas lista auxiliar, ya que luego nos interesa eliminar duplicados
## y así podemos llevar la cuenta de qué términos HPO7 genes han sido añadidos ya



## Rellenamos file_1 y construimos el fichero 2

for t in t_subt: # Iteramos por tejido
    file_2_hpo = [] # La lista de HPO se vacía al cambiar de tejido
    for ub in t_subt[t]: # Iteramos por término UBERON
        ub = ub.replace(':', '_')
        if ub in uberon_hpo.keys(): # Comprobamos que exista alguna relación UBERON-HPO
            for hpo in uberon_hpo[ub]: # Iteramos los términos HPO encontrados
                hpo = hpo.replace('_', ':')
                for hp in child_hpo(hpo): # Iteramos por el término HPO + sus hijos
                    if hp in hpo_genes.keys(): # Comprobamos que exista alguna relación HPO-gen
                        if hp not in file_1[t]: # Control de duplicados a nivel de HPO (fichero 1)
                            file_1[t].append(hp)
                        if hp not in file_2_hpo: # Control de duplicados a nivel de HPO (fichero 2)
                            file_2_hpo.append(hp)
                            file_2_genes = [] # La lista de genes usados se vacía cada vez que cambiamos de HPO
                            for g in range(len(hpo_genes[hp]["EntrezID"])): # Construimos el fichero 2
                                if hpo_genes[hp]["EntrezID"][g] not in file_2_genes: # Control de duplicados a nivel de gen (fichero 2)
                                    file_2_genes.append(hpo_genes[hp]["EntrezID"][g])
                                    hpo_gcode_gname.append(f'{hp}\t{hpo_genes[hp]["EntrezID"][g]}\t{hpo_genes[hp]["GeneName"][g]}\n')

## Construimos el fichero 1 a partir del diccionario file_1

for t in file_1.keys():
    hp_used = []
    line_ngenes = []
    for hp in file_1[t]:
        if hp not in hp_used:
            non_rep_genes = []
            for g in hpo_genes[hp]["EntrezID"]:
                if g not in non_rep_genes: # Control de duplicados a nivel de gen (fichero 1)
                    non_rep_genes.append(g)
            line_ngenes.append([f'{t}\t{hp}\t{hpo_genes[hp]["HPOname"]}\t{len(non_rep_genes)}\n', len(non_rep_genes)])
            hp_used.append(hp)
    line_ngenes.sort(key = lambda l: l[1], reverse = True) # Ordenamos de más a menos genes por cada tissue
    for i in line_ngenes:
        tissue_hpo_ngenes.append(i[0])



## Saving files

# date = datetime.now().strftime('%Y%m%d')


with open(f'tis_hpo_ngenes_{dataset}.tsv', 'w') as f:
    for line in tissue_hpo_ngenes: f.write(line)

print(f'File saved --> tis_hpo_ngenes_{dataset}.tsv')

with open(f'hpocode_genes_{dataset}.tsv', 'w') as f:
    for line in hpo_gcode_gname: f.write(line)

print(f'File saved --> hpocode_genes_{dataset}.tsv')
