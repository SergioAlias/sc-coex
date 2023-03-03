# Sergio Alías, 20221026
# Last modified 20221026

# Script para aplicar sobre el fichero tis_hpo_ngenes, de la forma: 
# tissue \t HPO code \t HPO name \t num genes
# Se obtiene otro fichero con igual estructura en el que se filtran los HPO con 20 genes o más
# los cuales no tengan ningún término hijo con 20 genes o más


import obonet
import networkx


dataset = "TabulaSapiens"


## Carga de datos

### Ontología HPO

hpo_obo = './HPO/hp.obo'

graph = obonet.read_obo(hpo_obo) # Graph creation from obo file

### Fichero de términos HPO y número de genes

filename = f'./tis_hpo_ngenes_{dataset}.tsv'

with open(filename) as f: f = f.readlines()


## Funciones

### Términos hijos HPO

def childHPO(hpo):
    """
    Genera una lista con el HPO y todos los hijos de un término HPO (pero NO hijos de hijos)

    Parámetros
    ----------
        hpo : str
            término HPO (ID)

    Returns
    -------
        [unnamed list] : list
            lista con los hijos del término hpo
    """
    return [parent for parent, child, key in graph.in_edges(hpo, keys=True)]


## Main code

header, f = f[0], f[1:]

lines = [line.strip('\n').split('\t') for line in f]

gene_dict = {line[1]: int(line[3]) for line in lines}

outfile = [header]

for i in range(0, len(lines)-1):
    hpo = lines[i][1]
    if gene_dict[hpo] >= 20:
        valid = True
        for child in childHPO(hpo):
            if child in gene_dict.keys() and gene_dict[child] >= 20:
                valid = False
        if valid:
            outfile.append(f[i])


with open(f'filtered_tis_hpo_ngenes_{dataset}.tsv', 'w') as f:
    for line in outfile: f.write(line)
