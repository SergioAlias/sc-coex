# -*- coding: utf-8 -*-
"""
@author: moran
"""
import csv
import numpy as np
import pandas as pd
from scipy import stats
#import os


tissues = ['liver-0-all',
'liver-1-all',
'liver-2-all',
'liver-3-all',
'liver-4-all',
'liver-5-all',
'liver-6-all',
'liver-7-all',
'liver-8-all',
'liver-9-all',
'liver-10-all',
'liver-11-all',
'liver-12-all',
'liver-13-all',
'liver-14-all',
'liver-15-all',
'liver-16-all']



# Create process and genes dictionary
# =====================================

processes_and_genes_file = csv.reader(open('process_and_associated_genes.csv', 'r'))  # processes and associated genes after filtering by GO domain ('biological process') and size 
proc_genes_dict = {}
for row in processes_and_genes_file:
    proc_genes_dict[row[0]] = row[1:]
    

def FC_tissue_v7(current_tissue):
    genes_FC_in_tissue_v7 = {}
    tis_file = csv.reader(open('tissues_FC_directory/' + current_tissue + '.csv', 'r'))  
    next(tis_file)
    for row in tis_file:
        gene = row[1]
        val = float(row[2])
        genes_FC_in_tissue_v7[gene] = val

    return genes_FC_in_tissue_v7


def TiPA_for_tissue(proc_genes_dict, current_tissue, genes_FC_tissue):
    
    ''' This function gets dictionary of processes and associated genes, 
    query tissue and FC vlaues for the tissue, and returns TiPA scores 
    alnog the size of the process'''
        
    scoring_and_length_per_tissue = {}

    for process in proc_genes_dict:
        list_of_tis_values = []
        for gene in proc_genes_dict[process]:
            if gene in genes_FC_tissue.keys():
                gene_value = float(genes_FC_tissue[gene])
                list_of_tis_values.append(gene_value)
        if len(list_of_tis_values) > 0:
            sorted_values = sorted(list_of_tis_values)
            tipa_score = stats.trim_mean(sorted_values, 0.1)  # calculation of TiPA score
            scoring_and_length_per_tissue[process] = (tipa_score)

    return scoring_and_length_per_tissue


def tipa_for_all_tissues():

    scoring_all_tissues = {}
    for tissue in tissues:
        genes_FC_tissue = FC_tissue_v7(tissue)
        scoring_all_tissues[tissue] = TiPA_for_tissue(proc_genes_dict, tissue, genes_FC_tissue)
    
    tipa_scores = {}
    for tissue in scoring_all_tissues:
        for process in scoring_all_tissues[tissue]:
            if process not in tipa_scores:
                tipa_scores[process] = {}
            tipa_scores[process][tissue] = scoring_all_tissues[tissue][process]
            
    scores_matrix = pd.DataFrame.from_dict(tipa_scores).transpose()
    scores_matrix.to_csv('tipa_scores_matrix_HPA_liver_all.csv')

    return tipa_scores  



tipa_for_all_tissues()



    
