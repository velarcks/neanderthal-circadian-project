#!/usr/bin/env python
# circadian_genes_candidate.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Imports genes that have evidence of circadian function from 5 different sources.
# Assigns a circadian function confident level to each of the genes (high, medium, or low).
# https://www.wikipathways.org/index.php/Pathway:WP3594
# https://pathcards.genecards.org/card/circadian_rythm_related_genes
# https://www.ncbi.nlm.nih.gov/biosystems/1458220
#
"""""""""""""""""""""""""""


# INPUT DATA
BIOSYS_DIR = '../data/genes_biosystems.txt'
CGDB_DIR = '../data/genes_cgdb_experimental.txt'
GO_DIR = '../data/genes_go.txt'
GWAS_DIR = '../data/genes_gwas.txt'
MCMAHON_DIR = '../data/genes_mcmahon.txt'
CONVERSION_DIR = '../data/gene_accession_conv.tsv'
# OUTPUT DATA
OUTPUT_DIR = '../data/circadian_genes_candidate.tsv'


import pandas as pd
import numpy as np
from functools import reduce

def load_file(file):
    f = [i.strip().split('\t')[0] for i in open(file).readlines()]
    return f

def list_to_column(input_list,source):
    df = pd.DataFrame(input_list, columns=[('GeneName')])
    df[source] = True
    return df


# IMPORT GENES FROM EVIDENCE SOURCES
biosystems = load_file(BIOSYS_DIR)
mcmahon = load_file(MCMAHON_DIR)
gwas = load_file(GWAS_DIR)
go = load_file(GO_DIR)
cgdb = load_file(CGDB_DIR)


# CREATE A DF FOR EACH SOURCE OF EVIDENCE
biosys_df = list_to_column(biosystems,'Biosystems')
mcmahon_df = list_to_column(mcmahon,'McMahon')
gwas_df = list_to_column(gwas,'GWAS')
go_df = list_to_column(go,'GO')
cgdb_df = list_to_column(cgdb,'CGDB')


# MERGE DATAFRAMES
dfs = [biosys_df, cgdb_df, gwas_df, go_df, mcmahon_df]
genes = reduce(lambda left,right: pd.merge(left,right,on=['GeneName'], 
                        how='outer'), dfs).drop_duplicates()


# ASSIGN CANDIDATE GENE A CONFIDENCE LEVEL BASED ON NUMBER OF EVIDENCE SOURCES
genes.loc[(genes[['GO', 'CGDB', 'GWAS', 'Biosystems']].isnull().sum(axis=1) == 1),'Confidence']='High'
genes.loc[(genes[['GO', 'CGDB', 'GWAS', 'Biosystems']].isnull().sum(axis=1) == 2),'Confidence']='Medium'
genes.loc[(genes[['GO', 'CGDB', 'GWAS', 'Biosystems']].isnull().sum(axis=1) >= 3),'Confidence']='Low'
genes.loc[(genes.McMahon.notnull()),'Confidence']='High'
genes.sort_values(by='GeneName', inplace=True)


# SAVE
genes = genes.fillna('NaN')
genes.to_csv(OUTPUT_DIR, sep='\t', index=False)
