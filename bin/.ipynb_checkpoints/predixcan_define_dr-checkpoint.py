#!/usr/bin/env python
# predixcan_define_dr.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Find list of divergently regulated genes by extracting genes 
# that have a p-value = 0
# 
"""""""""""""""""""""""""""


BASE_URL = 'https://raw.githubusercontent.com/colbrall/neanderthal_predixcan_manuscript/master/data/'
ALTAI_FILE = 'altai_update_pvalues_2sided.txt'
VINDIJA_FILE = 'vindija_pvalues_2sided.txt'
DENISOVA_FILE = 'denisovan_update_pvalues_2sided.txt'

CIRCADIAN_FILE = '../data/circadian_genes.list'
TISSUE_MAP = '../data/tissue_name_map.txt'


# --------------------------------------------------
import os
import pandas as pd
from functools import reduce


# --------------------------------------------------
def map_tissue_names(df):
    # Map tissue names to GTEx tissue name format
    tissue_name_map = pd.read_csv(TISSUE_MAP, sep='\t')
    # Fix tissue names
    df = pd.merge(df,tissue_name_map,left_on='tissue',
                  right_on='PublishedArchaics_Tissue')[['gene_id','GTEx_Tissue']]
    return df


def get_dr(base_url,path):
    # Load file
    url = os.path.join(base_url,path)
    df = pd.read_csv(url, sep='\t')
    
    # Convert columns to rows using melt
    df = df.melt(id_vars=['gene_id'], 
                        var_name='tissue', value_name='p-value')
    
    # Remove extra numbers at the end of each gene id
    df['gene_id'] = df['gene_id'].replace(to_replace ='\.\d+', value ='',regex=True)
    
    # Filter by p-value
    df = df[df['p-value']==0.0]
    
    return map_tissue_names(df)


def get_circadian(archaic,circadian,ind):
    #df = pd.read_csv(archaic, sep='\t')
    df = pd.merge(circadian,archaic,left_on='GeneID',right_on='gene_id')[['GeneID','GeneName','GTEx_Tissue']]
    df[ind] = 1
    return df
    

# --------------------------------------------------
# FIND DR GENES
altai_signif = get_dr(BASE_URL,ALTAI_FILE)
vindija_signif = get_dr(BASE_URL,VINDIJA_FILE)
denisova_signif = get_dr(BASE_URL,DENISOVA_FILE)


# LOAD LIST OF CIRCADIAN GENES
circadian_genes = pd.read_csv(CIRCADIAN_FILE, sep='\t')

# FIND CIRCADIAN DR
altai_cdr = get_circadian(altai_signif,circadian_genes,'Altai')
vindija_cdr = get_circadian(denisova_signif,circadian_genes,'Vindija')
denisova_cdr = get_circadian(vindija_signif,circadian_genes,'Denisova')

# MERGE ALL ARCHAIC DR DFS
DFS = [altai_cdr, vindija_cdr, denisova_cdr]
COLS = ['GeneID','GeneName','GTEx_Tissue']
archaics_cdr = reduce(lambda left,right: pd.merge(left,right,on=COLS,how='outer'), DFS)

# FILL NaN VALUES WITH 0
archaics_cdr.iloc[:,-3:] = archaics_cdr.iloc[:,-3:].fillna(0).astype(int)
archaics_cdr.sort_values(by='GeneName', inplace=True)

# FIND DR GENES IN COMMON BETWEEN ALL ARCHAICS
circ_dr = archaics_cdr[archaics_cdr.sum(axis=1)==3].iloc[:,:-3]
circ_dr.sort_values(by='GeneName', inplace=True)


# SAVE
altai_cdr.to_csv('../data/predixcan_dr_circadian_altai.tsv', index=False, sep='\t')
vindija_cdr.to_csv('../data/predixcan_dr_circadian_vindija.tsv', index=False, sep='\t')
denisova_cdr.to_csv('../data/predixcan_dr_circadian_denisova.tsv', index=False, sep='\t')
archaics_cdr.to_csv('../data/predixcan_dr_circadian_archaics.tsv', index=False, sep='\t')
circ_dr.to_csv('../data/predixcan_dr_circadian.tsv', index=False, sep='\t')
