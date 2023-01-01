#!/usr/bin/env python
# circadian_variants_gene.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Get variants inside circadian genes
# 
"""""""""""""""""""""""""""


# INPUT DATA DIRECTORY
KGP_DIR = '../data/circadian_variants_all.bed.gz'
GENES_DIR = '../data/circadian_genes.bed'
# OUTPUT DATA DIRECTORY
OUTPUT_FILE = '../data/circadian_variants_gene.bed'


import pandas as pd
import pybedtools


def get_variants(genomesprj,circadian):
    genomesprj = genomesprj.add_suffix('_wa')
    circadian = circadian.add_suffix('_wb')
    kgp_cols = genomesprj.columns.values.tolist()
    genes_cols = circadian.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(genomesprj)
    b = pybedtools.BedTool.from_dataframe(circadian)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=kgp_cols+genes_cols)
    return a_and_b_df


# LOAD DATA
genomesprj = pd.read_csv(KGP_DIR, sep='\t', compression='gzip')
circadian_genes = pd.read_csv(GENES_DIR, sep='\t').iloc[:,:-1]

# MERGE CIRCADIAN LOCI AND 1KGP VARIANTS
circadian_variants = get_variants(genomesprj,circadian_genes)

# FILTER COLUMNS
circadian_variants = circadian_variants.iloc[:,[0,1,2,3,4,-2,-1]].drop_duplicates()

# REMOVE COLUMN NAME SUFFIXES
circadian_variants.columns = circadian_variants.columns.str.replace('(_wa|_wb)', '', regex=True)


# SAVE
circadian_variants.to_csv(OUTPUT_FILE, index=False, sep='\t')
