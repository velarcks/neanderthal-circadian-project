#!/usr/bin/env python
# circadian_variants_introgressed.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Extract circadian variants in a set of introgressed variants
# from Browning et al., 2018
# 
"""""""""""""""""""""""""""


# INPUT DATA DIRECTORY
BROWNING = '../data/introgressed_variants_browning18_hg19.bed.gz'
CIRCADIAN_SNPS = '../data/circadian_variants_gene.bed'
CCRE_SNPS = '../data/circadian_variants_ccres.bed'
PROMOTER_SNPS = '../data/circadian_variants_promoters.bed'
# OUTPUT DATA DIRECTORY
OUTPUT_FILE = '../data/circadian_variants_introgressed.bed'


import pandas as pd
import pybedtools


def load_variant_sets(filename,reg):
    df = pd.read_csv(filename,sep='\t').iloc[:,[0,1,2,5,6]].drop_duplicates()
    df['Region'] = reg
    return df
    
def intersect_a_and_b(df_a,df_b):
    df_a = df_a.add_suffix('_wa')
    df_b = df_b.add_suffix('_wb')
    df_a_cols = df_a.columns.values.tolist()
    df_b_cols = df_b.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(df_a)
    b = pybedtools.BedTool.from_dataframe(df_b)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=df_a_cols+df_b_cols)
    return a_and_b_df

def collapse_columns(df,col):
    cols = df.columns.values.tolist()
    cols.remove(col)
    df[col] = df[col].astype(str)
    df = df.groupby(cols)[col].apply('|'.join).reset_index()
    return df


# LOAD FILES
browning18 = pd.read_csv(BROWNING, sep='\t', compression='gzip')
circ_genes = load_variant_sets(CIRCADIAN_SNPS, 'Gene')
circ_ccres = load_variant_sets(CCRE_SNPS, 'Regulatory')
circ_promoters = load_variant_sets(PROMOTER_SNPS, 'Promoter')

# MERGE ALL CIRCADIAN SNP DATAFRAMES
dfs = pd.concat([circ_genes,circ_ccres,circ_promoters])

# Collapse duplicated loci in multiple Regions
dfs = collapse_columns(dfs,'Region')

# INTERSECT INTROGRESSED AND CIRCADIAN VARIANTS
circ_introgressed = intersect_a_and_b(browning18,dfs)

# FILTER COLUMNS
circ_introgressed = circ_introgressed.iloc[:,[0,1,2,3,4,-3,-2,-1]]
circ_introgressed.columns = circ_introgressed.columns.str.replace('_wa|_wb', '', regex=True)


# SAVE
circ_introgressed.to_csv(OUTPUT_FILE, index=False, sep='\t')
