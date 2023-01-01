#!/usr/bin/env python
# circadian_variants_ccres.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Find cCREs that are flanking (1Mb) circadian genes.
# 
"""""""""""""""""""""""""""


# DATA DIRECTORY
CCRES_DIR = '../data/GRCh38-ccREs.liftOver.to.Hg19.bed.gz'
CIRCADIAN_DIR = '../data/circadian_variants_all.bed.gz'
OUTPUT_DIR = '../data/circadian_variants_ccres.bed'


import pybedtools
import pandas as pd


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


# LOAD DATA
ccres = pd.read_csv(CCRES_DIR, sep='\t', compression='gzip')
circadian_snps = pd.read_csv(CIRCADIAN_DIR, sep='\t', compression='gzip')#low_memory=False

# INTERSECT CIRCADIAN FLANKING SNPS AND CCRES
circadian_ccres = intersect_a_and_b(circadian_snps,ccres)

# FILTER COLUMNS
circadian_ccres = circadian_ccres.iloc[:,[0,1,2,3,4,5,6,-1]]

# REMOVE SUBSTRING FROM COLUMN NAMES
circadian_ccres.columns = circadian_ccres.columns.str.replace('_wa|_wb', '', regex=True)


# SAVE
circadian_ccres.to_csv(OUTPUT_DIR, index=False, sep='\t')
