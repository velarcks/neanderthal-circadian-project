#!/usr/bin/env python
# circadian_variants_promoter.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Find circadian variants in promoter regions using circadian promoter regions file
# 
"""""""""""""""""""""""""""


# DATA DIRECTORIRES
TSS_DIR = '../data/circadian_promoter_hg38.liftover.hg19.bed'
KGP_DIR = '../data/circadian_variants_all.bed.gz'
OUTPUT_FILE = '../data/circadian_variants_promoters.bed'


import pybedtools
import pandas as pd


def column_name_mapping(df):
    mapping = {df.columns[0]: 'Chr', 
               df.columns[1]: 'Start', 
               df.columns[2]: 'End',
               df.columns[3]: 'Ref',
               df.columns[4]: 'Alt'
              }
    df = df.rename(columns=mapping)
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


# LOAD PROMOTER LOCI
promoter = pd.read_csv(TSS_DIR, sep='\t')
genomesprj = pd.read_csv(KGP_DIR, sep='\t', compression='gzip')

# INTERSECT PROMOTER REGIONS AND 1KGP VARIANTS
snps = intersect_a_and_b(genomesprj,promoter)

# FILTER COLUMNS
snps = snps.iloc[:,[0,1,2,3,4,-2,-1]].drop_duplicates()

# REMOVE SUBSTRING FROM COLUMN NAMES
snps.columns = snps.columns.str.replace('_wa|_wb', '', regex=True)


# SAVE
snps.to_csv(OUTPUT_FILE, index=False, sep='\t')

