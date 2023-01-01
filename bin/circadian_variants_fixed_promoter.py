#!/usr/bin/env python
# circadian_variants_fixed_promoter.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Extract human specific and Neanderthal specific (Kuhlwilm et al. 2019) 
# circadian variants in promoter regions.
# 
"""""""""""""""""""""""""""


# INPUT DIRECTORIES
PROMOTER_DIR = '../data/circadian_promoter_hg38.liftover.hg19.bed'
AHMC_DIR = '../data/kuhlwilm19_archaic_fixed.tsv.gz'
HHMC_DIR = '../data/kuhlwilm19_human_fixed_tony.tsv.gz'


import pybedtools
import pandas as pd


def fix_locus_column(df):
    df[['Chr','End']] = df['POS'].str.split(':', expand=True)
    df['Start'] = df['End'].map(int)-1
    df['End'] = df['End'].map(int)
    df['Chr'] = df['Chr'].replace('^','chr',regex=True)
    df = df[['Chr', 'Start', 'End', 'REF', 'ALT']]
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

def load_fixed(filename):
    # KUHLWILM
    df = pd.read_csv(filename, sep='\t', compression='gzip')
    df = fix_locus_column(df)
    return df

def find_circadian(df1,df2):
    df = intersect_a_and_b(df1,df2).iloc[:,[0,1,2,3,4,-2,-1]]
    df.columns = df.columns.str.replace('_wa|_wb', '')
    return df


# LOAD CIRCADIAN PROMOTER REGIONS
promoter = pd.read_csv(PROMOTER_DIR, sep='\t')
circadian = pd.read_csv(CIRCADIAN, sep='\t')

# LOAD SETS OF FIXED VARIANTS
ahmc = load_fixed(AHMC_DIR)
hhmc = load_fixed(HHMC_DIR)

# FIND FIXED CIRCADIAN VARIANTS IN PROMOTER REGION
ahmc_circadian = find_circadian(ahmc,promoter)
hhmc_circadian = find_circadian(hhmc,promoter)


# SAVE
hhmc_circadian.to_csv('../data/circadian_variants_hhmc_promoters.bed', index=False, sep='\t')
ahmc_circadian.to_csv('../data/circadian_variants_ahmc_promoters.bed', index=False, sep='\t')
