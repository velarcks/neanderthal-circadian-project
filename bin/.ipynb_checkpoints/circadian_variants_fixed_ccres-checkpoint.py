#!/usr/bin/env python
# circadian_variants_fixed_ccres.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Extract human specific and Neanderthal specific (Kuhlwilm et al. 2019)
# circadian variants in regulatory regions.
# 
"""""""""""""""""""""""""""


# INPUT FILES
HHMC_FILE = '../data/circadian_variants_hhmc.bed'
AHMC_FILE = '../data/circadian_variants_ahmc.bed'
CCRE_FILE = '../data/GRCh38-ccREs.liftOver.to.Hg19.bed.gz'


# --------------------------------------------------
import pybedtools
import pandas as pd
#import warnings
#warnings.filterwarnings('ignore')


# --------------------------------------------------
def load_fixed(file):
    df = pd.read_csv(file, sep='\t', low_memory=False)
    return df

def get_fixed_ccres(fixed,ccre):
    df = intersect_a_and_b(fixed,ccre)
    df = df.iloc[:,:8]
    df.columns = df.columns.str.replace('_wa|_wb', '', regex=True)
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


# LOAD FILES
ahmc = load_fixed(AHMC_FILE)
hhmc = load_fixed(HHMC_FILE)
ccres = pd.read_csv(CCRE_FILE, sep='\t', compression='gzip')#.iloc[:,[0,1,2,-1]]

# FIND FIXED CIRCADIAN CCRES
ahmc_ccres = get_fixed_ccres(ahmc,ccres)
hhmc_ccres = get_fixed_ccres(hhmc,ccres)


# SAVE
hhmc_ccres.to_csv('../data/circadian_variants_hhmc_ccres.bed', index=False, sep='\t')
ahmc_ccres.to_csv('../data/circadian_variants_ahmc_ccres.bed', index=False, sep='\t')
