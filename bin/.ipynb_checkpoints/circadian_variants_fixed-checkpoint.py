#!/usr/bin/env python
# circadian_variants_fixed.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Extract circadian variants from human and archaic specific variants
# published in Kuhlwilm et al. 2019
# 
"""""""""""""""""""""""""""


# INPUT DATA DIRECTORY
HHMC_FILE = '../data/kuhlwilm19_human_fixed_tony.tsv.gz'
AHMC_FILE = '../data/kuhlwilm19_archaic_fixed.tsv.gz'
GENE_LOCI = '../data/circadian_genes.bed'
FLANK_LOCI = '../data/circadian_genes_flanking.bed'


import pandas as pd
import numpy as np
import pybedtools
#import warnings
#warnings.filterwarnings('ignore')


def load_kuhlwilm(filename):
    file = pd.read_csv(filename, sep='\t', compression='gzip')
    file = parse_dfs(file)
    return file

def parse_dfs(df):
    df[['Chr', 'End']] = df['POS'].str.split(':', expand=True)
    df['Start'] = df['End'].map(int) - 1
    df['End'] = df['End'].map(int)
    df['Chr'] = df['Chr'].replace('^', 'chr', regex=True)
    df = df[['Chr','Start','End', 'human_DAF', 'Gene_name', 'REF', 'ALT', 'CAnc', 
        'Altai_allele', 'Vindija_allele', 'Denisova_allele', 'consequence', 'dbSNP']]
    return df

def load_circadian_gene_loci(file):
    df = pd.read_csv(file, sep='\t', low_memory=False).iloc[:,:-1]
    #df.rename(columns={'Start_hg19':'Start','End_hg19':'End'},inplace=True)
    return df

def circadian_and_fixed(fixed):
    gene = pd.merge(circadian_snps,fixed, on=['Chr','Start','End'])
    flank = pd.merge(circadian_1m,fixed, on=['Chr','Start','End'])
    
    return pd.concat([gene,flank])

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

def edit_circadian_kuhlwilm_dfs(file,region):
    df = file.iloc[:,[0,1,2,3,4,-2,-1]]
    df.columns = df.columns.str.replace(r'(_wa|_wb)', '')
    df['Region'] = region
    return df


# Load Kuhlwilm variants
hhmc = load_kuhlwilm(HHMC_FILE).iloc[:,[0,1,2,5,6,11]]
ahmc = load_kuhlwilm(AHMC_FILE).iloc[:,[0,1,2,5,6,11]]

# Load circadian loci
circadian_genes = load_circadian_gene_loci(GENE_LOCI)
circadian_1M = load_circadian_gene_loci(FLANK_LOCI)

# INTERSECT KUHLWILM VARIANTS AND CIRCADIAN LOCI
hhmc_gene = intersect_a_and_b(hhmc,circadian_genes)
hhmc_flank = intersect_a_and_b(hhmc,circadian_1M)
ahmc_gene = intersect_a_and_b(ahmc,circadian_genes)
ahmc_flank = intersect_a_and_b(ahmc,circadian_1M)

# Select columns
hhmc_gene = edit_circadian_kuhlwilm_dfs(hhmc_gene, 'Gene')
hhmc_flank = edit_circadian_kuhlwilm_dfs(hhmc_flank, 'Flanking')
ahmc_gene = edit_circadian_kuhlwilm_dfs(ahmc_gene, 'Gene')
ahmc_flank = edit_circadian_kuhlwilm_dfs(ahmc_flank, 'Flanking')

# Concatenate hhmcs and ahmcs
hhmc_circ = pd.concat([hhmc_gene,hhmc_flank])
ahmc_circ = pd.concat([ahmc_gene,ahmc_flank])


# SAVE
hhmc_circ.to_csv('../data/circadian_variants_hhmc.bed', index=False, sep='\t')
ahmc_circ.to_csv('../data/circadian_variants_ahmc.bed', index=False, sep='\t')
