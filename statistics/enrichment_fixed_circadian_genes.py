#!/usr/bin/env python
# enrichment_fixed_circadian_genes.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Testing for enrichment of archaic and human specific variants (Kuhlwilm et al. 2019)
# inside circadian genes using Fisher's exact test.
# 
"""""""""""""""""""""""""""


# INPUT DATA DIRECTORY
CIRCADIAN_KUHLWILM = '../data/kuhlwilm19_circadian_and_flanking.bed.gz'
U = 4437803    # Total number of SNPs in Kuhlwilm19's set
HHMC = '../data/kuhlwilm19_human_fixed_tony.tsv.gz'
AHMC = '../data/kuhlwilm19_archaic_fixed.tsv.gz'
GENES_DIR = '../data/circadian_genes.bed'
# OUTPUT DATA DIRECTORY
OUTPUT_FILE = '../data/enrichment_fixed_circadian_genes.txt'


# --------------------------------------------------
from scipy import stats
import pybedtools
import numpy as np
import pandas as pd


def get_circadian():
    c = pd.read_csv(GENES_DIR, sep='\t')
    k = pd.read_csv(CIRCADIAN_KUHLWILM, sep='\t', compression='gzip')
    circadian_variants = intersect_a_and_b(k,c).iloc[:,:3].drop_duplicates()
    circadian_variants.columns = circadian_variants.columns.str.replace('_wa|_wb', '')
    return circadian_variants

def get_fixed(filename):
    file = pd.read_csv(filename, sep='\t', compression='gzip', low_memory=False)
    fixed = fix_locus_column(file)
    return fixed

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

def fix_locus_column(df):
    new_df = pd.DataFrame(columns=['Chr','End'])
    new_df[['Chr','End']] = df['POS'].str.split(':', expand=True)
    new_df['Start'] = new_df['End'].map(int)-1
    new_df['End'] = new_df['End'].map(int)
    new_df['Chr'] = new_df['Chr'].replace('^','chr',regex=True)
    new_df = new_df[['Chr', 'Start', 'End']]
    return new_df

    
# --------------------------------------------------
def fishers_exact(x,row_tot,col_tot,tot,group):
    a = len(x)
    b = len(row_tot) - a
    c = len(col_tot) - a
    d = tot - (a + b + c)

    #table = [[a, b],[c, d]]
    ar = np.array([[a, b],[c, d]])
    
    oddsratio, pvalue = stats.fisher_exact(ar)
    
    fishers_exact_results = "\n{} specific enrichment \nOR: {} \nP: {} \n{}\n".format(
                                                        group,oddsratio,pvalue,ar)
    return fishers_exact_results

def contingency_table(ar,group):
    df = pd.DataFrame(ar, columns=[f"{group} specific variants", f"Non-{group} specific variants"])
    df.index = ["Circadian variants", "Non-Circadian variants"] 
    df.loc['Column_Total']= df.sum(numeric_only=True, axis=0)
    df.loc[:,'Row_Total'] = df.sum(numeric_only=True, axis=1)
    return df


# --------------------------------------------------
def main():
    # FIND ALL CIRCADIAN VARIANTS IN KUHLWILM19'S SET
    circadian_variants = get_circadian()

    # FIXED VARIANTS
    archaic = get_fixed(AHMC)
    human = get_fixed(HHMC)

    # FIND Q2: Fixed circadian
    a_circadian = pd.merge(archaic,circadian_variants, on=['Chr','Start','End'])
    #
    h_circadian = pd.merge(human,circadian_variants, on=['Chr','Start','End'])

    # FISHER'S EXACT TEST
    # Archaic specific
    aar = fishers_exact(a_circadian,circadian_variants,archaic,U,'Archaic')
    print(aar)
    # contingency_table(aar,'Archaic')

    # Human specific
    har = fishers_exact(h_circadian,circadian_variants,human,U,'Human')
    print(har)
    # contingency_table(har,'Human')
    

    # SAVE RESULTS
    with open(OUTPUT_FILE, 'w') as output_file:
        output_file.write(aar + har)


# --------------------------------------------------
if __name__ == '__main__':
    main()

