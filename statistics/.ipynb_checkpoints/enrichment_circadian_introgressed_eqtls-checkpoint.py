#!/usr/bin/env python
# enrichment_circadian_introgressed_eqtls.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Hypergeometric distribution of introgressed variants in circadian genes 
#  with evidence of being eQTL in GTEx.
# 
"""""""""""""""""""""""""""


# INPUT DATA DIRECTORY
EQTLS = 4608446
CIRCADIAN_DIR = '../data/circadian_variants_gtex.bed.gz'
INTROGRESSED_DIR = '../data/gtex_v8_hg19_introgressed.bed.gz'
# OUTPUT DATA DIRECTORY
OUTPUT_FILE = '../data/enrichment_circadian_introgressed_eqtls.txt'


# --------------------------------------------------
import os, sys
import pandas as pd
import numpy as np
from scipy import stats
import pybedtools


def load_file(file):
    df = pd.read_csv(file, sep='\t', compression='gzip').iloc[:, 0:3].drop_duplicates()
    return df

def intersect_sets(df_a, df_b):
    df_a = df_a.add_suffix('_wa')
    df_b = df_b.add_suffix('_wb')
    a_cols = df_a.columns.values.tolist()
    b_cols = df_b.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(df_a)
    b = pybedtools.BedTool.from_dataframe(df_b)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=a_cols + b_cols)
    return a_and_b_df

def fishers_exact(U,right,down,x):
    # Fishers exact test
    a = len(x)
    b = len(right) - a
    c = len(down) - a
    d = U - (a + b + c)

    #table = [[a, b], [c, d]]
    ar = np.array([[a, b],[c, d]])

    oddsratio, pvalue = stats.fisher_exact(ar)

    fishers_results = '\nOR: {} \nP-VAL: {} \n{}'.format(oddsratio,pvalue,ar)

    
    # SAVE RESULTS
    with open(OUTPUT_FILE, 'w') as output_file:
        output_file.write(fishers_results)

    return fishers_results


# --------------------------------------------------
def main():
    # LOAD DATA
    circadian_df = load_file(CIRCADIAN_DIR)
    introgressed_df = load_file(INTROGRESSED_DIR)

    # FIND Q2
    q2 = intersect_sets(circadian_df,introgressed_df).iloc[:, 0:3].drop_duplicates()

    # FISHER'S EXACT TEST
    print(fishers_exact(EQTLS,circadian_df,introgressed_df,q2))


# --------------------------------------------------
if __name__ == '__main__':
    main()

