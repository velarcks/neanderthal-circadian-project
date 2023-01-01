#!/usr/bin/env python
# enrichment_fixed_circadian_ccres.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Fisher's exact test on Kuhlwilm fixed variants in humans and archaics
# 
"""""""""""""""""""""""""""

# INPUT DATA
U = 4437803
CIRCADIAN_KUHLWILM = '../data/kuhlwilm19_circadian_and_flanking.bed.gz'
HHMC = '../data/kuhlwilm19_human_fixed_tony.tsv.gz'
AHMC = '../data/kuhlwilm19_archaic_fixed.tsv.gz'
CCRE_DIR = '../data/GRCh38-ccREs.liftOver.to.Hg19.bed.gz'
CIRCADIAN_1M = '../data/circadian_genes_flanking.bed'
# OUTPUT DATA
OUTPUT_FILE = '../data/enrichment_fixed_circadian_ccres.txt'


# --------------------------------------------------
from scipy import stats
import pybedtools
import numpy as np
import pandas as pd


# --------------------------------------------------
def fix_locus_column(df):
    new_df = pd.DataFrame(columns=['Chr','End'])
    new_df[['Chr','End']] = df['POS'].str.split(':', expand=True)
    new_df['Start'] = new_df['End'].map(int)-1
    new_df['End'] = new_df['End'].map(int)
    new_df['Chr'] = new_df['Chr'].replace('^','chr',regex=True)
    new_df = new_df[['Chr', 'Start', 'End']]
    return new_df

def get_x_axis():
    x = pd.read_csv(CIRCADIAN_1M, sep='\t', low_memory=False)
    x.rename(columns={'Start_hg19':'Start','End_hg19':'End'},inplace=True)
    ccres = load_ccres()
    circadian_ccres = intersect_a_and_b(ccres,x).iloc[:,:3].drop_duplicates()
    circadian_ccres.columns = circadian_ccres.columns.str.replace('_wa|_wb', '', regex=True)
    return circadian_ccres

def load_ccres():
    ccres = pd.read_csv(CCRE_DIR, sep='\t')
    kuhlwilm = pd.read_csv(CIRCADIAN_KUHLWILM, sep='\t')
    ccres_kuhlwilm = intersect_a_and_b(kuhlwilm,ccres).iloc[:,:3].drop_duplicates()
    ccres_kuhlwilm.columns = ccres_kuhlwilm.columns.str.replace('_wa|_wb', '', regex=True)
    return ccres_kuhlwilm

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

def get_y_axis(filename):
    file = pd.read_csv(filename, sep='\t', low_memory=False)
    y = fix_locus_column(file)
    y_len = len(file['POS'].drop_duplicates())
    return y
    

# --------------------------------------------------
def fishers_exact(x,row_tot,col_tot,tot,group):
    a = len(x)
    b = len(row_tot) - a
    c = len(col_tot) - a
    d = tot - (a + b + c)

    #table = [[a, b],[c, d]]
    ar = np.array([[a, b],[c, d]])
    
    oddsratio, pvalue = stats.fisher_exact(ar)
    
    fishers_exact_results = "\n{} specific enrichment\nOR: {}, \nP: {}.\n{}\n".format(
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
    # PREPARE X AXIS: All circadian variants in universe
    circadian_ccres = get_x_axis()

    # PREPARE Y AXIS: Fixed human and archaic variants in universe
    y_a = get_y_axis(AHMC)
    #y_h = get_y('raw/kuhlwilm19_human_fixed.tsv.gz')
    y_h_tony = get_y_axis(HHMC)

    # FIND Q2: Fixed circadian variants
    # Archaic
    q2_archaic = pd.merge(circadian_ccres,y_a, on=['Chr','Start','End'])
    # Human
    q2_human = pd.merge(circadian_ccres,y_h_tony, on=['Chr','Start','End'])

    # FISHER'S: Archaic specific.
    aar = fishers_exact(q2_archaic,circadian_ccres,y_a,U,'Archaic')
    print(aar)
    #print(contingency_table(aar,'Archaic'))

    # FISHER'S: Human specific
    har = fishers_exact(q2_human,circadian_ccres,y_h_tony,U,'Human')
    print(har)
    #print(contingency_table(har,'Human'))

    
    # SAVE RESULTS
    with open(OUTPUT_FILE, 'w') as output_file:
        output_file.write(aar + har)


# --------------------------------------------------
if __name__ == '__main__':
    main()


