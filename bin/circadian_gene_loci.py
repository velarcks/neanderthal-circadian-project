#!/usr/bin/env python
# circadian_gene_loci.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
#
# Description: Retrieve loci for circadian genes and circadian flanking regions.
# Return a bed format dataframe containing medium to high confidence circadian genes.
# 
"""""""""""""""""""""""""""


# INPUT DATA DIRECTORY
CIRCADIAN_DIR = '../data/circadian_genes_candidate.tsv'
GENCODE_DIR = '../data/gencode.v29lift37.circadian.gtf.gz'
#GENCODE_DIR = '../data/temp_2.txt'
# OUTPUT DATA DIRECTORY
GENES_OUT = '../data/circadian_genes.bed'
FLANKING_OUT = '../data/circadian_genes_flanking.bed'
GENE_LIST = '../data/circadian_genes.list'


import pandas as pd
import gzip
import io
import warnings
warnings.filterwarnings('ignore')


def load_circadian_genes():
    df = pd.read_csv(CIRCADIAN_DIR, sep='\t', usecols = ['GeneName', 'Confidence'])
    # Remove low confidence genes
    df = df[~df.Confidence.isin(['Low'])]
    return df

def raw_gencode(gene_or_exon):
    with io.TextIOWrapper(gzip.open(GENCODE_DIR, 'r')) as f:
        data = f.readlines()
    df = pd.DataFrame(data)
    df.rename(columns={0:'Chr'},inplace=True)
    df = df.Chr.str.split('\t',expand=True)
    
    df = df[df[2] == gene_or_exon]
    df = df[[0,3,4,8]]
    
    convert_dict = {3: int, 4: int}
    df = df.astype(convert_dict)
    return df

def column_name_mapping(df):
    mapping = {df.columns[0]: 'Chr', 
               df.columns[1]: 'Start', 
               df.columns[2]: 'End', }
    df = df.rename(columns=mapping)
    return df

def replace_gene_strings(df):
    df['GeneID'] = df.GeneID.str.replace('"|gene_id |.[0-9]+_[0-9]+', '', regex=True)
    df['GeneName'] = df.GeneName.str.replace('"|gene_name ', '', regex=True)
    return df

def split_gene_columns(df,col):
    # Parse gencode file and reformat last column into different columns
    df = column_name_mapping(df)
    temp_df = df[8].str.split('; ', expand = True)[[0,col]] # , n = 2
    df['GeneID'] = temp_df[0]
    df['GeneName'] = temp_df[col]
    df.drop(columns=[8],inplace = True)
    df = replace_gene_strings(df)
    return df
    
def add_flanking_region(df):
    # GET FLANKING UPSTREAM region
    df['1M_Up'] = df.Start - 1000000
    # GET FLANKING DOWNSTREAM region
    df['1M_Down'] = df.End + 1000000
    # Change value to 0 if Upstream limit is a negative value
    mask = df['1M_Up'] < 0
    df.loc[mask,'1M_Up'] = 0
    return df

def get_flanking_loci(df):
    # UP
    up_1M = df[['Chr','1M_Up','Start','GeneID', 'GeneName']]
    up_1M['Direction'] = 'Flanking up'
    
    # DOWN
    down_1M = df[['Chr','End','1M_Down','GeneID', 'GeneName']]
    down_1M['Direction'] = 'Flanking down'
    
    up_1M = column_name_mapping(up_1M)
    down_1M = column_name_mapping(down_1M)
    
    flanking = pd.concat([up_1M,down_1M])
    return flanking


# LOAD GENCODE DATA
gencode_0 = raw_gencode('gene')

# PARSE GENCODE RAW FILE. SPLIT LAST COLUMN. 'GeneName' column index = 2
gencode_0 = split_gene_columns(gencode_0,2)

# REMOVE GENEIDS ENDING IN _PAR_Y
gencode_0 = gencode_0[~gencode_0['GeneID'].str.endswith('_Y')]

# ADD FLANKING REGION
gencode_0 = add_flanking_region(gencode_0)

# LOAD CIRCADIAN GENE SET
circadian_0 = load_circadian_genes()

# MERGE GENCODE AND CIRCADIAN GENES
circadian_bed = pd.merge(gencode_0,circadian_0, on='GeneName', how='inner')

# EXTRACT FLANKING REGIONS
flanking = get_flanking_loci(circadian_bed)
circadian_bed.drop(['1M_Up','1M_Down'],axis=1,inplace=True)

# EXTRACT LIST OF GENES
gene_list = circadian_bed[['GeneID','GeneName']].sort_values(by='GeneName')


# SAVE
circadian_bed.to_csv(GENES_OUT, index=False, sep='\t')
flanking.to_csv(FLANKING_OUT, index=False, sep='\t')
gene_list.to_csv(GENE_LIST, index=False, sep='\t')
