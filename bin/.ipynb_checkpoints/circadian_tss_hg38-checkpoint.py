#!/usr/bin/env python
# circadian_tss_hg38.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Transcription start sites in circadian genes were downloaded from Biomart
# Promoter regions are defined as: 5kb upstream and 1kb downstream of the TSS.
# This script creates a bed file containing the promoter region for each circadian gene.
# The output is hg38. Use intersectBed to liftOver to h19
# liftOver circadian_tss_hg38.bed /dors/capra_lab/data/ucsc/liftOver/hg38ToHg19.over.chain.gz \ 
#  raw/circadian_tss_hg38.liftover.hg19.bed raw/circadian_tss_hg38.liftover.hg19_unlifted.bed
# 
"""""""""""""""""""""""""""


TSS_DIR = '../data/circadian_tss_hg38.tsv'
CIRCADIAN_DIR = '../data/circadian_genes.bed'
CIRCADIAN_L = '../data/circadian_genes.list'
OUTPUT = '../data/circadian_promoter_hg38.bed'


import pandas as pd
import numpy as np


# LOAD CIRCADIAN GENES
circadian_bed = pd.read_csv(CIRCADIAN_DIR, sep='\t', usecols = ['Chr', 'GeneID'])
circadian_l = pd.read_csv(CIRCADIAN_L, sep='\t')

# TSS DOWNLOADED FROM BIOMART: GeneID and hg38 TSS loci
tss_hg38 = pd.read_csv(TSS_DIR, sep='\t')
tss_hg38 = pd.merge(tss_hg38,circadian_l, on='GeneID',how='inner').iloc[:,[0,2,1]]

# Select the first TSS site for each gene
tss_hg38 = tss_hg38.groupby(['GeneID','GeneName']).min().reset_index()

# MERGE CIRCADIAN AND TSS. Adds chromosome name
df = pd.merge(circadian_bed,tss_hg38, on='GeneID',how='inner')
# ADD EMPTY Start/End COLUMNS
df.insert(1, 'Start', np.NaN)
df.insert(2, 'End', np.NaN)

# FILL Start/End COLUMNS WITH -5K UP AND 1K DOWNSTREAM
df['Start'] = df['TSS'] - 5000
df['End'] = df['TSS'] + 1000

# RENAME Chr TO READY IT TO LIFTOVER
df.rename(columns={'Chr':'#Chr'},inplace=True)


# SAVE
df.to_csv(OUTPUT, index=False, sep='\t')
