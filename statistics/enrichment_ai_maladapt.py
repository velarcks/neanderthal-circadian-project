#!/usr/bin/env python
# enrichment_ai_maladapt.py


BASE_URL = 'https://raw.githubusercontent.com/xzhang-popgen/maladapt/main/empirical_scores/'
FILE_NAMES = '../data/maladapt_files.txt'
snps = '../data/circadian_variants_introgressed.bed'
OUTPUT = 'output/'


import pandas as pd
import os
from ai_functions import circadian_ai_enrichment


def column_name_mapping(df):
    mapping = {df.columns[0]: 'chr', 
               df.columns[1]: 'start', 
               df.columns[2]: 'end'}
    df = df.rename(columns=mapping)
    return df


#
files = [n.strip() for n in open(FILE_NAMES).readlines()]
files = [n for n in files if n.startswith('nea')]

# IMPORT INDIVIDUAL PREDICTIONS AND CONCATENATE
pred = pd.DataFrame()
for file in files:
    url = os.path.join(BASE_URL,file)
    f = pd.read_csv(url, sep=',')
    pred = pd.concat([pred,f])

# CHANGE COLUMN DATA TYPES
cols = pred.columns.values.tolist()
pred[cols[:-1]] = pred[cols[:-1]].astype('int')

# ADD CHR STRING
if not pred[cols[0]].map(str).values[0].startswith('chr'):
    pred[cols[0]] = pred[cols[0]].map(str).replace('^','chr',regex=True)

# 
pred = pred.groupby(cols[:-1], as_index=False)[cols[-1]].max()


result = circadian_ai_enrichment(pred,snps,0.9,1000)
print(result)


with open(OUTPUT, 'w') as f:
    f.write(result)




#import os
#os.system("python ./enrichment_ai.py ../data/Nea_to_CEU_af-0.25_w-0.02.tsv ../data/circadian_variants_introgressed.bed 0.5 5 -o ../data/enrichment_ai_genomatnn.txt")
