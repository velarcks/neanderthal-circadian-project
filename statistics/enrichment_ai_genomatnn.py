#!/usr/bin/env python
#enrichment_ai_genomatnn.py


base = 'https://raw.githubusercontent.com/grahamgower/genomatnn/main/predictions/'
file = 'Nea_to_CEU_af-0.25_w-0.02.tsv'
snps = '../data/circadian_variants_introgressed.bed'
OUTPUT = 'output/'


from ai_functions import circadian_ai_enrichment
import pandas as pd
import os


def column_name_mapping(df):
    mapping = {df.columns[0]: 'chr', 
               df.columns[1]: 'start', 
               df.columns[2]: 'end'}
    df = df.rename(columns=mapping)
    return df


url = os.path.join(base,file)
pred = pd.read_csv(url, sep='\t|,', engine='python')
pred = column_name_mapping(pred)
cols = pred.columns.values.tolist()
if not pred[cols[0]].map(str).values[0].startswith('chr'):
    pred[cols[0]] = pred[cols[0]].map(str).replace('^','chr',regex=True)


result = circadian_ai_enrichment(pred,snps,0.5,5)
print(result)


with open(OUTPUT, 'w') as f:
    f.write(result)
