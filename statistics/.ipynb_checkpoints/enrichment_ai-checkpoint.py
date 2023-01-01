#!/usr/bin/env python
# enrichment_ai.py
"""""""""""""""""""""""""""
# Author: Keila Velazquez-Arcelay
# 
# Description: Tests for overrepresentation of a set of variants in regions predicted 
# to contain adaptive introgression. Generates randum null sets of size equal to the 
# number of significant AI prediction regions.  Intersects each random set with the 
# set of observed variants to find the expected background values. Gets sum of regions 
# containing at least 1 circadian introgressed SNP in the set of circadian AI regions 
# and all the random iterations. Overlapping AI regions are merged into one region.
#
"""

import argparse
import pandas as pd
import numpy as np
import pybedtools


def genomatnn():
    pred = pd.read_csv(PRED, sep='\t|,', engine='python')
    pred = column_name_mapping(pred)
    cols = pred.columns.values.tolist()
    if not pred[cols[0]].map(str).values[0].startswith('chr'):
        pred[cols[0]] = pred[cols[0]].map(str).replace('^','chr',regex=True)
    
    
def maladapt():
    
    
    
def column_name_mapping(df):
    mapping = {df.columns[0]: 'chr', 
               df.columns[1]: 'start', 
               df.columns[2]: 'end'}
    df = df.rename(columns=mapping)
    return df


def merge_overlapping_regions(df):
    temp_df = pd.DataFrame()
    for i in df['chr'].drop_duplicates():
        temp = df.copy()
        temp = temp[temp['chr'].isin([i])]
        temp['group'] = ((temp['start']>temp['end'].shift())&\
                         (temp['end']>temp['start'].shift())).cumsum()
        # Join regions in each chromosome
        temp_df = pd.concat([temp_df, temp.groupby(['group']).agg({'start':'min', 'end': 'max'})])
    # Add chromosome names
    temp_df = pd.merge(sig_regions,temp_df,on='start').iloc[:,[0,1,-1]]
    temp_df.columns = temp_df.columns.str.replace(r'\_.*',r'')
    return temp_df


def intersect_a_and_b(df_a,df_b):
    df_a = df_a.add_suffix('_wa')
    df_b = df_b.add_suffix('_wb')
    df_a_cols = df_a.columns.values.tolist()
    df_b_cols = df_b.columns.values.tolist()
    a = pybedtools.BedTool.from_dataframe(df_a)
    b = pybedtools.BedTool.from_dataframe(df_b)
    a_and_b = a.intersect(b, wa=True, wb=True)
    a_and_b_df = pd.read_table(a_and_b.fn, names=df_a_cols+df_b_cols)
    #a_and_b_df = a_and_b_df[df_a_cols+[df_b_cols[-1]]]
    return a_and_b_df


# ============= PARSER =============

parser = argparse.ArgumentParser(description='.')
parser.add_argument('method', action='store', type=str, help='Method used to predict AI')
parser.add_argument('predictions', action='store', type=str, help='Input file')
parser.add_argument('snps', action='store', type=str, help='Input file')
parser.add_argument('threshold', action='store', type=float, help='Threshold')
parser.add_argument('iterations', action='store', type=int, help='Iterations')
parser.add_argument('-o', action='store', type=str, dest='output',
                    default='output.txt', help='Output file')


# =========== ARGUMENTS ===========
passed_args = parser.parse_args()
METHOD = passed_args.method
pred = passed_args.predictions
SNPS = passed_args.snps
THRESHOLD = passed_args.threshold
K = passed_args.iterations
OUTPUT = passed_args.output


# LOAD DATA
if METHOD == 'genomatnn'.lower():
    pred = genomatnn()
elif METHOD == 'maladapt'.lower():
    pred = maladapt()
else:
    break
    

snps = pd.read_csv(SNPS, sep='\t').iloc[:,[0,1,2]].drop_duplicates()
snps = column_name_mapping(snps)


# Extract AI regions
sig_regions = pred[pred.iloc[:,-1]>=THRESHOLD].iloc[:,:3].drop_duplicates()
sig_regions.sort_values(by=sig_regions.columns[[0,1]].values.tolist(),inplace=True)

ai_regions = merge_overlapping_regions(sig_regions)


# FIND NUMBER OF AI REGIONS CONTAINING AT LEAST 1 CIRCADIAN INTROGRESSED SNP
observed_ai_region_counts = 0
for i in range(len(ai_regions)):
    obs_intersection = intersect_a_and_b(snps,ai_regions.iloc[i:i+1,:]).iloc[:,:3].drop_duplicates()
    if len(obs_intersection) > 0: observed_ai_region_counts += 1
observed_mean = np.mean(observed_ai_region_counts)


# GENERATE RANDOM SETS
size = len(sig_regions)
rand_sets = [pred.iloc[:,:3].drop_duplicates().sample(n=size) for i in range(K)]


# INTERSECT EACH RANDOM SET WITH CIRCADIAN INTROGRESSED VARIANTS
regions_by_iter = []
for i in rand_sets:
    exp_intersection = intersect_a_and_b(snps,i)
    # Get circadian introgressed variants in set i
    expected_region_counts = 0
    for j in range(len(i)):
        # Get circadian introgressed variants in region j
        cols_intersect = exp_intersection.columns.values.tolist()
        cols_i = i.columns.values.tolist()
        merge = pd.merge(exp_intersection,i.iloc[j:j+1,:],
                         left_on=cols_intersect[3:6], 
                         right_on=cols_i[0:3]).iloc[:,0:3].drop_duplicates()
        if len(merge) > 0: expected_region_counts += 1
    regions_by_iter.append(expected_region_counts)
rand_mean = np.mean(regions_by_iter)


# PRINT MEANS
print('Observed regions: {}'.format(observed_ai_region_counts))
print('Expected average regions: {}'.format(rand_mean))


# ODDS RATIO
print('obs/exp: ', observed_ai_region_counts / rand_mean)


# P-VALUE
n = 0
for counts in regions_by_iter:
    if counts >= observed_ai_region_counts: n +=1
print('p-val: ', n/K)


# SAVE RESULTS
results = '\nRegions: {} \nObserved regions: {} \nExpected regions: {} \nOR: {} \nP: {} \n'.format(
    size,observed_ai_region_counts,rand_mean,observed_ai_region_counts/rand_mean,n/K)
with open(OUTPUT, 'w') as f:
    f.write(results)