import glob
import os
import pandas as pd
import re
import numpy as np
from matplotlib import pyplot as plt
import math
import seaborn as sns

# read in each _readcounts.txt as a dataframe
# store in a dict of dataframes w/ IDs as keys

root_path = '/labs/mbarna/users/tsusanto/ribosome_profiling/2_hiNd5-RNARFP-dep_20230511/stage9'

read_count_colnames = ['tx_id_v',
                       '5UTR',
                       'first15codons',
                       'CDS',
                       'last5codons',
                       '3UTR']

# import mart to label with gene names

mart = pd.read_csv('/labs/mbarna/users/adelexu/s27l-flag-riboprof/biomart/mart_jan2020_gene_transcript_name.txt', 
                   sep='\t', skiprows=1,
                   names=['gene_id', 'tx_id', 'tx_id_v', 'gene_name'])
name_for_tx = dict(zip(mart.tx_id, mart.gene_name))

# read in sample data

sample_dfs = {}

for source_path in glob.glob(os.path.join(root_path, '*', '*_readcounts.txt')):
    samp = os.path.basename(source_path).replace('_readcounts.txt', '').replace('20200829_lane', 'L')
    sample_dfs[samp] = pd.read_csv(source_path,
                                   sep='\t',
                                   header=None,
                                   names=read_count_colnames)
    sample_dfs[samp]['tx_id'] = [re.sub(r'\.\d*', '', tidv) for tidv in sample_dfs[samp]['tx_id_v']]
    sample_dfs[samp].set_index('tx_id', inplace=True)
    sample_dfs[samp].drop(columns='tx_id_v', inplace=True)
    sample_dfs[samp]['reads_per_tx'] = sample_dfs[samp].sum(axis=1)
    sample_dfs[samp]['rpm_per_tx'] = sample_dfs[samp]['reads_per_tx']/sample_dfs[samp]['reads_per_tx'].sum()*10**6

# CAUTION: _readcounts.txt rows cannot be compared directly by index -- different order of txs!!!
