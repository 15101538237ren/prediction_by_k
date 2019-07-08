# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os, re, pickle

#1. PCA Analysis
DATA_SET_DIR = '../DATA/DATA_SET'
FILE_END = '.bed.csv'
DATA_TYPE = {'MAX': '_K_max', 'MIN': '_K_min', 'METHY': '_methylation'}
K_LEN = [1, 3, 5]

INDEX_FILE = os.path.join(DATA_SET_DIR, 'dmr_merged_tile.bed.index')

index_array = pd.read_csv(INDEX_FILE, sep=',').values

non_zero_indexs = np.nonzero(index_array[:, 2])
non_zero_arr = index_array[non_zero_indexs]
non_zero_size = non_zero_indexs[0].size

zero_indexs = np.where(index_array[:, 2] == 0)[0]
sample_zero_indexs = np.random.choice(zero_indexs, non_zero_size, replace=False)
zero_arr= index_array[sample_zero_indexs]
def create_mini_data_set():
    pass
