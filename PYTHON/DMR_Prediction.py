# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os, re, pickle
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
#%matplotlib inline

#1. PCA Analysis
DATA_SET_DIR = '../DATA/DATA_SET'
FILE_END = '.bed.csv'
DATA_TYPE = {'MAX': '_K_max', 'MIN': '_K_min', 'METHY': '_methylation'}
K_LEN = [3]
INDEX_NAME = 'Indexs'
INDEX_DIR = os.path.join(DATA_SET_DIR, INDEX_NAME)
CHROMOSOMES = [str(i) for i in range(1, 23)]
CLASSES = [0,1,2,3]
CLASSES_LABELS = ['NON', 'Ecto DMR', 'Endo DMR', 'Meso DMR']
COLORS = ['r', 'g', 'b', 'k']


def create_mini_data_set(OUT_DIR, full_dataset=False):
    OUT_IDNEX_DIR = os.path.join(OUT_DIR, INDEX_NAME)
    if not os.path.exists(OUT_IDNEX_DIR):
        os.makedirs(OUT_IDNEX_DIR)

    for k_len in K_LEN:
        for key, item in DATA_TYPE.items():
            for chr_i in CHROMOSOMES:
                INDEX_FILE = os.path.join(INDEX_DIR, 'chr' + chr_i + ".csv")
                index_array = pd.read_csv(INDEX_FILE, sep=',', header=None).values.astype(int)
                if full_dataset:
                    DATA_SIZE = index_array.shape[0]
                    sampled_indexs = range(0, DATA_SIZE)
                else:
                    non_zero_indexs = np.nonzero(index_array[:, 2])[0]
                    zero_indexs = np.where(index_array[:, 2] == 0)[0]
                    DATA_SIZE = non_zero_indexs.size

                    sampled_nz_indexs =  np.random.choice(non_zero_indexs, DATA_SIZE, replace=False)
                    sampled_zero_indexs = np.random.choice(zero_indexs, DATA_SIZE, replace=False)
                    sampled_indexs = np.concatenate((sampled_nz_indexs, sampled_zero_indexs), axis=0)

                sampled_arr = index_array[sampled_indexs]
                out_indexs_fp = os.path.join(OUT_IDNEX_DIR , 'chr' + chr_i + ".csv")
                np.savetxt(out_indexs_fp, sampled_arr[:, :], fmt="%d,%d,%d", delimiter='\n')
                k_len_s = str(k_len)
                print("processing with %s %s chr %s " % (k_len_s, key, chr_i))

                input_fp = os.path.join(DATA_SET_DIR, k_len_s, k_len_s + item, 'chr' + chr_i + ".csv")
                data_pd = pd.read_csv(input_fp, sep=',', header=None).values.astype(float)
                sampled_data = data_pd[sampled_indexs]
                sampled_data = np.concatenate((sampled_data, sampled_arr[:,2].reshape(-1,1)), axis=1)
                OUT_DATA_DIR = os.path.join(OUT_DIR, k_len_s, k_len_s + item)
                if not os.path.exists(OUT_DATA_DIR):
                    os.makedirs(OUT_DATA_DIR)
                output_fp = os.path.join(OUT_DATA_DIR,'chr' + chr_i + ".csv")
                filed_count = 2 if key == 'METHY' else 3
                fmt_middle = ','.join(['%.4f' for i in range(filed_count * k_len)])
                np.savetxt(output_fp, sampled_data[:, :], fmt="%d,%d,"+fmt_middle+",%d", delimiter='\n')
def merged_diff_genomic_features(sep="\s+",line_end = "\n"):
    in_dir = '../DATA/Genomic_Features'
    for file_name in os.listdir(in_dir):
        if file_name.endswith(".bed"):
            ltws = []
            file_fp = os.path.join(in_dir, file_name)
            print("processing with %s" % file_name)
            with open(file_fp, "r") as file:
                lines = [line for line in file]

                for line_idx, line in enumerate(lines):
                    chr_i, start, end = re.split(sep, line.strip(line_end))[0:3]
                    ltw = "\t".join([chr_i, start, end])if int(start) < int(end) else "\t".join([chr_i, end, start])
                    ltws.append(ltw)

            with open(file_fp, "w") as file:
                file.write((line_end).join(ltws))
                file.write(line_end)

def pca_analysis(data_dir, k_len, data_type, chr_i):
    input_fp = os.path.join(data_dir, str(k_len), str(k_len) + DATA_TYPE[data_type], 'chr' + chr_i + ".csv")
    data_pd = pd.read_csv(input_fp, sep=',', header=None).values.astype(float)
    S_COL = 2
    COL_INDEXS = [i for i in range(S_COL, S_COL + 2* k_len)] if data_type == 'METHY' else list(np.array([[S_COL, S_COL + 2 * k_len] for i in range(0,k_len)]).flatten())
    x = data_pd[:, COL_INDEXS]
    x = StandardScaler().fit_transform(x)
    y = data_pd[:, -1]
    pca = PCA(n_components=2)
    pc = pca.fit_transform(x)
    return [pc, y]
def plot_pca(pc, y):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('PC1', fontsize=15)
    ax.set_ylabel('PC2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    for target, color in zip(CLASSES, COLORS):
        indexs = y == target
        ax.scatter(pc[indexs, 0], pc[indexs, 1], c=color, s=50)
    ax.legend(CLASSES_LABELS)
    ax.grid()
if __name__ == "__main__":
    OUT_DIR = '../DATA/MINI_DATA_SET'
    [pc, y] = pca_analysis(OUT_DIR, 1, 'MIN', '1')
    plot_pca(pc, y)