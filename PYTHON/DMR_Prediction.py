# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os, re, time
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier

from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc, accuracy_score, recall_score, precision_score, f1_score

import matplotlib.pyplot as plt
#%matplotlib inline

#1. PCA Analysis
DATA_SET_DIR = '../DATA/DATA_SET'
FILE_END = '.bed.csv'
DATA_TYPE = {'MAX': '_K_max', 'MIN': '_K_min', 'METHY': '_methylation'}
DATA_TYPE_ARR = ['MIN', 'MAX', 'METHY']
K_LEN = [1, 3, 5]
INDEX_NAME = 'Indexs'
INDEX_DIR = os.path.join(DATA_SET_DIR, INDEX_NAME)
CHROMOSOMES = [str(i) for i in range(1, 23)]
CLASSES = [0, 1, 2, 3]
CLASSES_LABELS = ['NON', 'Ecto DMR', 'Endo DMR', 'Meso DMR']
COLORS = ['r', 'g', 'b', 'k']
#COLORS = ['b', 'b', 'b', 'b']
RANDOM_STATE = 42
TEST_RATIO = 0.33

def create_mini_data_set(OUT_DIR, full_dataset=False):
    OUT_IDNEX_DIR = os.path.join(OUT_DIR, INDEX_NAME)
    if not os.path.exists(OUT_IDNEX_DIR):
        os.makedirs(OUT_IDNEX_DIR)

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

            sampled_nz_indexs = np.random.choice(non_zero_indexs, DATA_SIZE, replace=False)
            sampled_zero_indexs = np.random.choice(zero_indexs, DATA_SIZE, replace=False)
            sampled_indexs = np.concatenate((sampled_nz_indexs, sampled_zero_indexs), axis=0)

        for k_len in K_LEN:
            for key, item in DATA_TYPE.items():
                sampled_arr = index_array[sampled_indexs]
                out_indexs_fp = os.path.join(OUT_IDNEX_DIR , 'chr' + chr_i + ".csv")
                np.savetxt(out_indexs_fp, sampled_arr[:, :], fmt="%d,%d,%d,%d", delimiter='\n')
                k_len_s = str(k_len)
                print("processing with %s %s chr %s " % (k_len_s, key, chr_i))

                input_fp = os.path.join(DATA_SET_DIR, k_len_s, k_len_s + item, 'chr' + chr_i + ".csv")
                data_pd = pd.read_csv(input_fp, sep=',', header=None).values.astype(float)
                sampled_data = data_pd[sampled_indexs]
                sampled_data = np.concatenate((sampled_data[:, 0: -1], sampled_arr[:, 2:4].reshape(-1,2)), axis=1)
                OUT_DATA_DIR = os.path.join(OUT_DIR, k_len_s, k_len_s + item)
                if not os.path.exists(OUT_DATA_DIR):
                    os.makedirs(OUT_DATA_DIR)
                output_fp = os.path.join(OUT_DATA_DIR,'chr' + chr_i + ".csv")
                filed_count = 2 if key == 'METHY' else 3
                fmt_middle = ','.join(['%.4f' for i in range(filed_count * k_len)])
                np.savetxt(output_fp, sampled_data[:, :], fmt="%d,%d,"+fmt_middle+",%d,%d", delimiter='\n')

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

def preprocess_ChromeHMM(sep="\s+",line_end = "\n"):
    input_fp = '../DATA/Genomic_Features/HmmH1hescHMM.bed'
    ltws = []
    ChromeMM_DICT = {}
    with open(input_fp, "r") as file:
        lines = [line for line in file]

        for line_idx, line in enumerate(lines):
            chr_i, start, end, label = re.split(sep, line.strip(line_end))[0:4]
            label_id = label.split('_')[0]
            if str(label_id) not in ChromeMM_DICT.keys():
                ChromeMM_DICT[str(label_id)] = ' '.join(label.split('_')[1:])
            ltw = "\t".join([chr_i, start, end, label_id])if int(start) < int(end) else "\t".join([chr_i, end, start, label_id])
            ltws.append(ltw)

    # with open(input_fp, "w") as file:
    #     file.write((line_end).join(ltws))
    #     file.write(line_end)
    for label_id, label in ChromeMM_DICT.items():
        print("%s\t%s" %(label_id, label))

def read_single_csv_data(data_dir, k_len, data_type, chr_i):
    input_fp = os.path.join(data_dir, str(k_len), str(k_len) + DATA_TYPE[data_type], 'chr' + chr_i + ".csv")
    data_pd = pd.read_csv(input_fp, sep=',', header=None).values.astype(float)
    S_COL = 2
    COL_INDEXS = [i for i in range(S_COL, S_COL + 2* k_len)] if data_type == 'METHY' else list(np.array([[S_COL, S_COL + 2 * k_len] for i in range(0,k_len)]).flatten())
    COL_INDEXS += [-1]
    x = data_pd[:, COL_INDEXS].astype(float)
    y = data_pd[:, -2].astype(int)
    return [x, y]
def read_data(data_set_dir, k_len, data_type):
    c_x1, c_x2, c_y1, c_y2 = np.array([]), np.array([]), np.array([]), np.array([])
    for chr_i in CHROMOSOMES:
        x, y = read_single_csv_data(data_set_dir, k_len, data_type, chr_i)
        half_len = int(x.shape[0] / 2.0)
        train_end1 = int(half_len * (1.0 - TEST_RATIO))
        train_end2 = half_len + train_end1
        train_indexs = [i for i in range(0, train_end1)] + [i for i in range(half_len, train_end2)]
        test_indexs = [i for i in range(train_end1, half_len)] + [i for i in range(train_end2, x.shape[0])]

        c_x1 = np.concatenate((c_x1, x[train_indexs, :]), axis=0) if c_x1.size != 0 else x[train_indexs, :]
        c_y1 = np.concatenate((c_y1, y[train_indexs]), axis=0) if c_y1.size != 0 else y[train_indexs]

        c_x2 = np.concatenate((c_x2, x[test_indexs, :]), axis=0) if c_x2.size != 0 else x[test_indexs, :]
        c_y2 = np.concatenate((c_y2, y[test_indexs]), axis=0) if c_y2.size != 0 else y[test_indexs]
    return [c_x1, c_x2, c_y1, c_y2]
K_KEY = 'k'
D_KEY = 'd'
ND_KEY = 'nd'
ANNO_KEY = 'anno'
ANNO_ONE_HOT_KEY = 'anno_one_hot'
NTYPE_ANNOTATIONS = 16
def obtain_data_by_field_selection(data_set_dir, k_len, data_type, fields):
    data_list = read_data(data_set_dir, k_len, data_type)
    COL_INDEXS = []
    for item in fields:
        if item == K_KEY:
            for idx in range(0, k_len):
                COL_INDEXS.append(idx)
        elif item == D_KEY or item == ND_KEY:
            for idx in range(k_len, 2*k_len):
                COL_INDEXS.append(idx)
        elif item == ANNO_KEY or item == ANNO_ONE_HOT_KEY:
            COL_INDEXS.append(-1)
    for didx, data in enumerate(data_list[0: 2]): # do not change y label
        if ND_KEY in fields:
            target_cols = [idx for idx in range(k_len, 2*k_len)]
            data[:, target_cols] = 1.0/(1.0+np.absolute(data[:, target_cols]))
        data = data[:, COL_INDEXS]
        if ANNO_ONE_HOT_KEY in fields:
            target_cols = -1
            annotations = data[:, target_cols].astype(int)
            anno_len = annotations.shape[0]
            one_hot_annotations = np.zeros((anno_len, NTYPE_ANNOTATIONS))
            one_hot_annotations[np.arange(anno_len), annotations] = 1
            data_list[didx] = np.concatenate((data[:, 0 : -1], one_hot_annotations), axis=1)
    return data_list
# def read_data():
#     x_train = {}
#     y_train = {}
#     x_test = {}
#     y_test = {}
#     OUT_DIR = '../DATA/MINI_DATA_SET'
#     for k_i, klen in enumerate(K_LEN):
#         for d_i, datatype in enumerate(DATA_TYPE_ARR):
#
#             dict_key = str(klen) + '_' + datatype
#             x_train[dict_key] = c_x1
#             y_train[dict_key] = c_y1
#             x_test[dict_key] = c_x2
#             y_test[dict_key] = c_y2
#     return [x_train, y_train, x_test, y_test]
#
# [x_train, y_train, x_test, y_test] = read_data()
obtain_data_by_field_selection('../DATA/MINI_DATA_SET', 3, 'MIN', [K_KEY, ND_KEY, ANNO_KEY])

def plot_pca():
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    for k_i, klen in enumerate(K_LEN):
        for d_i, datatype in enumerate(DATA_TYPE_ARR):
            dict_key = str(klen) + '_' + datatype
            x = x_train[dict_key]
            x = StandardScaler().fit_transform(x)
            pca = PCA(n_components=2)
            pc = pca.fit_transform(x)
            y = y_train[dict_key]
            ax = axs[k_i, d_i]
            for target, color in zip(CLASSES, COLORS):
                indexs = y == target
                ax.scatter(pc[indexs, 0], pc[indexs, 1], c=color, s=10)

            ax.set_xlabel('PC1', fontsize=8)
            ax.set_ylabel('K,d:' + str(klen) + ' PC2', fontsize=8)
            ax.set_title(datatype, fontsize=8)
            ax.legend(CLASSES_LABELS)
            ax.grid()

    plt.savefig('../FIGURES/PCA.png', dpi=200)
    plt.show()

def plot_tsne():
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    for k_i, klen in enumerate(K_LEN):
        for d_i, datatype in enumerate(DATA_TYPE_ARR):
            dict_key = str(klen) + '_' + datatype
            x = x_train[dict_key]
            time_start = time.time()
            tsne = TSNE(n_components=2, perplexity=10, n_iter=250)
            pc = tsne.fit_transform(x)
            print('t-SNE in %.2f seconds'%(time.time() - time_start))
            y = y_train[dict_key]
            ax = axs[k_i, d_i]
            for target, color in zip(CLASSES, COLORS):
                indexs = y == target
                ax.scatter(pc[indexs, 0], pc[indexs, 1], c=color, s=10)

            ax.set_xlabel('TSNE 1', fontsize=8)
            ax.set_ylabel('K,d:' + str(klen) + ' TSNE 2', fontsize=8)
            ax.set_title(datatype, fontsize=8)
            ax.legend(CLASSES_LABELS)
            ax.grid()

    plt.savefig('../FIGURES/TSNE.png', dpi=200)
    plt.show()

def build_classifiers():
    classifier_names = [
        #"LR",
        #"KNN",
        "MLP",
        #"GBDT"
    ]

    classifiers = [
        #LogisticRegression(C=1, penalty='l1', max_iter=10000),
        #KNeighborsClassifier(n_neighbors=1)#,
        MLPClassifier(alpha=1e-5, hidden_layer_sizes=(100, 100, 20), random_state=RANDOM_STATE, max_iter=50000, learning_rate='adaptive')
        #GradientBoostingClassifier(n_estimators=5, learning_rate=.1, max_features=2, max_depth=2,random_state=RANDOM_STATE)
    ]
    return [classifier_names, classifiers]

def prepare_ml_data(dict_key):
    scaler = MinMaxScaler()
    X_train = x_train[dict_key]
    X_test = x_test[dict_key]

    Y_train = y_train[dict_key]
    Y_test = y_test[dict_key]
    Y_train[Y_train > 0] = 1
    Y_test[Y_test > 0] = 1

    X_train = scaler.fit_transform(X_train)
    X_test = scaler.fit_transform(X_test)
    return X_train, X_test, Y_train,Y_test
def prediction(klen, predict_type, X_train, X_test, Y_train,Y_test):
    [classifier_names, classifiers] = build_classifiers()
    for cidx, clf_name in enumerate(classifier_names):
        clf = classifiers[cidx].fit(X_train, Y_train)
        y_pred = clf.predict(X_test)
        if hasattr(clf, "decision_function"):
            Z = clf.decision_function(X_test)
        else:
            Z = clf.predict_proba(X_test)[:, 1]
        fpr_gb, tpr_gb, _ = roc_curve(Y_test, Z)
        roc = auc(fpr_gb, tpr_gb)

        print("%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (str(klen), predict_type,
                                                            clf_name, accuracy_score(Y_test, y_pred),
                                                            recall_score(Y_test, y_pred),
                                                            precision_score(Y_test, y_pred),
                                                            f1_score(Y_test, y_pred),
                                                            roc))
def machine_learning():
    X_train = {}
    X_test = {}
    Y_train = {}
    Y_test = {}

    for k_i, klen in enumerate(K_LEN):
        for d_i, datatype in enumerate(DATA_TYPE_ARR):
            dict_key = str(klen) + '_' + datatype
            X_train[dict_key], X_test[dict_key], Y_train[dict_key], Y_test[dict_key] = prepare_ml_data(dict_key)

    for k_i, klen in enumerate(K_LEN):
        pt1 = 'K'
        key1 = 'MIN'
        dk1 = str(klen) + '_' + key1
        prediction(klen, pt1, X_train[dk1], X_test[dk1], Y_train[dk1], Y_test[dk1])

        pt2 = 'METHY'
        key2 = pt2
        dk2 = str(klen) + '_' + key2
        prediction(klen, pt2, X_train[dk2], X_test[dk2], Y_train[dk2], Y_test[dk2])

        pt3 = 'K + METHY'
        X_train_c = np.concatenate((X_train[dk1], X_train[dk2]), axis=1)
        X_test_c = np.concatenate((X_test[dk1], X_test[dk2]), axis=1)
        prediction(klen, pt3, X_train_c, X_test_c, Y_train[dk1], Y_test[dk1])

if __name__ == "__main__":
    machine_learning()
    # OUT_DIR = '../DATA/MINI_DATA_SET'
    # create_mini_data_set(OUT_DIR)