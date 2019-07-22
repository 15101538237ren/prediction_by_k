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
DATA_TYPEs = ['MIN', 'MAX', 'METHY']
K_LEN = [1, 3, 5]
INDEX_NAME = 'Indexs'
INDEX_DIR = os.path.join(DATA_SET_DIR, INDEX_NAME)
CHROMOSOMES = [str(i) for i in range(1, 23)]

K_KEY = 'k'
D_KEY = 'd'
ND_KEY = 'nd'
ANNO_KEY = 'anno'
ANNO_ONE_HOT_KEY = 'anno_one_hot'
NTYPE_ANNOTATIONS = 16
ANNO_SELECTION = ANNO_KEY
FIELD_ARR = [ [K_KEY], [K_KEY, D_KEY], [K_KEY, ND_KEY], [K_KEY, ANNO_SELECTION],[K_KEY, D_KEY, ANNO_SELECTION], [K_KEY, ND_KEY, ANNO_SELECTION]]
FIELD_LABELS = ['k', 'k + d', 'k + nd', 'k + ChrMM', 'k + d + ChrMM', 'k + nd + ChrMM']
METHY_KEY = 'METHY'
K_METHY_KEY = K_KEY + "+" + METHY_KEY
DATA_TYPE_ARR = [K_KEY, METHY_KEY, K_METHY_KEY]

EACH_SUB_FIG_SIZE = 5
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

    with open(input_fp, "w") as file:
        file.write((line_end).join(ltws))
        file.write(line_end)
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

def plot_pca_or_tsne(DATA_DIR, klen, tsne = 0):
    FIG_TYPE = "TSNE" if tsne else "PCA"
    FIELD_ARR_HERE = FIELD_ARR[1 : ] if klen == 1 else FIELD_ARR
    FIELD_LABELS_HERE = FIELD_LABELS[1 : ] if klen == 1 else FIELD_LABELS
    hbins =len(FIELD_ARR_HERE)
    wbins = len(DATA_TYPE_ARR)

    fig, axs = plt.subplots(hbins, wbins, figsize=(hbins * EACH_SUB_FIG_SIZE, wbins * (EACH_SUB_FIG_SIZE - 2)))

    for f_i, filed in enumerate(FIELD_ARR_HERE):
        for d_i, type_name in enumerate(DATA_TYPE_ARR):
            if K_KEY == type_name:
                data_type = 'MIN'
                x_train, _, y_train, _ = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
            elif METHY_KEY == type_name:
                data_type = METHY_KEY
                x_train, _, y_train, _ = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
            elif K_METHY_KEY == type_name:
                data_type = 'MIN'
                x_train1, _, y_train1, _ = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
                data_type = METHY_KEY
                x_train2, _, y_train2, _ = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
                x_train = np.concatenate((x_train1, x_train2), axis=0)
                y_train = np.concatenate((y_train1, y_train2), axis=0)
            else:
                return
            if tsne:
                time_start = time.time()
                pca = TSNE(n_components=2, perplexity=10, n_iter=250)
                pc = pca.fit_transform(x_train)
                print('t-SNE in %.2f seconds' % (time.time() - time_start))
            else:
                x_train = StandardScaler().fit_transform(x_train)
                pca = PCA(n_components=2)
                pc = pca.fit_transform(x_train)
            ax = axs[f_i, d_i]
            for target, color in zip(CLASSES, COLORS):
                indexs = y_train == target
                ax.scatter(pc[indexs, 0], pc[indexs, 1], c=color, s=10)
            if not tsne:
                var_exp = [str(round(item, 2)) for item in pca.explained_variance_ratio_[0 : 2]]
                var_str = "pc1:"+ var_exp[0] + ", pc2:" + var_exp[1]
                ax.text(0, 0, var_str , fontsize=10, color="red")
            ax.set_xlabel('PC1', fontsize=12)
            if d_i == 0:
                ax.set_ylabel(FIELD_LABELS_HERE[f_i] , fontsize=12)
            if f_i == 0:
                ax.set_title(type_name, fontsize=12)
            if f_i ==0 and d_i == 0:
                ax.legend(CLASSES_LABELS)
            ax.grid()

    plt.savefig("../FIGURES/" + FIG_TYPE + "_Klen_" + str(klen) + ".png", dpi=200)
    plt.show()

def plot_figures():
    OUT_DATA_DIR = '../DATA/MINI_DATA_SET'
    for tsne in [0, 1]:
        for klen in K_LEN:
            plot_pca_or_tsne(OUT_DATA_DIR, klen, tsne=tsne)

def build_classifiers():
    classifier_names = [
        #"LR",
        #"KNN",
        "MLP",
        #"GBDT"
    ]

    classifiers = [
        #LogisticRegression(C=1, penalty='l1', max_iter=10000),
        # KNeighborsClassifier(n_neighbors=1)#,
        MLPClassifier(alpha=1e-5, hidden_layer_sizes=(100, 100, 20), random_state=RANDOM_STATE, max_iter=50000, learning_rate='adaptive')
        #GradientBoostingClassifier(n_estimators=5, learning_rate=.1, max_features=2, max_depth=2,random_state=RANDOM_STATE)
    ]
    return [classifier_names, classifiers]
def generate_table(out_dir, klen, performances):
    FIELD_LABELS_HERE = FIELD_LABELS[1:] if klen == 1 else FIELD_LABELS
    performances = np.array(performances)
    metrics = ['accuarcy', 'recall', 'precision', 'f1', 'roc']
    for c_i, classifier_name in enumerate(classifier_names):
        for m_i, metric in enumerate(metrics):
            out_fp = os.path.join(out_dir, str(klen) + '_' + classifier_name+ "_" + metric + ".csv")
            performance = performances[:,:, c_i, m_i]
            row_names = np.array(FIELD_LABELS_HERE).reshape((performance.shape[0]), 1)
            performance = np.concatenate((row_names, performance), axis=1)
            np.savetxt(out_fp, performance[:, :], fmt="%s,%s,%s,%s", delimiter='\n', header= ',' + ','.join(DATA_TYPE_ARR))

[classifier_names, classifiers] = build_classifiers()

def prediction(klen, predict_type, X_train, X_test, Y_train,Y_test):
    performances = []
    for cidx, clf_name in enumerate(classifier_names):
        clf = classifiers[cidx].fit(X_train, Y_train)
        y_pred = clf.predict(X_test)
        if hasattr(clf, "decision_function"):
            Z = clf.decision_function(X_test)
        else:
            Z = clf.predict_proba(X_test)[:, 1]
        fpr_gb, tpr_gb, _ = roc_curve(Y_test, Z)
        roc = auc(fpr_gb, tpr_gb)
        perf = [accuracy_score(Y_test, y_pred),
                       recall_score(Y_test, y_pred),
                       precision_score(Y_test, y_pred),
                       f1_score(Y_test, y_pred),
                       roc]
        print("%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (str(klen), predict_type,
                                                            clf_name, accuracy_score(Y_test, y_pred),
                                                            recall_score(Y_test, y_pred),
                                                            precision_score(Y_test, y_pred),
                                                            f1_score(Y_test, y_pred),
                                                            roc))
        performances.append(perf)
    return performances
def machine_learning(DATA_DIR, klen):
    FIELD_ARR_HERE = FIELD_ARR[1 : ] if klen == 1 else FIELD_ARR
    FIELD_LABELS_HERE = FIELD_LABELS[1 : ] if klen == 1 else FIELD_LABELS
    performances = [[] for item in range(len(FIELD_ARR_HERE))]

    for f_i, filed in enumerate(FIELD_ARR_HERE):
        for d_i, type_name in enumerate(DATA_TYPE_ARR):
            if K_KEY == type_name:
                data_type = 'MIN'
                x_train, x_test, y_train, y_test = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
            elif METHY_KEY == type_name:
                data_type = METHY_KEY
                x_train, x_test, y_train, y_test= obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
            elif K_METHY_KEY == type_name:
                data_type = 'MIN'
                x_train1, x_test1, y_train1, y_test1 = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
                data_type = METHY_KEY
                x_train2, x_test2, y_train2, y_test2 = obtain_data_by_field_selection(DATA_DIR, klen, data_type, filed)
                x_train = np.concatenate((x_train1, x_train2), axis=0)
                y_train = np.concatenate((y_train1, y_train2), axis=0)
                x_test = np.concatenate((x_test1, x_test2), axis=0)
                y_test = np.concatenate((y_test1, y_test2), axis=0)
            else:
                return
            y_train[y_train > 0] = 1
            y_test[y_test > 0] = 1
            predict_type = FIELD_LABELS_HERE[f_i] + " " + type_name
            perf = prediction(klen, predict_type, x_train, x_test, y_train, y_test)
            performances[f_i].append(perf)
    return performances
if __name__ == "__main__":
    # OUT_DIR = '../DATA/MINI_DATA_SET'
    # PERF_DIR='../DATA/ML_PERFORMANCE'
    # for klen in K_LEN:
    #     performances = machine_learning(OUT_DIR, klen)
    #     generate_table(PERF_DIR, klen, performances)
    # #
    # # create_mini_data_set(OUT_DIR)
    # preprocess_ChromeHMM()
    pass