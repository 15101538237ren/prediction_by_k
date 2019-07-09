# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os, re, pickle
from itertools import combinations

from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, recall_score, precision_score, f1_score

#Models for classifiers
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neural_network import MLPClassifier

TARGET_COLUMN_INDEX_OF_BED_FILE = 3
TEST_RATIO = 0.33
RANDOM_STATE = 42
WINDOW_SZS = [1]#1, 3, 5,7, 11, 23, 51,
modification_names = ["H3K4me1", "H3K4me3", "H3K9Me3", "H3K27ac", "H3K27Me3", "H3K36me3","Methylation"]#,,
differentiated_cell_types = ["dEC"]#, "dME", "dEN"
SCORING = ['f1_weighted']#'accuracy', 'precision', 'recall', ,'roc_auc'
SCORING_KEYS_FOR_PRINT = ['test_f1_weighted']#,'test_accuracy','test_precision', 'test_recall', 'test_f1_weighted','test_roc_auc'
DATA_DIR = '../DATA'

def read_tab_seperated_file_and_get_target_column(input_file_path, target_col_index, target_data_type=float, start_line= 1, sep="\s+",line_end = "\n"):
    """
    OBTAIN THE VALUES OF A TARGET COLUMN IN TSV FILE
    :param target_col_index: the target column index, starting from 0
    :param input_file_path: the path of the input tsv file
    :param start_line: the index of the start line in the file
    :param sep: column separator
    :param line_end: the line separator
    :return: the values of target column in tsv file
    """
    ret_value_list = []
    line_counter = 0

    with open(input_file_path, "r") as input_file:
        line = input_file.readline()
        while line:
            line_counter += 1
            if line_counter >= start_line:
                line_contents = re.split(sep, line.strip(line_end))
                try:
                    led = target_data_type(line_contents[target_col_index])
                    ret_value_list.append(led)
                except ValueError as e:
                    ret_value_list.append(0)
            line = input_file.readline()

    return ret_value_list

def combine_whole_dataset(dataset_dir, load=False):

    if load:
        with open(os.path.join(dataset_dir, 'DATA.pkl'), 'rb') as pkl_file:
            [k_data, histone_and_methy_data] = pickle.load(pkl_file)
    else:
        histone_and_methy_data = {}

        k_data = read_tab_seperated_file_and_get_target_column(os.path.join(dataset_dir, 'K_mean_in_each_peak.bed'),
                                                                   TARGET_COLUMN_INDEX_OF_BED_FILE)

        df = pd.read_csv(os.path.join(DATA_DIR, 'DATA_LABEL.txt'), sep='\s+', index_col=0)

        for idx in df.index:
            idx_s = str(idx)
            histone_and_methy_data[idx_s] = {}
            for col in df.columns:
                histone_and_methy_data[idx_s][col] = {}
                histone_and_methy_data[idx_s][col]['NAME'] = df.loc[idx, col]
                print(df.loc[idx, col])
                histone_and_methy_data[idx_s][col]['DATA'] = read_tab_seperated_file_and_get_target_column(os.path.join(dataset_dir, df.loc[idx, col] + '.bed'), TARGET_COLUMN_INDEX_OF_BED_FILE)
        with open(os.path.join(dataset_dir, 'DATA.pkl'), 'wb') as pkl_file:
            pickle.dump([k_data, histone_and_methy_data], pkl_file, -1)
    return [k_data, histone_and_methy_data]

def generate_the_combination_of_features():
    x_features = [] # the input feature for machine learning model
    y_features = []
    for target_modification_name in modification_names:

        left_features = [feature for feature in modification_names if feature != target_modification_name]
        for comb_feature_sz in range(0, 0):#len(left_features)
            comb_features = list(combinations(left_features, comb_feature_sz)) if comb_feature_sz else []

            if comb_feature_sz:
                for comb_feature in comb_features:
                    joined_features = [target_modification_name]
                    for feature in comb_feature:
                        joined_features.append(feature)
                    x_features.append(joined_features)
                    y_features.append([target_modification_name])
            else:
                x_features.append([target_modification_name])
                y_features.append([target_modification_name])

        x_features.append(["Methylation"])
        y_features.append([target_modification_name])
    return [x_features, y_features]

def sliding_window_of_data(x_data, y_data , window_sz):
    "Returns a sliding window (of width window_sz) over data (centered by (window_sz / 2 + 1)"
    mid = window_sz / 2

    sliced_data_length = len(x_data) - window_sz + 1
    sliced_x_data = np.zeros((sliced_data_length, window_sz * x_data.shape[1]))
    sliced_y_data = np.zeros((sliced_data_length, y_data.shape[1]))
    for i in range(sliced_data_length):
        sliced_x_data[i] = x_data[i : i + window_sz].flatten('F')
        sliced_y_data[i] = y_data[i + mid]
    return [sliced_x_data, sliced_y_data]
DATA_SET_DIR = os.path.join(DATA_DIR, 'DATA_FOR_PREDICTION')
original_dataset = combine_whole_dataset(DATA_SET_DIR, load=False)
k_data, histone_and_methy_data = original_dataset
hESC_data = histone_and_methy_data['hESC']

def sub_sampling(x_data, y_data, target_indexs, size = 1.0):
    x_under_represented = x_data[target_indexs]
    y_under_represented = y_data[target_indexs]

    x_over_represented_data = x_data[~target_indexs]
    y_over_represented_data = y_data[~target_indexs]

    selected_index = np.random.choice(len(x_over_represented_data), int(size * len(x_under_represented)))

    x_sampled_from_over_represented_data = x_over_represented_data[selected_index]
    y_sampled_from_over_represented_data = y_over_represented_data[selected_index]

    x_data = np.concatenate((x_under_represented, x_sampled_from_over_represented_data), axis=0)
    y_data = np.concatenate((y_under_represented, y_sampled_from_over_represented_data), axis=0)
    del x_under_represented, y_under_represented, x_over_represented_data, y_over_represented_data, x_sampled_from_over_represented_data, y_sampled_from_over_represented_data, selected_index
    x_data, y_data = shuffle(x_data, y_data, random_state=RANDOM_STATE)
    return [x_data, y_data]

def extract_target_data_by_features(differentiated_cell_type, x_features, y_features, include_k, window_sz):
    diff_cell_data = histone_and_methy_data[differentiated_cell_type]

    x_data = np.zeros((len(k_data), len(x_features) + include_k))
    for idx, feature in enumerate(x_features):
        x_data[:, idx] = hESC_data[feature]['DATA']
    if include_k:
        x_data[:, -1] = k_data

    y_data = np.zeros((len(k_data), len(y_features)))
    for idx, feature in enumerate(y_features):
        y_data[:, idx] = diff_cell_data[feature]['DATA']

    # Addressing the imbalanced data problem
    x_data, y_data = sub_sampling(x_data, y_data, np.logical_xor(x_data[:, 0], y_data[:, 0]), size = 1.5)
    x_data, y_data = sub_sampling(x_data, y_data, y_data.flatten() != 0, size = 1.0)

    # xor_indexs = np.logical_xor(x_data[:, 0], y_data[:, 0])
    #
    # print("ratio of changed now: %.2f" % (float(sum(xor_indexs)) / len(xor_indexs)))
    #
    # xor_indexs = y_data.flatten() != 0
    #
    # print("ratio of y post now: %.2f" % (float(sum(xor_indexs)) / len(xor_indexs)))

    if window_sz > 1:
        return sliding_window_of_data(x_data, y_data, window_sz)
    else:
        return [x_data, y_data]

def machine_learning_pipeline(x_features_list, y_features_list):

    print("include K\t diff cell type\tx features\twindow size\ttarget modification")
    for differentiated_cell_type in differentiated_cell_types:
        for idx, y_features in enumerate(y_features_list):
            for include_k in [1]:#0,
                for window_sz in WINDOW_SZS:
                    print("%d\t%s\t%s\t%d\t%s" %(include_k, differentiated_cell_type, ' '.join(x_features_list[idx]), window_sz, ' '.join(y_features)))
                    X, y = extract_target_data_by_features(differentiated_cell_type, x_features_list[idx], y_features, include_k, window_sz)
                    y = y.reshape(-1, )
                    [classifier_names, classifiers] = generate_classifiers(len(x_features_list[idx]), window_sz)
                    print("clf\tacc\trecall\tprec\tf1")
                    for cidx, cls_name in enumerate(classifier_names):
                        perform_machine_learning_prediction(cls_name, classifiers[cidx], X, y)
def generate_classifiers(feature_sz, window_sz):
    # data_sz = feature_sz * window_sz
    classifier_names = ["KNN"] # "Linear SVM", "RBF SVM", "Gaussian Process","Naive Bayes", "QDA","Decision Tree", "Neural Network",,, "AdaBoost", "Random Forest","LR",

    classifiers = [
        KNeighborsClassifier(n_neighbors=1),
    ]#SVC(kernel="linear", C=0.025), SVC(gamma=2, C=1), GaussianProcessClassifier(1.0 * RBF(1.0)),GaussianNB(),QuadraticDiscriminantAnalysis(), DecisionTreeClassifier(max_depth=5),
    # MLPClassifier(alpha=1e-5, hidden_layer_sizes= (50, 20), random_state = RANDOM_STATE, max_iter=10000, learning_rate='adaptive'),,AdaBoostClassifier()RandomForestClassifier(max_depth=5, n_estimators=10) LogisticRegression(C=10, penalty='l1', max_iter=10000),



    # classifier_names = ["LR"]
    # classifiers = [LogisticRegression(C=1, penalty='l1', max_iter=10000)]
    return [classifier_names, classifiers]
def perform_machine_learning_prediction(clf_name, clf, X, y):
    try:

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size= TEST_RATIO, random_state= RANDOM_STATE)
        clf = clf.fit(X, y)
        y_pred = clf.predict(X_test)
        print("%s\t%.2f\t%.2f\t%.2f\t%.2f" % (clf_name, accuracy_score(y_test, y_pred), recall_score(y_test, y_pred), precision_score(y_test, y_pred), f1_score(y_test, y_pred)))
        # print(classification_report(y_test, y_pred))
        print(confusion_matrix(y_test, y_pred))
    except ValueError as e:
        print(e)
if __name__ == "__main__":
    [x_features_list, y_features_list] = generate_the_combination_of_features()
    # machine_learning_pipeline(x_features_list, y_features_list)
    pass