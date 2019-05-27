# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from itertools import combinations
import os, re

TARGET_COLUMN_INDEX_OF_BED_FILE = 3
modification_names = ["H3K4me1", "H3K4me3", "H3K9Me3", "H3K27ac", "H3K27Me3", "H3K36me3", "Methylation"]
differentiated_cell_types = ["dEC", "dME", "dEN"]
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
                led = target_data_type(line_contents[target_col_index])
                ret_value_list.append(led)
            line = input_file.readline()
    return ret_value_list

def combine_whole_dataset(region='GENOME'):
    dataset_dir = os.path.join('DATA_SET', region)

    histone_and_methy_data = {}

    k_data = read_tab_seperated_file_and_get_target_column(os.path.join(dataset_dir, 'K_mean.bed'),
                                                               TARGET_COLUMN_INDEX_OF_BED_FILE)

    df = pd.read_csv('DATA/DATA_LABEL.txt', sep='\s+', index_col=0)

    for idx in df.index:
        idx_s = str(idx)
        histone_and_methy_data[idx_s] = {}
        for col in df.columns:
            histone_and_methy_data[idx_s][col] = {}
            histone_and_methy_data[idx_s][col]['NAME'] = df.loc[idx, col]
            print df.loc[idx, col]
            histone_and_methy_data[idx_s][col]['DATA'] = read_tab_seperated_file_and_get_target_column(os.path.join(dataset_dir, df.loc[idx, col] + '.bed'), TARGET_COLUMN_INDEX_OF_BED_FILE)
    return [k_data, histone_and_methy_data]

def generate_the_combination_of_features():
    x_features = [] # the input feature for machine learning model
    y_features = []
    for target_modification_name in modification_names:
        left_features = [feature for feature in modification_names if feature != target_modification_name]
        for comb_feature_sz in range(len(left_features)):
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

    return [x_features, y_features]

def extract_target_data_by_features(original_dataset, differentiated_cell_type, x_feature, y_feature, include_k):
    k_data, histone_and_methy_data = original_dataset

    hESC_data = histone_and_methy_data['hESC']
    diff_cell_data = histone_and_methy_data[differentiated_cell_type]

    x_data = np.zeros((len(k_data), len(x_feature) + include_k))
    for idx, feature in enumerate(x_feature):
        x_data[:, idx] = hESC_data[feature]['DATA']
    if include_k:
        x_data[:, -1] = k_data

    y_data = np.zeros((len(k_data), len(y_feature)))
    for idx, feature in enumerate(y_feature):
        y_data[:, idx] = diff_cell_data[feature]['DATA']

    return [x_data, y_data]
def machine_learning_pipeline(original_dataset, x_feature, y_feature):
    for include_k in [0, 1]:
        for differentiated_cell_type in differentiated_cell_types:
            x_data, y_data = extract_target_data_by_features(original_dataset, differentiated_cell_type, x_feature, y_feature, include_k)
if __name__ == "__main__":
    original_dataset = combine_whole_dataset()
    [x_feature, y_feature] = generate_the_combination_of_features()
    machine_learning_pipeline(original_dataset, x_feature, y_feature)