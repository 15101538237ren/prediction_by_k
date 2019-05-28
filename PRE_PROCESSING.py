# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from itertools import combinations, islice
import os, re, pickle

TARGET_COLUMN_INDEX_OF_BED_FILE = 3

modification_names = ["H3K4me1", "H3K4me3"]#, "H3K9Me3", "H3K27ac", "H3K27Me3", "H3K36me3", "Methylation"
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
            if line_counter > 1000:
                break
            line_counter += 1
            if line_counter >= start_line:
                line_contents = re.split(sep, line.strip(line_end))
                try:
                    led = target_data_type(line_contents[target_col_index])
                    ret_value_list.append(led)
                except ValueError, e:
                    ret_value_list.append(0)
            line = input_file.readline()

    return ret_value_list

def combine_whole_dataset(region='GENOME', load=False):
    dataset_dir = os.path.join('DATA_SET', region)
    if load:
        with open(os.path.join(dataset_dir, region + '.pkl'), 'rb') as pkl_file:
            [k_data, histone_and_methy_data] = pickle.load(pkl_file)
    else:
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
        with open(os.path.join(dataset_dir, region + '.pkl'), 'wb') as pkl_file:
            pickle.dump([k_data, histone_and_methy_data], pkl_file, -1)
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
        x_features.append(modification_names)
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
def extract_target_data_by_features(original_dataset, differentiated_cell_type, x_features, y_features, include_k, window_sz):
    k_data, histone_and_methy_data = original_dataset

    hESC_data = histone_and_methy_data['hESC']
    diff_cell_data = histone_and_methy_data[differentiated_cell_type]

    x_data = np.zeros((len(k_data), len(x_features) + include_k))
    for idx, feature in enumerate(x_features):
        x_data[:, idx] = hESC_data[feature]['DATA']
    if include_k:
        x_data[:, -1] = k_data

    y_data = np.zeros((len(k_data), len(y_features)))
    for idx, feature in enumerate(y_features):
        y_data[:, idx] = diff_cell_data[feature]['DATA']
    if window_sz > 1:
        return sliding_window_of_data(x_data, y_data, window_sz)
    else:
        return [x_data, y_data]
def machine_learning_pipeline(original_dataset, x_features_list, y_features_list):
    window_sz = 5 # must odd number
    for include_k in [0, 1]:
        for differentiated_cell_type in differentiated_cell_types:
            for idx, y_features in enumerate(y_features_list):
                x_data, y_data = extract_target_data_by_features(original_dataset, differentiated_cell_type, x_features_list[idx], y_features, include_k, window_sz)
                
if __name__ == "__main__":
    original_dataset = combine_whole_dataset(region='GENOME', load=True)
    [x_features_list, y_features_list] = generate_the_combination_of_features()
    machine_learning_pipeline(original_dataset, x_features_list, y_features_list)