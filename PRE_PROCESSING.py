# -*- coding: utf-8 -*-
import pandas as pd
import os

TARGET_COLUMN_INDEX_OF_BED_FILE = 3

def read_tab_seperated_file_and_get_target_column(input_file_path, target_col_index, start_line= 1, sep="\t",line_end = "\n"):
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
                line_contents = line.split(sep)
                led = line_contents[target_col_index].strip(line_end)
                ret_value_list.append(led)
            line = input_file.readline()
    return ret_value_list

def combine_whole_dataset(region='GENOME'):
    dataset_dir = os.path.join('DATA_SET', region)

    histone_and_methy_data = {}

    k_data = read_tab_seperated_file_and_get_target_column(os.path.join(dataset_dir, 'K_mean.bed'),
                                                               TARGET_COLUMN_INDEX_OF_BED_FILE)

    df = pd.read_csv('DATA/DATA_LABEL.txt', sep='\t', index_col=0)

    for idx in df.index:
        idx_s = str(idx)
        histone_and_methy_data[idx_s] = {}
        for col in df.columns:
            histone_and_methy_data[idx_s][col] = {}
            histone_and_methy_data[idx_s][col]['NAME'] = df.loc[idx, col]
            histone_and_methy_data[idx_s][col]['DATA'] = read_tab_seperated_file_and_get_target_column(os.path.join(dataset_dir, df.loc[idx, col] + '.bed'), TARGET_COLUMN_INDEX_OF_BED_FILE)

if __name__ == "__main__":
    combine_whole_dataset()