#!/usr/bin/env python
# coding: utf-8

# -*- coding: utf-8 -*-
import pandas as pd
from pandas.errors import ParserError
import numpy as np
import os, collections
import matplotlib.pyplot as plt

DATA_SET_DIR = '../DATA/Filtered_Feature_Analysis/'

GO_CLASS_NAMES = ["GO Molecular Function", "GO Biological Process", "GO Cellular Component"]
COL_NAMES = {'CATEGORY': '# Ontology', 'GO_ID': ' Term ID ', 'Q_Value': '  Hyper FDR Q-Val  '}

def underline_a_string(string_value):
    """
    Convert string value into its acronym, e.g What is Up
    :param string_value: input string
    :return: the acronym of string value
    """
    return "_".join(e for e in string_value.split(" "))

def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
def extract_GO_id_and_qvalue(class_name, input_fp, output_fp):
    """
    Extract the GO id in class name(e.g Biological Process) and the corresponding Hyper FDR Q-Value
    from input file path (input_fp) and output the result to output file path (output_fp)
    :param class_name: class name(e.g Biological Process)
    :param input_fp: input file path
    :param output_fp: output file path
    :return: None
    """
    print(input_fp)
    try:
        df = pd.read_csv(input_fp, sep='\t', header=1, index_col=False)
        extracted_df = df.loc[df[COL_NAMES['CATEGORY']] == class_name]
        if len(extracted_df):
            target_df = extracted_df[[COL_NAMES['GO_ID'], COL_NAMES['Q_Value']]]
            target_df.to_csv(output_fp, sep= "\t", header= False, index=False, line_terminator='\n')
    except ParserError:
        pass
def go_term_extraction(CLASS_NAME, INPUT_DIR, OUTPUT_DIR):
    OUT_SUBDIR = os.path.join(OUTPUT_DIR, underline_a_string(CLASS_NAME))
    mkdir(OUT_SUBDIR)
    for file_name in os.listdir(INPUT_DIR):
        if file_name.endswith(".tsv"):
            input_fp = os.path.join(INPUT_DIR, file_name)
            output_fp = os.path.join(OUT_SUBDIR, file_name)
            extract_GO_id_and_qvalue(CLASS_NAME, input_fp, output_fp)

def plot_distance_effect_on_go_term(INPUT_DIR, fig_dir):
    METHY_CLASS_NAMES = ["LM_K_0_5", "IM_K_0_5", "HM_K_0_5"]
    N_COL = 3
    N_ROW = 1
    EACH_SUB_FIG_SIZE = 5
    BASE_DISTANCE = 5
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
    for m_id, METHY_CLASS_NAME in enumerate(METHY_CLASS_NAMES):
        print("processing %s" % METHY_CLASS_NAME)
        ax = axs[m_id]
        dist_and_number_of_go_term_dict = {}
        input_dir = os.path.join(INPUT_DIR, METHY_CLASS_NAME)
        for file_name in os.listdir(input_dir):
            if file_name.endswith(".tsv"):
                input_fp = os.path.join(input_dir, file_name)
                dist = file_name.strip(".tsv").split("_")[-1]
                try:
                    df = pd.read_csv(input_fp, sep='\t', header=1, index_col=False)
                    df_categories = list(df[COL_NAMES['CATEGORY']].values)
                    num_of_go_term = len([item for item in df_categories if item in GO_CLASS_NAMES])
                except ParserError:
                    num_of_go_term = 0
                dist_and_number_of_go_term_dict[BASE_DISTANCE + int(dist)] = num_of_go_term
        sorted_dict = collections.OrderedDict(sorted(dist_and_number_of_go_term_dict.items()))
        dist_and_n = np.zeros((len(sorted_dict), 2))
        for kid, (key, value) in enumerate(sorted_dict.items()):
            dist_and_n[kid, 0] = key
            dist_and_n[kid, 1] = value
        ax.scatter(dist_and_n[:, 0], dist_and_n[:, 1], c='r', s= 5)
        ax.plot(dist_and_n[:, 0], dist_and_n[:, 1], c='b')
        max_n_x = dist_and_n[np.argmax(dist_and_n[:, 1]), 0]
        max_n = np.max(dist_and_n[:, 1])
        ax.vlines(x=max_n_x, ymin=0, ymax=max_n, color='r')
        ax.text(max_n_x + .5, 0.7 * (max_n), "Peak at %d kb to TSS" % max_n_x, fontsize=12)
        ax.set_xlabel('Dist to TSS')
        ax.set_ylabel('# of GO terms')
        ax.set_xlim([0, 220])
        ax.set_ylim([0, max_n])
        ax.set_title(METHY_CLASS_NAME)
    plt.savefig(os.path.join(fig_dir, "distance_effect.png"), dpi=200)
    # plt.show()

def plot_k_distance_to_tss(dist_dir, fig_dir):
    EACH_SUB_FIG_SIZE = 4
    methy_categories = ["LM", "IM", "HM"]
    N_ROW = len(methy_categories)

    k_categories = ["K_0_5", "K_5_15", "K_15_25", "K_25_35", "K_35_45", "K_45_55", "K_55_65",
                    "K_65_75", "K_75_85", "K_85_95", "K_95_100", ""]
    N_COL = len(k_categories)
    N_BIN = 100
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
    for i, methy_class in enumerate(methy_categories):
        for j, k_class in enumerate(k_categories):
            ax = axs[i][j]
            if k_class != "":
                input_fp = os.path.join(dist_dir, "%s_%s.tsv" %( methy_class, k_class))
            else:
                input_fp = os.path.join(dist_dir, "%s.tsv" % (methy_class))
            df = pd.read_csv(input_fp, sep="\t", header=None).values
            dists = df[:, -1]
            k_within_5kb_to_tss = dists[dists < 5000]
            proportion_of_k_out_5kb_to_tss = (1. - float(len(k_within_5kb_to_tss))/ len(dists)) * 100.0
            ax.hist(dists, N_BIN, density=True, facecolor='b', alpha=0.75)
            ax.text(50000, 0.00003, "%.0f%% k > 5kb" % proportion_of_k_out_5kb_to_tss)
            ax.set_xlim([-250000, 250000])
            ax.set_ylim([0, 0.000045])
            ax.set_xticks([])
            if i == 0:
                if k_class != "":
                    ax.set_title("%s" % k_class)
                else:
                    ax.set_title("%s" % "All")
            elif i == len(methy_categories) - 1:
                ax.set_xlabel("Dist to TSS")
                ax.set_xticks([-150000, 0, 150000])
            if j == 0:
                ax.set_ylabel('%s' % methy_class)
            else:
                ax.set_yticks([])
    plt.savefig(os.path.join(fig_dir, "K_distance_to_TSS.png"), dpi=200)
def annotation_stack_barplot(dist_dir, fig_dir, distance_category = 0): # distance category: 0: no restriction, 1: k <= 5kb, 2: k > 5kb
    EACH_SUB_FIG_SIZE = 4
    ChrmoHMM_LABELS = "Active Promoter", "Weak Promoter", "Poised Promoter", "Strong Enhancer", "Strong Enhancer", "Weak Enhancer", "Weak Enhancer", "Insulator", "Txn Transition", "Txn Elongation", "Weak Txn", "Repressed", "Heterochrom/lo", "Repetitive/CNV", "Repetitive/CNV"
    K_CATEGORY_LABEL = ["all_dist", "less_than_5kb", "greater_than_5kb"]
    NUMBER_OF_CHRMM_CLASS = len(ChrmoHMM_LABELS)
    chrHMMs = np.arange(NUMBER_OF_CHRMM_CLASS)

    methy_categories = ["LM", "IM", "HM"]
    N_ROW = len(methy_categories)

    k_categories = ["K_0_5", "K_5_15", "K_15_25", "K_25_35", "K_35_45", "K_45_55", "K_55_65",
                    "K_65_75", "K_75_85", "K_85_95", "K_95_100", ""]
    mod_k_categories = [k_categories[i] if i != len(k_categories) - 1 else "All" for i in range(len(k_categories)) ]
    ind = np.arange(len(k_categories))
    cm = plt.get_cmap('gist_rainbow')
    fig, axs = plt.subplots(N_ROW, 1, figsize=(EACH_SUB_FIG_SIZE * 5, N_ROW * EACH_SUB_FIG_SIZE))
    for i, methy_class in enumerate(methy_categories):
        ax = axs[i]
        class_data = np.zeros((NUMBER_OF_CHRMM_CLASS, len(k_categories)))
        for j, k_class in enumerate(k_categories):
            if k_class != "":
                input_fp = os.path.join(dist_dir, "%s_%s.bed" %( methy_class, k_class))
            else:
                input_fp = os.path.join(dist_dir, "%s.bed" % (methy_class))
            df = pd.read_csv(input_fp, sep="\t", header=None).values
            if distance_category == 1:
                dists = np.abs(df[:, -2]) <= 5000
                classes = df[dists, -1]
            elif distance_category == 2:
                dists =np.abs(df[:, -2]) > 5000
                classes = df[dists, -1]
            else:
                classes = df[:, -1]
            unique, counts = np.unique(classes, return_counts=True)
            normed_counts = counts / float(counts.sum())
            unique_dc = dict(zip(unique, normed_counts))
            for cls in chrHMMs:
                if cls + 1 in unique_dc.keys():
                    class_data[cls][j] = unique_dc[cls + 1]
                else:
                    class_data[cls][j] = 0.
        pls = []
        for cls in chrHMMs:
            if cls == 0:
                pl = ax.bar(ind, tuple(list(class_data[cls, :])), 0.35)
            else:
                sum_prev = list(np.sum(class_data[0: cls, :], axis=0))
                pl = ax.bar(ind, tuple(list(class_data[cls, :])), 0.35, bottom=tuple(sum_prev))
            for item in pl:
                item.set_color(cm(1. * cls / NUMBER_OF_CHRMM_CLASS))
            pls.append(pl)
        ax.set_ylabel('Percentage in %s' % methy_class)
        ax.set_xticks(ind)
        ax.set_xticklabels(mod_k_categories)
        ax.set_yticks(np.arange(0.2, 1.1, 0.2))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        if i == 0:
            ax.legend(tuple(pls), tuple(ChrmoHMM_LABELS), loc='center left', bbox_to_anchor=(1, 0.5))
        elif i == len(methy_class) - 1:
            ax.set_xlabel('ChromHMM')
    plt.savefig(os.path.join(fig_dir, "ChromoHMM_%s.png" % K_CATEGORY_LABEL[distance_category]), dpi=200)
if __name__ == "__main__":

    fig_dir = "../FIGURES/Functional_Analysis"
    mkdir(fig_dir)
    INPUT_DIR = os.path.join(DATA_SET_DIR, "DATA", "DISTANCE_EXPLORATION_OF_GO")
    # plot_distance_effect_on_go_term(INPUT_DIR, fig_dir)
    dist_dir = os.path.join(DATA_SET_DIR, "DATA", "K_DISTANCE_TO_TSS")
    annotation_stack_barplot(dist_dir, fig_dir, distance_category=1)
    annotation_stack_barplot(dist_dir, fig_dir, distance_category=2)