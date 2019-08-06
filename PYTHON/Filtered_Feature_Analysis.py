#!/usr/bin/env python
# coding: utf-8

# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os, re, math, time
import matplotlib.pyplot as plt
from sklearn.metrics import mutual_info_score
from scipy.spatial.distance import jensenshannon
from sklearn import manifold
from matplotlib.ticker import NullFormatter
from scipy.stats import variation

# %matplotlib inline


def mkdirs(dir_paths):
    for dir_path in dir_paths:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)


DATA_SET_DIR = '/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis/'
K_Cluster_dir = os.path.join(DATA_SET_DIR, "K_Cluster_Analysis")
ChrmoHMM_LABELS = "Active Promoter", "Weak Promoter", "Poised Promoter", "Strong Enhancer", "Strong Enhancer", "Weak Enhancer", "Weak Enhancer", "Insulator", "Txn Transition", "Txn Elongation", "Weak Txn", "Repressed", "Heterochrom/lo", "Repetitive/CNV", "Repetitive/CNV"
NUMBER_OF_CHRMM_CLASS = len(ChrmoHMM_LABELS)

EACH_SUB_FIG_SIZE = 5
CHROMOSOMES = [i for i in range(1, 23)]

LOG_SCALE = False
JSD = True  # Using JSD or Mutual information
HISTOGRAM = False  # plot histogram or the distribution by pdf
MI_BIN = 10
N_CLASSES = 5
NBIN = 20 if LOG_SCALE else 20
METRIC = "JSD" if JSD else "1-MI"


def get_percentile_index_arr(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX):
    df = pd.read_csv(DATA_FP, sep='\t', header=None).values
    ks = df[:, COL_INDEX]
    if K_analyze:
        if LOG_SCALE:
            ks[ks <= 0] = 0.001
            ks = np.log10(ks)
            df[:, COL_INDEX] = ks
            percentiles = [-2.0]
        else:
            percentiles = [0.0]
    else:
        if METHY_K_LESS:
            percentiles = [0.0]
        else:
            percentiles = [0.8]

    percentiles_to_query = np.array([i * 20 for i in range(1, N_CLASSES + 1)]) if K_analyze else np.array(
        [25, 50, 65, 100])
    percentiles_results = list(np.percentile(ks, percentiles_to_query))
    for item in percentiles_results:
        percentiles.append(float(item))
    print(percentiles)
    df_indexs_arr = []
    for p_idx in range(len(percentiles) - 1):
        df_indexs = np.logical_and(ks > percentiles[p_idx], ks < percentiles[p_idx + 1])
        df_indexs_arr.append(df_indexs)
    return df_indexs_arr, df, percentiles


def hist_K_with_5_percentiles(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX, FIG_DIR):
    fig_dir = os.path.join(FIG_DIR, "Others")
    mkdirs([fig_dir])
    df_indexs_arr, df, percentiles = get_percentile_index_arr(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    Ks = df[:, COL_INDEX]
    xhist, _, _ = ax.hist(Ks, 20, density=True, facecolor='b', alpha=0.75)
    pK = xhist / xhist.sum()
    entropy = -np.sum(pK * np.log10(pK))
    for percentile in percentiles:
        ax.axvline(x=percentile, color='r')
    if K_analyze:
        if LOG_SCALE:
            xlabel = 'log10(k)'
            xlim = [-1, 1]
        else:
            xlabel = 'k'
            xlim = [0, 10]
    else:
        xlabel = 'methylation'
        if METHY_K_LESS:
            xlim = [0, 0.8]
        else:
            xlim = [0.8, 1]
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Probability')
    ax.set_title("Histogram of %s, Entropy: %.2f" % (xlabel, entropy))
    ax.set_xlim(xlim)
    plt.savefig(os.path.join(fig_dir, "Hist.png"), dpi=200)


def PREPARE_FIG_DIR(FIGURE_DIR, K_analyze):
    if K_analyze:
        if LOG_SCALE:
            sub_fig_dir = "K_Log10"
        else:
            sub_fig_dir = "K"
    else:
        sub_fig_dir = "Methy"
    print(sub_fig_dir + "_" + METRIC)
    FIG_DIR = os.path.join(FIGURE_DIR, sub_fig_dir)
    mkdirs([FIG_DIR])

    return FIG_DIR

# METHY_K_LESS = True
# FILE_NAME = "METHY_K_LESS" if METHY_K_LESS else "METHY_K_GREATER"
# DATA_FP = os.path.join(DATA_SET_DIR, "DATA", FILE_NAME + ".tsv")
# FIGURE_DIR = os.path.join("../FIGURES/Filtered_Feature_Analysis/", FILE_NAME)
# mkdirs([FIGURE_DIR])
# for K_analyze in [False, True]:
#     FIG_DIR = PREPARE_FIG_DIR(FIGURE_DIR, K_analyze)
#     COL_INDEX = 4 if K_analyze else 3
#     hist_K_with_5_percentiles(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX, FIG_DIR)

def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi


def plot_ChrMM(DATA_FP, K_analyze, COL_INDEX, METHY_K_LESS, METHY_LESS_NAME, FIG_DIR, CLASS_LABELS):
    dist_dict = {}
    entropy_dict = {}
    fig_dir = os.path.join(FIG_DIR, "ChromHMM")
    mkdirs([fig_dir])
    cm = plt.get_cmap('gist_rainbow')
    ind = np.arange(len(CLASS_LABELS))
    chrHMMs = np.arange(NUMBER_OF_CHRMM_CLASS)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    df_indexs_arr, df, percentiles = get_percentile_index_arr(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX)

    class_data = np.zeros((len(chrHMMs), len(ind)))
    class_data_non_normed = np.zeros((len(chrHMMs), len(ind)))
    class_fp = os.path.join(DATA_SET_DIR, "ChrmoHMM", METHY_LESS_NAME, "HmmH1hescHMM.csv")
    categories = pd.read_csv(class_fp, sep="\t", header=None).values
    for cid, class_label in enumerate(CLASS_LABELS):
        classes = categories[df_indexs_arr[cid], -1].astype(int)
        unique, counts = np.unique(classes, return_counts=True)
        unique = unique.astype(int)
        normed_counts = counts / float(counts.sum())
        unique_dc = dict(zip(unique, normed_counts))
        unique_dc_non_normed = dict(zip(unique, counts))
        for cls in chrHMMs:
            class_data[cls][cid] = unique_dc[cls + 1]
            class_data_non_normed[cls][cid] = unique_dc_non_normed[cls + 1]
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
    ax.set_xlabel('ChromHMM Class')
    ax.set_ylabel('Percentage')
    ax.set_xticks(ind)
    ax.set_xticklabels(CLASS_LABELS)
    ax.set_yticks(np.arange(0.2, 1.1, 0.2))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(tuple(pls), tuple(ChrmoHMM_LABELS), loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(os.path.join(fig_dir, "ChromHMM_PERCENTAGE_OF_K.png"), dpi=200)

    N_ROW = 3
    N_COL = 5
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))

    for cls in chrHMMs:
        row = int(cls / N_COL)
        col = int(cls % N_COL)
        ax = axs[row][col]
        y_max = np.max(class_data[cls, :]) + 0.08
        pl = ax.bar(ind, tuple(list(class_data[cls, :])), 0.35)
        for item in pl:
            item.set_color(cm(1. * cls / NUMBER_OF_CHRMM_CLASS))
        for cid, class_label in enumerate(CLASS_LABELS):
            ax.text(ind[cid] - 0.2, class_data[cls, cid] + 0.002, str(round(class_data[cls, cid] * 100.0, 2)) + "%",
                    fontsize=14, weight="bold")
        if col == 0:
            ax.set_ylabel('Percentage')
        ax.set_xticks(ind)
        ax.set_xticklabels(CLASS_LABELS, fontsize=10)
        ax.set_ylim(0, y_max)
        ax.set_title(ChrmoHMM_LABELS[cls], fontsize=18)
    plt.savefig(os.path.join(fig_dir, "ChromHMM_SEP_PERCENTAGE_OF_K.png"), dpi=200)

    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
    for cls in chrHMMs:
        row = int(cls / N_COL)
        col = int(cls % N_COL)
        ax = axs[row][col]
        class_data_non_normed[cls, :] = class_data_non_normed[cls, :] / np.sum(class_data_non_normed[cls, :])
        pl = ax.bar(ind, tuple(list(class_data_non_normed[cls, :])), 0.35)
        for item in pl:
            item.set_color(cm(1. * cls / NUMBER_OF_CHRMM_CLASS))
        for cid, class_label in enumerate(CLASS_LABELS):
            ax.text(ind[cid] - 0.2, class_data_non_normed[cls, cid] + 0.002,
                    str(round(class_data_non_normed[cls, cid] * 100.0, 2)) + "%", fontsize=14, weight="bold")
        if col == 0:
            ax.set_ylabel('Percentage')
        ax.set_xticks(ind)
        ax.set_xticklabels(CLASS_LABELS, fontsize=10)
        y_max = np.max(class_data_non_normed[cls, :]) + 0.08
        ax.set_ylim(0, y_max)
        ax.set_title(ChrmoHMM_LABELS[cls], fontsize=18)
    plt.savefig(os.path.join(fig_dir, "ChromHMM_SEP_PERCENTAGE_OF_ChroMM.png"), dpi=200)

    # Plot the pdf
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
    if K_analyze:
        if LOG_SCALE:
            xlim = [-1, 1]
        else:
            xlim = [0, 10]
    else:
        if METHY_K_LESS:
            xlim = [0, 0.8]
        else:
            xlim = [0.8, 1]
    for cls in chrHMMs:
        row = int(cls / N_COL)
        col = int(cls % N_COL)
        ax = axs[row][col]
        ks_indexs = categories[:, -1].astype(int) == (cls + 1)
        ks = df[ks_indexs, COL_INDEX]
        xhist, _, _ = ax.hist(df[:, COL_INDEX], bins=NBIN, density=True, color='black', edgecolor='black', alpha=0.5,
                              linewidth=0.5)  #
        yhist, _, _ = ax.hist(ks, bins=NBIN, density=True, color=cm(1. * cls / NUMBER_OF_CHRMM_CLASS),
                              edgecolor='black', alpha=0.5, linewidth=0.5)  #
        if JSD:
            distance = jensenshannon(xhist, yhist)
        else:
            max_mi = calc_MI(xhist, xhist, MI_BIN)
            mi = calc_MI(xhist, yhist, MI_BIN) / max_mi
            distance = 1 - mi
        dist_dict[ChrmoHMM_LABELS[cls]] = distance
        pK = yhist / yhist.sum()
        entropy = -np.sum(pK * np.log10(pK))
        entropy_dict[ChrmoHMM_LABELS[cls]] = entropy
        print(ChrmoHMM_LABELS[cls] + "\t" + str(round(float(distance), 2)))
        ax.text(float(np.mean(np.array(xlim))), max([float(np.max(xhist)), float(np.max(yhist))]) * 0.75,
                "%s:%.2f" % (METRIC, distance), color='black', fontsize=18)

        ax.text(float(np.mean(np.array(xlim))), max([float(np.max(xhist)), float(np.max(yhist))]) * 0.6,
                "%s:%.2f" % ("Entropy", entropy), color='black', fontsize=18)
        ax.set_ylabel('Probability')
        ax.set_xlim(xlim)
        ax.set_title(ChrmoHMM_LABELS[cls], fontsize=18)
    plt.savefig(os.path.join(fig_dir, "ChromHMM_DIST_ChroMM_" + METRIC + ".png"), dpi=200)
    return [dist_dict, entropy_dict]


def plot_K_hist_in_different_markers(DATA_FP, BASE_DIR, SUB_DIR_NAME, K_analyze, COL_INDEX, METHY_K_LESS, FIG_DIR, CLASS_LABELS):
    dist_dict = {}
    entropy_dict = {}
    df_indexs_arr, df, _ = get_percentile_index_arr(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX)
    ind = np.arange(len(CLASS_LABELS))
    cm = plt.get_cmap('gist_rainbow')
    input_fps = []
    file_names = []
    for file_name in os.listdir(BASE_DIR):
        if file_name.endswith(".csv"):
            input_fps.append(os.path.join(BASE_DIR, file_name))
            file_names.append(file_name.strip(".csv"))
    N_FILES = len(input_fps)
    N_COL = 5
    Max_NROW = 3
    N_SUBFIG_PER_FIG = Max_NROW * N_COL
    NFIG = int(math.ceil(N_FILES / N_SUBFIG_PER_FIG))
    fig_dir = os.path.join(FIG_DIR, SUB_DIR_NAME)
    mkdirs([fig_dir])
    for i in range(NFIG):
        if NFIG > 1:
            fig_fp = os.path.join(fig_dir, SUB_DIR_NAME + "_" + str(i) + ".png")
        else:
            fig_fp = os.path.join(fig_dir, SUB_DIR_NAME + ".png")
        base_index = i * N_SUBFIG_PER_FIG
        N_remaining_files = N_FILES - base_index
        N_ROW = int(math.ceil((N_remaining_files) / N_COL)) if N_remaining_files <= N_SUBFIG_PER_FIG else Max_NROW

        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
        SUB_FIG_RANGE = N_SUBFIG_PER_FIG if N_remaining_files > N_SUBFIG_PER_FIG else N_remaining_files

        for j in range(SUB_FIG_RANGE):
            row = j // N_COL
            col = j % N_COL
            if N_ROW == 1:
                ax = axs[col]
            else:
                ax = axs[row][col]
            index = base_index + j
            class_fp = input_fps[index]
            categories = pd.read_csv(class_fp, sep="\t", header=None).values
            class_data = []
            cat_indexs = categories[:, -1].astype(int) == 1
            for cid, class_label in enumerate(CLASS_LABELS):
                df_indexs = df_indexs_arr[cid]
                intersect = np.logical_and(df_indexs, cat_indexs)
                # int_data = df[intersect]
                cls_count = np.sum(categories[intersect, -1])
                class_data.append(cls_count)
            class_data = np.array(class_data)
            class_data = class_data[:] / np.sum(class_data)
            y_max = np.max(class_data) + 0.08
            pl = ax.bar(ind, tuple(list(class_data)), 0.35)
            for item in pl:
                item.set_color(cm(1. * j / SUB_FIG_RANGE))
            for cid, class_label in enumerate(CLASS_LABELS):
                ax.text(ind[cid] - 0.2, class_data[cid] + 0.002,
                        str(round(class_data[cid] * 100.0, 2)) + "%", fontsize=10, weight="bold")
            if col == 0:
                ax.set_ylabel("Percentage", fontsize=18)
            ax.set_xticks(ind)
            ax.set_xticklabels(CLASS_LABELS, fontsize=12)
            ax.set_xlim(-0.5, len(CLASS_LABELS) - 0.5)
            ax.set_ylim(0, y_max)
            ax.set_title(file_names[index], fontsize=18)
        plt.savefig(fig_fp, dpi=200)
        if NFIG > 1:
            fig_fp = os.path.join(fig_dir, SUB_DIR_NAME + "_DIST_" + str(i) + "_" + METRIC + ".png")
        else:
            fig_fp = os.path.join(fig_dir, SUB_DIR_NAME + "_DIST_" + METRIC + ".png")
        # Plot the pdf
        if K_analyze:
            if LOG_SCALE:
                xlim = [-1, 1]
            else:
                xlim = [0, 10]
        else:
            if METHY_K_LESS:
                xlim = [0, 0.8]
            else:
                xlim = [0.8, 1]
        fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
        for j in range(SUB_FIG_RANGE):
            row = j // N_COL
            col = j % N_COL
            if N_ROW == 1:
                ax = axs[col]
            else:
                ax = axs[row][col]
            index = base_index + j
            class_fp = input_fps[index]
            categories = pd.read_csv(class_fp, sep="\t", header=None).values
            class_data = []
            ks_indexs = categories[:, -1].astype(int) == 1
            ks = df[ks_indexs, COL_INDEX]
            xhist, _, _ = ax.hist(df[:, COL_INDEX], bins=NBIN, density=True, color='black', edgecolor='black',
                                  alpha=0.5,
                                  linewidth=0.5)  #
            yhist, _, _ = ax.hist(ks, bins=NBIN, density=True, color=cm(1. * j / SUB_FIG_RANGE), edgecolor='black',
                                  alpha=0.5, linewidth=0.5)  #
            if JSD:
                distance = jensenshannon(xhist, yhist)
            else:
                max_mi = calc_MI(xhist, xhist, MI_BIN)
                mi = calc_MI(xhist, yhist, MI_BIN) / max_mi
                distance = 1 - mi
            print(file_names[index] + "\t" + str(round(float(distance), 2)))
            dist_dict[file_names[index]] = distance
            pK = yhist / yhist.sum()
            entropy = -np.sum(pK * np.log10(pK))
            entropy_dict[file_names[index]] = entropy
            ax.text(float(np.mean(np.array(xlim))), max([float(np.max(xhist)), float(np.max(yhist))]) * 0.75,
                    "%s:%.2f" % (METRIC, distance), color='black', fontsize=18)

            ax.text(float(np.mean(np.array(xlim))), max([float(np.max(xhist)), float(np.max(yhist))]) * 0.6,
                    "%s:%.2f" % ("Entropy", entropy), color='black', fontsize=18)
            ax.set_ylabel('Probability')
            ax.set_xlim(xlim)
            ax.set_title(file_names[index], fontsize=18)
        plt.savefig(fig_fp, dpi=200)
    return [dist_dict, entropy_dict]


def plot_k_histogram_versus_diff_methylation_interval():
    data_fp = os.path.join(DATA_SET_DIR, "DATA", "METHY_AND_K.tsv")
    FIGURE_DIR = "../FIGURES/Filtered_Feature_Analysis/"
    df = pd.read_csv(data_fp, sep='\t', header=None).values
    METHY_COL_IDNEX = 3
    K_COL_INDEX = 4
    BINS = 20
    N_COL = 5
    N_ROW = 5
    fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
    METHY_LOWER_BOUNDS = [i*0.1 for i in range(0, 10)]
    TITLE_LABELS = ["%.1f-%.1f" %(item, item + .1) for item in METHY_LOWER_BOUNDS]

    for m_idx, methy_lower_bound in enumerate(METHY_LOWER_BOUNDS):
        row = int(m_idx / N_COL)
        col = int(m_idx % N_COL)
        ax = axs[row][col]
        Methys = df[:, METHY_COL_IDNEX]
        target_methy_rows = np.logical_and(Methys > methy_lower_bound, Methys < (methy_lower_bound + 0.1))
        Ks = df[target_methy_rows , K_COL_INDEX]
        xhist, _, _ = ax.hist(Ks, BINS, density=True, facecolor='b', alpha=0.75)
        pK = xhist / xhist.sum()
        entropy = -np.sum(pK * np.log10(pK))
        ax.set_xlabel('K')
        ax.set_ylabel('Probability')
        ax.set_title("%s in %s, Entropy: %.2f" % ('K', TITLE_LABELS[m_idx], entropy))
        ax.set_xlim([0, 10])
    plt.savefig(os.path.join(FIGURE_DIR, "Hist_OF_K_IN_DIFF_METHY_INTERVAL.png"), dpi=200)
def run_plot_pipeline():
    DATA_RECORDS = []
    for K_analyze in [True, False]:
        for METHY_K_LESS in [False, True]:  #
            TMP_ARR = []
            FILE_NAME = "METHY_K_LESS" if METHY_K_LESS else "METHY_K_GREATER"
            DATA_FP = os.path.join(DATA_SET_DIR, "DATA", FILE_NAME + ".tsv")
            FIGURE_DIR = os.path.join("../FIGURES/Filtered_Feature_Analysis/", FILE_NAME)
            mkdirs([FIGURE_DIR])
            CLASS_LABELS = ['slowest', 'slow', 'middle', 'fast', 'fastest'] if K_analyze else ['<25%', '25%-50%',
                                                                                               '50%-65%',
                                                                                               '>65%']
            N_CLASSES = len(CLASS_LABELS)

            FIG_DIR = PREPARE_FIG_DIR(FIGURE_DIR, K_analyze)
            COL_INDEX = 4 if K_analyze else 3
            hist_K_with_5_percentiles(DATA_FP, METHY_K_LESS, K_analyze, COL_INDEX, FIG_DIR)
            element = plot_ChrMM(DATA_FP, K_analyze, COL_INDEX, METHY_K_LESS, FILE_NAME, FIG_DIR, CLASS_LABELS)
            TMP_ARR.append(element)
            regions = ["General", "Genomic_Regions", "Histone_Modification", "Histone_Modification_Enzyme",
                       "TFBS"]  # []  #,
            for REGION in regions:
                element = plot_K_hist_in_different_markers(DATA_FP, os.path.join(DATA_SET_DIR, REGION, FILE_NAME), REGION,
                                                           K_analyze,
                                                           COL_INDEX, METHY_K_LESS, FIG_DIR, CLASS_LABELS)
                TMP_ARR.append(element)
# The input should be a csv file with columns: chr_i, start_site, k
def cluster_k_by_distance(data_fp, out_fp, maximum_distance_winthin_cluster=300, minimum_number_of_k_per_cluster=5):
    df = pd.read_csv(data_fp, sep='\t', header=None).values
    chrs = df[:, 0].astype(int)
    cords = df[:, 1].astype(int)
    methys = df[:, 3]
    ks = df[:, 4]
    N_SITES = cords.shape[0] #number of the CpG sites
    cluster_ids = np.arange(0, N_SITES) #initiate the cluster ids
    ind_current_site = 0
    cluster_ks = {}
    cluster_methys = {}
    cluster_ds = {}
    cluster_regions = {}
    while ind_current_site < N_SITES - 1:
        curr_loc = cords[ind_current_site]
        ind_next_site = ind_current_site + 1
        for ind_nxt in range(ind_next_site, N_SITES):
            distance= cords[ind_nxt] - curr_loc
            if distance > maximum_distance_winthin_cluster or (ind_nxt == N_SITES - 1 and distance <= maximum_distance_winthin_cluster):
                cluster_ids[ind_current_site: ind_nxt] = ind_current_site
                cluster_regions[ind_current_site] = [chrs[ind_current_site], curr_loc, cords[ind_nxt - 1]]
                cluster_methys[ind_current_site] = methys[ind_current_site: ind_nxt]
                cluster_ks[ind_current_site] = ks[ind_current_site: ind_nxt]
                cluster_ds[ind_current_site] = cords[ind_current_site + 1: ind_nxt] - cords[ind_current_site: ind_nxt - 1]
                ind_current_site = ind_nxt # update the index of the nex site which has less than maximum_distance_winthin_cluster with current index site
                break
    unique, counts = np.unique(cluster_ids, return_counts=True)
    # unique_dc = dict(zip(unique, counts))
    clusters_ids_interested = unique[counts >= minimum_number_of_k_per_cluster]

    data = {}
    with open(out_fp, "w") as out_file:
        for cid in np.arange(0, clusters_ids_interested.shape[0]):
            cluster_id = clusters_ids_interested[cid]
            data[cid + 1] = {}
            k_in_this_cluster = cluster_ks[cluster_id]
            methys_in_this_cluster = cluster_methys[cluster_id]
            d_in_this_cluster = cluster_ds[cluster_id]

            data[cid + 1]['region'] = cluster_regions[cluster_id]
            data[cid + 1]['ks'] = k_in_this_cluster
            data[cid + 1]['num_k'] = k_in_this_cluster.shape[0]
            data[cid + 1]['median_k'] = np.median(k_in_this_cluster)
            data[cid + 1]['coef_var_k'] = variation(k_in_this_cluster)

            data[cid + 1]['methys'] = methys_in_this_cluster
            data[cid + 1]['median_methy'] = np.median(methys_in_this_cluster)
            data[cid + 1]['coef_var_methy'] = variation(methys_in_this_cluster)

            data[cid + 1]['ds'] = cluster_ds[cluster_id]
            data[cid + 1]['mean_d'] = np.mean(d_in_this_cluster)
            data[cid + 1]['coef_var_d'] = variation(d_in_this_cluster)
            region = data[cid + 1]['region']
            ltw = "chr%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % ( region[0], region[1], region[2], data[cid + 1]['num_k'],
                                                                                data[cid + 1]['median_k'], data[cid + 1]['coef_var_k'], data[cid + 1]['median_methy'],
                                                                                data[cid + 1]['coef_var_methy'], data[cid + 1]['mean_d'], data[cid + 1]['coef_var_d'])
            out_file.write(ltw)
        print("write %s successful!" % out_fp)

def prepare_cluster_data():
    K_Cluster_dir = os.path.join(DATA_SET_DIR, "K_Cluster_Analysis")
    out_dir = os.path.join(K_Cluster_dir, "CLUSTER_DATA")
    mkdirs([out_dir])
    for chr_i in range(2, 12):
        input_fp = os.path.join(K_Cluster_dir, 'INPUT_DATA','chr%d.bed' % chr_i)
        out_fp = os.path.join(out_dir, 'chr%d.bed' % chr_i)
        cluster_k_by_distance(input_fp, out_fp)

def label_dmr():
    DMR_DIR = os.path.join(DATA_SET_DIR, "DMR")
    input_fp = os.path.join(DMR_DIR, "DATA", "DMR_ORIGIN_TABLE.bed")
    out_fp = os.path.join(DMR_DIR, "DMR_Labeled.bed")
    df = pd.read_csv(input_fp, sep='\t', header=None).values
    out_data = np.zeros((df.shape[0], 4))
    out_data[:, 0:3] = df[:, 0:3]
    col_indexs = [np.where(r == 1)[0][0] for r in df[:, -3:].astype(int)]
    hesc_methy_levels = df[:, 3]
    other_data_levels = df[:, 4: 4+3]
    labels = []
    for idx, col in enumerate(col_indexs):
        hesc_lev = hesc_methy_levels[idx]
        other_lev = other_data_levels[idx, col]
        label = -1 * (col + 1) if hesc_lev > other_lev else (col + 1)
        labels.append(label)
    labels = np.array(labels)
    out_data[:, -1] = labels
    np.savetxt(out_fp, out_data[:, :], fmt="chr%d\t%d\t%d\t%d", delimiter='\n')
    print("finish %s" % out_fp)
def plot_cluster():
    K_Cluster_dir = os.path.join(DATA_SET_DIR, "K_Cluster_Analysis")
    data_fp = os.path.join(K_Cluster_dir, "Cluster_ChromoHMM_chr1.tsv")
    df = pd.read_csv(data_fp, sep='\t', header=None).values
    data_to_cluster = df[:, 3: -1]
    print(data_to_cluster.shape)
    n_neighbors = 10

    fig = plt.figure(figsize=(15, 8))
    plt.suptitle("Manifold Learning with %i points, %i neighbors"
                 % (1000, n_neighbors), fontsize=14)

    # Perform Isomap Manifold learning.
    t0 = time.time()
    trans_data = manifold.Isomap(n_neighbors, n_components=2) \
        .fit_transform(data_to_cluster).T
    t1 = time.time()
    print("%s: %.2g sec" % ('ISO', t1 - t0))

    ax = fig.add_subplot(251)
    plt.scatter(trans_data[0], trans_data[1], cmap=plt.cm.rainbow)
    plt.title("%s (%.2g sec)" % ('Isomap', t1 - t0))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis('tight')
def plot_ChrMM_of_Cluster(cluster_fps, CLASS_LABELS):
    fig_dir = "../FIGURES/Filtered_Feature_Analysis/ChromHMM_of_Cluster"
    mkdirs([fig_dir])
    N_FILES = len(cluster_fps)
    chrHMMs = np.arange(NUMBER_OF_CHRMM_CLASS)
    class_data = np.zeros((len(chrHMMs), N_FILES))
    for fid, cluster_fp in enumerate(cluster_fps):
        df = pd.read_csv(cluster_fp, sep='\t', header=None).values
        data_HMMs = df[:, -1].astype(int)
        unique, counts = np.unique(data_HMMs, return_counts=True)
        normed_counts = counts / float(counts.sum())
        unique_dc = dict(zip(unique, normed_counts))
        for cls in chrHMMs:
            class_data[cls][fid] = unique_dc[cls + 1]
    cm = plt.get_cmap('gist_rainbow')
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ind = np.arange(N_FILES)
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
    ax.set_ylabel('Percentage')
    ax.set_xticks(ind)
    ax.set_xticklabels(CLASS_LABELS, fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.arange(0.2, 1.1, 0.2))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(tuple(pls), tuple(ChrmoHMM_LABELS), loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title('ChromHMM_High_K_CV_AND_Global', fontsize=18)
    plt.savefig(os.path.join(fig_dir, "ChromHMM_High_K_CV_AND_Global.png"), dpi=200)


def random_select_bed_lines(input_fp, output_fp):
    file_sz = os.path.getsize(input_fp) / 1e6 # get the file size in mb
    MAX_FILE_SZ = 4.8
    if file_sz > MAX_FILE_SZ:
        ratio = MAX_FILE_SZ / file_sz
    else:
        ratio = 1.
    with open(input_fp, "r") as input_f:
        lines = input_f.readlines()
        N_LINES = len(lines)
        SELECT_LINE_COUNTS = int(ratio * N_LINES)
        selected_indexs = list(np.random.choice(N_LINES, SELECT_LINE_COUNTS,replace=False))
        lines = np.array(lines)
        selected_lines = lines[selected_indexs]
    with open(output_fp, "w") as out_f:
        for line in selected_lines:
            out_f.write(line)
        print("write %s sucessful" % output_fp)
def filter_the_bed_file_to_max_size(INPUT_DIR, OUTPUT_DIR):
    mkdirs([OUTPUT_DIR])
    for file_name in os.listdir(INPUT_DIR):
        if file_name.endswith(".bed"):
            input_fp = os.path.join(INPUT_DIR, file_name)
            output_fp = os.path.join(OUTPUT_DIR, file_name)
            random_select_bed_lines(input_fp, output_fp)
if __name__ == "__main__":

    # cluster_f1 = os.path.join(K_Cluster_dir, "Methylated_And_High_K_Coef_Var_Cluster.tsv")
    # cluster_f2 = os.path.join(K_Cluster_dir, "Cluster_ChromoHMM.tsv")
    # names = ['Methylated_And_High_K_Coef_Var_Cluster', 'Cluster_ChromoHMM']
    # plot_ChrMM_of_Cluster([cluster_f1, cluster_f2], names)
    INPUT_DIR = os.path.join(DATA_SET_DIR, "DATA", "PERCENTILED_K")
    OUTPUT_DIR = os.path.join(DATA_SET_DIR, "DATA", "FILTERED_PERCENTILED_K")
    filter_the_bed_file_to_max_size(INPUT_DIR, OUTPUT_DIR)