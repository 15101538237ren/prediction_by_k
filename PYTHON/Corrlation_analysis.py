# -*- coding:utf-8 -*-
import math,os,collections
import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True)
GENOMIC_FEATURES_DIR = '/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Genomic_Features/'
EACH_SUB_FIG_SIZE = 5
def mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def _workers_count():
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_count = os.cpu_count()
    return cpu_count

def read_bed_file_and_store_pos_to_a_struct(bed_tsv_file_path):
    """
    Read tsv format bed file (No chr prefix at 1st colum), and store the chr, pos, k or methylation level info into dict
    :param bed_tsv_file_path: tsv format bed file input path
    :return: the dict with chr, pos, k or methylation level info
    """
    df = pd.read_csv(bed_tsv_file_path, sep='\t', header=None).values
    chrs = df[:, 0]
    unique_chr, _ = np.unique(chrs, return_counts=True)
    dict_to_store = {chr_i: {} for chr_i in unique_chr}
    for item in df:
        chr_i, pos, value, _ = item
        dict_to_store[chr_i][int(pos)] = value
    return [dict_to_store]

def filter_d_length_to_generate_CpG_pairs(Chr_CpG_pos_and_methy_dict, D_MAX):
    """
    Extract the CpG pairs with distance of d, and there can be CpGs between CpG pairs with d distance
    :param Chr_CpG_pos_and_methy_dict: The dict with chr, pos, k or methylation level info provided
    :param d: the distance that required for CpG pairs
    :return: An array that stores all CpG pairs with distance of d
    """
    CpG_Pairs_Dict = {}
    D_RANGE = range(2, D_MAX)
    Chr_CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[0]
    for d in D_RANGE:
        if d % 10 == 0:
            print("Calculating CpG pairs with %d distance" % d)
        for chr_i in Chr_CpG_pos_and_methy_dict.keys():
            CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[chr_i]
            for key in CpG_pos_and_methy_dict.keys():
                pos = key
                methy_level_1 = CpG_pos_and_methy_dict[key]
                pos_off_by_d = pos + d
                if pos_off_by_d in CpG_pos_and_methy_dict.keys():
                    methy_level_2 = CpG_pos_and_methy_dict[pos_off_by_d]
                    methy_levels_of_CpG_pairs =[methy_level_1, methy_level_2]
                    if d not in CpG_Pairs_Dict.keys():
                        CpG_Pairs_Dict[d] = []
                    CpG_Pairs_Dict[d].append(methy_levels_of_CpG_pairs)
    return [CpG_Pairs_Dict]

def filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg(Chr_CpG_pos_and_methy_dict, D_MAX):
    """
    Extract the CpG pairs with distance of d, and required no CpGs between CpG pairs with d distance
    :param d: the distance that required for CpG pairs
    :param length_of_od:
    :param od_keys:
    :param od_vals:
    :return:
    """
    CpG_Pairs_list = {}
    Chr_CpG_pos_and_methy_dict=Chr_CpG_pos_and_methy_dict[0]
    for chr_i in Chr_CpG_pos_and_methy_dict.keys():
        CpG_pos_and_methy_dict = Chr_CpG_pos_and_methy_dict[chr_i]
        sorted_struct = collections.OrderedDict(sorted(CpG_pos_and_methy_dict.items()))
        od_keys = list(sorted_struct.keys())
        od_vals = list(sorted_struct.values())
        for i in range(0, len(sorted_struct) - 1):
            # 后一个位点与前一个位点之间的距离为d
            pre_od_key = od_keys[i]
            post_od_key = od_keys[i + 1]
            dist = post_od_key - pre_od_key
            if 0 < dist <D_MAX:
                methy_levels_of_CpG_pairs = [od_vals[i], od_vals[i + 1]]
                if dist not in CpG_Pairs_list.keys():
                    CpG_Pairs_list[dist] = []
                CpG_Pairs_list[dist].append(methy_levels_of_CpG_pairs)
    return [CpG_Pairs_list]

def calc_C_d_by_pearson_correlation(CpG_Pairs_Dict):
    """
    calc the Pearson correlation of CpG pairs list
    :param CpG_pairs: The input CpG pairs list
    :return: The dict with {d: correlation}
    """
    RD_dict = {}
    CpG_Pairs_Dict = CpG_Pairs_Dict[0]
    for DISTANCE, CpG_pairs in CpG_Pairs_Dict.items():
        if len(CpG_pairs) > 5:
            CpG_arr = np.array(CpG_pairs)
            r_d = np.corrcoef(CpG_arr[:, 0], CpG_arr[:, 1])[0][1]
            RD_dict[DISTANCE] = r_d
    return [RD_dict] if RD_dict else []

def write_RD_into(RD_dict, out_fp):
    if RD_dict:
        RD_dict = RD_dict[0]
        sorted_struct = collections.OrderedDict(sorted(RD_dict.items()))
        od_keys = list(sorted_struct.keys())
        ltw = "\n".join(["%d\t%.3f" %(d, RD_dict[d]) for d in od_keys]) + "\n"
        with open(out_fp, 'w') as RD_file:
            RD_file.write(ltw)
def calc_correlation(bed_tsv_file_path, out_R_d_correlation_path,  d_max, is_inter_with_other_cpg):
    [Chr_CpG_pos_and_methy_dict] = read_bed_file_and_store_pos_to_a_struct(bed_tsv_file_path)
    if is_inter_with_other_cpg:
        CpG_pairs = filter_d_length_to_generate_CpG_pairs([Chr_CpG_pos_and_methy_dict], d_max)
    else:
        CpG_pairs = filter_d_length_to_generate_CpG_pairs_not_inter_with_other_cpg([Chr_CpG_pos_and_methy_dict], d_max)
    RD_dict = calc_C_d_by_pearson_correlation(CpG_pairs)
    write_RD_into(RD_dict, out_R_d_correlation_path)
def prepare_config_file(CONFIG_FP):
    regions = ["ChrmoHMM", "ADDED", "General", "Genomic_Regions", "Histone_Modification", "Histone_Modification_Enzyme"]  # ,"TFBS"
    k_or_methy = ['K', 'Methy']
    ltws = []
    for REGION in regions:
        for km in k_or_methy:
            base_dir = os.path.join(GENOMIC_FEATURES_DIR, REGION, km + "_intersected_extracted")
            out_rd_corr_dir = os.path.join(GENOMIC_FEATURES_DIR, REGION, km + "_RD")
            mkdir(out_rd_corr_dir)
            for file_name in os.listdir(base_dir):
                if file_name.endswith(".csv"):
                    bed_tsv_file_path = os.path.join(base_dir, file_name)
                    out_R_d_correlation_path = os.path.join(out_rd_corr_dir, file_name)
                    ltws.append('\t'.join([bed_tsv_file_path, out_R_d_correlation_path]))
    with open(CONFIG_FP, 'w') as CONFIG_file:
        ltw = '\n'.join(ltws) + '\n'
        CONFIG_file.write(ltw)
def plot_local_regression_and_RD():
    regions = ["ChrmoHMM", "ADDED", "General", "Genomic_Regions", "Histone_Modification", "Histone_Modification_Enzyme"]
    k_or_methy = ['K', 'Methy']
    N_COL = 5
    MAX_D = 2000
    cm = plt.get_cmap('gist_rainbow')
    FIG_DIR = os.path.join("../FIGURES/Correlation_Analysis/")
    for REGION in regions:
        for km in k_or_methy:
            fig_dir = os.path.join(FIG_DIR, km)
            mkdir(fig_dir)
            out_rd_corr_dir = os.path.join(GENOMIC_FEATURES_DIR, REGION, km + "_RD")
            file_labels = []
            file_paths = []
            for file_name in os.listdir(out_rd_corr_dir):
                if file_name.endswith(".csv"):
                    file_labels.append(file_name.strip(".csv"))
                    RD_path = os.path.join(out_rd_corr_dir, file_name)
                    file_paths.append(RD_path)
            fig_fp = os.path.join(fig_dir,  REGION + ".png")

            N_FILES = len(file_labels)
            N_ROW = int(math.ceil((N_FILES) / N_COL))
            fig, axs = plt.subplots(N_ROW, N_COL, figsize=(N_COL * EACH_SUB_FIG_SIZE, N_ROW * EACH_SUB_FIG_SIZE))
            for j in range(N_FILES):
                row = j // N_COL
                col = j % N_COL
                if N_ROW == 1:
                    ax = axs[col]
                else:
                    ax = axs[row][col]
                file_path = file_paths[j]
                RD_df = pd.read_csv(file_path, sep="\t", header=None).values
                sns.regplot(x=RD_df[:, 0], y=RD_df[:, 1], ax=ax, lowess=True, ci=95, color=cm(1. * j / N_FILES))
                ax.set_xticks(range(0, MAX_D + 1, 200))
                ax.set_xlim(0, MAX_D)
                ax.set_ylim(0, 1.0)
                ax.set_title(file_labels[j], fontsize=18)
            plt.savefig(fig_fp, dpi=200)
def wrapper_for_correlation_analysis(row):
    bed_tsv_file_path = str(row[0])
    out_R_d_correlation_path = str(row[1])
    print("%s\t%s" % (bed_tsv_file_path, out_R_d_correlation_path))
    D_MAX = 2000
    is_inter_with_other_cpg = True
    calc_correlation(bed_tsv_file_path, out_R_d_correlation_path, D_MAX, is_inter_with_other_cpg)
def call_for_correlation_analysis(CONFIG_FP):
    num_workers = _workers_count()
    print('Starting call_for_correlation_analysis with', num_workers, 'workers')
    parametersDF = pd.read_csv(CONFIG_FP, index_col=0, sep='\t', lineterminator='\n')
    tups = parametersDF.itertuples(name=None)
    pool = mp.Pool(processes=num_workers)
    pool.map_async(wrapper_for_correlation_analysis, tups).get()

if __name__ == "__main__":
    CONFIG_FP = '../DATA/RD_CONFIG.tsv'
    call_for_correlation_analysis(CONFIG_FP)