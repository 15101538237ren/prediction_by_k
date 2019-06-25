# -*- coding: utf-8 -*-

import os, re
import numpy as np
def read_bed_file_and_convert(input_file_path, out_file_path, start_line= 1, sep="\s+",line_end = "\n"):
    '''
    :param input_file_path:
    :param out_file_path:
    :param start_line:
    :param sep:
    :param line_end:
    :return:
    '''
    line_counter = 0
    ltws = []
    with open(input_file_path, "r") as input_file:
        print("processing %s" % input_file_path)
        line = input_file.readline()
        while line:
            line_counter += 1
            if line_counter >= start_line:
                if line_counter % 500000 == 0:
                    print( "%d lines processed" % line_counter)
                line_contents = re.split(sep, line.strip(line_end))
                try:
                    chr_i, start, end, fraction, _ = line_contents
                    start = int(start) + 1
                    end = start + 1
                    methy_reads, tot_reads = re.split('/', fraction.strip('\''))
                    ltws.append('\t'.join([chr_i, str(start), str(end), methy_reads, str(int(tot_reads) - int(methy_reads))]) )
                except ValueError as e:
                    pass
            line = input_file.readline()

    with open(out_file_path, "w") as out_file:
        out_file.write((line_end).join(ltws))
        out_file.write(line_end)
    print('convert %s successful' % input_file_path)

def convert_file_format_pipeline():
    base_dir = '/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/REPLI_BS/'
    output_dir = '/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Output/'
    for file in [f for f in os.listdir(base_dir)]:
        if file.endswith('.bed'):
            read_bed_file_and_convert(os.path.join(base_dir, file),
                                      os.path.join(output_dir, file))

def filter_bed_by_methylation_level(input_fp, output_fp, methy_threshold = 0.5):
    line_counter = 0
    with open(output_fp, "w") as output_file:
        with open(input_fp, "r") as input_file:
            line = input_file.readline()
            while line:
                line_counter += 1
                if line_counter % 200000 == 0:
                    print( "%d lines processed" % line_counter)
                line_contents = re.split(r"\s+", line.strip("\n"))
                try:
                    chr_no, start_site, end_site, cpg_id, fraction = line_contents
                    if float(fraction) > methy_threshold:
                        output_file.write(line)
                except ValueError as e:
                    pass
                line = input_file.readline()
def calculate_methy_cpg_reads_ratio(bed_file_path, window_sz= 0):
    tot_cpg_reads = np.array([0.])
    methylated_cpg_reads = np.array([0.])

    cpg_counter = 0
    window_idx = 0
    with open(bed_file_path, "r") as input_file:
        line = input_file.readline()
        while line:
            line_contents = re.split(r"\s+", line.strip("\n"))
            try:
                _, _, _, m_reads, tot_reads, _ = line_contents
                if window_sz:
                    if cpg_counter and (cpg_counter % window_sz == 0):
                        window_idx += 1
                        tot_cpg_reads = np.append(tot_cpg_reads, 0)
                        methylated_cpg_reads = np.append(methylated_cpg_reads, 0)

                    cpg_counter += 1
                methylated_cpg_reads[window_idx] += int(m_reads)
                tot_cpg_reads[window_idx] += int(tot_reads)
            except ValueError as e:
                pass
            line = input_file.readline()
    return float(np.mean(methylated_cpg_reads / tot_cpg_reads))

def calculate_methy_cpg_reads_ratio_by_fraction(bed_file_path, window_sz= 0):
    cpg_counter = 0
    window_idx = 0

    fraction_array = []
    tmp_arr = []
    with open(bed_file_path, "r") as input_file:
        line = input_file.readline()
        while line:
            line_contents = re.split(r"\s+", line.strip("\n"))
            try:
                fraction = float(line_contents[4])
                if window_sz:
                    if cpg_counter and (cpg_counter % window_sz == 0):
                        window_idx += 1
                        fraction_array.append(np.mean(np.array(tmp_arr)))
                        tmp_arr = []

                    cpg_counter += 1
                    tmp_arr.append(fraction)
                else:
                    fraction_array.append(fraction)
            except ValueError as e:
                pass
            line = input_file.readline()
    return float(np.mean(np.array(fraction_array)))
def calc_ratio():
    base_dir = '/Users/Ren/PycharmProjects/prediction_by_k/DATA_SET/REPLI_BS-METHYLATION_ANALYSIS/'
    out_fp = os.path.join(base_dir, 'methylated_reads_ratio.tsv')
    with open(out_fp, 'w') as ouf_file:
        ouf_file.write("Px\tWindow size\tRepli-RRBS\tRepli-BS\tGEO\n")
        for item in ['P11', 'P12', 'P13']:
            for window_sz in [10, 100, 1000, 0]:

                repli_rrbs_ratio = calculate_methy_cpg_reads_ratio(os.path.join(base_dir, item, 'RRBS_its.bed'),
                                                             window_sz=window_sz)
                repli_bs_ratio = calculate_methy_cpg_reads_ratio(os.path.join(base_dir, item, item + '_its.bed'),
                                                           window_sz=window_sz)
                wgbs_bulk_ratio = calculate_methy_cpg_reads_ratio_by_fraction(os.path.join(base_dir, item, 'Methy_its.bed'),
                                                             window_sz=window_sz)

                ouf_file.write("%s\t%d\t%.3f\t%.3f\t%.3f\n" % (item, window_sz, repli_rrbs_ratio , repli_bs_ratio, wgbs_bulk_ratio))
                print( "%s\t%d\t%.3f\t%.3f\t%.3f" % (item, window_sz, repli_rrbs_ratio , repli_bs_ratio, wgbs_bulk_ratio))
if __name__ == "__main__":
    convert_file_format_pipeline()
    # calc_ratio()
    # filter_bed_by_methylation_level('/Users/Ren/PycharmProjects/prediction_by_k/DATA_SET/REPLI_BS-METHYLATION_ANALYSIS/GSM1112841.bed','/Users/Ren/PycharmProjects/prediction_by_k/DATA_SET/REPLI_BS-METHYLATION_ANALYSIS/GSM1112841_filtered.bed')