# -*- coding: utf-8 -*-

import os, re
def extract_promoter(bed_fp, out_fp, upstream = 2000, after_tss = 0):
    ltws = []
    with open(bed_fp, "r") as file:
        lines = [line for line in file]

        for line_idx, line in enumerate(lines):
            chr_i, start, end, strand = re.split("\s+", line.strip("\n"))[0:4]
            if strand =='+':
                start = int(start)
                end = start + after_tss
                start = start - upstream
            else:
                end = int(end)
                start = end - after_tss
                end = end + upstream

            ltw = "\t".join([chr_i, str(start), str(end), strand])
            ltws.append(ltw)

    with open(out_fp, "w") as file:
        file.write(("\n").join(ltws))
        file.write("\n")
        print("write %s successful!" % out_fp)

# Separate TFBS_Cluster file into separated file by TF name.
def separate_TFBS_by_name(bed_fp, out_dir):
    TF_dict = {}
    with open(bed_fp, "r") as file:
        lines = [line for line in file]

        for line_idx, line in enumerate(lines):
            _, _, _, tf_name = re.split("\s+", line.strip("\n"))[0:4]
            if tf_name not in TF_dict.keys():
                TF_dict[tf_name] =[]
            TF_dict[tf_name].append(line)
    for key, lines in TF_dict.items():
        out_fp = os.path.join(out_dir, key + ".bed")
        with open(out_fp, "w") as file:
            file.write(("").join(lines))
            print("write %s successful!" % out_fp)
