#!/bin/bash
#$ -N Calc_methy_ratio
#$ -q pub*,ionode,rxn
#$ -m beas

# 1. filter the bulk CpGs with 0.5 methylation level
# bulk_hesc_methy_data=/pub/hongleir/data/Methy-data/GSM1112841.bed
# use REPLI_BS_AND_RRBS.py: call the filter_bed_by_methylation_level($bulk_hesc_methy_data, filtered_bulk_hesc_methy_data)

# 2. packaging repli-rrbs data

cp /pub/hongleir/data/Methy-data/GSM916051.bed bulk_rrbs/

packaging_dir=/pub/hongleir/DATA_FOR_ANALYSIS/repli_rrbs

mkdir $packaging_dir
cd $packaging_dir

cp /pub/hongleir/DATA_FOR_ANALYSIS/P1*/results/bismark_methylation_calls/methylation_coverage/P1*_READ1_val_1_bismark_bt2_pe.bismark.cov.gz .

gunzip *.gz

# do some formatting by regular expression with sublime to convert cov format into standard bed format

cd ..

tar -jcvf 'repli_rrbs_bed.tar.bz2' $packaging_dir

# Right now, we should have:
# 1. All  Repli_rrbs data files: P11 (20.2MB), P12 (357.3MB), P13 (233.5MB) 
# 2. RRBS-CpG-filtered Repli-BS data: RRBS_locations_S1toS6fraction_1hBrdu_0Chase_nascent.bed
# 3. Bulk RRBS data: GSM1112841_filtered.bed

# Then for each P1* data:
# we intesect P1* with RRBS_locations_S1toS6fraction_1hBrdu_0Chase_nascent.bed to get overlapping CpGs> INTERSECT.bed, P1*_overlapped.bed, Repli_bs_overlapped.bed
# then intersect Bulk RRBS data with reduced overlapped CpGs> Methy_its.bed
# Finally, by applying the calc_ratio() in REPLI_BS_AND_RRBS.py, we got the methy fraction for these three different data types.


