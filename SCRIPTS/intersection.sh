#!/bin/bash
#$ -N Intersection
#$ -q pub*,ionode,rxn
#$ -m beas

#hESC, dEN, dEC, dME intersection with RRBS CpGs

module load bedtools/2.25.0
module load bedops/2.4.14

work_dir=/pub/hongleir/DATA_FOR_ANALYSIS
cd $work_dir

f1=repli_bs/RRBS_locations_S1toS6fraction_1hBrdu_0Chase_nascent.bed
sort -k 1,1 -k2,2n $f1> $f1.sorted
rm -f $f1
mv $f1.sorted $f1
f2=bulk_rrbs/GSM1112820.bed
awk 'NR%2!=0' $f2|sort -k 1,1 -k2,2n> $f2.sorted
rm -f $f2
mv $f2.sorted $f2
f3=bulk_rrbs/GSM1112839.bed
awk 'NR%2!=0' $f3|sort -k 1,1 -k2,2n> $f3.sorted
rm -f $f3
mv $f3.sorted $f3
f4=bulk_rrbs/GSM1112841.bed
#awk 'NR%2!=0' $f4|sort -k 1,1 -k2,2n> $f4.sorted
rm -f $f4
mv $f4.sorted $f4
f5=bulk_rrbs/GSM916051.bed
#awk 'NR%2!=0' $f5|sort -k 1,1 -k2,2n> $f5.sorted
rm -f $f5
mv $f5.sorted $f5

intersect_out=bulk_rrbs/rrbs_cpgs_intersect.bed

bedtools intersect -a $f1 -b $f2 -sorted> $f2.f1.intersect
bedtools intersect -a $f2.f1.intersect -b $f3 -sorted>$f3.f2.intersect
bedtools intersect -a $f3.f2.intersect -b $f4 -sorted>$f4.f3.intersect
bedtools intersect -a $f4.f3.intersect -b $f5 -sorted>$intersect_out
rm -f $f2.f1.intersect $f3.f2.intersect $f4.f3.intersect

bedtools intersect -a $f2 -b $intersect_out -wa> $f2.intersected
bedtools intersect -a $f3 -b $intersect_out -wa> $f3.intersected
bedtools intersect -a $f4 -b $intersect_out -wa> $f4.intersected
bedtools intersect -a $f5 -b $intersect_out -wa> $f5.intersected

max_file=bulk_rrbs/max_hESC_dME_dEC_dEN.bed
max_file_filtered=bulk_rrbs/max_hESC_dME_dEC_dEN_0.5_filtered.bed

cat $f2.intersected $f3.intersected $f4.intersected $f5.intersected|sort -k 1,1 -k2,2n|bedtools merge -c 5 -o max> $max_file

awk '$4>0.5' $max_file> $max_file_filtered

bedtools intersect -a $f2.intersected -b $max_file_filtered -wa> $f2.inter_filtered
bedtools intersect -a $f3.intersected -b $max_file_filtered -wa> $f3.inter_filtered
bedtools intersect -a $f4.intersected -b $max_file_filtered -wa> $f4.inter_filtered
bedtools intersect -a $f5.intersected -b $max_file_filtered -wa> $f5.inter_filtered


mkdir K

cp ../prediction_by_k/DATA/K.bed K

bedtools intersect -a K/K.bed -b $max_file -wa> K/K_not_filtered.bed
bedtools intersect -a K/K.bed -b $max_file_filtered -wa> K/K_filtered_by_0.5.bed

intersect_dir=hESC_dME_dEC_dEN_interesected

cp K/* $max_file $max_file_filtered bulk_rrbs/*.intersected bulk_rrbs/*.inter_filtered $intersect_dir

tar -jcvf $intersect_dir.tar.bz2 $intersect_dir


scp hongleir@hpc.oit.uci.edu:/pub/hongleir/DATA_FOR_ANALYSIS/hESC_dME_dEC_dEN_interesected.tar.bz2 /Users/Ren/PycharmProjects/prediction_by_k/DATA/DATA_FOR_ANALYSIS/


