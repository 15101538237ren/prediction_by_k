#!/bin/bash
#$ -N PEAK_CALLING
#$ -q pub*,ionode,rxn
#$ -m beas

module load bedtools/2.25.0
module load bedops/2.4.14

indir=/pub/hongleir/data/ChIP-Seq_Peaks
merged_filtered=$indir/MERGED'_FILTERED_BY_METHYLATION.bed'
out_file=$indir/MERGED'.bed'
chopped_file=$indir/MERGED'_200bp.bed'
merged_chopped=$indir/MERGED'_Chopped.bed'

cat $indir/*island'.bed'| sort -k 1,1 -k2,2n | bedtools merge> $out_file
bedops --chop 200 $out_file >$chopped_file
bedtools intersect -a prediction_by_k/data/K.bed -b $chopped_file | sort -k 1,1 -k2,2n>prediction_by_k/data/K_filtered.bed
bedtools intersect -a $chopped_file -b prediction_by_k/data/K_filtered.bed -wa | sort -k 1,1 -k2,2n| bedtools merge>$merged_filtered
bedops --chop 200 $merged_filtered >$merged_chopped
rm $out_file $chopped_file

bedtools map -a $merged_chopped -b prediction_by_k/data/K_filtered.bed -c 4 -o mean> prediction_by_k/data/K_mean.bed

indir=/pub/hongleir/data/ChIP-Seq_Peaks
merged_filtered=$indir/MERGED'_FILTERED_BY_METHYLATION.bed'
cd $indir

for f in *_island.bed
do
  filename=${f%%_*}
  bedtools intersect -a $merged_chopped -b $f -wa | sed 's/^\([^ ]*\)/\1 1/' >$filename'_overlap.bed'
  bedtools intersect -a $merged_chopped -b $f -v | sed 's/^\([^ ]*\)/\1 0/' >$filename'_no_overlap.bed'
  cat $filename'_overlap.bed' $filename'_no_overlap.bed' | sort -k 1,1 -k2,2n >$filename'_binarized.bed'
  rm -f $filename'_overlap.bed'
  rm -f $filename'_no_overlap.bed'
  echo $filename
done

methy_data_dir=/pub/hongleir/data/Methy-data
cd $methy_data_dir
mkdir -p output

for f in *.bed
do
    filename=${f%%\.*}
    bedtools intersect -a $f -b $merged_chopped | sort -k 1,1 -k2,2n >$filename'.tmp'
    bedtools map -a $merged_chopped -b $filename'.tmp' -c 5 -o mean>output/$filename'_mean.bed'
    rm -f $filename'.tmp'
    echo $filename
done

module purge