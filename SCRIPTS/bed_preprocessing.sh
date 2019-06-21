#!/bin/bash
#$ -N PEAK_CALLING
#$ -q pub*,ionode,rxn
#$ -m beas

module load bedtools/2.25.0
module load bedops/2.4.14

indir=/pub/hongleir/data/ChIP-Seq_Peaks
K_dir=/pub/hongleir/data/K

peak_merged=$indir/peak_merged.bed
peaks_overlapped_with_K_filtered=$indir/peaks_overlapped_with_K_filtered.bed
original_K_data=$K_dir/K.bed
filtered_K_data=$K_dir/K_filtered_by_merged_peaks.bed
K_mean_in_each_peak=$K_dir/K_mean_in_each_peak.bed


# Chip-Seq Data data Processing
# 1. Merge all identified peaks in data
cat $indir/*island'.bed'| sort -k 1,1 -k2,2n | bedtools merge>$peak_merged
# 2. Obtain the overlapped CpGs in K data and Peaks
bedtools intersect -a $original_K_data -b $peak_merged -wa | sort -k 1,1 -k2,2n>$filtered_K_data
bedtools intersect -a $peak_merged -b $filtered_K_data -wa | sort -k 1,1 -k2,2n| bedtools merge>$peaks_overlapped_with_K_filtered
# 3. Calculate the average K in each peak region
bedtools map -a $peaks_overlapped_with_K_filtered -b $filtered_K_data -c 4 -o mean> $K_mean_in_each_peak

cd $indir

for f in *_island.bed
do
  filename=${f%%_*}
  bedtools intersect -a $peaks_overlapped_with_K_filtered -b $f -wa | sort -k 1,1 -k2,2n| bedtools merge | sed 's/^\([^ ]*\)/\1 1/' >$filename'_overlap.bed'
  bedtools intersect -a $peaks_overlapped_with_K_filtered -b $f -v | sort -k 1,1 -k2,2n| bedtools merge | sed 's/^\([^ ]*\)/\1 0/' >$filename'_no_overlap.bed'
  cat $filename'_overlap.bed' $filename'_no_overlap.bed' | sort -k 1,1 -k2,2n >$filename'_binarized.bed'
  rm -f $filename'_overlap.bed' $filename'_no_overlap.bed'
  echo $filename
done


cd $indir
DATA_FOR_PREDICTION=/pub/hongleir/DATA_FOR_PREDICTION
mkdir -p $DATA_FOR_PREDICTION
for f in *'_binarized.bed'
do
  filename=${f%%_*}
  cp $f $DATA_FOR_PREDICTION
  mv $DATA_FOR_PREDICTION/$filename'_binarized.bed' $DATA_FOR_PREDICTION/$filename'.bed'
  echo $filename
done

# Methylation BED data Processing
methy_data_dir=/pub/hongleir/data/Methy-data
cd $methy_data_dir

for f in *.bed
do
    filename=${f%%\.*}
    echo $filename
    bedtools intersect -a $f -b $peaks_overlapped_with_K_filtered | sort -k 1,1 -k2,2n >$filename'.tmp'
    bedtools map -a $peaks_overlapped_with_K_filtered -b $filename'.tmp' -c 5 -o mean>$DATA_FOR_PREDICTION/$filename'.bed'
    rm -f $filename'.tmp'
done

cp $K_mean_in_each_peak $DATA_FOR_PREDICTION

tar -jcvf $DATA_FOR_PREDICTION'.tar.bz2' $DATA_FOR_PREDICTION
module purge