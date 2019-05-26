#!/bin/bash
#$ -N PEAK_CALLING
#$ -q pub*,ionode,rxn
#$ -m beas

module load bedtools/2.25.0
module load bedops/2.4.14

indir=/pub/hongleir/data/ChIP-Seq_Peaks
K_dir=prediction_by_k/data

chopped_file=$indir/MERGED'_200bp.bed'
merged_chopped=$indir/MERGED'_Chopped.bed'

# Chip-Seq Data data Processing

cat $indir/*island'.bed'| sort -k 1,1 -k2,2n | bedtools merge>$indir/MERGED'.bed'
bedops --chop 200 $indir/MERGED'.bed' >$chopped_file
bedtools intersect -a $K_dir/K.bed -b $chopped_file | sort -k 1,1 -k2,2n>$K_dir/K_filtered.bed
bedtools intersect -a $chopped_file -b $K_dir/K_filtered.bed -wa | sort -k 1,1 -k2,2n| bedtools merge>$indir/MERGED'_FILTERED_BY_METHYLATION.bed'
bedops --chop 200 $indir/MERGED'_FILTERED_BY_METHYLATION.bed' >$merged_chopped
rm $indir/MERGED'.bed' $chopped_file

bedtools map -a $merged_chopped -b $K_dir/K_filtered.bed -c 4 -o mean> $K_dir/K_mean.bed

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

# Methylation BED data Processing
for f in *.bed
do
    filename=${f%%\.*}
    bedtools intersect -a $f -b $merged_chopped | sort -k 1,1 -k2,2n >$filename'.tmp'
    bedtools map -a $merged_chopped -b $filename'.tmp' -c 5 -o mean>output/$filename'_mean.bed'
    rm -f $filename'.tmp'
    echo $filename
done

module purge