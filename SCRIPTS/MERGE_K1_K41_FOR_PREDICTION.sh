#!/bin/bash

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/
PREDICT_PROJ_DIR=/Users/emmanueldollinger/PycharmProjects/prediction_by_k
PREPROCESS_DIR=$PREDICT_PROJ_DIR/DATA/Preprocessing
DMR_DIR=$PREDICT_PROJ_DIR/DATA/Genomic_Features/DMR
DATA_DIR=$BASE_DIR/DATA/Repli_BS
K_RATES_DIR=$DATA_DIR/K_RATES

f1=$K_RATES_DIR/1.bed
f2=$K_RATES_DIR/41.bed
f3=$PREPROCESS_DIR/Merged_k1_k41.bed
f4=/pub/hongleir/data/Methy-data
f5=$PREPROCESS_DIR/Merged_k1_k41_max.bed
f6=$PREPROCESS_DIR/Merged_k1_k41_min.bed
f7=$PREPROCESS_DIR/GENOME_100bp_TILES.bed
f10=$PREPROCESS_DIR/merged_dmr.bed
f11=$PREPROCESS_DIR/dmr_ovlp_tile.bed
f12=$PREPROCESS_DIR/non_dmr_ovlp_tiles.bed
f13=$PREPROCESS_DIR/dmr_merged_tile.bed
f14=$PREPROCESS_DIR/GSM1112841_intersected_with_merged_k.bed
# 1. Merge K1 K41 together and sort the merged file: 
cat $f1 $f2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f3
# Get the min and max file
bedtools merge -c 4,5 -o max -i $f3>$f5
bedtools merge -c 4,5 -o min -i $f3>$f6

# 2. Interesect the merged k sites with the methylation file in hESC
scp $f3 hongleir@hpc.oit.uci.edu:$f4
# Server side : operation
qrsh
module load bedtools/2.25.0
cd $f4
bedtools intersect -a GSM1112841.bed -b Merged_k1_k41.bed -wa| sort -k 1,1 -k2,2n|bedtools merge -c 5 -o max>GSM1112841_intersected_with_merged_k.bed
# download the intersected sites of methylation data
scp hongleir@hpc.oit.uci.edu:$f4/GSM1112841_intersected_with_merged_k.bed $PREPROCESS_DIR

# now we have 3 files: Merged_k1_k41_max.bed Merged_k1_k41_min.bed GSM1112841_intersected_with_merged_k.bed


# 3. merge DMR file and label each record by cell type class
cd $DMR_DIR
FILE_COUNTER=0
for f in *.bed
do
  FILE_COUNTER=$((FILE_COUNTER+1))
  echo $FILE_COUNTER
  awk -v fc=$FILE_COUNTER 'BEGIN {FS="\t"; OFS=","} {print $0"\t"fc}' $f| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f.tmp
done

cat *.tmp | gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$DMR_DIR/merged_dmr.bed
rm *.tmp


#now have merged_dmr.bed

# Stat (histogram) the dmr region length
awk 'BEGIN {FS="\t"; OFS=","} {print $3-$2}' $DMR_DIR/merged_dmr.bed>$DMR_DIR/region_length_of_merged_dmr.txt

# After histogram, we found the most frequent DMR size is 100. 
# Therefore, we partition the Genome by 100bp
cd $PREPROCESS_DIR

#wget https://raw.githubusercontent.com/arq5x/bedtools/master/genomes/human.hg19.genome
#bedtools makewindows -g human.hg19.genome -w 100 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$f7
# make the tile file without overlap
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3-1}' $f7>$f7".tmp"
rm -f $f7
mv $f7".tmp" $f7
#Filter the tiles by max neighbor distance at 10000
bedtools closest -a $f7 -b $f6 -D a -t first -k 5|gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f7".closest5"
awk -F "\t" '{ if($9<10000) { print } }' $f7".closest5"| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$9}'| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4 -o max>$f7".filtered"
rm $f7".closest5"
awk -F "\t" '{ if($4>-10000) { print } }' $f7".filtered"| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f7".filtered_n"
rm $f7".filtered"
mv $f7".filtered_n" $f7".filtered"


#cp $DMR_DIR/merged_dmr.bed $f10
# -a tile file, -b merged_dmr, the minimum fraction be satisfied at least 50% of A (50bp) is covered
# if the overlapping between the genome tile and DMR are over 50bp (50 percent), then the tile is labeled as DMR region.

bedtools intersect -a $f7".filtered" -b $f10 -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f11
#remove the useless columns
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8}' $f11>$f11".tmp"
rm -f $f11
mv $f11".tmp" $f11

#obtain the non-dmr overlapping tile file

bedtools intersect -a $f7".filtered" -b $f10 -f 0.5 -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f12
# add 0 to last column to indicate the class is 0 relative to class label in dmr
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"0}' $f12>$f12".tmp"
rm -f $f12
mv $f12".tmp" $f12

cat $f11 $f12| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4 -o max>$f13
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f13>$f13".index"


#can parallel with the last one
# For each tile, find the nearest 1, 5, 10, 20 Ks and the corresponding distance by bedtools closest.
# Overlap the genome tiles with DMR
prefix='chr'

BASE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/
PREDICT_PROJ_DIR=/Users/emmanueldollinger/PycharmProjects/prediction_by_k
PREPROCESS_DIR=$PREDICT_PROJ_DIR/DATA/Preprocessing
DMR_DIR=$PREDICT_PROJ_DIR/DATA/Genomic_Features/DMR
DATA_DIR=$BASE_DIR/DATA/Repli_BS
K_RATES_DIR=$DATA_DIR/K_RATES

f1=$K_RATES_DIR/1.bed
f2=$K_RATES_DIR/41.bed
f3=$PREPROCESS_DIR/Merged_k1_k41.bed
f4=/pub/hongleir/data/Methy-data
f5=$PREPROCESS_DIR/Merged_k1_k41_max.bed
f6=$PREPROCESS_DIR/Merged_k1_k41_min.bed
f7=$PREPROCESS_DIR/GENOME_100bp_TILES.bed
f10=$PREPROCESS_DIR/merged_dmr.bed
f11=$PREPROCESS_DIR/dmr_ovlp_tile.bed
f12=$PREPROCESS_DIR/non_dmr_ovlp_tiles.bed
f13=$PREPROCESS_DIR/dmr_merged_tile.bed
f14=$PREPROCESS_DIR/GSM1112841_intersected_with_merged_k.bed
for N_neighbors in 5; # 3 5
do
	OUT_NDIR=$PREPROCESS_DIR/$N_neighbors
	f8=$OUT_NDIR/$N_neighbors"_K_min.bed"
	f9=$OUT_NDIR/$N_neighbors"_K_max.bed"
	f92=$OUT_NDIR/$N_neighbors"_methylation.bed"
	mkdir -p $OUT_NDIR
	#-a tile file, -b K min file
	bedtools closest -a $f7".filtered" -b $f6 -D a -t first -k $N_neighbors|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10/10000.0}'|gsort -k 1,1 -k2,2n --parallel=8  -S 50%| bedtools merge -c 4,5,6 -o collapse>$f8
	#awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5"\t"$6}' $f8 | gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$f8".csv"
	#echo $N_neighbors" K min is done"
	
	bedtools closest -a $f7".filtered" -b $f5 -D a -t first -k $N_neighbors|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10/10000.0}'|gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6 -o collapse>$f9
	#awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5"\t"$6}' $f9 | gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$f9".csv"
	#echo $N_neighbors" K max is done"
	

	bedtools closest -a $f7".filtered" -b $f14 -D a -t first -k $N_neighbors |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9/10000.0}'|gsort -k 1,1 -k2,2n --parallel=8  -S 50%| bedtools merge -c 4,5 -o collapse>$f92
	#awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}' $f92 | gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$f92".csv"
	echo $N_neighbors" methy is done"
done

f15=$PREPROCESS_DIR/3/3_K_min.bed
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f15| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f15".tmp"

# Check the order of chr index are consistent between the k neighbor version and dmr region version
diff $f15".tmp" $f13".index"
# No differenece Checked! Now can use $f13 as index
rm $f15".tmp"

prefix='chr'
awk 'BEGIN {FS="\t"; OFS=","} {print NR-1"\t"$1"\t"$2"\t"$4}' $f13 | gsed 's/\t/,/g'| gsed -e "s/^$prefix//" >$f13".index"


DATA_SET_DIR=$PREDICT_PROJ_DIR/DATA/DATA_SET

mv 1 3 5 $f13".index" $DATA_SET_DIR









