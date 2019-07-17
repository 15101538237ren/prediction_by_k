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
DATA_SET_DIR=$PREDICT_PROJ_DIR/DATA/DATA_SET
Genomic_Features_DIR=$PREDICT_PROJ_DIR/DATA/Genomic_Features
f16=$DATA_SET_DIR/Genomic_Features
ChromoHMM=$Genomic_Features_DIR/HmmH1hescHMM.bed
f17=$ChromoHMM".interesected"
f18=$ChromoHMM".non_interesected"
f19=$ChromoHMM".merged"
f20=$PREPROCESS_DIR/dmr_merged_tile_with_ChromeMM.bed
index_dir=$PREPROCESS_DIR/Indexs
GENOME_DIR=$PREPROCESS_DIR/hg19
TILE_METHY_DATA=$PREPROCESS_DIR/TILE_METHY_REP
prefix='chr'
declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

HISTONE_DIR=$PREDICT_PROJ_DIR/DATA/DATA_FOR_PREDICTION
declare -a HISTONES=('GSM772800' 'GSM772752' 'GSM772756' 'GSM997249' 'GSM772750' 'GSM772751')

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
  echo $f" "$FILE_COUNTER
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


mkdir -p $index_dir
for chr in "${CHRS[@]}"; 
do
    bedextract $chr $f13 | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4}'|gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $index_dir/$chr.csv; 
done

#awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $index_dir/chr1.bed>$index_dir/chr1.bed".index"


for N_neighbors in 1; # 3 5
do
	OUT_NDIR=$PREPROCESS_DIR/$N_neighbors
	f8=$OUT_NDIR/$N_neighbors"_K_min"
	f9=$OUT_NDIR/$N_neighbors"_K_max"
	f92=$OUT_NDIR/$N_neighbors"_methylation"
	mkdir -p $OUT_NDIR
	mkdir -p $f8 $f9 $f92
	mkdir -p $f8"_ind"
	#-a tile file, -b K min file
	#bedtools closest -a $f7".filtered" -b $f6 -D a -t first -k $N_neighbors|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10/10000.0}'|gsort -k 1,1 -k2,2n --parallel=8  -S 50%| bedtools merge -c 4,5,6 -o collapse>$f8".bed"
	for chr in "${CHRS[@]}"; 
	do
	    bedextract $chr $f8".bed" > $f8"_ind"/$chr".bed" #| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5"\t"$6}' | gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $f8/$chr".csv"; 
	done
	echo $N_neighbors" K min is done"
	
	#bedtools closest -a $f7".filtered" -b $f5 -D a -t first -k $N_neighbors|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10/10000.0}'|gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6 -o collapse>$f9".bed"
	
	for chr in "${CHRS[@]}"; 
	do
	    bedextract $chr $f9".bed" |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5"\t"$6}' |gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $f9/$chr".csv"; 
	done
	echo $N_neighbors" K max is done"

	# #bedtools closest -a $f7".filtered" -b $f14 -D a -t first -k $N_neighbors |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9/10000.0}'|gsort -k 1,1 -k2,2n --parallel=8  -S 50%| bedtools merge -c 4,5 -o collapse>$f92".bed"
	for chr in "${CHRS[@]}"; 
	do
	    bedextract $chr $f92".bed" | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}' | gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$f92/$chr".csv"; 
	done
	echo $N_neighbors" methy is done"
done

f15=$PREPROCESS_DIR/1/1_K_min_ind/chr1.bed
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f15| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f15".tmp"

# Check the order of chr index are consistent between the k neighbor version and dmr region version
diff $f15".tmp" $index_dir/chr1.bed".index"
# No differenece Checked! Now can use $f13 as index
rm $f15".tmp"

prefix='chr'
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4}' $f13 | gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$f13".index" 




cp -r $index_dir $DATA_SET_DIR/Indexs

# Process ChromoHMM, first call preprocess_ChromeHMM() in Python file to format it into 4 col .bed
gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $ChromoHMM>$ChromoHMM".sorted"
rm $ChromoHMM
mv $ChromoHMM".sorted" $ChromoHMM

bedtools intersect -a $f7".filtered" -b $ChromoHMM -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f17
#remove the useless columns
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8}' $f17>$f17".tmp"
rm -f $f17
mv $f17".tmp" $f17

#obtain the non-dmr overlapping tile file

bedtools intersect -a $f7".filtered" -b $ChromoHMM -f 0.5 -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f18
# add 0 to last column to indicate the class is 0 relative to class label in dmr
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"0}' $f18>$f18".tmp"
rm -f $f18
mv $f18".tmp" $f18

cat $f17 $f18| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4 -o max>$f19
rm $f17 $f18

awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f19>$f19".tmp"

for chr in "${CHRS[@]}"; 
do
    bedextract $chr $f13 | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4}'|gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $index_dir/$chr.csv; 
done

awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f13>$f13".tmp"

diff $f13".tmp" $f19".tmp"
rm -f $f13".tmp" $f19".tmp"


#Process Tile K data for ML: copied from K_Estimation/DATA/Repli_BS/K_RATES/41/100/41.bed and 1.bed
for ti in 1 2;
do
	f21=$TILE_METHY_DATA$ti".bed"
	echo $f21
	#gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f21>$f21".sorted"
	#rm $f21
	#mv $f21".sorted" $f21
	bedtools intersect -a $f7".filtered" -b $f21 -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f21".tmp"
	#remove the useless columns
	awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"\t"$9"\t"$10}' $f21".tmp">$f21".intersected"

	#obtain the non-dmr overlapping tile file

	bedtools intersect -a $f7".filtered" -b $f21 -f 0.5 -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f21".tmp"
	# add 0 to last column to indicate the class is 0 relative to class label in dmr
	awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"0"\t"0"\t"0}' $f21".tmp">$f21".non_intersected"
	rm -f $f21".tmp"
	cat $f21".intersected" $f21".non_intersected"| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6 -o max>$f21".merged"
	rm $f21".non_intersected" $f21".intersected"
done
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f21".merged">$f21".tmp"
diff $f13".tmp" $f21".tmp" 


for his in "${HISTONES[@]}"; 
do
	f22=$HISTONE_DIR/$his".bed"
	echo $f22

	bedtools intersect -a $f7".filtered" -b $f22 -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f22".tmp"
	#remove the useless columns
	awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$8"}' $f22".tmp">$f22".intersected"

	#obtain the non-dmr overlapping tile file

	bedtools intersect -a $f7".filtered" -b $f22 -f 0.5 -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f22".tmp"
	# add 0 to last column to indicate the class is 0 relative to class label in dmr
	awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"0}' $f22".tmp">$f22".non_intersected"
	rm -f $f22".tmp"
	cat $f22".intersected" $f22".non_intersected"| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4 -o max>$f22".merged"
	rm $f22".non_intersected" $f22".intersected"    
done


paste $f13 $f19 | cut -f 1,2,3,4,8 >$f20
for chr in "${CHRS[@]}"; 
do
    bedextract $chr $f20 | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}'|gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $index_dir/$chr.csv; 
done





# download hg19
cd $GENOME_DIR
FTP_LINK=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
for chr in "${CHRS[@]}"; 
do
	file_name=$chr".fa.gz"
    wget $FTP_LINK$file_name .
    gunzip $file_name
done

f1=$K_RATES_DIR/55.bed
ffn=$PREPROCESS_DIR/merged_k_and_methy.bed

# gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f1>$f1".sorted"
# rm -f $f1
# mv $f1".sorted" $f1
# gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f14>$f14".sorted"
# rm -f $f14
# mv $f14".sorted" $f14
bedtools intersect -a $f1 -b $f14 -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$ffn
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$9}' $ffn|gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$ffn".csv"

