#!/bin/bash

KRATE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS/K_RATES
FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis
DATA_DIR=$FEATURE_ANALYSIS/DATA

mkdir -p $FEATURE_ANALYSIS $DATA_DIR
# Steps
# 0. Prepare the methylation and K files
METHY_FILE=$DATA_DIR/GSM1112841_55_intersected.bed
K_FILE=$DATA_DIR/55.bed

#intersect the K with methy to get the overlapped K sites.
bedtools intersect -a $K_FILE -b $METHY_FILE -wa -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$K_FILE".tmp"
rm $K_FILE
mv $K_FILE".tmp" $K_FILE

# Check the order of data
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $METHY_FILE> $METHY_FILE".tmp"
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $K_FILE> $K_FILE".tmp"
diff $METHY_FILE".tmp" $K_FILE".tmp"

#if consistent with each other
rm $METHY_FILE".tmp" $K_FILE".tmp"

# 1. Separate K and methylation files by 0.8 methylation threshold.
METHY_K=$DATA_DIR/METHY_AND_K.bed
METHY_K_LESS=$DATA_DIR/METHY_K_LESS.bed #FILES THAT HAVE < 0.8 METHYLATION LEVEL AND ITS CORRESPONDING Ks
METHY_K_INTER=$DATA_DIR/METHY_K_INTER.bed
METHY_K_GREATER=$DATA_DIR/METHY_K_GREATER.bed
# Data ordering checked
paste $METHY_FILE $K_FILE | cut -f 1,2,3,4,8,9 >$METHY_K
# Process Methylation > 0.8
awk  -F "\t" '{ if($4 < 0.4) { print } }' $METHY_K > $METHY_K_LESS
awk  -F "\t" '{ if($4 > 0.4 && $4 < 0.8) { print } }' $METHY_K > $METHY_K_INTER
awk  -F "\t" '{ if($4 > 0.8) { print } }' $METHY_K > $METHY_K_GREATER
METHY_K_CSV=$DATA_DIR/METHY_AND_K.tsv
METHY_K_LESS_CSV=$DATA_DIR/METHY_K_LESS.tsv #FILES THAT HAVE < 0.8 METHYLATION LEVEL AND ITS CORRESPONDING Ks
METHY_K_INTER_CSV=$DATA_DIR/METHY_K_INTER.tsv
METHY_K_GREATER_CSV=$DATA_DIR/METHY_K_GREATER.tsv
gsed -e "s/^chr//" $METHY_K> $METHY_K_CSV
gsed -e "s/^chr//" $METHY_K_LESS> $METHY_K_LESS_CSV
gsed -e "s/^chr//" $METHY_K_INTER> $METHY_K_INTER_CSV
gsed -e "s/^chr//" $METHY_K_GREATER> $METHY_K_GREATER_CSV

# Copy features bed files
OLD_DIR=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Genomic_Features
cd $OLD_DIR

for dir_name in *;
do
	echo $dir_name
	dist=$FEATURE_ANALYSIS/$dir_name
	mkdir -p $dist
	cp $dir_name/*.bed $dist
done


FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis
DATA_DIR=$FEATURE_ANALYSIS/DATA
METHY_K_LESS=$DATA_DIR/METHY_K_LESS.bed #FILES THAT HAVE < 0.8 METHYLATION LEVEL AND ITS CORRESPONDING Ks
METHY_K_INTER=$DATA_DIR/METHY_K_INTER.bed
METHY_K_GREATER=$DATA_DIR/METHY_K_GREATER.bed
# WORK_DIR=$FEATURE_ANALYSIS/ADDED
# TMP_DIR=$FEATURE_ANALYSIS/TMP
# declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')
# Process for ADDED directories
# cd $WORK_DIR
# for f in *.bed;
# do
# 	filename=${f%%.*}
#   	echo $filename
#   	f1=$filename."tmp"
#   	awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f | gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$f1
# 	for chr in "${CHRS[@]}"; 
# 	do
# 	    bedextract $chr $f1 > $TMP_DIR/$chr.bed; 
# 	done
# 	cat $TMP_DIR/*.bed |gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$f
# 	rm $TMP_DIR/*.bed $f1
# done


declare -a MARKERS=('ADDED') #'General' 'Genomic_Regions' 'Histone_Modification' 'Histone_Modification_Enzyme' 'TFBS')
#declare -a MARKERS=('ChrmoHMM')
ChromoHMM=0
for region in "${MARKERS[@]}"; #
do
    f1=$METHY_K_GREATER
	WORK_DIR=$FEATURE_ANALYSIS/$region
	OUT_DIR=$WORK_DIR/METHY_K_GREATER
	mkdir -p $OUT_DIR
	cd $WORK_DIR

	#Do intersection
	for f in *.bed;
	do
		filename=${f%%.*}
	  	echo $filename
		f2=$OUT_DIR/$filename"_i.bed"
		f3=$OUT_DIR/$filename"_v.bed"
		f4=$OUT_DIR/$filename".csv"
		if [ $ChromoHMM -eq 1 ]
		then
			bedtools intersect -a $f1 -b $f -wa -wb -sorted | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}'| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6,7 -o max >$f2
		else
			bedtools intersect -a $f1 -b $f -wa -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6 -o max | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"1}'>$f2
		fi
		bedtools intersect -a $f1 -b $f -v -sorted|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"0}'>$f3
		cat $f2 $f3 |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' | gsed -e "s/^chr//"| gsort -k 1,1 -k2,2n --parallel=8  -S 50%> $f4
		rm $f2 $f3
	done
done






#CG content analysis
OLD_FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/FEATURE_ANALYSIS
hg19=$OLD_FEATURE_ANALYSIS/hg19.fa
CG_Analysis=$FEATURE_ANALYSIS/CG_Analysis
RADIUS=25
mkdir -p $CG_Analysis
declare -a DATA_TYPES=('METHY_K_LESS' 'METHY_K_INTER' 'METHY_K_GREATER')

for DT in "${DATA_TYPES[@]}"; 
do
	f1=$DATA_DIR/$DT".bed"
	f2=$CG_Analysis/$DT"_flank_$RADIUS.bed"
	f3=$CG_Analysis/$DT"_CG.bed"
	echo $f2
	awk -v radius=$RADIUS 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2-radius"\t"$3+radius-1"\t"$4"\t"$5"\t"$6}' $f1 >$f2

	echo $f3
	bedtools nuc -fi $hg19 -bed $f2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {if (NR!=1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8} }'|gsed -e "s/^chr//">$f3
	rm $f2
done

# K Cluster Analysis
FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis
DATA_DIR=$FEATURE_ANALYSIS/DATA
METHY_K=$DATA_DIR/METHY_AND_K.bed
CLUSTER_DIR=$FEATURE_ANALYSIS/K_Cluster_Analysis
declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

#Prepare data for cluster
cd $CLUSTER_DIR
INPUT_DATA=$CLUSTER_DIR/INPUT_DATA
mkdir -p $INPUT_DATA
for chr in "${CHRS[@]}";
do
    bedextract $chr $METHY_K | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|gsed -e "s/^chr//">$INPUT_DATA/$chr.bed; 
done

CLUSTER_DATA=$CLUSTER_DIR/CLUSTER_DATA
MERGED_CLUSTER_DATA=$CLUSTER_DIR/Cluster_merged.bed

cat $CLUSTER_DATA/*.bed |gsort -k 1,1 -k2,2n --parallel=8  -S 50%| sed 's/nan/0.0/'>$MERGED_CLUSTER_DATA

fHMM=$FEATURE_ANALYSIS/ChrmoHMM/HmmH1hescHMM.bed
ChromoHMM_LABELED_CLUSTER=$CLUSTER_DIR/Cluster_ChromoHMM.bed
bedtools intersect -a $MERGED_CLUSTER_DATA -b $fHMM -f 0.5 -wa -wb -sorted|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$14}' | gsort -k 1,1 -k2,2n --parallel=8  -S 50% |bedtools merge -c 4,5,6,7,8,9,10,11 -o max>$ChromoHMM_LABELED_CLUSTER
rm -f $MERGED_CLUSTER_DATA

#CG content_ANNOTATION
OLD_FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/FEATURE_ANALYSIS
hg19=$OLD_FEATURE_ANALYSIS/hg19.fa
Cluster_CG=$CLUSTER_DIR/Cluster_CG.bed
bedtools nuc -fi $hg19 -bed $ChromoHMM_LABELED_CLUSTER| gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {if (NR!=1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13} }' >$Cluster_CG


#Intersect cluster data with other annotations
f1=$ChromoHMM_LABELED_CLUSTER
OUT_DIR=$CLUSTER_DIR/ANNOTATION_DATA
mkdir -p $OUT_DIR

declare -a MARKERS=('ADDED' 'General' 'Genomic_Regions' 'Histone_Modification' 'Histone_Modification_Enzyme' 'TFBS' ) # 'DMR'
DMR=0
for region in "${MARKERS[@]}"; #
do
	WORK_DIR=$FEATURE_ANALYSIS/$region
	cd $WORK_DIR
	OUT_DIR_HERE=$OUT_DIR/$region
	mkdir -p $OUT_DIR_HERE
	#Do intersection
	for f in *.bed;
	do
		filename=${f%%.*}
	  	echo $filename
		f2=$OUT_DIR_HERE/$filename"_i.bed"
		f3=$OUT_DIR_HERE/$filename"_v.bed"
		f4=$OUT_DIR_HERE/$filename".tsv"
		if [ $DMR -eq 1 ]
		then
			bedtools intersect -a $f1 -b $f -wa -wb -sorted | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$15}'| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4 -o max>$f2 #
		else
			bedtools intersect -a $f1 -b $f -wa -sorted | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"1}' | gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4 -o max>$f2 #
		fi
		bedtools intersect -a $f1 -b $f -v -sorted| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"0}'>$f3
		cat $f2 $f3 |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4}' | gsort -k 1,1 -k2,2n --parallel=8  -S 50%| gsed -e "s/^chr//"> $f4
		rm $f2 $f3
	done
done

awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f1| gsed -e "s/^chr//"> $f1".tmp"
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}' $f4> $f4".tmp"
diff $f1".tmp" $f4".tmp"
rm $f1".tmp" $f4".tmp"

declare -a CHRS=('chr1')
for chr in "${CHRS[@]}"; 
do
    bedextract $chr $f1 | gsed -e "s/^chr//" > $CLUSTER_DIR/Cluster_ChromoHMM_$chr.tsv; 
done

gsed -e "s/^chr//" $f1 >$CLUSTER_DIR/Cluster_ChromoHMM.tsv

f1=$CLUSTER_DIR/Cluster_ChromoHMM.bed
f2=$CLUSTER_DIR/Methylated_And_Low_K_Cluster.bed #_Full
f3=$CLUSTER_DIR/Cluster_ChromoHMM_cols_reduced.bed
f4=$CLUSTER_DIR/Methylated_And_High_K_Cluster.bed
#Extract Methylated (>0.7) but low K's clusters
awk  -F "\t" '{ if($5 < 1 && $7 > 0.5) { print $1"\t"$2"\t"$3} }' $f1 > $f2
awk  -F "\t" '{ if($5 > 7 && $7 > 0.5) { print $1"\t"$2"\t"$3} }' $f1 > $f4

awk  -F "\t" '{ print $1"\t"$2"\t"$3}' $f1 > $f3

f1=$FEATURE_ANALYSIS/ADDED/METHY_K_GREATER/NANOG.csv
f2=$FEATURE_ANALYSIS/ADDED/METHY_K_GREATER/Methylated_And_Low_K_NANOG.bed
f3=$FEATURE_ANALYSIS/ADDED/METHY_K_GREATER/K_and_NANOG.bed
awk  -F "\t" '{ if($4 < 2 && $3 > 0.5 && $6 == 1) { print "chr"$1"\t"$2"\t"$2+1} }' $f1 > $f2
awk  -F "\t" '{ if($6 == 1) { print "chr"$1"\t"$2"\t"$2+1} }' $f1 > $f3


f1=$FEATURE_ANALYSIS/DATA/METHY_AND_K.bed
f2=$FEATURE_ANALYSIS/DATA/Methylated_And_Low_K.bed
#Extract Methylated (>0.7) but low K's clusters
awk  -F "\t" '{ if($5 < 1 && $4 > 0.5) { print $1"\t"$2"\t"$3} }' $f1 > $f2
f3=$FEATURE_ANALYSIS/DATA/METHY_AND_K_cols_reduced.bed
awk  -F "\t" '{ print $1"\t"$2"\t"$3}' $f1 > $f3


f1=$CLUSTER_DIR/Cluster_ChromoHMM.bed
f2=$CLUSTER_DIR/Methylated_And_Low_K_Cluster.bed
#Extract Methylated (>0.7) but low K's clusters
awk  -F "\t" '{ if($5 < 1 && $7 > 0.5) { print $1"\t"$2"\t"$3} }' $f1 > $f2

f1=$CLUSTER_DIR/Cluster_ChromoHMM.tsv
f2=$CLUSTER_DIR/Methylated_And_High_K_Coef_Var_Cluster.tsv
#Extract Methylated (>0.7) but low K's clusters
awk  -F "\t" '{ if($6 > 0.7 && $7 > 0.5) { print} }' $f1 > $f2
f1=$CLUSTER_DIR/Cluster_ChromoHMM.bed
f3=$CLUSTER_DIR/Methylated_And_High_K_Coef_Var_Cluster.bed
awk  -F "\t" '{ if($6 > 0.7 && $7 > 0.5) { print $1"\t"$2"\t"$3} }' $f1 > $f3

FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis
DATA_DIR=$FEATURE_ANALYSIS/DATA
#site specific merthylated and low k
f1=$DATA_DIR/METHY_AND_K.bed
f2=$DATA_DIR/Inter_Methylated_And_Low_K.bed
f3=$DATA_DIR/Unmethylated_And_Low_K.bed
# f4=$DATA_DIR/Methylated_And_Low_K.bed
# f5=$DATA_DIR/Methylated_And_Low_K_sorted.bed
#Extract Methylated (>0.7) but low K's clusters
# awk  -F "\t" '{ if($5 < 0.4 && $4 < 0.3) { print $1"\t"$2"\t"$3} }' $f1  > $f3
awk  -F "\t" '{ if($5 < 0.5 && $4 > 0.5 && $4 < 0.8) { print $1"\t"$2"\t"$3} }' $f1  > $f2 #| gsed -e "s/^chr//"
gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f4 > $f5

f3=$DATA_DIR/Methylated_And_High_K.tsv
awk  -F "\t" '{ if($5 > 7 && $4 > 0.9) { print $1"\t"$2"\t"$3} }' $f1 | gsed -e "s/^chr//"> $f3 #$1"\t"$2"\t"$3

