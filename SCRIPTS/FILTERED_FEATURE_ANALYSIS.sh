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
METHY_K_GREATER=$DATA_DIR/METHY_K_GREATER.bed
# Data ordering checked
paste $METHY_FILE $K_FILE | cut -f 1,2,3,4,8,9 >$METHY_K
# Process Methylation > 0.8
awk  -F "\t" '{ if($4 < 0.8) { print } }' $METHY_K > $METHY_K_LESS
awk  -F "\t" '{ if($4 > 0.8) { print } }' $METHY_K > $METHY_K_GREATER

METHY_K_LESS_CSV=$DATA_DIR/METHY_K_LESS.tsv #FILES THAT HAVE < 0.8 METHYLATION LEVEL AND ITS CORRESPONDING Ks
METHY_K_GREATER_CSV=$DATA_DIR/METHY_K_GREATER.tsv
gsed -e "s/^chr//" $METHY_K_LESS> $METHY_K_LESS_CSV
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
METHY_K_GREATER=$DATA_DIR/METHY_K_GREATER.bed

# declare -a MARKERS=('General' 'Genomic_Regions' 'Histone_Modification' 'Histone_Modification_Enzyme' 'TFBS')
declare -a MARKERS=
ChromoHMM=1
for region in "${MARKERS[@]}"; #
do
    f1=$METHY_K_LESS
	WORK_DIR=$FEATURE_ANALYSIS/$region
	OUT_DIR=$WORK_DIR/METHY_K_LESS
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



