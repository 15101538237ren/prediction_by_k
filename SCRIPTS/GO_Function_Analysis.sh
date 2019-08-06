#!/bin/bash

# Generate methy and k percentile partitioned K data
FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis
DATA_DIR=$FEATURE_ANALYSIS/DATA
METHY_AND_K=$DATA_DIR/METHY_AND_K.bed

OUT_DIR=$DATA_DIR/PERCENTILED_K
mkdir -p $OUT_DIR
FILTER_K=$DATA_DIR/FILTERED_PERCENTILED_K
declare -a methy_percentiles=(0 0.5 0.8 1.0) # Methylation percentiles
declare -a methy_labels=("LM" "IM" "HM") # Methylation interval labels
declare -a k_percentiles=(0.00 0.54 1.07 1.66 2.40 3.24 3.98 4.68 5.50 6.31 7.41 10.00) # K percentiles
declare -a k_labels=("K_0_5" "K_5_15" "K_15_25" "K_25_35" "K_35_45" "K_45_55" "K_55_65" "K_65_75" "K_75_85" "K_85_95" "K_95_100") # Methylation interval labels

for (( i=0; i<${#methy_percentiles[@]} - 1; i++ ));
do
	for (( j=0; j<${#k_percentiles[@]} - 1; j++ ));
	do
		file_name="${methy_labels[$i]}_${k_labels[$j]}.bed"
		echo $file_name
		f=$OUT_DIR/$file_name
		#echo "methy: ${methy_percentiles[$i]} - ${methy_percentiles[(($i + 1))]}, k: ${k_percentiles[$j]} - ${k_percentiles[(($j + 1))]}"
		awk -v lm=${methy_percentiles[$i]} -v hm=${methy_percentiles[(($i + 1))]} -v lk=${k_percentiles[$j]} -v hk=${k_percentiles[(($j + 1))]} -F "\t" '{ if( $4 > lm && $4 < hm && $5 > lk && $5 < hk) { print $1"\t"$2"\t"$3 } }' $METHY_AND_K  > $f
		# Use python to filter the file
		f2=$FILTER_K/$file_name
		gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f2> $f2".tmp"
		rm $f2
		mv $f2".tmp" $f2
	done
	file_name="${methy_labels[$i]}.bed"
	echo $file_name
	f=$OUT_DIR/$file_name
	# awk -v lm=${methy_percentiles[$i]} -v hm=${methy_percentiles[(($i + 1))]} -F "\t" '{ if( $4 > lm && $4 < hm ) { print $1"\t"$2"\t"$3 } }' $METHY_AND_K  > $f
	# Use python to filter the file
	f2=$FILTER_K/$file_name
	gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f2> $f2".tmp"
	rm $f2
	mv $f2".tmp" $f2
done

# Distance to nearest TSS
FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Filtered_Feature_Analysis
DATA_DIR=$FEATURE_ANALYSIS/DATA
GENE_FP=$DATA_DIR/GENES.bed

awk -F "\t" '{ if( $4 == "+") { print $1"\t"$2"\t"$2"\t+\t"$5} else { print $1"\t"$3"\t"$3"\t-\t"$5} }' $GENE_FP| gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$GENE_FP".tmp"
mv $GENE_FP $GENE_FP".bk"
mv $GENE_FP".tmp" $GENE_FP
DIST_DIR=$DATA_DIR/K_DISTANCE_TO_TSS
mkdir -p $DIST_DIR
cd $OUT_DIR
fHMM=$DATA_DIR/HmmH1hescHMM.bed
for f in *.bed;
do
	filename=${f%%.*}
	echo $filename
	f2=$DIST_DIR/$filename".tsv"
	#bedtools closest -a $f -b $GENE_FP -D b  | sort -k 1,1 -k2,2n --parallel=8  -S 50%| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$9}' | gsed -e "s/^chr//" > $f2
	f3=$DIST_DIR/$filename".bed"
	#awk 'BEGIN {FS="\t"; OFS=","} {print "chr"$1"\t"$2"\t"$2+1"\t"$3}' $f2>$f3
	bedtools intersect -a $f3 -b $fHMM -f 0.5 -wa -wb -sorted |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$f3".tmp"
	rm $f3
	mv $f3".tmp" $f3
done 

