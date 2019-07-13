KRATE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS/K_RATES
FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/FEATURE_ANALYSIS
ORIGIN_DATA=$FEATURE_ANALYSIS/ORIGIN_DATA
CLASSIFIED_DATA=$FEATURE_ANALYSIS/CLASSIFIED_DATA
ChrmoHMM_INTERSECTED=$FEATURE_ANALYSIS/ChrmoHMM_INTERSECTED
fHMM=$FEATURE_ANALYSIS/HmmH1hescHMM.bed

declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')
declare -a CLASSES=('slowest' 'slow' 'middle' 'fast' 'fastest')

prefix='chr'

# Prepare the bed data and sort them
for rep in 1 41;
do
	f1_dir=$ORIGIN_DATA/Rep$rep
	mkdir -p $f1_dir
	f1=$KRATE_DIR/$rep.bed
	f1_dest=$f1_dir/tile_2.bed
	cp $f1 $f1_dest
	echo $f1_dest
	for TileSize in 100 200 500 1000 2000;
	do
		f1=$KRATE_DIR/$rep/$TileSize/$rep.bed
		f1_dest=$f1_dir/tile_$TileSize.bed
		cp $f1 $f1_dest
		echo $f1_dest
	done

	cd $f1_dir
	for f in *.bed;
	do
		filename=${f%%.*}
		gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f>$f".sorted"
		echo $f
		rm -f $f
		mv $f".sorted" $f
		mkdir -p $filename

		for chr in "${CHRS[@]}"; 
		do
		    bedextract $chr $f >$filename/$chr.bed; 
		done 

		cat $filename/*.bed |gsort -k 1,1 -k2,2n --parallel=8  -S 50%| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4}'|gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$filename.csv
		rm -rf $filename
	done
done

#Intersect the data with ChromoHMM annotation
mkdir -p $ChrmoHMM_INTERSECTED
for rep in 1 41;
do
	DEST_DIR=$ChrmoHMM_INTERSECTED/Rep$rep
	mkdir -p $DEST_DIR
	for TileSize in 2 100 200 500 1000 2000;
	do
		for cls in "${CLASSES[@]}"; 
		do
			f1=$CLASSIFIED_DATA/Rep$rep/"tile_"$TileSize"_"$cls".bed"
			f2=$DEST_DIR/"tile_"$TileSize"_"$cls
			echo $f2
			bedtools intersect -a $f1 -b $fHMM -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$8}' >$f2".bed"
			awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}' $f2".bed"|gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $f2".csv"; 
		done 
	done
done

# fast and slow K in different histone modification context
declare -a HISTONE_GEOIDS=('GSM772800' 'GSM772752' 'GSM772756' 'GSM997249' 'GSM772750' 'GSM772751')
declare -a HISTONE_LABELS=('H3K4me1' 'H3K4me3' 'H3K9Me3' 'H3K27ac' 'H3K27Me3' 'H3K36me3')
# HPC server side
indir=/pub/hongleir/data/ChIP-Seq_Peaks
histone_esc=$indir/histone_esc
mkdir -p $histone_esc
for i in "${!HISTONE_GEOIDS[@]}"; 
do
	gid=${HISTONE_GEOIDS[$i]}
	f1=$indir/$gid"_island.bed"
	f2=$histone_esc/${HISTONE_LABELS[$i]}".bed"
	echo $f2
	cp $f1 $f2
done
awk '$4 > 0  {print ;}' $indir/GSM772751_continuous.bed >$indir/H3K36me3.bed

cd $indir
f1='histone_esc.tar.bz2'
tar -jcvf $f1  histone_esc

#local side: download data from server
indir=/pub/hongleir/data/ChIP-Seq_Peaks
f1='histone_esc.tar.bz2'
FEATURE_ANALYSIS=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/FEATURE_ANALYSIS
scp hongleir@hpc.oit.uci.edu:$indir/$f1 $FEATURE_ANALYSIS

#Intersect the data with Histone Modification annotation
HISTONE_INTERSECTED=$FEATURE_ANALYSIS/HISTONE_INTERSECTED
mkdir -p $HISTONE_INTERSECTED

for histone_type in H3K36me3; #"${HISTONE_LABELS[@]}"; 
do
	OUT_DIR=$HISTONE_INTERSECTED/$histone_type
	HISTONE_BED_FILE=$FEATURE_ANALYSIS/histone_esc/$histone_type".bed"
	gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $HISTONE_BED_FILE>$HISTONE_BED_FILE".tmp"
	rm -f $HISTONE_BED_FILE
	mv $HISTONE_BED_FILE".tmp" $HISTONE_BED_FILE
	
	mkdir -p $OUT_DIR
	for rep in 1 41;
	do
		DEST_DIR=$OUT_DIR/Rep$rep
		mkdir -p $DEST_DIR
		for TileSize in 100 200 500 1000 2000;
		do
			for cls in "${CLASSES[@]}"; 
			do
				f1=$CLASSIFIED_DATA/Rep$rep/"tile_"$TileSize"_"$cls".bed"
				f2=$DEST_DIR/"tile_"$TileSize"_"$cls
				echo $f2
				bedtools intersect -a $f1 -b $HISTONE_BED_FILE -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$8}' >$f2".bed"
				awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}' $f2".bed"|gsed 's/\t/,/g'| gsed -e "s/^$prefix//"> $f2".csv";
				rm -f $f2".bed"
			done 
		done
	done
done


