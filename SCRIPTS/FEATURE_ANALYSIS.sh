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
for rep in 55; #1 41
do
	f1_dir=$ORIGIN_DATA/Rep$rep
	mkdir -p $f1_dir
	f1=$KRATE_DIR/$rep.bed
	# awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$3"\t"$2"\t"$4"\t"$5}' $f1> $f1".sorted"
	# rm $f1
	# mv $f1".sorted" $f1

	f1_dest=$f1_dir/tile_2.bed
	cp $f1 $f1_dest
	# echo $f1_dest
	# for TileSize in 100 200 500 1000 2000;
	# do
	# 	f1=$KRATE_DIR/$rep/$TileSize/$rep.bed
	# 	f1_dest=$f1_dir/tile_$TileSize.bed
	# 	cp $f1 $f1_dest
	# 	echo $f1_dest
	# done

	cd $f1_dir
	for f in *.bed;
	do
		filename=${f%%.*}
		gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f>$f".sorted"
		echo $filename
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
for rep in 55; # 1 41
do
	DEST_DIR=$ChrmoHMM_INTERSECTED/Rep$rep
	mkdir -p $DEST_DIR
	for TileSize in 2; # 100 200 500 1000 2000
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
GENOMIC_FEATURES=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Genomic_Features
HISTONE_INTERSECTED=$FEATURE_ANALYSIS/HISTONE_INTERSECTED
mkdir -p $HISTONE_INTERSECTED

for histone_type in "${HISTONE_LABELS[@]}"; #; 
do
	OUT_DIR=$HISTONE_INTERSECTED/$histone_type
	HISTONE_BED_FILE=$GENOMIC_FEATURES/Histone_Modification/$histone_type".bed"
	gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $HISTONE_BED_FILE>$HISTONE_BED_FILE".tmp"
	rm -f $HISTONE_BED_FILE
	mv $HISTONE_BED_FILE".tmp" $HISTONE_BED_FILE
	
	mkdir -p $OUT_DIR
	for rep in 55; #1 41
	do
		DEST_DIR=$OUT_DIR/Rep$rep
		mkdir -p $DEST_DIR
		for TileSize in 2; #100 200 500 1000 2000
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

#CpG density/CG content analysis
hg19=$FEATURE_ANALYSIS/hg19.fa
CG_Analysis=$FEATURE_ANALYSIS/CG_Analysis
rep=55
f1=$KRATE_DIR/$rep.bed
f1_dest=$CG_Analysis/$rep
cp $f1 $f1_dest".bed"
for lr in 25 50 75; #lr : local radius  
do
	f2=$f1_dest"_flank_"$lr".bed"
	echo $f2
	# awk -v radius=$lr 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2-radius"\t"$3+radius-1"\t"$4"\t"$5}' $f1_dest".bed">$f2
	f3=$f1_dest"_flank_"$lr"_GC.bed"
	echo $f3
	bedtools nuc -fi $hg19 -bed $f2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {if (NR!=1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7} }' >$f3 #|gsed 's/\t/,/g'| gsed -e "s/^chr//">$f3
done
# Promoter downloaded from UCSC table browser, using 2000bp upstream the TSS for each knowngene
Promoter=$CG_Analysis/Promoter.bed
gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $Promoter | bedtools merge -s -c 6 -o distinct>$Promoter".filtered"
rm -f $Promoter
mv $Promoter".filtered" $Promoter

declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')
# After Python Promoter calculation
for chr in "${CHRS[@]}"; 
do
    bedextract $chr $Promoter>$CG_Analysis/$chr.bed; 
done

cat chr*.bed | gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$Promoter".filtered"
rm -f $Promoter chr*.bed
mv $Promoter".filtered" $Promoter

f1_dest=$CG_Analysis/$rep
CGI=$CG_Analysis/CGI.bed

gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $CGI|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}'>$CGI".filtered"
rm -f $CGI
mv $CGI".filtered" $CGI

# Annotate by Promoter
for lr in 25; #lr : local radius   50 75
do
	f3=$f1_dest"_flank_"$lr"_GC"
	f4=$f3"_Promoter_annotated.bed"
	echo $f3
	bedtools intersect -a $f3".bed" -b $Promoter -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% |bedtools merge -c 4,5,6 -o max | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"1}'>$f3"_intersected.bed"
	bedtools intersect -a $f3".bed" -b $Promoter -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6 -o max | awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"0}'>$f3"_non_ovlp.bed"
	cat $f3"_intersected.bed" $f3"_non_ovlp.bed" | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > $f4
	rm $f3"_intersected.bed" $f3"_non_ovlp.bed"
	f3=$f1_dest"_flank_"$lr"_GC_Promoter_annotated"
	f4=$f1_dest"_flank_"$lr"_GC_Promoter_and_CGI_annotated.bed"
	f5=$f1_dest"_flank_"$lr"_Promoter_CGI_ChrHMM_annotated.csv"
	echo $f3
	bedtools intersect -a $f3".bed" -b $CGI -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% |bedtools merge -c 4,5,6,7 -o max | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"1}'>$f3"_i.bed"
	bedtools intersect -a $f3".bed" -b $CGI -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%|bedtools merge -c 4,5,6,7 -o max | awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"0}'>$f3"_v.bed"
	cat $f3"_i.bed" $f3"_v.bed" | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > $f4
	rm $f3"_i.bed" $f3"_v.bed"
	
	bedtools intersect -a $f4 -b $fHMM -f 0.5 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$12}'|gsed 's/\t/,/g'| gsed -e "s/^chr//"> $f5
done



# Methylation and K analysis
# Intersect K by merged_esimated_sites
f1=$K_RATES_DIR/55.bed
f4=/pub/hongleir/data/Methy-data
K_Methylation_Analysis=$FEATURE_ANALYSIS/K_Methylation_Analysis

scp $f1 hongleir@hpc.oit.uci.edu:$f4
module load bedtools/2.25.0
bedtools intersect -a $f4/GSM1112841.bed -b $f4/55.bed -wa| sort -k 1,1 -k2,2n|bedtools merge -c 5 -o max>GSM1112841_55_intersected.bed

scp hongleir@hpc.oit.uci.edu:$f4/GSM1112841_55_intersected.bed $K_Methylation_Analysis

ffn=$K_Methylation_Analysis/K_and_methy_intersected
f14=$K_Methylation_Analysis/GSM1112841_55_intersected.bed

gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f1>$f1".sorted"
rm -f $f1
mv $f1".sorted" $f1

gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f14>$f14".sorted"
rm -f $f14
mv $f14".sorted" $f14

bedtools intersect -a $f1 -b $f14 -wa -wb -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%>$ffn".bed"
awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$9}' $ffn".bed" |gsed 's/\t/,/g'| gsed -e "s/^$prefix//">$ffn".csv"

#CONSERVATION ANALYSIS
CONSERVATION=$FEATURE_ANALYSIS/CONSERVATION
# Download data from ucsc
#rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw $CONSERVATION

# Download DNaseI Hypersensitive Site Master List (125 cell types) from ENCODE/Analysis
# rsync -avz --progress rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseMasterSites/wgEncodeAwgDnaseMasterSites.bed.gz $OTHER_MARKERS

# Download Transcription Factor ChIP-seq Clusters (161 factors) from ENCODE with Factorbook Motifs
declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

GENOMIC_FEATURES=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Genomic_Features
# TFBS=$GENOMIC_FEATURES/TFBS
# rsync -avz --progress rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz $TFBS

TMP=$GENOMIC_FEATURES/TMP
mkdir -p $TMP

for f in $GENOMIC_FEATURES/*/*.bed;
do
	echo $f
	gsort -k 1,1 -k2,2n --parallel=8  -S 50% -i $f | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}'>$f".filtered"
	rm -f $f
	mv $f".filtered" $f
	# After Python Promoter calculation
	for chr in "${CHRS[@]}"; #
	do
	    bedextract $chr $f >$TMP/$chr.txt; 
	done
	cat $TMP/chr*.txt | gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3}'>$f".filtered"
	rm -f $f $TMP/chr*.txt
	mv $f".filtered" $f
done

KRATE_DIR=/Users/emmanueldollinger/MatlabProjects/K_Estimation/DATA/Repli_BS/K_RATES
GENOMIC_FEATURES=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/Genomic_Features

f1=$KRATE_DIR/55.bed
WORK_DIR=$GENOMIC_FEATURES/TFBS14
OUT_DIR=$WORK_DIR/K_intersected
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
	bedtools intersect -a $f1 -b $f -wa -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% |bedtools merge -c 4 -o max| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"1}'>$f2
	bedtools intersect -a $f1 -b $f -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"0}'>$f3
	cat $f2 $f3 |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}' | gsed -e "s/^chr//"| gsort -k 1,1 -k2,2n --parallel=8  -S 50%> $f4
	rm $f2 $f3
done
#gsed 's/\t/,/g'| gsed -e "s/^chr//">


WORK_DIR=$GENOMIC_FEATURES/ChrmoHMM
mkdir -p $WORK_DIR
cp $fHMM $WORK_DIR
OUT_DIR=$WORK_DIR/K_intersected
mkdir -p $OUT_DIR
#Do intersection
cd $WORK_DIR
for f in *.bed;
do
	filename=${f%%.*}
  	echo $filename
  	f2=$OUT_DIR/$filename"_i.bed"
	f3=$OUT_DIR/$filename"_v.bed"
	f4=$OUT_DIR/$filename".csv"
	bedtools intersect -a $f1 -b $f -wa -wb -sorted |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$9}' >$f2
	bedtools intersect -a $f1 -b $f -v -sorted|awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"0}'>$f3
	cat $f2 $f3 |awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$4"\t"$5}' | gsed -e "s/^chr//"| gsort -k 1,1 -k2,2n --parallel=8  -S 50%> $f4
	rm $f2 $f3
done










