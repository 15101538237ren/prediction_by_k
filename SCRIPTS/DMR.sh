WORK_DIR=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/WORK_DIR
cd $WORK_DIR
hesc=GSM1112841.bed
ec=GSM1112820.bed
merged_hesc_ec=MERGED.bed
dmr=EC_DMR.bed

declare -a CHRS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

for chr in "${CHRS[@]}"; 
do
    bedextract $chr $dmr>$chr.bed; 
done

cat chr*.bed | gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$dmr".filtered"
rm -f $dmr chr*.bed
mv $dmr".filtered" $dmr


bedtools intersect -a $hesc -b $ec -wa -wb -sorted | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$8}' | gsort -k 1,1 -k2,2n --parallel=8  -S 50% >$merged_hesc_ec

bedtools intersect -a $merged_hesc_ec -b $dmr -wa -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%  |gsed 's/\t/,/g'| gsed -e "s/^chr//"> mdmr_intersected.bed 
bedtools intersect -a $merged_hesc_ec -b $dmr -v -sorted | gsort -k 1,1 -k2,2n --parallel=8  -S 50%  |gsed 's/\t/,/g'| gsed -e "s/^chr//"> mdmr_unintersected.bed 

cat