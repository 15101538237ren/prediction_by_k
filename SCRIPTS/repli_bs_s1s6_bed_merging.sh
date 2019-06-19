#!/bin/bash
#$ -N REPLI_BS_S1S6_BED_MERGING
#$ -q rxn
#$ -m beas

module load bedtools/2.25.0
module load bedops/2.4.14

cd '/pub/hongleir/DATA_FOR_ANALYSIS/repli_bs'

cat *.bed |sort -k 1,1 -k2,2n | bedtools merge -c 4,5 -o sum >'S1_S6.bed'