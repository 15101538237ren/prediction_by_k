#!/bin/bash
#$ -N DATA_DOWNLOAD_BED_CONVERSION
#$ -q pub64,ionode,rxn
#$ -m beas
function usage {
    echo "usage: DATA_DOWNLOAD.sh [-f format] [-i infile] [-o outdir]"
    echo "  -h          display help"
    echo "  -f format   specify the format of the files you are downloading from geo, without .gz ending"
    echo "  -i infile   specify the input file with geo ids line by line"
    echo "  -o outdir   specify the target download directory"
    exit 1
}
if [  $# -lt 3 ] 
then 
	usage
	exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
        -o | --outdir )         shift
                                data_dir=$1
                                ;;
        -i | --infile )    	shift
				inputfp="$(pwd)/$1"
                                ;;
        -f | --format )         shift
                                format=$1
                                ;;
	-h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

mkdir -p $data_dir
cd $data_dir
while read p; do
    SUBP=$(echo $p|sed 's/...$/nnn/')
    RECEVE_FILE_NAME="$p.$format.gz"
    wget -O "$RECEVE_FILE_NAME" "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/$SUBP/$p/suppl/*.$format.gz"
    gunzip $RECEVE_FILE_NAME
done <$inputfp
