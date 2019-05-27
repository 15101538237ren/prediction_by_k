#!/bin/bash
#$ -N PEAK_CALLING
#$ -q pub*,ionode,rxn
#$ -m beas

function usage {
    echo "usage: PEAK_CALLING.sh [-c control] [-f infile] [-i indir] [-o outdir]"
    echo "  -h          display help"
    echo "  -s The sicer path"
    echo "  -c With control or not, e.g 1 with control, 0 without control"
    echo "  -d sh executing dir"
    echo "  -a bed file, specify the path of bed file"
    echo "  -b control file, specify the path of bed file"
    echo "  -i input dir   specify the input data directory"
    echo "  -p broad peak: 1 otherwise 0"
    echo "  -o output dir   specify the output data directory"
    exit 1
}

if [  $# -lt 6 ]
then
        usage
        exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
        -s | --sicer_path )        shift
                                sicer_path=$1
                                ;;
        -o | --outdir )         shift
                                outdir=$1
                                ;;
        -i | --indir )          shift
                                indir=$1
                                ;;
        -c | --control )        shift
                                control=$1
                                ;;
        -d | --sh_exe_dir )     shift
                                sh_exe_dir=$1
                                ;;
        -a | --bedfile )        shift
                                bed_file=$1
                                ;;
        -b | --controlfile )    shift
                                control_file=$1
                                ;;
        -p | --broadpeak )      shift
                                broad_peak=$1
                                ;;
        -f | --inputfile )      shift
                                inputfile=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

SICER=$sicer_path

# Default parameters for peak calling
species="hg19"
redundancy_threshold=1
window_size=200
fragment_size=150
effective_genome_fraction=0.74
FDR=0.01
gap_size=200
e_value=100
cmd="SICER.sh"

mkdir -p $outdir

if [ $control -eq 1 ]; then
	# With control
	cmd='SICER.sh'
	value=$FDR
else
	# Without control data
	cmd='SICER-rb.sh'
	value=$e_value
fi

if [ $broad_peak -eq 1 ]; then
	#broad peak
	gap_size=$((3*$window_size))
else
    gap_size=$window_size
fi

cd "$sh_exe_dir"

sh "$SICER/$cmd" $indir $bed_file $control_file $outdir $species $redundancy_threshold $window_size $fragment_size $effective_genome_fraction $gap_size $value


