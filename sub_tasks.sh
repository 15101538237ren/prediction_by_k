#!/bin/bash
function usage {
    echo "usage: sub_tasks.sh [-i input_file]"
    echo "  -h          display help"
    echo "  -s  the program is running on server or not"
    echo "  -i input file   specify the input file"
    exit 1
}

if [  $# -lt 2 ]
then
        usage
        exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
        -i | --input_file )     shift
                                input_file=$1
                                ;;
        -s | --server )         shift
                                server=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done
base_dir="/data/users/hongleir"
data_dir="/pub/hongleir/data"
sicer_path="/pub/hongleir/SICER_V1.1/SICER"
if [ $server -ne 1 ]; then
    # With control
    base_dir="/Users/Ren/PycharmProjects/prediction_by_k"
    data_dir="$base_dir/data"
    sicer_path="$base_dir/SICER_V1.1/SICER"
fi

peak_dir="$data_dir/ChIP-Seq_Peaks"
chip_seq_data_dir="$data_dir/ChIP-Seq-data"
while IFS=$'\t' read -r -a cols
	do
	dir_path="$base_dir/peak_calling/${cols[0]}"
	mkdir -p $dir_path
	cp "$base_dir/PEAK_CALLING.sh" $dir_path
	bed_file="${cols[1]}.bed"
	control_file="${cols[2]}.bed"
	control=${cols[4]}
	broad_peak=${cols[3]}
	if [ $control -eq 1 ]; then
		# With control
		sh "$dir_path/PEAK_CALLING.sh" -s "$sicer_path" -d "$dir_path" -c 1 -i "$chip_seq_data_dir" -o "$peak_dir" -a "$bed_file" -b "$control_file" -p $broad_peak
	else
		# Without control data
		sh "$dir_path/PEAK_CALLING.sh" -s "$sicer_path" -d "$dir_path" -c 0 -i "$chip_seq_data_dir" -o "$peak_dir" -a "$bed_file" -p $broad_peak
	fi
done < "$input_file"
