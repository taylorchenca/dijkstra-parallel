#!/bin/bash
job_name=dijkstra_serial
#
#SBATCH --job-name=$job_name
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --nodes=1
#SBATCH --output=$job_name-%j.out 
#SBATCH --time=23:59:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hchen7@scu.edu #

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=job_name-$current_time.log
file_path=/WAVE/projects/COEN-319-Sp20/hchen7/hw3/ 
for size in 100 500 1000 5000 10000; do
    for k in `seq 1 20`; do   
        source=$RANDOM
        let "source %= $size"
        printf $0, 1, $size , $k, $source ,
        $file_path/dijkstra-parallel/dijkstra "$file_path/test_data/${size}.graph" $source "$file_path/test_data/${size}_serial.out" >> "$file_path/logs/$log_file_name"
    done 
done
