#!/bin/bash
#
#SBATCH --job-name=dijkstra_mpi
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1
#SBATCH --output=dijkstra_mpi-%j.out 
#SBATCH --time=23:59:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hchen7@scu.edu #

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=dijkstra_mpi-$current_time.log
file_path=/WAVE/projects/COEN-319-Sp20/hchen7/hw3 
for size in 100 500 1000 5000 10000; do
    for thread_num in 1 2 4 8 12 14 16 20 24 28; do
        for k in `seq 1 20`; do   
            source=$RANDOM
            let "source %= $size"
            mpirun -n $thread_num $file_path/dijkstra-parallel/dijkstra_mpi "$file_path/test_data/${size}.graph" $source "$file_path/test_data/${size}_serial.out" "$file_path/logs/$log_file_name"
        done 
    done
done
