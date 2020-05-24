#!/bin/bash
job_name=dijkstra_mpi
#
#SBATCH --job-name=$job_name
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=64G
#SBATCH --nodes=$0
#SBATCH --output=$job_name-%j.out 
#SBATCH --time=23:59:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hchen7@scu.edu #

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=job_name-$current_time.log

for size in 100 500 1000 5000 10000; do
    for thread_num in 1 2 4 8 12 14 16 20 24 28; do
        for k in `seq 1 20`; do   
            source=$RANDOM
            let "source %= $size"
            printf $0, $thread_num, $size , $k, $source ,
            mpirun -n $thread_num /WAVE/projects/COEN-319-Sp20/hchen7/hw3/dijkstra-parallel/dijkstra_mpi "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}.graph" $source "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}_serial.out" >> "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /logs/$log_file_name"
        done 
    done
done
