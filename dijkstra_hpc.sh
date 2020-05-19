#!/bin/bash
#
#SBATCH --job-name=matmult
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --nodes=1
#SBATCH --output=matmult-%j.out #SBATCH --time=10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hchen7@scu.edu #
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=matmult-$current_time.log

/WAVE/projects/COEN-319-Sp20/hchen7/hw3/dijkstra-parallel/dijkstra 1000 1000 1000  >> /WAVE/projects/COEN-319-Sp20/hchen7/hw2/logs/$log_file_name

for size in 100 500 1000 5000 10000
do
  for source in
  /WAVE/projects/COEN-319-Sp20/hchen7/hw3/dijkstra-parallel/dijkstra "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}.graph" 0 "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}_serial.out" >> "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /logs/$log_file_name"
done