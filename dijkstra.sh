#!/bin/bash
#
#SBATCH --job-name=dijkstra
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --nodes=1
#SBATCH --output=dijkstra-%j.out #SBATCH --time=10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hchen7@scu.edu #
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=dijkstra-$current_time.log

/WAVE/projects/COEN-319-Sp20/hchen7/hw2/dijkstra/dijkstra 1000 1000 1000  >> /WAVE/projects/COEN-319-Sp20/hchen7/hw2/logs/$log_file_name
/WAVE/projects/COEN-319-Sp20/hchen7/hw2/dijkstra/dijkstra 5000 5000 5000  >> /WAVE/projects/COEN-319-Sp20/hchen7/hw2/logs/$log_file_name
/WAVE/projects/COEN-319-Sp20/hchen7/hw2/dijkstra/dijkstra 1000 1000 2000  >> /WAVE/projects/COEN-319-Sp20/hchen7/hw2/logs/$log_file_name
/WAVE/projects/COEN-319-Sp20/hchen7/hw2/dijkstra/dijkstra 1000 2000 5000  >> /WAVE/projects/COEN-319-Sp20/hchen7/hw2/logs/$log_file_name
/WAVE/projects/COEN-319-Sp20/hchen7/hw2/dijkstra/dijkstra 9000 2500 3750  >> /WAVE/projects/COEN-319-Sp20/hchen7/hw2/logs/$log_file_name