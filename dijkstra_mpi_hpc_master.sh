for node_count in `seq 1 3`; do   
    sbatch ./dijkstra_mpi_hpc.sh $node_count dijkstra_mpi
done