current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=dijkstra_mpi_local-$current_time.log
file_path="/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 "
export PMIX_MCA_gds=hash
for size in 100 500 1000 5000 10000; do
# for size in 100
  # echo -n $size , >> "$file_path/logs/$log_file_name"
  mpirun -n 4 "$file_path/pr2/dijkstra_mpi" "$file_path/test_data/${size}.graph" 0 "$file_path/test_data/${size}_mpi.out" "$file_path/logs/$log_file_name"
done