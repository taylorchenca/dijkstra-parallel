current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=dijkstra_omp_local-$current_time.log
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
for size in 100 500 1000 5000 10000
do
  "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /pr2/dijkstra_omp" "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}.graph" 0 "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}_omp.out" >> "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /logs/$log_file_name"
done