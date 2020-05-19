current_time=$(date "+%Y.%m.%d-%H.%M.%S")
log_file_name=dijkstra_serial_local-$current_time.log
for size in 100 500 1000 5000 10000
do
  "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /pr2/dijkstra" "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}.graph" 0 "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /test_data/${size}_serial.out" >> "/Users/kexinchen/SCU Google Drive/COEN 319/Homework/Homework 3 /logs/$log_file_name"
done