#!/bin/bash

# Define your lists here
#i_list=(1 2 3 4 5 6 7 8)          # list of values for i
#j_list=(3 4 5 6 7 8)          # list of values for j

source lw_lc_lists.sh

main_folder=$(pwd)
# Optional: uncomment to see what was loaded
# echo "lw_list: ${lw_list[@]}"
# echo "lc_list: ${lc_list[@]}"

read -r N < n_samples.txt
# Define your list dynamically
i_list=($(seq 1 "$N"))

# Create folders
for lw in "${lw_list[@]}"; do
    for lc in "${lc_list[@]}"; do
	for i in "${i_list[@]}"; do
		folder="lw_${lw}__lc_${lc}/sample_${i}"
		# Copy files/folders into the new folder
		cp -R new_file_to_add/* "$folder"/
        done
    done
done
sleep 5

# Create folders
for lw in "${lw_list[@]}"; do
    for lc in "${lc_list[@]}"; do
	for i in "${i_list[@]}"; do
		folder="lw_${lw}__lc_${lc}/sample_${i}"
		# Copy files/folders into the new folder
		cd "$folder"
		python3 cost_mat_calc.py
		cd "$main_folder"
        done
    done
done
sleep 5

