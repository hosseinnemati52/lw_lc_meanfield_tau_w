#!/bin/bash
set -euo pipefail

# ---- CONFIGURATION ----
# Define your list (can be any values)
read -r N < ../n_samples.txt

# Define your list dynamically
i_list=($(seq 1 "$N"))

# ---- MAIN LOOP ----
for i in "${i_list[@]}"; do
    folder="sample_${i}"
    printf "Creating folder: %s\n" "$folder"
    mkdir -p -- "$folder"
done

