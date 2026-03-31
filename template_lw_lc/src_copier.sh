#!/bin/bash
set -euo pipefail

# ---- CONFIGURATION ----
read -r N < ../n_samples.txt

# Define your list dynamically
i_list=($(seq 1 "$N"))

# ---- MAIN LOOP ----
for i in "${i_list[@]}"; do
    folder="sample_${i}"
    
    # Copy files/folders into the new folder
    cp -R src/* "$folder"/
done

