#!/bin/bash

# --- 1. specify base paths ---
DATA_BASE="/mnt/volume/tacmem_workspace/derivatives/data"
OUTPUT_BASE="/mnt/volume/tacmem_workspace/derivatives/ASHS_seg/greedy_reg"

# --- 2. greedy registration ---
for sub_dir in ${DATA_BASE}/sub-*; do
    
    #extract the subject ID from the full path
    SUB_ID=$(basename "$sub_dir")
    
    #define file paths for this specific subject
    T1="${sub_dir}/anat/${SUB_ID}_T1w.nii"
    T2="${sub_dir}/anat/${SUB_ID}_T2w.nii"
    SUB_OUT_DIR="${OUTPUT_BASE}/${SUB_ID}"
    OUT_MAT="${SUB_OUT_DIR}/rigid.mat"

    #check if both T1 and T2 images exist before running
    if [[ -f "$T1" && -f "$T2" ]]; then
        echo "Processing ${SUB_ID}..."

        #create the output directory if it doesn't exist
        mkdir -p "$SUB_OUT_DIR"

        #run the greedy command
        greedy -d 3 -a \
          -i "$T2" "$T1" \
          -ia-image-centers \
          -dof 6 \
          -o "$OUT_MAT" \
          -n 100x50x10 \
          -m NMI
          
        echo "finished ${SUB_ID}. Matrix saved to ${OUT_MAT}"
        echo "-----------------------------------------------"
    else
        echo "skipping ${SUB_ID}: Missing T1 or T2 files."
    fi
done
