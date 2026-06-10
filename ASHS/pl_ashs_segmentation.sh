#!/bin/bash

# --- 1. specify system path ---
ASHS_ROOT="/home/pengliu/Software/ASHS/ashs-fastashs"
ATLAS="/home/pengliu/Software/ASHS/atlas"

# --- 2. specify base path ---
DATA_BASE="/mnt/volume/tacmem_workspace/derivatives/data"
OUTPUT_BASE="/mnt/volume/tacmem_workspace/derivatives/ASHS_seg/hippocampus_seg"
GREEDY_BASE="/mnt/volume/tacmem_workspace/derivatives/ASHS_seg/greedy_reg"

# --- 3. specify subject ---
subjects=("sub-03" "sub-05" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-20" "sub-24" "sub-27" "sub-28" "sub-29" "sub-30" "sub-32" "sub-33" "sub-34" "sub-35" "sub-36" "sub-37" "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-47")

# --- 4. ASHS segmentation ---
for subject_id in "${subjects[@]}"; do
    
    echo "processing ${subject_id}..."

    # specify paths for this subject
    SUBJECT_OUT_DIR="${OUTPUT_BASE}/${subject_id}"
    INPUT_T1="${DATA_BASE}/${subject_id}/anat/${subject_id}_T1w.nii"
    INPUT_T2="${DATA_BASE}/${subject_id}/anat/${subject_id}_T2w.nii"
    RIGID_MATRIX="${GREEDY_BASE}/${subject_id}/rigid.mat"

    # check if the greedy matrix exists before running
    if [ ! -f "$RIGID_MATRIX" ]; then
        echo "ERROR: Matrix file for ${subject_id} not found at ${RIGID_MATRIX}. Skipping."
        continue
    fi

    # create the output folder (now uses the correct path)
    mkdir -p "$SUBJECT_OUT_DIR"

    # ASHS segmentation
    ${ASHS_ROOT}/bin/ashs_main.sh -P -I ${subject_id} \
        -a ${ATLAS}/ashs_atlas_abc \
        -g "${INPUT_T1}" \
        -f "${INPUT_T2}" \
        -w "${SUBJECT_OUT_DIR}" \
        -N -M -m "${RIGID_MATRIX}"

    echo "finished ${subject_id}"
    echo "------------------------------------------------"
    
done
