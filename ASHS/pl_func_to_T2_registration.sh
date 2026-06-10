#!/bin/bash

# --- 1. specify base path ---
BASE_DIR="/mnt/volume/tacmem_workspace/derivatives"

# --- 2. specify subject and analysis lists ---
subjects=("sub-03" "sub-05" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-20" "sub-24" "sub-27" "sub-28" "sub-29" "sub-30" "sub-32" "sub-33" "sub-34" "sub-35" "sub-36" "sub-37" "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-47")

analyses=("TNT" "recall")

# --- 3. func to T2 registration ---
for subject_id in "${subjects[@]}"; do
    echo "================================================"
    echo "Processing Registration for: ${subject_id}"

    # define paths constant for the subject
    # RF_REF: use the ASHS segmentation as the reference for the high-res T2* grid
    RF_REF="${BASE_DIR}/ASHS_seg/hippocampus_seg/${subject_id}/final/${subject_id}_left_lfseg_corr_usegray.nii.gz"
    REG_MAT="${BASE_DIR}/ASHS_seg/greedy_reg/${subject_id}/rigid.mat"

    # pre-check: ensure the transformation matrix and reference exist
    if [[ ! -f "$REG_MAT" ]] || [[ ! -f "$RF_REF" ]]; then
        echo "SKIPPING ${subject_id}: Missing rigid.mat or reference T2* segmentation."
        continue
    fi

    # loop through the T-map result folders (TNT, recall, etc.)
    for analysis_type in "${analyses[@]}"; do
        SUBJ_DIR="${BASE_DIR}/results/individual_analyses/${analysis_type}/${subject_id}"

        # check if this specific analysis folder exists for the subject
        if [[ ! -d "$SUBJ_DIR" ]]; then
            echo "  [!] Folder not found: $analysis_type for $subject_id"
            continue
        fi

        echo "  --- Analysis: $analysis_type ---"

        # find and process every T-map inside the analysis folder
        for MOV_TMAP in "${SUBJ_DIR}"/spmT_*.nii; do
            
            # ensure the wildcard actually found files
            [[ -e "$MOV_TMAP" ]] || { echo "    No spmT_*.nii files found in $analysis_type"; continue; }
            
            # prepare output filename
            T_BASE=$(basename "$MOV_TMAP" .nii)
            OUT_TMAP="${SUBJ_DIR}/${T_BASE}_in_T2star_space.nii.gz"

            echo "    -> Reslicing $T_BASE..."
            
            # apply the rigid transformation using Greedy
            greedy -d 3 \
                -rf "$RF_REF" \
                -ri LINEAR \
                -r "$REG_MAT" \
                -rm "$MOV_TMAP" \
                "$OUT_TMAP"
        done
    done
done

echo "================================================"
echo "All T-map registration tasks complete!"
