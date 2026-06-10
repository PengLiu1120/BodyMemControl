#!/bin/bash
# Extract Max and Mean T-values from vlPFC and insula ROI masks
# Both masks and T-maps are in MNI space

# --- 1. specify base paths ---
BASE_DIR="/mnt/volume/tacmem_workspace/derivatives"
RESULT_DIR="${BASE_DIR}/results/vlPFC_insula_T_map_values"
mkdir -p "${RESULT_DIR}"

# --- 2. specify subject list ---
subjects=("sub-03" "sub-05" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-20" "sub-24" "sub-27" "sub-28" "sub-29" "sub-30" "sub-32" "sub-33" "sub-34" "sub-35" "sub-36" "sub-37" "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-47")

# --- 3. specify which contrast maps to extract ---
# key = label used in output, value = spmT filename
declare -A contrast_maps
contrast_maps=(["spmT_nothink_think"]="spmT_0002")

# --- 4. define ROI masks in MNI space ---
declare -A roi_masks
roi_masks=(
    ["vlPFC_right"]="${BASE_DIR}/results/group_analyses/TNT/nothink-think/ROI_r_vlPFC.nii"
    ["insula_right"]="${BASE_DIR}/results/group_analyses/TNT/nothink-think/ROI_r_insula.nii"
    ["insula_left"]="${BASE_DIR}/results/group_analyses/TNT/nothink-think/ROI_l_insula.nii"
)

export LC_NUMERIC="C"

# create temp dir, clean up on exit
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# --- 5. verify all masks exist ---
echo "Checking ROI masks..."
for roi_name in "${!roi_masks[@]}"; do
    mask_path="${roi_masks[$roi_name]}"
    if [ ! -f "$mask_path" ]; then
        echo "  [!] ERROR: Missing mask: ${roi_name} -> ${mask_path}"
        exit 1
    else
        echo "  [✓] Found: ${roi_name} -> ${mask_path}"
    fi
done
echo ""

# --- 6. write CSV header ---
OUTPUT_CSV="${RESULT_DIR}/vlPFC_insula_T_values.csv"
echo "Subject,Contrast_Map,ROI,Mean_T,Max_T" > "${OUTPUT_CSV}"

# --- 7. T map value extraction ---
for subject_id in "${subjects[@]}"; do
    echo "Processing ${subject_id}..."

    SUBJ_DIR="${BASE_DIR}/results/subject_analyses/TNT/${subject_id}"

    if [ ! -d "$SUBJ_DIR" ]; then
        echo "  [!] Missing directory for ${subject_id}, skipping."
        continue
    fi

    for contrast_label in "${!contrast_maps[@]}"; do
        contrast_file="${contrast_maps[$contrast_label]}"
        TMAP="${SUBJ_DIR}/${contrast_file}.nii"

        # also check for .nii.gz
        if [ ! -f "$TMAP" ]; then
            TMAP="${SUBJ_DIR}/${contrast_file}.nii.gz"
        fi

        if [ ! -f "$TMAP" ]; then
            echo "  [!] Missing T-map: ${contrast_file} for ${subject_id}, skipping."
            continue
        fi

        for roi_name in "${!roi_masks[@]}"; do
            mask_path="${roi_masks[$roi_name]}"

            # resample group mask to match subject T-map voxel dimensions
            # since both are in MNI space, only voxel size may differ
            RESAMPLED_MASK="${TMPDIR}/${subject_id}_${roi_name}_resampled.nii.gz"

            flirt \
                -in "$mask_path" \
                -ref "$TMAP" \
                -out "$RESAMPLED_MASK" \
                -applyxfm \
                -init "${FSLDIR}/etc/flirtsch/ident.mat" \
                -interp nearestneighbour 2>/dev/null

            # binarise resampled mask to remove interpolation artefacts
            BINARY_MASK="${TMPDIR}/${subject_id}_${roi_name}_binary.nii.gz"
            fslmaths "$RESAMPLED_MASK" -thr 0.5 -bin "$BINARY_MASK"

            # check mask has voxels after resampling
            NVOX=$(fslstats "$BINARY_MASK" -V | awk '{print $1}')

            if [ "${NVOX}" -eq 0 ] 2>/dev/null; then
                echo "  [!] Warning: Empty mask after resampling for ${roi_name}, skipping."
                MEAN_T="NA"
                MAX_T="NA"
            else
                # extract mean and max T-values within binary mask
                MEAN_T=$(fslstats "$TMAP" -k "$BINARY_MASK" -M | xargs)
                MAX_T=$(fslstats  "$TMAP" -k "$BINARY_MASK" -R | awk '{print $2}' | xargs)
            fi

            # write to CSV
            echo "${subject_id},${contrast_label},${roi_name},${MEAN_T},${MAX_T}" \
                >> "${OUTPUT_CSV}"

            echo "  [✓] ${roi_name} | Mean=${MEAN_T} | Max=${MAX_T}"
        done
    done
done

echo ""
echo "Done! Results saved to: ${OUTPUT_CSV}"
