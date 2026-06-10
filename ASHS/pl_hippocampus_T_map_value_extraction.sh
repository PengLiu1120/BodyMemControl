#!/bin/bash
# --- 1. base paths ---
BASE_DIR="/mnt/volume/tacmem_workspace/derivatives"
RESULT_DIR="${BASE_DIR}/results/hippocampus_T_map_values"

# --- 2. specify subject and analysis lists ---
subjects=("sub-03" "sub-05" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-20" "sub-24" "sub-27" "sub-28" "sub-29" "sub-30" "sub-32" "sub-33" "sub-34" "sub-35" "sub-36" "sub-37" "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-47")

analyses=("TNT" "recall")

# define label names
declare -A names
names=([1]="CA1" [2]="CA2" [3]="DG" [4]="CA3" [5]="Tail" [6]="Sulcus" [7]="MISC" [8]="SUB" [9]="ERC" [10]="BA35" [11]="BA36" [12]="PHC" [13]="MISC2" [14]="CS")

# fix for decimal point interpretation in some environments
export LC_NUMERIC="C"

# create a temp dir for binary masks, clean up on exit
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# --- 3. T map value extraction ---
for subject_id in "${subjects[@]}"; do
    echo "Processing ${subject_id}"
    for analysis in "${analyses[@]}"; do
        SUBJ_ANAL_DIR="${BASE_DIR}/results/individual_analyses/${analysis}/${subject_id}"

        if [ ! -d "$SUBJ_ANAL_DIR" ]; then
            echo "  [!] Missing analysis directory for ${subject_id} / ${analysis}, skipping."
            continue
        fi

        # ensure output directory exists
        mkdir -p "${RESULT_DIR}/${analysis}"

        for TMAP in "${SUBJ_ANAL_DIR}"/spmT_*_in_T2star_space.nii.gz; do
            [ -e "$TMAP" ] || continue

            MAP_NAME=$(basename "$TMAP" _in_T2star_space.nii.gz)

            for hemi in left right; do
                SEG="${BASE_DIR}/ASHS_seg/hippocampus_seg/${subject_id}/final/${subject_id}_${hemi}_lfseg_corr_usegray.nii.gz"
                OUTPUT_FILE="${RESULT_DIR}/${analysis}/${subject_id}_${MAP_NAME}_${hemi}_stats.txt"

                if [ ! -f "$SEG" ]; then
                    echo "  [!] Missing $hemi segmentation for ${subject_id}, skipping."
                    continue
                fi

                {
                    echo "Subject: ${subject_id} | Analysis: ${analysis} | Map: ${MAP_NAME} | Hemisphere: ${hemi}"
                    echo "------------------------------------------------------------"
                    printf "%-3s | %-20s | %-8s | %-8s\n" "ID" "Subfield Name" "Mean_T" "Max_T"
                    echo "------------------------------------------------------------"

                    for id in {1..14}; do
                        # create a binary mask for this specific label
                        LABEL_MASK="${TMPDIR}/mask_${subject_id}_${hemi}_label${id}.nii.gz"
                        fslmaths "$SEG" -thr $id -uthr $id -bin "$LABEL_MASK"

                        # check if the mask has any voxels
                        NVOX=$(fslstats "$LABEL_MASK" -V | awk '{print $1}')

                        if [ "${NVOX}" -eq 0 ] 2>/dev/null; then
                            CLEAN_MEAN="0.0000"
                            CLEAN_MAX="0.0000"
                        else
                            # compute stats using binary mask — no intensity thresholds, so negatives are preserved
                            MEAN=$(fslstats "$TMAP" -k "$LABEL_MASK" -M)
                            MAX=$(fslstats "$TMAP" -k "$LABEL_MASK" -R | awk '{print $2}')
                            CLEAN_MEAN=$(echo "$MEAN" | xargs)
                            CLEAN_MAX=$(echo "$MAX" | xargs)
                        fi

                        printf "%02d  | %-20s | %-8.4f | %-8.4f\n" "$id" "${names[$id]}" "$CLEAN_MEAN" "$CLEAN_MAX"
                    done

                } > "$OUTPUT_FILE"

                echo "  [✓] Written: $(basename $OUTPUT_FILE)"
            done
        done
    done
done

echo ""
echo "Done! All stats files are in $RESULT_DIR"
