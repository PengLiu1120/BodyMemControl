import os
import csv

# --- specify base paths ---
base_dir = "/mnt/volume/tacmem_workspace/derivatives/results/hippocampus_T_map_values"
analyses = ["TNT", "recall"]

print(f"starting collation from: {base_dir}")

for analysis in analyses:
    analysis_path = os.path.join(base_dir, analysis)
    output_csv = os.path.join(base_dir, f"hippocampus_stats_{analysis}.csv")
    
    if not os.path.isdir(analysis_path):
        print(f"Skipping {analysis}: Directory not found.")
        continue

    analysis_data = []

    # --- CHANGE 1: SORT THE FILES ---
    # get all .txt files produced by your Bash script and sort them
    files = [f for f in os.listdir(analysis_path) if f.endswith("_stats.txt")]
    files.sort() # <--- This ensures sub-03 comes before sub-30
    
    for filename in files:
        # expected filename: sub-03_spmT_0001_left_stats.txt
        parts = filename.replace("_stats.txt", "").split("_")
        
        if len(parts) < 4:
            continue
            
        # --- CHANGE 2: REMOVE 'sub-' ---
        # Before: subject = parts[0]
        # After: strip 'sub-' from the string
        subject = parts[0].replace("sub-", "") # <--- Result: "03" instead of "sub-03"
        
        tmap_id = f"{parts[1]}_{parts[2]}" 
        hemi = parts[3]         

        file_path = os.path.join(analysis_path, filename)
        
        with open(file_path, "r") as f:
            lines = f.readlines()
            # parse data rows (skipping headers)
            for line in lines[4:]:
                if "|" in line:
                    cols = [c.strip() for c in line.split("|")]
                    
                    if len(cols) == 4:
                        analysis_data.append({
                            "Subject": subject,
                            "Analysis": analysis,
                            "Hemisphere": hemi,
                            "Condition_Map": tmap_id,
                            "Subfield_ID": cols[0],
                            "Subfield_Name": cols[1],
                            "Mean_T": cols[2],
                            "Max_T": cols[3]
                        })

    # --- write to specific CSV for this analysis ---
    headers = ["Subject", "Analysis", "Hemisphere", "Condition_Map", "Subfield_ID", "Subfield_Name", "Mean_T", "Max_T"]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        writer.writerows(analysis_data)
        
    print(f"successfully created: {output_csv} ({len(analysis_data)} rows)")

print("processing completed.")
