#!/bin/bash
#
# GROMACS MD Analysis Script (RMSD, Rg, SASA, RMSF, HBond)
#
# This script automates the calculation of six key molecular dynamics metrics
# using GROMACS tools and converts the resulting .xvg files to .csv format.
# It also generates a Python script (ana.py) for easy data loading.
#
# ASSUMPTIONS FOR AUTOMATION:
# 1. The script assumes the input files are named 'topol.tpr' and 'prod.xtc'.
# 2. It assumes GROMACS index group '1' is 'Protein' (used for Rg, SASA, RMSF, HBond).
# 3. It assumes GROMACS index group '3' is 'C-alpha' (used for C-alpha RMSF).
# 4. It assumes GROMACS index group '4' is 'Backbone' (used for RMSD fitting/calculation).
#
# To run: chmod +x gromacs_analysis.sh && ./gromacs_analysis.sh

# --- Configuration ---
TPR_FILE="step5_1.tpr"
XTC_FILE="nopbc.xtc"
OUTPUT_DIR="analysis_results"

# --- Function to Convert XVG to CSV ---
# Removes GROMACS comments (@, #) and replaces spaces with commas.
convert_xvg_to_csv() {
    local xvg_file="$1"
    local csv_file="$2"

    if [ -f "$xvg_file" ]; then
        # Use awk to filter comments and replace delimiters
        # $0 !~ /^[#@]/ : Skips lines starting with # or @
        # gsub(/[[:space:]]+/, ",") : Replaces multiple spaces/tabs with a single comma
        # $1=$1 : Forces awk to re-evaluate the line, applying the gsub substitution
        # The result might contain a leading comma if the first column is time/index, which is fine for CSV parsing.
        awk '$0 !~ /^[#@]/ {gsub(/[[:space:]]+/, ","); $1=$1; print}' "$xvg_file" > "$csv_file"
        rm "$xvg_file" # Clean up temporary XVG file
        echo "Successfully converted to $csv_file"
    else
        echo "Error: Temporary file $xvg_file not found for conversion."
    fi
}
# --------------------------------------


# --- Setup and Checks ---

echo "--- Starting GROMACS Analysis Script ---"

# Check for GROMACS availability
if ! command -v gmx &> /dev/null
then
    echo "Error: GROMACS (gmx) could not be found. Please ensure it is in your PATH."
    exit 1
fi

# Check for input files
if [ ! -f "$TPR_FILE" ] || [ ! -f "$XTC_FILE" ]; then
    echo "Error: Required files not found."
    echo "Expected: $TPR_FILE and $XTC_FILE in the current directory."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "Output directory created: $OUTPUT_DIR"

# --- 1. Root Mean Square Deviation (RMSD) Calculation ---
echo ""
echo "--- 1. Calculating RMSD (Fit on Backbone, calculate on Backbone) ---"
echo "Piping selections '4' (Backbone) for reference and '4' (Backbone) for RMSD."

TEMP_RMSD_XVG="$OUTPUT_DIR/rmsd.temp.xvg"
FINAL_RMSD_CSV="$OUTPUT_DIR/rmsd.csv"

# Run GROMACS command to create temporary XVG
echo -e "4\n4" | gmx rms -s "$TPR_FILE" -f "$XTC_FILE" -o "$TEMP_RMSD_XVG" -fit rot+trans -tu ns 2> /dev/null

if [ $? -eq 0 ]; then
    convert_xvg_to_csv "$TEMP_RMSD_XVG" "$FINAL_RMSD_CSV"
    echo "RMSD calculated and converted successfully: $FINAL_RMSD_CSV"
else
    echo "Warning: RMSD calculation may have failed or had non-critical warnings. CSV file not generated."
fi

# --- 2. Radius of Gyration (Rg) Calculation ---
echo ""
echo "--- 2. Calculating Radius of Gyration (Rg) for Protein ---"
echo "Piping selection '1' (Protein) for Rg calculation."

TEMP_RG_XVG="$OUTPUT_DIR/rg.temp.xvg"
FINAL_RG_CSV="$OUTPUT_DIR/rg.csv"

# Run GROMACS command to create temporary XVG
echo "1" | gmx gyrate -s "$TPR_FILE" -f "$XTC_FILE" -o "$TEMP_RG_XVG" -tu ns 2> /dev/null

if [ $? -eq 0 ]; then
    convert_xvg_to_csv "$TEMP_RG_XVG" "$FINAL_RG_CSV"
    echo "Rg calculated and converted successfully: $FINAL_RG_CSV"
else
    echo "Warning: Rg calculation may have failed or had non-critical warnings. CSV file not generated."
fi


# --- 3. Solvent Accessible Surface Area (SASA) Calculation ---
echo ""
echo "--- 3. Calculating SASA for Protein ---"
echo "Piping selection '1' (Protein) for SASA calculation."

TEMP_SASA_XVG="$OUTPUT_DIR/sasa.temp.xvg"
FINAL_SASA_CSV="$OUTPUT_DIR/sasa.csv"

# Run GROMACS command to create temporary XVG
echo "1" | gmx sasa -s "$TPR_FILE" -f "$XTC_FILE" -o "$TEMP_SASA_XVG" -tu ns 2> /dev/null

if [ $? -eq 0 ]; then
    convert_xvg_to_csv "$TEMP_SASA_XVG" "$FINAL_SASA_CSV"
    echo "SASA calculated and converted successfully: $FINAL_SASA_CSV"
else
    echo "Warning: SASA calculation may have failed or had non-critical warnings. CSV file not generated."
fi

# --- 4. Root Mean Square Fluctuation (RMSF) Calculation (Protein) ---
echo ""
echo "--- 4. Calculating RMSF (Fluctuation per residue) for entire Protein ---"
echo "Piping selection '1' (Protein) for RMSF calculation."

TEMP_RMSF_PROT_XVG="$OUTPUT_DIR/rmsf_protein.temp.xvg"
FINAL_RMSF_PROT_CSV="$OUTPUT_DIR/rmsf_protein.csv"

# Run GROMACS command to create temporary XVG
# Note: gmx rmsf implicitly fits the structure to the average coordinates.
echo "1" | gmx rmsf -s "$TPR_FILE" -f "$XTC_FILE" -o "$TEMP_RMSF_PROT_XVG" -tu ns 2> /dev/null

if [ $? -eq 0 ]; then
    convert_xvg_to_csv "$TEMP_RMSF_PROT_XVG" "$FINAL_RMSF_PROT_CSV"
    echo "RMSF (Protein) calculated and converted successfully: $FINAL_RMSF_PROT_CSV (Residue Index vs Fluctuation)"
else
    echo "Warning: RMSF (Protein) calculation may have failed or had non-critical warnings. CSV file not generated."
fi

# --- 5. Hydrogen Bond (HBond) Calculation ---
echo ""
echo "--- 5. Calculating Intra-Molecular Hydrogen Bonds (Protein-Protein) ---"
echo "Piping selections '1' (Donor) and '1' (Acceptor)."

TEMP_HBOND_XVG="$OUTPUT_DIR/hbond.temp.xvg"
FINAL_HBOND_CSV="$OUTPUT_DIR/hbond.csv"

# Run GROMACS command to create temporary XVG
# -num: Output file for the total number of hydrogen bonds over time.
echo -e "1\n1" | gmx hbond -s "$TPR_FILE" -f "$XTC_FILE" -num "$TEMP_HBOND_XVG" -tu ns 2> /dev/null

if [ $? -eq 0 ]; then
    convert_xvg_to_csv "$TEMP_HBOND_XVG" "$FINAL_HBOND_CSV"
    echo "Hydrogen Bond count calculated and converted successfully: $FINAL_HBOND_CSV (Time vs Total H-Bonds)"
else
    echo "Warning: Hydrogen Bond calculation may have failed or had non-critical warnings. CSV file not generated."
fi

# --- 6. Root Mean Square Fluctuation (RMSF) Calculation (C-Alpha) ---
echo ""
echo "--- 6. Calculating RMSF (Fluctuation per residue) for C-Alpha atoms ---"
echo "Piping selection '3' (C-Alpha) for RMSF calculation."

TEMP_RMSF_CALPHA_XVG="$OUTPUT_DIR/rmsf_calpha.temp.xvg"
FINAL_RMSF_CALPHA_CSV="$OUTPUT_DIR/rmsf_calpha.csv"

# Run GROMACS command to create temporary XVG
# Note: gmx rmsf implicitly fits the structure to the average coordinates.
echo "3" | gmx rmsf -s "$TPR_FILE" -f "$XTC_FILE" -o "$TEMP_RMSF_CALPHA_XVG" -tu ns 2> /dev/null

if [ $? -eq 0 ]; then
    convert_xvg_to_csv "$TEMP_RMSF_CALPHA_XVG" "$FINAL_RMSF_CALPHA_CSV"
    echo "RMSF (C-Alpha) calculated and converted successfully: $FINAL_RMSF_CALPHA_CSV (Residue Index vs Fluctuation)"
else
    echo "Warning: RMSF (C-Alpha) calculation may have failed or had non-critical warnings. CSV file not generated."
fi


# --- 7. Create Python Loader Script (ana.py) ---
echo ""
echo "--- 7. Creating Python utility script: ana.py ---"

# Use a heredoc to write the contents of ana.py to a file
cat << EOF > $OUTPUT_DIR/ana.py
import pandas as pd
import os
import glob
from typing import Dict
import pickle

def load_csv_files(folder_path: str) -> Dict[str, pd.DataFrame]:
    """
    Loads all .csv files from a specified folder into a dictionary of pandas DataFrames.

    Args:
        folder_path: The path to the directory containing the CSV files.

    Returns:
        A dictionary where keys are the filenames (without extension) and
        values are the corresponding pandas DataFrames.
    """
    if not os.path.isdir(folder_path):
        print(f"Error: Folder not found at '{folder_path}'")
        return {}

    # Pattern to match all CSV files in the folder
    search_pattern = os.path.join(folder_path, "*.csv")
    csv_files = glob.glob(search_pattern)

    if not csv_files:
        print(f"Warning: No CSV files found in '{folder_path}'")
        return {}

    data_frames = {}

    for file_path in csv_files:
        try:
            # Extract the filename without the path and extension to use as the dictionary key
            base_name = os.path.splitext(os.path.basename(file_path))[0]

            # Read the CSV file. The GROMACS-derived CSVs have no header.
            df = pd.read_csv(file_path, header=None)
            data_frames[base_name] = df
            print(f"Loaded: {base_name} ({len(df)} rows)")

        except Exception as e:
            print(f"Error loading file {file_path}: {e}")

    return data_frames

# --- Example Usage ---
if __name__ == '__main__':
    print("\n--- Python Loader Script Demonstration ---")

    # Check if the analysis_results folder exists (created by the main Bash script)
    actual_results_dir = "analysis_results"
    if os.path.isdir(actual_results_dir):
        print(f"Attempting to load data from: {actual_results_dir}")
        md_data = load_csv_files(actual_results_dir)

        if md_data:
            print("\nSuccessfully loaded DataFrames. Dictionary keys:")
            print(list(md_data.keys()))
            with open('ana.pickle', 'wb') as file:
                pickle.dump(md_data, file)

            

        else:
            print("No data loaded. Run the main script first to generate data.")
    else:
        print(f"Directory '{actual_results_dir}' not found. Run the main analysis script first.")

EOF
echo "Python script 'ana.py' created in the current directory."


# --- Summary ---
echo ""
echo "--- Analysis Complete ---"
echo "Results saved in the '$OUTPUT_DIR/' directory in CSV format:"
ls -l "$OUTPUT_DIR"/*.csv
echo ""
echo "Python utility script 'ana.py' has been generated. You can now load your data easily:"
echo "    python ana.py"
echo "or inside a Jupyter notebook:"
echo "    from ana import load_csv_files"
echo "    data = load_csv_files('analysis_results')"
echo "--------------------------"
echo "You now have 6 key analysis files and a data loader."
