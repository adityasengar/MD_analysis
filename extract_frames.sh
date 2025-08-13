#!/bin/bash
#
# This script extracts a specified percentage of frames from a GROMACS trajectory file (.xtc).
# It is designed to be run from within a simulation directory (e.g., 'run1', 'run2', etc.)
# where 'step7_1.xtc' and 'step7_1.tpr' files are present.
#
# How it works:
# 1. Sets a variable for the desired percentage of frames to extract.
# 2. Loads the necessary GROMACS modules for execution.
# 3. Checks for the existence of 'step7_1.xtc' and 'step7_1.tpr' in the current directory.
# 4. Uses 'gmx check' to determine the total number of frames in 'step7_1.xtc'.
# 5. Calculates an 'INTERVAL' to skip frames, ensuring the desired percentage is extracted.
# 6. Executes 'gmx trjconv' to perform the extraction. It automatically selects the 'System' group (which includes protein, ligand, solvent, etc.) by piping '0' to the prompt.
# 7. Saves the reduced trajectory to a new .xtc file in the same directory.
#

# Percentage of frames to extract (e.g., 5 for 5%)
PERCENTAGE_TO_EXTRACT=5

# Load relevant modules
module purge
module load gcc cuda hwloc cmake python fftw openblas intel-oneapi-mkl gromacs

# Define files relative to the current directory
XTC_FILE="step7_1.xtc"
TPR_FILE="step7_1.tpr"
OUTPUT_XTC="step7_prot_SOLU_${PERCENTAGE_TO_EXTRACT}percent.xtc"

echo "Processing current directory..."

if [ -f "${XTC_FILE}" ] && [ -f "${TPR_FILE}" ]; then
    # Get total number of frames
    TOTAL_FRAMES=$(gmx check -f "${XTC_FILE}" 2>&1 | grep "^Step" | awk '{print $2}')

    if [ -z "${TOTAL_FRAMES}" ]; then
        echo "Error: Could not determine total frames for ${XTC_FILE}. Skipping."
    else
        # Calculate the interval to extract the desired percentage of frames
        # For example, if PERCENTAGE_TO_EXTRACT is 5, we want 1 frame every 20 frames (100/5)
        INTERVAL=$(( 100 / PERCENTAGE_TO_EXTRACT ))

        if [ "${INTERVAL}" -eq 0 ]; then
            echo "Error: Invalid percentage or interval calculation. Skipping."
        else
            echo "Total frames in ${XTC_FILE}: ${TOTAL_FRAMES}"
            echo "Extracting 1 frame every ${INTERVAL} frames to get approximately ${PERCENTAGE_TO_EXTRACT}% of data."

            # Extract frames, piping '0' for SOLU selection
            # gmx trjconv will ask for a group, '0' corresponds to 'SOLU' based on common GROMACS index files.
            echo 0 | gmx trjconv -f "${XTC_FILE}" -s "${TPR_FILE}" -o "${OUTPUT_XTC}" -skip "${INTERVAL}"

            if [ $? -eq 0 ]; then
                echo "Successfully extracted frames to ${OUTPUT_XTC}"
            else
                echo "Error extracting frames for ${XTC_FILE}"
            fi
	fi
    fi
else
    echo "Required files (${XTC_FILE} or ${TPR_FILE}) not found in the current directory. Skipping."
fi

echo "Script finished."

