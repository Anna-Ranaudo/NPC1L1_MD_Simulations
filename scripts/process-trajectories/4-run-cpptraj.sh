#!/bin/bash
# =======================================================================================
# Script to process MD trajectories using Cpptraj
# Usage: ./4-run-cpptraj.sh
# =======================================================================================
# ---------------------------------------------------------------------------------------
# 1. Configuration Variables
# ---------------------------------------------------------------------------------------
# System details
SYSTEM="6v3f"           # Directory name for system (e.g., 6v3f)
SYSTEM_UPPER="6V3F"     # System name in filenames (e.g., 6V3F)
STATE="bound"           # State (e.g., bound, apo)
REPLICA="rep1"          # Replica directory (e.g., rep1)

# Processing parameters
STRIP_MASK=":1-1262"    # Mask for the protein/ligand to keep (strips everything else)

# Base paths
BASE_DIR="../NPC1L1_MD_Simulations"
SCRIPTS_DIR="${BASE_DIR}/scripts/process-trajectories"

# Input/Output paths
TOPOLOGY_FILE="${BASE_DIR}/input_structures/${SYSTEM}/${STATE}/bilayer_${SYSTEM_UPPER}-x-MD.amber_lipid.top"
TRAJ_DIR="${BASE_DIR}/run/${SYSTEM}/${STATE}/${REPLICA}"
OUTPUT_DIR="${TRAJ_DIR}/processing"  # Output directory for processed files

# ---------------------------------------------------------------------------------------
# 2. Setup and Checks
# ---------------------------------------------------------------------------------------
echo "Starting trajectory processing for ${SYSTEM} - ${STATE} - ${REPLICA}..."

# Check if directories exist
if [ ! -d "$TRAJ_DIR" ]; then
    echo "Error: Trajectory directory not found: $TRAJ_DIR"
    exit 1
fi

if [ ! -f "$TOPOLOGY_FILE" ]; then
    echo "Error: Topology file not found: $TOPOLOGY_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# ---------------------------------------------------------------------------------------
# 3. Generate Trajectory List
# ---------------------------------------------------------------------------------------
# Find all .nc files in the trajectory directory
# We sort them to ensure correct order (assuming filenames are sortable, e.g., 01_..., 02_...)
TRAJ_FILES=$(ls -v "$TRAJ_DIR"/07*.nc "$TRAJ_DIR"/08*.nc 2>/dev/null)

if [ -z "$TRAJ_FILES" ]; then
    echo "Error: No .nc files found in $TRAJ_DIR"
    exit 1
fi

# Create a temporary file with trajin commands
TRAJIN_CMD_FILE="${OUTPUT_DIR}/trajin_commands.tmp"
> "$TRAJIN_CMD_FILE" # Clear file

for traj in $TRAJ_FILES; do
    echo "trajin $traj" >> "$TRAJIN_CMD_FILE"
done

echo "Found $(wc -l < "$TRAJIN_CMD_FILE") trajectory files."

# ---------------------------------------------------------------------------------------
# 4. Prepare Cpptraj Input Files
# ---------------------------------------------------------------------------------------

# Function to process template
# Usage: process_template <template_file> <output_file>
process_template() {
    local template=$1
    local output=$2
    
    # Copy template to output
    cp "$template" "$output"
    
    # Replace placeholders using sed
    # We use | as delimiter to avoid issues with / in paths
    sed -i "s|AMBER_TOPOLOGY|$TOPOLOGY_FILE|g" "$output"
    sed -i "s|OUTPUT_DIR|$OUTPUT_DIR|g" "$output"
    sed -i "s|STRIP_MASK|$STRIP_MASK|g" "$output"
    
    # Special handling for TRAJIN_FILES
    if grep -q "TRAJIN_FILES" "$output"; then
        # Use sed to read the content of TRAJIN_CMD_FILE and replace the marker
        # The 'r' command appends file content *after* the match, so we delete the match line
        sed -i "/TRAJIN_FILES/r $TRAJIN_CMD_FILE" "$output"
        sed -i "/TRAJIN_FILES/d" "$output"
    fi
}

echo "Generating input files..."

# 1. Concatenate
process_template "${SCRIPTS_DIR}/1-concatenate.in" "${OUTPUT_DIR}/1-concatenate.run.in"

# 2. Strip Water/Membrane
process_template "${SCRIPTS_DIR}/2-stripwater-membrane.in" "${OUTPUT_DIR}/2-stripwater-membrane.run.in"

# 3. Treat PBCs
process_template "${SCRIPTS_DIR}/3-treat-pbcs.in" "${OUTPUT_DIR}/3-treat-pbcs.run.in"

# ---------------------------------------------------------------------------------------
# 5. Run Cpptraj
# ---------------------------------------------------------------------------------------

echo "Running Step 1: Concatenation..."
cpptraj -i "${OUTPUT_DIR}/1-concatenate.run.in" > "${OUTPUT_DIR}/1-concatenate.log"
if [ $? -ne 0 ]; then echo "Step 1 failed! Check log."; exit 1; fi

echo "Running Step 2: Stripping Water/Membrane..."
cpptraj -i "${OUTPUT_DIR}/2-stripwater-membrane.run.in" > "${OUTPUT_DIR}/2-stripwater-membrane.log"
if [ $? -ne 0 ]; then echo "Step 2 failed! Check log."; exit 1; fi

echo "Running Step 3: Treating PBCs..."
cpptraj -i "${OUTPUT_DIR}/3-treat-pbcs.run.in" > "${OUTPUT_DIR}/3-treat-pbcs.log"
if [ $? -ne 0 ]; then echo "Step 3 failed! Check log."; exit 1; fi

# ---------------------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------------------

rm "$TRAJIN_CMD_FILE"
# Optional: remove run input files
# rm "${OUTPUT_DIR}"/*.run.in

echo "Processing complete! Files are in $OUTPUT_DIR"
