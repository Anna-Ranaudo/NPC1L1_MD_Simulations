#!/bin/bash
# PCA with GROMACS on NPC1L1 MD trajectories.
# Prerequisites: .nc → .xtc conversion done; concatenated trajectory already built with gmx trjcat.
# To concatenate: gmx trjcat -f run-md1/<traj>.xtc run-md2/<traj>.xtc ... -o <Nt><traj>.xtc

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR="/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25"
REPLICAS=("run-md1" "run-md2" "run-md3" "run-md4" "run-md5")
N_REP=${#REPLICAS[@]}

FIT_GROUP=1    # transmembrane region (group 1 in index.ndx)
ANAL_GROUP=0   # whole protein      (group 0 in index.ndx)

# Each entry: "system_id|pdb_code|condition|traj_name_no_ext"
# traj_name_no_ext: filename without .xtc (same for nc and xtc versions)
SYSTEMS=(
    "open_bound|6v3f|col|07-08-prot-lig-pbc"
    "open_apo|6v3f|no-col|07-08-prot-pbc"
    "closed_bound|6v3h|col|07-08-prot-lig-pbc"
    "closed_apo|6v3h|no-col|07-08-prot-pbc"
)

# =============================================================================
# MAIN LOOP OVER SYSTEMS
# =============================================================================
for SYSTEM in "${SYSTEMS[@]}"; do
    IFS='|' read -r SYS_ID PDB COND TRAJ_NAME <<< "$SYSTEM"

    DATA_DIR="${BASE_DIR}/${PDB}/data/${COND}"
    CAT_TRAJ="${DATA_DIR}/${N_REP}t${TRAJ_NAME}.xtc"   # e.g. 5t07-08-prot-lig-pbc.xtc
    REF_GRO="${DATA_DIR}/${REPLICAS[0]}/f0-${TRAJ_NAME}.gro"  # first frame of replica 1
    NDX="${DATA_DIR}/index.ndx"
    PCA_DIR="${DATA_DIR}/pca-gromacs"

    echo ""
    echo "============================================================"
    echo "  PCA – system: ${SYS_ID}"
    echo "============================================================"

    mkdir -p "${PCA_DIR}"
    # gmx covar/anaeig write output files (eigenvec.trr, average.pdb, …) to the cwd
    cd "${PCA_DIR}" || { echo "ERROR: cannot enter ${PCA_DIR}"; continue; }

    # -------------------------------------------------------------------------
    # STEP 1 – covar: compute cumulative PCA (fit on TM region, analyse whole protein)
    # -------------------------------------------------------------------------
    echo "${FIT_GROUP} ${ANAL_GROUP}" | gmx covar \
        -f "${CAT_TRAJ}" \
        -s "${REF_GRO}" \
        -n "${NDX}"
    # Outputs written here: eigenvec.trr  eigenval.xvg  average.pdb  covar.log  covar.xvg

    # -------------------------------------------------------------------------
    # STEP 2 – anaeig: project each replica onto the cumulative PC1-PC2 space
    # All replicas use the SAME cumulative average.pdb so projections are comparable.
    # -------------------------------------------------------------------------
    for i in "${!REPLICAS[@]}"; do
        REP="${REPLICAS[$i]}"
        REP_NUM=$((i + 1))
        echo "${FIT_GROUP} ${ANAL_GROUP}" | gmx anaeig \
            -f "${DATA_DIR}/${REP}/${TRAJ_NAME}.xtc" \
            -s average.pdb \
            -v eigenvec.trr \
            -eig eigenval.xvg \
            -first 1 -last 2 \
            -2d "plot12-md${REP_NUM}.xvg" \
            -n "${NDX}"
    done

    # -------------------------------------------------------------------------
    # STEP 3 – anaeig: RMSF along PC1 and PC2 (which residues drive each PC)
    # -------------------------------------------------------------------------
    for PC in 1 2; do
        echo "${ANAL_GROUP}" | gmx anaeig \
            -f "${CAT_TRAJ}" \
            -s "${REF_GRO}" \
            -v eigenvec.trr \
            -first ${PC} -last ${PC} \
            -rmsf "eigrmsf-${SYS_ID}-pc${PC}.xvg" \
            -n "${NDX}"
    done

    # -------------------------------------------------------------------------
    # STEP 4 – anaeig: extreme structures along PC1 and PC2 (100 interpolated frames)
    # see  https://manual.gromacs.org/2024.2/onlinehelp/gmx-anaeig.html
    # -------------------------------------------------------------------------
    for PC in 1 2; do
        echo "${FIT_GROUP} ${ANAL_GROUP}" | gmx anaeig \
            -f "${CAT_TRAJ}" \
            -s "${REF_GRO}" \
            -v eigenvec.trr \
            -first ${PC} -last ${PC} \
            -extr "${SYS_ID}-extr-pc${PC}.pdb" \
            -nframes 100 \
            -n "${NDX}"
    done

done
