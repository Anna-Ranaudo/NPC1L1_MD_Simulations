#!/bin/bash
# Joint PCA with GROMACS: builds a PC space shared by two systems and projects each of them
# into it, producing the 2D projections read by plot-2d-projections-pca_compare-systems.py.
# Companion to run-pca-with-gromacs.sh, which runs the PCA of each system separately.
#
# Prerequisites:
#   - per-system PCA already run: its average.pdb files are used as fit references here
#   - joint trajectories already built with gmx trjcat (2 x 5 replicas), e.g.
#     gmx trjcat -f <sysA>/run-md1/<traj>.xtc ... <sysB>/run-md5/<traj>.xtc -o 10t-<name>.xtc
#   - index.ndx present in the joint directory of each pair, on BASE_DIR
#
# Trajectories, fit references and index files are read from BASE_DIR; everything is written
# to final_data/ in the repository, next to the projections it produces. Note that eigenvec.trr
# and covar.log land there too but are excluded by .gitignore (*.trr, *.log), so they stay as
# local files and are not versioned.

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR="/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

FIT_GROUP=1    # transmembrane region, 'r_612-782_r_1092-1242' in index.ndx
ANAL_GROUP=0   # whole protein, 'Protein' in index.ndx

PER_SYS_DIR="all-pca-gromacs-07-08"       # per-system PCA outputs
TMFIT_DIR="pca-fit-transmemb-calc-all"    # transmembrane-fit generation of the PCA

N_REP=5                  # replicas per system  -> per-system trajectories are ${N_REP}t<base>.xtc
N_JOINT=$((2 * N_REP))   # replicas in a joint trajectory -> ${N_JOINT}t-<name>.xtc

# trajectory base name per condition (col = cholesterol-bound, no-col = apo)
traj_base() {
    case "$1" in
        col)    echo "07-08-prot-lig-pbc" ;;
        no-col) echo "07-08-prot-pbc" ;;
    esac
}

# the final_data subdirectories the two conformations and the two conditions map to
state_of() {
    case "$1" in
        6v3f) echo "open" ;;
        6v3h) echo "closed" ;;
    esac
}

binding_of() {
    case "$1" in
        col)    echo "bound" ;;
        no-col) echo "apo" ;;
    esac
}

# Both helpers read the index file from ${NDX}, set by the loop, and work in the cwd, which the
# loop sets to the final_data directory the projections belong to.

# Build the joint PC space from the concatenated trajectory of the two systems.
# Outputs into the cwd: eigenvec.trr  eigenval.xvg  average.pdb  covar.log
# Usage: build_space <joint_traj> <fit_reference>
build_space() {
    echo "${FIT_GROUP} ${ANAL_GROUP}" | gmx covar \
        -f "$1" \
        -s "$2" \
        -n "${NDX}"
}

# Project one system's cumulative trajectory onto the PC1-PC2 space built in the cwd.
# Usage: project <traj> <fit_reference> <output_xvg>
project() {
    echo "${FIT_GROUP} ${ANAL_GROUP}" | gmx anaeig \
        -f "$1" \
        -s "$2" \
        -v eigenvec.trr \
        -eig eigenval.xvg \
        -first 1 -last 2 \
        -2d "$3" \
        -n "${NDX}"
}

# =============================================================================
# PART 1 – bound vs apo: joint space of the col + no-col trajectories of one conformation
# =============================================================================
for PDB in 6v3f 6v3h; do
    CONF_DIR="${BASE_DIR}/${PDB}/data"
    JOINT_DIR="${CONF_DIR}/col-and-no-col"
    NDX="${JOINT_DIR}/${TMFIT_DIR}/index.ndx"
    WORK_DIR="${REPO_DIR}/final_data/$(state_of "${PDB}")/bound-vs-apo/pca-2d-proj"

    echo ""
    echo "============================================================"
    echo "  joint PCA – bound vs apo, ${PDB}"
    echo "============================================================"

    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}" || { echo "ERROR: cannot enter ${WORK_DIR}"; continue; }

    # the joint space is fit on the bound (col) per-system average
    build_space "${JOINT_DIR}/${N_JOINT}t-07-08.xtc" \
                "${CONF_DIR}/col/${PER_SYS_DIR}/${TMFIT_DIR}/average.pdb"

    # fit reference of each projection: the system's top-level per-system average,
    # i.e. not the TM-fit one used in PART 2 — kept as it was originally run
    for COND in col no-col; do
        SYS_DIR="${CONF_DIR}/${COND}/${PER_SYS_DIR}"
        project "${SYS_DIR}/${N_REP}t$(traj_base "${COND}").xtc" \
                "${SYS_DIR}/average.pdb" \
                "plot-${COND}-${N_REP}t.xvg"
    done
done

# =============================================================================
# PART 2 – open vs closed: joint space of the 6v3f + 6v3h trajectories of one condition
# =============================================================================
for COND in col no-col; do
    JOINT_DIR="${BASE_DIR}/compare-6v3f-6v3h/${COND}-6v3f-6v3h"
    NDX="${JOINT_DIR}/${TMFIT_DIR}/index.ndx"
    WORK_DIR="${REPO_DIR}/final_data/comparison_open-vs-closed/$(binding_of "${COND}")/pca-2d-proj"
    TRAJ="${N_REP}t$(traj_base "${COND}").xtc"

    echo ""
    echo "============================================================"
    echo "  joint PCA – open vs closed, ${COND}"
    echo "============================================================"

    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}" || { echo "ERROR: cannot enter ${WORK_DIR}"; continue; }

    # the joint space is fit on the open (6v3f) per-system average
    build_space "${JOINT_DIR}/${N_JOINT}t-${COND}.xtc" \
                "${BASE_DIR}/6v3f/data/${COND}/${PER_SYS_DIR}/${TMFIT_DIR}/average.pdb"

    # fit reference of each projection: the system's TM-fit average
    for PDB in 6v3f 6v3h; do
        SYS_DIR="${BASE_DIR}/${PDB}/data/${COND}/${PER_SYS_DIR}"
        project "${SYS_DIR}/${TRAJ}" \
                "${SYS_DIR}/${TMFIT_DIR}/average.pdb" \
                "plot-${COND}-${PDB}.xvg"
    done
done
