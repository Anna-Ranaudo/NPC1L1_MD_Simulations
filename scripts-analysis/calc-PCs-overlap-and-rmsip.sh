#!/bin/bash
# Compares the essential subspaces of the NPC1L1 MD replicas, for every pair of replicas of
# each system, with gmx anaeig -over on the eigenvectors of the per-replica PCA.
#
# Prerequisite: per-replica PCA present in <system>/run-md<N>/pca-gromacs-07-08/
#
# Two metrics come out of the same gmx call, and they behave differently with respect to
# N_EIGENVECTORS, so they are written to two different directories:
#
# final_data/<state>/<binding>/overlap-<N>PCs-single-replica/   DEPENDS on N_EIGENVECTORS
#   eigenvectors-overlap-md<i>-md<j>.xvg   the curve written by gmx: overlap between the first
#                                     N_EIGENVECTORS eigenvectors of replica i and the first x
#                                     eigenvectors of replica j, for x = 1 ... 3783. Note that
#                                     only the eigenvectors of replica i are restricted to N:
#                                     the curve itself spans the whole set. The point at
#                                     x = N_EIGENVECTORS is RMSIP^2, the last point is 1.000 by
#                                     construction.
#   rmsip-matrix.csv                  RMSIP of every pair, i.e. sqrt of the curve read at
#                                     x = N_EIGENVECTORS, as a replica x replica matrix
#
# final_data/<state>/<binding>/covariance-overlap-single-replica/   does NOT depend on N
#   covariance-overlap-md<i>-md<j>.txt  the gmx anaeig output, which carries the overlap of the
#                                     covariance matrices ('normalized' and 'shape'), the
#                                     traces and the command line of the run
#   covariance-overlap-matrix.csv     the 'normalized' values as a replica x replica matrix
#
# The covariance overlap uses all 3783 dimensions and the eigenvalues, hence it is the same
# whatever N_EIGENVECTORS is: it lives outside the per-N directory so that running this script
# at a different N does not produce a second identical copy of it.
#
# Both matrices are extracted from the files above, so they can be refreshed without re-running
# gmx: a pair whose curve is already present is skipped.

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR="/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

N_REP=5                 # replicas per system, in run-md1 ... run-md<N_REP>
N_EIGENVECTORS=2        # number of considered eigenvectors of replica i

PER_REP_DIR="pca-gromacs-07-08"   # per-replica PCA, inside run-md<N>

# Each entry: "system_id|pdb_code|condition"; system_id is <state>_<binding>, which also gives
# the final_data subdirectory the results are written to.
SYSTEMS=(
    "open_bound|6v3f|col"
    "open_apo|6v3f|no-col"
    "closed_bound|6v3h|col"
    "closed_apo|6v3h|no-col"
)

# Optional arguments: <N_EIGENVECTORS> [system_id ...], to run a different subspace dimension
# or a subset of the systems, e.g.   calc-PCs-overlap-and-rmsip.sh 5 open_apo
# With no arguments the values set above are used, for all four systems.
if [ $# -ge 1 ]; then
    N_EIGENVECTORS="$1"
    shift
fi

if [ $# -ge 1 ]; then
    SELECTED=()
    for SYSTEM in "${SYSTEMS[@]}"; do
        for WANTED in "$@"; do
            [ "${SYSTEM%%|*}" = "${WANTED}" ] && SELECTED+=("${SYSTEM}")
        done
    done
    if [ ${#SELECTED[@]} -eq 0 ]; then
        echo "ERROR: no system matches '$*'"
        exit 1
    fi
    SYSTEMS=("${SELECTED[@]}")
fi

# Write a replica x replica matrix from the values held in the VALUE associative array.
# Usage: write_matrix <output_csv>
write_matrix() {
    local header="" row
    for j in $(seq 1 ${N_REP}); do header="${header},md${j}"; done
    echo "${header}" > "$1"

    for i in $(seq 1 ${N_REP}); do
        row="md${i}"
        for j in $(seq 1 ${N_REP}); do
            if [ "${i}" -eq "${j}" ]; then
                row="${row},1.000"
            else
                row="${row},${VALUE[${i},${j}]}"
            fi
        done
        echo "${row}" >> "$1"
    done
}

# =============================================================================
# MAIN LOOP OVER SYSTEMS
# =============================================================================
for SYSTEM in "${SYSTEMS[@]}"; do
    IFS='|' read -r SYS_ID PDB COND <<< "$SYSTEM"

    DATA_DIR="${BASE_DIR}/${PDB}/data/${COND}"
    SYS_OUT="${REPO_DIR}/final_data/${SYS_ID%_*}/${SYS_ID#*_}"
    OVERLAP_DIR="${SYS_OUT}/overlap-${N_EIGENVECTORS}PCs-single-replica"
    COV_DIR="${SYS_OUT}/covariance-overlap-single-replica"

    echo ""
    echo "============================================================"
    echo "  subspace overlap – system: ${SYS_ID}"
    echo "============================================================"

    mkdir -p "${OVERLAP_DIR}" "${COV_DIR}"

    # -------------------------------------------------------------------------
    # STEP 1 – gmx anaeig on every pair of replicas
    # -------------------------------------------------------------------------
    for i in $(seq 1 $((N_REP - 1))); do
        for j in $(seq $((i + 1)) ${N_REP}); do
            PAIR="md${i}-md${j}"
            CURVE="${OVERLAP_DIR}/eigenvectors-overlap-${PAIR}.xvg"
            COV_FILE="${COV_DIR}/covariance-overlap-${PAIR}.txt"

            if [ -s "${CURVE}" ]; then
                echo "  ${PAIR}: curve already present, gmx skipped"
                continue
            fi

            # the covariance overlap does not depend on N_EIGENVECTORS, so an already captured
            # one is kept as it is and this run only produces the curve
            [ -s "${COV_FILE}" ] && CAPTURE="/dev/null" || CAPTURE="${COV_FILE}"

            echo "  comparing ${PAIR}"
            gmx anaeig \
                -v "${DATA_DIR}/run-md${i}/${PER_REP_DIR}/eigenvec.trr" \
                -v2 "${DATA_DIR}/run-md${j}/${PER_REP_DIR}/eigenvec.trr" \
                -first 1 -last ${N_EIGENVECTORS} \
                -over "${CURVE}" \
                > "${CAPTURE}" 2>&1
        done
    done

    # -------------------------------------------------------------------------
    # STEP 2 – RMSIP: the curve at x = N_EIGENVECTORS is RMSIP^2
    # -------------------------------------------------------------------------
    unset VALUE; declare -A VALUE
    for i in $(seq 1 $((N_REP - 1))); do
        for j in $(seq $((i + 1)) ${N_REP}); do
            CURVE="${OVERLAP_DIR}/eigenvectors-overlap-md${i}-md${j}.xvg"

            if [ ! -s "${CURVE}" ]; then
                echo "  WARNING: $(basename "${CURVE}") missing, left out of rmsip-matrix.csv"
                continue
            fi

            OVERLAP=$(awk -v n="${N_EIGENVECTORS}" '$1 == n {print $2}' "${CURVE}")
            VALUE[${i},${j}]=$(awk -v o="${OVERLAP}" 'BEGIN {printf "%.3f", sqrt(o)}')
            VALUE[${j},${i}]="${VALUE[${i},${j}]}"
        done
    done
    write_matrix "${OVERLAP_DIR}/rmsip-matrix.csv"

    # -------------------------------------------------------------------------
    # STEP 3 – covariance overlap: the 'normalized' value of each pair
    # -------------------------------------------------------------------------
    unset VALUE; declare -A VALUE
    for i in $(seq 1 $((N_REP - 1))); do
        for j in $(seq $((i + 1)) ${N_REP}); do
            COV_FILE="${COV_DIR}/covariance-overlap-md${i}-md${j}.txt"

            if [ ! -s "${COV_FILE}" ]; then
                echo "  WARNING: $(basename "${COV_FILE}") missing, left out of the matrix"
                continue
            fi

            VALUE[${i},${j}]=$(awk '/normalized:/ {printf "%.3f", $2}' "${COV_FILE}")
            VALUE[${j},${i}]="${VALUE[${i},${j}]}"
        done
    done
    write_matrix "${COV_DIR}/covariance-overlap-matrix.csv"

    echo "  written: rmsip-matrix.csv and covariance-overlap-matrix.csv"
done
