#!/usr/bin/env bash

# USAGE
# bash bbsplit_index.sh genomes.csv <outpath>
set -euo pipefail

FASTA_LIST="${1:-}"
OUTDIR="${2:-results/bbsplit_index}"
MEM_GB="${MEM_GB:-8}"
THREADS="${THREADS:-4}"

if [[ -z "$FASTA_LIST" ]]; then
    echo "Usage: $0 refs.csv [outdir]"
    exit 1
fi

if [[ ! -f "$FASTA_LIST" ]]; then
    echo "Error: CSV file not found: $FASTA_LIST"
    exit 1
fi

mkdir -p "$OUTDIR"

declare -a REF_ARGS=()

while IFS=, read -r short_name fasta_path; do
    # skip empty lines
    [[ -z "${short_name// }" ]] && continue
    [[ -z "${fasta_path// }" ]] && continue

    # trim whitespace
    short_name="$(echo "$short_name" | xargs)"
    fasta_path="$(echo "$fasta_path" | xargs)"

    if [[ ! -f "$fasta_path" ]]; then
        echo "Error: FASTA not found for '$short_name': $fasta_path"
        exit 1
    fi

    REF_ARGS+=("ref_${short_name}=${fasta_path}")
done < "$FASTA_LIST"

if [[ ${#REF_ARGS[@]} -eq 0 ]]; then
    echo "Error: No valid references found in $FASTA_LIST"
    exit 1
fi

echo "Building BBSplit index in: $OUTDIR"
echo "Using memory: ${MEM_GB}g"
echo "Using threads: ${THREADS}"

echo "bbsplit.sh -Xmx${MEM_GB}g ${REF_ARGS[@]} path=${OUTDIR} threads=${THREADS} ${OUTDIR}/index_build.log" 

bbsplit.sh \
    "-Xmx${MEM_GB}g" \
    "${REF_ARGS[@]}" \
    "path=${OUTDIR}" \
    "threads=${THREADS}" \
    > "${OUTDIR}/index_build.log" 2>&1

echo "BBSplit index build completed."
echo "Index directory: $OUTDIR"
echo "Log file: ${OUTDIR}/index_build.log"
