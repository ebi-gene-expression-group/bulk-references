#!/bin/bash
set -euo pipefail
#SBATCH --job-name=nf-core-references
#SBATCH --time=07-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=12000
#SBATCH --output=nf-core-references-%j-%x.out
#SBATCH --error=nf-core-references-%j-%x.err
#SBATCH --partition=production


timestamp=$(date +"%Y_%m_%d_%I_%M_%p")
NXF_ANSI_LOG=false
export NXF_TTY_WIDTH=999


# Slurm executes a submitted script from a temporary spool copy, so
# BASH_SOURCE[0] is not the repository path. Prefer the directory from which
# sbatch was submitted and retain the local-script fallback for direct runs.
SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
cd "${SCRIPT_DIR}"

# Preserve the historical reference build by default. Set this to true when
# generating the separate STAR index required by Parabricks rna_fq2bam.
PARABRICKS_STAR_INDEX=${PARABRICKS_STAR_INDEX:-false}

if [[ $# -ne 1 ]]; then
	echo "Usage: $0 <species-yaml>" >&2
	exit 2
fi

export SPECIES=$1

export WORKSUBDIR=$(basename $SPECIES .yaml)

nextflow_config_args=()
if [[ "${PARABRICKS_STAR_INDEX,,}" == "true" ]]; then
	nextflow_config_args=(-c "${SCRIPT_DIR}/parabricks.config")
fi

nextflow run nf-core-references/main.nf \
	--input ${SPECIES} \
	--outdir ${BULK_REFERENCES_DIR} \
	--tools "star,salmon,kallisto,faidx,createsequencedictionary,intervals,sizes,tabix" \
	"${nextflow_config_args[@]}" \
	-with-trace ${BULK_REFERENCES_DIR}/nf-core-references_trace_${SPECIES}_$timestamp.txt 
