#!/bin/bash
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


export SPECIES=$1

export WORKSUBDIR=$(basename $SPECIES .yaml)

nextflow run nf-core-references/main.nf \
	--input ${SPECIES} \
	--outdir ${BULK_REFERENCES_DIR} \
	--tools "star,salmon,kallisto,faidx,createsequencedictionary,intervals,sizes,tabix" \
	-with-trace ${BULK_REFERENCES_DIR}/nf-core-references_trace_${SPECIES}_$timestamp.txt 

