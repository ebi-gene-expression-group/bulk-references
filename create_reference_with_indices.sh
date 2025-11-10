#!/bin/bash
#SBATCH --job-name=nf-core-references
#SBATCH --time=07-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=12000
#SBATCH --output=nf-core-references-%j-%x.out
#SBATCH --error=nf-core-references-%j-%x.err
#SBATCH --partition=production

source activate /hps/software/users/ma/service/isl/conda/config/envs/nextflow=25.10.0
timestamp=$(date +"%Y_%m_%d_%I_%M_%p")
NXF_ANSI_LOG=false
export NXF_TTY_WIDTH=999


i=$1

nextflow run nf-core-references/main.nf \
	--input ${i} \
	--outdir $BULK_REFERENCES_DIR \
	--tools "star,salmon,kallisto,faidx,createsequencedictionary,intervals,sizes,tabix" \
	-with-trace $BULK_REFERENCES_DIR/nf-core-references_trace_$i_$timestamp.txt 

