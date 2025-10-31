#!/bin/bash
#SBATCH --job-name=nf-core-references
#SBATCH --time=01-23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=12000
#SBATCH --output=nf-core-references-%j-%x.out
#SBATCH --error=nf-core-references-%j-%x.err
#SBATCH --partition=production

source activate /hps/software/users/ma/service/isl/conda/config/envs/nextflow=25.10.0

NXF_ANSI_LOG=false

nextflow run nf-core-references/main.nf \
	--input datasheet.yaml \
	--outdir ./references_indices 

