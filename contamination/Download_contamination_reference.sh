#!/usr/bin/env bash +x
set -euo pipefail

##############################################################################
# Download contamination references for nf-core/rnaseq
# Supports both BBSplit (alignment-based) and Kraken2 (k-mer-based) approaches
##############################################################################

: "${BULK_REFERENCES_DIR:?BULK_REFERENCES_DIR must be set}"

# Configuration
OUTDIR="$BULK_REFERENCES_DIR/contamination"
FUNGI_DIR="$OUTDIR/fungi_refseq"
BACT_DIR="$OUTDIR/bacteria_refseq"
META_DIR="$OUTDIR/metadata"
MODE="${MODE:-bbsplit}"  # Options: bbsplit, kraken2, both

mkdir -p "$OUTDIR" "$FUNGI_DIR" "$BACT_DIR" "$META_DIR"

# Helper function for downloads
download() {
    local url="$1"
    local dest="$2"
    echo "  Downloading $(basename "$dest")..."
    curl -LfsS "$url" -o "$dest"
}

# Download common bacterial contaminant
echo "==> Downloading E. coli K-12 reference..."
download \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" \
  "$OUTDIR/ecoli_k12.fa.gz"

gunzip -f "$OUTDIR/ecoli_k12.fa.gz"

# Animal-specific contamination (minimal - just E. coli)
echo "==> Setting up animal contamination references..."
printf "ecoli_k12,%s\n" "$OUTDIR/ecoli_k12.fa" > "$OUTDIR/bbsplit_animal.txt"

echo "==> Animal contamination setup complete."
echo "    Output: $OUTDIR/bbsplit_animal.txt"


# Plant-specific contamination (comprehensive)
echo "==> Setting up plant contamination references..."

# Download UniVec_Core
echo "==> Downloading UniVec_Core (vector contamination)..."
download \
  "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core" \
  "$OUTDIR/univec_core.fa"

# Download assembly summaries
echo "==> Downloading assembly summary files..."
download \
  "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt" \
  "$META_DIR/fungi_assembly_summary.txt"

download \
  "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" \
  "$META_DIR/bacteria_assembly_summary.txt"

# Generate fungi manifest (FILTERED to common/important genera)
echo "==> Generating filtered fungi manifest (common plant pathogens)..."
awk -F '\t' '
BEGIN { OFS="\n" }
$1 !~ /^#/ && 
$11 == "latest" && 
$12 == "Complete Genome" &&
$20 != "na" &&
($8 ~ /^Saccharomyces([[:space:]]|$)/ ||    # Model organism/contaminant
 $8 ~ /^Candida([[:space:]]|$)/ ||          # Human/lab contaminant
 $8 ~ /^Aspergillus([[:space:]]|$)/ ||      # Common mold
 $8 ~ /^Fusarium([[:space:]]|$)/ ||         # Plant pathogen
 $8 ~ /^Botrytis([[:space:]]|$)/ ||         # Plant pathogen
 $8 ~ /^Magnaporthe([[:space:]]|$)/ ||      # Rice blast
 $8 ~ /^Ustilago([[:space:]]|$)/ ||         # Smut fungi
 $8 ~ /^Puccinia([[:space:]]|$)/) {         # Rust fungi
    path = $20
    sub(/\/$/, "", path)
    n = split(path, a, "/")
    asm = a[n]
    print path "/" asm "_genomic.fna.gz"
}' "$META_DIR/fungi_assembly_summary.txt" > "$META_DIR/fungi_fna_urls.txt"

# Generate plant pathogen bacteria manifest
echo "==> Generating plant pathogen bacteria manifest..."
awk -F '\t' '
$1 !~ /^#/ &&
$11 == "latest" &&
$12 == "Complete Genome" &&
($8 ~ /^Pseudomonas syringae([[:space:]]|$)/ || 
 $8 ~ /^Agrobacterium tumefaciens([[:space:]]|$)/ ||
 $8 ~ /^Xanthomonas([[:space:]]|$)/ ||
 $8 ~ /^Erwinia([[:space:]]|$)/) &&
$20 != "na" {
    path = $20
    sub(/\/$/, "", path)
    n = split(path, a, "/")
    asm = a[n]
    print path "/" asm "_genomic.fna.gz"
}' "$META_DIR/bacteria_assembly_summary.txt" > "$META_DIR/plant_pathogen_fna_urls.txt"

# Report counts
FUNGI_COUNT=$(wc -l < "$META_DIR/fungi_fna_urls.txt")
BACT_COUNT=$(wc -l < "$META_DIR/plant_pathogen_fna_urls.txt")
echo "    Found $FUNGI_COUNT fungal genomes"
echo "    Found $BACT_COUNT bacterial genomes"

# Download genomes
echo "==> Downloading fungal genomes..."
wget -nv -nc -P "$FUNGI_DIR" -i "$META_DIR/fungi_fna_urls.txt" || {
    echo "WARNING: Some fungal downloads failed, continuing..."
}

echo "==> Downloading bacterial plant pathogen genomes..."
wget -nv -nc -P "$BACT_DIR" -i "$META_DIR/plant_pathogen_fna_urls.txt" || {
    echo "WARNING: Some bacterial downloads failed, continuing..."
}

# Combine FASTA files
echo "==> Combining fungal FASTA files..."
find "$FUNGI_DIR" -name '*.fna.gz' -print0 | sort -z | xargs -0 zcat > "$OUTDIR/fungi_combined.fa"
FUNGI_SIZE=$(stat -f%z "$OUTDIR/fungi_combined.fa" 2>/dev/null || stat -c%s "$OUTDIR/fungi_combined.fa" 2>/dev/null || echo "unknown")
echo "    Combined fungi size: $FUNGI_SIZE bytes"

echo "==> Combining bacterial FASTA files..."
find "$BACT_DIR" -name '*.fna.gz' -print0 | sort -z | xargs -0 zcat > "$OUTDIR/bacteria_combined.fa"
BACT_SIZE=$(stat -f%z "$OUTDIR/bacteria_combined.fa" 2>/dev/null || stat -c%s "$OUTDIR/bacteria_combined.fa" 2>/dev/null || echo "unknown")
echo "    Combined bacteria size: $BACT_SIZE bytes"

# Generate BBSplit manifest 
echo "==> Generating BBSplit manifest..." 
{
  printf "ecoli,%s\n" "$OUTDIR/ecoli_k12.fa" 
  printf "fungi,%s\n" "$OUTDIR/fungi_combined.fa" 
  printf "bacteria,%s\n" "$OUTDIR/bacteria_combined.fa" 
  printf "univec,%s\n" "$OUTDIR/univec_core.fa" 
} > "$OUTDIR/bbsplit_plant.txt"

# Generate checksums
echo "==> Recording checksums..."
(
    cd "$OUTDIR"
    sha256sum \
        ecoli_k12.fa \
        univec_core.fa \
        fungi_combined.fa \
        bacteria_combined.fa \
        bbsplit_plant.txt > "$META_DIR/reference_checksums.sha256"
)

# Summary
echo ""
echo "==> Setup complete!"
echo "    Contamination type: $CONTAMINATION_TYPE"
echo "    Output directory: $OUTDIR"
echo "    BBSplit manifest: $OUTDIR/bbsplit_plant.txt"
echo ""
echo "==> Next steps for nf-core/rnaseq:"
echo "    1. Build BBSplit index:"
echo "       bbsplit.sh ref=fungi_combined.fa,bacteria_combined.fa,univec_core.fa \\"
echo "         path=$OUTDIR/bbsplit_index/ -Xmx64g"
echo ""
echo "    2. Use in nf-core/rnaseq pipeline:"
echo "       nextflow run nf-core/rnaseq \\"
echo "         --remove_ribo_rna true \\"
echo "         --bbsplit_fasta_list $OUTDIR/bbsplit_plant.txt \\"
echo "         ..."
echo ""
