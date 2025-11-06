#!/usr/bin/env bash
# Usage: ./format_genome_entries.sh input_file.txt ENSEMBL_RELNO ENSEMBLGENOMES_RELNO
# Example: ./format_genome_entries.sh genome_references.conf 112 60
#
# Reads each line, skips comments or empty lines, replaces RELNO with proper release numbers,
# and outputs YAML for each entry.
# This script converts atlas-annotations/sh/gtf/genome_references.conf to YAML for nf-core/references piepline 

input_file="$1"
ensembl_rel="$2"
ensemblgenomes_rel="$3"

if [[ -z "$input_file" || -z "$ensembl_rel" || -z "$ensemblgenomes_rel" ]]; then
  echo "Usage: $0 input_file.txt ENSEMBL_RELNO ENSEMBLGENOMES_RELNO"
  exit 1
fi

VERIFY_GZIP="${VERIFY_GZIP:-0}"

url_exists () {
  # 1) HEAD request, follow redirects; ensure 2xx, non-empty content-length if present
  local u="$1"

  # -I: HEAD, -L: follow redirects, -s: silent, -S: show errors, --fail: non-2xx => error
  if ! curl -I -L -sS --fail --retry 3 --retry-connrefused --max-time 30 "$u" > /tmp/headers.$$ 2>/dev/null; then
    rm -f /tmp/headers.$$
    return 1
  fi
  # Optionally ensure not zero-length if server returns Content-Length
  if grep -i '^Content-Length:' /tmp/headers.$$ >/dev/null 2>&1; then
    local sz
    sz=$(awk 'BEGIN{IGNORECASE=1} /^Content-Length:/ {print $2}' /tmp/headers.$$ | tr -d '\r')
    if [[ -n "$sz" && "$sz" -eq 0 ]]; then
      rm -f /tmp/headers.$$
      return 1
    fi
  fi
  rm -f /tmp/headers.$$
  return 0
}

gzip_ok () {
  # 2) Stream a small byte-range and ask gzip to test; saves bandwidth vs full file
  # Many servers support Range; if not, curl falls back to full – still works, but slower.
  local u="$1"
  # Try first 1MB; adjust if needed
  if curl -L -sS --fail --retry 2 --retry-connrefused --max-time 120 -r 0-1048575 "$u" | gzip -t 2>/dev/null; then
    return 0
  else
    return 1
  fi
}

validate_three () {
  # args: fasta_url gtf_url cdna_url
  local fasta="$1" gtf="$2" cdna="$3"
  local ok=0

  for label in fasta gtf cdna; do
    ok=0
    local u
    eval u=\$$label
    [[ -z "$u" ]] && continue
    if ! url_exists "$u"; then
      echo "Warning: $label URL not reachable: $u" >&2
      ok=1
      continue
    fi
    if [[ "$VERIFY_GZIP" == "1" ]]; then
      if ! gzip_ok "$u"; then
        echo "Warning: $label URL is reachable but gzip test failed: $u" >&2
        ok=1
      fi
    fi
  done

  return $ok
}

while IFS= read -r line; do
  # Skip empty lines or comments
  [[ -z "$line" || "$line" =~ ^# ]] && continue

  # Split into fields
  species=$(echo "$line" | awk '{print $1}')
  taxid=$(echo "$line" | awk '{print $2}')
  source=$(echo "$line" | awk '{print $3}')
  fasta_url=$(echo "$line" | awk '{print $4}')
  cdna_url=$(echo "$line" | awk '{print $5}')
  gtf_url=$(echo "$line" | awk '{print $6}')
  genome=$(echo "$line" | awk '{print $7}')

  # Determine which release to use based on FTP domain
  if [[ "$fasta_url" == *"ftp.ensemblgenomes.org"* ]]; then
    rel="$ensemblgenomes_rel"
  elif [[ "$fasta_url" == *"ftp.ensembl.org"* ]]; then
    rel="$ensembl_rel"
  else
    echo "Warning: Unknown FTP source for $species — skipping."
    continue
  fi

  source_name="Ensembl"
  # convert ftp to https and replace RELNO
  fasta_url_mod=$(echo "$fasta_url" | sed "s|ftp://|http://|" | sed "s|RELNO|$rel|g")
  gtf_url_mod=$(echo "$gtf_url" | sed "s|ftp://|http://|" | sed "s|RELNO|$rel|g")
  cdna_url_mod=$(echo "$cdna_url" | sed "s|ftp://|http://|" | sed "s|RELNO|$rel|g")

  # Validate URLs before emitting YAML
  if validate_three "$fasta_url_mod" "$gtf_url_mod" "$cdna_url_mod"; then
    : # all good
  else
    echo "Skipping ${species} ($genome) due to URL validation failures." >&2
    continue
  fi

  

  # Capitalize species nicely (e.g. saccharomyces_cerevisiae → Saccharomyces_cerevisiae)
  species_cap=$(echo "$species" | awk -F'_' '{print toupper(substr($1,1,1)) substr($1,2) "_" toupper(substr($2,1,1)) substr($2,2)}')

  # Print YAML entry
  cat <<EOF
- genome: $genome
  fasta: "$fasta_url_mod"
  gtf: "$gtf_url_mod"
  source_version: "${source_name}_${rel}"
  species: "$species_cap"
  source: "$source_name"
  # Add these fields to ensure index generation
  mito_name: "MT"
  macs_gsize: 12100000
EOF

done < "$input_file"
