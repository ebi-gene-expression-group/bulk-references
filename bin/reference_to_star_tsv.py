#!/usr/bin/env python3

import argparse
import gzip
import math
import statistics
import sys
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def fasta_stats(fasta_path):
    lengths = []
    seq_len = 0

    with open_text(fasta_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_len > 0:
            lengths.append(seq_len)

    if not lengths:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")

    total_len = sum(lengths)
    n_contigs = len(lengths)

    lengths_sorted = sorted(lengths, reverse=True)
    running = 0
    n50 = 0
    half = total_len / 2
    for length in lengths_sorted:
        running += length
        if running >= half:
            n50 = length
            break

    return {
        "genome_size_bp": total_len,
        "n_contigs": n_contigs,
        "n50_bp": n50,
    }


def parse_gtf_attributes(attr_string):
    attrs = {}
    for field in attr_string.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def intron_stats_from_gtf(gtf_path):
    # transcript_id -> list of exon intervals
    transcripts = defaultdict(list)

    with open_text(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = parts

            if feature != "exon":
                continue

            a = parse_gtf_attributes(attrs)
            transcript_id = a.get("transcript_id")

            # Skip entries without transcript_id
            if transcript_id is None:
                continue

            start = int(start)
            end = int(end)
            if end < start:
                continue

            transcripts[transcript_id].append((chrom, strand, start, end))

    intron_lengths = []

    for transcript_id, exons in transcripts.items():
        if len(exons) < 2:
            continue

        # Keep only coordinates; assume transcript exons are on same chrom/strand
        exons_sorted = sorted(exons, key=lambda x: x[2])

        for i in range(len(exons_sorted) - 1):
            _, _, exon1_start, exon1_end = exons_sorted[i]
            _, _, exon2_start, exon2_end = exons_sorted[i + 1]

            intron_start = exon1_end + 1
            intron_end = exon2_start - 1
            intron_len = intron_end - intron_start + 1

            if intron_len > 0:
                intron_lengths.append(intron_len)

    if not intron_lengths:
        return {
            "max_intron_bp": 0,
            "p99_intron_bp": 0,
        }

    intron_lengths.sort()

    def percentile(values, p):
        if not values:
            return 0
        idx = max(0, min(len(values) - 1, math.ceil((p / 100.0) * len(values)) - 1))
        return values[idx]

    return {
        "max_intron_bp": intron_lengths[-1],
        "p99_intron_bp": percentile(intron_lengths, 99),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compute STAR helper TSV from genome FASTA and annotation GTF"
    )
    parser.add_argument("--dataset", required=True, help="Dataset/species name for TSV row")
    parser.add_argument("--fasta", required=True, help="Genome FASTA (.fa/.fasta, optionally .gz)")
    parser.add_argument("--gtf", required=True, help="Annotation GTF (optionally .gz)")
    parser.add_argument(
        "--out",
        default="-",
        help="Output TSV file path, or '-' for stdout"
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append row to existing TSV instead of overwriting"
    )
    args = parser.parse_args()

    fasta = fasta_stats(args.fasta)
    introns = intron_stats_from_gtf(args.gtf)

    header = [
        "dataset",
        "genome_size_bp",
        "n_contigs",
        "n50_bp",
        "max_intron_bp",
        "p99_intron_bp",
    ]

    row = [
        args.dataset,
        str(fasta["genome_size_bp"]),
        str(fasta["n_contigs"]),
        str(fasta["n50_bp"]),
        str(introns["max_intron_bp"]),
        str(introns["p99_intron_bp"]),
    ]

    output_line = "\t".join(row) + "\n"
    header_line = "\t".join(header) + "\n"

    if args.out == "-":
        sys.stdout.write(header_line)
        sys.stdout.write(output_line)
        return

    mode = "a" if args.append else "w"

    write_header = True
    if args.append:
        try:
            with open(args.out, "r") as fh:
                first = fh.readline()
                if first.strip():
                    write_header = False
        except FileNotFoundError:
            write_header = True

    with open(args.out, mode) as out_fh:
        if write_header:
            out_fh.write(header_line)
        out_fh.write(output_line)


if __name__ == "__main__":
    main()
