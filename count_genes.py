#!/usr/bin/env python3
from pathlib import Path
import csv
from collections import defaultdict


def extract_locus_tag(header_line: str) -> str:
    """Extract locus_tag from an FFN/FASTA header line."""
    header = header_line[1:].strip() if header_line.startswith(">") else header_line.strip()

    for token in header.split(";"):
        token = token.strip()
        if token.startswith("locus_tag="):
            return token.split("=", 1)[1].strip()

    # Fallback to first whitespace-delimited token if locus_tag is absent.
    return header.split()[0] if header else ""


def load_sequence_stats_by_locus(ffn_path: Path) -> tuple[dict, dict]:
    """Return mappings of locus_tag -> GC count and locus_tag -> total nucleotide count."""
    gc_counts = {}
    total_nt_counts = {}
    current_locus = ""
    current_gc = 0
    current_total = 0

    with ffn_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_locus:
                    gc_counts[current_locus] = current_gc
                    total_nt_counts[current_locus] = current_total

                current_locus = extract_locus_tag(line)
                current_gc = 0
                current_total = 0
                continue

            sequence = line.upper()
            current_gc += sequence.count("G") + sequence.count("C")
            current_total += len(sequence)

    if current_locus:
        gc_counts[current_locus] = current_gc
        total_nt_counts[current_locus] = current_total

    return gc_counts, total_nt_counts


def main():
    blast_results_dir = Path("Blast_Results")
    ffn_path = Path("Pseudomonas_aeruginosa_PAO1_107.ffn")
    output_csv = Path("locus_counts.csv")

    if not blast_results_dir.exists():
        print(f"Error: {blast_results_dir} directory not found")
        return

    if not ffn_path.exists():
        print(f"Error: {ffn_path} file not found")
        return
    
    print(f"{blast_results_dir} exists: {blast_results_dir.exists()}")
    print(f"Loading sequence stats from {ffn_path}...")
    gc_counts, total_nt_counts = load_sequence_stats_by_locus(ffn_path)
    print(f"Loaded sequence stats for {len(gc_counts)} loci")

    results_files = [p for p in sorted(blast_results_dir.glob("*")) if p.is_file()]
    print(f"Processing {len(results_files)} file(s) from {blast_results_dir}")

    # Count every occurrence of each locus across all files.
    total_counts = defaultdict(int)
    # Count locus presence once per genome file.
    per_genome_counts = defaultdict(int)

    for file_path in results_files:
        filename = file_path.name
        print(f"Processing {filename}...")

        try:
            loci_in_this_file = set()
            with file_path.open("r", encoding="utf-8") as handle:
                for line in handle:
                    line = line.strip()
                    if not line:
                        continue

                    fields = line.split("\t")
                    if len(fields) < 2:
                        # Skip malformed lines that do not include a locus tag column.
                        continue

                    locus = fields[1].strip()
                    if not locus:
                        continue

                    total_counts[locus] += 1
                    loci_in_this_file.add(locus)

            for locus in loci_in_this_file:
                per_genome_counts[locus] += 1
        except Exception as e:
            print(f"  Error processing {filename}: {e}")

    loci = sorted(set(total_counts) | set(per_genome_counts))
    if not loci:
        print("No locus tags were found in the BLAST result files.")
        return

    with output_csv.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(loci)
        writer.writerow([total_counts[locus] for locus in loci])
        writer.writerow([per_genome_counts[locus] for locus in loci])
        writer.writerow([gc_counts.get(locus, 0) for locus in loci])
        writer.writerow([total_nt_counts.get(locus, 0) for locus in loci])
        writer.writerow([
            round((gc_counts.get(locus, 0) / total_nt_counts[locus]) * 100, 2)
            if total_nt_counts.get(locus, 0)
            else 0.0
            for locus in loci
        ])

    print(f"Wrote locus counts to {output_csv}")
    print("Row 1: locus tags")
    print("Row 2: total occurrences across all files")
    print("Row 3: unique-file occurrences (max 1 count per file)")
    print("Row 4: GC nucleotide counts from Pseudomonas_aeruginosa_PAO1_107.ffn")
    print("Row 5: total nucleotide counts from Pseudomonas_aeruginosa_PAO1_107.ffn")
    print("Row 6: GC percentage of total nucleotides")


if __name__ == "__main__":
    main()