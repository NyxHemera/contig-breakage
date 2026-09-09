#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import subprocess
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Process genome FASTA files, run blastn, and filter BLAST results."
        )
    )
    parser.add_argument(
        '--genomes_dir',
        default='Genomes',
        help='Directory containing raw genome files (default: Genomes)',
    )
    parser.add_argument(
        '--processed_genomes_dir',
        default='Processed_Genomes',
        help='Directory containing processed FASTA files (default: Processed_Genomes)',
    )
    parser.add_argument(
        '--blast_results_dir',
        default='Blast_Results',
        help='Directory to write BLAST results (default: Blast_Results)',
    )
    parser.add_argument(
        '--skip_genome_processing',
        action='store_true',
        help='Skip processing genome files and use existing files in processed_genomes_dir',
    )
    parser.add_argument(
        '--max_files',
        type=int,
        default=None,
        help='Optional maximum number of files to process (default: process all)',
    )
    parser.add_argument(
        '--log_file',
        default='extract_genome_sequences.log',
        help='Path to log file for status/warning/error output (default: extract_genome_sequences.log)',
    )
    return parser.parse_args()


def log_message(message, log_file_path=None):
    print(message)
    if log_file_path is None:
        return

    log_file_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().isoformat(timespec='seconds')
    with open(log_file_path, 'a') as log_f:
        log_f.write(f"[{timestamp}] {message}\n")


def iter_files(directory):
    if not directory.exists():
        raise FileNotFoundError(f"{directory} was not found")

    if not directory.is_dir():
        raise NotADirectoryError(f"{directory} is not a directory")

    with os.scandir(directory) as entries:
        for entry in entries:
            if entry.is_file():
                yield Path(entry.path)


def iter_processed_fasta_files(processed_dir):
    for file_path in iter_files(processed_dir):
        if file_path.suffix.lower() == '.fasta':
            yield file_path

def extract_sequence_ends(sequence, length=250):
    if len(sequence) == 0:
        return {'first_250': '', 'last_250': ''}
    elif len(sequence) <= length * 2:
        # If sequence is shorter than 500 chars, return the whole thing
        return {'first_250': sequence, 'last_250': sequence}
    else:
        return {
            'first_250': sequence[:length],
            'last_250': sequence[-length:]
        }

def process_genome_file(file_path):
    sequences = {}
    
    with open(file_path, 'r') as f:
        current_header = None
        current_sequence = ""
        
        for line in f:
            line = line.rstrip('\n')
            
            if line.startswith('>'):
                # Save the previous sequence
                if current_header is not None:
                    sequences[current_header] = extract_sequence_ends(current_sequence)
                
                # Start a new sequence
                current_header = line
                current_sequence = ""
            else:
                # Append to current sequence
                current_sequence += line.strip()
        
        # Save the last sequence
        if current_header is not None:
            sequences[current_header] = extract_sequence_ends(current_sequence)
    
    return sequences

def write_processed_fasta(output_file_path, sequences):
    with open(output_file_path, 'w') as out_f:
        for header, ends in sequences.items():
            clean_header = header[1:] if header.startswith('>') else header
            header_parts = clean_header.split(' ', 1)
            header_prefix = header_parts[0]
            header_suffix = f" {header_parts[1]}" if len(header_parts) > 1 else ''

            out_f.write(f">{header_prefix}_beg{header_suffix}\n")
            out_f.write(f"{ends['first_250']}\n")

            out_f.write(f">{header_prefix}_end{header_suffix}\n")
            out_f.write(f"{ends['last_250']}\n")

def filter_blast_result_file(result_file_path, min_identity=70.0, log_file_path=None):
    best_rows = {}
    expected_columns = 14
    pending_fragment = None
    malformed_row_count = 0
    repaired_row_count = 0

    with open(result_file_path, 'r') as in_f:
        for raw_line in in_f:
            line = raw_line.rstrip('\n')
            if not line:
                continue

            if pending_fragment is not None:
                # Some files occasionally contain a line split mid-record; join once.
                line = pending_fragment + line
                pending_fragment = None
                repaired_row_count += 1

            columns = line.split('\t')
            if len(columns) < expected_columns:
                pending_fragment = line
                continue

            if len(columns) != expected_columns:
                malformed_row_count += 1
                continue

            header = columns[0]
            try:
                percent_identity = float(columns[3])
            except ValueError:
                continue

            if percent_identity < min_identity:
                continue

            existing = best_rows.get(header)
            # Keep the first row encountered when percent identities are equal.
            if existing is None or percent_identity > existing[0]:
                best_rows[header] = (percent_identity, line)

    if pending_fragment is not None:
        malformed_row_count += 1

    if repaired_row_count > 0:
        log_message(
            f"    Repaired {repaired_row_count} split BLAST row(s) in {result_file_path}",
            log_file_path=log_file_path,
        )

    if malformed_row_count > 0:
        log_message(
            f"    Skipped {malformed_row_count} malformed BLAST row(s) in {result_file_path}",
            log_file_path=log_file_path,
        )

    with open(result_file_path, 'w') as out_f:
        for _, best_line in best_rows.values():
            out_f.write(f"{best_line}\n")

def filter_all_blast_results(blast_results_dir, min_identity=70.0, log_file_path=None):
    for result_file in sorted(iter_files(blast_results_dir)):
        if result_file.suffixes[-2:] != ['.blastn', '.txt']:
            continue
        filter_blast_result_file(
            result_file,
            min_identity=min_identity,
            log_file_path=log_file_path,
        )
        log_message(f"  Filtered BLAST results: {result_file}", log_file_path=log_file_path)

def run_blastn(
    processed_dir,
    blast_results_dir,
    blast_db_path='pseudodb',
    max_files=None,
    log_file_path=None,
):
    log_message(f"\nRunning blastn on files in {processed_dir}...", log_file_path=log_file_path)
    blast_results_dir.mkdir(parents=True, exist_ok=True)

    processed_count = 0

    for fasta_file in sorted(iter_processed_fasta_files(processed_dir)):
        if max_files is not None and processed_count >= max_files:
            break

        processed_count += 1
        output_file_name = f"{fasta_file.stem}.blastn.txt"
        output_file_path = blast_results_dir / output_file_name
        
        log_message(f"  Running blastn for {fasta_file.name}...", log_file_path=log_file_path)
        try:
            command = [
                'blastn',
                '-query', str(fasta_file),
                '-db', blast_db_path,
                '-out', str(output_file_path),
                '-outfmt', '6 qseqid sseqid sacc pident nident qlen length evalue slen qstart qend sstart send sstrand qstrand'
            ]

            log_message(str(command), log_file_path=log_file_path)
            
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            log_message(
                f"    blastn completed for {fasta_file.name}. Results saved to {output_file_path}",
                log_file_path=log_file_path,
            )
            filter_blast_result_file(output_file_path, log_file_path=log_file_path)
            log_message(
                f"    Filtered to best hit per header in {output_file_path}",
                log_file_path=log_file_path,
            )
            if result.stderr:
                log_message(
                    f"    blastn stderr for {fasta_file.name}:\n{result.stderr}",
                    log_file_path=log_file_path,
                )

        except subprocess.CalledProcessError as e:
            log_message(
                f"  Error running blastn for {fasta_file.name}: {e}",
                log_file_path=log_file_path,
            )
            log_message(f"    stdout: {e.stdout}", log_file_path=log_file_path)
            log_message(f"    stderr: {e.stderr}", log_file_path=log_file_path)
        except Exception as e:
            log_message(
                f"  An unexpected error occurred for {fasta_file.name}: {e}",
                log_file_path=log_file_path,
            )

def main():
    args = parse_args()
    genomes_dir = Path(args.genomes_dir)
    processed_dir = Path(args.processed_genomes_dir)
    blast_results_dir = Path(args.blast_results_dir)
    max_files = args.max_files
    log_file_path = Path(args.log_file) if args.log_file else None

    if max_files is not None and max_files < 1:
        log_message("Error: --max_files must be >= 1", log_file_path=log_file_path)
        return

    if not args.skip_genome_processing:
        if not genomes_dir.exists():
            log_message(
                f"Error: {genomes_dir} directory not found",
                log_file_path=log_file_path,
            )
            return

        processed_dir.mkdir(parents=True, exist_ok=True)

        processed_count = 0
        for file_path in sorted(iter_files(genomes_dir)):
            if max_files is not None and processed_count >= max_files:
                break

            filename = file_path.name
            log_message(f"Processing {filename}...", log_file_path=log_file_path)

            try:
                sequences = process_genome_file(file_path)
                log_message(
                    f"  Found {len(sequences)} sequence sections",
                    log_file_path=log_file_path,
                )

                output_filename = f"{file_path.stem}.fasta"
                output_path = processed_dir / output_filename
                write_processed_fasta(output_path, sequences)
                log_message(
                    f"  Wrote processed FASTA: {output_path}",
                    log_file_path=log_file_path,
                )
                processed_count += 1
            except Exception as e:
                log_message(
                    f"  Error processing {filename}: {e}",
                    log_file_path=log_file_path,
                )
    else:
        if not processed_dir.exists() or not processed_dir.is_dir():
            log_message(
                f"Error: processed genomes directory not found or invalid: {processed_dir}",
                log_file_path=log_file_path,
            )
            return

    run_blastn(
        processed_dir,
        blast_results_dir,
        max_files=max_files,
        log_file_path=log_file_path,
    )
    # print("\nFiltering all existing BLAST result files...")
    # filter_all_blast_results(blast_results_dir, log_file_path=log_file_path)


if __name__ == '__main__':
    main()