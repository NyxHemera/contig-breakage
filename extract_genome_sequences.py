#!/usr/bin/env python3
from pathlib import Path
import subprocess

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

def filter_blast_result_file(result_file_path, min_identity=70.0):
    best_rows = {}

    with open(result_file_path, 'r') as in_f:
        for raw_line in in_f:
            line = raw_line.rstrip('\n')
            if not line:
                continue

            columns = line.split('\t')
            if len(columns) < 3:
                continue

            header = columns[0]
            try:
                percent_identity = float(columns[2])
            except ValueError:
                continue

            if percent_identity < min_identity:
                continue

            existing = best_rows.get(header)
            # Keep the first row encountered when percent identities are equal.
            if existing is None or percent_identity > existing[0]:
                best_rows[header] = (percent_identity, line)

    with open(result_file_path, 'w') as out_f:
        for _, best_line in best_rows.values():
            out_f.write(f"{best_line}\n")

def filter_all_blast_results(blast_results_dir, min_identity=70.0):
    for result_file in sorted(blast_results_dir.glob('*.blastn.txt')):
        filter_blast_result_file(result_file, min_identity=min_identity)
        print(f"  Filtered BLAST results: {result_file}")

def run_blastn(processed_dir, blast_results_dir, blast_db_path='pseudodb'):
    print(f"\nRunning blastn on files in {processed_dir}...")
    blast_results_dir.mkdir(exist_ok=True)

    for fasta_file in processed_dir.glob('*.fasta'):
        output_file_name = f"{fasta_file.stem}.blastn.txt"
        output_file_path = blast_results_dir / output_file_name
        
        print(f"  Running blastn for {fasta_file.name}...")
        try:
            command = [
                'blastn',
                '-query', str(fasta_file),
                '-db', blast_db_path,
                '-out', str(output_file_path),
                '-outfmt', '6'
            ]
            
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            print(f"    blastn completed for {fasta_file.name}. Results saved to {output_file_path}")
            filter_blast_result_file(output_file_path)
            print(f"    Filtered to best hit per header in {output_file_path}")
            if result.stderr:
                print(f"    blastn stderr for {fasta_file.name}:\n{result.stderr}")

        except subprocess.CalledProcessError as e:
            print(f"  Error running blastn for {fasta_file.name}: {e}")
            print(f"    stdout: {e.stdout}")
            print(f"    stderr: {e.stderr}")
        except Exception as e:
            print(f"  An unexpected error occurred for {fasta_file.name}: {e}")

def main():
    genomes_dir = Path('Genomes')
    processed_dir = Path('Processed_Genomes')
    blast_results_dir = Path('Blast_Results')
    
    if not genomes_dir.exists():
        print(f"Error: {genomes_dir} directory not found")
        return

    processed_dir.mkdir(exist_ok=True)
    
    all_results = {}
    
    for file_path in sorted(genomes_dir.glob('*')):
        if file_path.is_file():
            filename = file_path.name
            print(f"Processing {filename}...")
            
            try:
                sequences = process_genome_file(file_path)
                all_results[filename] = sequences
                print(f"  Found {len(sequences)} sequence sections")

                output_filename = f"{file_path.stem}.fasta"
                output_path = processed_dir / output_filename
                write_processed_fasta(output_path, sequences)
                print(f"  Wrote processed FASTA: {output_path}")
            except Exception as e:
                print(f"  Error processing {filename}: {e}")

    run_blastn(processed_dir, blast_results_dir)
    print("\nFiltering all existing BLAST result files...")
    filter_all_blast_results(blast_results_dir)


if __name__ == '__main__':
    main()