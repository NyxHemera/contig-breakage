# Contig Break Identification Standard Workflow

This repository contains a set of scripts for processing genome FASTA files, creating local BLAST databases, and comparing sequence fragments to identify matches. The main workflow starts with `extract_genome_sequences.py`, which prepares genome sequences for BLAST-based screening.

## Step 1 - Extract Genome Sequences:

`extract_genome_sequences.py` does the following:

1. Reads one or more raw genome files from a source directory.
2. Splits each FASTA record into a leading sequence segment and a trailing sequence segment.
3. Writes those segments to a processed FASTA directory.
4. Runs `blastn` against a local BLAST database.
5. Filters BLAST output to keep only the strongest hit for each query header.

### Command-line usage

Run the script from the project root:

```bash
python extract_genome_sequences.py
```

Example:

```bash
python extract_genome_sequences.py --processed_genomes_dir Processed_Genomes/ --blast_results_dir Blast_Results/ --log_file blast_logs.txt --skip_genome_processing
```

### Available arguments

```text
--genomes_dir
    Directory containing raw genome FASTA files
    Default: Genomes

--processed_genomes_dir
    Directory for processed FASTA outputs. Serves as output if --skip_genome_processing is disabled, and input if it is enabled.
    Default: Processed_Genomes

--blast_results_dir
    Directory to store BLAST result files
    Default: Blast_Results

--skip_genome_processing
    Skip the genome processing step and use files already present in processed_genomes_dir

--max_files
    Limit how many genome files are processed

--log_file
    Path for the script log
    Default: extract_genome_sequences.log
```

## Step 2 - Count Genes:

`count_genes.py` does the following:

1. Loads the BLAST result files from the `Blast_Results/` directory.
2. Extracts each locus tag from the second column of each BLAST row.
3. Counts how many times each locus appears across all files.
4. Counts how many unique files contain each locus at least once.
5. Reads the reference `.ffn` file and calculates GC counts and total nucleotide counts per locus.
6. Writes a summary CSV for downstream comparison and visualization.

### Command-line usage

Run the script from the project root:

```bash
python count_genes.py
```

### What the script reads

The script expects:

- a directory named `Blast_Results/` containing BLAST output files
- a reference genome FASTA file named `Pseudomonas_aeruginosa_PAO1_107.ffn`

It checks that both exist before processing.


### Output CSV

The script writes `locus_counts.csv` with the following rows:

1. locus tags
2. total occurrences across all files
3. unique-file occurrences
4. GC nucleotide counts from the reference `.ffn` file
5. total nucleotide counts from the reference `.ffn` file
6. GC percentage for each locus

## Step 3 - Organize Blast Results:

`organize_blast_results.py` does the following:

1. Scans the `Blast_Results/` directory for BLAST output files.
2. Reads each result row and keeps only rows associated with sequence break positions.
3. Groups rows by locus tag into separate TSV files.
4. Writes the organized output into a dedicated directory such as `Organized_Blast_Results/`.
5. Adds metadata about the original genome and the break position for downstream analysis.

### Command-line usage

Run the script from the project root:

```bash
python organize_blast_results.py --blast_results_dir Blast_Results/ --output_dir Organized_Blast_Results/ --progress_every 50 --max_open_files 128
```

### Available arguments

```text
--blast_results_dir
    Directory containing BLAST result files
    Default: Blast_Results

--output_dir
    Directory where organized locus-specific TSV files will be written
    Default: Organized_Blast_Results

--progress_every
    Print progress every N files while processing
    Default: 100

--max_open_files
    Maximum number of output files to keep open at once
    Default: 256

--max_files
    Maximum number of files to process from the sorted input directory
    Default: 0, which means process all files
```

### How it works

For each row, it reads values such as:

- `subject_length`
- `query_start`
- `query_end`
- `subject_start`
- `subject_end`

Then it determines a `break_position`:

- For a `beg` record, if `query_start == 1`, `break_position = subject_start`; otherwise it sets `break_position = -1`.
- For an `end` record, if `query_end == 250`, `break_position = subject_end`; otherwise it sets `break_position = -1`.

This is meant to record whether the sequence break landed at the expected location for the beginning or end of the extracted segment.

## Step 4 - Get Largest Values:

`get_largest_values.py` does the following:

1. Reads a CSV file and selects the row indicated by `--row_to_target`.
2. Pairs each locus tag with the corresponding value from that row.
3. Sorts those pairs by value in descending order.
4. Prints the top `--num_results` locus/value pairs.

### Command-line usage

Run the script from the project root:

```bash
python get_largest_values.py --csv_path locus_counts.csv --row_to_target 1 --num_results 10
```

### Available arguments

```text
--csv_path
    Path to the input CSV file

--row_to_target
    0-based row index containing values to analyze

--num_results
    Number of top values to return along with their corresponding locus tags
```

## Step 5 - Generate Visuals:

`generate_visuals.py` does the following:

1. Reads a CSV file and selects the row indicated by `--row_to_target`.
2. Validates that the first row contains locus tags and the target row has matching values.
3. Creates a histogram from the selected row across all loci.
4. Saves the plot as a PNG image with a name derived from the provided `--title` argument.
5. Supports either a linear histogram or a circular polar histogram for circular-genome views.

### Command-line usage

Run the script from the project root:

Example:

```bash
python generate_visuals.py --csv_path locus_counts.csv --row_to_target 1 --title TotalCounts
```

Circular example:

```bash
python generate_visuals.py --csv_path locus_counts.csv --row_to_target 1 --title UniqueCounts --circular
```

### Available arguments

```text
--csv_path
    Path to the input CSV file

--row_to_target
    0-based row index containing values to plot

--title
    Plot title and y-axis label prefix; the PNG is written as <title>_histogram.png

--circular
    Render a circular histogram instead of the default linear bar plot
```

## Step 6 - Check Break Consistency:

`check_break_consistency.py` can be run directly with `--locus_tag`, or it can read locus tags from standard input. This makes it easy to pipe the output of `get_largest_values.py` into a focused consistency check for the most important loci.

### Example pipeline

```bash
python get_largest_values.py --csv_path locus_counts_new.csv --row_to_target 1 --num_results 5 | python check_break_consistency.py --blast_results_dir Organized_Blast_Results/ --exclude_no_break
```

This command does the following:

1. `get_largest_values.py` reads `locus_counts_new.csv` and ranks the values in row 1.
2. It prints the top 5 loci as `locus,value` pairs.
3. `check_break_consistency.py` reads those lines from stdin.
4. It extracts the locus tag before the comma and checks the matching TSV files in `Organized_Blast_Results/`.
5. It excludes `break_position == -1` from the printed counts and histogram output when `--exclude_no_break` is set.

### Command-line usage

```bash
python check_break_consistency.py --blast_results_dir Organized_Blast_Results/ --locus_tag locus_1 locus_2 locus_3
```

or:

```bash
python check_break_consistency.py --blast_results_dir Organized_Blast_Results/ --exclude_no_break
```

with input via stdin:

```bash
python get_largest_values.py --csv_path locus_counts_new.csv --row_to_target 1 --num_results 5 | python check_break_consistency.py --blast_results_dir Organized_Blast_Results/ --exclude_no_break
```

### Available arguments

```text
--locus_tag
    One or more locus tags to inspect directly from the command line

--blast_results_dir
    Directory containing per-locus TSV files, for example Organized_Blast_Results

--exclude_no_break
    Exclude break_position = -1 from printed counts and histogram output
```

### What the script reads

The script expects one or more per-locus TSV files in the organized BLAST directory. It reads the `break_position` column from each file and counts how often each position occurs.

The script also accepts stdin lines in either of these forms:

```text
locus_tag
locus_tag,value
```

This is why the output from `get_largest_values.py` works cleanly in a pipe: the script strips everything after the first comma and keeps the locus identifier.

### How it works

For each selected locus:

- the script builds the path `Organized_Blast_Results/<locus_tag>.tsv`
- it opens the TSV and reads the `break_position` values
- it counts each position using a frequency counter
- it optionally removes `-1` entries when `--exclude_no_break` is set
- it generates a histogram PNG in `Break_Position_Histograms/`
- it prints the counts and the image output path to the terminal

### Output

The script prints break-position counts such as:

```text
Break position counts for PA14_0004
1	12
2	4
45	6
Histogram saved to Break_Position_Histograms/PA14_0004_break_position_histogram.png
```

It also creates a histogram image for each locus in the `Break_Position_Histograms/` directory.


# Contig Break Identification Alternative Workflow

## Step 1 - Create Genome Databases:

`create_genome_dbs.py` does the following:

1. Scans a folder of genome sequence files.
2. Builds a BLAST database for each genome in that folder.
3. Writes each database into its own subdirectory under the output folder.
4. Uses `makeblastdb` so the database can later be used by `blastn`.

### Command-line usage

```bash
python create_genome_dbs.py --genome_folder /storage/project/r-sbrown365-0/shared/00_BackupData/E_Mehlferber/P_aeruginosa_30K_Genome_Data/30K_Genomes_passedFilter --output_folder Genome_DBs
```

### Available arguments

```text
--genome_folder
    Folder containing genome sequence files

--output_folder
    Folder where per-genome BLAST database folders will be created

--dbtype
    BLAST DB type to create: nucl or prot
    Default: nucl

--limit
    Only build databases for the first N genome files in sorted order
    Default: all files
```

## Step 2 - Run Filtered Reference Against Genomes:

`run_filtered_against_genomes.py` does the following:

1. Reads a root directory containing BLAST database folders.
2. Finds each database prefix in the subdirectories under `--db_root_dir`.
3. Runs `blastn` using a single query genome file against each database.
4. Writes one output TSV per database folder into the specified output directory.
5. Limits the run to a configurable number of databases when `--max_databases` is provided.

### Command-line usage

```bash
python run_filtered_against_genomes.py --db_root_dir Genome_DBs/ --query_genome_file PAO1_filtered.ffn --output_dir Filtered_Blast_Run_Results/ --max_databases 20
```

### Available arguments

```text
--db_root_dir
    Directory containing one subfolder per genome database

--query_genome_file
    FASTA/FFN query genome file to compare against all databases

--output_dir
    Directory where BLAST result files will be written

--max_databases
    Optional cap on how many database folders to process
    Default: all folders
```

## Step 3 - Compare Lengths:

`compare_lengths.py` does the following:

1. Reads a directory of BLAST result files or a single BLAST output file.
2. Parses the reference `.ffn` file to determine expected gene lengths.
3. Groups BLAST rows by gene.
4. Classifies each gene result as `Complete`, `Incomplete`, `Missing`, or `Split` based on alignment length, gene length, and contig-break heuristics.
5. Saves the output files into `Compare_length_output/<gene>/Blast_gene_status/...`.

### Command-line usage

```bash
python compare_lengths.py Filtered_Blast_Run_Results/ PAO1_filtered.ffn
```

## Step 4 - Organize Filter Compare Lengths:

`organize_filter_compare_lengths_output.py` does the following:

1. Reads a gene-specific folder created by `compare_lengths.py`.
2. Looks in `Blast_gene_status/Split` for the split BLAST files to process.
3. Keeps only the highest-PID hit per contig within each file.
4. Builds the best sequential query-chain using a configurable nucleotide gap threshold.
5. Writes the filtered results into a new output directory under the gene name.

### Command-line usage

```bash
python3 organize_filter_compare_lengths_output.py --input_gene_folder Compare_length_output/PA4212/ --output_directory Filtered_Sequential_Compared_Lengths_Output/ --max_gap_nt 50
```

### Available arguments

```text
--input_gene_folder
    Gene folder containing the output generated by compare_lengths.py

--output_directory
    Root directory where filtered output files will be written

--max_files
    Optional limit on how many split files to process

--min_pid
    Optional minimum PID threshold before chain selection

--max_gap_nt
    Maximum allowed gap between adjacent rows in the query sequence chain
    Default: 10
```

### What the script reads

The script expects a gene folder with a `Blast_gene_status/Split` directory, which is produced by `compare_lengths.py`.

It reads each split-file TSV, removes redundant contig hits by keeping the highest-PID match per contig, and then selects the best sequential chain of query coordinates.

### How it works

For each split BLAST file:

- it loads the TSV into a pandas DataFrame
- it filters to the single highest-PID hit per contig
- it sorts rows by query coordinates
- it searches for the best chain where adjacent rows are separated by at most `--max_gap_nt`
- it writes the resulting filtered rows to the output directory under the gene name



