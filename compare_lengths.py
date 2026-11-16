#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os, glob
import sys
import math
import re

if __name__ == "__main__":
    # length = sys.argv[1]
    # length = 10608
    blast_output = sys.argv[1]
    # blast_output = "comparelengthstestcases/testPA0041.fna.tsv"
    # gene = sys.argv[3]
    # gene = "PA0041"
    ffn_file = sys.argv[2]
    # ffn_file = "PAO1_filtered.ffn"
    # gene_proportion = sys.argv[4]
    gene_proportion = 0.05 
    # gene_difference = sys.argv[5]
    gene_difference = 0.4
    # tolerance = sys.argv[6]
    tolerance = 50
    # split_tolerance = sys.argv[7]
    split_tolerance = 5

def parse_ffn_genes(ffn_file):
    """Return a dict mapping gene names to their lengths from a .ffn file."""
    gene_lengths = {}
    current_gene = None
    current_seq = []

    with open(ffn_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_gene is not None:
                    gene_lengths[current_gene] = len("".join(current_seq).replace("\n", ""))
                current_gene = line[1:].split()[0]
                current_seq = []
            else:
                if current_gene is not None:
                    current_seq.append(line)

        if current_gene is not None:
            gene_lengths[current_gene] = len("".join(current_seq).replace("\n", ""))

    return gene_lengths

def get_gene_folders(gene, type_st):
    return f"Compare_length_output/{gene}/Blast_gene_status/{type_st}"


def split_blast_by_gene(blast_output):
    """Return a dict mapping each gene name to its corresponding BLAST rows."""
    blast_table = pd.read_csv(blast_output, sep="\t", names=['Gene', 'Prot_acc', 'Prot_acc_stripped', 'PID', 'nident', 'Gene_Length', 'Aln_length', 'evalue', 'slen', 'qstart', 'qend', 'Start', 'End', 'Strand', 'QStrand'])
    gene_tables = {}

    for gene_name in sorted(blast_table['Gene'].unique()):
        gene_tables[gene_name] = blast_table[blast_table['Gene'] == gene_name].copy()

    return gene_tables

def list_blast_files(blast_output):
    if os.path.isdir(blast_output):
        return sorted(os.path.join(blast_output, file) for file in os.listdir(blast_output) if os.path.isfile(os.path.join(blast_output, file)))
    if os.path.isfile(blast_output):
        return [blast_output]
    raise FileNotFoundError(f"BLAST output path not found: {blast_output}")

def make_folders():
    for gene in ffn_gene_lengths:
        complete=f"Compare_length_output/{gene}/Blast_gene_status/Complete"
        incomplete=f"Compare_length_output/{gene}/Blast_gene_status/Incomplete"
        missing=f"Compare_length_output/{gene}/Blast_gene_status/Missing"
        split=f"Compare_length_output/{gene}/Blast_gene_status/Split"

        for folder in [complete, incomplete, missing, split]:
            os.makedirs(folder, exist_ok=True)

def saving_funct(type_st, gene, table):
    folder = get_gene_folders(gene, type_st)
    # print(f"Saving to folder: {folder}")
    """Function to save blast output to defined folder"""
    table.to_csv(f"{folder}/{filename}", index = False, header = False, sep = "\t")

    with open(f"Compare_length_output/{gene}/{type_st}_output.txt", "a") as f:
        f.write(filename_no_ext_clean + "\n")

# # Print gene length and what gene proportion equates to
# print(length, gene_proportion)

ffn_gene_lengths = parse_ffn_genes(ffn_file)
# print(ffn_gene_lengths)

blast_files = list_blast_files(blast_output)
make_folders()

total_files = len(blast_files)
for completed_count, blast_file in enumerate(blast_files, start=1):
    filename = os.path.basename(blast_file)
    filename_no_ext = os.path.splitext(filename)[0]
    gene_tables = split_blast_by_gene(blast_file)

    for gene, table in gene_tables.items():
        filename_no_ext_clean = re.sub(f"_{gene}$", "", filename_no_ext)

        length = ffn_gene_lengths[gene]

        # If blast output is empty, assign to missing
        if len(table) == 0:
            saving_funct("Missing", gene, table)

        # If blast output yeilds one result, assign to complete
        elif len(table) == 1:
            saving_funct("Complete", gene, table)
        else:
            # If gene length in the first line matches alignment length in the first line, save to complete
            if table.iloc[0]['Gene_Length'] == table.iloc[0]['Aln_length']:
                saving_funct("Complete", gene, table)
            else:
                # Calculate length threshold
                length_threshold = math.floor(float(length) - (float(length) * float(gene_proportion)))

                # Calculate split threshold
                split_threshold = table.iloc[0]['slen'] - float(split_tolerance)

                # Establishes contig breaks (if a contig starts or ends in the middle of the gene of interest)
                if table.iloc[0]['Start'] == 1:
                    saving_funct("Split", gene, table)
                elif table.iloc[0]['End'] == 1:
                    saving_funct("Split", gene, table)
                elif table.iloc[0]['Start'] >= split_threshold:
                    saving_funct("Split", gene, table)
                elif table.iloc[0]['End'] >= split_threshold:
                    saving_funct("Split", gene, table)
                else:
                    # If alignment length of the first line is greater than or equal to the length threshold, assign to complete
                    if table.iloc[0]['Aln_length'] >= length_threshold:
                        saving_funct("Complete", gene, table)
                    else:
                        # Calculate gene difference threshold
                        gene_diff_threshold = math.floor(float(length) * float(gene_difference))

                        # If alignment length of the first line is less than the length threshold, assign to missing
                        if table.iloc[0]['Aln_length'] < gene_diff_threshold:
                            saving_funct("Missing", gene, table)
                        else:
                            if table.iloc[0]['qstart'] == 1:
                                # If the query start in the first line equals one
                                # Calculate the differences between query start and end of the first two lines
                                diff_qstart = abs(table.iloc[1]['qstart'] - table.iloc[0]['qend'])
                                diff_qend = abs(table.iloc[1]['qend'] - float(length))

                                # If the differences are less than or equal to the tolerances and the contig accessions match, assign as incomplete (biologically split)
                                # If the contig accessions do not match, assign to contig split
                                if (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] == table.iloc[1]['Prot_acc']:
                                    saving_funct("Incomplete", gene, table)
                                elif (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] != table.iloc[1]['Prot_acc']:
                                    saving_funct("Split", gene, table)
                                else:
                                    saving_funct("Complete", gene, table)
                            elif table.iloc[1]['qstart'] == 1:
                                # If the query start in the second line equals one
                                # Calculate the differences between query start and end of the first two lines
                                diff_qstart = abs(table.iloc[0]['qstart'] - table.iloc[1]['qend'])
                                diff_qend = abs(table.iloc[0]['qend'] - float(length))

                                # If the differences are less than or equal to the tolerances and the contig accessions match, assign as incomplete (biologically split)
                                # If the contig accessions do not match, assign to contig split
                                if (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] == table.iloc[1]['Prot_acc']:
                                    saving_funct("Incomplete", gene, table)
                                elif (diff_qstart <= float(tolerance)) and (diff_qend <= float(tolerance)) and table.iloc[0]['Prot_acc'] != table.iloc[1]['Prot_acc']:
                                    saving_funct("Split", gene, table)
                                else:
                                    saving_funct("Complete", gene, table)
                            else:
                                saving_funct("Complete", gene, table)

    print(f"Finished {filename}; fraction complete: {completed_count / total_files:.4f} ({completed_count}/{total_files})")