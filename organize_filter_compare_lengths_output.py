#!/usr/bin/env python3

import argparse
import os
from dataclasses import dataclass

import pandas as pd


DEFAULT_MAX_GAP_NT = 10


@dataclass
class BlastRow:
	idx: int
	q_lo: int
	q_hi: int
	pid: float


def chain_score(chain, rows, gene_length):
	first = rows[chain[0]]
	last = rows[chain[-1]]
	covered = sum(rows[i].q_hi - rows[i].q_lo + 1 for i in chain)
	boundary_gap = abs(first.q_lo - 1) + abs(gene_length - last.q_hi)
	return (
		covered,
		-boundary_gap,
		len(chain),
		-first.q_lo,
		last.q_hi,
	)


def keep_highest_pid_per_contig(df):
	# Column 2 is contig ID, column 4 is PID in the TSV schema.
	contig_col = df.iloc[:, 1]
	pid_col = pd.to_numeric(df.iloc[:, 3], errors="coerce")

	work = df.copy()
	work["_contig_id"] = contig_col
	work["_pid"] = pid_col

	# Stable sort keeps deterministic tie-breaking by original row order.
	work = work.sort_values(by="_pid", ascending=False, kind="mergesort")
	work = work.drop_duplicates(subset="_contig_id", keep="first")

	return work.drop(columns=["_contig_id", "_pid"])


def best_sequential_chain(df, min_pid=None, max_gap_nt=DEFAULT_MAX_GAP_NT):
	if min_pid is not None:
		df = df[df.iloc[:, 3].astype(float) >= min_pid]
	if df.empty:
		return []

	gene_length = int(float(df.iloc[0, 5]))
	rows = []
	for idx, row in df.iterrows():
		qstart = int(float(row.iloc[9]))
		qend = int(float(row.iloc[10]))
		q_lo = min(qstart, qend)
		q_hi = max(qstart, qend)
		pid = float(row.iloc[3])
		rows.append(BlastRow(idx=idx, q_lo=q_lo, q_hi=q_hi, pid=pid))

	rows.sort(key=lambda r: (r.q_lo, r.q_hi, -r.pid))
	n = len(rows)

	best_ending_at = [[i] for i in range(n)]

	for i in range(n):
		for j in range(i):
			# Allow up to a configurable gap between adjacent segments in query space.
			diff = rows[i].q_lo - rows[j].q_hi
			if 1 <= diff <= (max_gap_nt + 1):
				candidate = best_ending_at[j] + [i]
				if chain_score(candidate, rows, gene_length) > chain_score(best_ending_at[i], rows, gene_length):
					best_ending_at[i] = candidate

	best_chain = max(best_ending_at, key=lambda c: chain_score(c, rows, gene_length))
	ordered_original_indices = [rows[i].idx for i in best_chain]
	return ordered_original_indices


def process_split_files(input_gene_folder, output_root, max_files=None, min_pid=None, max_gap_nt=DEFAULT_MAX_GAP_NT):
	gene_name = os.path.basename(os.path.normpath(input_gene_folder))
	split_dir = os.path.join(input_gene_folder, "Blast_gene_status", "Split")

	if not os.path.isdir(split_dir):
		raise FileNotFoundError(f"Split directory not found: {split_dir}")

	output_gene_dir = os.path.join(output_root, gene_name)
	os.makedirs(output_gene_dir, exist_ok=True)

	files = sorted(
		f for f in os.listdir(split_dir)
		if os.path.isfile(os.path.join(split_dir, f))
	)
	if max_files is not None:
		files = files[:max_files]

	for idx, filename in enumerate(files, start=1):
		in_path = os.path.join(split_dir, filename)
		out_path = os.path.join(output_gene_dir, filename)

		df = pd.read_csv(in_path, sep="\t", header=None)
		df = keep_highest_pid_per_contig(df)
		chain_indices = best_sequential_chain(df, min_pid=min_pid, max_gap_nt=max_gap_nt)
		filtered_df = df.loc[chain_indices]

		# Keep chain order as computed from qstart/qend sequence.
		filtered_df.to_csv(out_path, sep="\t", header=False, index=False)

		completed_fraction = idx / len(files)
		completed_percent = completed_fraction * 100
		print(f"Processed {idx}/{len(files)} files ({completed_fraction:.2%}, {completed_percent:.1f}%)")


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Filter Split BLAST files to the best sequential qstart/qend chain "
			"and write results under output_directory/gene_name."
		)
	)
	parser.add_argument(
		"--input_gene_folder",
		type=str,
		required=True,
		help="Gene folder containing Blast_gene_status/Split files.",
	)
	parser.add_argument(
		"--output_directory",
		type=str,
		required=True,
		help="Root output directory for filtered Split files.",
	)
	parser.add_argument(
		"--max_files",
		type=int,
		default=None,
		help="Only process the first N Split files in sorted order (default: all files).",
	)
	parser.add_argument(
		"--min_pid",
		type=float,
		default=None,
		help="Only retain rows with PID >= X before sequential alignment (default: no PID filter).",
	)
	parser.add_argument(
		"--max_gap_nt",
		type=int,
		default=DEFAULT_MAX_GAP_NT,
		help="Maximum nucleotide gap allowed between adjacent rows in the query sequence chain (default: 10).",
	)
	return parser.parse_args()


def main():
	args = parse_args()

	if args.max_files is not None and args.max_files < 1:
		raise SystemExit("max_files must be a positive integer")
	if args.min_pid is not None and args.min_pid < 0:
		raise SystemExit("min_pid must be non-negative")
	if args.max_gap_nt < 0:
		raise SystemExit("max_gap_nt must be non-negative")

	process_split_files(
		args.input_gene_folder,
		args.output_directory,
		max_files=args.max_files,
		min_pid=args.min_pid,
		max_gap_nt=args.max_gap_nt,
	)


if __name__ == "__main__":
	main()
