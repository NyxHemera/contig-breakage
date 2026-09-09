#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from collections import OrderedDict
from pathlib import Path
from typing import IO
from typing import Iterator

def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Organize BLAST result rows into one file per locus tag."
	)
	parser.add_argument(
		"--blast_results_dir",
		default="Blast_Results",
		help="Directory containing BLAST result files (default: Blast_Results)",
	)
	parser.add_argument(
		"--output_dir",
		default="Organized_Blast_Results",
		help="Directory to save organized BLAST result files (default: Organized_Blast_Results)",
	)
	parser.add_argument(
		"--progress_every",
		type=int,
		default=100,
		help=(
			"Print progress every N files while processing (default: 100, set to 0 to disable)"
		),
	)
	parser.add_argument(
		"--max_open_files",
		type=int,
		default=256,
		help="Maximum number of output files to keep open at once (default: 256)",
	)
	return parser.parse_args()


def iter_blast_files(blast_results_dir: Path) -> Iterator[Path]:
	if not blast_results_dir.exists():
		raise FileNotFoundError(f"{blast_results_dir} was not found")

	if not blast_results_dir.is_dir():
		raise NotADirectoryError(f"{blast_results_dir} is not a directory")

	with os.scandir(blast_results_dir) as entries:
		for entry in entries:
			if entry.is_file():
				yield Path(entry.path)


def organize_blast_results(
	blast_results_dir: Path,
	output_dir: Path,
	progress_every: int,
	max_open_files: int,
) -> int:
	output_dir.mkdir(parents=True, exist_ok=True)
	if max_open_files < 1:
		raise ValueError("--max_open_files must be >= 1")

	output_handles: OrderedDict[str, IO[str]] = OrderedDict()
	output_writers: dict[str, csv.writer] = {}
	result_count = 0
	file_count = 0
	header = [
		"query_id",
		"locus_tag",
		"percent_identity",
		"alignment_length",
		"mismatches",
		"gap_opens",
		"query_start",
		"query_end",
		"subject_start",
		"subject_end",
		"evalue",
		"bit_score",
		"break_type",
		"break_position",
	]

	def get_writer(locus_tag: str) -> csv.writer:
		if locus_tag in output_writers:
			output_handles.move_to_end(locus_tag)
			return output_writers[locus_tag]

		if len(output_handles) >= max_open_files:
			old_locus_tag, old_handle = output_handles.popitem(last=False)
			old_handle.close()
			del output_writers[old_locus_tag]

		output_path = output_dir / f"{locus_tag}.tsv"
		file_exists = output_path.exists()
		is_empty = (not file_exists) or output_path.stat().st_size == 0
		handle = output_path.open("a", encoding="utf-8", newline="")
		writer = csv.writer(handle, delimiter="\t")

		if is_empty:
			writer.writerow(header)

		output_handles[locus_tag] = handle
		output_writers[locus_tag] = writer
		return writer

	try:
		for blast_file in iter_blast_files(blast_results_dir):
			file_count += 1
			with blast_file.open("r", encoding="utf-8") as handle:
				for raw_line in handle:
					line = raw_line.strip()
					if not line:
						continue

					columns = line.split()
					if len(columns) < 12:
						continue

					break_type = columns[0][-3:]
					if break_type not in {"beg", "end"}:
						continue

					locus_tag = columns[1]
					subject_start = columns[8]
					subject_end = columns[9]
					break_position = subject_start if break_type == "beg" else subject_end

					writer = get_writer(locus_tag)
					writer.writerow(columns + [break_type, break_position])
					result_count += 1

			if progress_every > 0 and file_count % progress_every == 0:
				print(f"Processed {file_count} files, mapped {result_count} rows...", flush=True)
	finally:
		for handle in output_handles.values():
			handle.close()

	return result_count


def main() -> None:
	args = parse_args()
	blast_results_dir = Path(args.blast_results_dir)
	output_dir = Path(args.output_dir)
	result_count = organize_blast_results(
		blast_results_dir,
		output_dir,
		args.progress_every,
		args.max_open_files,
	)
	print(f"Saved {result_count} rows into {output_dir}")


if __name__ == "__main__":
	main()
