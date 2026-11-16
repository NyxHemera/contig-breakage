#!/usr/bin/env python3
"""Create BLAST databases for each genome file in a folder.

Usage:
	python create_genome_dbs.py /path/to/genomes /path/to/output
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Build a BLAST database for each file in a genome folder, "
			"writing each database to its own output subfolder."
		)
	)
	parser.add_argument(
		"--genome_folder",
		type=Path,
		help="Folder containing genome sequence files.",
	)
	parser.add_argument(
		"--output_folder",
		type=Path,
		help="Folder where per-genome BLAST database folders will be created.",
	)
	parser.add_argument(
		"--dbtype",
		choices=("nucl", "prot"),
		default="nucl",
		help="BLAST DB type to create (default: nucl).",
	)
	parser.add_argument(
		"--limit",
		type=int,
		default=None,
		help="Only build databases for the first N genome files in sorted order (default: all files).",
	)
	return parser.parse_args()


def build_db_for_file(genome_file: Path, output_root: Path, dbtype: str) -> None:
	per_genome_out = output_root / genome_file.name
	per_genome_out.mkdir(parents=True, exist_ok=True)

	db_prefix = per_genome_out / "db"
	cmd = [
		"makeblastdb",
		"-in",
		str(genome_file),
		"-dbtype",
		dbtype,
		"-out",
		str(db_prefix),
	]

	print(f"[RUN] {' '.join(cmd)}")
	subprocess.run(cmd, check=True)
	print(f"[OK] Created DB for {genome_file.name} in {per_genome_out}")


def main() -> int:
	args = parse_args()

	if shutil.which("makeblastdb") is None:
		print("Error: 'makeblastdb' was not found in PATH.", file=sys.stderr)
		return 1

	genome_folder = args.genome_folder.resolve()
	output_folder = args.output_folder.resolve()

	if not genome_folder.exists() or not genome_folder.is_dir():
		print(f"Error: genome_folder does not exist or is not a directory: {genome_folder}", file=sys.stderr)
		return 1

	output_folder.mkdir(parents=True, exist_ok=True)

	genome_files = sorted(path for path in genome_folder.iterdir() if path.is_file())
	if args.limit is not None:
		if args.limit < 1:
			print(f"Error: --limit must be a positive integer, got {args.limit}.", file=sys.stderr)
			return 2
		genome_files = genome_files[: args.limit]
		print(f"Building databases for the first {len(genome_files)} genome file(s) only.")

	if not genome_files:
		print(f"No files found in genome_folder: {genome_folder}")
		return 0

	failures = 0
	for genome_file in genome_files:
		try:
			build_db_for_file(genome_file, output_folder, args.dbtype)
		except subprocess.CalledProcessError as exc:
			failures += 1
			print(f"[FAIL] makeblastdb failed for {genome_file}: exit code {exc.returncode}", file=sys.stderr)

	if failures:
		print(f"Completed with {failures} failure(s).", file=sys.stderr)
		return 1

	print("Completed successfully.")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
