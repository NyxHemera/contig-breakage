#!/usr/bin/env python3

import argparse
import shlex
import subprocess
import sys
from pathlib import Path


OUTFMT = (
	"6 qseqid sseqid sacc pident nident qlen length evalue slen "
	"qstart qend sstart send sstrand qstrand"
)


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Run blastn for one query genome against every BLAST database found "
			"in subdirectories of a database root directory."
		)
	)
	parser.add_argument(
		"--db_root_dir",
		required=True,
		type=Path,
		help="Directory containing one subfolder per genome database.",
	)
	parser.add_argument(
		"--query_genome_file",
		required=True,
		type=Path,
		help="FASTA/FFN query genome file to blast against all databases.",
	)
	parser.add_argument(
		"--output_dir",
		required=True,
		type=Path,
		help="Directory where BLAST result files will be written.",
	)
	parser.add_argument(
		"--max_databases",
		type=int,
		default=None,
		help=(
			"Optional cap on how many database folders to process. "
			"Defaults to all folders."
		),
	)
	return parser.parse_args()


def validate_paths(db_root_dir: Path, query_genome_file: Path) -> None:
	if not db_root_dir.is_dir():
		raise FileNotFoundError(f"Database root directory not found: {db_root_dir}")
	if not query_genome_file.is_file():
		raise FileNotFoundError(f"Query genome file not found: {query_genome_file}")


def find_db_prefix(db_folder: Path) -> Path:
	# BLAST nucleotide DBs include a .nin file for each database prefix.
	nin_files = sorted(db_folder.glob("*.nin"))
	if not nin_files:
		raise FileNotFoundError(f"No .nin database file found in {db_folder}")
	if len(nin_files) > 1:
		prefixes = ", ".join(str(p.with_suffix("")) for p in nin_files)
		raise RuntimeError(
			f"Multiple BLAST DB prefixes found in {db_folder}: {prefixes}"
		)
	return nin_files[0].with_suffix("")


def run_blastn(query_file: Path, db_prefix: Path, out_file: Path) -> None:
	command = [
		"blastn",
		"-query",
		str(query_file),
		"-db",
		str(db_prefix),
		"-outfmt",
		OUTFMT,
		"-out",
		str(out_file),
	]
	print(f"COMMAND: {shlex.join(command)}", flush=True)
	subprocess.run(command, check=True)


def main() -> int:
	args = parse_args()

	db_root_dir = args.db_root_dir.resolve()
	query_genome_file = args.query_genome_file.resolve()
	output_dir = args.output_dir.resolve()

	validate_paths(db_root_dir, query_genome_file)
	if args.max_databases is not None and args.max_databases <= 0:
		raise RuntimeError("--max_databases must be a positive integer")

	output_dir.mkdir(parents=True, exist_ok=True)

	db_folders = sorted(path for path in db_root_dir.iterdir() if path.is_dir())
	if not db_folders:
		raise RuntimeError(f"No database folders found in {db_root_dir}")

	if args.max_databases is not None:
		db_folders = db_folders[: args.max_databases]

	for db_folder in db_folders:
		db_prefix = find_db_prefix(db_folder)
		out_file = output_dir / f"{db_folder.name}.tsv"

		print(f"Running blastn for {db_folder.name} -> {out_file.name}")
		run_blastn(query_genome_file, db_prefix, out_file)

	print(f"Completed blastn runs for {len(db_folders)} database folders.")
	return 0


if __name__ == "__main__":
	try:
		raise SystemExit(main())
	except (FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as exc:
		print(f"Error: {exc}", file=sys.stderr)
		if isinstance(exc, subprocess.CalledProcessError):
			cmd = exc.cmd if isinstance(exc.cmd, list) else [str(exc.cmd)]
			print(f"FAILED_COMMAND: {shlex.join(cmd)}", file=sys.stderr)
		raise SystemExit(1)
