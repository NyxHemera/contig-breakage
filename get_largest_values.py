#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_locus_and_target_values(csv_path: Path, row_to_target: int) -> tuple[list[str], list[int]]:
	with csv_path.open("r", encoding="utf-8", newline="") as handle:
		rows = list(csv.reader(handle))

	if not rows:
		raise ValueError(f"{csv_path} is empty")

	if len(rows) < row_to_target + 1:
		raise ValueError(f"Expected at least {row_to_target + 1} rows in {csv_path}, found {len(rows)}")

	loci = [value.strip() for value in rows[0] if value.strip()]
	target_row = [int(value) for value in rows[row_to_target] if value.strip()]

	if len(loci) != len(target_row):
		raise ValueError(
			"Locus tag row and target row have different lengths: "
			f"{len(loci)} vs {len(target_row)}"
		)

	return loci, target_row


def get_largest_pairs(loci: list[str], values: list[int], num_results: int) -> list[tuple[str, int]]:
	if num_results < 1:
		raise ValueError("num_results must be at least 1")

	paired = list(zip(loci, values))
	sorted_pairs = sorted(paired, key=lambda item: item[1], reverse=True)
	return sorted_pairs[:num_results]


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Return largest values from a target CSV row with their corresponding locus tags."
	)
	parser.add_argument("--csv_path", required=True, help="Path to the input CSV file")
	parser.add_argument(
		"--row_to_target",
		type=int,
		required=True,
		help="0-based row index containing values to analyze",
	)
	parser.add_argument(
		"--num_results",
		type=int,
		required=True,
		help="Number of top values to return",
	)
	return parser.parse_args()


def main() -> None:
	args = parse_args()
	csv_path = Path(args.csv_path)

	if not csv_path.exists():
		raise FileNotFoundError(f"{csv_path} was not found")

	loci, target_values = load_locus_and_target_values(csv_path, args.row_to_target)
	top_results = get_largest_pairs(loci, target_values, args.num_results)

	for locus, value in top_results:
		print(f"{locus},{value}")


if __name__ == "__main__":
	main()
