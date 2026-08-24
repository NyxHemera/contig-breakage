#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys


def load_locus_counts(csv_path: Path, row_to_target: int) -> tuple[list[str], list[int]]:
	with csv_path.open("r", encoding="utf-8", newline="") as handle:
		rows = list(csv.reader(handle))

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


def create_histogram(
	loci: list[str],
	values: list[int],
	output_path: Path,
	title: str,
	y_label: str,
) -> None:
	figure_width = 18
	figure_height = 8
	fig, ax = plt.subplots(figsize=(figure_width, figure_height))

	positions = list(range(len(loci)))
	ax.bar(positions, values, color="#2b6cb0", edgecolor="#1a365d", linewidth=0.5)
	ax.set_xlabel("Locus tag")
	ax.set_ylabel(y_label)
	ax.set_title(title)

	label_step = max(1, len(loci) // 25)
	label_positions = positions[::label_step]
	if label_positions and label_positions[-1] != positions[-1]:
		label_positions.append(positions[-1])
	ax.set_xticks(label_positions)
	ax.set_xticklabels([loci[index] for index in label_positions], rotation=90, fontsize=6)
	ax.tick_params(axis="x", length=0)

	ax.grid(axis="y", linestyle="--", alpha=0.3)
	fig.tight_layout()
	fig.savefig(output_path, dpi=300, bbox_inches="tight")
	plt.close(fig)


def main() -> None:
	csv_path_str, row_to_target, title = sys.argv[1], sys.argv[2], sys.argv[3]
	csv_path = Path(csv_path_str)

	if not csv_path.exists():
		raise FileNotFoundError(f"{csv_path} was not found")
	
	loci, target_row = load_locus_counts(csv_path, int(row_to_target))

	filename = f"{title}_histogram.png"

	create_histogram(
		loci,
		target_row,
		Path(filename),
		f"{title} by locus tag",
		title,
	)

	print(f"Saved histogram to {filename}")


if __name__ == "__main__":
	main()
