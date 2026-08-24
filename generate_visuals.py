#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_locus_counts(csv_path: Path) -> tuple[list[str], list[int], list[int]]:
	with csv_path.open("r", encoding="utf-8", newline="") as handle:
		rows = list(csv.reader(handle))

	if len(rows) < 3:
		raise ValueError(f"Expected at least 3 rows in {csv_path}, found {len(rows)}")

	loci = [value.strip() for value in rows[0] if value.strip()]
	total_counts = [int(value) for value in rows[1] if value.strip()]

	unique_row_index = 3 if len(rows) > 3 else 2
	unique_counts = [int(value) for value in rows[unique_row_index] if value.strip()]

	if len(loci) != len(total_counts):
		raise ValueError(
			"Locus tag row and total-count row have different lengths: "
			f"{len(loci)} vs {len(total_counts)}"
		)

	if len(loci) != len(unique_counts):
		raise ValueError(
			"Locus tag row and unique-count row have different lengths: "
			f"{len(loci)} vs {len(unique_counts)}"
		)

	return loci, total_counts, unique_counts


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
	csv_path = Path("locus_counts.csv")

	if not csv_path.exists():
		raise FileNotFoundError(f"{csv_path} was not found")

	loci, total_counts, unique_counts = load_locus_counts(csv_path)

	create_histogram(
		loci,
		total_counts,
		Path("total_occurrences_histogram.png"),
		"Total occurrences across all files by locus tag",
		"Total occurrences across all files",
	)
	create_histogram(
		loci,
		unique_counts,
		Path("unique_file_occurrences_histogram.png"),
		"Unique file occurrences by locus tag",
		"Unique file occurrences",
	)

	print("Saved histogram to total_occurrences_histogram.png")
	print("Saved histogram to unique_file_occurrences_histogram.png")


if __name__ == "__main__":
	main()
