#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
	circular: bool = False,
) -> None:
	if circular:
		fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
		positions = list(range(len(loci)))
		n_loci = len(loci)
		theta = [2 * math.pi * idx / n_loci for idx in positions]
		bar_width = (2 * math.pi / n_loci) * 0.95

		ax.bar(theta, values, width=bar_width, bottom=0, color="#2b6cb0", edgecolor="#1a365d", linewidth=0.5)
		ax.set_theta_offset(math.pi / 2)
		ax.set_theta_direction(-1)
		ax.set_title(title, pad=24)
		ax.set_ylabel(y_label, labelpad=24)

		label_step = max(1, n_loci // 25)
		label_positions = positions[::label_step]
		if label_positions and label_positions[-1] != positions[-1]:
			label_positions.append(positions[-1])
		ax.set_xticks([theta[index] for index in label_positions])
		ax.set_xticklabels([loci[index] for index in label_positions], fontsize=6)
		ax.grid(linestyle="--", alpha=0.3)
	else:
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


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Generate locus histograms from a CSV count table.")
	parser.add_argument("--csv_path", required=True, help="Path to the input CSV file")
	parser.add_argument("--row_to_target", type=int, required=True, help="0-based row index containing values to plot")
	parser.add_argument("--title", required=True, help="Plot title and y-axis label prefix")
	parser.add_argument(
		"--circular",
		action="store_true",
		help="Render a circular histogram suitable for circular genomes",
	)
	return parser.parse_args()


def main() -> None:
	args = parse_args()
	csv_path = Path(args.csv_path)

	if not csv_path.exists():
		raise FileNotFoundError(f"{csv_path} was not found")
	
	loci, target_row = load_locus_counts(csv_path, args.row_to_target)

	filename = f"{args.title}_histogram.png"

	create_histogram(
		loci,
		target_row,
		Path(filename),
		f"{args.title} by locus tag",
		args.title,
		circular=args.circular,
	)

	print(f"Saved histogram to {filename}")


if __name__ == "__main__":
	main()
