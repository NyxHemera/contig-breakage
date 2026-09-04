#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt


HISTOGRAM_DIR = Path("Break_Position_Histograms")


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Count break_position values for one or more locus tags and generate histograms."
		)
	)
	parser.add_argument(
		"--locus_tag",
		nargs="+",
		help=(
			"One or more locus tags to inspect (maps to <locus_tag>.tsv). "
			"If omitted, locus tags are read from stdin."
		),
	)
	parser.add_argument(
		"--blast_results_dir",
		required=True,
		help="Directory containing per-locus TSV files (for example Organized_Blast_Results)",
	)
	return parser.parse_args()


def parse_locus_from_line(raw_line: str) -> str | None:
	line = raw_line.strip()
	if not line:
		return None
	if "," in line:
		return line.split(",", 1)[0].strip() or None
	return line


def collect_locus_tags(cli_locus_tags: list[str] | None) -> list[str]:
	collected: list[str] = []
	seen: set[str] = set()

	if cli_locus_tags:
		for locus_tag in cli_locus_tags:
			cleaned = locus_tag.strip()
			if cleaned and cleaned not in seen:
				seen.add(cleaned)
				collected.append(cleaned)

	if not sys.stdin.isatty():
		for raw_line in sys.stdin:
			parsed = parse_locus_from_line(raw_line)
			if parsed and parsed not in seen:
				seen.add(parsed)
				collected.append(parsed)

	return collected


def parse_break_position(value: str) -> int | str:
	try:
		return int(value)
	except ValueError:
		return value


def break_sort_key(value: int | str) -> tuple[int, int | str]:
	if isinstance(value, int):
		return (0, value)
	return (1, value)


def count_break_positions(tsv_path: Path) -> Counter[int | str]:
	with tsv_path.open("r", encoding="utf-8", newline="") as handle:
		reader = csv.DictReader(handle, delimiter="\t")
		if not reader.fieldnames or "break_position" not in reader.fieldnames:
			raise ValueError(f"Missing break_position column in {tsv_path}")

		counts: Counter[int | str] = Counter()
		try:
			for row in reader:
				raw_break_position = (row.get("break_position") or "").strip()
				if not raw_break_position:
					continue
				counts[parse_break_position(raw_break_position)] += 1
		except csv.Error as err:
			line_num = reader.line_num
			raise ValueError(
				f"Malformed TSV data in {tsv_path} at line {line_num}: {err}"
			) from err

	return counts


def plot_histogram(locus_tag: str, counts: Counter[int | str], output_dir: Path) -> Path:
	output_dir.mkdir(parents=True, exist_ok=True)

	ordered_positions = sorted(counts, key=break_sort_key)
	x_labels = [str(position) for position in ordered_positions]
	y_values = [counts[position] for position in ordered_positions]
	x_positions = list(range(len(x_labels)))

	fig_width = max(8, min(30, len(x_labels) * 0.35))
	fig, ax = plt.subplots(figsize=(fig_width, 6))
	ax.bar(x_positions, y_values)
	ax.set_xlim(-0.5, len(x_positions) - 0.5)
	ax.margins(x=0)
	ax.set_xlabel("break_position")
	ax.set_ylabel("count")
	ax.set_title(f"Break position frequency for {locus_tag}")

	if len(x_labels) <= 25:
		ax.set_xticks(x_positions)
		ax.set_xticklabels(x_labels)
		if len(x_labels) > 12:
			ax.tick_params(axis="x", labelrotation=90)
	else:
		max_tick_labels = 20
		step = math.ceil(len(x_labels) / max_tick_labels)
		tick_indices = list(range(0, len(x_labels), step))
		if tick_indices[-1] != len(x_labels) - 1:
			tick_indices.append(len(x_labels) - 1)
		ax.set_xticks(tick_indices)
		ax.set_xticklabels([x_labels[index] for index in tick_indices], rotation=90)

	plt.tight_layout()
	histogram_path = output_dir / f"{locus_tag}_break_position_histogram.png"
	fig.savefig(histogram_path, dpi=200)
	plt.close(fig)

	return histogram_path


def print_counts(locus_tag: str, counts: Counter[int | str]) -> None:
	print(f"Break position counts for {locus_tag}")
	# for break_position in sorted(counts, key=break_sort_key):
	# 	print(f"{break_position}\t{counts[break_position]}")


def main() -> None:
	args = parse_args()
	blast_results_dir = Path(args.blast_results_dir)
	locus_tags = collect_locus_tags(args.locus_tag)

	if not blast_results_dir.exists():
		raise FileNotFoundError(f"{blast_results_dir} was not found")

	if not blast_results_dir.is_dir():
		raise NotADirectoryError(f"{blast_results_dir} is not a directory")

	if not locus_tags:
		raise ValueError(
			"No locus tags provided. Use --locus_tag or pipe tags via stdin (locus or locus,count)."
		)

	missing_locus_tags: list[str] = []
	data_error_locus_tags: list[str] = []

	for index, locus_tag in enumerate(locus_tags):
		tsv_path = blast_results_dir / f"{locus_tag}.tsv"
		if not tsv_path.exists():
			missing_locus_tags.append(locus_tag)
			print(f"No TSV found for locus tag {locus_tag}: {tsv_path}", file=sys.stderr)
			continue

		try:
			counts = count_break_positions(tsv_path)
		except ValueError as err:
			data_error_locus_tags.append(locus_tag)
			print(str(err), file=sys.stderr)
			continue

		if not counts:
			print(f"No break_position values found in {tsv_path}")
			continue

		histogram_path = plot_histogram(locus_tag, counts, HISTOGRAM_DIR)
		print_counts(locus_tag, counts)
		print(f"Histogram saved to {histogram_path}")

		if index != len(locus_tags) - 1:
			print()

	if missing_locus_tags or data_error_locus_tags:
		raise SystemExit(1)


if __name__ == "__main__":
	main()
