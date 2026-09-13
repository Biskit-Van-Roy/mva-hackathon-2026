"""Build and validate Track 1 submission CSV files before using a submission."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path


COLUMNS = [
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
]


@dataclass
class ValidationResult:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def _is_positive_integer(value: str) -> bool:
    try:
        return int(value) > 0 and str(int(value)) == value.strip()
    except (TypeError, ValueError):
        return False


def _validate_variant(row: dict[str, str], suffix: str, row_number: int, result: ValidationResult) -> None:
    chrom = row[f"chrom_{suffix}"].strip()
    pos = row[f"pos_{suffix}"].strip()
    ref = row[f"ref_{suffix}"].strip()
    alt = row[f"alt_{suffix}"].strip()
    values = [chrom, pos, ref, alt]

    if suffix == "2" and not any(values):
        return
    if not all(values):
        result.errors.append(f"Row {row_number}: variant {suffix} is only partially populated.")
        return
    if not chrom.startswith("chr"):
        result.warnings.append(
            f"Row {row_number}: '{chrom}' does not use the official example convention 'chrN'; verify the VCF spelling."
        )
    if not _is_positive_integer(pos):
        result.errors.append(f"Row {row_number}: pos_{suffix} must be a positive integer.")
    if not ref or not alt:
        result.errors.append(f"Row {row_number}: ref_{suffix} and alt_{suffix} are required.")
    if ref != ref.upper() or alt != alt.upper():
        result.warnings.append(f"Row {row_number}: alleles should be uppercase.")


def validate_rows(fieldnames: list[str] | None, rows: list[dict[str, str]]) -> ValidationResult:
    result = ValidationResult()
    if fieldnames != COLUMNS:
        result.errors.append(
            "CSV columns or order do not match the official format. Expected: " + ",".join(COLUMNS)
        )
        return result
    if not 1 <= len(rows) <= 10:
        result.errors.append(f"Submission must contain 1 to 10 rows; found {len(rows)}.")

    previous_epcr = math.inf
    seen: set[tuple[str, ...]] = set()
    for index, row in enumerate(rows, start=2):
        if row["proband_id"].strip() != "PROBAND01":
            result.errors.append(f"Row {index}: only proband_id PROBAND01 is accepted.")

        _validate_variant(row, "1", index, result)
        _validate_variant(row, "2", index, result)

        try:
            epcr = float(row["epcr"])
            if not math.isfinite(epcr) or not 0 < epcr <= 1:
                raise ValueError
            if epcr > previous_epcr:
                result.warnings.append(
                    f"Row {index}: EPCR order is not descending; the server will re-sort it."
                )
            previous_epcr = epcr
        except ValueError:
            result.errors.append(f"Row {index}: epcr must be a finite number in (0, 1].")

        finding_type = row["finding_type"].strip().lower()
        if finding_type not in {"primary", "secondary"}:
            result.errors.append(f"Row {index}: finding_type must be primary or secondary.")

        key = tuple(row[column].strip().upper() for column in COLUMNS[1:9])
        if key in seen:
            result.errors.append(f"Row {index}: duplicate variant or pair.")
        seen.add(key)

    return result


def validate_csv(path: Path) -> ValidationResult:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return validate_rows(reader.fieldnames, rows)


def build_csv(input_path: Path, output_path: Path) -> ValidationResult:
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    candidates = payload.get("candidates") if isinstance(payload, dict) else payload
    if not isinstance(candidates, list):
        raise ValueError("JSON must be a list or an object with a 'candidates' list.")

    normalized: list[dict[str, str]] = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ValueError("Every candidate must be a JSON object.")
        row = {column: str(candidate.get(column, "")) for column in COLUMNS}
        row["proband_id"] = row["proband_id"] or "PROBAND01"
        row["finding_type"] = row["finding_type"] or "primary"
        for allele_column in ("ref_1", "alt_1", "ref_2", "alt_2"):
            row[allele_column] = row[allele_column].upper()
        normalized.append(row)

    normalized.sort(key=lambda row: float(row["epcr"]), reverse=True)
    result = validate_rows(COLUMNS, normalized)
    if not result.valid:
        return result

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)
    return result


def _print_result(result: ValidationResult) -> None:
    for warning in result.warnings:
        print(f"WARNING: {warning}")
    for error in result.errors:
        print(f"ERROR: {error}", file=sys.stderr)
    print("VALID" if result.valid else "INVALID")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser("validate", help="Validate an existing CSV")
    validate_parser.add_argument("csv_path", type=Path)

    build_parser = subparsers.add_parser("build", help="Build a CSV from candidate JSON")
    build_parser.add_argument("input_json", type=Path)
    build_parser.add_argument("output_csv", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "validate":
            result = validate_csv(args.csv_path)
        else:
            result = build_csv(args.input_json, args.output_csv)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    _print_result(result)
    return 0 if result.valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
