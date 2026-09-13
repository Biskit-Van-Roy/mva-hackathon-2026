#!/usr/bin/env python3
"""Add non-destructive genotype QC labels to a private candidate TSV.

This script never removes candidate records. It appends normalized genotype,
allele-balance, and heuristic review fields, then writes an aggregate JSON
summary that contains no loci, alleles, genes, or patient-identifying rows.
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable


PATHOGENIC_CLIN_SIG = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic_low_penetrance",
    "likely_pathogenic_low_penetrance",
}
CONFLICTING_CLIN_SIG = {"conflicting_interpretations_of_pathogenicity"}
APPENDED_FIELDS = ["GT_NORMALIZED", "GT_CLASS", "ALLELE_BALANCE", "QC_TIER", "QC_REASONS"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add heuristic, non-destructive genotype QC labels to candidate TSV rows."
    )
    parser.add_argument("--input-tsv", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--min-dp", type=float, default=10.0)
    parser.add_argument("--min-gq", type=float, default=20.0)
    parser.add_argument("--het-ab-min", type=float, default=0.20)
    parser.add_argument("--het-ab-max", type=float, default=0.80)
    parser.add_argument("--hom-alt-ab-min", type=float, default=0.80)
    args = parser.parse_args()
    if args.min_dp < 0 or args.min_gq < 0:
        parser.error("DP and GQ thresholds must be non-negative")
    if not 0 <= args.het_ab_min <= args.het_ab_max <= 1:
        parser.error("heterozygous allele-balance bounds must satisfy 0 <= min <= max <= 1")
    if not 0 <= args.hom_alt_ab_min <= 1:
        parser.error("--hom-alt-ab-min must be between 0 and 1")
    return args


def partial_path(path: Path) -> Path:
    return path.with_name(path.name + ".partial")


def parse_number(value: str) -> float | None:
    value = (value or "").strip()
    if not value or value == ".":
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def parse_max_af(value: str) -> float | None:
    values: list[float] = []
    for token in (value or "").split("&"):
        number = parse_number(token)
        if number is not None:
            values.append(number)
    return max(values) if values else None


def clin_sig_tokens(value: str) -> set[str]:
    tokens: set[str] = set()
    for token in re.split(r"[&,|/]", value or ""):
        normalized = re.sub(r"[\s-]+", "_", token.strip().lower())
        if normalized:
            tokens.add(normalized)
    return tokens


def normalize_gt(value: str) -> str:
    value = (value or "").strip().replace("|", "/")
    if not value or value in {".", "./."}:
        return "./."
    return value


def classify_gt(value: str) -> str:
    gt = normalize_gt(value)
    if gt == "./.":
        return "NO_CALL"
    alleles = gt.split("/")
    if any(not allele.isdigit() for allele in alleles):
        return "OTHER"
    numeric = [int(allele) for allele in alleles]
    if len(numeric) == 1:
        return "HEMI_ALT" if numeric[0] > 0 else "HEMI_REF"
    if all(allele == 0 for allele in numeric):
        return "HOM_REF"
    if len(set(numeric)) == 1 and numeric[0] > 0:
        return "HOM_ALT"
    if any(allele > 0 for allele in numeric):
        return "HET"
    return "OTHER"


def allele_balance(value: str) -> float | None:
    tokens = (value or "").split(",")
    if len(tokens) < 2:
        return None
    depths: list[float] = []
    for token in tokens:
        number = parse_number(token)
        if number is None or number < 0:
            return None
        depths.append(number)
    total = sum(depths)
    if total <= 0:
        return None
    return sum(depths[1:]) / total


def quantiles(values: Iterable[float]) -> dict[str, float | int | None]:
    ordered = sorted(values)
    if not ordered:
        return {"n": 0, "min": None, "p25": None, "median": None, "p75": None, "max": None}

    def q(fraction: float) -> float:
        position = (len(ordered) - 1) * fraction
        lower = math.floor(position)
        upper = math.ceil(position)
        if lower == upper:
            return ordered[lower]
        weight = position - lower
        return ordered[lower] * (1 - weight) + ordered[upper] * weight

    return {
        "n": len(ordered),
        "min": round(ordered[0], 3),
        "p25": round(q(0.25), 3),
        "median": round(q(0.50), 3),
        "p75": round(q(0.75), 3),
        "max": round(ordered[-1], 3),
    }


def frequency_bin(value: float | None) -> str:
    if value is None:
        return "missing"
    if value <= 0.00001:
        return "at_most_0.00001"
    if value <= 0.0001:
        return "0.00001_to_0.0001"
    if value <= 0.001:
        return "0.0001_to_0.001"
    if value <= 0.01:
        return "0.001_to_0.01"
    return "above_0.01_clinical_override"


def qc_row(row: dict[str, str], args: argparse.Namespace) -> tuple[dict[str, str], str, list[str], float | None]:
    gt_normalized = normalize_gt(row.get("GT", ""))
    gt_class = classify_gt(gt_normalized)
    dp = parse_number(row.get("DP", ""))
    gq = parse_number(row.get("GQ", ""))
    ab = allele_balance(row.get("AD", ""))

    review_reasons: list[str] = []
    missing_reasons: list[str] = []

    if row.get("FILTER", "") not in {"PASS", "."}:
        review_reasons.append("VARIANT_FILTER_NOT_PASS")
    if gt_class in {"NO_CALL", "HOM_REF", "HEMI_REF", "OTHER"}:
        review_reasons.append(f"GT_{gt_class}")
    if dp is None:
        missing_reasons.append("DP_MISSING")
    elif dp < args.min_dp:
        review_reasons.append("DP_BELOW_THRESHOLD")
    if gq is None:
        missing_reasons.append("GQ_MISSING")
    elif gq < args.min_gq:
        review_reasons.append("GQ_BELOW_THRESHOLD")
    if ab is None:
        missing_reasons.append("AD_OR_ALLELE_BALANCE_MISSING")
    elif gt_class == "HET" and not args.het_ab_min <= ab <= args.het_ab_max:
        review_reasons.append("HET_ALLELE_BALANCE_OUTSIDE_RANGE")
    elif gt_class in {"HOM_ALT", "HEMI_ALT"} and ab < args.hom_alt_ab_min:
        review_reasons.append("ALT_ALLELE_BALANCE_BELOW_THRESHOLD")

    if review_reasons:
        tier = "REVIEW_LOW_CONFIDENCE"
    elif missing_reasons:
        tier = "REVIEW_MISSING_METRICS"
    else:
        tier = "PASS_HEURISTIC_QC"

    reasons = review_reasons + missing_reasons
    enriched = dict(row)
    enriched.update(
        {
            "GT_NORMALIZED": gt_normalized,
            "GT_CLASS": gt_class,
            "ALLELE_BALANCE": "" if ab is None else f"{ab:.4f}",
            "QC_TIER": tier,
            "QC_REASONS": ";".join(reasons),
        }
    )
    return enriched, tier, reasons, ab


def main() -> int:
    args = parse_args()
    if not args.input_tsv.is_file():
        raise FileNotFoundError(f"Input TSV does not exist: {args.input_tsv}")
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    output_partial = partial_path(args.output_tsv)
    summary_partial = partial_path(args.summary)
    output_partial.unlink(missing_ok=True)
    summary_partial.unlink(missing_ok=True)

    counters: dict[str, Counter[str]] = {
        "filter_counts": Counter(),
        "gt_class_counts": Counter(),
        "qc_tier_counts": Counter(),
        "qc_reason_counts": Counter(),
        "impact_counts": Counter(),
        "max_af_bins": Counter(),
    }
    dp_values: list[float] = []
    gq_values: list[float] = []
    het_ab_values: list[float] = []
    gene_counts: Counter[str] = Counter()
    total_records = 0
    pathogenic_count = 0
    conflicting_count = 0

    csv.field_size_limit(sys.maxsize)
    try:
        with args.input_tsv.open("rt", encoding="utf-8", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError("Input TSV has no header")
            required = {"FILTER", "GT", "AD", "DP", "GQ", "IMPACT", "MAX_AF", "CLIN_SIG", "Gene"}
            missing = sorted(required - set(reader.fieldnames))
            if missing:
                raise ValueError(f"Input TSV is missing required fields: {', '.join(missing)}")
            output_fields = list(reader.fieldnames) + [
                field for field in APPENDED_FIELDS if field not in reader.fieldnames
            ]
            with output_partial.open("wt", encoding="utf-8", newline="") as output:
                writer = csv.DictWriter(
                    output,
                    fieldnames=output_fields,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="ignore",
                )
                writer.writeheader()
                for row in reader:
                    total_records += 1
                    enriched, tier, reasons, ab = qc_row(row, args)
                    writer.writerow(enriched)

                    gt_class = enriched["GT_CLASS"]
                    counters["filter_counts"][row.get("FILTER", "") or "missing"] += 1
                    counters["gt_class_counts"][gt_class] += 1
                    counters["qc_tier_counts"][tier] += 1
                    counters["qc_reason_counts"].update(reasons)
                    counters["impact_counts"][row.get("IMPACT", "") or "missing"] += 1
                    counters["max_af_bins"][frequency_bin(parse_max_af(row.get("MAX_AF", "")))] += 1

                    dp = parse_number(row.get("DP", ""))
                    gq = parse_number(row.get("GQ", ""))
                    if dp is not None:
                        dp_values.append(dp)
                    if gq is not None:
                        gq_values.append(gq)
                    if gt_class == "HET" and ab is not None:
                        het_ab_values.append(ab)

                    gene = (row.get("Gene", "") or "").strip()
                    if gene:
                        gene_counts[gene] += 1
                    clin_tokens = clin_sig_tokens(row.get("CLIN_SIG", ""))
                    if clin_tokens & PATHOGENIC_CLIN_SIG:
                        pathogenic_count += 1
                    if clin_tokens & CONFLICTING_CLIN_SIG:
                        conflicting_count += 1

        multi_gene_counts = [count for count in gene_counts.values() if count >= 2]
        summary = {
            "summary_version": "CandidateGenotypeQCV1",
            "total_records": total_records,
            "thresholds_are_heuristic_not_diagnostic": True,
            "thresholds": {
                "min_dp": args.min_dp,
                "min_gq": args.min_gq,
                "het_ab_min": args.het_ab_min,
                "het_ab_max": args.het_ab_max,
                "hom_alt_ab_min": args.hom_alt_ab_min,
            },
            "filter_counts": dict(sorted(counters["filter_counts"].items())),
            "gt_class_counts": dict(sorted(counters["gt_class_counts"].items())),
            "qc_tier_counts": dict(sorted(counters["qc_tier_counts"].items())),
            "qc_reason_counts": dict(sorted(counters["qc_reason_counts"].items())),
            "impact_counts": dict(sorted(counters["impact_counts"].items())),
            "max_af_bins": dict(sorted(counters["max_af_bins"].items())),
            "clin_sig_pathogenic_or_likely_pathogenic": pathogenic_count,
            "clin_sig_conflicting": conflicting_count,
            "distinct_ensembl_gene_ids": len(gene_counts),
            "genes_with_at_least_two_candidates": len(multi_gene_counts),
            "records_in_genes_with_at_least_two_candidates": sum(multi_gene_counts),
            "dp_distribution": quantiles(dp_values),
            "gq_distribution": quantiles(gq_values),
            "het_allele_balance_distribution": quantiles(het_ab_values),
        }
        with summary_partial.open("wt", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        os.replace(output_partial, args.output_tsv)
        os.replace(summary_partial, args.summary)
    except Exception:
        output_partial.unlink(missing_ok=True)
        summary_partial.unlink(missing_ok=True)
        raise

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
