#!/usr/bin/env python3
"""Filter a VEP-annotated single-sample VCF without bcftools split-vep.

The script parses the INFO/CSQ schema from the VCF header, retains broad
rare-disease candidates, copies retained VCF records verbatim, and writes one
representative passing CSQ annotation per retained record to a private TSV.
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
from contextlib import ExitStack
from pathlib import Path
from typing import IO, Iterable


CSQ_HEADER_RE = re.compile(r"^##INFO=<ID=CSQ,.*?Format: ([^\"]+)")

PATHOGENIC_CLIN_SIG = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic_low_penetrance",
    "likely_pathogenic_low_penetrance",
}
CONFLICTING_CLIN_SIG = {"conflicting_interpretations_of_pathogenicity"}

ANNOTATION_FIELDS = [
    "SYMBOL",
    "Gene",
    "Feature",
    "BIOTYPE",
    "Consequence",
    "IMPACT",
    "HGVSc",
    "HGVSp",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "VARIANT_CLASS",
    "CANONICAL",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "ENSP",
    "SWISSPROT",
    "GENE_PHENO",
    "SIFT",
    "PolyPhen",
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
]

TSV_FIELDS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "GT",
    "AD",
    "DP",
    "GQ",
    *ANNOTATION_FIELDS,
    "FILTER_REASONS",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a broad rare-disease candidate set from VEP CSQ annotations."
    )
    parser.add_argument("--input", required=True, type=Path, help="VEP VCF or VCF.GZ")
    parser.add_argument(
        "--output-vcf", required=True, type=Path, help="Uncompressed candidate VCF"
    )
    parser.add_argument("--output-tsv", required=True, type=Path, help="Candidate TSV")
    parser.add_argument("--summary", required=True, type=Path, help="Summary JSON")
    parser.add_argument(
        "--max-af",
        type=float,
        default=0.01,
        help="Maximum population AF for functional candidates (default: 0.01)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Progress interval in input records (default: 1000000)",
    )
    args = parser.parse_args()
    if not 0 <= args.max_af <= 1:
        parser.error("--max-af must be between 0 and 1")
    return args


def open_text(path: Path) -> IO[str]:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
    return path.open("rt", encoding="utf-8", errors="strict", newline="")


def get_info_value(info: str, key: str) -> str | None:
    prefix = f"{key}="
    for item in info.split(";"):
        if item.startswith(prefix):
            return item[len(prefix) :]
    return None


def parse_float(value: str) -> float | None:
    if not value or value == ".":
        return None
    parsed: list[float] = []
    for token in value.split("&"):
        if token and token != ".":
            try:
                parsed.append(float(token))
            except ValueError:
                continue
    return max(parsed) if parsed else None


def clean(value: object) -> str:
    return str(value if value is not None else "").replace("\t", " ").replace("\n", " ")


def clin_sig_tokens(value: str) -> set[str]:
    """Return normalized, exact ClinVar significance labels from a CSQ field."""
    tokens: set[str] = set()
    for token in re.split(r"[&,|/]", value):
        normalized = re.sub(r"[\s-]+", "_", token.strip().lower())
        if normalized:
            tokens.add(normalized)
    return tokens


def has_pathogenic_clin_sig(value: str) -> bool:
    return bool(clin_sig_tokens(value) & PATHOGENIC_CLIN_SIG)


def has_conflicting_clin_sig(value: str) -> bool:
    return bool(clin_sig_tokens(value) & CONFLICTING_CLIN_SIG)


def annotation_score(annotation: dict[str, str]) -> tuple[int, int, int, float]:
    clinical = int(has_pathogenic_clin_sig(annotation.get("CLIN_SIG", "")))
    impact = annotation.get("IMPACT", "").upper()
    severity = 3 if impact == "HIGH" else 2 if impact == "MODERATE" else 0
    splice = int("splice_region_variant" in annotation.get("Consequence", ""))
    max_af = parse_float(annotation.get("MAX_AF", ""))
    rarity = 1.0 if max_af is None else 1.0 - max_af
    return clinical, severity, splice, rarity


def evaluate_annotation(
    annotation: dict[str, str], max_af_threshold: float
) -> tuple[bool, list[str]]:
    impact = annotation.get("IMPACT", "").upper()
    consequence = annotation.get("Consequence", "")
    clin_sig = annotation.get("CLIN_SIG", "")
    max_af = parse_float(annotation.get("MAX_AF", ""))

    high = impact == "HIGH"
    moderate = impact == "MODERATE"
    splice_region = "splice_region_variant" in consequence
    functional = high or moderate or splice_region
    rare = max_af is None or max_af <= max_af_threshold
    clinical = has_pathogenic_clin_sig(clin_sig)

    reasons: list[str] = []
    if high:
        reasons.append("IMPACT_HIGH")
    if moderate:
        reasons.append("IMPACT_MODERATE")
    if splice_region:
        reasons.append("SPLICE_REGION")
    if max_af is None:
        reasons.append("MAX_AF_MISSING")
    elif max_af <= max_af_threshold:
        reasons.append("MAX_AF_WITHIN_THRESHOLD")
    if clinical:
        reasons.append("CLIN_SIG_PATHOGENIC_OR_LIKELY_PATHOGENIC")

    return (functional and rare) or clinical, reasons


def parse_annotations(csq_value: str, csq_fields: list[str]) -> Iterable[dict[str, str]]:
    for raw_annotation in csq_value.split(","):
        values = raw_annotation.split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        yield dict(zip(csq_fields, values))


def format_sample(format_text: str, sample_text: str) -> dict[str, str]:
    if not format_text or format_text == "." or not sample_text:
        return {}
    keys = format_text.split(":")
    values = sample_text.split(":")
    return dict(zip(keys, values))


def partial_path(path: Path) -> Path:
    return path.with_name(path.name + ".partial")


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(f"Input does not exist: {args.input}")

    for path in (args.output_vcf, args.output_tsv, args.summary):
        path.parent.mkdir(parents=True, exist_ok=True)

    vcf_partial = partial_path(args.output_vcf)
    tsv_partial = partial_path(args.output_tsv)
    summary_partial = partial_path(args.summary)
    for path in (vcf_partial, tsv_partial, summary_partial):
        path.unlink(missing_ok=True)

    counters: dict[str, int | float | str] = {
        "filter_version": "BroadRareDiseaseV2",
        "input_records": 0,
        "records_with_csq": 0,
        "candidate_records": 0,
        "impact_high": 0,
        "impact_moderate": 0,
        "splice_region": 0,
        "clin_sig_pathogenic_or_likely_pathogenic": 0,
        "clin_sig_conflicting": 0,
        "max_af_missing": 0,
        "max_af_at_most_0_001": 0,
        "max_af_threshold": args.max_af,
    }

    csq_fields: list[str] | None = None
    sample_count = 0

    try:
        with ExitStack() as stack:
            source = stack.enter_context(open_text(args.input))
            vcf_out = stack.enter_context(
                vcf_partial.open("wt", encoding="utf-8", newline="")
            )
            tsv_handle = stack.enter_context(
                tsv_partial.open("wt", encoding="utf-8", newline="")
            )
            tsv_out = csv.DictWriter(
                tsv_handle,
                fieldnames=TSV_FIELDS,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            tsv_out.writeheader()

            for line_number, line in enumerate(source, start=1):
                if line.startswith("##"):
                    match = CSQ_HEADER_RE.match(line)
                    if match:
                        csq_fields = match.group(1).rstrip('\">').split("|")
                    vcf_out.write(line)
                    continue

                if line.startswith("#CHROM"):
                    columns = line.rstrip("\r\n").split("\t")
                    sample_count = max(0, len(columns) - 9)
                    vcf_out.write(
                        "##candidate_filter=<ID=BroadRareDiseaseV2,"
                        f"MaxPopulationAF={args.max_af},"
                        'Description="Retain HIGH/MODERATE/splice-region annotations '
                        "when MAX_AF is missing or within threshold, plus annotations "
                        "whose CLIN_SIG has an exact pathogenic or likely-pathogenic "
                        'classification; conflicting classifications are not a clinical override">\n'
                    )
                    vcf_out.write(line)
                    continue

                if line.startswith("#") or not line.strip():
                    vcf_out.write(line)
                    continue

                if csq_fields is None:
                    raise ValueError("The VCF header does not define INFO/CSQ Format fields")

                columns = line.rstrip("\r\n").split("\t")
                if len(columns) < 8:
                    raise ValueError(f"Malformed VCF record at line {line_number}")

                counters["input_records"] = int(counters["input_records"]) + 1
                input_records = int(counters["input_records"])
                if args.progress_every and input_records % args.progress_every == 0:
                    print(
                        f"Processed {input_records:,} input records; "
                        f"retained {int(counters['candidate_records']):,}",
                        file=sys.stderr,
                        flush=True,
                    )

                csq_value = get_info_value(columns[7], "CSQ")
                if not csq_value or csq_value == ".":
                    continue

                counters["records_with_csq"] = int(counters["records_with_csq"]) + 1
                passing: list[tuple[dict[str, str], list[str]]] = []
                for annotation in parse_annotations(csq_value, csq_fields):
                    keep, reasons = evaluate_annotation(annotation, args.max_af)
                    if keep:
                        passing.append((annotation, reasons))

                if not passing:
                    continue

                selected, reasons = max(passing, key=lambda item: annotation_score(item[0]))
                counters["candidate_records"] = int(counters["candidate_records"]) + 1

                impact = selected.get("IMPACT", "").upper()
                consequence = selected.get("Consequence", "")
                clin_sig = selected.get("CLIN_SIG", "")
                max_af = parse_float(selected.get("MAX_AF", ""))
                if impact == "HIGH":
                    counters["impact_high"] = int(counters["impact_high"]) + 1
                if impact == "MODERATE":
                    counters["impact_moderate"] = int(counters["impact_moderate"]) + 1
                if "splice_region_variant" in consequence:
                    counters["splice_region"] = int(counters["splice_region"]) + 1
                if has_pathogenic_clin_sig(clin_sig):
                    counters["clin_sig_pathogenic_or_likely_pathogenic"] = int(
                        counters["clin_sig_pathogenic_or_likely_pathogenic"]
                    ) + 1
                if has_conflicting_clin_sig(clin_sig):
                    counters["clin_sig_conflicting"] = int(
                        counters["clin_sig_conflicting"]
                    ) + 1
                if max_af is None:
                    counters["max_af_missing"] = int(counters["max_af_missing"]) + 1
                elif max_af <= 0.001:
                    counters["max_af_at_most_0_001"] = int(
                        counters["max_af_at_most_0_001"]
                    ) + 1

                sample = format_sample(
                    columns[8] if len(columns) > 8 else "",
                    columns[9] if len(columns) > 9 else "",
                )
                row = {
                    "CHROM": columns[0],
                    "POS": columns[1],
                    "ID": columns[2],
                    "REF": columns[3],
                    "ALT": columns[4],
                    "QUAL": columns[5],
                    "FILTER": columns[6],
                    "GT": sample.get("GT", ""),
                    "AD": sample.get("AD", ""),
                    "DP": sample.get("DP", ""),
                    "GQ": sample.get("GQ", ""),
                    **{field: selected.get(field, "") for field in ANNOTATION_FIELDS},
                    "FILTER_REASONS": ";".join(reasons),
                }
                tsv_out.writerow({key: clean(value) for key, value in row.items()})
                vcf_out.write(line)

        counters["sample_count"] = sample_count
        counters["csq_field_count"] = len(csq_fields or [])
        with summary_partial.open("wt", encoding="utf-8") as summary_handle:
            json.dump(counters, summary_handle, indent=2, sort_keys=True)
            summary_handle.write("\n")

        os.replace(vcf_partial, args.output_vcf)
        os.replace(tsv_partial, args.output_tsv)
        os.replace(summary_partial, args.summary)
    except Exception:
        for path in (vcf_partial, tsv_partial, summary_partial):
            path.unlink(missing_ok=True)
        raise

    print(json.dumps(counters, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
