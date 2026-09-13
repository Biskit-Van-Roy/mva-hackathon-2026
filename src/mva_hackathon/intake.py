"""Inspect the gated inputs without exporting patient-level variant rows or prose."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO


HPO_PATTERN = re.compile(r"HP:\d{7}", re.IGNORECASE)
GRCH38_CHR1_LENGTH = 248_956_422


def _open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_meta_identifier(line: str) -> str | None:
    match = re.search(r"ID=([^,>]+)", line)
    return match.group(1) if match else None


def _infer_build(reference: str | None, contig_lengths: dict[str, int]) -> str:
    reference_lower = (reference or "").lower()
    if any(token in reference_lower for token in ("grch38", "hg38")):
        return "GRCh38"
    chr1_length = contig_lengths.get("chr1") or contig_lengths.get("1")
    if chr1_length == GRCH38_CHR1_LENGTH:
        return "GRCh38"
    if chr1_length:
        return f"unknown (chr1 length={chr1_length})"
    return "unknown"


def _variant_type(ref: str, alts: list[str]) -> str:
    if len(alts) > 1:
        return "multiallelic"
    alt = alts[0]
    if alt.startswith("<") or "[" in alt or "]" in alt or alt == "*":
        return "symbolic_or_breakend"
    if len(ref) == len(alt) == 1:
        return "snv"
    return "indel_or_mnv"


def _classify_gt(sample_field: str, format_field: str) -> str:
    keys = format_field.split(":")
    try:
        gt_index = keys.index("GT")
    except ValueError:
        return "GT_missing_from_FORMAT"
    values = sample_field.split(":")
    if gt_index >= len(values):
        return "GT_missing_value"
    gt = values[gt_index].replace("|", "/")
    if gt in {".", "./."}:
        return "no_call"
    alleles = gt.split("/")
    if all(allele == "0" for allele in alleles):
        return "hom_ref"
    if all(allele == alleles[0] and allele not in {"0", "."} for allele in alleles):
        return "hom_alt"
    if any(allele not in {"0", "."} for allele in alleles):
        return "het_or_multiallelic_alt"
    return "other"


def inspect_vcf(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"VCF not found: {path}")

    info_ids: set[str] = set()
    format_ids: set[str] = set()
    filter_definitions: set[str] = set()
    contig_lengths: dict[str, int] = {}
    reference: str | None = None
    sample_count = 0
    sample_fingerprint: str | None = None
    variant_count = 0
    chrom_counts: Counter[str] = Counter()
    filter_counts: Counter[str] = Counter()
    type_counts: Counter[str] = Counter()
    genotype_counts: Counter[str] = Counter()
    malformed_rows = 0

    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("##reference="):
                reference = line.rstrip().split("=", 1)[1]
                continue
            if line.startswith("##contig=<"):
                contig_id = _parse_meta_identifier(line)
                length_match = re.search(r"length=(\d+)", line, re.IGNORECASE)
                if contig_id and length_match:
                    contig_lengths[contig_id] = int(length_match.group(1))
                continue
            if line.startswith("##INFO=<"):
                identifier = _parse_meta_identifier(line)
                if identifier:
                    info_ids.add(identifier)
                continue
            if line.startswith("##FORMAT=<"):
                identifier = _parse_meta_identifier(line)
                if identifier:
                    format_ids.add(identifier)
                continue
            if line.startswith("##FILTER=<"):
                identifier = _parse_meta_identifier(line)
                if identifier:
                    filter_definitions.add(identifier)
                continue
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                sample_names = columns[9:]
                sample_count = len(sample_names)
                if sample_names:
                    joined = "\0".join(sample_names).encode("utf-8")
                    sample_fingerprint = hashlib.sha256(joined).hexdigest()[:12]
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                malformed_rows += 1
                continue

            chrom, _pos, _id, ref, alt, _qual, filter_value, _info = fields[:8]
            alts = alt.split(",")
            variant_count += 1
            chrom_counts[chrom] += 1
            type_counts[_variant_type(ref, alts)] += 1
            for current_filter in filter_value.split(";"):
                filter_counts[current_filter or "EMPTY"] += 1
            if len(fields) >= 10:
                genotype_counts[_classify_gt(fields[9], fields[8])] += 1

    build = _infer_build(reference, contig_lengths)
    annotations = {
        "vep_csq": "CSQ" in info_ids,
        "snpeff_ann": "ANN" in info_ids,
        "clinvar_like": any(key in info_ids for key in ("CLNSIG", "CLNREVSTAT", "CLNDN")),
        "population_af_like": sorted(
            key for key in info_ids if key.upper() in {"AF", "GNOMAD_AF", "GNOMAD_EXOMES_AF", "MAX_AF"}
        ),
    }
    return {
        "path_name": path.name,
        "compressed_size_bytes": path.stat().st_size,
        "reference_header": reference,
        "inferred_build": build,
        "grch38_ready": build == "GRCh38",
        "sample_count": sample_count,
        "sample_name_fingerprint": sample_fingerprint,
        "contig_count": len(contig_lengths),
        "variant_count": variant_count,
        "malformed_rows": malformed_rows,
        "chromosome_counts": dict(chrom_counts),
        "variant_type_counts": dict(type_counts),
        "filter_counts": dict(filter_counts),
        "genotype_counts_first_sample": dict(genotype_counts),
        "info_fields": sorted(info_ids),
        "format_fields": sorted(format_ids),
        "filter_definitions": sorted(filter_definitions),
        "annotation_flags": annotations,
    }


def extract_hpo_ids_from_docx(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Phenotype document not found: {path}")
    with zipfile.ZipFile(path) as archive:
        document_xml = archive.read("word/document.xml")
    root = ET.fromstring(document_xml)
    text = " ".join(node.text or "" for node in root.iter() if node.tag.endswith("}t"))
    return sorted({match.upper() for match in HPO_PATTERN.findall(text)})


def build_report(vcf_path: Path, phenotype_path: Path | None) -> dict:
    hpo_ids = extract_hpo_ids_from_docx(phenotype_path) if phenotype_path else []
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "privacy_note": "Aggregate technical metadata only; no clinical prose or individual variant rows.",
        "vcf": inspect_vcf(vcf_path),
        "phenotype": {
            "document_present": phenotype_path is not None,
            "hpo_ids": hpo_ids,
            "hpo_count": len(hpo_ids),
        },
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True, help="Path to .vcf or .vcf.gz")
    parser.add_argument("--phenotype", type=Path, help="Optional clinical phenotype .docx")
    parser.add_argument("--output", type=Path, required=True, help="Output aggregate JSON")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        report = build_report(args.vcf, args.phenotype)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (FileNotFoundError, KeyError, OSError, ET.ParseError, zipfile.BadZipFile) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    vcf = report["vcf"]
    print(f"Intake complete: {vcf['variant_count']:,} VCF rows; build={vcf['inferred_build']}")
    print(f"HPO terms detected: {report['phenotype']['hpo_count']}")
    print(f"Report: {args.output}")
    if not vcf["grch38_ready"]:
        print("WARNING: GRCh38 could not be confirmed from the VCF header.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
