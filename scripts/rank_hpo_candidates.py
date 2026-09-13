#!/usr/bin/env python3
"""Rank private rare-disease candidates using transparent HPO and variant evidence.

The script produces two rankings:
1. an unbiased ranking based on phenotype, variant evidence, QC, and a
   non-confirmatory recessive-model signal;
2. an MVA-aware ranking that adds an explicit, auditable mechanism-gene bonus.

Scores are prioritization heuristics, not probabilities or clinical
classifications. No candidate is removed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path


RANKER_VERSION = "HPOVariantRankerV2"
PATHOGENIC_CLIN_SIG = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic_low_penetrance",
    "likely_pathogenic_low_penetrance",
}
CONFLICTING_CLIN_SIG = {"conflicting_interpretations_of_pathogenicity"}

VARIANT_OUTPUT_FIELDS = [
    "HPO_SCORE",
    "FAMILY_HISTORY_HPO_SCORE",
    "VARIANT_EVIDENCE_SCORE",
    "INHERITANCE_SCORE",
    "INHERITANCE_HYPOTHESIS",
    "MVA_MECHANISM_GENE",
    "SCORE_UNBIASED",
    "RANK_UNBIASED",
    "SCORE_MVA_AWARE",
    "RANK_MVA_AWARE",
]

GENE_OUTPUT_FIELDS = [
    "GENE_GROUP_KEY",
    "GENE_SYMBOL",
    "CANDIDATE_COUNT",
    "PASS_QC_COUNT",
    "HET_PASS_QC_COUNT",
    "HOM_ALT_PASS_QC_COUNT",
    "HPO_ANNOTATION_COUNT",
    "HPO_SCORE",
    "FAMILY_HISTORY_HPO_SCORE",
    "TOP_VARIANT_EVIDENCE_SCORE",
    "SECOND_VARIANT_EVIDENCE_SCORE",
    "INHERITANCE_SCORE",
    "INHERITANCE_HYPOTHESIS",
    "MVA_MECHANISM_GENE",
    "GENE_SCORE_UNBIASED",
    "GENE_RANK_UNBIASED",
    "GENE_SCORE_MVA_AWARE",
    "GENE_RANK_MVA_AWARE",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank VEP candidates with HPO semantic similarity and transparent evidence components."
    )
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--hp-obo", required=True, type=Path)
    parser.add_argument("--genes-to-phenotype", required=True, type=Path)
    parser.add_argument("--proband-hpo", required=True, type=Path)
    parser.add_argument("--family-hpo", type=Path)
    parser.add_argument("--output-variants", required=True, type=Path)
    parser.add_argument("--output-genes", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument(
        "--mechanism-gene",
        action="append",
        default=[],
        help="Gene symbol receiving an explicit bonus in the MVA-aware ranking; repeat as needed.",
    )
    parser.add_argument("--phenotype-weight", type=float, default=10.0)
    parser.add_argument("--family-weight", type=float, default=1.0)
    parser.add_argument("--mechanism-bonus", type=float, default=3.0)
    args = parser.parse_args()
    for path in (
        args.candidates,
        args.hp_obo,
        args.genes_to_phenotype,
        args.proband_hpo,
    ):
        if not path.is_file():
            parser.error(f"required input does not exist: {path}")
    if args.family_hpo is not None and not args.family_hpo.is_file():
        parser.error(f"family HPO input does not exist: {args.family_hpo}")
    if min(args.phenotype_weight, args.family_weight, args.mechanism_bonus) < 0:
        parser.error("score weights must be non-negative")
    return args


def partial_path(path: Path) -> Path:
    return path.with_name(path.name + ".partial")


def clean(value: object) -> str:
    return str(value if value is not None else "").replace("\t", " ").replace("\n", " ")


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
    values = [parse_number(token) for token in (value or "").split("&")]
    present = [value for value in values if value is not None]
    return max(present) if present else None


def clin_sig_tokens(value: str) -> set[str]:
    result: set[str] = set()
    for token in re.split(r"[&,|/]", value or ""):
        normalized = re.sub(r"[\s-]+", "_", token.strip().lower())
        if normalized:
            result.add(normalized)
    return result


def read_hpo_ids(path: Path | None) -> list[str]:
    if path is None:
        return []
    result: list[str] = []
    seen: set[str] = set()
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            value = raw.split("#", 1)[0].strip().split("\t", 1)[0]
            if not value:
                continue
            if not re.fullmatch(r"HP:\d{7}", value):
                raise ValueError(f"invalid HPO identifier in {path}: {value}")
            if value not in seen:
                result.append(value)
                seen.add(value)
    if path is not None and not result:
        raise ValueError(f"HPO input is empty: {path}")
    return result


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str], str | None]:
    parents: dict[str, set[str]] = defaultdict(set)
    alt_to_primary: dict[str, str] = {}
    data_version: str | None = None
    current_id: str | None = None
    current_parents: set[str] = set()
    current_alt_ids: list[str] = []
    obsolete = False

    def commit() -> None:
        if current_id and not obsolete:
            parents[current_id].update(current_parents)
            for alt_id in current_alt_ids:
                alt_to_primary[alt_id] = current_id

    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if data_version is None and line.startswith("data-version:"):
                data_version = line.split(":", 1)[1].strip()
            if line == "[Term]":
                commit()
                current_id = None
                current_parents = set()
                current_alt_ids = []
                obsolete = False
            elif line.startswith("id: HP:"):
                current_id = line.split("id:", 1)[1].strip()
            elif line.startswith("alt_id: HP:"):
                current_alt_ids.append(line.split("alt_id:", 1)[1].strip())
            elif line.startswith("is_a: HP:"):
                current_parents.add(line.split()[1])
            elif line == "is_obsolete: true":
                obsolete = True
    commit()
    return dict(parents), alt_to_primary, data_version


def read_gene_annotations(path: Path, canonicalize) -> dict[str, set[str]]:
    gene_terms: dict[str, set[str]] = defaultdict(set)
    csv.field_size_limit(sys.maxsize)
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"gene_symbol", "hpo_id"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError(
                "genes_to_phenotype header must contain gene_symbol and hpo_id"
            )
        for row in reader:
            symbol = (row.get("gene_symbol") or "").strip().upper()
            hpo_id = canonicalize((row.get("hpo_id") or "").strip())
            if symbol and hpo_id:
                gene_terms[symbol].add(hpo_id)
    return dict(gene_terms)


def read_candidates(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    csv.field_size_limit(sys.maxsize)
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("candidate TSV has no header")
        required = {
            "SYMBOL",
            "Consequence",
            "IMPACT",
            "MAX_AF",
            "CLIN_SIG",
            "SIFT",
            "PolyPhen",
            "QC_TIER",
            "GT_CLASS",
        }
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"candidate TSV is missing: {', '.join(missing)}")
        return [dict(row) for row in reader], list(reader.fieldnames)


def gene_identity(row: dict[str, str], row_index: int) -> tuple[str, str]:
    """Return a grouping key and display symbol without merging unrelated blanks."""
    symbol = (row.get("SYMBOL") or "").strip().upper()
    if symbol:
        return f"SYMBOL:{symbol}", symbol
    ensembl_gene = (row.get("Gene") or "").strip()
    if ensembl_gene:
        return f"ENSEMBL:{ensembl_gene}", "UNMAPPED_SYMBOL"
    chrom = (row.get("CHROM") or "NA").strip()
    pos = (row.get("POS") or str(row_index)).strip()
    return f"INTERGENIC:{chrom}:{pos}:{row_index}", "UNMAPPED"


def consequence_score(value: str) -> float:
    weights = {
        "transcript_ablation": 5.0,
        "splice_acceptor_variant": 5.0,
        "splice_donor_variant": 5.0,
        "stop_gained": 5.0,
        "frameshift_variant": 5.0,
        "start_lost": 4.0,
        "stop_lost": 4.0,
        "transcript_amplification": 3.5,
        "inframe_insertion": 3.0,
        "inframe_deletion": 3.0,
        "missense_variant": 2.5,
        "protein_altering_variant": 2.0,
        "splice_region_variant": 2.0,
    }
    return max((weights.get(term, 0.0) for term in (value or "").split("&")), default=0.0)


def rarity_score(value: str) -> float:
    af = parse_max_af(value)
    if af is None:
        return 1.0
    if af <= 0.00001:
        return 4.0
    if af <= 0.0001:
        return 3.5
    if af <= 0.001:
        return 2.5
    if af <= 0.01:
        return 1.0
    return -2.0


def variant_evidence_score(row: dict[str, str]) -> float:
    impact_scores = {"HIGH": 4.0, "MODERATE": 2.5, "LOW": 1.0, "MODIFIER": 0.0}
    score = impact_scores.get((row.get("IMPACT") or "").upper(), 0.0)
    score += consequence_score(row.get("Consequence", ""))
    score += rarity_score(row.get("MAX_AF", ""))

    clin_tokens = clin_sig_tokens(row.get("CLIN_SIG", ""))
    if clin_tokens & PATHOGENIC_CLIN_SIG:
        score += 3.0
    if clin_tokens & CONFLICTING_CLIN_SIG:
        score -= 1.0
    if "deleterious" in (row.get("SIFT") or "").lower():
        score += 0.5
    if "damaging" in (row.get("PolyPhen") or "").lower():
        score += 0.5

    qc_tier = row.get("QC_TIER", "")
    if qc_tier == "PASS_HEURISTIC_QC":
        score += 1.0
    elif qc_tier == "REVIEW_LOW_CONFIDENCE":
        score -= 1.5
    return round(score, 4)


def inheritance_for_rows(rows: list[dict[str, str]]) -> tuple[float, str]:
    pass_rows = [row for row in rows if row.get("QC_TIER") == "PASS_HEURISTIC_QC"]
    hom_alt = sum(row.get("GT_CLASS") == "HOM_ALT" for row in pass_rows)
    het = sum(row.get("GT_CLASS") == "HET" for row in pass_rows)
    if hom_alt:
        return 3.0, "HOMOZYGOUS_RECESSIVE_CANDIDATE"
    if het >= 2:
        return 2.5, "UNPHASED_COMPOUND_HET_HYPOTHESIS"
    if het == 1:
        return 0.5, "SINGLE_HETEROZYGOUS_CANDIDATE"
    return 0.0, "NO_STRONG_RECESSIVE_PATTERN"


def assign_ranks(items: list[dict[str, object]], score_field: str, rank_field: str, tie_field: str) -> None:
    ordered = sorted(
        range(len(items)),
        key=lambda index: (-float(items[index][score_field]), str(items[index].get(tie_field, "")), index),
    )
    for rank, index in enumerate(ordered, start=1):
        items[index][rank_field] = rank


def main() -> int:
    args = parse_args()
    rows, input_fields = read_candidates(args.candidates)
    if not rows:
        raise ValueError("candidate TSV contains no rows")

    parents, alt_to_primary, ontology_version = parse_obo(args.hp_obo)

    def canonicalize(term: str) -> str:
        if not term:
            return ""
        return alt_to_primary.get(term, term) if term in parents or term in alt_to_primary else ""

    @lru_cache(maxsize=None)
    def ancestors(term: str) -> frozenset[str]:
        canonical = canonicalize(term)
        if not canonical:
            return frozenset()
        result = {canonical}
        for parent in parents.get(canonical, set()):
            result.update(ancestors(parent))
        return frozenset(result)

    proband_terms = [canonicalize(term) for term in read_hpo_ids(args.proband_hpo)]
    family_terms = [canonicalize(term) for term in read_hpo_ids(args.family_hpo)]
    if any(not term for term in proband_terms + family_terms):
        raise ValueError("one or more query HPO terms are absent from the ontology")

    gene_terms = read_gene_annotations(args.genes_to_phenotype, canonicalize)
    total_reference_genes = len(gene_terms)
    term_gene_counts: Counter[str] = Counter()
    candidate_symbols = {(row.get("SYMBOL") or "").strip().upper() for row in rows}
    candidate_symbols.discard("")
    candidate_propagated: dict[str, set[str]] = {}

    for symbol, direct_terms in gene_terms.items():
        propagated: set[str] = set()
        for term in direct_terms:
            propagated.update(ancestors(term))
        term_gene_counts.update(propagated)
        if symbol in candidate_symbols:
            candidate_propagated[symbol] = propagated

    def information_content(term: str) -> float:
        return -math.log((term_gene_counts.get(term, 0) + 1) / (total_reference_genes + 1))

    def similarity(symbol: str, query_terms: list[str]) -> float:
        if not query_terms:
            return 0.0
        gene_ancestors = candidate_propagated.get(symbol, set())
        if not gene_ancestors:
            return 0.0
        term_scores: list[float] = []
        for query in query_terms:
            query_ancestors = ancestors(query)
            denominator = information_content(query)
            common = query_ancestors & gene_ancestors
            best = max((information_content(term) for term in common), default=0.0)
            term_scores.append(0.0 if denominator <= 0 else min(1.0, best / denominator))
        return sum(term_scores) / len(term_scores)

    rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    group_symbols: dict[str, str] = {}
    for row_index, row in enumerate(rows, start=1):
        group_key, display_symbol = gene_identity(row, row_index)
        row["__GENE_GROUP_KEY"] = group_key
        rows_by_gene[group_key].append(row)
        group_symbols[group_key] = display_symbol

    mechanism_genes = {gene.strip().upper() for gene in args.mechanism_gene if gene.strip()}
    gene_records: list[dict[str, object]] = []
    gene_components: dict[str, dict[str, object]] = {}

    for group_key, gene_rows in rows_by_gene.items():
        symbol = group_symbols[group_key]
        hpo_symbol = symbol if symbol not in {"UNMAPPED", "UNMAPPED_SYMBOL"} else ""
        evidence = sorted((variant_evidence_score(row) for row in gene_rows), reverse=True)
        top = evidence[0]
        second = evidence[1] if len(evidence) > 1 else 0.0
        hpo_score = similarity(hpo_symbol, proband_terms)
        family_score = similarity(hpo_symbol, family_terms)
        inheritance_score, inheritance_hypothesis = inheritance_for_rows(gene_rows)
        mechanism = hpo_symbol in mechanism_genes
        unbiased = (
            args.phenotype_weight * hpo_score
            + args.family_weight * family_score
            + top
            + 0.5 * second
            + inheritance_score
        )
        aware = unbiased + (args.mechanism_bonus if mechanism else 0.0)
        pass_rows = [row for row in gene_rows if row.get("QC_TIER") == "PASS_HEURISTIC_QC"]
        record: dict[str, object] = {
            "GENE_GROUP_KEY": group_key,
            "GENE_SYMBOL": symbol,
            "CANDIDATE_COUNT": len(gene_rows),
            "PASS_QC_COUNT": len(pass_rows),
            "HET_PASS_QC_COUNT": sum(row.get("GT_CLASS") == "HET" for row in pass_rows),
            "HOM_ALT_PASS_QC_COUNT": sum(row.get("GT_CLASS") == "HOM_ALT" for row in pass_rows),
            "HPO_ANNOTATION_COUNT": len(gene_terms.get(hpo_symbol, set())),
            "HPO_SCORE": round(hpo_score, 6),
            "FAMILY_HISTORY_HPO_SCORE": round(family_score, 6),
            "TOP_VARIANT_EVIDENCE_SCORE": round(top, 4),
            "SECOND_VARIANT_EVIDENCE_SCORE": round(second, 4),
            "INHERITANCE_SCORE": inheritance_score,
            "INHERITANCE_HYPOTHESIS": inheritance_hypothesis,
            "MVA_MECHANISM_GENE": "YES" if mechanism else "NO",
            "GENE_SCORE_UNBIASED": round(unbiased, 6),
            "GENE_SCORE_MVA_AWARE": round(aware, 6),
        }
        gene_records.append(record)
        gene_components[group_key] = record

    assign_ranks(gene_records, "GENE_SCORE_UNBIASED", "GENE_RANK_UNBIASED", "GENE_GROUP_KEY")
    assign_ranks(gene_records, "GENE_SCORE_MVA_AWARE", "GENE_RANK_MVA_AWARE", "GENE_GROUP_KEY")

    ranked_rows: list[dict[str, object]] = []
    for row in rows:
        group_key = row["__GENE_GROUP_KEY"]
        gene = gene_components[group_key]
        evidence = variant_evidence_score(row)
        inheritance_score = float(gene["INHERITANCE_SCORE"])
        unbiased = (
            args.phenotype_weight * float(gene["HPO_SCORE"])
            + args.family_weight * float(gene["FAMILY_HISTORY_HPO_SCORE"])
            + evidence
            + inheritance_score
        )
        mechanism = gene["MVA_MECHANISM_GENE"] == "YES"
        aware = unbiased + (args.mechanism_bonus if mechanism else 0.0)
        enriched: dict[str, object] = dict(row)
        enriched.update(
            {
                "HPO_SCORE": gene["HPO_SCORE"],
                "FAMILY_HISTORY_HPO_SCORE": gene["FAMILY_HISTORY_HPO_SCORE"],
                "VARIANT_EVIDENCE_SCORE": evidence,
                "INHERITANCE_SCORE": inheritance_score,
                "INHERITANCE_HYPOTHESIS": gene["INHERITANCE_HYPOTHESIS"],
                "MVA_MECHANISM_GENE": gene["MVA_MECHANISM_GENE"],
                "SCORE_UNBIASED": round(unbiased, 6),
                "SCORE_MVA_AWARE": round(aware, 6),
            }
        )
        ranked_rows.append(enriched)

    assign_ranks(ranked_rows, "SCORE_UNBIASED", "RANK_UNBIASED", "SYMBOL")
    assign_ranks(ranked_rows, "SCORE_MVA_AWARE", "RANK_MVA_AWARE", "SYMBOL")

    output_paths = (args.output_variants, args.output_genes, args.summary)
    for path in output_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    partials = [partial_path(path) for path in output_paths]
    for path in partials:
        path.unlink(missing_ok=True)

    try:
        variant_fields = input_fields + [field for field in VARIANT_OUTPUT_FIELDS if field not in input_fields]
        with partials[0].open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=variant_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in sorted(ranked_rows, key=lambda item: int(item["RANK_MVA_AWARE"])):
                writer.writerow({key: clean(row.get(key, "")) for key in variant_fields})

        with partials[1].open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=GENE_OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for record in sorted(gene_records, key=lambda item: int(item["GENE_RANK_MVA_AWARE"])):
                writer.writerow({key: clean(record.get(key, "")) for key in GENE_OUTPUT_FIELDS})

        summary = {
            "ranker_version": RANKER_VERSION,
            "scores_are_heuristics_not_probabilities_or_clinical_classifications": True,
            "candidate_records": len(rows),
            "candidate_gene_symbols": len(candidate_symbols),
            "candidate_gene_groups": len(rows_by_gene),
            "unmapped_symbol_records": sum(not (row.get("SYMBOL") or "").strip() for row in rows),
            "candidate_genes_with_hpo_annotations": sum(symbol in gene_terms for symbol in candidate_symbols),
            "reference_gene_count": total_reference_genes,
            "ontology_term_count": len(parents),
            "ontology_data_version": ontology_version,
            "proband_hpo_count": len(proband_terms),
            "family_history_hpo_count": len(family_terms),
            "mechanism_gene_count_configured": len(mechanism_genes),
            "candidate_records_in_mechanism_genes": sum(
                1 for row in rows if (row.get("SYMBOL") or "").strip().upper() in mechanism_genes
            ),
            "candidate_mechanism_genes_present": sum(gene in candidate_symbols for gene in mechanism_genes),
            "inheritance_hypothesis_counts": dict(
                sorted(Counter(str(record["INHERITANCE_HYPOTHESIS"]) for record in gene_records).items())
            ),
            "weights": {
                "phenotype": args.phenotype_weight,
                "family_history": args.family_weight,
                "mechanism_bonus": args.mechanism_bonus,
                "second_variant_gene_weight": 0.5,
            },
            "ranking_lanes": {
                "unbiased": "HPO + family-history HPO + variant evidence + unconfirmed recessive-pattern signal",
                "mva_aware": "unbiased score plus explicit mechanism-gene bonus",
            },
        }
        with partials[2].open("wt", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        for partial, final in zip(partials, output_paths):
            os.replace(partial, final)
    except Exception:
        for path in partials:
            path.unlink(missing_ok=True)
        raise

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
