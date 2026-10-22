"""Transparent, research-only ranking for Track 2 drug-repurposing evidence.

The score is an auditable prioritisation heuristic. It is not a probability,
clinical classification, treatment recommendation, or prescribing tool.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "track2-drug-evidence-v1"
REQUIRED_EVIDENCE_TYPES = {"primary_mechanistic", "regulatory", "safety"}
POSITIVE_COMPONENTS = {
    "mechanistic_fit": 5,
    "evidence_strength": 5,
    "regulatory_fit": 3,
    "experimental_testability": 4,
    "scalability": 4,
}
PENALTY_COMPONENTS = {
    "safety_risk": 5,
    "chromosomal_instability_risk": 5,
    "patient_specific_conflict": 5,
}
CSV_COLUMNS = [
    "RANK",
    "CANDIDATE_ID",
    "GENERIC_NAME",
    "ELIGIBLE_FOR_SHORTLIST",
    "BENEFIT_SCORE",
    "RISK_PENALTY",
    "FINAL_SCORE",
    *[name.upper() for name in POSITIVE_COMPONENTS],
    *[name.upper() for name in PENALTY_COMPONENTS],
    "TARGET_INTERVENTION",
    "EVIDENCE_COUNT",
    "EXCLUSION_REASONS",
    "RESEARCH_ONLY",
]


@dataclass
class ValidationResult:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


@dataclass(frozen=True)
class RankedCandidate:
    source: dict[str, Any]
    eligible: bool
    benefit_score: float
    risk_penalty: float
    final_score: float
    exclusion_reasons: tuple[str, ...]


def load_weights(path: Path) -> dict[str, dict[str, float]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"weights schema_version must be {SCHEMA_VERSION}")

    weights: dict[str, dict[str, float]] = {}
    for group_name, expected_components in (
        ("positive_weights", POSITIVE_COMPONENTS),
        ("penalty_weights", PENALTY_COMPONENTS),
    ):
        group = payload.get(group_name)
        if not isinstance(group, dict) or set(group) != set(expected_components):
            raise ValueError(
                f"{group_name} must contain exactly: {', '.join(expected_components)}"
            )
        normalized: dict[str, float] = {}
        for name, value in group.items():
            try:
                numeric = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"weight {group_name}.{name} must be numeric") from exc
            if not math.isfinite(numeric) or numeric < 0:
                raise ValueError(f"weight {group_name}.{name} must be finite and non-negative")
            normalized[name] = numeric
        weights[group_name] = normalized
    return weights


def _bounded_score(
    candidate: dict[str, Any],
    group_name: str,
    component_name: str,
    maximum: int,
    result: ValidationResult,
) -> float:
    scores = candidate.get(group_name)
    candidate_id = candidate.get("candidate_id", "<missing>")
    if not isinstance(scores, dict) or component_name not in scores:
        result.errors.append(
            f"{candidate_id}: missing {group_name}.{component_name}"
        )
        return 0.0
    try:
        value = float(scores[component_name])
    except (TypeError, ValueError):
        result.errors.append(
            f"{candidate_id}: {group_name}.{component_name} must be numeric"
        )
        return 0.0
    if not math.isfinite(value) or not 0 <= value <= maximum:
        result.errors.append(
            f"{candidate_id}: {group_name}.{component_name} must be between 0 and {maximum}"
        )
        return 0.0
    return value


def validate_candidate(candidate: dict[str, Any]) -> ValidationResult:
    result = ValidationResult()
    candidate_id = candidate.get("candidate_id")
    if not isinstance(candidate_id, str) or not candidate_id.strip():
        result.errors.append("candidate_id must be a non-empty string")
        candidate_id = "<missing>"

    for field_name in ("generic_name", "target_intervention", "mechanism_summary"):
        value = candidate.get(field_name)
        if not isinstance(value, str) or not value.strip():
            result.errors.append(f"{candidate_id}: {field_name} must be a non-empty string")

    if candidate.get("research_only") is not True:
        result.errors.append(f"{candidate_id}: research_only must be true")
    if not isinstance(candidate.get("market_approved"), bool):
        result.errors.append(f"{candidate_id}: market_approved must be boolean")

    jurisdictions = candidate.get("approval_jurisdictions")
    if not isinstance(jurisdictions, list) or not all(
        isinstance(value, str) and value.strip() for value in jurisdictions
    ):
        result.errors.append(
            f"{candidate_id}: approval_jurisdictions must be a list of non-empty strings"
        )

    evidence = candidate.get("evidence")
    evidence_types: set[str] = set()
    if not isinstance(evidence, list) or not evidence:
        result.errors.append(f"{candidate_id}: evidence must be a non-empty list")
    else:
        for index, record in enumerate(evidence, start=1):
            if not isinstance(record, dict):
                result.errors.append(f"{candidate_id}: evidence {index} must be an object")
                continue
            evidence_type = record.get("evidence_type")
            citation = record.get("citation")
            supports = record.get("supports")
            if isinstance(evidence_type, str):
                evidence_types.add(evidence_type)
            else:
                result.errors.append(
                    f"{candidate_id}: evidence {index} requires evidence_type"
                )
            if not isinstance(citation, str) or not citation.strip():
                result.errors.append(f"{candidate_id}: evidence {index} requires citation")
            if not isinstance(supports, str) or not supports.strip():
                result.errors.append(f"{candidate_id}: evidence {index} requires supports")

    missing_types = REQUIRED_EVIDENCE_TYPES - evidence_types
    if missing_types:
        result.warnings.append(
            f"{candidate_id}: missing shortlist evidence types: {', '.join(sorted(missing_types))}"
        )

    risks = candidate.get("risks")
    if not isinstance(risks, list) or not risks or not all(
        isinstance(value, str) and value.strip() for value in risks
    ):
        result.errors.append(f"{candidate_id}: risks must be a non-empty list")

    for name, maximum in POSITIVE_COMPONENTS.items():
        _bounded_score(candidate, "positive_scores", name, maximum, result)
    for name, maximum in PENALTY_COMPONENTS.items():
        _bounded_score(candidate, "penalty_scores", name, maximum, result)
    return result


def _eligibility(candidate: dict[str, Any]) -> tuple[bool, tuple[str, ...]]:
    reasons: list[str] = []
    if candidate.get("market_approved") is not True:
        reasons.append("NOT_MARKET_APPROVED")
    if not candidate.get("approval_jurisdictions"):
        reasons.append("NO_APPROVAL_JURISDICTION")
    evidence_types = {
        record.get("evidence_type")
        for record in candidate.get("evidence", [])
        if isinstance(record, dict)
    }
    for evidence_type in sorted(REQUIRED_EVIDENCE_TYPES - evidence_types):
        reasons.append(f"MISSING_{evidence_type.upper()}_EVIDENCE")
    explicit_reason = str(candidate.get("exclusion_reason") or "").strip()
    if explicit_reason:
        reasons.append(explicit_reason)
    return not reasons, tuple(reasons)


def rank_candidates(
    candidates: list[dict[str, Any]],
    weights: dict[str, dict[str, float]],
) -> tuple[list[RankedCandidate], ValidationResult]:
    validation = ValidationResult()
    seen_ids: set[str] = set()
    ranked: list[RankedCandidate] = []

    for candidate in candidates:
        if not isinstance(candidate, dict):
            validation.errors.append("every candidate must be an object")
            continue
        candidate_result = validate_candidate(candidate)
        validation.errors.extend(candidate_result.errors)
        validation.warnings.extend(candidate_result.warnings)
        candidate_id = str(candidate.get("candidate_id") or "")
        if candidate_id in seen_ids:
            validation.errors.append(f"duplicate candidate_id: {candidate_id}")
        seen_ids.add(candidate_id)
        if candidate_result.errors:
            continue

        benefit = sum(
            float(candidate["positive_scores"][name])
            * weights["positive_weights"][name]
            for name in POSITIVE_COMPONENTS
        )
        penalty = sum(
            float(candidate["penalty_scores"][name])
            * weights["penalty_weights"][name]
            for name in PENALTY_COMPONENTS
        )
        eligible, reasons = _eligibility(candidate)
        final_score = max(0.0, min(100.0, benefit - penalty))
        ranked.append(
            RankedCandidate(
                source=candidate,
                eligible=eligible,
                benefit_score=round(benefit, 4),
                risk_penalty=round(penalty, 4),
                final_score=round(final_score, 4),
                exclusion_reasons=reasons,
            )
        )

    ranked.sort(
        key=lambda item: (
            not item.eligible,
            -item.final_score,
            str(item.source["generic_name"]).casefold(),
            str(item.source["candidate_id"]),
        )
    )
    return ranked, validation


def _candidate_row(item: RankedCandidate, rank: int | None) -> dict[str, str]:
    source = item.source
    row = {
        "RANK": str(rank or ""),
        "CANDIDATE_ID": str(source["candidate_id"]),
        "GENERIC_NAME": str(source["generic_name"]),
        "ELIGIBLE_FOR_SHORTLIST": "YES" if item.eligible else "NO",
        "BENEFIT_SCORE": f"{item.benefit_score:.4f}",
        "RISK_PENALTY": f"{item.risk_penalty:.4f}",
        "FINAL_SCORE": f"{item.final_score:.4f}",
        "TARGET_INTERVENTION": str(source["target_intervention"]),
        "EVIDENCE_COUNT": str(len(source["evidence"])),
        "EXCLUSION_REASONS": ";".join(item.exclusion_reasons) or ".",
        "RESEARCH_ONLY": "YES",
    }
    for name in POSITIVE_COMPONENTS:
        row[name.upper()] = str(source["positive_scores"][name])
    for name in PENALTY_COMPONENTS:
        row[name.upper()] = str(source["penalty_scores"][name])
    return row


def build_outputs(
    input_path: Path,
    weights_path: Path,
    output_csv: Path,
    summary_path: Path,
) -> ValidationResult:
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"input schema_version must be {SCHEMA_VERSION}")
    candidates = payload.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("input must contain a non-empty candidates list")

    weights = load_weights(weights_path)
    ranked, validation = rank_candidates(candidates, weights)
    if not validation.valid:
        return validation

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    shortlist_rank = 0
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in ranked:
            rank: int | None = None
            if item.eligible:
                shortlist_rank += 1
                rank = shortlist_rank
            writer.writerow(_candidate_row(item, rank))

    summary = {
        "schema_version": SCHEMA_VERSION,
        "research_only": True,
        "score_is_not_probability_or_clinical_recommendation": True,
        "candidate_count": len(ranked),
        "eligible_count": sum(item.eligible for item in ranked),
        "excluded_count": sum(not item.eligible for item in ranked),
        "required_evidence_types": sorted(REQUIRED_EVIDENCE_TYPES),
        "positive_weights": weights["positive_weights"],
        "penalty_weights": weights["penalty_weights"],
        "formula": "clamp(sum(positive_score*weight)-sum(penalty_score*weight),0,100)",
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return validation


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json", type=Path)
    parser.add_argument("weights_json", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("summary_json", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = build_outputs(
            args.input_json,
            args.weights_json,
            args.output_tsv,
            args.summary_json,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    for warning in result.warnings:
        print(f"WARNING: {warning}")
    for error in result.errors:
        print(f"ERROR: {error}", file=sys.stderr)
    print("VALID" if result.valid else "INVALID")
    return 0 if result.valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
