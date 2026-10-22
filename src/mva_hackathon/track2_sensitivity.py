"""Deterministic weight sensitivity analysis for the Track 2 ranking."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from collections import Counter
from pathlib import Path

from .drug_ranking import load_weights, rank_candidates


def analyze(candidates_path: Path, weights_path: Path, factors: tuple[float, ...]) -> dict:
    payload = json.loads(candidates_path.read_text(encoding="utf-8"))
    candidates = payload["candidates"]
    baseline = load_weights(weights_path)
    component_keys = [
        *(('positive_weights', key) for key in baseline['positive_weights']),
        *(('penalty_weights', key) for key in baseline['penalty_weights']),
    ]
    winner_counts: Counter[str] = Counter()
    order_counts: Counter[tuple[str, ...]] = Counter()
    scenario_count = 0

    for multipliers in itertools.product(factors, repeat=len(component_keys)):
        scenario_count += 1
        perturbed = {
            "positive_weights": dict(baseline["positive_weights"]),
            "penalty_weights": dict(baseline["penalty_weights"]),
        }
        for (group, component), multiplier in zip(component_keys, multipliers, strict=True):
            perturbed[group][component] *= multiplier
        ranked, validation = rank_candidates(candidates, perturbed)
        if not validation.valid:
            raise ValueError("candidate validation failed during sensitivity analysis")
        eligible_order = tuple(
            item.source["candidate_id"] for item in ranked if item.eligible
        )
        if eligible_order:
            winner_counts[eligible_order[0]] += 1
            order_counts[eligible_order] += 1

    baseline_ranked, baseline_validation = rank_candidates(candidates, baseline)
    if not baseline_validation.valid:
        raise ValueError("candidate validation failed for baseline weights")
    baseline_rows = [
        {
            "candidate_id": item.source["candidate_id"],
            "generic_name": item.source["generic_name"],
            "eligible": item.eligible,
            "benefit_score": item.benefit_score,
            "risk_penalty": item.risk_penalty,
            "final_score": item.final_score,
            "exclusion_reasons": list(item.exclusion_reasons),
        }
        for item in baseline_ranked
    ]
    return {
        "research_only": True,
        "factors": list(factors),
        "scenario_count": scenario_count,
        "baseline_ranking": baseline_rows,
        "winner_counts": dict(sorted(winner_counts.items())),
        "winner_fractions": {
            key: round(value / scenario_count, 6)
            for key, value in sorted(winner_counts.items())
        },
        "distinct_eligible_orders": len(order_counts),
        "most_common_eligible_orders": [
            {"order": list(order), "count": count}
            for order, count in order_counts.most_common(10)
        ],
    }


def write_outputs(summary: dict, json_path: Path, tsv_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["candidate_id", "winner_count", "winner_fraction"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for candidate_id, count in summary["winner_counts"].items():
            writer.writerow({
                "candidate_id": candidate_id,
                "winner_count": count,
                "winner_fraction": summary["winner_fractions"][candidate_id],
            })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidates", type=Path)
    parser.add_argument("weights", type=Path)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("winners_tsv", type=Path)
    parser.add_argument("--factors", type=float, nargs="+", default=(0.75, 1.0, 1.25))
    args = parser.parse_args()
    factors = tuple(args.factors)
    if not factors or any(value <= 0 for value in factors):
        parser.error("all sensitivity factors must be positive")
    summary = analyze(args.candidates, args.weights, factors)
    write_outputs(summary, args.summary_json, args.winners_tsv)
    print(f"VALID scenarios={summary['scenario_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
