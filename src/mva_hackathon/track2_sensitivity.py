"""Track 2 synthetic sensitivity analysis; never interpret scores clinically."""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

from .drug_ranking import SCHEMA_VERSION, load_weights, rank_candidates

def analyze(candidates, weights, factors=(0.8, 1.0, 1.2)):
    keys = [(group, key) for group in ("positive_weights", "penalty_weights")
            for key in weights[group]]
    baseline, validation = rank_candidates(candidates, weights)
    if not validation.valid:
        raise ValueError("; ".join(validation.errors))
    base_eligible = [x.source["candidate_id"] for x in baseline if x.eligible]
    eligible_set = set(base_eligible)
    winner_counts = Counter()
    ordering_counts = Counter()
    eligible_changed = 0
    score_ranges = {x.source["candidate_id"]: [float("inf"), float("-inf")] for x in baseline}
    scenarios = 0
    for scale_values in itertools.product(factors, repeat=len(keys)):
        varied = {group: dict(values) for group, values in weights.items()}
        for (group, key), scale in zip(keys, scale_values):
            varied[group][key] = weights[group][key] * scale
        ranked, check = rank_candidates(candidates, varied)
        if not check.valid:
            raise ValueError("; ".join(check.errors))
        ids = tuple(x.source["candidate_id"] for x in ranked if x.eligible)
        if set(ids) != eligible_set:
            eligible_changed += 1
        ordering_counts[ids] += 1
        if ids:
            winner_counts[ids[0]] += 1
        for item in ranked:
            bounds = score_ranges[item.source["candidate_id"]]
            bounds[0] = min(bounds[0], item.final_score)
            bounds[1] = max(bounds[1], item.final_score)
        scenarios += 1
    return {
        "schema_version": SCHEMA_VERSION,
        "research_only": True,
        "synthetic_input_only": True,
        "scenario_count": scenarios,
        "relative_weight_factors": list(factors),
        "vary_all_weights_independently": True,
        "baseline_eligible_order": base_eligible,
        "baseline_excluded": [x.source["candidate_id"] for x in baseline if not x.eligible],
        "eligibility_change_scenarios": eligible_changed,
        "baseline_order_frequency": ordering_counts[tuple(base_eligible)],
        "baseline_order_fraction": ordering_counts[tuple(base_eligible)] / scenarios,
        "winner_counts": dict(winner_counts),
        "final_score_ranges": {k: {"min": v[0], "max": v[1]} for k, v in score_ranges.items()},
        "limitations": [
            "Synthetic candidates only; not a ranking of actual treatments.",
            "Perturbation factors are sensitivity assumptions, not biological uncertainties.",
            "Eligibility rules do not depend on weights by design.",
            "Clamping to 0-100 can hide score differences; investigate ties.",
        ],
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("candidates_json", type=Path)
    parser.add_argument("weights_json", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()
    data = json.loads(args.candidates_json.read_text(encoding="utf-8"))
    if data.get("schema_version") != SCHEMA_VERSION or not isinstance(data.get("candidates"), list):
        parser.error("invalid candidate schema")
    if not all(str(c.get("candidate_id", "")).startswith("synthetic-") for c in data["candidates"]):
        parser.error("this public demo accepts synthetic-* candidate IDs only")
    summary = analyze(data["candidates"], load_weights(args.weights_json))
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"VALID: {summary['scenario_count']} synthetic scenarios")
    print(f"Baseline order fraction: {summary['baseline_order_fraction']:.4f}")

if __name__ == "__main__":
    main()
