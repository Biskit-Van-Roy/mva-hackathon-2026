"""Tests for research-only synthetic Track 2 sensitivity results."""
import copy
import json
from pathlib import Path

from mva_hackathon.drug_ranking import load_weights
from mva_hackathon.track2_sensitivity import analyze

ROOT = Path(__file__).resolve().parents[1]


def example():
    candidates = json.loads((ROOT / "config/track2_candidates.example.json").read_text())["candidates"]
    weights = load_weights(ROOT / "config/track2_scoring.json")
    return candidates, weights


def test_full_factorial_reproducible():
    candidates, weights = example()
    result = analyze(candidates, weights)
    assert result["scenario_count"] == 3 ** 8 == 6561
    assert result["eligibility_change_scenarios"] == 0
    assert result["baseline_eligible_order"] == ["synthetic-readthrough-probe"]
    assert result["winner_counts"]["synthetic-readthrough-probe"] == 6561


def test_two_eligible_candidates_can_reverse_order():
    candidates, weights = example()
    competitor = copy.deepcopy(candidates[0])
    competitor["candidate_id"] = "synthetic-competitor"
    competitor["generic_name"] = "SYNTHETIC_COMPETITOR"
    # Different trade-offs create score/rank variation under changed weights.
    competitor["positive_scores"] = {
        "mechanistic_fit": 3, "evidence_strength": 5, "regulatory_fit": 3,
        "experimental_testability": 4, "scalability": 3
    }
    competitor["penalty_scores"] = {
        "safety_risk": 1, "chromosomal_instability_risk": 0,
        "patient_specific_conflict": 2
    }
    result = analyze([candidates[0], competitor], weights)
    assert result["scenario_count"] == 6561
    assert result["eligibility_change_scenarios"] == 0
    assert len(result["baseline_eligible_order"]) == 2
    assert sum(result["winner_counts"].values()) == 6561
    assert 0 <= result["baseline_order_fraction"] <= 1


def test_score_ranges_are_ordered():
    candidates, weights = example()
    result = analyze(candidates, weights, factors=(0.9, 1.1))
    assert result["scenario_count"] == 2 ** 8
    for bounds in result["final_score_ranges"].values():
        assert 0 <= bounds["min"] <= bounds["max"] <= 100
