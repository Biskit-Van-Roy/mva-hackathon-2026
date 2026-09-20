from __future__ import annotations

import csv
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

from mva_hackathon.drug_ranking import (
    SCHEMA_VERSION,
    build_outputs,
    load_weights,
    rank_candidates,
)


WEIGHTS = {
    "schema_version": SCHEMA_VERSION,
    "positive_weights": {
        "mechanistic_fit": 5,
        "evidence_strength": 4,
        "regulatory_fit": 5,
        "experimental_testability": 5,
        "scalability": 5,
    },
    "penalty_weights": {
        "safety_risk": 4,
        "chromosomal_instability_risk": 6,
        "patient_specific_conflict": 4,
    },
}


def candidate(candidate_id: str, name: str) -> dict:
    return {
        "candidate_id": candidate_id,
        "generic_name": name,
        "market_approved": True,
        "approval_jurisdictions": ["TEST"],
        "target_intervention": "test mechanism",
        "mechanism_summary": "Synthetic unit-test record.",
        "research_only": True,
        "positive_scores": {
            "mechanistic_fit": 5,
            "evidence_strength": 4,
            "regulatory_fit": 3,
            "experimental_testability": 4,
            "scalability": 3,
        },
        "penalty_scores": {
            "safety_risk": 1,
            "chromosomal_instability_risk": 0,
            "patient_specific_conflict": 0,
        },
        "evidence": [
            {"evidence_type": "primary_mechanistic", "citation": "synthetic://p", "supports": "test"},
            {"evidence_type": "regulatory", "citation": "synthetic://r", "supports": "test"},
            {"evidence_type": "safety", "citation": "synthetic://s", "supports": "test"},
        ],
        "risks": ["synthetic risk"],
        "exclusion_reason": "",
    }


class DrugRankingTests(unittest.TestCase):
    def test_penalty_can_reverse_benefit_order(self) -> None:
        safe = candidate("safe", "Safe")
        risky = candidate("risky", "Risky")
        risky["positive_scores"]["evidence_strength"] = 5
        risky["penalty_scores"]["chromosomal_instability_risk"] = 5

        ranked, validation = rank_candidates([risky, safe], {
            "positive_weights": WEIGHTS["positive_weights"],
            "penalty_weights": WEIGHTS["penalty_weights"],
        })

        self.assertTrue(validation.valid)
        self.assertEqual(ranked[0].source["candidate_id"], "safe")
        self.assertGreater(ranked[0].final_score, ranked[1].final_score)

    def test_ineligible_candidate_remains_visible_but_unranked(self) -> None:
        eligible = candidate("eligible", "Eligible")
        excluded = candidate("excluded", "Excluded")
        excluded["market_approved"] = False
        excluded["approval_jurisdictions"] = []
        excluded["exclusion_reason"] = "NEGATIVE_CONTROL"

        ranked, validation = rank_candidates([excluded, eligible], {
            "positive_weights": WEIGHTS["positive_weights"],
            "penalty_weights": WEIGHTS["penalty_weights"],
        })

        self.assertTrue(validation.valid)
        self.assertEqual(len(ranked), 2)
        self.assertTrue(ranked[0].eligible)
        self.assertFalse(ranked[1].eligible)
        self.assertIn("NOT_MARKET_APPROVED", ranked[1].exclusion_reasons)

    def test_build_outputs_is_deterministic_and_auditable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            input_path = root / "candidates.json"
            weights_path = root / "weights.json"
            output_path = root / "ranking.tsv"
            summary_path = root / "summary.json"
            input_path.write_text(json.dumps({
                "schema_version": SCHEMA_VERSION,
                "candidates": [candidate("b", "Beta"), candidate("a", "Alpha")],
            }), encoding="utf-8")
            weights_path.write_text(json.dumps(WEIGHTS), encoding="utf-8")

            result = build_outputs(input_path, weights_path, output_path, summary_path)

            self.assertTrue(result.valid)
            with output_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["GENERIC_NAME"] for row in rows], ["Alpha", "Beta"])
            self.assertEqual([row["RANK"] for row in rows], ["1", "2"])
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertTrue(summary["score_is_not_probability_or_clinical_recommendation"])
            self.assertEqual(summary["eligible_count"], 2)

    def test_rejects_out_of_range_component(self) -> None:
        invalid = candidate("invalid", "Invalid")
        invalid["positive_scores"]["mechanistic_fit"] = 6
        ranked, validation = rank_candidates([invalid], {
            "positive_weights": WEIGHTS["positive_weights"],
            "penalty_weights": WEIGHTS["penalty_weights"],
        })
        self.assertFalse(validation.valid)
        self.assertEqual(ranked, [])
        self.assertTrue(any("between 0 and 5" in error for error in validation.errors))

    def test_weight_schema_must_match_exactly(self) -> None:
        malformed = deepcopy(WEIGHTS)
        del malformed["positive_weights"]["scalability"]
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "weights.json"
            path.write_text(json.dumps(malformed), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must contain exactly"):
                load_weights(path)


if __name__ == "__main__":
    unittest.main()
