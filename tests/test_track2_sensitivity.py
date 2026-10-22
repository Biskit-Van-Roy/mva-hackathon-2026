from __future__ import annotations

import json
import tempfile
from pathlib import Path

from mva_hackathon.track2_sensitivity import analyze


def test_sensitivity_counts_every_scenario() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = analyze(
        root / "config/track2_candidates.curated.json",
        root / "config/track2_scoring.json",
        (0.75, 1.0, 1.25),
    )
    assert summary["scenario_count"] == 3**8
    assert sum(summary["winner_counts"].values()) == summary["scenario_count"]
    assert summary["baseline_ranking"][0]["candidate_id"] == "gentamicin-readthrough-probe"


def test_sensitivity_output_is_json_serializable() -> None:
    root = Path(__file__).resolve().parents[1]
    summary = analyze(
        root / "config/track2_candidates.curated.json",
        root / "config/track2_scoring.json",
        (1.0,),
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        path = Path(temporary_directory) / "summary.json"
        path.write_text(json.dumps(summary), encoding="utf-8")
        assert json.loads(path.read_text(encoding="utf-8"))["scenario_count"] == 1
