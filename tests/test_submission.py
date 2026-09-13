from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from mva_hackathon.submission import COLUMNS, build_csv, validate_csv


class SubmissionTests(unittest.TestCase):
    def test_build_sorts_by_epcr_and_validates(self) -> None:
        candidates = {
            "candidates": [
                {
                    "chrom_1": "chr2",
                    "pos_1": 200,
                    "ref_1": "a",
                    "alt_1": "g",
                    "epcr": 0.2,
                    "finding_type": "secondary",
                },
                {
                    "chrom_1": "chr1",
                    "pos_1": 100,
                    "ref_1": "c",
                    "alt_1": "t",
                    "chrom_2": "chr1",
                    "pos_2": 150,
                    "ref_2": "g",
                    "alt_2": "a",
                    "epcr": 0.9,
                    "finding_type": "primary",
                },
            ]
        }
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            input_path = root / "candidates.json"
            output_path = root / "predictions.csv"
            input_path.write_text(json.dumps(candidates), encoding="utf-8")

            result = build_csv(input_path, output_path)

            self.assertTrue(result.valid)
            self.assertTrue(validate_csv(output_path).valid)
            with output_path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(list(rows[0]), COLUMNS)
            self.assertEqual(rows[0]["epcr"], "0.9")
            self.assertEqual(rows[0]["ref_1"], "C")
            self.assertEqual(rows[0]["proband_id"], "PROBAND01")

    def test_rejects_partial_pair(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "invalid.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=COLUMNS)
                writer.writeheader()
                writer.writerow(
                    {
                        "proband_id": "PROBAND01",
                        "chrom_1": "chr1",
                        "pos_1": "100",
                        "ref_1": "A",
                        "alt_1": "G",
                        "chrom_2": "chr1",
                        "epcr": "0.8",
                        "finding_type": "primary",
                    }
                )

            result = validate_csv(path)

            self.assertFalse(result.valid)
            self.assertTrue(any("partially populated" in error for error in result.errors))


if __name__ == "__main__":
    unittest.main()
