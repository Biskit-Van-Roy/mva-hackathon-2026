from __future__ import annotations

import gzip
import tempfile
import unittest
import zipfile
from pathlib import Path

from mva_hackathon.intake import build_report


class IntakeTests(unittest.TestCase):
    def test_aggregate_report_detects_grch38_and_hpo(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            vcf_path = root / "synthetic.vcf.gz"
            phenotype_path = root / "phenotype.docx"

            vcf_text = "\n".join(
                [
                    "##fileformat=VCFv4.2",
                    "##reference=GRCh38",
                    "##contig=<ID=chr1,length=248956422>",
                    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic">',
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND01",
                    "chr1\t100\t.\tA\tG\t99\tPASS\tCSQ=G|missense_variant\tGT\t0/1",
                    "",
                ]
            )
            with gzip.open(vcf_path, "wt", encoding="utf-8") as handle:
                handle.write(vcf_text)

            document_xml = (
                '<?xml version="1.0" encoding="UTF-8"?>'
                '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
                "<w:body><w:p><w:r><w:t>HP:0001250 and hp:0001263</w:t></w:r></w:p></w:body>"
                "</w:document>"
            )
            with zipfile.ZipFile(phenotype_path, "w") as archive:
                archive.writestr("word/document.xml", document_xml)

            report = build_report(vcf_path, phenotype_path)

            self.assertTrue(report["vcf"]["grch38_ready"])
            self.assertEqual(report["vcf"]["variant_count"], 1)
            self.assertEqual(report["vcf"]["genotype_counts_first_sample"]["het_or_multiallelic_alt"], 1)
            self.assertTrue(report["vcf"]["annotation_flags"]["vep_csq"])
            self.assertEqual(report["phenotype"]["hpo_ids"], ["HP:0001250", "HP:0001263"])


if __name__ == "__main__":
    unittest.main()
