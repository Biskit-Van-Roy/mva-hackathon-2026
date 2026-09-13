# MVA Hackathon 2026 — Track 1 Variant Identification Report

**Participant:** Roy-Hidalgo  
**Proband:** PROBAND01  
**Reference assembly:** GRCh38  
**Report date:** 13 September 2026  
**Submission file:** `Roy-Hidalgo_hpo-mva-v2.csv`  
**Reproducibility repository:** [Biskit-Van-Roy/mva-hackathon-2026](https://github.com/Biskit-Van-Roy/mva-hackathon-2026)  
**Audited code commit:** [`05d67bfc71fac8d30bc9faa64ddbe4b458f5d7df`](https://github.com/Biskit-Van-Roy/mva-hackathon-2026/commit/05d67bfc71fac8d30bc9faa64ddbe4b458f5d7df)

> This is a research hackathon result, not a clinical diagnosis. Clinical use would require orthogonal confirmation, professional variant interpretation, and segregation or molecular phasing.

## Executive summary

The top-ranked hypothesis is an unphased compound-heterozygous pair in **BUB1B**, a gene with an established autosomal-recessive relationship to mosaic variegated aneuploidy syndrome 1 (MVA1). The pair comprises:

1. **chr15:40209701 T>G**, `NM_001211.6:c.2210T>G`, `p.Leu737Ter`; and
2. **chr15:40220612 T>G**, `NM_001211.6:c.3006T>G`, `p.Asn1002Lys`.

Both variants were heterozygous and passed the prespecified heuristic genotype-quality checks. The nonsense variant is classified by ClinVar as pathogenic/likely pathogenic. The missense variant was absent from the population and clinical fields available in the local VEP annotation and was predicted deleterious by SIFT and probably damaging by PolyPhen. BUB1B ranked first both in an unbiased lane and in a separate MVA-aware lane. Thus, the result was not created by the explicit MVA-gene bonus.

The pair remains **unphased**: the single-sample VCF reports `GT=0/1` for both records, with missing `PGT` and `PID`. No parental samples, BAM, CRAM, or FASTQ files were available for phasing. The submission therefore represents a strong causal hypothesis rather than confirmed biallelic inheritance.

## Data and privacy

Analysis used the gated single-proband VCF and a structured phenotype profile. Raw sequencing reads were not used. No raw VCF, FASTQ, BAM, CRAM, or clinical document was uploaded to the public repository. Patient-level inputs and generated candidate, ranking, and submission files remained under Git-ignored `data/` and `artifacts/` paths.

Before the initial push, tracked paths were audited for genomic files and restricted-data directories. No sensitive paths were tracked. The final submission CSV and its checksum were independently confirmed as ignored by Git.

The repository contains reproducible code, tests, configuration examples, and methodological documentation. It must be publicly accessible when the Track 1 entry is submitted, as required by the official interface.

## Phenotype representation

Seven proband HPO terms and one family-history term were manually curated and validated against Human Phenotype Ontology release `v2026-09-01`:

| Role           | HPO term   | Label                          |
| -------------- | ---------- | ------------------------------ |
| Proband        | HP:0002859 | Rhabdomyosarcoma               |
| Proband        | HP:0000121 | Nephrocalcinosis               |
| Proband        | HP:0004322 | Short stature                  |
| Proband        | HP:0001508 | Failure to thrive              |
| Proband        | HP:0003202 | Skeletal muscle atrophy        |
| Proband        | HP:0001622 | Premature birth                |
| Proband        | HP:0001518 | Small for gestational age      |
| Family history | HP:0200067 | Recurrent spontaneous abortion |

Reference files were obtained from the [HPO v2026-09-01 release](https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01):

- `hp.obo`: SHA-256 `93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22`
- `genes_to_phenotype.txt`: SHA-256 `507a17bff9c49e6329fbd88b1f91734fa7e9fdea5c06304044c0892eb6ab248c`

The ranking used ontology ancestor propagation and an information-content similarity score. Of 490 candidate gene symbols, 202 had HPO annotations in the selected reference.

## Variant processing

The supplied GRCh38 VCF was functionally annotated with Ensembl Variant Effect Predictor release 116. The retained annotation exposed canonical/MANE transcript selection, predicted consequence and impact, HGVS expressions, population-frequency fields, ClinVar significance, SIFT, and PolyPhen.

The broad rare-disease filter retained:

- `HIGH` or `MODERATE` impact annotations;
- splice-region annotations;
- functional candidates with missing maximum population allele frequency or `MAX_AF <= 0.01`; and
- exact pathogenic or likely-pathogenic ClinVar annotations as a clinical override.

Conflicting ClinVar assertions were recorded but did not act as a pathogenic override. One best passing annotation was selected per VCF record without deleting the underlying source record. The process yielded 607 candidate records.

## Genotype quality control

Quality control was non-destructive: all 607 candidates were retained and labeled. Default thresholds were:

- minimum depth: `DP >= 10`;
- minimum genotype quality: `GQ >= 20`;
- heterozygous allele balance: `0.20–0.80`; and
- homozygous-alternate allele balance: `>= 0.80`.

Records were labeled `PASS_HEURISTIC_QC`, `REVIEW_LOW_CONFIDENCE`, or `REVIEW_MISSING_METRICS`. These labels are prioritization aids, not diagnostic classifications.

## Inheritance modeling and ranking

Variants were grouped by gene symbol. Blank symbols were separated by Ensembl gene identifier or locus so that unrelated unannotated records could not form a false compound-heterozygous hypothesis.

For each gene, the unbiased score combined:

`10 × HPO similarity + 1 × family-history similarity + top variant evidence + 0.5 × second variant evidence + inheritance signal`

Variant evidence integrated consequence severity, VEP impact, rarity, ClinVar status, SIFT/PolyPhen predictions, and genotype-QC tier. The inheritance component assigned a transparent heuristic signal to a passing homozygous-alternate candidate, two or more passing heterozygous candidates, or a single passing heterozygous candidate.

A second MVA-aware lane added a fixed three-point bonus to configured mechanism genes (`BUB1B`, `CEP57`, and `TRIP13`). Both lanes retained every candidate and exposed their component scores. BUB1B ranked first in both lanes:

| Metric                      |                           Result |
| --------------------------- | -------------------------------: |
| Candidate variants in BUB1B |                                2 |
| Passing QC                  |                                2 |
| Passing heterozygous        |                                2 |
| HPO similarity              |                         0.722994 |
| Top variant evidence        |                             16.5 |
| Second variant evidence     |                              8.0 |
| Unbiased gene score         |                        30.229935 |
| Unbiased rank               |                                1 |
| MVA-aware rank              |                                1 |
| Inheritance label           | UNPHASED_COMPOUND_HET_HYPOTHESIS |

The submitted EPCR of `0.90` was assigned from the convergence of phenotype fit, the known recessive disease mechanism, a pathogenic/likely-pathogenic loss-of-function allele, a rare predicted-damaging second allele, strong genotype metrics, and first place in both ranking lanes. It is a transparent expert heuristic, not an empirically calibrated posterior probability.

## Ranked finding

| Rank | Gene  | Variant 1          | Variant 2          | Finding type | EPCR |
| ---: | ----- | ------------------ | ------------------ | ------------ | ---: |
|    1 | BUB1B | chr15:40209701 T>G | chr15:40220612 T>G | primary      | 0.90 |

### Variant 1: BUB1B p.Leu737Ter

- GRCh38: `chr15:40209701 T>G`
- MANE Select: `NM_001211.6:c.2210T>G`; `NP_001202.5:p.Leu737Ter`
- Consequence: stop gained; VEP impact `HIGH`
- Genotype: `0/1`; `AD=21,25`; `DP=46`; `GQ=99`; allele balance `0.5435`
- Population evidence: local maximum AF `9.982×10⁻5`; gnomAD exome AF `7.867×10⁻5`; gnomAD genome AF `3.286×10⁻5`
- Clinical evidence: ClinVar Variation ID 533901, aggregate pathogenic/likely pathogenic with multiple submitters and no conflicts

The [ClinVar record](https://www.ncbi.nlm.nih.gov/clinvar/variation/533901/) confirms the GRCh38 coordinate, MANE HGVS, nonsense consequence, and pathogenic/likely-pathogenic aggregate classification. ClinVar also cautions that population frequency at this site may have poor-quality metrics; frequency was therefore supportive but not decisive.

### Variant 2: BUB1B p.Asn1002Lys

- GRCh38: `chr15:40220612 T>G`
- MANE Select: `NM_001211.6:c.3006T>G`; `NP_001202.5:p.Asn1002Lys`
- Consequence: missense; VEP impact `MODERATE`
- Genotype: `0/1`; `AD=15,13`; `DP=28`; `GQ=99`; allele balance `0.4643`
- Population and clinical fields: no existing-variant identifier, population AF, or ClinVar classification in the local annotation
- In-silico evidence: SIFT `deleterious (0.01)`; PolyPhen `probably damaging (0.997)`

This second allele is a candidate variant of uncertain clinical significance in the present analysis. Computational predictions alone are insufficient to establish pathogenicity. Its contribution depends critically on independent confirmation, segregation, phase, and additional clinical or functional interpretation.

### Gene–disease and inheritance evidence

[ClinGen's BUB1B curation](https://search.clinicalgenome.org/kb/gene-dosage/HGNC:1149) identifies BUB1B as associated with an autosomal-recessive phenotype and summarizes biallelic BUB1B findings in MVA1. It also lists `NM_001211.6` as the MANE Select transcript. This supports compound-recessive modeling but does not prove that this particular missense allele is disease-causing.

## Phase assessment

The VCF contains phase-related FORMAT definitions (`PGT`, `PID`), but both candidate records have `PGT=.` and `PID=.`. Their genotypes are unphased `0/1`, not phased `0|1` or `1|0`. A single-sample VCF cannot establish parental origin. The known hypomorphic BUB1B allele `rs576524605` was not recorded in the variant-only VCF; this absence cannot be interpreted as a confident homozygous-reference call without coverage or gVCF evidence.

The current evidence does not distinguish:

- **trans:** one candidate allele on each homolog, consistent with compound recessive inheritance; from
- **cis:** both candidate alleles on the same homolog, leaving the other copy unaffected.

The next diagnostic step would be targeted confirmation and parental segregation. Read-backed phasing could be attempted if aligned or raw reads become available and the physical distance/read structure permits it. The [ClinGen PM3 guidance](https://www.clinicalgenome.org/docs/pm3-recommendation-for-in-trans-criterion-pm3-version-1.0/) is relevant when evaluating observations in trans for recessive disease.

## Alternative-candidate review

Competing high-ranked genes were manually reviewed to test whether BUB1B's lead resulted from a ranking artifact.

- A homozygous PEX5 intronic/microsatellite deletion initially received a high VEP consequence. Exact ClinVar review showed predominantly benign/likely-benign evidence with conflicting aggregate classification, so it was strongly deprioritized.
- A two-variant HERC2 hypothesis contained one variant of uncertain significance and one likely-benign, relatively more frequent allele on a non-MANE/NMD transcript; it was deprioritized.
- A heterozygous loss-of-function finding in another disease-associated gene was treated as a possible incidental/secondary finding, not as an explanation for the recessive MVA phenotype. It was not included in the Track 1 submission pending consent-aware clinical handling.

Only the primary BUB1B pair was submitted. This avoided adding low-confidence rows merely to fill the ten-row allowance.

## Validation and reproducibility

The following controls were completed:

- all 607 candidate records were present in the V2 variant ranking;
- 490 mapped symbols and 521 independent gene groups were represented;
- 36 records without a symbol were no longer collapsed into one false `UNMAPPED` gene;
- BUB1B remained rank 1 before and after the explicit MVA bonus;
- SHA-256 verification passed for the ranked-variant TSV, ranked-gene TSV, and summary JSON;
- SHA-256 verification passed for the final submission CSV;
- submission and intake unit tests passed (`3 passed` total); and
- the public Git tree was audited to exclude restricted genomic and derived patient files.

The auditable implementation is in:

- `scripts/filter_vep_candidates.py`
- `scripts/summarize_candidate_qc.py`
- `scripts/rank_hpo_candidates.py`
- `src/mva_hackathon/submission.py`
- `tests/test_submission.py`
- `tests/test_intake.py`

The repository commit recorded at the top of this report contains the code state used for this submission. Private derived outputs are intentionally excluded, but their schemas, checks, and generation logic are documented.

## Limitations

1. Phase is unknown, and no parental genotypes are available.
2. No BAM, CRAM, or FASTQ data were available for read-level confirmation, local realignment, coverage assessment, mosaicism analysis, or structural-variant rescue.
3. The second BUB1B allele lacks a local ClinVar classification and is supported partly by in-silico predictions.
4. The rare-disease filter and QC thresholds are heuristic and may miss low-quality true variants, deep intronic variants, regulatory variants, repeat-associated variants, CNVs, or complex SVs.
5. HPO similarity depends on database annotations and the completeness of the provided phenotype.
6. EPCR is a prioritization estimate rather than a calibrated clinical probability.
7. A negative result for the known hypomorphic allele in a variant-only VCF is not equivalent to a confident reference call.

## Artificial-intelligence assistance disclosure

OpenAI ChatGPT Plus was used to assist with code drafting, debugging, command construction, privacy auditing, and synthesis of the methods report. The ChatGPT data control **“Improve the model for everyone” was disabled** during this work. No raw genomic file or clinical document was uploaded to ChatGPT; selected derived summaries and candidate-level outputs were reviewed interactively. All commands were run by the participant in the local environment, and key outputs were checked with unit tests, row counts, manual database review, and SHA-256 verification. AI output was treated as assistance, not as clinical interpretation.

## Required acknowledgement

_"This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA_
_Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize_
_sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed_
_their data and their story to advance research into this rare disease. We acknowledge their trust in making this_
_Hackathon possible."_

## Official submission references

- [Track 1 submission specification](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py)
- [Official evaluator](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)
- [Official Hackathon Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py)
- [MVA Hackathon gated dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)

## Pre-submission checklist

- [x] Re-open the official rules and insert the exact required acknowledgement above.
- [ ] Make the GitHub repository public and verify its URL without authentication.
- [ ] Confirm the report's commit hash matches the final public code revision.
- [ ] Validate `Roy-Hidalgo_hpo-mva-v2.csv` locally one final time.
- [ ] Upload the CSV and this Markdown report through the official Track 1 interface.
- [ ] Record the submission receipt without publishing restricted source data.
- [ ] Delete restricted data and derived datasets within the required post-hackathon window and email the organizers to confirm deletion.
