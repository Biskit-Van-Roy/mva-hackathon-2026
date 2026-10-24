# Track 2 Report — MVA: evidence-constrained readthrough hypothesis
**Participant:** Biskit-Van-Roy (GitHub handle; replace with Hugging Face username in submitted filename)
**Date:** 2026-10-23
**Status:** pre-submission draft; research only; not medical advice

## Executive summary
This project evaluates approved aminoglycoside antibiotics as *hypothesis-generating* candidates to investigate premature termination codon readthrough in the setting of mosaic variegated aneuploidy (MVA). The gene-level disease mechanism centers on BUB1B/BUBR1 and spindle-assembly-checkpoint dysfunction. No evidence of variant-specific BUB1B functional rescue, MVA efficacy, or acceptable therapeutic exposure is established by this analysis. We do not propose administering these antibiotics to a patient.

## 1. Disease mechanism and intervention rationale
BUB1B encodes BUBR1, a key participant in mitotic checkpoint function. Pathogenic biallelic BUB1B variation can produce MVA type 1 and defective chromosome segregation. A proposed nonsense allele could potentially serve as a readthrough substrate, subject to confirmation of the relevant transcript, stop codon, local nucleotide context, transcript abundance and nonsense-mediated decay. A missense allele is mechanistically distinct and cannot itself be rescued by premature-stop readthrough. No case-specific coordinates or sensitive clinical records are included here.

Aminoglycosides have been experimentally observed to promote premature-stop readthrough in certain other genes and model systems, but these observations cannot be generalized to BUB1B or translated into clinical effectiveness without direct work.

## 2. Candidate evidence and critical uncertainties

| Candidate | Repurposing hypothesis | Evidence basis | Missing validation | Main concerns |
|---|---|---|---|---|
| Gentamicin | Investigate readthrough in a BUB1B nonsense-variant model | Primary experimental readthrough in other genes | BUB1B-specific readthrough and functional checkpoint rescue | Renal and auditory/vestibular toxicity; exposure feasibility |
| Amikacin | Same readthrough hypothesis, mechanistically complementary comparator | Primary experimental readthrough in other genes, context dependent | Same | Renal and auditory/vestibular toxicity; exposure feasibility |
| Tobramycin | Comparator for relative readthrough behavior | Primary experimental readthrough observations in other models | Same | Renal and auditory/vestibular toxicity; exposure feasibility |

These medicines are approved for specified bacterial infections in particular jurisdictions, **not for MVA**. Approval status, formulations, authoritative regulatory labels, and safety evidence must be documented for each candidate before describing it as eligible under the software's evidence gate.

## 3. Methods and reproducible software
Public source: https://github.com/Biskit-Van-Roy/mva-hackathon-2026

The Track 2 program validates candidate records and evidence types, checks market approval/jurisdiction, retains exclusion reasons, and calculates an auditable heuristic:

`score = clamp(sum(positive_component * weight) - sum(risk_component * penalty_weight), 0, 100)`

The weight configuration is in `config/track2_scoring.json`. The method does not estimate drug response probabilities or clinical benefit.

Local testing on 2026-10-23: **11 automated tests passed**. The synthetic demonstration included two candidates (one eligible; one excluded). The baseline eligible candidate's score was 67/100.

Sensitivity analysis uses a full factorial adjustment of eight weights by multipliers {0.8, 1.0, 1.2}, for **3^8 = 6,561** scenarios. With this *synthetic* dataset, the sole eligible candidate remained first in 6,561 scenarios. Its final score ranged from 45.6 to 88.4; the excluded example remained at 0. This is primarily a software reproducibility/eligibility check, **not comparative robustness between actual drugs**, because only one candidate is eligible. A test with two eligible candidates was also executed, but its assertions do not establish actual rank reversal.

No patient-specific numerical drug ranking is claimed in this report. No individualized genomic data were shared publicly.

## 4. Proposed validation experiments (future work, not performed)
1. Verify molecular transcript and effect of candidate variants and whether a premature-stop transcript is present.
2. Evaluate nonsense-mediated decay and the local termination-codon sequence context.
3. In controlled preclinical assays, compare protein readthrough and BUBR1 abundance with appropriate positive, negative and toxicity controls.
4. Assess BUBR1 activity, spindle checkpoint integrity, chromosome missegregation and cell viability.
5. Assess dose–response and toxicity in a clinically relevant model without extrapolating laboratory concentrations into a recommended patient regimen.
6. Reassess alternatives if no meaningful biological rescue or acceptable safety margin is observed.

## 5. Strengths and limitations
**Strengths:** explicit evidence-gating; transparent weighting; separate risk penalties; visibility of exclusions; version-controlled code; synthetic public examples; reproducible test results.

**Limitations:** no BUB1B-specific drug-response data; uncertain variant/transcript relevance; indirect mechanistic evidence from other genes; no patient safety/efficacy evidence; no authoritative per-candidate regulatory dossiers incorporated into a real scoring dataset; synthetic sensitivity analysis cannot rank medicines. A computational score must not substitute for expert genetic, pharmacologic, or clinical review.

## 6. References
- ClinGen, BUB1B gene–disease validity curation: https://search.clinicalgenome.org/CCID:004299
- Keeling KM, Bedwell DM (2002), experimental aminoglycoside-mediated premature-stop suppression in other genes: https://pubmed.ncbi.nlm.nih.gov/12072912/
- Howard MT et al. (2004), variation in nonsense-mutation readthrough response: https://pubmed.ncbi.nlm.nih.gov/14991821/
- Floquet et al. (2012), sequence-context considerations: https://pubmed.ncbi.nlm.nih.gov/22479203/

## 7. AI-assisted methods disclosure — MUST COMPLETE BEFORE SUBMISSION
Drafting/code review assistance: OpenAI ChatGPT. **Participant must specify actual product plan/tier, relevant data-handling/training setting, what material was shared, and which claims/code were independently verified.** Do not assert a privacy or retention setting unless confirmed.

## Submission blockers / mandatory review
- [ ] Replace report filename with the actual Hugging Face username, as requested by organizers.
- [ ] Check each cited publication's title, year and relevance against its primary text.
- [ ] Document per-drug authoritative jurisdiction-specific regulatory and safety sources.
- [ ] Have a scientifically qualified reviewer inspect the variant/mechanism reasoning and candidates.
- [ ] Record and upload a 3-minute YouTube/Vimeo pitch; obtain working URL.
- [ ] Verify the public GitHub repository reflects the final report, and ensure no restricted data is exposed.
- [ ] Upload this Markdown or exported PDF to the Hugging Face Submit – Track 2 page with GitHub and video URLs.
- [ ] Save the submission confirmation. A pushed Git commit is **not** a Track 2 submission.
