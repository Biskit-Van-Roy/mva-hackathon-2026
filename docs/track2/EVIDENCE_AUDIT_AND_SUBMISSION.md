# Track 2 — Evidence audit and pre-submission plan

Status: research hypothesis; not a treatment recommendation.
As of 2026-10-23. Keep case-specific genomic and clinical data outside this public repository.

## Supported mechanistic interpretation
- BUB1B encodes BUBR1, a spindle-assembly-checkpoint protein. Gene–disease association does not itself validate a particular variant or therapeutic rescue.
- Readthrough of premature termination codons is sequence-context dependent.
- Keeling & Bedwell (2002; PMID 12072912) tested gentamicin, amikacin and tobramycin using constructs involving IDUA and P53; this is **indirect evidence**, not BUB1B-specific evidence.
- Howard et al. (2004; PMID 14991821) tested mutations in dystrophin; some responded and others were refractory. Do not extrapolate percentages to BUB1B.

## Candidate audit

| Candidate | What is supported | Crucial missing data | Interpretation |
|---|---|---|---|
| Gentamicin | Experimental readthrough in unrelated nonsense variants | BUB1B allele-specific readthrough, functional rescue, exposure and toxicity feasibility | Hypothesis only |
| Amikacin | Context-dependent readthrough, sometimes comparable to gentamicin in experimental systems | Same missing BUB1B evidence | Hypothesis only |
| Tobramycin | Readthrough observed in model systems, often lower than gentamicin in Keeling & Bedwell | Same missing BUB1B evidence | Comparator, not validated therapy |

All three are aminoglycosides; substantial toxicity considerations include kidney and auditory/vestibular damage. Approval as an antibiotic must not be described as approval for MVA.

## Critical scientific gates before final ranking
1. Independently confirm reference genome, reference transcript and effect of each relevant variant. Do not publish restricted coordinates.
2. Determine if the purported nonsense variant yields an mRNA substrate amenable to readthrough; consider nonsense-mediated decay.
3. A missense variant cannot be repaired by readthrough of a premature stop; distinguish independent allele mechanisms.
4. Check any drug-specific primary evidence, regulatory approval jurisdiction, and authoritative safety information separately.
5. Define a level-by-level ordinal rubric for each scoring component before entering real scores.
6. Score unknown evidence as unknown, not zero proven risk. Consider explicit exclusion or missing-data reporting.
7. Include full sensitivity analysis (weight perturbations, eligibility stability and rank reversal) and a no-evidence baseline.
8. Keep risk scores separate from benefit; do not present 0–100 heuristic scores as response probabilities.
9. Propose **preclinical** assays: protein abundance, checkpoint function, chromosome missegregation, viability and toxicity. No patient dosing or clinical recommendation.
10. Distinguish hypotheses, data-derived observations and confirmed findings in final report.

## Sources for the audit
- Keeling KM, Bedwell DM (2002). https://pubmed.ncbi.nlm.nih.gov/12072912/
- Howard MT et al. (2004). https://pubmed.ncbi.nlm.nih.gov/14991821/
- ClinGen BUB1B gene–disease curation: https://search.clinicalgenome.org/CCID:004299

## Reproducibility and delivery checklist
- [x] Track 2 ranking implementation merged into main (PR #1).
- [x] 8 automated tests passed locally, as reported on 2026-10-23.
- [x] Synthetic ranking generated and validated locally.
- [x] Preliminary scientific comparison documented.
- [ ] Collect authoritative per-drug regulatory and safety documents and record jurisdiction.
- [ ] Curate primary experimental evidence and explicit scoring rubric.
- [ ] Run and save sensitivity analyses; report whether candidate order changes.
- [ ] Produce PDF report, at most the official required format/length.
- [ ] Record a <=3-minute pitch if required by submission form.
- [ ] Verify the exact Track 2 submission fields in the live portal, submit, and retain receipt/confirmation.

Deadline listed on official page: 2026-10-24 23:59 UTC.
