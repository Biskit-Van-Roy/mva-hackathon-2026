# Track 2 drug-repurposing evidence ranking

This module prioritises approved medicines for laboratory follow-up. It does not
estimate treatment efficacy and must not be used for prescribing or clinical
decision-making.

## Design principles

1. Every component and weight is visible in version-controlled JSON.
2. Benefit and risk are calculated separately.
3. Candidates that fail eligibility stay in the output with explicit reasons.
4. A high mechanism score cannot erase a chromosome-instability or patient-specific risk.
5. Shortlisting requires primary mechanistic, regulatory and safety evidence.
6. Synthetic public examples exercise the pipeline without exposing patient data.

## Eligibility gate

A candidate can enter the shortlist only when it:

- is recorded as market approved;
- has at least one approval jurisdiction;
- includes primary mechanistic evidence;
- includes an authoritative regulatory source;
- includes an authoritative safety source; and
- has no explicit exclusion reason.

Eligibility and score are deliberately separate. An investigational compound can
have a high experimental score but cannot enter an approved-drug repurposing shortlist.

## Score

The default configuration produces a 0–100 research-prioritisation score:

```text
benefit = sum(positive component × positive weight)
risk penalty = sum(risk component × penalty weight)
final = clamp(benefit - risk penalty, 0, 100)
```

Positive components cover mechanistic fit, evidence strength, regulatory fit,
experimental testability and scalability. Penalties cover general safety,
possible worsening of chromosomal instability and patient-specific conflicts.

The component scales are ordinal rubrics, not biological measurements. Before
the final submission, the report must define every score level and include a
sensitivity analysis showing whether conclusions change under reasonable weights.

## Public demonstration

Run the synthetic example:

```bash
python -m mva_hackathon.drug_ranking \
  config/track2_candidates.example.json \
  config/track2_scoring.json \
  artifacts/track2/example_ranking.tsv \
  artifacts/track2/example_summary.json
```

The example names and citations are explicitly synthetic. Real candidate evidence
will be curated separately and must cite primary publications and authoritative
regulatory sources.
