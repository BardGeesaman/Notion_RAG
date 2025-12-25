# Structural Alerts Guide

## Overview
Structural alert filtering is a fast, rule-based screen to flag common substructures associated with assay interference, reactivity, or general “medchem liabilities”.

Common uses:
- **Compound registration**: catch obvious liabilities early.
- **HTS triage**: prioritize follow-up on cleaner chemotypes.
- **Lead optimization**: track whether changes introduce new liabilities.

## Filter Types

### PAINS (≈480 patterns)
PAINS (Pan-Assay INterference compoundS) patterns are associated with promiscuous assay behavior (e.g., aggregation, redox cycling, covalent reactivity). A PAINS hit does not guarantee a compound is “bad”, but it is a strong warning to validate with orthogonal assays and careful controls.

### Brenk (≈105 patterns)
Brenk filters capture “unwanted substructures” that frequently correlate with poor developability, reactivity, or other liabilities. These are medium-severity flags that typically require context-aware review.

### Lilly (≈20 patterns, MVP)
The Lilly filter here is an MVP, curated SMARTS subset targeting common liabilities such as reactive electrophiles (acyl halides, Michael acceptors), toxicophores (nitro aromatics, anilines), and unstable motifs (peroxides, epoxides). It is intended as a quick triage list rather than a complete representation of internal Lilly rules.

## Traffic Light System
- **RED**: Any PAINS match (high severity). Treat as likely assay-interfering; validate with controls/orthogonal assays.
- **YELLOW**: Brenk and/or Lilly matches only (medium severity). Proceed with caution; context matters.
- **GREEN**: No alerts detected by the enabled filters.

## API Usage

### Single compound
`POST /api/alerts/check`

Request:
- `smiles`: string
- `filters`: optional list of filter names (`["pains","brenk","lilly"]`). Omit/None runs all.

### Batch (up to 100)
`POST /api/alerts/batch`

Request:
- `smiles_list`: list of SMILES (max 100 per request)
- `filters`: optional list of filter names

### List available filters
`GET /api/alerts/filters`

Response includes:
- `name`, `pattern_count`, `description`

## Limitations
- **Flags are not verdicts**: Many “flagged” compounds can be valid in context (assay type, concentration, controls, mechanism).
- **Green is not a guarantee**: Filters are incomplete by design and do not cover all liabilities.
- **ML toxicity prediction is deferred**: This MVP is structural/rule-based only; learned toxicity models are planned for a later phase.



