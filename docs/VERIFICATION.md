# Verification record

## Executed checks

Environment: R 4.5.2, Linux.

Command executed from the repository root:

```sh
Rscript tests/test_core.R
```

Result: **18 checks passed; all 14 R files parsed successfully.**

Checks cover:

- Rejection of fractional/negative count input and duplicate sample IDs.
- Reference level and expected-group-size guards.
- Training-only feature ranking/centering/scaling, including a test that
  perturbing held-out feature values cannot change fitted preprocessing.
- Stratified inner-fold construction.
- Separation of all-common-gene and overlap-only comparison scopes.
- Arithmetic on the supplied overlap file: 47 unique genes, one concordant
  gene ANGPTL4, and Pearson r approximately -0.872 over those rounded rows.

The three supplied PNGs were downloaded; the LFC and ML figures were visually
inspected for interpretation. Original supplied bytes are preserved and their
SHA-256 checksums recorded in `SOURCE_MANIFEST.tsv`.

The public GEO series matrices were inspected to establish title-based sample
patterns and the GSE104704 supplementary normalized-count description.

## Not executed / not established

- No end-to-end DESeq2 or limma analysis was run.
- No new GO or KEGG enrichment results were computed.
- No ML model was fit or performance estimate recomputed.
- CRAN/Bioconductor model and data APIs have not been integration-tested in this
  environment; only base R packages were available during the recorded tests.
- The full input for the historical 14,962-gene correlation was not reproduced.
- A verified raw-count matrix for GSE104704 was not supplied. The RNA pipeline
  therefore deliberately requires a user-provided documented raw-count source.
- Clinical validation, independent external validation, causal regional effects
  and biomarker robustness are not established.

Passing parsing and unit checks establishes only the properties actually
tested. Future biological runs must record package versions, inputs,
provenance, sample selection, diagnostics and outputs before upgrading status.
