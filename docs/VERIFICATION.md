# Verification record

## Executed checks

Environment: R 4.5.2, Linux.

Command executed from the repository root:

```sh
Rscript tests/test_core.R
```

Initial publication: **18 checks passed; all 14 R files parsed successfully.**
The supplied-folder publication adds original archived scripts and a public-data
downloader. **All 23 R files parse; the same 18 core checks pass.**
All 35 copied source files were checked against their Drive MD5 checksums,
and all local Markdown links were checked for valid targets.
The main report was then added from a checksum-matching chat attachment,
bringing the total to 36 original files. All 11 newly attached platform plots
were byte-identical to the plots already published in `supplied_project/`.
They are now embedded directly in the root README and platform galleries.
The full RNA-seq DE CSV was subsequently uploaded from a checksum-matching
attachment: 19,130 data rows, seven columns, 2,244,278 bytes. The total is now
37 original files. The other three accompanying RNA-seq CSV attachments matched
the existing copies byte-for-byte and were not duplicated.
The public-input downloader was executed successfully: both GEO inputs
downloaded, matched the supplied Drive MD5 checksums, and the duplicate RNA
input path was restored. The large input bytes are not committed to GitHub.

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
