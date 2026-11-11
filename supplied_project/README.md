# Original supplied project materials

This directory preserves the user-supplied Drive project, separately from the
corrected analysis code at the repository root.

**Source:** [Google Drive project folder](https://drive.google.com/drive/folders/13qSBHAletCzn71uysY7_O-ToyOJwi7D8).

## Contents and status

- **36 original files copied byte-for-byte**, with their Drive MD5 checksums
  verified before publication.
- Original `Scripts/` includes all eight R scripts and the supplied ZIP bundle.
  These are historical source code, **not the corrected entry points**.
- `RNA_seq/` contains six original figures, three smaller CSVs, and the pipeline
  explanation PDF. Its input metadata are under `GSE104704/`.
- `micro_array/` contains six original figures, the model-results CSV and the
  full supplied GO/KEGG enrichment tables.
- The main project report PDF was added from the uploaded attachment and
  checksum-matched to its Drive original.
- Root-level cross-platform figures, two summary/overlap CSVs and
  `summary_comparison.txt` retain their original names and contents.
- `.DS_Store` files are excluded. The Drive folder's trailing space in
  `micro_array ` is removed from its GitHub directory name.

The complete source inventory, checksums, links and per-file publication status
are recorded in [FILE_MANIFEST.tsv](FILE_MANIFEST.tsv).

## Main report: included directly

[Read the full project report PDF](Transcriptomic%20Analysis%20of%20Alzheimer%27s%20Disease%20A%20Cross-Platform%20Comparative%20Stud.pdf)

The report is preserved unchanged. Its methods describe training-only ML
selection, whereas the archived source scripts use global selection. Some
figure captions also differ from the report's narrative. Inclusion is source
preservation, not endorsement of every numerical or scientific claim.

## Large private tables: linked, not copied

The authenticated Drive media connection rejected these files because its
response limit is 1,000,000 bytes. Anonymous download links returned login HTML,
not the file contents. No Drive sharing permissions were changed.

| Supplied file | Size | Original source |
|---|---:|---|
| RNA-seq `all_DE_results_rnaseq.csv` | 2,244,278 bytes | [Open full RNA-seq DE table](https://drive.google.com/file/d/1trEFMMjy6qBOh42OX-t7Da5jjaeQlBv4/view?usp=drivesdk) |
| Microarray `de_results_microarray.csv` | 5,065,409 bytes | [Open full microarray DE table](https://drive.google.com/file/d/1y3Lwwq5G4A3zQ2CaSsODNTRUI5CckOSc/view?usp=drivesdk) |

These links may require the source owner's Drive permission. **The two tables
are not present as GitHub blobs.** To include them directly later, upload those
exact files here on GitHub or supply them as chat attachments. No fake result
CSV or placeholder PDF is used in their place.

## Large public GEO inputs: reproducible downloads

The supplied microarray matrix and both copies of the RNA-seq supplementary
table match public GEO files by MD5 checksum. To retrieve the same bytes under
their original folder paths, run from the repository root:

```sh
Rscript supplied_project/download_public_inputs.R
```

This downloads the 41,956,695-byte microarray matrix and 5,452,297-byte RNA-seq
supplementary table, checks MD5, and places a copy of the RNA table at each of its
two original locations. These large files are **downloadable inputs, not
committed GitHub blobs**.

The RNA supplementary table remains normalized/library-size-adjusted data.
Publishing or downloading it does not make it valid raw-count input for the
corrected DESeq2 pipeline.

## Interpretation warning

These are historical, supplied artifacts—not newly rerun results. The original
ML scripts perform global feature selection before validation, while the report
and plot labels describe training-only selection. The exact executable provenance
for those figures remains unresolved. Their scores, predictions and feature
rankings are not independently verified performance or biomarker evidence.

The source prose and report can contain claims that the corrected README does
not endorse. In particular, sign reversals alone do not establish a causal
regional explanation. Preserving source files is not validation of their claims.

## Figure galleries

- [RNA-seq supplied figures and tables](RNA_seq/README.md)
- [Microarray supplied figures and tables](micro_array/README.md)

Use the repository-root analysis scripts for corrected code. Do not run the
archived scripts and expect the corrected methodology or paths.
