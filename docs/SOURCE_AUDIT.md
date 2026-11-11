# Source audit and corrections

## Source

User-provided Drive folder:
https://drive.google.com/drive/folders/13qSBHAletCzn71uysY7_O-ToyOJwi7D8

All eight standalone R scripts in `Scripts/` were read. Duplicate basenames
ending in `_1.R` are the microarray enrichment/ML versions, as confirmed by
their dataset headers. Original checksums and file IDs are recorded in
`SOURCE_MANIFEST.tsv`; the revised scripts are not claimed to be exact copies.

Only the three root PNGs and two requested root CSVs are included as historical
artifacts. Large platform inputs, PDFs, `.DS_Store`, ZIP duplicates and other
unrequested files are not published.

## Problems found

1. Proposed README called GSE132903 Affymetrix and claimed limma-voom. The supplied
   summary identifies Illumina HumanHT-12 v4; the source script uses ordinary
   limma. Corrected code checks GPL10558.
2. Proposed README DEG totals and ML rankings differ from the supplied summary.
   No invented reconciliation is made. Historical totals are attributed and
   unverified model-performance claims are excluded from conclusions.
3. Original RNA script integer-coerced a GEO expression matrix. The series matrix
   is metadata-only, and GEO describes the supplementary table as DESeq
   library-size-adjusted counts. Corrected code requires separately verified raw
   counts and records provenance instead of silently coercing input.
4. Original sample selection guessed the first matching metadata column and
   used broad substring matching. Corrected code uses exact observed title
   patterns and checks expected group sizes; Young RNA samples are excluded.
5. Both original ML scripts selected DEGs using the full cohort before validation.
   The RNA script also used full-cohort fitted VST values. Corrected code starts
   from all raw RNA features or source processed microarray features, with
   training-only feature selection/scaling, and nested LASSO tuning.
6. Enrichment lacked an assay-specific universe. Corrected code supplies mapped
   tested genes and handles insufficient/empty results.
7. Root run commands and comparison relative paths did not agree. All corrected
   entry points now require the repository root and use platform-specific paths.
8. The overlap summary uses |LFC| > 0.5, whereas original DE scripts used > 1.
   One shared config uses > 0.5; both fresh analyses must be rerun together.
9. The supplied 47-gene table has r approximately -0.872, not -0.438. The scatter
   labels -0.438 for 14,962 genes. Corrected comparison exports both scopes
   separately. Full-gene correlation remains unverified.
10. The supplied prose asserted that region explains sign reversals and that
    genes are robust biomarkers. These causal/validation claims are not supported
    by the supplied observational comparison and have been removed.

## What “corrected” means here

Code has been revised to remove the identified unsafe operations and document
assumptions. Base-R regression tests and parsing are not a substitute for running
DESeq2, limma, enrichment and model fitting on verified inputs. No replacement
results have been fabricated. The historical ML artifact is preserved with
explicit warnings, not endorsed.
