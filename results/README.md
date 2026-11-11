# Historical supplied tables

`cross_platform_summary.csv` and `overlap_DE_genes.csv` are unchanged source
artifacts, not outputs of a corrected run.

The overlap CSV has 47 unique genes, 46 sign reversals and one same-sign gene
(ANGPTL4). Pearson r calculated from these 47 rounded rows is approximately
-0.872. The summary's -0.438 refers to a different, larger gene set.

The summary's model names and accuracy fields are historical and must not be
used as verified performance. The archived scripts select features globally,
whereas the report and plot labels claim training-only selection; their exact
executable provenance remains unresolved.

New cross-platform outputs are written under `results/rerun/`; they do not
overwrite these artifacts. Full historical platform DE tables and raw count
inputs are not included.
