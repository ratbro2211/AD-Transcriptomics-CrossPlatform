# Microarray: supplied historical outputs

These files are preserved from Drive; **they have not been regenerated using
the corrected pipeline**. The report/plot labels state training-only selection,
but the archived script uses global DEGs; supplied ML scores are not independently verified.

## Differential expression

![Supplied microarray volcano](volcano_microarray.png)

[Full supplied DE table—42,179 rows](de_results_microarray.csv)

The full CSV is included directly and matches the original Drive file's checksum.
It is preserved unchanged; the number of result rows is not a unique-gene count.

## Supplied enrichment results

These are source outputs, not newly confirmed enrichment results.

![Supplied microarray GO plot](go_bp_microarray.png)

[Full supplied GO enrichment table](go_results_microarray.csv)

![Supplied microarray KEGG plot](kegg_microarray.png)

[Full supplied KEGG enrichment table](kegg_results_microarray.csv)

## Supplied ML—not independently verified

![Original microarray model comparison](ml_comparison_microarray.png)

![Original microarray feature importance](rf_importance_microarray.png)

[Original score table](ml_model_results_microarray.csv)

The feature-importance rankings are exploratory and do not establish biomarkers.
