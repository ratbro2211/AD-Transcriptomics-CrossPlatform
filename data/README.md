# Input data and provenance

Large/private working inputs are not committed. GEO datasets remain available
from [GSE104704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104704)
and [GSE132903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132903).

## RNA-seq (required before a fresh DESeq2 run)

The GEO series matrix provides metadata, not a usable raw-count matrix.
The submitter's `GSE104704_RNA-Seq_Table.txt.gz` contains DESeq
library-size-adjusted tag counts. **Do not round, truncate, or coerce those
normalized values into raw counts.**

Obtain bona fide gene-level raw read counts from a documented processing of
the study's sequencing reads, or a verified raw-count source. NCBI-computed
RNA-seq counts, if available for this accession, require checking their
processing and exact sample IDs. No availability is assumed by this repository.
FASTQ alignment/quantification is outside this pipeline's scope.

Provide `data/GSE104704/raw_counts.tsv`:

- A tab-delimited matrix, first column named `gene_id`, one unique gene per row.
- Remaining column names must be exact GEO GSM sample accessions.
- Values must be nonnegative integer raw counts, not CPM/TPM/FPKM/normalized counts.
- At minimum include the 12 AD and 10 Old GSM samples; Young samples are excluded.
- Stable Ensembl IDs, Entrez IDs, or valid HGNC symbols are accepted for annotation.
  Do not mix identifier systems. Ambiguous/unmapped symbols are excluded from
  enrichment and cross-platform comparisons, not silently re-labelled.

Also provide `data/GSE104704/counts_provenance.tsv`, one row with columns:

```text
dataset	input_type	source_url	processing
GSE104704	raw_gene_counts	<actual upstream source URL>	<actual reference, quantifier and counting procedure>
```

Replace the bracketed descriptions with actual provenance; do not copy them
as assertions. The script checks the declaration and count values, but **a
declaration plus integer values is not proof of raw-count origin**.

The pipeline downloads GEO metadata and checks exactly 12 AD and 10 Old controls.
It stops before DESeq2 if raw counts/provenance are missing.

## Microarray

The pipeline downloads the processed GSE132903 series matrix using GEOquery,
checks platform GPL10558 and exactly 97 AD / 98 ND samples, and maps probes
through `illuminaHumanv4.db`. It expects log-scale values and performs an
explicit additional quantile normalization for the exploratory DE analysis.
It does not claim to reproduce raw-intensity processing or run `neqc()`.

The ML input uses the source processed log matrix, before this additional
full-cohort quantile normalization. The study submitter's preprocessing may
already have used the full cohort; that limitation is not removable from
processed data and must be reported.
