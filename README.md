# RNAseq_MCF7_JQ1_E2

RNA-seq differential-expression analysis of MCF7 cells treated with the BET
inhibitor JQ1, with and without estradiol (E2).

## Data source and sample selection

The eight samples are from [GEO GSE55922](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55922),
“Bromodomain protein BRD4 is required for estrogen receptor-dependent
transcription and enhancer activation” (GSM1348541–GSM1348544 and
GSM1348549–GSM1348552). Their sequencing reads are part of
[SRA study SRP040136](https://www.ncbi.nlm.nih.gov/sra?term=SRP040136) and
[BioProject PRJNA241317](https://www.ebi.ac.uk/ena/browser/view/PRJNA241317).

The selected groups are JQ1 plus E2, JQ1 plus vehicle, control siRNA plus E2,
and control siRNA plus vehicle, with two replicates per group. The study also
contains BRD4-siRNA samples; those are not included in this eight-sample
comparison. Sample labels and run accessions are recorded in
[`data/sample_metadata.csv`](data/sample_metadata.csv).

## Reproduce the analysis

Install R and the required R packages:

```sh
sudo apt-get install r-base-core
Rscript scripts/install_R_packages.R
```

Prepare the gene-count matrix from the recount3 data archive, then run the
analysis:

```sh
Rscript scripts/prepare_counts.R
Rscript analysis.R
```

The source files are from the recount3 SRP040136 archive and use GENCODE v29.
They contain base-pair coverage sums, not ordinary integer read counts.
`prepare_counts.R` applies recount3's documented `targetSize / AUC` scaling
(`targetSize = 4e7`) and rounds the result, matching recount3's
`transform_counts(by = "auc")` method. The small compressed source files and
derived matrix are included under `data/`; the much larger FASTQ files are not
needed. The GEO-supplied normalized-count table is not used as edgeR input.

The analysis filters genes at CPM > 1 in at least two samples, applies TMM
normalization, performs MDS/PCA and expression-distribution QC, estimates
dispersion, and runs an edgeR likelihood-ratio test for the JQ1 main effect
while adjusting for E2/vehicle status. This compares JQ1 with control siRNA
across both estrogen conditions; it does not test a JQ1-by-E2 interaction.
KEGG over-representation uses human Ensembl-to-Entrez mappings and KEGG
pathway membership. No batch correction is applied because the sample metadata
does not provide a batch variable.

## Outputs

`data/gene_counts.csv` contains the derived gene-level count matrix, and
`data/library_scaling_factors.csv` records the per-sample AUC scaling factors.
Analysis tables, session information, and QC/differential-expression plots are
written under `results/`.
