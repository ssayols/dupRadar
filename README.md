[![](https://img.shields.io/badge/Altmetric-9-blue.svg)](https://www.altmetric.com/details/12855878)

# dupRadar

Assessment of duplication rates in RNA-Seq datasets.

## Introduction to dupRadar

PCR clonal artefacts originating from NGS library preparation can affect both 
genomic as well as RNA-Seq applications when protocols are pushed to their 
limits. In RNA-Seq however the artifactual reads are not easy to tell apart 
from normal read duplication due to natural over-sequencing of highly expressed 
genes. Especially when working with little input material or single cells 
assessing the fraction of duplicate reads is an important quality control step 
for NGS data sets. Up to now there are only tools to calculate the global 
duplication rates that do not take into account the effect of gene expression 
levels which leaves them of limited use for RNA-Seq data.

## Installation

To install this package, start R and enter:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dupRadar")
```

## Documentation

To view documentation for the version of this package installed in your system, start R and enter:

```r
browseVignettes("dupRadar")
```

or access the [pkgdown documentation](https://ssayols.github.io/dupRadar/index.html).

## Dependencies

### Imports (mandatory for core functionality)

* [Rsubread](http://bioconductor.org/packages/Rsubread/): Mapping, quantification and variant analysis of sequencing data.

## Citation

Sayols S, Scherzinger D, Klein H (2016). "dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data." BMC Bioinformatics, 17, 428. [doi: 10.1186/s12859-016-1276-2](http://dx.doi.org/10.1186/s12859-016-1276-2).
