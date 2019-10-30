Variation Analysis
================
Marc Vaudel
2019-10-28

# 1\. Variation Analysis

In the [previous chapter](novel_peptides.md), novel peptides were
identified by Johansson *et al.* [(1)](#references) using a
proteogenomic search strategy. In addition, this search allowed
identifying peptides presenting single amino acid variants (SAAV).

## Libraries

We will need the following libraries, please make sure that they are
installed.

``` r
library(tidyr)
library(dplyr)
library(ggplot2)
library(scico)

theme_set(theme_bw(base_size = 11))
```

## Genotype - peptide abundance relationship

Figure 7h displays the intensity of peptides produced by the reference
and alternative alleles of corresponding SNPs for each genotype.

![Figure\_7h](resources/images/Fig7h.png?raw=true
"Johansson et al. Fig 7h")

> Figure 7h in Johansson *et al.*
[(1)](#references).

##### [:thought\_balloon:](answers.md#thought_balloon-if-we-assume-a-linear-relationship-between-number-of-alleles-and-peptide-abundance-what-should-be-the-peptide-distribution-for-each-genotype) *If we assume a linear relationship between number of alleles and peptide abundance, what should be the peptide distribution for each genotype?*

Note that the genotyping data corresponding to the peptide intensities
were not made available, but more examples can be found in
[Supplementary Figure
15](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf).

##### :speech\_balloon: Do the results presented follow a linear trend? What can affect the linearity of the relationship between number of alleles and intensity distribution?

## Abundance of variant peptides relative to the proteins

In [Supplementary Table
1](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape)
and [Supplementary Table
6](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape),
the authors provide the results of their proteogenomic analysis for
genes and for peptides carrying alternative alleles of SAAVs,
respectively. In addition, [Supplementary Table
1](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape)
contains experimental design and sample characterization information.
The tables were extracted to an R-friendly text format for this
tutorial, and are available in
[resources/data/tumor.gz](resources/data/tumor.gz),
[resources/data/proteins.gz](resources/data/proteins.gz), and
[resources/data/saav.gz](resources/data/saav.gz).

##### :pencil2: Load the data in R as in the code below.

``` r
tumorDF <- read.table(
    file = "resources/data/tumor.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)

genesDF <- read.table(
    file = "resources/data/genes.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
) %>%
    select(
        -contains("nspsm")
    )

saavDF <- read.table(
    file = "resources/data/saav.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)
```

Note that for the sake of speed, the last columns of the genes table
were
skipped.

##### [:thought\_balloon:](answers.md#thought_balloon-if-we-assume-a-linear-relationship-between-number-of-alleles-and-peptide-abundance-what-should-be-the-peptide-distribution-for-each-genotype) *How do we need to transform the tables to compare SAAV peptide-level intensities to gene-level intensities?*

For the sake of time, this data transformation was run for you already.
The proposed solution is available in
[scripts/merge\_table\_1\_6.R](resources/data/saav.gz), and the
reformated SAAV peptides table is available in
[resources/data/table16.gz](resources/data/table16.gz).

##### :pencil2: Load the reformated SAAV peptides table.

``` r
saavDF <- read.table(
    file = "resources/data/table16.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)
```

## References

1)  [Breast cancer quantitative proteome and proteogenomic
    landscape](https://www.ncbi.nlm.nih.gov/pubmed/30962452)
