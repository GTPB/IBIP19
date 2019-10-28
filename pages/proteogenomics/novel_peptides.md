Novel Peptides
================
Marc Vaudel
2019-10-28

# 1\. Novel Peptides

Mutations in the genome can alter the sequence so that non-coding
sections become coding. *e.g.* through the introduction of a start codon
or the alteration of a stop codon. In order to identify these
non-canonical genomic products, protein databases that capture genetic
variation and non-canonical genomic products are generated either by
enriching canonical protein sequences or by running six reading frame
translation of the entire genome
[(1)](#references).

##### [:thought\_balloon:](answers.md#thought_balloon-based-on-you-knowledge-of-peptide-and-protein-identification-can-you-anticipate-challenges-posed-by-these-proteogenomic-databases) *Based on you knowledge of peptide and protein identification, can you anticipate challenges posed by these proteogenomic databases?*

## Processing

This tutorial is a notebook that contains [R](r-project.org) code that
can be run directly from the *Rmd* file. It assumes that the R working
directory is the proteogenomics folder of the repository, *e.g.*
`/myfolder/IBIP19/pages/proteogenomics`. We recommend using RStudio to
run this tutorial.

## Libraries

You will need the following libraries, please make sure that they are
installed.

``` r
library(dplyr)
```

## Data set

In this tutorial, we will analyze the non-canonical genomic products
identified in breast cancer by Johansson *et al.* [(2)](#references).
Note that this tutorial does not cover the database generation, search,
and validation of identification results. These bioinformatic procedures
are very demanding and we strongly advise to make sure that they are in
place at your lab or at the facility processing the data before
conducting any proteogenomic experiment. The proteogenomic
identification results by Johansson *et al.* [(2)](#references) are
reported in Supplementary Data 6, available
[here](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape)
in the course
repository.

##### [:thought\_balloon:](answers.md#thought_balloon-what-do-the-different-columns-in-the-table-represent) *What do the different columns in the table represent?*

For this tutorial, the *Novel Peptides* table was extracted to an
R-friendly text format, and is available
[here](resources/data/novel_peptides.txt).

##### :pencil2: Load the data in R as in the code below.

``` r
novelPeptidesDF <- read.table(
    file = "resources/data/novel_peptides.txt",
    header = T,
    sep = "\t",
    comment.char = "",
    stringsAsFactors = F
)
```

##### :pencil2: Find the different classes of loci represented.

``` r
knitr::kable(
    as.data.frame(
        table(
            novelPeptidesDF$class
            )
        )
    )
```

| Var1                     | Freq |
| :----------------------- | ---: |
| exonic                   |   17 |
| exonic-intronic          |    4 |
| exonic-ncRNA\_exonic     |    2 |
| exonic-splicing          |    1 |
| exonic-upstream          |    1 |
| exonic-UTR5              |   12 |
| intergenic               |  172 |
| intergenic-ncRNA\_exonic |    1 |
| intronic                 |   91 |
| intronic-exonic          |   10 |
| ncRNA\_exonic            |   18 |
| ncRNA\_exonic-exonic     |    1 |
| ncRNA\_intronic          |   22 |
| splicing-intronic        |    1 |
| upstream                 |    3 |
| UTR3                     |    1 |
| UTR5                     |   16 |
| UTR5-exonic              |   14 |
| UTR5-upstream            |    1 |

##### [:thought\_balloon:](answers.md#thought_balloon-based-on-you-knowledge-of-peptide-and-protein-identification-can-you-anticipate-challenges-posed-by-these-proteogenomic-databases) *How can we speculate on these different classes of loci yielding detectable peptides?*

## References

1)  [Proteogenomics: concepts, applications and computational
    strategies](https://www.ncbi.nlm.nih.gov/pubmed/25357241)

2)  [Breast cancer quantitative proteome and proteogenomic
    landscape](https://www.ncbi.nlm.nih.gov/pubmed/30962452)
