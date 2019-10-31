CNA-Protein
================
Marc Vaudel
2019-10-28

# CNA-Protein

In the previous chapters, we investigated the presence of novel or
variant peptides in Breast cancer tumors as done by Johansson *et al.*
[(1)](#references) using a proteogenomic search strategy. Here, we look
at the effect of copy number alteration (CNA). Copy number variation is
associated with increased or decreased gene expression, but the
association with protein abundances is less pronounced, one says that
the effect of CNAs is post-translationally attenuated
[(1,2)](#references).

##### :speech\_balloon: In you opinion, what mechanisms can affect the association between CNA, gene expression, and protein abundances?

## Libraries

We will need the following libraries, please make sure that they are
installed.

``` r
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtable)
library(grid)
library(conflicted)

theme_set(theme_bw(base_size = 11))
```

## CNA - RNA - Protein abundance relationship

Johansson *et al.* [(1)](#references) conduct a correlation analysis
similar to Gonçalves *et al.* [(3)](#references), where the correlation
of protein abundances with CNA is compared to the correlation of RNA
levels with CNA. These values are provided in [Supplementary Table
5](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape).
The correlation values were extracted in R-friendly format and are
availabe in
[resources/data/cna-rna-protein](resources/data/cna-rna-protein).

##### :pencil2: Import the correlation results.

``` r
# Load the file, select the columns of interest, format the columns and content

cnaCorDF <- read.table(
    file = "resources/data/cna-rna-protein",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)
```

##### [:thought\_balloon:](answers.md#thought_balloon-what-do-the-columns-represent-what-is-the-difference-between-pearson-and-spearman-correlations) What do the columns represent? What is the difference between Pearson and Spearman correlations?

##### :pencil2: Plot the correlation results for proteins against mRNA as in Figure 6 of Johansson *et al.* [(1)](#references) and Figure 1 of Gonçalves *et al.* [(3)](#references).

``` r
# Build the scatter plot

scatterPlot <- ggplot(
    data = cnaCorDF
) +
    geom_point(
        mapping = aes(
            x = mRNA_Spearman_correlation,
            y = protein_Spearman_correlation
        ),
        col = "black",
        alpha = 0.2
    ) +
    geom_density_2d(
        mapping = aes(
            x = mRNA_Spearman_correlation,
            y = protein_Spearman_correlation
        ),
        col = "white",
    ) +
    scale_x_continuous(
        name = "RNA-CNA Correlation (Spearman)"
    ) +
    scale_y_continuous(
        name = "Protein-CNA Correlation (Spearman)"
    )


# Build the density plots

rnaDensityPlot <- ggplot(
    data = cnaCorDF
) + theme_minimal() + 
    geom_density(
        mapping = aes(
            x = mRNA_Spearman_correlation
        ),
        fill = "black",
        alpha = 0.1
    ) +
    scale_x_continuous(
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        expand = c(0, 0)
    ) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()
    )

proteinDensityPlot <- ggplot(
    data = cnaCorDF
) + theme_minimal() + 
    geom_density(
        mapping = aes(
            x = protein_Spearman_correlation
        ),
        fill = "black",
        alpha = 0.1
    ) +
    scale_x_continuous(
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        expand = c(0, 0)
    ) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()
    ) + 
    coord_flip()


# Make grobs from plots

scatterGrob <- ggplotGrob(scatterPlot)
rnaDensityGrob <- ggplotGrob(rnaDensityPlot)
proteinDensityGrob <- ggplotGrob(proteinDensityPlot)


# Insert the densities as new row and column in the scatter grob

mergedGrob <- rbind(scatterGrob[1:4, ], rnaDensityGrob[7, ], scatterGrob[5:nrow(scatterGrob), ], size = "last")
mergedGrob$heights[5] <- unit(0.2, "null")

proteinDensityGrob <- gtable_add_rows(
        x = proteinDensityGrob, 
        heights = unit(1, "null"), 
        pos = 0
    )

mergedGrob <- cbind(mergedGrob[, 1:5], proteinDensityGrob[, 5], mergedGrob[, 6:ncol(mergedGrob)], size = "first")
mergedGrob$widths[6] <- unit(0.2, "null")


# Plot

grid.draw(mergedGrob)
```

![](cna-protein_files/figure-gfm/correlation_plot-1.png)<!-- -->

## References

1)  [Breast cancer quantitative proteome and proteogenomic
    landscape](https://www.ncbi.nlm.nih.gov/pubmed/30962452)

2)  [The genomic and transcriptomic architecture of 2,000 breast tumours
    reveals novel
    subgroups](https://www.ncbi.nlm.nih.gov/pubmed/22522925)

3)  [Widespread Post-transcriptional Attenuation of Genomic Copy-Number
    Variation in Cancer](https://www.ncbi.nlm.nih.gov/pubmed/29032074)
