##
#
# Code used to generate the figures in the answers.
# 
##

# Peptide intensity vs allele distribution

peptideAlleleDF <- data.frame(
    genotype = factor(rep(c("AA", "AB", "BB"), 2)),
    intensity = c(0, 50, 100, -100, -50, -0),
    peptide = factor(c(rep("Alt", 3), rep("Ref", 3)), levels = c("Ref", "Alt"))
)

peptideAllelePlot <- ggplot(
    data = peptideAlleleDF
) +
    theme_bw(
        base_size = 16
    ) +
    geom_col(
        mapping = aes(
            x = genotype,
            y = intensity,
            fill = peptide
        )
    ) +
    scale_x_discrete(
        name = "Genotype"
    ) + 
    scale_y_continuous(
        name = "Peptide Abundance [%]",
        breaks = (-2:2)*50,
        labels = abs((-2:2)*50)
    ) +
    geom_hline(
        yintercept = 0
    ) +
    scale_fill_manual(
        name = "SAAV",
        values = scico(
            n = 2,
            palette = "cork",
            begin = 0.2,
            end = 0.8,
            direction = -1
        ),
        breaks = c("Alt", "Ref")
    ) +
    theme(
        legend.position = "top",
        panel.grid.minor.y = element_blank()
    )

png("pages/proteogenomics/resources/images/peptideAllelePlot.png", width = 800, height = 600)
peptideAllelePlot
dummy <- dev.off()


# Logged peptide and protein intensity distributions

saavDF <- read.table(
    file = "resources/data/table16.gz",
    header = T,
    sep = "\t",
    comment.char = "",
    quote = "",
    stringsAsFactors = F
)

tumorColumns <- c("OSL.53E_saavPeptide", "OSL.53E_gene", "OSL.567_saavPeptide", "OSL.567_gene", "OSL.3FF_saavPeptide", "OSL.3FF_gene", "OSL.55F_saavPeptide", "OSL.55F_gene", "OSL.46A_saavPeptide", "OSL.46A_gene", "OSL.4B0_saavPeptide", "OSL.4B0_gene", "OSL.4D6_saavPeptide", "OSL.4D6_gene", "OSL.485_saavPeptide", "OSL.485_gene", "OSL.41B_saavPeptide", "OSL.41B_gene", "OSL.4AF_saavPeptide", "OSL.4AF_gene", "OSL.46E_saavPeptide", "OSL.46E_gene", "OSL.494_saavPeptide", "OSL.494_gene", "OSL.457_saavPeptide", "OSL.457_gene", "OSL.48B_saavPeptide", "OSL.48B_gene", "OSL.4B4_saavPeptide", "OSL.4B4_gene", "OSL.449_saavPeptide", "OSL.449_gene", "OSL.44E_saavPeptide", "OSL.44E_gene", "OSL.3EB_saavPeptide", "OSL.3EB_gene", "OSL.43C_saavPeptide", "OSL.43C_gene", "OSL.493_saavPeptide", "OSL.493_gene", "OSL.4D9_saavPeptide", "OSL.4D9_gene", "OSL.56F_saavPeptide", "OSL.56F_gene", "OSL.405_saavPeptide", "OSL.405_gene", "OSL.49E_saavPeptide", "OSL.49E_gene", "OSL.441_saavPeptide", "OSL.441_gene", "OSL.430_saavPeptide", "OSL.430_gene", "OSL.4FA_saavPeptide", "OSL.4FA_gene", "OSL.43A_saavPeptide", "OSL.43A_gene", "OSL.406_saavPeptide", "OSL.406_gene", "OSL.47C_saavPeptide", "OSL.47C_gene", "OSL.524_saavPeptide", "OSL.524_gene", "OSL.443_saavPeptide", "OSL.443_gene", "OSL.458_saavPeptide", "OSL.458_gene", "OSL.53D_saavPeptide", "OSL.53D_gene", "OSL.540_saavPeptide", "OSL.540_gene", "OSL.42E_saavPeptide", "OSL.42E_gene", "OSL.40A_saavPeptide", "OSL.40A_gene", "OSL.40E_saavPeptide", "OSL.40E_gene", "OSL.3FA_saavPeptide", "OSL.3FA_gene", "OSL.521_saavPeptide", "OSL.521_gene", "OSL.46D_saavPeptide", "OSL.46D_gene", "OSL.54D_saavPeptide", "OSL.54D_gene", "OSL.4BA_saavPeptide", "OSL.4BA_gene", "OSL.579_saavPeptide", "OSL.579_gene", "OSL.57B_saavPeptide", "OSL.57B_gene")

intensitiesDF <- saavDF %>%
    gather(
        !!tumorColumns,
        key = "sample",
        value = "intensity"
    ) %>%
    filter(
        !is.na(intensity) & !is.infinite(intensity)
    ) %>%
    separate(
        col = sample,
        into = c("tumor", "species"),
        sep = "_"
    )

intensitiesDF$species <- factor(
    intensitiesDF$species, 
    levels = c("gene", "saavPeptide")
)
levels(intensitiesDF$species) <- c("Gene", "SAAV Peptide")

plot <- ggplot(
    data = intensitiesDF
) +
    theme_bw(
        base_size = 16
    ) +
    geom_violin(
        mapping = aes(
            x = tumor,
            y = log10(intensity),
            col = tumor
        ),
        fill = NA
    ) +
    geom_boxplot(
        mapping = aes(
            x = tumor,
            y = log10(intensity),
            col = tumor
        ),
        fill = NA,
        outlier.shape = NA
    ) +
    facet_grid(
        species ~ .
    ) + 
    scale_fill_manual(
        values = scico(
            n = length(unique(intensitiesDF$tumor)),
            palette = "batlow",
            begin = 0.25,
            end = 0.75
        )
    ) +
    scale_color_manual(
        values = scico(
            n = length(unique(intensitiesDF$tumor)),
            palette = "batlow",
            begin = 0.25,
            end = 0.75
        )
    ) + 
    scale_y_continuous(
        name = "Intensity [log10]"
    ) + 
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
            
        ),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

png("pages/proteogenomics/resources/images/peptideProteinDistribution.png", width = 800, height = 1200)
plot
dummy <- dev.off()
