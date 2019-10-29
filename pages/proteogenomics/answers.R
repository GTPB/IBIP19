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
