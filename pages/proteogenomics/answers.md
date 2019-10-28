---
layout: page
title: Answers-Proteogenomics
---

# Answers - Proteogenomics

##### [:speech_balloon:](novel_peptides.md#speech_balloon-based-on-you-knowledge-of-peptide-and-protein-identification-can-you-anticipate-challenges-posed-by-these-proteogenomic-databases) _Based on you knowledge of peptide and protein identification, can you anticipate challenges posed by these proteogenomic databases?_

The first challenge posed by proteogenomic databases is that they are very large, and therefore require much more computational power. The inflation in the collection of possible sequences, called search space, also increases the probability that two peptide sequences are similar, hence reducing the discrimination power of search engines. Thus reduced discrimination power results in a reduced identification rate at fixed error rate [(1)](#references). The estimation of error rates in proteomics usually relies on the target-decoy approach [(2)](#references), which is very reliable to track random matches, but not the errors due to partial matches [(3)](#references). Consequently, error rate estimation is very challenging in proteogenomic databases and requires spectrum-level inspection of the matches [(4)](#references). Finally, The increased number of similar sequences reduces the probability to find peptides unique to a protein, hence complexifying protein inference [(5)](#references).

##### [:speech_balloon:](novel_peptides.md#speech_balloon-what-do-the-different-columns-in-the-table-represent) _What do the different columns in the table represent?_

This document contains two tables, a list of novel peptides and a list of peptides containing single amino acid variats (SAAV).

- Novel Peptides

| Column Name | Description |
| ----------- | ----------- |
| Peptide Sequence | The amino acid sequence of the identified peptide. |
| MSGF+ SpecEva | The score (e-value) produced by the search engine (ms-gf+) when evaluating the match between the peptide and the fragmentation spectrum. The lower, the better. |
| # PSMs | The number of spectra where this peptide was matched. The more, the better. |
| Annotation - GRCh37 - hg19 | The annotation of the genetic locus coding the peptide. _GRCh37 - hg19_ corresponds to the build of the genome used for the analysis. |
| Chromosome, start, end, strand, locus number | Genomic coordinates of the locus coding the peptide in the GRCh37 build. |
| Category and sequence similarity to known proteins | Quality control report on possible mismatch with other proteins. Note that this is provided at two different steps of the bioinformatic pipeline. See _Proteogenomics search and Class-specific FDR_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |
| Peptide explained by nsSNP in CanProVar 2.0 | Known variation that could explain the peptide. |
| hg38_coordinates | Genomic coordinates of the locus coding the peptide in the GRCh38 build. |
| (closest) matched protein | Protein best matched with this peptide. |
| Nterm.seq.3aa. - Aligned sequence - Cterm.seq.3aa. | Sequence, upstream, and downstream amino acids in the mached protein. |
| Identity | Identity score between the peptide and its matched counterpart. The higher the better. |
| Peptide length - Alignment length - # mismatches - # gaps | Summary information on the identified and matched peptides. |
| top 0.2% MHCflurry | class I peptide/MHC binding affinity prediction: see _NSearch neoantigen candidates in draft proteomes data and MHC binding prediction_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |
| Identified in Kim et al draft proteome | Indicates whether the peptide was identified in [(6)](#references). |
| Class | Class of the locus coding this peptide. |
| Associated gene (closest genes in the genome) | Nearest gene for the given locus. |
| Category | Category of the nearest gene. |
| TMT quantification | Abundance of the peptide in the different tumors screened. |
| Orthogonal data | Additional data on the peptide: see _Orthogonal validation data of proteogenomics searches_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |

- SAAV

| Column Name | Description |
| ----------- | ----------- |
| Peptide Sequence | The amino acid sequence of the identified peptide. |
| Position of amino acid substitution | The position of the SAAV on the peptide sequence. |
| Matched protein(s) | The proteins matched. |
| Peptide sequence with modifications | The peptide sequence annotated with the masses of the modifications. |
| Peptide length | The sequence length in number of amino acids. |
| Fragmentation Method | The method used to fragment the peptide during mass spectrometry analysis. |
| Precursor m/z | The mass over charge ratio of the precursor selected for fragmentation. |
| IsotopeError | The isotopic difference between the measured m/z and the m/z of the peptide. |
| Precursor Error (ppm) | The relative difference in parts per million (ppm) between the measured m/z and the m/z of the peptide after isotopic correction. |
| Charge | The charge of the precursor. |
| DeNovoScore, MSGFScore | Score of the match between the fragmentation spectrum and the peptide. The higher, the better. |
| SpecEValue, EValue | E-value of the score of the match between the fragmentation spectrum and the peptide. The lower, the better. |
| SAAV Only - Peptide FDR | FDR estimation at the given score for the SAAV peptides only. See _Proteogenomics search and Class-specific FDR_ in the [Supplementary Information](../resources/Johansson_et_al_breast_cancer_quantitative_proteome_and_proteogenomic_landscape/supplementary_information.pdf). |
| SpectraFile, TMT set, ScanNum | File, TMT set, and scan number for the spectrum where the peptide was identified. |
| SpectrumAI verification | Summary information of the spectrum-level validation of the SAAV using flanking fragment ions [(4)](#references). |
| TMT intensity | Intensity of the fragment ions from the different TMT channels. |



## References

(1) [Anatomy and evolution of database search engines-a central component of mass spectrometry based proteomic workflows](https://www.ncbi.nlm.nih.gov/pubmed/28902424)

(2) [Target‐decoy search strategy for increased confidence in large‐scale protein identifications by mass spectrometyr](https://www.ncbi.nlm.nih.gov/pubmed/17327847)

(3) [Analysis of the resolution limitations of peptide identification algorithms](https://www.ncbi.nlm.nih.gov/pubmed/21995378)

(4) [Discovery of coding regions in the human genome by integrated proteogenomics analysis workflow](https://www.ncbi.nlm.nih.gov/pubmed/29500430)

(5) [Interpretation of shotgun proteomic data: the protein inference problem](https://www.ncbi.nlm.nih.gov/pubmed/16009968)

(6) [A draft map of the human proteome](https://www.ncbi.nlm.nih.gov/pubmed/24870542)
