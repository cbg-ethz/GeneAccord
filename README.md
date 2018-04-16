# GeneAccord #

This R package has functions to detect clonal exclusivity or co-occurrence patterns of altered genes or pathways in a cohort of cancer patients. Such clonally exclusive gene pairs may uncover synergies between co-existing clones.

### Dependencies ###

GeneAccord has been developed with R version 3.4.3. 

### Functionality ###

<p align="center">
	<img src="inst/ext/Cartoon_GeneAccordAlgorithm_Rpackage.png?raw=true" alt="Schematic overview of the procedure of GeneAccord"/>
</p>

The input data are the mutated gene-clone assignments from a cohort of cancer patients and from a collection of tree inferences to take into account uncertainty in the phylogenetic tree inference. The first step is to compute the overall rates of clonal exclusivity for all gene pairs for each patient separately. These rates reflect the expected prevalence of the clonal exclusivity pattern in each patient and for each gene or pathway pair. From these, the distribution of the test statistic under the null hypothesis of the likelihood ratio test is computed. For each gene pair, a paramater is computed that indicates whether the pair tends to be mutated in different clones more often than expected. Such pairs can be selected and tested for significance. Significant pairs may indicate synergistic effects between the co-existing clones.


Please refer to the vignette and to the R documentation of methods and classes.

### Install ###

Open R and input:

```{r}
install.packages("devtools")

# to install its dependencies
install.packages(c("biomaRt", "caTools", "dplyr", "ggplot2", "graphics", "grDevices", "gtools", "ggpubr", "magrittr", "maxLik", "RColorBrewer", "reshape2", "stats", "tibble", "utils"))

library(devtools)

install_github("cbg-ethz/GeneAccord")

library(GeneAccord)
```

### Contacts ###

Ariane L. Moore ( ariane.moore __ bsse.ethz.ch )

### Citation ###

Ariane L. Moore, Jack Kuipers, Jochen Singer, Elodie Burcklen, Peter Schraml, Christian Beisel, Holger Moch, Niko Beerenwinkel
Intra-tumor heterogeneity and co-occurring clones in renal cell carcinoma (2018)
bioRxiv number; doi: doi


