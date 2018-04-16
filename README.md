# GeneAccord #

This R package has functions to detect clonal exclusivity or co-occurrence patterns of altered genes or pathways in a cohort of cancer patients.

### Dependencies ###

GeneAccord has been developed with R version 3.4.3. 

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

### Functionality ###

Please refer to the vignette and to the R documentation of methods and classes.


### Contacts ###

Ariane L. Moore ( ariane.moore __ bsse.ethz.ch )

### Citation ###

Ariane L. Moore, Jack Kuipers, Jochen Singer, Elodie Burcklen, Peter Schraml, Christian Beisel, Holger Moch, Niko Beerenwinkel
Intra-tumor heterogeneity and co-occurring clones in renal cell carcinoma (2018)
bioRxiv number; doi: doi


