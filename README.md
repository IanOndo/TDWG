
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TDWG: cleaning occurrence records using the World Geographical Scheme for Recording Plant Distributions (WGSRPD) <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

`TDWG` (Ondo 2024) contains multiple functions to select or spatially
filter occurrence records for plant species. It combines botanical
information from hundreds of thousands of plant species from the Kew
World Checklist of Vascular Plants (WCVP) (Govaerts et al. 2021)
database with geographical information derived from The International
Working Group on Taxonomic Databases for Plant Sciences (TDWG) at
approximately “country” level and upwards.  The package was used to help
identifying plants geographic ranges at broad scale and cleaning
occurrence records for 36,687 distribution models in [*The global
distribution of plants used by
humans*](https://www.science.org/doi/10.1126/science.adg8028) (Pironon
and Ondo et al. 2024). 

## *Installation*

Make sure to have [*R*](https://cloud.r-project.org/ "R") or
[*Rstudio*](https://rstudio.com/products/rstudio/download/ "Rstudio")
installed on your machine. Some R packages need to be compiled from
source, so if you are on Windows, you need to install
[*Rtools*](http://cran.r-project.org/bin/windows/Rtools/) too.  

Install *TDWG* with the following instructions:

``` r
devtools::install_github("IanOndo/TDWG")
library(TDWG)
```

## Important note

The package is currently using the 2022 version of the WCVP database, so
no the lastest updates of the taxonomy. Current developments aim at
incorporating functions from the R package
[*rWCVP*](https://github.com/matildabrown/rWCVP) (Brown et al. 2023).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-rwcvp" class="csl-entry">

Brown, Matilda J. M., Barnaby E. Walker, Nicholas Black, Rafaël
Govaerts, Ian Ondo, Robert Turner, and Eimear Nic Lughadha. 2023.
“rWCVP: A Companion r Package to the World Checklist of Vascular
Plants.” *New Phytologist*.

</div>

<div id="ref-WCVP" class="csl-entry">

Govaerts, R., E. Nic Lughadha, N. Black, R. Turner, and A. Paton. 2021.
“The World Checklist of Vascular Plants, a Continuously Updated Resource
for Exploring Global Plant Diversity.” Journal Article. *Sci Data* 8
(1): 215. <https://doi.org/10.1038/s41597-021-00997-6>.

</div>

<div id="ref-TDWG" class="csl-entry">

Ondo, Ian. 2024. *TDWG: Cleaning Occurrence Records Using the World
Geographical Scheme for Recording Plant Distributions (WGSRPD)*.
<https://github.com/IanOndo/TDWG>.

</div>

<div id="ref-UsefulPlants" class="csl-entry">

Pironon and Ondo, M. Diazgranados, R. Allkin, A. C. Baquero, R.
Cámara-Leret, C. Canteiro, Z. Dennehy-Carr, et al. 2024. “The Global
Distribution of Plants Used by Humans.” *Science* 383 (6680): 293–97.
<https://doi.org/10.1126/science.adg8028>.

</div>

</div>
