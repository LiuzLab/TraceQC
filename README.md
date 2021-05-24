# TraceQC <img src="man/figures/hexsticker.png" align="right" height="140"/>

TraceQC is a R package for quality control (QC) of CRISPR Lineage Tracing Sequence Data. With simple programming, users can create an HTML QC report page and plots for the QC. Here is an example QC report page: [Click here](https://htmlpreview.github.io/?https://github.com/LiuzLab/TraceQC-Supplementary/blob/master/docs/index.html)!

## Installation

Dependencies:
- reticulate which is an R interface to Python: https://rstudio.github.io/reticulate/
- ComplexHeatmap:http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
- DECIPHER: https://bioconductor.org/packages/release/bioc/html/DECIPHER.html

Currently, installing TraceQC is only available using `devtools`, and a personal Github token is required to install TraceQC. Please follow the following scripts for the installation inside an R session.

```
BiocManager::install(c("ComplexHeatmap", "DECIPHER"))
install.packages(c("fastqcr", "reticulate", "tictoc"))
```

```r
if(!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("LiuzLab/TraceQC")
```
After installation of TraceQC, run `install_external_packages()` to install required external tools/packages needed by TraceQC.
