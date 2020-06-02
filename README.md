# TraceQC <img src="man/figures/hexsticker.png" align="right" height="140"/>

TraceQC is a package for quality control of CRISPR Lineage Tracing Seqence Data.

## Installation

Dependencies:
- reticulate which is an R interface to Python: https://rstudio.github.io/reticulate/
- ComplexHeatmap:http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
- DECIPHER: https://bioconductor.org/packages/release/bioc/html/DECIPHER.html

Currently, installation of TraceQC is only available using `devtools`, and a personal Github token is required to install TraceQC. Please follow the following scripts for the installation inside a R session.

```
BiocManager::install(c("circlize", "ComplexHeatmap", "DECIPHER"))
install.packages(c("fastqcr", "reticulate", "tictoc"))
```

```r
if(!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
AUTH_TOKEN <- "YOUR TOKEN can be found at https://github.com/settings/tokens"
devtools::install_github("LiuzLab/TraceQC", auth_token=AUTH_TOKEN)
```
After installation of TraceQC, run `install_external_packages()` to install required external tools/packages needed by TraceQC.


## Example

```r
library(TraceQC)
library(fastqcr)
input_file <- system.file("extdata", "test_data",
                          "fastq", "example.fastq", package="TraceQC")
ref_file <- system.file("extdata", "test_data", "ref",
                        "ref.txt", package="TraceQC")
qc_dir <- tempdir()
fastqc(system.file("extdata", "test_data",
                   "fastq", package = "TraceQC"),
       qc.dir=qc_dir)

input_qc_path <- get_qcpath(input_file, qc_dir)

example_obj <- TraceQC(input_file = input_file,
                       ref_file = ref_file,
                       fastqc_file = input_qc_path)

save(example_obj, file="data/example_obj.rda")
```
