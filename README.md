# TraceQC <img src="man/figures/hexsticker.png" align="right" height="140"/>

TraceQC is a R package for quality control (QC) of CRISPR Lineage Tracing Seqence Data. With a simple programming, users can create a HTML QC report page and plots for the QC. Here is an example QC report page: [Click here](https://htmlpreview.github.io/?https://github.com/LiuzLab/TraceQC-Supplementary/blob/master/docs/index.html)!

## Installation

Dependencies:
- reticulate which is an R interface to Python: https://rstudio.github.io/reticulate/
- ComplexHeatmap:http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
- DECIPHER: https://bioconductor.org/packages/release/bioc/html/DECIPHER.html
- 
Currently, installation of TraceQC is only available using `devtools`, and a personal Github token is required to install TraceQC. Please follow the following scripts for the installation inside a R session.

```
BiocManager::install(c("ComplexHeatmap", "DECIPHER"))
install.packages(c("fastqcr", "reticulate", "tictoc"))
```

```r
if(!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("LiuzLab/TraceQC")
```
After installation of TraceQC, run `install_external_packages()` to install required external tools/packages needed by TraceQC.


## Examples

A FASTQ file and a reference file are required to use TraceQC. The reference is a text file which contains information as follows:

```
ATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGTAGACGCACCTCCACCCCACAGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAA
target 23 140
spacer 87 107
PAM 107 110
```

The first line of the reference file represents a construct sequence. The other lines indicates target, spacer, and PAM regions of the construct. In these lines, two numbers next to a region name specify the start and end locations of the region. Be aware that locations are 0-based.

`generate_qc_report` is used to create a QC HTML report. The following script shows an example to create.

```r
library(TraceQC)
obj <- generate_qc_report(
 input_file = system.file("extdata", "test_data",
                          "fastq", "example_14d.fastq.gz", package="TraceQC"),
 ref_file = system.file("extdata", "test_data",
                        "ref", "ref.txt", package="TraceQC"),
 preview = T,
 title = "TraceQC report for example_14d.fastq",
 ncores=1
 )
summary(obj)
```

The following example shows how to create an TraceQC object.

```r
library(TraceQC)
library(fastqcr)
input_file <- system.file("extdata", "test_data",
                          "fastq", "example.fastq.gz", package="TraceQC")
ref_file <- system.file("extdata", "test_data", "ref",
                        "ref.txt", package="TraceQC")
qc_dir <- tempdir()
fastqc(system.file("extdata", "test_data",
                   "fastq", package = "TraceQC"),
       qc.dir=qc_dir)

input_qc_path <- get_qcpath(input_file, qc_dir)

obj <- TraceQC(input_file = input_file, ref_file = ref_file, fastqc_file = input_qc_path)
```

