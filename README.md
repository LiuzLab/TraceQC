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


## Examples

A FASTQ file and a reference file are required to use TraceQC. The reference is a text file which contains information as follows:

```
ATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGTAGACGCACCTCCACCCCACAGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAA
target 24 140
```

The first line of the reference file represents a construct sequence. The second line indicates the target region of the construct. In the lines, two numbers next to a region name specify the region's start and end locations. Locations should be 1-based, i.e., the first location is 1. If users want to add additional regions like the spacer region or PAM region, they can add more lines containing the additional regions. The format of the regions is the same as the target region. Here is an example of the reference file with other regions:


```
ATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGTAGACGCACCTCCACCCCACAGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAA
target 24 140
spacer 88 107
PAM 108 110
```


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

The following example shows how to create a TraceQC object.

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

