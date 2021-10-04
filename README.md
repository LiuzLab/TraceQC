# TraceQC <img src="man/figures/hexsticker.png" align="right" height="140"/>

TraceQC is a R package for quality control (QC) of CRISPR Lineage Tracing Sequence Data. 

## Installation

```r
if(!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("LiuzLab/TraceQC")
```
To install the Python packages for traceQC, run:
```
pip install biopython pandas tqdm
```

## Tutorial
The tutorial of TraceQC pipeline for bulk DNA sequencing is [here](https://liuzlab.github.io/TraceQC/tutorial). The dataset is sampled from [hgRNA dataset](https://www.nature.com/articles/nmeth.4108).

The tutorial of TraceQC pipeline for single-cell RNA sequencing is [here](https://liuzlab.github.io/TraceQC/tutorial_sc). The dataset is sampled from [Carlin dataset](https://www.sciencedirect.com/science/article/pii/S0092867420305547).

## TraceQC annotated reference format
The reference is a text file which contains information as follows:

    ATGGACTATCATATGCTTACCG...CCGGTAGACGCACCTCCACCCCACAGTGGGGTTAGAGCTAGAAATA
    target 23 140

The first line of the reference file is required should be the construct
sequence. The second line is also required should be the target barcode
region of the construct. In this lines, two numbers next to a region
name specify the start and end locations of the region. Locations should
be 1-based, i.e. the first location is indicated as 1. Users can
optionally add additional regions such as spacer region or PAM region in
the same format. Here is an example of the refenence file with
additional regions:

    ATGGACTATCATATGCTTACCG...CCGGTAGACGCACCTCCACCCCACAGTGGGGTTAGAGCTAGAAATA
    target 24 140
    spacer 88 107
    PAM 108 110

The examples of annotated hgRNA reference sequence is aviable [here](https://github.com/LiuzLab/TraceQC/blob/master/inst/extdata/test_data/ref/ref_hgRNA_invitro.txt). 
The examples of annotated Carlin reference sequence is aviable [here](https://github.com/LiuzLab/TraceQC/blob/master/inst/extdata/test_data/ref/ref_carlin.txt). 

## TraceQC output format
| Column      | Description |
| ----------- | ----------- |
| character      | A mutation identification string       |
| type   | The type of mutation (deletion, insertion and substitution).        |
| start      | The starting positioin of mutation.       |
| length   | The length of mutation.        |
| alt   | The altered sequence.        |
| count   | The read count of mutation.        |
| cell   | The cell IDs that contain this mutation.        |

## Documentation
The full documentation of TraceQC functions is available [here](https://github.com/LiuzLab/TraceQC/blob/master/TraceQC_1.0.1.pdf).

References
----------
Kalhor, R., Mali, P., & Church, G. M. (2017). Rapidly evolving homing CRISPR barcodes. Nature methods, 14(2), 195-200.

Bowling, S., Sritharan, D., Osorio, F. G., Nguyen, M., Cheung, P., Rodriguez-Fraticelli, A., ... & Camargo, F. D. (2020). An engineered CRISPR-Cas9 mouse line for simultaneous readout of lineage histories and gene expression profiles in single cells. Cell, 181(6), 1410-1422.

