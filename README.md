# TraceQC

Trace is a package for quality control of CRISPR Lineage Tracing Seqence Data.

## Installation

Currently, installation of TraceQC is only available using `devtools`, and a personal Github token is required to install TraceQC. Please follow the following script for the installation.

```
if(!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
AUTH_TOKEN <- "YOUR TOKEN can be found at https://github.com/settings/tokens"
devtools::install_github("LiuzLab/TraceQC", auth_token=AUTH_TOKEN)
```

