# Repository for the paper "Quantile balancing inverse probability weighting for non-probability samples"

## Paper

+ Accepted to Survey Methodology journal (planned for December 2025).
+ [Arxiv version](https://arxiv.org/abs/2403.09726)
+ Latest version [here](paper/paper-nonprob-qcal.pdf).

## Requirements

R packages

``` r
install.packages("nonprobsvy")
install.packages(c("jointCalib", "sampling", "laeken", "survey", "data.table", "ggplot2")) ## statistical
install.packages(c("data.table", "ggplot2", "xtable", "stringi")) ## processing
install.packages(c("doSNOW", "progress", "foreach")) ## paralell computing
```

Codes were developed under the following R version

``` r
> sessionInfo()
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.3
```

## Structure

Structure of the repo:

-   `codes/`:
    -   `sim-code-yang-2020.R` -- function with all estimators discussed in the paper
    -   `sim-code-run.R` -- script to run function `yang_sim()` and the main simulation
    -   `sim-processing-results.qmd` -- notebook with codes processing results
-   `figs` -- figures for the plot
-   `results/` -- simulation results in the `RDS` format
-   `paper/` -- pdf file with the paper and source files

## Financing

Work on this paper was supported by the National Science Centre, OPUS 20 grant no. 2020/39/B/HS4/00941.
