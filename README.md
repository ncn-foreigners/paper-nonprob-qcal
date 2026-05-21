# Repository for the paper "Quantile balancing inverse probability weighting for non-probability samples"

## Paper

+ Paper: Beręsewicz, M., Szymkowiak, M. and Chlebicki, P. (2025). Quantile balancing inverse probability weighting for non-probability samples. Survey Methodology, 51(2), 533-559. Paper available at http://www.statcan.gc.ca/pub/12-001-x/2025002/article/00005-eng.pdf.
+ [Arxiv version](https://arxiv.org/abs/2403.09726)

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

## License

- This repository is licensed under Creative Commons Attribution 4.0 International (CC BY 4.0). See `LICENSE`.
- Scope and third-party-material note: see `LICENSE-content`.
