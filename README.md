# cashr

`cashr`: Solve the Empirical Bayes Normal Means problem with correlated noise, with applications to large-scale multiple testing problems and the False Discovery Rate (FDR) control

## Installation

1. Install [`R`](https://cran.r-project.org).
2. Install the `Rmosek` package according to online instructions such as
- https://docs.mosek.com/8.1/rmosek/install-interface.html
- https://gist.github.com/mikelove/67ea44d5be5a053e599257fe357483dc
- https://rdrr.io/cran/ashr/f/inst/rmosek-mac.md
3. Once `Rmosek` is intalled, install `cashr` by running the following commands in the `R` interactive environment:
``` r
install.packages("devtools")
devtools::install_github("LSun/cashr")
```
