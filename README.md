# cashr

Empirical Bayes normal means with correlated noise, for large-scale multiple testing with
false discovery rate (FDR) control.

> Lei Sun and Matthew Stephens. *Solving the Empirical Bayes Normal Means Problem with Correlated
> Noise.* [arXiv:1812.07488](https://arxiv.org/abs/1812.07488).

You have estimates `x` with standard errors `s`, for example one per gene, and the noise
`(x - θ)/s` is N(0, 1) for each gene but correlated across genes. `cash()` estimates the
distribution of the effects θ together with the realized distribution of the correlated noise.
It returns local FDRs, q-values, posterior means and local false sign rates.

The same method in Python: [cashpy](https://github.com/LSun/cashpy).

## Install

```r
install.packages("remotes")
remotes::install_github("LSun/cashr")
```

No dependencies beyond base R, and no MOSEK licence.

## Quick start

The leukemia data of the paper (Table 1): 7,128 genes.

```r
library(cashr)
d <- read.csv("https://raw.githubusercontent.com/LSun/cashr_paper/master/examples/data/leukemia.csv")
fit <- cash(d$x, d$s)
fit
#> cash fit (newton): n = 7128, K = 25, L = 10, pi0 = 0.8147, penloglik = 2339.5795, 9 iterations, 3.06s
sum(fit$qvalue <= 0.1)            # discoveries at FDR 0.1
#> [1] 387
```

## What you get

| | |
|---|---|
| `fit$qvalue` | q-values: call `qvalue <= α` discoveries at FDR α |
| `fit$lfdr` | local false discovery rates |
| `fit$pm` | posterior means of θ |
| `fit$lfsr`, `fit$svalue` | local false sign rates and s-values |
| `fit$pi0` | estimated proportion of nulls |
| `fit$fitted_g` | prior on θ: mixture of N(0, sd²), with sd = 0 a point mass |
| `fit$omega` | coefficients of the fitted noise distribution |

For z-scores that are known to be all null, `ecn(z)` fits only the noise distribution, and
`plot()` shows it.

Main options of `cash()` (all documented in `?cash`):
- `s`: one standard error for all, or one per estimate. The default is 1, for z-scores.
- `gd.order = 10`: the number of terms in the noise model.
- `gd.priority = TRUE`: start from "everything is null" and update the noise model before the
  signal. This is the conservative choice.
- `solver = "mosek"` (optional): the 2018 algorithm, kept for comparison. It needs REBayes,
  SQUAREM, a [MOSEK licence](https://www.mosek.com/products/academic-licenses/) (free for
  academics) and Rmosek ≥ 10, which is not on CRAN: see
  [the Rmosek install guide](https://docs.mosek.com/12.0/rmosek/install-interface.html).

## More

- **Examples:** fits of three datasets from the paper, in
  [cashr_paper/examples](https://github.com/LSun/cashr_paper/tree/master/examples).
- **Reproducing the paper:** every figure, Table 1 and the numerical claims, in
  [cashr_paper/reproduce](https://github.com/LSun/cashr_paper/tree/master/reproduce).
- **Changes from 0.2.0:** see [NEWS.md](NEWS.md). The 2018 version is still available with
  `remotes::install_github("LSun/cashr@v0.2.0")`.

## License

GPL (≥ 3)
