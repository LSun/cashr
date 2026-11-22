L <- 6

noise_density <- function(z, omega) as.vector(std.gaussDeriv(z, length(omega) - 1) %*% omega)

test_that("orthonormal Hermite recursion matches the explicit polynomials", {
  y <- seq(-4, 4, length = 17)
  H <- hermite_orthonormal(y, 4)
  He <- cbind(1, y, y^2 - 1, y^3 - 3 * y, y^4 - 6 * y^2 + 3)
  expect_equal(H, sweep(He, 2, sqrt(factorial(0:4)), "/"), ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("basis columns are scaled derivatives of phi", {
  z <- seq(-3, 3, length = 13); h <- 1e-5
  B <- std.gaussDeriv(z, L); Bp <- std.gaussDeriv(z + h, L); Bm <- std.gaussDeriv(z - h, L)
  for (l in 1:L) expect_equal(B[, l + 1], (Bp[, l] - Bm[, l]) / (2 * h) / sqrt(l), tolerance = 1e-6)
})

small <- function() {
  set.seed(0)
  list(x = c(-2.5, -0.3, 0, 0.8, 3.1), s = c(0.7, 1, 1.3, 0.9, 1.1), sd = c(0, 0.5, 1.5),
       omega = c(1, 0.1 * rnorm(L)))
}
quad <- function(f, lo = -40, hi = 40) stats::integrate(Vectorize(f), lo, hi, rel.tol = 1e-10, subdivisions = 1000L)$value

test_that("likelihood, posterior-mean and P(theta > 0) arrays match numerical integration", {
  p <- small()
  a <- cashr:::cash_arrays(p$x, p$s, p$sd, L)
  sc <- exp(a$log_scale)
  PM <- cashr:::posterior_mean_array(a); PP <- cashr:::positive_prob_array(a)
  for (j in seq_along(p$x)) for (k in 2:3) {
    dens <- function(th) stats::dnorm(th, 0, p$sd[k]) * noise_density((p$x[j] - th) / p$s[j], p$omega) / p$s[j]
    expect_equal(sum(a$F[j, k, ] * p$omega) * sc[j], quad(dens), tolerance = 1e-6)
    expect_equal(sum(PM[j, k, ] * p$omega) * sc[j], quad(function(th) th * dens(th)), tolerance = 1e-6)
    expect_equal(sum(PP[j, k, ] * p$omega) * sc[j], quad(dens, 0, 40), tolerance = 1e-6)
  }
})

sim <- function(n, seed) {
  set.seed(seed)
  theta <- ifelse(runif(n) < 0.7, 0, rnorm(n, 0, 2))
  s <- runif(n, 0.6, 1.4)
  list(x = theta + s * rnorm(n), s = s, theta = theta)
}

test_that("cash is calibrated on independent noise", {
  d <- sim(3000, 3)
  fit <- cash(d$x, d$s)
  expect_true(fit$converged)
  expect_lt(abs(fit$pi0 - 0.7), 0.06)
  disc <- fit$qvalue < 0.1
  expect_lt(mean(d$theta[disc] == 0), 0.15)
  expect_true(all(fit$lfdr >= 0 & fit$lfdr <= 1))
  expect_true(all(fit$lfsr >= fit$lfdr - 1e-12))
  expect_true(min(1 + cashr:::gauss_deriv_ratio(seq(-10, 10, by = 0.001), 10)[, -1] %*% fit$omega[-1]) > 0)
})

test_that("ecn recovers N(0,1) from iid noise", {
  set.seed(5)
  z <- rnorm(3000)
  fit <- ecn(z)
  expect_lt(max(abs(fit$omega[-1])), 0.1)
  expect_equal(ecn_moments(z, 4)[1], 1)
})

test_that("newton solver matches the original MOSEK algorithm", {
  skip_if(Sys.getenv("CASHR_TEST_MOSEK") != "1", "set CASHR_TEST_MOSEK=1 (needs Rmosek, REBayes, SQUAREM and a licence)")
  d <- sim(1000, 1)
  a <- cash(d$x, d$s)
  b <- cash(d$x, d$s, solver = "mosek", maxiter = 100)
  expect_gte(a$penloglik, b$penloglik - 1e-3)
  expect_equal(a$pi0, b$pi0, tolerance = 1e-2)
})

## Frozen fit of the original cashr 0.2.0 (MOSEK algorithm) on real correlated noise
## (GTEx liver null z-scores, n = 10,000); the Python package runs the same test.
test_that("cash and ecn match the original R + MOSEK reference fit", {
  ref <- function(f) read.csv(test_path("r_reference", f))
  d <- ref("ref_fit.csv"); g <- ref("ref_g.csv"); om <- ref("ref_omega.csv")
  obj <- ref("ref_obj.csv")[[1]]; z <- ref("z_null.csv")
  fit <- cash(d$x, d$s)
  expect_true(fit$converged)
  expect_gte(fit$penloglik, obj - 1e-6)
  expect_lt(abs(fit$pi0 - g$pi[1]), 1e-3)
  expect_lt(max(abs(fit$omega - om$omega)), 5e-3)
  expect_lt(max(abs(fit$lfdr - d$lfdr)), 2e-3)
  expect_lt(max(abs(fit$lfsr - d$lfsr)), 2e-3)
  expect_lt(max(abs(fit$pm - d$pm)), 1e-2)
  expect_equal(sum(fit$qvalue < 0.1), sum(d$q < 0.1))
  e <- ecn(z[[2]])
  expect_lt(max(abs(e$omega - om$ecn)), 1e-3)
})

## Fig. 5 simulation (g1 = 4, sim 234), where the omega constraint is nearly active at z ~ 9.5.
## A development version before 1.0.0 kept the cutting-plane working set across iterations;
## near-duplicate tail rows stalled the barrier and the fit stopped after a decrease
## (penloglik -18675.2962).
test_that("no stall when the noise constraint is active in the tail", {
  d <- read.csv(test_path("stall_case", "stall_case.csv"))
  fit <- cash(d$x, d$s)
  expect_true(fit$converged)
  expect_gt(fit$penloglik, -18674.955)
})

## Fig. 5 simulation (g1 = 2, sim 206) has |z| up to 28. At the all-null start the first
## omega-step sees null densities ~1e-170 for those rows; without row scaling 1 / r^2 overflowed
## the barrier Hessian. Best known penloglik: -17004.540428.
test_that("far-tail outliers at the all-null start do not break the omega-step", {
  d <- read.csv(test_path("tail_outlier", "tail_outlier.csv"))
  fit <- expect_no_warning(cash(d$x, d$s))
  expect_true(fit$converged)
  expect_gt(fit$penloglik, -17004.5406)
})

## Fig. 5 simulation (g1 = 4, sim 56). Late in a barrier centering a grid slack of ~1e-16 stayed
## positive when updated along the Newton direction but rounded to 0 when recomputed, giving
## Inf / NaN in the Hessian (qr.default: NA/NaN/Inf). Python fits it to -18070.248718.
test_that("barrier accepts only points that are feasible when recomputed", {
  d <- read.csv(test_path("slack_rounding", "slack_rounding.csv"))
  fit <- cash(d$x, d$s)
  expect_true(fit$converged)
  expect_equal(fit$penloglik, -18070.248718, tolerance = 1e-9)
})

## Real correlated null z-scores (paper Fig. 1, z.sel column 3; sd 1.70, no signal). Updating pi
## first (cashr 0.2.0's default order) stops at a local optimum with pi0 = 0.009 and every
## z-score called signal; the default omega-first order reaches pi0 = 1 (penloglik +424).
test_that("conservative omega-first start on inflated all-null data", {
  z <- read.csv(test_path("null_inflated", "null_inflated.csv"))$z
  fit <- cash(z, 1)
  expect_gt(fit$pi0, 0.99)
  expect_equal(sum(fit$qvalue <= 0.1), 0)
  old <- cash(z, 1, gd.priority = FALSE)
  expect_lt(old$pi0, 0.05)
  expect_gt(fit$penloglik, old$penloglik + 400)
})
