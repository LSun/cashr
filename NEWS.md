# cashr 1.0.0

* New default solver (`solver = "newton"`): log-barrier Newton for both convex blocks with exact
  cutting planes for f >= 0. No MOSEK licence needed; fits converge fully and are 10-30x faster.
* The 2018 algorithm remains available as `solver = "mosek"` (SQUAREM + REBayes + MOSEK), ported
  from Rmosek's removed `scopt` interface to exponential cones (Rmosek >= 10).
* `cash()` now returns `pi0`, `lfsr`, `svalue`, `loglik` and timing; `cdf_g()`, `density_ecn()`,
  `ecn_moments()`, `omega_penalty()` and `hermite_orthonormal()` are exported.
* Fixed: `omega.pen`; `lfsr` for `gd.order < 3`; numerically stable Hermite basis.
* Tests (testthat), including a frozen reference fit of cashr 0.2.0 shared with the Python package.
* Conservative start and update order. `cash()` starts from "everything is null" (pi = point
  mass at 0, N(0, 1) noise) and by default updates the noise model omega before pi in each
  iteration (`gd.priority = TRUE`). Updating pi first, the 0.2.0 order (`gd.priority = FALSE`),
  can stop at a local optimum that calls all-null correlated data signal. On the paper's
  simulations both orders reach the same optimum.
* The block-coordinate ascent is monotone: a block update that would lower the objective is
  rejected.
* Numerical robustness of the solvers, found by auditing all 7,000 paper simulations against
  the Python package:
  - The cutting-plane working set is rebuilt on every omega-step. Carrying it across iterations
    accumulated near-duplicate tail constraints, stalled the barrier and could end a fit early.
  - Rows are rescaled to max 1, which leaves the optimum unchanged: observations with null
    density ~1e-170 at the all-null start overflowed the barrier Hessian.
  - The barrier line search accepts a point only if it is strictly feasible with the residuals
    recomputed exactly. Before, a slack of ~1e-16 could round to 0 and produce NaN.

# cashr 0.2.0

* Version used for the 2018 paper (arXiv:1812.07488v1/v2).
