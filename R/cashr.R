#' @title Fit Empirical Bayes Normal Means with Correlated Noise
#'
#' @description This is the main interface for fitting correlated EBNM models
#'   based on algorithms proposed by Sun and Stephens.  The default
#'   behaviour is simply to run the biconvex optimization and return
#'   the result and posterior calculations.
#'
#' @param x A p vector of observations
#' @param s A scalar or a p vector of standard deviations.
#' @param deltaAt0 Logical, indicating whether to use a point mass at zero as one of components for a mixture distribution of the prior.
#' @param gd.order
#' @param omega.lambda
#' @param omega.rho
#' @param omega.pen
#' @param mixsd.mult
#' @param gd.priority Logical, indicating whether to optimizer over prior
#' @param control
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom utils capture.output modifyList
#'
cash = function (x, s = 1,
                 deltaAt0 = TRUE,
                 gd.order = 10,
                 omega.lambda = 10, omega.rho = 0.5, omega.pen = NULL,
                 mixsd.mult = sqrt(2),
                 gd.priority = FALSE,
                 control = list(maxiter = 50)) {
  if (s > 0 & length(s) == 1L) {
    L <- gd.order

    s = rep(s, length(x))
  } else if (length(x) != length(s) | !all(s > 0)) {
    stop("s should either be a positive number or a vector of positive numbers with the same length of x")
  }

  ## setting a dense grid of sd for gaussian mixture prior
  if (deltaAt0) {
    sd = c(0, autoselect.mixsd(x, s, mult = mixsd.mult))
    pi_prior = c(10, rep(1, length(sd) - 1))
  } else {
    sd = autoselect.mixsd(x, s, mult = mixsd.mult)
    pi_prior = rep(1, length(sd))
  }

  array_F = array_f(x, s, sd, L, mixcompdist = 'normal', gd.normalized = TRUE)
  array_F = aperm(array_F, c(2, 3, 1))

  if (is.null(omega.pen)) {
    if (is.null(omega.lambda)) {
      omega_prior = rep(0, L)
    } else if (is.null(omega.rho)) {
      omega_prior = rep(omega.lambda, L)
      omega_prior[seq(1, L, by = 2)] = 0
    } else {
      omega_prior = omega.lambda / sqrt(omega.rho^(1 : L))
      omega_prior[seq(1, L, by = 2)] = 0
    }
  } else if (length(omega_pen) == length(L)) {
    omega_prior = omega.pen
  } else {
    stop('The length of penalty parameters of omega should be L.')
  }

  hermite = Hermite(L)
  z.extra = seq(-10, 10, by = 0.001)
  gd0.std = dnorm(z.extra)
  matrix_lik_z = cbind(gd0.std)
  for (i in 1 : L) {
    gd.std = (-1)^i * hermite[[i]](z.extra) * gd0.std / sqrt(factorial(i))
    matrix_lik_z = cbind(matrix_lik_z, gd.std)
  }

  res = biopt(array_F, matrix_lik_z, pi_prior, omega_prior, control, primal = FALSE, gd.priority)
  pihat = res$pihat
  what = res$what
  fitted_g = normalmix(pi = pihat, mean = 0, sd = sd)
  penloglik = -res$B

  array_PM = array_pm(x, s, sd, L, gd.normalized = TRUE)
  array_PM = aperm(array_PM, c(2, 3, 1))
  theta_pm = colSums(t(apply(pihat * array_PM, 2, colSums)) * what) / colSums(t(apply(pihat * array_F, 2, colSums)) * what)

  lfdr = lfdr_top(pihat[1], what, x, s, L) / colSums(t(apply(pihat * array_F, 2, colSums)) * what)
  qvalue = qval.from.lfdr(lfdr)

  output <- list(fitted_g = fitted_g,
                 gd.order = L,
                 omega = what,
                 penloglik = penloglik,
                 niter = res$niter,
                 converged = res$converged,
                 pm = theta_pm,
                 lfdr = as.vector(lfdr),
                 qvalue = qvalue,
                 x = x,
                 s = s,
                 deltaAt0 = deltaAt0,
                 array_F = array_F,
                 array_PM = array_PM
                 )
  class(output) <- "cash"

  return(output)
}

#' @export

summary.cash <- function (output, ...) {
  output[1 : 6]
}

#' @export

print.cash <- function (output, ...) {
  print(summary.cash(output, ...))
}

#' @export

get_svalue <- function (output) {
  array_PP <- array_PosProb(output$x, output$s, output$deltaAt0, output$fitted_g$sd, gd.ord = output$gd.order, gd.normalized = TRUE)
  array_PP = aperm(array_PP, c(2, 3, 1))
  theta_PosProb <- colSums(t(apply(output$fitted_g$pi * array_PP, 2, colSums)) * output$omega) / colSums(t(apply(output$fitted_g$pi * output$array_F, 2, colSums)) * output$omega)
  lfsr <- compute_lfsr(1 - output$lfdr - theta_PosProb, output$lfdr)
  svalue <- qval.from.lfdr(lfsr)
  return(list(lfsr = lfsr,
              svalue = svalue))
}

compute_lfsr <- function (NegativeProb, ZeroProb) {
  ifelse(NegativeProb > 0.5 * (1 - ZeroProb), 1 - NegativeProb, NegativeProb + ZeroProb)
}

bifixpoint = function(pinw, array_F, matrix_lik_z, pi_prior, w_prior, primal, gd.priority){
  Kpi = dim(array_F)[1]
  Lw = dim(array_F)[2]
  pi = pinw[1 : Kpi]
  w = pinw[-(1 : Kpi)]
  if (gd.priority) {
    matrix_lik_w = apply(array_F * pi, 2, colSums)
    if (primal) {
      g_current = matrix_lik_w %*% w
      w_new = c(1, w.mosek.primal(matrix_lik_w, w_prior, w.init = c(g_current, w[-1]))$w)
    } else {
      w_new = c(1, w.mosek(matrix_lik_w, matrix_lik_z, w_prior, w.init = w)$w)
    }
    matrix_lik_pi = apply(aperm(array_F, c(2, 1, 3)) * w_new, 2, colSums)
    pi_new = mixIP(matrix_lik_pi, pi_prior, pi)$pihat
  } else {
    matrix_lik_pi = apply(aperm(array_F, c(2, 1, 3)) * w, 2, colSums)
    pi_new = mixIP(matrix_lik_pi, pi_prior, pi)$pihat
    matrix_lik_w = apply(array_F * pi_new, 2, colSums)
    if (primal) {
      g_current = matrix_lik_w %*% w
      w_new = c(1, w.mosek.primal(matrix_lik_w, w_prior, w.init = c(g_current, w[-1]))$w)
    } else {
      w_new = c(1, w.mosek(matrix_lik_w, matrix_lik_z, w_prior, w.init = w)$w)
    }
  }
  # w_new = c(1, w.cvxr.uncns(matrix_lik_w, w.init = w)$primal_values[[1]])
  return(c(pi_new, w_new))
}

# w.cvxr.uncns = function (matrix_lik_w, w.init = NULL) {
#   FF <- matrix_lik_w[, -1]
#   f <- matrix_lik_w[, 1]
#   p <- ncol(FF)
#   w <- CVXR::Variable(p)
#   objective <- CVXR::Maximize(CVXR::SumEntries(CVXR::Log(FF %*% w + f)))
#   prob <- CVXR::Problem(objective)
#   if (is.null(w.init)) {
#     capture.output(result <- solve(prob), file = "/dev/null")
#   } else {
#     capture.output(result <- solve(prob, warm_start = w.init[-1]), file = "/dev/null")
#   }
#   return(result)
# }

w.mosek = function (matrix_lik_w, matrix_lik_z, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[, 1]
  B = matrix_lik_z[, -1]
  b = matrix_lik_z[, 1]
  m = ncol(A)
  nA = nrow(A)
  nB = nrow(B)
  AB <- rbind(A, B)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    g.init <- as.vector(matrix_lik_w %*% w.init)
    v.init <- c(1 / g.init, rep(0, nB))
    v.init.list <- list(xx = v.init)
    P$sol <- list(itr = v.init.list, bas = v.init.list)
  }
  P$c <- c(a, b)
  P$A <- Matrix::Matrix(t(AB), sparse = TRUE)
  if (is.null(w_prior) | all(w_prior == 0) | missing(w_prior)) {
    P$bc <- rbind(rep(0, m), rep(0, m))
  } else {
    P$bc <- rbind(-w_prior, w_prior)
  }
  P$bx <- rbind(rep(0, nA + nB), rep(Inf, nA + nB))
  opro <- matrix(list(), nrow = 5, ncol = nA)
  rownames(opro) <- c("type", "j", "f", "g", "h")
  opro[1, ] <- as.list(rep("log", nA))
  opro[2, ] <- as.list(1 : nA)
  opro[3, ] <- as.list(rep(-1, nA))
  opro[4, ] <- as.list(rep(1, nA))
  opro[5, ] <- as.list(rep(0, nA))
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$suc - z$sol$itr$slc
  list(w = w, status = status)
}

w.mosek.primal = function (matrix_lik_w, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[,1]
  m = ncol(A)
  n = nrow(A)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    w.init.list <- list(xx = w.init)
    P$sol <- list(bas = w.init.list)
  }
  P$c <- rep(0, n + m)
  P$A <- Matrix::Matrix(cbind(diag(n), -A), sparse = TRUE)
  P$bc <- rbind(a, a)
  P$bx <- rbind(c(rep(0, n), rep(-Inf, m)),
                c(rep(Inf, n), rep(Inf, m)))
  if (missing(w_prior) | is.null(w_prior) | all(w_prior == 0)) {
    opro <- matrix(list(), nrow = 5, ncol = n)
    rownames(opro) <- c("type", "j", "f", "g", "h")
    opro[1, ] <- as.list(rep("log", n))
    opro[2, ] <- as.list(1 : n)
    opro[3, ] <- as.list(rep(-1, n))
    opro[4, ] <- as.list(rep(1, n))
    opro[5, ] <- as.list(rep(0, n))
  } else {
    opro <- matrix(list(), nrow = 5, ncol = n + m)
    rownames(opro) <- c("type", "j", "f", "g", "h")
    opro[1, ] <- as.list(c(rep("log", n), rep("pow", m)))
    opro[2, ] <- as.list(1 : (n + m))
    opro[3, ] <- as.list(c(rep(-1, n), w_prior))
    opro[4, ] <- as.list(c(rep(1, n), rep(2, m)))
    opro[5, ] <- as.list(rep(0, n + m))
  }
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$xx[-(1 : n)]
  list(w = w, status = status)
}


mixIP = function (matrix_lik, prior, pi_init = NULL, control = list()) {
  if(!requireNamespace("REBayes", quietly = TRUE)) {
    stop("mixIP requires installation of package REBayes")}
  control = set_control_mixIP(control)
  n = nrow(matrix_lik)
  k = ncol(matrix_lik)
  A = rbind(diag(length(prior)),matrix_lik) # add in observations corresponding to prior
  w = c(prior-1,rep(1,n))
  A = A[w!=0,] #remove zero weight entries, as these otherwise cause errors
  w = w[w!=0]
  #w = rep(1,n+k)
  res = REBayes::KWDual(A, rep(1,k), normalize(w), control=control)
  return(list(pihat = normalize(res$f), niter = NULL, converged=(res$status=="OPTIMAL"), control=control))
}

set_control_squarem=function(control,nobs){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if (nobs > 50000) control.default$trace = TRUE
  control.default$tol = min(0.1/nobs,1.e-7) # set default convergence criteria to be more stringent for larger samples
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

set_control_mixIP=function(control){
  control.default=list(rtol=1e-6)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

biopt = function (array_F, matrix_lik_z, pi_prior, w_prior, control, primal, gd.priority) {
  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                         tol = 1e-07, maxiter = 5000, trace = FALSE)
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)
  Kpi = dim(array_F)[1]
  Lw = dim(array_F)[2]
  pinw_init = c(1, rep(0, Kpi - 1), 1, rep(0, Lw - 1))
  res = SQUAREM::squarem(par = pinw_init, fixptfn = bifixpoint, objfn = binegpenloglik,
                array_F = array_F, matrix_lik_z = matrix_lik_z, pi_prior = pi_prior, w_prior = w_prior, primal = primal, gd.priority = gd.priority, control = controlinput)
  return(list(pihat = normalize(pmax(0, res$par[1 : Kpi])),
              what = res$par[-(1 : Kpi)],
              B = res$value.objfn,
              niter = res$fpevals,
              converged = res$convergence))
}

normalmix = function (pi, mean, sd) {
  structure(data.frame(pi, mean, sd), class = "normalmix")
}

binegpenloglik = function (pinw, array_F, matrix_lik_z, pi_prior, w_prior, primal, gd.priority)
{
  return(-bipenloglik(pinw, array_F, pi_prior, w_prior))
}

bipenloglik = function (pinw, array_F, pi_prior, w_prior) {
  K = dim(array_F)[1]
  pi = pinw[1 : K]
  w = pinw[-(1 : K)]
  loglik = sum(log(pmax(0, colSums(t(apply(pi * array_F, 2, colSums)) * w))))
  subset = (pi_prior != 1)
  priordens = sum((pi_prior - 1)[subset] * log(pi[subset]))
  return(loglik + priordens - sum(abs(w[-1]) * w_prior))
}

normalize = function (x) {
  return(x/sum(x))
}


array_f = function (betahat, sebetahat, sd, gd.ord, mixcompdist, gd.normalized) {
  if (mixcompdist == "normal") {
    array_f = array_f.normal(betahat, sebetahat, sd, gd.ord, gd.normalized)
  } else {
    stop ("invalid prior mixture")
  }
  return(array_f)
}

## this function is vectorized for x
## more efficient if let it run for x at once
gauss.deriv = function(x, ord) {
  return(dnorm(x) * EQL::hermite(x, ord))
}

Hermite = function (gd.ord) {
  x <- PolynomF::polynom()
  H <- PolynomF::polylist(x, - 1 + x^2)
  if (gd.ord >= 3) {
    for(n in 2 : (gd.ord - 1))
      H[[n+1]] <- x * H[[n]] - n * H[[n-1]]
  }
  return(H)
}

array_f.normal = function (betahat, sebetahat, sd, gd.ord, gd.normalized) {
  sd.mat = sqrt(outer(sebetahat^2, sd^2, FUN = "+"))
  beta.std.mat = betahat / sd.mat
  temp2 = array(dim = c(dim(beta.std.mat), gd.ord + 1))
  temp2[, , 1] = dnorm(beta.std.mat)
  hermite = Hermite(gd.ord)
  if (gd.normalized) {
    for (i in 1 : gd.ord) {
      temp2[, , i + 1] = temp2[, , 1] * hermite[[i]](beta.std.mat) * (-1)^i / sqrt(factorial(i))
    }
  } else {
    for (i in 1 : gd.ord) {
      temp2[, , i + 1] = temp2[, , 1] * hermite[[i]](beta.std.mat) * (-1)^i
    }
  }
  # temp2.test = outer(beta.std.mat, 0:gd.ord, FUN = gauss.deriv)
  se.std.mat = sebetahat / sd.mat
  temp1 = exp(outer(log(se.std.mat), 0 : gd.ord + 1, FUN = "*"))
  array_f = temp1 * temp2 / sebetahat
  rm(temp1)
  rm(temp2)
  return(array_f)
}

array_pm = function (betahat, sebetahat, sd, gd.ord, gd.normalized) {
  sd.mat = sqrt(outer(sebetahat^2, sd^2, FUN = "+"))
  beta.std.mat = betahat / sd.mat
  temp2 = array(dim = c(dim(beta.std.mat), gd.ord + 1))
  temp2_0 = dnorm(beta.std.mat)
  hermite = Hermite(gd.ord + 1)
  if (gd.normalized) {
    for (i in 0 : gd.ord) {
      temp2[, , i + 1] = temp2_0 * hermite[[i + 1]](beta.std.mat) * (-1)^(i + 1) / sqrt(factorial(i))
    }
  } else {
    for (i in 1 : (gd.ord + 1)) {
      temp2[, , i + 1] = temp2_0 * hermite[[i + 1]](beta.std.mat) * (-1)^(i + 1)
    }
  }
  # temp2.test = outer(beta.std.mat, 0:gd.ord, FUN = gauss.deriv)
  se.std.mat = sebetahat / sd.mat
  sd.std.mat2 = t(sd / t(sd.mat))^2
  temp1 = exp(outer(log(se.std.mat), 0 : gd.ord, FUN = "*"))
  for (ii in 0 : gd.ord) {
    temp1[, , ii + 1] = temp1[, , ii + 1] * sd.std.mat2
  }
  array_pm = (-1) * temp1 * temp2
  rm(temp1)
  rm(temp2)
  return(array_pm)
}

lfdr_top = function (pi0, w, betahat, sebetahat, gd.ord) {
  hermite = Hermite(gd.ord)
  gd.mat = matrix(0, ncol = gd.ord + 1, nrow = length(betahat))
  gd.mat[, 1] = dnorm(betahat / sebetahat)
  for (l in 1 : gd.ord) {
    gd.mat[, l + 1] = (-1)^l * hermite[[l]](betahat / sebetahat) * gd.mat[, 1] / sqrt(factorial(l))
  }
  gd.std.mat = gd.mat / sebetahat
  return(gd.std.mat %*% w * pi0)
}

array_PosProb = function (betahat, sebetahat, deltaAt0, sd, gd.ord, gd.normalized) {
  sd.mat = sqrt(outer(sebetahat^2, sd^2, FUN = "+"))
  beta.std.mat = betahat / sd.mat
  se.std.mat <- sebetahat / sd.mat
  sd.std.mat <- t(sd / t(sd.mat))
  sd.se.mat <- 1 / outer(sebetahat, sd, FUN = "/")
  beta.std.sd.se.mat <- beta.std.mat * sd.se.mat
  pdf.beta.std.mat <- dnorm(beta.std.mat)
  pdf.beta.std.sd.se.mat <- dnorm(beta.std.sd.se.mat)
  cdf.beta.std.sd.se.mat <- pnorm(beta.std.sd.se.mat)
  temp2 = array(dim = c(dim(beta.std.mat), gd.ord + 1))
  hermite = Hermite(gd.ord)
  temp2[, , 1] <- cdf.beta.std.sd.se.mat * pdf.beta.std.mat / sd.mat
  temp2[, , 2] <- cdf.beta.std.sd.se.mat * hermite[[1]](beta.std.mat) * (-1) * pdf.beta.std.mat * se.std.mat^2 / sebetahat +
    sd.se.mat * pdf.beta.std.sd.se.mat * pdf.beta.std.mat * se.std.mat^2 / sebetahat
  if (gd.normalized) {
    temp2[, , 3] <- (cdf.beta.std.sd.se.mat * hermite[[2]](beta.std.mat) * (-1)^2 * pdf.beta.std.mat * se.std.mat^3 / sebetahat +
        2 * sd.se.mat * se.std.mat^3 / sebetahat * pdf.beta.std.sd.se.mat * hermite[[1]](beta.std.mat) * (-1) * pdf.beta.std.mat +
        sd.se.mat^2 * se.std.mat^3 / sebetahat * hermite[[1]](beta.std.sd.se.mat) * (-1) * pdf.beta.std.sd.se.mat * pdf.beta.std.mat) / sqrt(factorial(2))
    for (i in 3 : gd.ord) {
      temp2[, , (i + 1)] <- cdf.beta.std.sd.se.mat * hermite[[i]](beta.std.mat) / sqrt(factorial(i)) * (-1)^i * pdf.beta.std.mat * se.std.mat^(i + 1) / sebetahat +
        i * sd.se.mat * se.std.mat^(i + 1) / sebetahat * pdf.beta.std.sd.se.mat * hermite[[i - 1]](beta.std.mat) / sqrt(factorial(i)) * (-1)^(i - 1) * pdf.beta.std.mat +
        sd.se.mat^i * se.std.mat^(i + 1) / sebetahat * hermite[[i - 1]](beta.std.sd.se.mat) / sqrt(factorial(i)) * (-1)^(i - 1) * pdf.beta.std.sd.se.mat * pdf.beta.std.mat +
        Reduce('+', lapply(2 : (i - 1), function (m) {
          sqrt(choose(i, m)) * hermite[[m - 1]](beta.std.sd.se.mat) / sqrt(factorial(m)) * pdf.beta.std.sd.se.mat * sd.std.mat^m *
            hermite[[i - m]](beta.std.mat) / sqrt(factorial(i - m)) * pdf.beta.std.mat * se.std.mat^(i - m) / sd.mat * (-1)^(i - 1)
        }))
    }
  } else {
    temp2[, , 3] <- temp3_1 * hermite[[2]](beta.std.mat) * (-1)^2 * temp2_0 +
      2 * temp3 * temp3_0 * hermite[[1]](beta.std.mat) * (-1) * temp2_0 +
      temp3^2 * hermite[[1]](beta.std.mat * temp3) * (-1) * temp3_0 * temp2_0
    for (i in 3 : gd.ord) {
      temp2[, , i + 1] <- ((temp3_1 * hermite[[i]](beta.std.mat) * (-1)^i +
                              i * temp3 * temp3_0 * hermite[[i - 1]](beta.std.mat) * (-1)^(i - 1) +
                              temp3^i * hermite[[i - 1]](beta.std.mat * temp3) * temp3_0) * temp2_0 +
                             sum(sapply(2 : (i - 1), function (m) {
                               exp(lchoose(i, m) + m * temp3) * hermite[[m - 1]](beta.std.mat * temp3) *
                                 hermite[[i - m]](beta.std.mat)
                             })) * temp3_0 * temp2_0) / sqrt(factorial(i))
    }
  }
  # temp2.test = outer(beta.std.mat, 0:gd.ord, FUN = gauss.deriv)
  if (deltaAt0) {
    array_PosProb = temp2
    array_PosProb[, 1, ] <- 0
  } else {
    array_PosProb = temp2
  }
  rm(temp2)
  return(array_PosProb)
}

autoselect.mixsd = function (betahat, sebetahat, mult)
{
  sebetahat = sebetahat[sebetahat != 0]
  sigmaamin = min(sebetahat)/10
  if (all(betahat^2 <= sebetahat^2)) {
    sigmaamax = 8 * sigmaamin
  }
  else {
    sigmaamax = 2 * sqrt(max(betahat^2 - sebetahat^2))
  }
  if (mult == 0) {
    return(c(0, sigmaamax/2))
  }
  else {
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

plot.ghat = function (fitted.g, mixcompdist = "normal", xlim = c(-10, 10)) {
  pi = fitted.g$pi
  mean = fitted.g$mean
  sd = fitted.g$sd
  x = seq(xlim[1], xlim[2], 0.01)
  if (mixcompdist == "normal") {
    y = sum(pi * pnorm(x, mean, sd))
  }
}

ghat.cdf = function (x, fitted.g) {
  pi = fitted.g$pi
  mean = fitted.g$mean
  sd = fitted.g$sd
  return(sum(pi * pnorm(x, mean, sd)))
}

qval.from.lfdr <- function (lfdr) {
  if (sum(!is.na(lfdr)) == 0) {
    return(rep(NA, length(lfdr)))
  }
  o = order(lfdr)
  qvalue = rep(NA, length(lfdr))
  qvalue[o] = (cumsum(sort(lfdr))/(1 : sum(!is.na(lfdr))))
  return(qvalue)
}
