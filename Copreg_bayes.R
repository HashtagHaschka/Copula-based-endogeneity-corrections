## =============================================================================
##  Copreg_bayes.R
##
##  Bayesian joint estimation of the Gaussian copula endogeneity correction,
##  Haschka (2025), Oxford Bulletin of Economics and Statistics.
##
##  Requires copreg-core.R to be sourced first.  Base R only.
## =============================================================================

if (!exists(".copreg_model"))
  stop("Source copreg-core.R before this file.")


## -----------------------------------------------------------------------------
##  Copula-based endogeneity corrections in R
##  https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  Copyright (C) 2026 Rouven E. Haschka
##  ORCID: https://orcid.org/0000-0002-2916-9745
##
##  If this code contributes to work you publish, please cite the software
##
##    Haschka, R. E. (2026). Copula-based endogeneity corrections in R.
##    https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  and the paper this estimator implements
##
##    Haschka, R. E. (2025). Bayesian inference for joint estimation models
##      using copulas to handle endogenous regressors. Oxford Bulletin of
##      Economics and Statistics.
##
##  This program is free software: you can redistribute it and/or modify it
##  under the terms of the GNU General Public License as published by the Free
##  Software Foundation, either version 3 of the License, or (at your option)
##  any later version.
##
##  This program is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
##  more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program.  If not, see <https://www.gnu.org/licenses/>.
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##  What this estimator does, and why it is its own class
## -----------------------------------------------------------------------------
##  Every other estimator in the toolbox computes marginal CDFs first, builds a
##  copula term from them, runs a regression and bootstraps the whole thing.
##  This one does none of that.  The CDFs of the regressors are not plugged in
##  but are unknown parameters: each regressor gets probability masses on its
##  uniquely observed values (Section 3.1), and those masses are drawn along
##  with the regression coefficients, the error variance and the full copula
##  correlation matrix in one MCMC run.  There is no first stage and nothing is
##  estimated a priori, so there are no plug-in estimates whose uncertainty
##  would have to be bootstrapped back in; the posterior carries all of it.
##
##  Consequences for the interface:
##
##    * cdf and ties do not exist.  The CDF is a parameter, and ties are the
##      normal case -- the masses sit on unique values, so a binary regressor
##      simply has two of them (the application in the paper has 121 weekly
##      dummies of exactly that kind).
##    * nboots is replaced by iterations, burnin and thin.
##    * rho IS reported here, unlike in 2sCOPE-np: the copula correlation
##      matrix is a parameter of the model, so its posterior comes for free.
##    * confint() returns credible intervals, not confidence intervals, and no
##      asymptotic argument is involved anywhere.
##
##  Ordering convention throughout, taken from the paper: z (endogenous), then
##  x (exogenous), then e.  The intercept is a regressor but not part of the
##  copula -- it has one unique value, so it carries no distributional
##  information (Web Appendix E does the same).


## =============================================================================
##  1.  Base R replacements for invgamma, mvtnorm, copula and LaplacesDemon
## =============================================================================
##  The reference implementation pulls in four packages for draws and densities
##  that are a few lines each.  Writing them out keeps the dependency list at
##  Formula, and the copula density in particular gets much faster: the
##  Gaussian copula log density is given in closed form in Equation (3), so
##  there is no reason to route N observations through a generic dCopula().

.bayes_rinvgamma <- function(n, shape, rate) 1 / stats::rgamma(n, shape, rate = rate)

## Wishart by the Bartlett decomposition.  V = LL' lower Cholesky, A lower
## triangular with chi distributed diagonal, then LA(LA)' ~ W(nu, V).
.bayes_rwish <- function(nu, V) {
  p <- ncol(V)
  L <- t(chol(V))
  A <- matrix(0, p, p)
  diag(A) <- sqrt(stats::rchisq(p, nu - seq_len(p) + 1))
  if (p > 1L) A[lower.tri(A)] <- stats::rnorm(p * (p - 1L) / 2L)
  tcrossprod(L %*% A)
}

## X ~ IW(nu, S)  <=>  X^{-1} ~ W(nu, S^{-1}).  Scale parameterisation, i.e.
## E(X) = S / (nu - p - 1), which is what the paper's Equation (8) uses.
.bayes_riwish <- function(nu, S) {
  Si <- chol2inv(chol(S))
  chol2inv(chol(.bayes_rwish(nu, Si)))
}

## Multivariate normal draw and log density, both from one Cholesky factor.
.bayes_rmvn <- function(mu, R) as.vector(mu + crossprod(R, stats::rnorm(ncol(R))))
.bayes_dmvn <- function(x, mu, R) {
  z <- backsolve(R, x - mu, transpose = TRUE)
  -0.5 * length(x) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z * z)
}


## =============================================================================
##  2.  Margins: probability masses on the unique values of a regressor
## =============================================================================
##  Section 3.1.  For a generic regressor w the density is treated as a
##  nonparametric function putting mass lambda_j on the j-th smallest unique
##  value, and the margin entering the copula is the cumulative sum.
##
##  Two details that the reference implementations disagree on, resolved here
##  in favour of the version that survives ties:
##
##    * cells are the unique VALUES, not the observations.  With m unique
##      values out of N observations the mass vector has length m, and the
##      counts v enter the Dirichlet full conditional.  Indexing by rank()
##      instead breaks as soon as two observations share a value, which for
##      the dummies of the empirical application is every observation.
##    * the margin is evaluated at the midpoint of the mass, cumsum(l) - l/2,
##      not at its upper edge.  For a continuous regressor with all masses of
##      order 1/N the difference is negligible; for a binary one it is the
##      difference between a sensible normal score and Phi^{-1}(1) = Inf.

.bayes_margin <- function(v) {
  uv  <- sort(unique(v))
  grp <- match(v, uv)
  m   <- length(uv)
  list(uv = uv, grp = grp, m = m, cnt = tabulate(grp, nbins = m),
       f = factor(grp, levels = seq_len(m)))
}

## normal scores of the margin implied by the current masses
.bayes_xi <- function(lambda, mg) {
  Fmid <- (cumsum(lambda) - lambda / 2) * (mg$m / (mg$m + 0.01))
  stats::qnorm(Fmid[mg$grp])
}

## Gibbs draw of lambda, Web Appendix B following Ng et al. (2011): a copula
## correlated uniform per observation, pushed through the Gamma(1,1) quantile
## function, aggregated within value cells, plus one Gamma(1,1) prior
## pseudo-count per cell for the Dir(1) prior.  Normalising gives a draw from
## Dir(1 + v).
.bayes_draw_lambda <- function(u, mg) {
  g <- stats::qgamma(u, shape = 1, rate = 1)
  G <- as.numeric(tapply(g, mg$f, sum))
  G[is.na(G)] <- 0
  G <- G + stats::rgamma(mg$m, shape = 1, rate = 1)
  G / sum(G)
}


## =============================================================================
##  3.  Log posterior, score and working weight
## =============================================================================
##  Only the parts that depend on (alpha, beta, delta) and sigma^2 are needed:
##  the two Metropolis-Hastings steps condition on Sigma and lambda, so the
##  categorical likelihood terms and the Dirichlet priors of Equation (4) are
##  constant and cancel in the acceptance ratio.  Dropping them is not an
##  approximation, it is the same ratio computed with fewer operations.
##
##  Writing A = Sigma^{-1} - I and splitting xi = (xi_zx, xi_e), the copula
##  exponent is
##
##      xi' A xi  =  xi_zx' A_zx,zx xi_zx  +  2 xi_e (A_zx,e' xi_zx)  +  A_ee xi_e^2
##
##  whose first term is constant within the step.  Precomputing q = A_zx,e' xi_zx
##  once turns each evaluation into O(N) instead of O(N d^2).

.bayes_lpost <- function(cf, s2, y, X, q, aee, pr, horseshoe, phi2) {
  if (s2 <= 0) return(-Inf)
  e   <- y - as.vector(X %*% cf)
  xe  <- e / sqrt(s2)
  q2  <- sum(xe * xe)
  lp  <- -0.5 * (2 * sum(q * xe) + aee * q2) -
    0.5 * length(e) * log(s2) - 0.5 * q2
  ## inverse Gamma prior on sigma^2, Equation (7)
  lp <- lp - (pr$a + 1) * log(s2) - pr$b / s2
  ## prior on the regression coefficients, intercept included
  lp + if (horseshoe) {
    ## The horseshoe of Equation (5) with phi, tau and iota marginalised out.
    ## That marginal has no closed form; Carvalho, Polson & Scott (2010,
    ## Thm. 1) bracket it by
    ##
    ##   (K/2) log(1 + 4/gamma^2)  <  p(gamma)  <  K log(1 + 2/gamma^2),
    ##
    ## and the upper bound is the usual working form.  Both ends are right:
    ## p behaves like 1/gamma^2 far out, which is the Cauchy tail Equation (5)
    ## puts there, and diverges at zero, which is what does the shrinking.
    ##
    ## Reading Equation (6) literally, as p(gamma) proportional to 1/|gamma|,
    ## does not work.  That density is not integrable at zero, so the spike
    ## there is infinitely deep rather than merely tall.  The proposal moves
    ## the whole coefficient vector at once, so as soon as one coefficient
    ## lands near zero the prior term explodes, the move is accepted, and
    ## every later move away from it is rejected: acceptance collapses to a
    ## fraction of a percent and the trace becomes a flat line.  The log-log
    ## divergence below shrinks just as hard and is integrable.
    sum(log(log1p(2 / pmax(cf * cf, 1e-300))))
  } else {
    sum(stats::dnorm(cf, 0, sqrt(phi2), log = TRUE))
  }
}

## First and second derivative of the log full conditional of tau = log sigma^2.
## Web Appendix C, plus the Jacobian of the transformation: the target in tau
## is p(sigma^2) * sigma^2, so the prior contributes -(a+1) + b/sigma^2 + 1,
## which is the -a + b/sigma^2 below.  The acceptance ratio must carry the
## matching (tau_new - tau_cur); the reference implementation omits it there
## while using this score, which biases sigma^2.
.bayes_tau_derivs <- function(cf, s2, y, X, q, aee, pr) {
  e   <- y - as.vector(X %*% cf)
  xe  <- e / sqrt(s2)
  qxe <- sum(q * xe)
  q2  <- sum(xe * xe)
  list(score = 0.5 * qxe + 0.5 * (aee + 1) * q2 - 0.5 * length(e) -
         pr$a + pr$b / s2,
       f2 = -0.25 * qxe - 0.5 * (aee + 1) * q2 - pr$b / s2)
}


## =============================================================================
##  4.  Block-recursive Gibbs step for the copula covariance matrix
## =============================================================================
##  Exogeneity of x means Cov(xi_x, xi_e) = 0.  Web Appendix A of the paper
##  draws W from a plain inverse Wishart full conditional, which does not
##  respect those zeros: the inverse Wishart family has no way to hold entries
##  of a covariance matrix at zero.
##
##  Instead W is reparameterised hierarchically,
##
##      xi_z = B_x xi_x + beta_e xi_e + eps,
##      xi_x ~ MN(0, Sigma_xx),  xi_e ~ N(0, sigma_e^2),  eps ~ MN(0, Omega),
##
##  with the three components independent.  The implied W has the zero blocks
##  by construction, and the map from (B_x, beta_e, Omega, Sigma_xx, sigma_e^2)
##  to W is bijective given the restrictions, so nothing is lost.  With
##  C = [B_x | beta_e] the four conjugate full conditionals are
##
##      Sigma_xx | .  ~  IW(nu1 + N, Psi1 + S_xx)
##      sigma_e^2 | . ~  IG(a + N/2, b + S_ee/2)
##      C | Omega, .  ~  MN(Mbar, Omega, Vbar)
##      Omega | C, .  ~  IW(nu2 + N, Psi2 + S_epseps)
##
##  drawn in that order, C using the Omega of the previous cycle.  Sigma is
##  then the correlation matrix of the reconstructed W.

.bayes_draw_W <- function(xi_z, xi_x, xi_e, Omega, pr) {

  N <- length(xi_e)
  K <- ncol(xi_z)
  L <- if (is.null(xi_x)) 0L else ncol(xi_x)

  Xt <- if (L > 0L) cbind(xi_x, xi_e) else matrix(xi_e, N, 1L)
  St <- crossprod(Xt)
  Sz <- crossprod(xi_z, Xt)

  Sxx <- if (L > 0L) .bayes_riwish(pr$nu1 + N, pr$Psi1 + crossprod(xi_x)) else NULL
  se2 <- .bayes_rinvgamma(1, pr$a_e + N / 2, pr$b_e + sum(xi_e^2) / 2)

  Vb <- chol2inv(chol(pr$V0inv + St))
  Vb <- (Vb + t(Vb)) / 2
  Mb <- (pr$M0 %*% pr$V0inv + Sz) %*% Vb
  C  <- Mb + t(chol(Omega)) %*% matrix(stats::rnorm(K * (L + 1L)), K, L + 1L) %*%
    chol(Vb)

  R     <- xi_z - Xt %*% t(C)
  Omega <- .bayes_riwish(pr$nu2 + N, pr$Psi2 + crossprod(R))

  Bx <- C[, seq_len(L), drop = FALSE]
  be <- C[, L + 1L, drop = FALSE]

  W11 <- Omega + se2 * tcrossprod(be)
  if (L > 0L) W11 <- W11 + Bx %*% Sxx %*% t(Bx)
  W13 <- se2 * be
  W <- if (L > 0L) {
    W12 <- Bx %*% Sxx
    rbind(cbind(W11, W12, W13),
          cbind(t(W12), Sxx, matrix(0, L, 1L)),
          cbind(t(W13), matrix(0, 1L, L), se2))
  } else {
    rbind(cbind(W11, W13), cbind(t(W13), se2))
  }

  s <- sqrt(diag(W))
  list(Sigma = W / outer(s, s), Omega = Omega)
}


## =============================================================================
##  5.  The sampler
## =============================================================================
##  Algorithm 1 of Web Appendix D: one Metropolis-Hastings step for the
##  regression coefficients using the iteratively weighted least squares
##  proposal of Gamerman (1997), one for log sigma^2 using a Newton proposal,
##  then Gibbs steps for W and for the masses.
##
##  Draws are thinned as they are produced rather than afterwards.  The
##  reference implementation keeps every iterate and thins at the end, which
##  for the empirical application of the paper -- 1,002,000 iterations and some
##  2,600 mass columns -- would need around 21 GB.  Keeping only the retained
##  iterates brings that down to a few MB.

.bayes_sampler <- function(y, X, cop, mgs, iterations, burnin, thin,
                           horseshoe, pr, start, verbose) {

  N  <- length(y)
  p  <- ncol(X)
  K  <- length(cop$endo)
  L  <- length(cop$exog)
  d  <- K + L + 1L
  ci <- seq_len(d - 1L)                       # copula columns other than e

  ndraw <- length(seq.int(burnin + 1L, iterations, by = thin))
  out <- list(
    coefficients = matrix(NA_real_, ndraw, p, dimnames = list(NULL, colnames(X))),
    sigma2       = numeric(ndraw),
    Sigma        = matrix(NA_real_, ndraw, d * (d - 1L) / 2L),
    lambda       = lapply(mgs, function(m) matrix(NA_real_, ndraw, m$m)),
    accept       = c(coefficients = 0, sigma2 = 0))
  names(out$lambda) <- names(mgs)

  cf   <- start$coefficients
  s2   <- start$sigma2
  Sig  <- start$Sigma
  lam  <- start$lambda
  phi2 <- start$phi2
  Om   <- diag(K)

  ## The same progress bar the bootstrap uses.  It cannot be redrawn on every
  ## iteration here: at 102,000 of them that is 102,000 writes to the console
  ## for a bar some 190 characters wide, which costs more than it shows.
  ## Redrawing at most a thousand times is beyond the resolution of the bar
  ## and free.
  keep <- 0L
  nxt  <- burnin + 1L
  step <- max(1L, iterations %/% 1000L)
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (it in seq_len(iterations)) {

    ## --- normal scores of the margins under the current masses --------------
    Xi <- matrix(0, N, d)
    for (j in seq_len(d - 1L)) Xi[, j] <- .bayes_xi(lam[[j]], mgs[[j]])

    Sinv <- chol2inv(chol(Sig))
    A    <- Sinv - diag(d)
    q    <- as.vector(Xi[, ci, drop = FALSE] %*% A[ci, d])
    aee  <- A[d, d]

    ## --- Metropolis-Hastings for the regression coefficients ---------------
    ## IWLS proposal.  The working weight is the exact observed information
    ## per observation, (Sigma^{-1})_{dd} / sigma^2, so the proposal covariance
    ## does not depend on the coefficients and reusing it in the reverse
    ## density is exact rather than approximate.
    e  <- y - as.vector(X %*% cf)
    sg <- sqrt(s2)
    nu <- (q + aee * e / sg) / sg + e / s2
    Rp <- chol(chol2inv(chol((Sinv[d, d] / s2) * crossprod(X))))

    mu  <- cf + as.vector(crossprod(Rp, Rp %*% crossprod(X, nu)))
    cfp <- .bayes_rmvn(mu, Rp)

    ep  <- y - as.vector(X %*% cfp)
    nup <- (q + aee * ep / sg) / sg + ep / s2
    mup <- cfp + as.vector(crossprod(Rp, Rp %*% crossprod(X, nup)))

    la <- .bayes_lpost(cfp, s2, y, X, q, aee, pr, horseshoe, phi2) -
      .bayes_lpost(cf, s2, y, X, q, aee, pr, horseshoe, phi2) +
      .bayes_dmvn(cf, mup, Rp) - .bayes_dmvn(cfp, mu, Rp)
    if (is.finite(la) && log(stats::runif(1)) < la) {
      cf <- cfp
      out$accept[["coefficients"]] <- out$accept[["coefficients"]] + 1
    }

    ## --- Gibbs for the coefficient variances, normal prior only ------------
    if (!horseshoe)
      phi2 <- .bayes_rinvgamma(p, pr$a_phi + 0.5, pr$b_phi + cf^2 / 2)

    ## --- Metropolis-Hastings for log sigma^2 -------------------------------
    tau <- log(s2)
    dv  <- .bayes_tau_derivs(cf, s2, y, X, q, aee, pr)
    Pc  <- if (dv$f2 < -1e-12) -1 / dv$f2 else 1
    mc  <- tau + Pc * dv$score
    tp  <- stats::rnorm(1, mc, sqrt(Pc))
    s2p <- exp(tp)

    dvp <- .bayes_tau_derivs(cf, s2p, y, X, q, aee, pr)
    Pp  <- if (dvp$f2 < -1e-12) -1 / dvp$f2 else 1
    mp  <- tp + Pp * dvp$score

    la <- .bayes_lpost(cf, s2p, y, X, q, aee, pr, horseshoe, phi2) -
      .bayes_lpost(cf, s2, y, X, q, aee, pr, horseshoe, phi2) +
      stats::dnorm(tau, mp, sqrt(Pp), log = TRUE) -
      stats::dnorm(tp, mc, sqrt(Pc), log = TRUE) +
      (tp - tau)                                  # Jacobian of tau = log sigma^2
    if (is.finite(la) && log(stats::runif(1)) < la) {
      s2 <- s2p
      out$accept[["sigma2"]] <- out$accept[["sigma2"]] + 1
    }

    ## --- Gibbs for W, hence Sigma ------------------------------------------
    ## xi_e is e / sigma directly.  Routing it through qnorm(pnorm(e, sd = s))
    ## is the same number in exact arithmetic but returns +-Inf once |e/sigma|
    ## exceeds about 8.3, where pnorm saturates at 0 or 1.
    Xi[, d] <- (y - as.vector(X %*% cf)) / sqrt(s2)
    dw  <- .bayes_draw_W(Xi[, seq_len(K), drop = FALSE],
                         if (L > 0L) Xi[, K + seq_len(L), drop = FALSE] else NULL,
                         Xi[, d], Om, pr)
    Sig <- dw$Sigma
    Om  <- dw$Omega

    ## --- Gibbs for the masses ----------------------------------------------
    Rs <- chol(Sig)
    U  <- stats::pnorm(matrix(stats::rnorm(N * d), N, d) %*% Rs)
    for (j in seq_len(d - 1L)) lam[[j]] <- .bayes_draw_lambda(U[, j], mgs[[j]])

    ## --- store -------------------------------------------------------------
    if (it == nxt) {
      keep <- keep + 1L
      out$coefficients[keep, ] <- cf
      out$sigma2[keep]         <- s2
      out$Sigma[keep, ]        <- Sig[lower.tri(Sig)]
      for (j in seq_len(d - 1L)) out$lambda[[j]][keep, ] <- lam[[j]]
      nxt <- nxt + thin
    }

    if (verbose && (it %% step == 0L || it == iterations))
      utils::setTxtProgressBar(pb, it)
  }

  out$accept <- out$accept / iterations
  out
}


## =============================================================================
##  6.  CopRegBAYES
## =============================================================================
##
##  formula        y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##                 as everywhere else in the toolbox: factors, interactions,
##                 transformations and "-1" all work, and the position around
##                 the "|" decides what is treated as endogenous.  Every
##                 non-intercept column of the design matrix enters the copula,
##                 so a four level factor contributes three columns with two
##                 unique values each, exactly as the weekly dummies of the
##                 paper's application do.
##
##  iterations,    102,000 iterations, the first 2,000 discarded and every
##  burnin, thin   100th of the rest retained, which is the setting of the
##                 paper's simulations and leaves 1,000 draws.  Draws are
##                 thinned while sampling, so only the retained ones are ever
##                 held in memory.
##
##  horseshoe      TRUE uses the horseshoe prior of Equation (5) in the
##                 marginalised form of Equation (6), p(gamma) proportional to
##                 1/|gamma|, on every coefficient including the intercept.
##                 FALSE puts an independent N(0, phi^2) on each coefficient
##                 with phi^2 ~ IG, drawn by a conjugate Gibbs step.
##
##  prior.args     list overriding the hyperparameters, all defaulting to the
##                 paper: a and b for sigma^2 (0.001 each, Section 3.2.2),
##                 Psi1, nu1, Psi2, nu2, a_e, b_e, M0 and V0 for the blocks of
##                 W, and a_phi, b_phi for the normal prior.
##
##  start          list with any of coefficients, sigma2, Sigma, lambda.
##                 Defaults follow Web Appendix D: OLS for the coefficients and
##                 the error variance, the identity for Sigma, and one draw
##                 from Dir(1) per regressor for the masses.
##
##  Use validity() for the identification checks and the convergence
##  diagnostics, including Gelman-Rubin.

CopRegBAYES <- function(formula, data,
                        iterations = 102000,
                        burnin = 2000,
                        thin = 100,
                        horseshoe = TRUE,
                        prior.args = list(),
                        start = NULL,
                        subset = NULL,
                        contrasts = NULL,
                        verbose = interactive(),
                        cdf, ties, nboots) {

  cl <- match.call()

  if (!missing(cdf) || !missing(ties))
    stop("The Bayesian approach estimates the marginal CDFs rather than ",
         "plugging them in: the\n  probability masses defining them are ",
         "parameters, drawn along with everything else.\n  The 'cdf' and ",
         "'ties' arguments therefore do not apply here.", call. = FALSE)
  if (!missing(nboots))
    stop("Inference here is posterior rather than bootstrap. Use 'iterations',",
         " 'burnin' and\n  'thin' instead of 'nboots'.", call. = FALSE)

  iterations <- as.integer(iterations)
  burnin     <- as.integer(burnin)
  thin       <- as.integer(thin)
  if (burnin >= iterations)
    stop("'burnin' must be smaller than 'iterations'.", call. = FALSE)
  if (thin < 1L) stop("'thin' must be at least 1.", call. = FALSE)

  info <- .copreg_model(formula, data, subset, contrasts)
  X <- info$X; y <- info$y
  N <- length(y); p <- ncol(X)

  ## --- which columns enter the copula, in the order z, x, e -----------------
  endo <- info$endo_cols
  exog <- setdiff(info$exo_cols, if (info$has_intercept) 1L else integer(0))
  K <- length(endo); L <- length(exog)
  if (K < 1L) stop("No endogenous regressor found.", call. = FALSE)

  const <- vapply(c(endo, exog), function(j) length(unique(X[, j])) < 2L,
                  logical(1))
  if (any(const))
    stop("Regressor(s) with a single unique value cannot enter the copula: ",
         paste(colnames(X)[c(endo, exog)][const], collapse = ", "),
         ".\n  A constant carries no distributional information.", call. = FALSE)

  ## P for the endogenous regressors, W for the exogenous ones and xi for the
  ## structural error, as everywhere else in the toolbox; the paper writes
  ## z, x and e for the same three things.  The stars mark normal scores:
  ## these are correlations of Phi^{-1}(F(.)), which is what the frequentist
  ## estimators report as rho(P*, xi*) as well.
  cop <- list(endo = endo, exog = exog,
              names = c(colnames(X)[endo], colnames(X)[exog], "xi"))
  mgs <- lapply(c(endo, exog), function(j) .bayes_margin(X[, j]))
  names(mgs) <- colnames(X)[c(endo, exog)]
  d <- K + L + 1L

  ## --- hyperparameters, defaults as in the paper ---------------------------
  pr <- list(a = 0.001, b = 0.001,                       # sigma^2, Equation (7)
             a_phi = 0.001, b_phi = 0.001,               # normal prior variance
             nu1 = L + 2, Psi1 = diag(max(L, 1L)),       # Sigma_xx
             nu2 = K + 2, Psi2 = diag(K),                # Omega
             a_e = 0.001, b_e = 0.001,                   # sigma_e^2
             M0 = matrix(0, K, L + 1L), V0 = diag(L + 1L))
  ## An explicit NULL means "leave the default", not "delete it": assigning
  ## NULL into a list drops the element, which would leave the sampler without
  ## a hyperparameter it needs.
  prior.args <- prior.args[!vapply(prior.args, is.null, logical(1))]
  if (length(prior.args) > 0L) {
    bad <- setdiff(names(prior.args), names(pr))
    if (length(bad) > 0L)
      stop("Unknown entries in 'prior.args': ", paste(bad, collapse = ", "),
           ".\n  Available: ", paste(names(pr), collapse = ", "), ".",
           call. = FALSE)
    pr[names(prior.args)] <- prior.args
  }
  pr$V0inv <- chol2inv(chol(pr$V0))

  ## --- starting values, Web Appendix D --------------------------------------
  ols <- stats::lm.fit(X, y)
  st <- list(coefficients = unname(ols$coefficients),
             sigma2 = sum(ols$residuals^2) / max(1L, N - p),
             Sigma = diag(d),
             lambda = lapply(mgs, function(m)
               { g <- stats::rgamma(m$m, 1); g / sum(g) }),
             phi2 = rep(1000, p))
  start <- start[!vapply(start, is.null, logical(1))]
  if (length(start) > 0L) {
    bad <- setdiff(names(start), names(st))
    if (length(bad) > 0L)
      stop("Unknown entries in 'start': ", paste(bad, collapse = ", "), ".",
           call. = FALSE)
    st[names(start)] <- start
  }
  if (anyNA(st$coefficients))
    stop("The design matrix is rank deficient, so OLS gives no starting ",
         "values.\n  Drop the collinear column(s) or supply 'start'.",
         call. = FALSE)

  if (verbose)
    message("Sampling ", format(iterations, big.mark = ","), " iterations, ",
            "keeping ", length(seq.int(burnin + 1L, iterations, by = thin)),
            " draws ...")

  fit <- .bayes_sampler(y, X, cop, mgs, iterations, burnin, thin,
                        horseshoe, pr, st, verbose)

  ## --- labels for the free elements of Sigma --------------------------------
  ij <- which(lower.tri(diag(d)), arr.ind = TRUE)
  colnames(fit$Sigma) <- paste0("rho(", cop$names[ij[, 2L]], "*, ",
                                cop$names[ij[, 1L]], "*)")
  free <- !(ij[, 2L] > K & ij[, 2L] < d & ij[, 1L] == d)

  draws <- cbind(fit$coefficients, sigma2 = fit$sigma2,
                 fit$Sigma[, free, drop = FALSE])

  ## A stuck parameter is a warning sign about the run, so it is raised here
  ## rather than waiting for validity() to be called.  An acceptance rate near
  ## zero says the same thing about the whole vector at once: the Metropolis
  ## step is proposing moves that are never taken, and the draws describe a
  ## point rather than a posterior.
  frozen <- colnames(draws)[apply(draws, 2, function(v) all(v == v[1L]))]
  if (length(frozen) > 0L)
    warning("The chain never moved for: ", paste(frozen, collapse = ", "),
            ".\n  These draws describe a single point, not a posterior. ",
            "Check the trace with\n  plot(fit, which = ", deparse(frozen[1L]),
            ", type = \"trace\") before reading anything off them.",
            call. = FALSE)
  if (fit$accept[["coefficients"]] < 0.01)
    warning("The coefficients were accepted in only ",
            formatC(100 * fit$accept[["coefficients"]], format = "f",
                    digits = 2), "% of iterations.\n  The sampler is barely ",
            "moving; the draws are unlikely to represent the posterior.",
            call. = FALSE)

  cf <- colMeans(fit$coefficients)
  names(cf) <- colnames(X)
  fitted <- as.vector(X %*% cf)

  structure(list(
    coefficients = cf,
    posterior.median = apply(fit$coefficients, 2, stats::median),
    std.error = apply(fit$coefficients, 2, stats::sd),
    sigma2 = mean(fit$sigma2),
    vcov = stats::var(fit$coefficients),
    draws = draws,
    coefficient.draws = fit$coefficients,
    sigma2.draws = fit$sigma2,
    Sigma.draws = fit$Sigma,
    Sigma.free = free,
    lambda.draws = fit$lambda,
    margins = mgs,
    copula.names = cop$names,
    acceptance = fit$accept,
    fitted.values = fitted,
    residuals = y - fitted,
    y = y, X = X, n = N,
    formula = formula, terms = info$terms, mf = info$mf,
    xlevels = info$xlevels, contrasts = info$contrasts,
    na.action = info$na.action,
    endo.names = colnames(X)[endo], exog.names = colnames(X)[exog],
    rho.names = paste0("rho(", colnames(X)[endo], "*, xi*)"),
    horseshoe = horseshoe, prior.args = pr,
    iterations = iterations, burnin = burnin, thin = thin,
    ndraws = nrow(fit$coefficients),
    method = "Bayesian copula correction (Haschka 2025)",
    call = cl), class = "copregbayes")
}


## =============================================================================
##  7.  Methods
## =============================================================================

print.copregbayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nPosterior means:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n", x$ndraws, " draws from ", format(x$iterations, big.mark = ","),
      " iterations.\n", sep = "")
  invisible(x)
}


summary.copregbayes <- function(object, level = 0.95, ...) {
  a <- (1 - level) / 2
  tab <- t(apply(object$draws, 2, function(v)
    c(`P. Mean` = mean(v), `P. Median` = stats::median(v), `Sd` = stats::sd(v),
      stats::quantile(v, c(a, 1 - a)))))
  ## Only the correlations between an endogenous regressor and the structural
  ## error are shown.  Those are the endogeneity, and the whole point of the
  ## model; the correlations among the regressors themselves are nuisance
  ## parameters that the joint estimation happens to produce as well, and with
  ## twelve regressors there are sixty-six of them.  They stay on the object.
  p   <- length(object$coefficients)
  rho <- tab[rownames(tab) %in% object$rho.names, , drop = FALSE]
  structure(list(call = object$call, method = object$method,
                 coefficients = tab[seq_len(p), , drop = FALSE],
                 sigma2 = tab[p + 1L, , drop = FALSE],
                 rho = rho,
                 n.other = nrow(tab) - p - 1L - nrow(rho),
                 restricted = colnames(object$Sigma.draws)[!object$Sigma.free],
                 acceptance = object$acceptance, horseshoe = object$horseshoe,
                 n = object$n, ndraws = object$ndraws, level = level,
                 iterations = object$iterations, burnin = object$burnin,
                 thin = object$thin),
            class = "summary.copregbayes")
}


print.summary.copregbayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nRegression coefficients:\n")
  print(format(x$coefficients, digits = digits), quote = FALSE)
  cat("\nStructural error variance:\n")
  print(format(x$sigma2, digits = digits), quote = FALSE)
  cat("\nEndogeneity: rho(P*, xi*) is the correlation between the normal score",
      "\n  of an endogenous regressor and that of the structural error.\n")
  print(format(x$rho, digits = digits), quote = FALSE)
  if (x$n.other > 0L)
    cat("  ", x$n.other, " further correlations among the regressors are ",
        "estimated jointly and\n  sit in $Sigma.draws; ", length(x$restricted),
        " more are held at zero, which is what makes\n  the exogenous ",
        "regressors exogenous.\n", sep = "")
  cat("\nThe two quantile columns are ", format(100 * x$level), "% credible ",
      "intervals: no asymptotic\n  argument is involved, and the marginal ",
      "CDFs are estimated jointly rather\n  than plugged in.\n", sep = "")
  cat("\n", x$n, " observations. ", x$ndraws, " draws kept from ",
      format(x$iterations, big.mark = ","), " iterations (burn-in ",
      format(x$burnin, big.mark = ","), ", thinning ", x$thin, ").\n", sep = "")
  cat("Prior on the coefficients: ",
      if (x$horseshoe) "horseshoe" else "normal with inverse Gamma variance",
      ". Acceptance ", paste0(names(x$acceptance), " ",
                              formatC(x$acceptance, format = "f", digits = 3),
                              collapse = ", "), ".\n", sep = "")
  invisible(x)
}


coef.copregbayes    <- function(object, ...) object$coefficients
vcov.copregbayes    <- function(object, ...) object$vcov
nobs.copregbayes    <- function(object, ...) object$n
formula.copregbayes <- function(x, ...) x$formula
fitted.copregbayes  <- function(object, ...) object$fitted.values

residuals.copregbayes <- function(object, ...) object$residuals

## Credible intervals, the posterior quantiles themselves rather than a normal
## approximation to them.
confint.copregbayes <- function(object, parm, level = 0.95, ...) {
  a <- (1 - level) / 2
  ci <- t(apply(object$draws, 2, stats::quantile, c(a, 1 - a)))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

predict.copregbayes <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  mt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(mt, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X  <- stats::model.matrix(mt, mf, contrasts.arg = object$contrasts)
  as.vector(X[, names(object$coefficients), drop = FALSE] %*% object$coefficients)
}


## -----------------------------------------------------------------------------
##  plot
## -----------------------------------------------------------------------------
##  'which' selects what to look at, by name or by position, because a model of
##  any size has far more parameters than fit on a screen; the default is the
##  endogenous regressors, which are the ones the correction is about.
##
##  type = "cdf" is different in kind from the other four.  It shows the
##  marginal distribution of a regressor as estimated by the model -- the
##  cumulative sums of the probability masses, with a pointwise credible band
##  across the draws.  That object simply does not exist in the frequentist
##  estimators, where the CDF is a fixed plug-in with no uncertainty attached.

plot.copregbayes <- function(x, which = NULL,
                             type = c("trace", "acf", "pacf", "density", "cdf"),
                             level = 0.95, ask = NULL, ...) {

  type <- match.arg(type)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (type == "cdf") {
    nm <- names(x$lambda.draws)
    if (is.null(which)) which <- x$endo.names
    if (is.numeric(which)) which <- nm[which]
    bad <- setdiff(which, nm)
    if (length(bad) > 0L)
      stop("No margin is estimated for: ", paste(bad, collapse = ", "),
           ".\n  Available: ", paste(nm, collapse = ", "), ".", call. = FALSE)
    a <- (1 - level) / 2
    if (is.null(ask)) ask <- length(which) > 1L && grDevices::dev.interactive()
    graphics::par(ask = ask)
    for (v in which) {
      L <- x$lambda.draws[[v]]
      F <- t(apply(L, 1, cumsum))
      q <- apply(F, 2, stats::quantile, c(a, 0.5, 1 - a))
      u <- x$margins[[v]]$uv
      graphics::plot(u, q[2, ], type = "s", ylim = c(0, 1), xlab = v,
                     ylab = "posterior CDF",
                     main = paste0("Estimated marginal CDF of ", v), ...)
      graphics::lines(u, q[1, ], type = "s", lty = 2, col = "grey40")
      graphics::lines(u, q[3, ], type = "s", lty = 2, col = "grey40")
      graphics::rug(x$margins[[v]]$uv, col = "grey60")
      graphics::legend("bottomright", bty = "n", lty = c(1, 2),
                       col = c("black", "grey40"),
                       legend = c("posterior median",
                                  paste0(format(100 * level), "% band")))
    }
    return(invisible(x))
  }

  nm <- colnames(x$draws)
  if (is.null(which)) which <- x$endo.names
  if (is.numeric(which)) which <- nm[which]
  bad <- setdiff(which, nm)
  if (length(bad) > 0L)
    stop("Not a parameter of this model: ", paste(bad, collapse = ", "),
         ".\n  Available: ", paste(nm, collapse = ", "), ".", call. = FALSE)

  if (is.null(ask)) ask <- length(which) > 1L && grDevices::dev.interactive()
  graphics::par(ask = ask)
  for (v in which) {
    d <- x$draws[, v]
    switch(type,
           trace   = graphics::plot(seq_along(d), d, type = "l",
                                    xlab = "draw", ylab = v,
                                    main = paste("Trace of", v), ...),
           acf     = stats::acf(d, main = paste("ACF of", v), ...),
           pacf    = stats::pacf(d, main = paste("PACF of", v), ...),
           density = { dd <- stats::density(d)
                       graphics::plot(dd, xlab = v,
                                      main = paste("Posterior of", v), ...)
                       graphics::abline(v = mean(d), lty = 2) })
  }
  invisible(x)
}


## -----------------------------------------------------------------------------
##  validity
## -----------------------------------------------------------------------------
##  Two things at once, both of which the paper treats as prerequisites rather
##  than results.
##
##  Identification, Section 4.2.  The three known failure modes are normality
##  of the endogenous regressor, non-Gaussian regressor-error dependence, and a
##  misspecified structural error distribution.  Only the first is testable
##  from the data alone, and it is the one that matters most: if z is normal
##  then E(e|z) is linear and the copula cannot separate it from the regression
##  itself.  The paper's own reassurance is that the Bayesian version degrades
##  gracefully -- the posterior of rho(z,e) concentrates near zero rather than
##  producing confident nonsense -- so that posterior is reported alongside.
##
##  Convergence.  Geweke, effective sample size and the lag-one autocorrelation
##  come out of the chain that is already there.  Gelman-Rubin does not: it
##  needs several chains from dispersed starts, which means running the sampler
##  again.  Web Appendix E draws those starts from N(OLS, OLS se) for the
##  coefficients, LKJ for the correlation matrix and Dir(1) for the masses, and
##  that is what chains = TRUE does.  It costs as many runs as chains.

.bayes_geweke <- function(v, first = 0.1, last = 0.5) {
  n <- length(v)
  a <- v[seq_len(floor(first * n))]
  b <- v[(n - floor(last * n) + 1L):n]
  s <- function(u) {
    m <- max(1L, floor(10 * log10(length(u))))
    r <- stats::acf(u, lag.max = m, plot = FALSE)$acf[-1L]
    stats::var(u) * (1 + 2 * sum(r * (1 - seq_along(r) / (m + 1)))) / length(u)
  }
  se <- s(a) + s(b)
  if (!is.finite(se) || se <= 0) return(NA_real_)
  (mean(a) - mean(b)) / sqrt(se)
}

## Effective sample size by the initial positive sequence: sum the sample
## autocorrelations up to the first negative one and stop there.
##
## The truncation has to fire when the FIRST autocorrelation is already
## negative, which is the hardest truncation of all -- no lags at all, so
## ESS = n.  Excluding that case leaves the whole sequence out to lag
## floor(10 log10 n) in the sum, noise terms included, and 1 + 2 sum(r) can
## then approach zero from above and send the estimate far past n.  For a well
## mixing chain the sign of the lag-one estimate is close to a coin flip, so
## this is the common case rather than a corner: the diagnostic was least
## reliable exactly where the chain was best behaved, and it erred towards
## claiming more precision than there was.
.bayes_ess <- function(v) {
  n <- length(v)
  m <- max(1L, floor(10 * log10(n)))
  r <- stats::acf(v, lag.max = m, plot = FALSE)$acf[-1L]
  if (!all(is.finite(r))) return(NA_real_)      # a series that never moved
  k <- which(r < 0)[1L]
  if (!is.na(k)) r <- r[seq_len(k - 1L)]
  max(1, n / (1 + 2 * sum(r)))
}

## LKJ(1), the uniform distribution over correlation matrices, by the onion
## method; used only for dispersed starting values.
.bayes_rlkj <- function(d) {
  if (d < 2L) return(diag(d))
  b <- (d - 1) / 2
  r <- 2 * stats::rbeta(1, b, b) - 1
  R <- matrix(c(1, r, r, 1), 2L, 2L)
  if (d > 2L) for (m in 2:(d - 1L)) {
    b <- b - 0.5
    y <- stats::rbeta(1, m / 2, b)
    u <- stats::rnorm(m); u <- u / sqrt(sum(u * u))
    z <- t(chol(R)) %*% (sqrt(y) * u)
    R <- rbind(cbind(R, z), c(as.vector(z), 1))
  }
  R
}

validity.copregbayes <- function(object, chains = FALSE, power = 0.8,
                                 verbose = interactive(), ...) {

  nz <- .validity_nonnormality(object$X[, object$endo.names, drop = FALSE],
                               object$endo.names, object$n, power)

  ## posterior of the endogeneity correlations, one per endogenous regressor
  ze <- object$rho.names[object$rho.names %in% colnames(object$draws)]
  rho <- t(apply(object$draws[, ze, drop = FALSE], 2, function(v)
    c(`P. Mean` = mean(v), `2.5%` = unname(stats::quantile(v, 0.025)),
      `97.5%` = unname(stats::quantile(v, 0.975)),
      `P(rho > 0)` = mean(v > 0))))

  conv <- data.frame(
    Geweke = apply(object$draws, 2, .bayes_geweke),
    ESS    = apply(object$draws, 2, .bayes_ess),
    `AC(1)` = apply(object$draws, 2, function(v) {
      a <- stats::acf(v, lag.max = 1, plot = FALSE)$acf[2L]
      if (is.finite(a)) a else NA_real_
    }),
    check.names = FALSE)

  ## A parameter whose draws never move is not a missing diagnostic, it is a
  ## finding: the Metropolis step never accepted for it.  The three columns
  ## above are all NA for such a parameter -- acf() of a constant is NaN --
  ## so without naming it the row would read as though something merely could
  ## not be computed.
  stuck <- colnames(object$draws)[
    apply(object$draws, 2, function(v) all(v == v[1L]))]

  gr <- NULL
  if (!isFALSE(chains)) {
    nc <- if (isTRUE(chains)) 4L else as.integer(chains)
    if (nc < 2L) stop("Gelman-Rubin needs at least two chains.", call. = FALSE)
    ## The caller's frame has to be taken here, in validity()'s own body.
    ## Inside the lapply() below, parent.frame() would be the frame of the
    ## lapply() call instead, and the refit would look for the model's data
    ## anywhere but where it lives -- so it worked only when the data happened
    ## to be a global variable, and silently used the wrong object when a
    ## global of the same name existed.
    where <- parent.frame()
    cl <- object$call
    cl$verbose <- verbose
    se <- object$std.error
    if (verbose)
      message("Gelman-Rubin needs one further run of the sampler per chain; ",
              nc, " to go.")
    reps <- lapply(seq_len(nc), function(i) {
      if (verbose) message("Chain ", i + 1L, " of ", nc + 1L, ":")
      cl$start <- list(
        coefficients = stats::rnorm(length(object$coefficients),
                                    object$coefficients, se),
        Sigma = .bayes_rlkj(length(object$copula.names)),
        lambda = lapply(object$margins, function(m)
          { g <- stats::rgamma(m$m, 1); g / sum(g) }))
      eval(cl, where)$draws
    })
    reps <- c(list(object$draws), reps)
    n <- nrow(reps[[1L]]); M <- length(reps)
    gr <- vapply(seq_len(ncol(reps[[1L]])), function(j) {
      x <- vapply(reps, function(r) r[, j], numeric(n))
      B <- n * stats::var(colMeans(x))
      W <- mean(apply(x, 2, stats::var))
      if (W <= 0) return(NA_real_)
      sqrt(((n - 1) / n * W + B / n) / W)
    }, numeric(1))
    names(gr) <- colnames(object$draws)
  }

  structure(list(nonnormality = nz$table, thresholds = nz$thresholds,
                 endogeneity = rho, convergence = conv, gelman.rubin = gr,
                 stuck = stuck,
                 acceptance = object$acceptance, ndraws = object$ndraws,
                 endo.names = object$endo.names),
            class = "copregbayes.validity")
}


print.copregbayes.validity <- function(x, digits = max(3L, getOption("digits") - 3L),
                                       ...) {
  i <- 0L; nxt <- function() paste0("[", i <<- i + 1L, "] ")
  cat("\nChecks for the Bayesian copula correction\n")

  cat("\n", nxt(), "Nonnormality of the endogenous regressors\n", sep = "")
  cat("    Identification comes from the nonlinearity of E(e|z), which is\n",
      "    present only when z is nonnormal; under normality the copula cannot\n",
      "    tell regressor variation from error variation.\n", sep = "")
  print(format(x$nonnormality, digits = digits))

  cat("\n", nxt(), "Posterior of the endogeneity correlations\n", sep = "")
  print(format(x$endogeneity, digits = digits), quote = FALSE)
  cat("    An interval covering zero says the data carry no evidence of\n",
      "    endogeneity, which is the honest reading; it is not a test decision.\n",
      sep = "")

  cat("\n", nxt(), "Convergence of the chain\n", sep = "")
  cat("    Geweke compares the first tenth with the last half of the draws and\n",
      "    is a z statistic, so |Geweke| > 2 is a warning sign. ESS is the\n",
      "    effective number of independent draws behind ", x$ndraws,
      " retained ones.\n", sep = "")
  print(format(x$convergence, digits = digits))
  if (length(x$stuck) > 0L)
    cat("    The chain never moved for: ", paste(x$stuck, collapse = ", "),
        ".\n    That is why those rows are NA, and it is a finding rather than",
        " a gap:\n    the sampler is not exploring the posterior of these",
        " parameters at all.\n", sep = "")
  cat("    Acceptance rates: ",
      paste0(names(x$acceptance), " ",
             formatC(x$acceptance, format = "f", digits = 3), collapse = ", "),
      "\n", sep = "")

  if (!is.null(x$gelman.rubin)) {
    cat("\n", nxt(), "Gelman-Rubin from dispersed starts\n", sep = "")
    print(format(x$gelman.rubin, digits = digits), quote = FALSE)
    cat("    Gelman et al. read convergence as GR <= 1.1.\n")
  } else {
    cat("\n", nxt(), "Gelman-Rubin was not computed\n", sep = "")
    cat("    It needs several chains from dispersed starts, so it costs one\n",
        "    further run of the sampler per chain. Pass chains = TRUE, or a\n",
        "    number, to compute it.\n", sep = "")
  }

  cat("\nSource: Haschka (2025), Oxford Bulletin of Economics and Statistics\n")
  invisible(x)
}


### REFERENCES ----------------------------------------------------------------
##
## Carvalho, C. M., N. G. Polson, and J. G. Scott (2010). The horseshoe
##   estimator for sparse signals. Biometrika 97(2), 465-480.
##
## Gamerman, D. (1997). Sampling from the posterior distribution in generalized
##   linear mixed models. Statistics and Computing 7(1), 57-68.
##
## Haschka, R. E. (2025). Bayesian inference for joint estimation models using
##   copulas to handle endogenous regressors. Oxford Bulletin of Economics and
##   Statistics.
##
## Ng, K. W., G. Tian, and M. Tang (2011). Dirichlet and Related Distributions.
##   Chichester: Wiley.
##
## Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##   estimation using copulas. Marketing Science 31(4), 567-586.
##
## Zhang, X., W. J. Boscardin, and T. R. Belin (2006). Sampling correlation
##   matrices in Bayesian models with correlated latent variables. Journal of
##   Computational and Graphical Statistics 15(4), 880-896.
