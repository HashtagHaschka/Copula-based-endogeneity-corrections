## =============================================================================
##  CopRegPANEL.R
##
##  Copula correction for linear panel models with fixed effects and correlated
##  regressors, Haschka (2022), Journal of Marketing Research 59(4), 861-880.
##
##  This is a different model class from the cross-sectional estimators: it is
##  estimated by maximum likelihood rather than augmented least squares, and it
##  is deliberately not registered with copreg(), just as plm is separate from
##  lm.  It does share the CDF estimators and the general conventions of
##  copreg-core.R, which has to be sourced first.
## =============================================================================

if (!exists(".copreg_fit"))
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
##    Haschka, R. E. (2022). Handling endogenous regressors using copulas: A
##      generalization to linear panel models with fixed effects and correlated
##      regressors. Journal of Marketing Research 59(4), 861-880.
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


## =============================================================================
##  1.  Panel infrastructure
## =============================================================================

## Panel-aware lag, lead and diff.  These are placed in the environment the
## formula is evaluated in, so that lag(y), lag(y, 2), lag(x):w and the like
## work inside the formula exactly as any other function would, without turning
## the data into a special class first.  Values never cross a panel boundary.
.panel_env <- function(id, parent = parent.frame()) {
  
  e <- new.env(parent = parent)
  
  e$lag <- function(x, k = 1L) {
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 0L)
      stop("lag(): 'k' must be a single non-negative integer.", call. = FALSE)
    if (k == 0L) return(x)
    n <- length(x)
    out <- rep(NA_real_, n)
    i <- seq_len(n)
    ok <- i > k
    ok[ok] <- id[i[ok]] == id[i[ok] - k]
    out[ok] <- x[i[ok] - k]
    out
  }
  
  e$lead <- function(x, k = 1L) {
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 0L)
      stop("lead(): 'k' must be a single non-negative integer.", call. = FALSE)
    if (k == 0L) return(x)
    n <- length(x)
    out <- rep(NA_real_, n)
    i <- seq_len(n)
    ok <- i + k <= n
    ok[ok] <- id[i[ok]] == id[i[ok] + k]
    out[ok] <- x[i[ok] + k]
    out
  }
  
  e$diff <- function(x, k = 1L) x - e$lag(x, k)
  
  e
}


## Forward orthogonal deviations (Arellano 1993).  Footnote 4 of the paper:
## combining a fixed-effects transformation with the GLS transformation that
## removes the resulting nonspherical errors yields exactly FOD, so applying
## FOD to Equation 3 lands directly on Equation 10.
##
##   w*_it = sqrt(k / (k + 1)) * (w_it - mean of the remaining w_i,t+1..T),
##   k = T_i - t,  t = 1, ..., T_i - 1
##
## The original code built this as chol(solve(D D')) D with D the differencing
## matrix, which is the same transformation up to an overall sign -- the sign
## cancels because y and every regressor are transformed alike -- but costs
## O(T^3) per panel instead of O(T).
##
## The last observation of every panel is consumed, and panels with a single
## observation drop out entirely.
.fod <- function(M, id) {
  
  M  <- as.matrix(M)
  sp <- split(seq_len(nrow(M)), factor(id, levels = unique(id)))
  
  keep <- unlist(lapply(sp, function(ix)
    if (length(ix) >= 2L) ix[-length(ix)] else integer(0)), use.names = FALSE)
  
  out <- matrix(NA_real_, length(keep), ncol(M),
                dimnames = list(NULL, colnames(M)))
  
  pos <- 0L
  for (ix in sp) {
    Ti <- length(ix)
    if (Ti < 2L) next
    B <- M[ix, , drop = FALSE]
    ## reverse cumulative sums: S[t, ] = sum over s >= t
    S <- apply(B, 2L, function(v) rev(cumsum(rev(v))))
    if (is.null(dim(S))) S <- matrix(S, Ti, ncol(B))
    tt <- seq_len(Ti - 1L)
    k  <- Ti - tt
    out[pos + tt, ] <- sqrt(k / (k + 1)) *
      (B[tt, , drop = FALSE] - S[tt + 1L, , drop = FALSE] / k)
    pos <- pos + Ti - 1L
  }
  
  list(M = out, keep = keep)
}


## Formula and data processing.  Reuses .copreg_model() so that the formula
## semantics, the endogenous/exogenous split and the diagnostics are identical
## to the cross-sectional estimators; everything panel specific happens around
## it.
## The 'intercept' argument of CopRegPANEL() is not passed down here: it acts
## on the transformed model, which does not exist yet at this point.  What
## happens here is the opposite -- the intercept column is built and then
## dropped again, see below.
.panel_model <- function(formula, data, index, subset = NULL,
                         contrasts = NULL) {
  
  if (!is.data.frame(data)) data <- as.data.frame(data)
  
  ## --- index ---------------------------------------------------------------
  if (missing(index) || is.null(index))
    stop("'index' must name the panel identifier, and normally also the time ",
         "variable: index = c(\"firm\", \"year\").", call. = FALSE)
  if (!is.character(index) || length(index) > 2L)
    stop("'index' must be one or two variable names.", call. = FALSE)
  miss <- setdiff(index, names(data))
  if (length(miss) > 0L)
    stop("Index variable(s) not found in the data: ",
         paste(miss, collapse = ", "), call. = FALSE)
  
  idv <- index[1L]
  tmv <- if (length(index) == 2L) index[2L] else NA_character_
  
  if (!is.null(subset)) data <- data[subset, , drop = FALSE]
  
  ## --- sort by (id, time) --------------------------------------------------
  ## FOD and lag() both read the rows in order, so this is not cosmetic.  With
  ## no time variable the given row order within a panel is taken as the time
  ## order and the user is told so.
  if (is.na(tmv)) {
    message("No time variable in 'index'. The existing row order within each ",
            "panel is\n  taken as the time order.")
    ord <- order(data[[idv]])
  } else {
    if (anyDuplicated(data[, c(idv, tmv)]))
      stop("The index does not identify the rows uniquely: some (",
           idv, ", ", tmv, ") pair occurs more than once.", call. = FALSE)
    ord <- order(data[[idv]], data[[tmv]])
  }
  data <- data[ord, , drop = FALSE]
  
  ## --- panel-aware lag/lead/diff in the formula ----------------------------
  F <- Formula::as.Formula(formula)
  environment(F) <- .panel_env(data[[idv]], environment(formula))
  
  ## --- intercept -----------------------------------------------------------
  ## The model matrix is always built WITH an intercept, and its column is
  ## dropped again below.  This is not cosmetic: without an intercept
  ## model.matrix codes the first factor with one dummy per level, whose sum is
  ## the constant, and the transformation sends the constant to zero, leaving a
  ## rank-deficient design.  Building with an intercept gives the factors
  ## treatment contrasts, which is the correct coding once the individual
  ## effects absorb the constant.
  ##
  ## A "- 1" in the formula therefore has no effect here; the 'intercept'
  ## argument controls the transformed model instead, see CopRegPANEL().
  Fi <- if (length(F)[2L] >= 2L) stats::update(F, . ~ . + 1 | .)
  else stats::update(F, . ~ . + 1)
  environment(Fi) <- environment(F)
  
  ## .copreg_model checks the endogenous regressors for ties, but on this
  ## model the copula is applied to the transformed variable, not the raw one.
  ## Its warning would therefore describe the wrong scale, so it is muffled
  ## here and CopRegPANEL() repeats the check after the transformation.
  info <- withCallingHandlers(
    .copreg_model(Fi, data, subset = NULL, contrasts = contrasts),
    warning = function(w) {
      if (grepl("substantial ties", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  
  ic <- match("(Intercept)", colnames(info$X))
  if (!is.na(ic)) {
    info$X <- info$X[, -ic, drop = FALSE]
    info$endo_cols <- info$endo_cols - (info$endo_cols > ic)
    info$exo_cols  <- setdiff(seq_len(ncol(info$X)), info$endo_cols)
    info$part <- info$part[-ic]; info$ord <- info$ord[-ic]
    info$factor_cols <- info$factor_cols[info$factor_cols != ic]
    info$factor_cols <- info$factor_cols - (info$factor_cols > ic)
    info$has_intercept <- FALSE
  }
  
  ## --- rows that survived na.omit ------------------------------------------
  keep <- if (is.null(info$na.action)) seq_len(nrow(data))
  else seq_len(nrow(data))[-unclass(info$na.action)]
  id   <- data[[idv]][keep]
  tm   <- if (is.na(tmv)) seq_along(keep) else data[[tmv]][keep]
  
  c(info, list(id = id, time = tm, index = index, rows = keep,
               data = data[keep, , drop = FALSE], Formula = F))
}


## =============================================================================
##  2.  Likelihood
## =============================================================================
##
##  Equation 12.  With K endogenous regressors x, L exogenous regressors z and
##  the error e of the FOD-transformed model,
##
##      L(theta) = prod_it det(Xi)^(-1/2)
##                 exp[ -1/2 (xi_x, xi_z, xi_e)' (Xi^-1 - I) (xi_x, xi_z, xi_e) ]
##                 x phi(e_it, sigma^2),
##
##  where xi = Phi^-1(Fhat(.)) are the normal scores of the regressors and
##  xi_e = Phi^-1(Phi(e, sigma^2)).  The first line is the Gaussian copula
##  density; the second is the marginal density of the error, which is what
##  ties beta into the likelihood once e is replaced by y - x'beta - z'delta.
##
##  Xi is assembled rather than estimated as a whole:
##
##    - the block among the regressors is observed, so the empirical
##      correlation of their normal scores is plugged in
##    - the correlations between the exogenous regressors and the error are
##      restricted to zero, that being what exogeneity means here
##    - only the K correlations between the endogenous regressors and the error
##      are free parameters
##
##  That is what keeps the model parsimonious: adding exogenous regressors to
##  the copula costs no additional parameters, and it is exactly the step that
##  removes the omitted-variable bias the paper's appendix describes.
##
##  Parameterisation for the optimiser: rho = tanh(par) and sigma^2 = exp(par),
##  so the search is unconstrained.

## 'cop' indexes the columns of X that are random variables and therefore enter
## the copula; an intercept is excluded.  'endo_pos' gives the positions of the
## endogenous regressors within cop, i.e. the rows of Xi carrying a free rho.
##
## Two simplifications make this exact and cheap.
##
## First, the normal score of the error needs no transformation at all:
##
##   xi_e = Phi^-1( Phi(e, sigma^2) ) = Phi^-1( Phi(e / sigma) ) = e / sigma,
##
## because Phi(., sigma^2) in Equation 11 is the CDF of N(0, sigma^2) while
## Phi^-1 is the standard normal quantile.  Computing it as a round trip
## through pnorm and qnorm, as the original code did, needs the probabilities
## clipped away from 0 and 1, which truncates the likelihood at |e/sigma| of
## about 7 and leaves a kink there.  The direct form has neither problem.
##
## Second, the regressor block S of Xi is the empirical correlation of the
## normal scores and does not depend on any parameter, so with
##
##   Xi = [ S  r ]        v = S^-1 r,   k = 1 - r'v,
##        [ r' 1 ]
##
##   det(Xi) = det(S) k,
##   z'(Xi^-1 - I)z = q0 + (zr'v)^2 / k - 2 xi_e (zr'v) / k + xi_e^2 (1/k - 1),
##   q0 = zr' S^-1 zr - zr'zr,
##
## the parts that cost O(n m^2) are computed once and every likelihood
## evaluation is left with an O(n m) matrix-vector product.  Positive
## definiteness of Xi reduces to S being positive definite, checked once, and
## k > 0.
##
## .panel_precompute() does the fixed work; .panel_negll() is what the
## optimiser calls.

.panel_precompute <- function(Zreg, Xi_reg) {
  m  <- ncol(Zreg)
  ch <- tryCatch(chol(Xi_reg), error = function(e) NULL)
  if (is.null(ch))
    stop("The empirical correlation matrix of the copula data is not positive ",
         "definite.", call. = FALSE)
  Sinv <- chol2inv(ch)
  list(Zreg = Zreg, Sinv = Sinv, logdetS = 2 * sum(log(diag(ch))),
       q0 = rowSums((Zreg %*% Sinv) * Zreg) - rowSums(Zreg^2), m = m)
}

.panel_negll <- function(par, y, X, pre, endo_pos, fix_rho = FALSE) {
  
  p <- ncol(X)
  m <- pre$m
  K <- length(endo_pos)
  n <- length(y)
  
  beta <- par[seq_len(p)]
  rho  <- if (fix_rho) numeric(K) else tanh(par[p + seq_len(K)])
  ls2  <- par[length(par)]
  if (!is.finite(ls2) || ls2 > 700 || ls2 < -700) return(1e10)
  s2 <- exp(ls2)
  
  e <- y - as.vector(X %*% beta)
  if (!all(is.finite(e))) return(1e10)
  
  r <- numeric(m); r[endo_pos] <- rho
  v <- as.vector(pre$Sinv %*% r)
  k <- 1 - sum(r * v)
  if (!is.finite(k) || k <= 1e-10) return(1e10)      # Xi not positive definite
  
  ze <- e / sqrt(s2)
  g  <- as.vector(pre$Zreg %*% v)
  
  quad <- pre$q0 + g * g / k - 2 * ze * g / k + ze * ze * (1 / k - 1)
  
  ll <- -0.5 * (n * (pre$logdetS + log(k)) + sum(quad)) -
    0.5 * n * (log(2 * pi) + ls2) - 0.5 * sum(ze * ze)
  
  if (!is.finite(ll)) return(1e10)
  -ll
}


## Optimiser.  BFGS is quick when it works but leans on numerical gradients of
## a surface that is not always smooth, which is why the original code reached
## for Nelder-Mead.  The default tries BFGS and silently falls back to
## Nelder-Mead whenever BFGS errors, returns a non-finite value or fails to
## converge; which one was used is recorded in the fitted object.
.panel_optim <- function(start, method = c("BFGS", "Nelder-Mead"),
                         maxit = 5000L, ...) {
  
  method <- match.arg(method)
  
  run <- function(m, st, mx)
    tryCatch(stats::optim(st, .panel_negll, method = m,
                          control = list(maxit = mx, reltol = 1e-10), ...),
             error = function(e) NULL)
  
  o <- run(method, start, maxit)
  ok <- !is.null(o) && is.finite(o$value) && o$value < 1e9 && o$convergence == 0L
  
  if (!ok && method == "BFGS") {
    o2 <- run("Nelder-Mead", start, max(maxit, 20000L))
    if (!is.null(o2) && is.finite(o2$value) && o2$value < 1e9) {
      o2$method <- "Nelder-Mead"
      o2$fallback <- TRUE
      return(o2)
    }
  }
  
  if (is.null(o) || !is.finite(o$value) || o$value >= 1e9)
    return(NULL)
  
  o$method <- method
  o$fallback <- FALSE
  o
}


## One full estimation pass on a given (already FOD-transformed) sample.
.panel_estimate <- function(y, X, endo, cdf, ties, method, start = NULL,
                            maxit = 5000L, restricted = TRUE, cop = NULL) {
  
  p <- ncol(X); K <- length(endo)
  if (is.null(cop)) cop <- seq_len(p)
  endo_pos <- match(endo, cop)
  m <- length(cop)
  
  U <- matrix(NA_real_, nrow(X), m,
              dimnames = list(NULL, colnames(X)[cop]))
  for (j in seq_len(m)) U[, j] <- .cdf_estimate(X[, cop[j]], cdf, ties)
  Zreg <- stats::qnorm(pmin(pmax(U, 1e-12), 1 - 1e-12))
  Xi_reg <- if (m == 1L) matrix(1, 1, 1) else stats::cor(Zreg)
  if (m > 1L) {
    ev <- eigen(Xi_reg, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < 1e-8)
      stop("The empirical correlation matrix of the copula data is singular. ",
           "Two regressors carry the same information after the ",
           "transformation; this happens in particular with a full set of ",
           "time dummies.", call. = FALSE)
  }
  pre <- .panel_precompute(Zreg, Xi_reg)
  
  if (is.null(start)) {
    f0 <- .lm.fit(X, y)
    if (f0$rank < p) return(NULL)
    start <- c(f0$coefficients, numeric(K),
               log(max(sum(f0$residuals^2) / max(length(y) - p, 1L), 1e-8)))
  }
  
  o <- .panel_optim(start, method = method, maxit = maxit,
                    y = y, X = X, pre = pre, endo_pos = endo_pos)
  if (is.null(o)) return(NULL)
  
  ## restricted fit, all rho = 0, for the likelihood ratio test.  Skipped in
  ## the bootstrap, where only the point estimates are needed.
  o0 <- if (!restricted) NULL else
    .panel_optim(c(o$par[seq_len(p)], o$par[length(o$par)]),
                 method = method, maxit = maxit,
                 y = y, X = X, pre = pre, endo_pos = endo_pos, fix_rho = TRUE)
  
  list(beta   = o$par[seq_len(p)],
       rho    = tanh(o$par[p + seq_len(K)]),
       sigma2 = exp(o$par[length(o$par)]),
       logLik = -o$value,
       logLik0 = if (is.null(o0)) NA_real_ else -o0$value,
       convergence = o$convergence, method = o$method,
       fallback = o$fallback, par = o$par,
       Zreg = Zreg, Xi_reg = Xi_reg)
}


## =============================================================================
##  3.  Bootstrap
## =============================================================================
##
##  Equation 26.  Cross-sectional units are resampled, not rows: draw
##  I_1, ..., I_N i.i.d. uniform on {1, ..., N} and take the whole block
##  (y_i, x_i, z_i) of each drawn unit.  Keeping the units intact is what
##  preserves the within-panel dependence, and it is what makes the standard
##  errors robust to serial correlation and to heteroskedasticity
##  (Goncalves & Kilian 2004, cited in the paper).  Resampling rows instead
##  breaks exactly that and produces standard errors that are too small.
##
##  The asymptotics run in N, so a small cross section makes these imprecise.
##
##  Each replicate refits the whole model, re-estimating the marginal CDFs on
##  the resampled data: they are generated regressors, and treating them as
##  given is the reason the information matrix cannot be used here.

## Which columns enter the copula?
##
## Not every regressor is a random variable whose marginal distribution the
## copula can describe.  An intercept is constant.  Time dummies are worse:
## they are deterministic functions of the period and identical for every
## panel, so the transformation maps them to the same values throughout, and
## once centred a full set of them spans one dimension less than it has
## columns.  The regression is unaffected -- the design matrix keeps full rank
## -- but the correlation matrix of the copula data becomes singular.
##
## The redundant columns are therefore dropped from the copula while staying in
## the model.  The endogenous regressors are ordered first so that they are
## never the ones dropped; if one of them is redundant the model is not
## identified and that is reported instead.
##
## The selection is made once, on the full sample, and reused in every
## bootstrap replicate so that the parameter vector keeps its length.
.panel_copula_cols <- function(Xv, endo, cdf, ties, intercept) {
  
  p    <- ncol(Xv)
  cand <- if (isTRUE(intercept)) seq_len(p)[-1L] else seq_len(p)
  if (length(cand) == 0L)
    stop("No regressor left to build the copula from.", call. = FALSE)
  
  Z <- vapply(cand, function(j)
    stats::qnorm(pmin(pmax(.cdf_estimate(Xv[, j], cdf, ties), 1e-12),
                      1 - 1e-12)), numeric(nrow(Xv)))
  Zc <- sweep(Z, 2L, colMeans(Z), "-")
  
  ep   <- match(endo, cand)
  ordv <- c(ep, setdiff(seq_along(cand), ep))
  qz   <- qr(Zc[, ordv, drop = FALSE], tol = 1e-7)
  sel  <- sort(ordv[qz$pivot[seq_len(qz$rank)]])
  
  if (!all(ep %in% sel))
    stop("An endogenous regressor carries no independent variation in the ",
         "copula data after the transformation, so its dependence with the ",
         "error is not identified.", call. = FALSE)
  
  drop <- setdiff(seq_along(cand), sel)
  if (length(drop) > 0L)
    message("Kept in the model but left out of the copula, because they carry ",
            "no independent\n  variation after the transformation (typically ",
            "a full set of time dummies): ",
            paste(colnames(Xv)[cand[drop]], collapse = ", "), ".")
  
  cand[sel]
}


.panel_boot <- function(y, X, id, endo, cdf, ties, method, nboots, start,
                        cop = NULL, verbose = interactive(), parallel = FALSE,
                        ncores = NULL, max_tries = 50L) {
  
  sp <- split(seq_along(y), factor(id, levels = unique(id)))
  N  <- length(sp)
  p  <- ncol(X); K <- length(endo)
  
  one <- function(...) {
    for (a in seq_len(max_tries)) {
      ix <- unlist(sp[sample.int(N, N, replace = TRUE)], use.names = FALSE)
      Xb <- X[ix, , drop = FALSE]
      if (qr(Xb)$rank < p) next
      r <- tryCatch(.panel_estimate(y[ix], Xb, endo, cdf, ties, method,
                                    start = start, restricted = FALSE,
                                    cop = cop),
                    error = function(e) NULL)
      if (is.null(r)) next
      ## the same resample without the copula, i.e. the plain within estimator,
      ## so that the standard errors are directly comparable
      fo <- tryCatch(.lm.fit(Xb, y[ix]), error = function(e) NULL)
      if (is.null(fo) || fo$rank < p) next
      return(c(r$beta, r$rho, r$sigma2, fo$coefficients))
    }
    NULL
  }
  
  ptype <- if (isFALSE(parallel) || is.null(parallel)) "no"
  else if (isTRUE(parallel))
    if (.Platform$OS.type == "unix") "multicore" else "snow"
  else match.arg(parallel, c("no", "multicore", "snow"))
  
  if (ptype != "no") {
    if (is.null(ncores)) {
      nc <- parallel::detectCores()
      ncores <- if (is.na(nc)) 1L else max(1L, nc - 1L)
    }
    ncores <- max(1L, min(as.integer(ncores), nboots))
    if (ncores == 1L) ptype <- "no"
  }
  
  if (ptype != "no") {
    
    if (verbose)
      message("Bootstrapping standard errors (", nboots, " replicates on ",
              ncores, " cores, ", ptype, ") ...")
    
    if (ptype == "multicore") {
      reps <- parallel::mclapply(seq_len(nboots), one, mc.cores = ncores)
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      env <- environment(.panel_boot)
      parallel::clusterExport(cl, .copreg_export_names(env), envir = env)
      parallel::clusterSetRNGStream(cl)
      reps <- parallel::parLapply(cl, seq_len(nboots), one)
    }
    bad  <- vapply(reps, function(r) is.null(r) || inherits(r, "try-error"),
                   logical(1))
    fail <- sum(bad)
    out  <- if (all(bad)) NULL else do.call(rbind, reps[!bad])
    
  } else {
    
    if (verbose) {
      message("Bootstrapping standard errors (", nboots, " replicates) ...")
      pb <- utils::txtProgressBar(min = 0, max = nboots, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    out <- matrix(NA_real_, nboots, 2L * p + K + 1L)
    fail <- 0L
    for (b in seq_len(nboots)) {
      r <- one()
      if (is.null(r)) fail <- fail + 1L else out[b, ] <- r
      if (verbose) utils::setTxtProgressBar(pb, b)
    }
    out <- out[stats::complete.cases(out), , drop = FALSE]
  }
  
  if (is.null(out) || nrow(out) == 0L)
    stop("No bootstrap replicate could be estimated. With very few panels or ",
         "a weakly identified model the likelihood often fails to converge on ",
         "resampled data.", call. = FALSE)
  if (fail > 0L)
    warning(fail, " of ", nboots, " bootstrap replicates did not converge and ",
            "were dropped.", call. = FALSE)
  
  colnames(out) <- c(colnames(X), paste0("rho.", colnames(X)[endo]), "sigma2",
                     paste0("within.", colnames(X)))
  out
}



## =============================================================================
##  4.  CopRegPANEL
## =============================================================================
##
##  formula     y ~ endogenous | exogenous, the same syntax as the other
##              estimators.  Everything lm() accepts is allowed, plus the
##              panel-aware lag(), lead() and diff():
##
##                y ~ lag(y) + price | promo + as.factor(year)
##                y ~ price + lag(price, 2) | lag(promo) + x:lag(x)
##
##              A lagged dependent variable belongs in the endogenous part.
##              FOD does not make it predetermined: FOD(y_{t-1}) averages over
##              future values of the lagged series, which contain the current
##              error, so the Nickell bias survives the transformation intact.
##              That is the situation the paper's dynamic panel case addresses.
##
##  data        a data.frame
##
##  index       c("id", "time"), as in plm.  The time variable is used to sort;
##              gaps in it are accepted and the observed sequence is treated as
##              consecutive.
##
##  intercept   FALSE by default.  A fixed-effects transformation removes
##              anything time invariant, so an intercept is identified only
##              beside time dummies, where it stands for the reference period.
##
##  cdf         marginal CDF estimator, default "kde.plugin": a Gaussian kernel
##              with the Polansky & Baker plug-in bandwidth, the estimator the
##              original implementation used.  The paper itself uses the density
##              estimator of Ngom, Moussa & De Dieu Nkurunziza (2018) integrated
##              by quadrature, and its footnote 5 explicitly allows conventional
##              kernel estimation or the empirical CDF instead; any choice from
##              copreg-core.R is accepted.
##
##  method      "BFGS" with an automatic Nelder-Mead fallback, or "Nelder-Mead".
##
##  Returns an object of class "copregpanel".

CopRegPANEL <- function(formula, data, index,
                        cdf = "kde.plugin",
                        ties = "max",
                        intercept = FALSE,
                        nboots = 199,
                        method = c("BFGS", "Nelder-Mead"),
                        start = NULL,
                        maxit = 5000L,
                        subset = NULL,
                        contrasts = NULL,
                        parallel = FALSE,
                        ncores = NULL,
                        verbose = interactive()) {
  
  cl     <- match.call()
  cdf    <- match.arg(cdf, .copreg_cdf_choices)
  ties   <- match.arg(ties, .copreg_ties_choices)
  method <- match.arg(method)
  if (!is.numeric(nboots) || length(nboots) != 1L || nboots < 2)
    stop("'nboots' must be a single number >= 2.", call. = FALSE)
  nboots <- as.integer(nboots)
  
  info <- .panel_model(formula, data, index, subset, contrasts)
  y <- info$y; X <- info$X; id <- info$id
  
  ## --- FOD ------------------------------------------------------------------
  tr <- .fod(cbind(y, X), id)
  yv <- tr$M[, 1L]
  Xv <- tr$M[, -1L, drop = FALSE]
  colnames(Xv) <- colnames(X)
  idv <- id[tr$keep]
  
  ## An intercept in the transformed model is a free constant of the pooled
  ## regression, which is what the original implementation offered.  It is not
  ## the intercept of the structural model: that one is time invariant and the
  ## transformation removes it along with the individual effects, exactly as
  ## footnote 2 of the paper describes.  It is collinear with a full set of
  ## time dummies, which the rank check below catches.
  if (isTRUE(intercept)) Xv <- cbind(`(Intercept)` = 1, Xv)
  
  Ti  <- table(id)
  if (any(Ti < 2L))
    message(sum(Ti < 2L), " panel(s) hold a single observation and drop out ",
            "of the estimation:\n  the transformation consumes one ",
            "observation per panel.")
  if (nrow(Xv) <= ncol(Xv))
    stop("Only ", nrow(Xv), " observations survive the transformation, for ",
         ncol(Xv), " regressors.", call. = FALSE)
  if (qr(Xv)$rank < ncol(Xv))
    stop("The transformed design matrix is rank deficient. Either a regressor ",
         "is time invariant within panels, or two of them are collinear",
         if (isTRUE(intercept))
           " -- note that an intercept is collinear with a full set of time dummies"
         else "", ".", call. = FALSE)
  
  ## ties, now on the scale the copula actually sees
  tie_msg <- character(0)
  for (j in match(info$endo_names, colnames(Xv))) {
    v  <- Xv[, j]
    tb <- tabulate(match(v, unique(v)))
    mx <- max(tb) / length(v)
    if (length(tb) <= 20L || mx > 0.05)
      tie_msg <- c(tie_msg, sprintf("%s (%d distinct values, largest tie group %.1f%%)",
                                    colnames(Xv)[j], length(tb), 100 * mx))
  }
  if (length(tie_msg) > 0L)
    warning("After the transformation the endogenous regressor(s) still carry ",
            "substantial ties: ", paste(tie_msg, collapse = "; "),
            ".\n  The copula inverts an estimated CDF, and a point mass leaves ",
            "a plateau whose inverse\n  is not unique (Qian & Xie 2024).",
            call. = FALSE)
  
  N <- length(unique(idv))
  if (N < 30L)
    warning("Only ", N, " cross-sectional units. Both the estimator and the ",
            "bootstrap rely on asymptotics in N, so the standard errors may ",
            "be imprecise.", call. = FALSE)
  
  endo <- match(info$endo_names, colnames(Xv))
  if (anyNA(endo))
    stop("Internal error locating the endogenous columns.", call. = FALSE)
  
  ## --- maximum likelihood ---------------------------------------------------
  if (verbose) message("Maximising the likelihood ...")
  cop <- .panel_copula_cols(Xv, endo, cdf, ties, intercept)
  fit <- .panel_estimate(yv, Xv, endo, cdf, ties, method, start = start,
                         maxit = maxit, cop = cop)
  if (is.null(fit))
    stop("The likelihood could not be maximised. Try method = \"Nelder-Mead\", ",
         "different starting values, or check the identification diagnostics ",
         "of an OLS fit first.", call. = FALSE)
  if (isTRUE(fit$fallback))
    message("BFGS did not converge; the reported fit comes from Nelder-Mead.")
  if (!identical(fit$convergence, 0L))
    warning("The optimiser reported convergence code ", fit$convergence,
            ". Treat the estimates with care: refit with ",
            "method = \"Nelder-Mead\", a larger 'maxit', or different ",
            "starting values and check that the answer is the same.",
            call. = FALSE)
  
  beta <- fit$beta; names(beta) <- colnames(Xv)
  rho  <- fit$rho;  names(rho)  <- colnames(Xv)[endo]
  
  ## --- bootstrap ------------------------------------------------------------
  B <- .panel_boot(yv, Xv, idv, endo, cdf, ties, method, nboots,
                   start = fit$par, cop = cop, verbose = verbose,
                   parallel = parallel, ncores = ncores)
  
  p  <- ncol(Xv); K <- length(endo)
  V  <- stats::cov(B[, seq_len(p), drop = FALSE])
  dimnames(V) <- list(names(beta), names(beta))
  se     <- sqrt(diag(V))
  rho.se <- apply(B[, p + seq_len(K), drop = FALSE], 2L, stats::sd)
  names(rho.se) <- names(rho)
  s2.se  <- stats::sd(B[, p + K + 1L])
  
  ## the uncorrected benchmark: the within estimator on the same resamples.
  ## Haschka (2022) uses the ratio of the two standard errors as the warning
  ## sign for a lack of identification.
  wb <- .lm.fit(Xv, yv)$coefficients
  names(wb) <- colnames(Xv)
  w.se <- apply(B[, p + K + 1L + seq_len(p), drop = FALSE], 2L, stats::sd)
  names(w.se) <- colnames(Xv)
  se.ratio <- se / w.se
  
  ## --- fixed effects and the two scales ------------------------------------
  ## alpha_i = mean(y_i) - mean(x_i)'beta - mean(z_i)'delta, as in the paper
  ## and as plm's fixef() does it
  ## as in lm(), keep the row names so that observations stay identifiable
  nmr <- rownames(info$mf)
  nmt <- nmr[tr$keep]
  
  bs <- beta[colnames(X)]              # structural part, intercept dropped
  fit_s_all <- as.vector(X %*% bs)
  alpha <- tapply(y - fit_s_all, id, mean)
  alpha <- stats::setNames(as.vector(alpha), names(alpha))
  aid   <- as.vector(alpha[as.character(id)])
  
  fit_t <- stats::setNames(as.vector(Xv %*% beta), nmt)   # transformed model
  res_t <- stats::setNames(as.vector(yv), nmt) - fit_t
  y     <- stats::setNames(as.vector(y), nmr)
  fit_s <- stats::setNames(aid + fit_s_all, nmr)   # structural, original rows
  res_s <- y - fit_s
  
  ## --- tests on rho ---------------------------------------------------------
  lr <- if (is.na(fit$logLik0)) NULL else {
    st <- 2 * (fit$logLik - fit$logLik0)
    c(chisq = st, df = K, p = stats::pchisq(st, K, lower.tail = FALSE))
  }
  Vr <- stats::cov(B[, p + seq_len(K), drop = FALSE])
  dimnames(Vr) <- list(names(rho), names(rho))
  wald <- .wald(rho, Vr, seq_len(K))
  
  ## --- fit statistics -------------------------------------------------------
  wthn <- function(v) v - ave(v, id, FUN = mean)
  yw  <- wthn(y); Xw <- apply(X, 2L, wthn)
  ybar <- tapply(y, id, mean)
  Xbar <- do.call(rbind, lapply(split(as.data.frame(X), id), colMeans))
  Xbar <- Xbar[, colnames(X), drop = FALSE]
  r2 <- c(
    transformed = 1 - sum(res_t^2) / sum(yv^2),
    structural  = 1 - sum(res_s^2) / sum((y - mean(y))^2),
    within      = suppressWarnings(stats::cor(yw, as.vector(Xw %*% bs))^2),
    between     = suppressWarnings(stats::cor(ybar,
                                              as.vector(Xbar %*% bs))^2),
    overall     = suppressWarnings(stats::cor(y, fit_s_all)^2))
  
  k <- p + K + 1L
  structure(list(
    call = cl, method = "Panel copula MLE (Haschka 2022)",
    cdf = cdf, ties = ties, intercept = isTRUE(intercept),
    nboots = nrow(B),
    coefficients = beta, std.error = se, vcov = V, boot = B,
    rho = rho, rho.se = rho.se, rho.vcov = Vr,
    within.coefficients = wb, within.std.error = w.se, se.ratio = se.ratio, sigma2 = fit$sigma2, sigma2.se = s2.se,
    logLik = fit$logLik, logLik.restricted = fit$logLik0, npar = k,
    AIC = -2 * fit$logLik + 2 * k,
    BIC = -2 * fit$logLik + log(nrow(Xv)) * k,
    lr.test = lr, wald.test = wald,
    convergence = fit$convergence, optim.method = fit$method,
    fallback = fit$fallback,
    fitted.structural = fit_s, residuals.structural = res_s,
    fitted.transformed = fit_t, residuals.transformed = res_t,
    fixef = alpha, r.squared = r2,
    y = y, X = X, y.transformed = yv, X.transformed = Xv,
    id = id, time = info$time, id.transformed = idv,
    N = N, Tbar = mean(Ti), Trange = range(Ti),
    n = nrow(X), n.transformed = nrow(Xv),
    endogenous = info$endo_names,
    exogenous = setdiff(colnames(X), info$endo_names),
    endo.cols = endo, index = info$index,
    terms = info$terms, Formula = info$Formula,
    xlevels = info$xlevels, contrasts = info$contrasts,
    diagnostics = local({
      Zk <- fit$Zreg[, match(endo, cop), drop = FALSE]
      colnames(Zk) <- names(rho)
      .copreg_diagnostics(yv, Xv, Zk, Zk,
                          list(endo_cols = endo,
                               exo_cols = setdiff(seq_len(p), endo),
                               endo_names = names(rho),
                               has_intercept = FALSE), res_t)
    })
  ), class = "copregpanel")
}


## =============================================================================
##  5.  S3 methods
## =============================================================================
##
##  Two scales exist side by side and every extractor takes 'structural':
##
##    structural = TRUE   the model as written down, y = alpha_i + x'beta +
##                        z'delta + e, on the original rows.  This is the
##                        default, matching the other estimators in the toolbox.
##    structural = FALSE  the FOD-transformed model that maximum likelihood
##                        actually fits, on the transformed rows.  This is the
##                        plm convention, and the residuals here are the ones
##                        the likelihood treats as N(0, sigma^2), so they are
##                        what to inspect when checking the error assumption.

coef.copregpanel  <- function(object, ...) object$coefficients
vcov.copregpanel  <- function(object, ...) object$vcov
nobs.copregpanel  <- function(object, ...) object$n
formula.copregpanel <- function(x, ...) x$Formula

logLik.copregpanel <- function(object, ...)
  structure(object$logLik, df = object$npar,
            nobs = object$n.transformed, class = "logLik")

residuals.copregpanel <- function(object, structural = TRUE, ...) {
  if (isTRUE(structural)) object$residuals.structural
  else object$residuals.transformed
}

fitted.copregpanel <- function(object, structural = TRUE, ...) {
  if (isTRUE(structural)) object$fitted.structural
  else object$fitted.transformed
}

## The individual effects, alpha_i = mean(y_i) - mean(x_i)'beta, exactly as in
## the paper and as plm::fixef computes them.
if (!exists("fixef")) fixef <- function(object, ...) UseMethod("fixef")
fixef.copregpanel <- function(object, ...) object$fixef

## Predictions use the structural model.  With structural = TRUE the panel
## effect is added, which needs the unit to have been in the estimation sample;
## unknown units get NA rather than a silently substituted average.
predict.copregpanel <- function(object, newdata = NULL, structural = TRUE, ...) {
  
  if (is.null(newdata)) {
    return(if (isTRUE(structural)) object$fitted.structural
           else object$fitted.transformed)
  }
  
  if (!isTRUE(structural))
    stop("Predictions on the transformed scale need the whole panel of a unit; ",
         "use structural = TRUE.", call. = FALSE)
  
  idv <- object$index[1L]
  if (!idv %in% names(newdata))
    stop("newdata must carry the panel identifier '", idv, "'.", call. = FALSE)
  
  Terms <- stats::delete.response(object$terms)
  environment(Terms) <- .panel_env(newdata[[idv]], environment(Terms))
  mf <- stats::model.frame(Terms, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X  <- stats::model.matrix(Terms, mf, contrasts.arg = object$contrasts)
  X  <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  b  <- object$coefficients[colnames(X)]
  if (anyNA(b))
    stop("newdata produced columns that are not in the fitted model.",
         call. = FALSE)
  
  out <- drop(X %*% b)
  a <- object$fixef[as.character(newdata[[idv]])]
  if (anyNA(a))
    warning("Panel(s) in newdata were not in the estimation sample; their ",
            "individual effect is unknown and the prediction is NA.",
            call. = FALSE)
  out + as.vector(a)
}

confint.copregpanel <- function(object, parm, level = 0.95,
                                type = c("normal", "percentile"), ...) {
  type <- match.arg(type)
  cf <- object$coefficients
  if (missing(parm)) parm <- names(cf)
  a  <- (1 - level) / 2
  if (type == "normal") {
    q  <- stats::qnorm(1 - a)
    ci <- cbind(cf - q * object$std.error, cf + q * object$std.error)
  } else {
    ci <- base::t(apply(object$boot[, seq_along(cf), drop = FALSE], 2,
                        stats::quantile, probs = c(a, 1 - a), names = FALSE))
  }
  dimnames(ci) <- list(names(cf), paste0(round(100 * c(a, 1 - a), 1), " %"))
  ci[parm, , drop = FALSE]
}

print.copregpanel <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("\n", x$method, "\n", sep = "")
  if (!is.null(x$call)) { cat("\nCall:\n"); print(x$call) }
  cat("\nBalanced panel: n = ", x$N, ", T = ", x$Trange[1],
      if (x$Trange[1] != x$Trange[2]) paste0("-", x$Trange[2]),
      ", observations = ", x$n, "\n", sep = "")
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

summary.copregpanel <- function(object, ...) {
  z <- object$coefficients / object$std.error
  tab <- cbind(Estimate = object$coefficients,
               `Std. Error` = object$std.error,
               `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  
  est <- c(object$rho, sigma2 = object$sigma2)
  ses <- c(object$rho.se, object$sigma2.se)
  names(est) <- c(paste0("rho(", names(object$rho), "*, xi*)"), "sigma2")
  zr  <- est / ses
  rtab <- cbind(Estimate = est, `Std. Error` = ses, `z value` = zr,
                `Pr(>|z|)` = 2 * stats::pnorm(-abs(zr)))
  rownames(rtab) <- names(est)
  ## sigma^2 is bounded below by zero, so a two-sided test against zero is not
  ## meaningful; only rho is tested
  rtab[nrow(rtab), 3:4] <- NA_real_
  
  structure(list(call = object$call, method = object$method,
                 cdf = object$cdf, ties = object$ties, nboots = object$nboots,
                 coefficients = tab, rho = rtab,
                 residuals.structural = object$residuals.structural,
                 residuals.transformed = object$residuals.transformed,
                 r.squared = object$r.squared,
                 logLik = object$logLik, AIC = object$AIC, BIC = object$BIC,
                 npar = object$npar,
                 lr.test = object$lr.test, wald.test = object$wald.test,
                 N = object$N, Trange = object$Trange, Tbar = object$Tbar,
                 n = object$n, n.transformed = object$n.transformed,
                 index = object$index, endogenous = object$endogenous,
                 optim.method = object$optim.method,
                 convergence = object$convergence, fallback = object$fallback,
                 intercept = object$intercept,
                 diagnostics = object$diagnostics),
            class = "summary.copregpanel")
}

print.summary.copregpanel <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      signif.stars = getOption("show.signif.stars"),
                                      ...) {
  
  cat("\n", x$method, "\n", sep = "")
  if (!is.null(x$call)) { cat("\nCall:\n"); print(x$call) }
  
  cat("\nPanel: index = (", paste(x$index, collapse = ", "), ")\n", sep = "")
  cat("  ", x$N, " cross-sectional units, T = ", x$Trange[1],
      if (x$Trange[1] != x$Trange[2]) paste0(" to ", x$Trange[2]) else "",
      if (x$Trange[1] != x$Trange[2])
        paste0(" (mean ", formatC(x$Tbar, digits = 3), ")") else "",
      "\n", sep = "")
  cat("  ", x$n, " observations, ", x$n.transformed,
      " after the forward orthogonal deviations transformation\n", sep = "")
  cat(if (isTRUE(x$intercept))
    "  the transformed regression carries a free constant\n"
    else paste0("  no constant in the transformed regression: the structural",
                " intercept is\n  time invariant and goes with the",
                " individual effects\n"))
  
  qn <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("\nResiduals of the structural model (y - alpha_i - x'beta - z'delta):\n")
  print(structure(stats::quantile(x$residuals.structural), names = qn),
        digits = digits)
  cat("\nResiduals of the transformed model, which the likelihood treats as",
      "normal:\n")
  print(structure(stats::quantile(x$residuals.transformed), names = qn),
        digits = digits)
  
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, has.Pvalue = TRUE,
                      P.values = TRUE, na.print = "NA")
  
  cat("\nDependence parameters: rho(P*, xi*) is the correlation between the",
      "normal\n  score of an endogenous regressor and that of the error, and",
      "sigma^2 is the\n  variance of the error of the transformed model.",
      "rho = 0 means no endogeneity.\n")
  stats::printCoefmat(x$rho, digits = digits, signif.stars = signif.stars,
                      has.Pvalue = TRUE, P.values = TRUE, na.print = "")
  
  if (!is.null(x$lr.test))
    cat("\nNo endogeneity, all rho = 0:\n  likelihood ratio  chi-squared = ",
        formatC(x$lr.test["chisq"], digits = digits), " on ", x$lr.test["df"],
        " df, p = ", format.pval(x$lr.test["p"], digits = digits), "\n",
        sep = "")
  if (!is.null(x$wald.test))
    cat("  bootstrap Wald   chi-squared = ",
        formatC(x$wald.test["chisq"], digits = digits), " on ",
        x$wald.test["df"], " df, p = ",
        format.pval(x$wald.test["p"], digits = digits), "\n", sep = "")
  
  r2 <- x$r.squared
  cat("\nR-squared:\n")
  cat("  transformed model ", formatC(r2["transformed"], digits = digits),
      "   structural model ", formatC(r2["structural"], digits = digits),
      "\n", sep = "")
  cat("  within ", formatC(r2["within"], digits = digits),
      "   between ", formatC(r2["between"], digits = digits),
      "   overall ", formatC(r2["overall"], digits = digits), "\n", sep = "")
  cat("  within, between and overall are squared correlations excluding the\n",
      "  individual effects.\n", sep = "")
  
  cat("\nLog-likelihood ", formatC(x$logLik, digits = digits),
      " on ", x$npar, " parameters;  AIC ", formatC(x$AIC, digits = digits),
      ",  BIC ", formatC(x$BIC, digits = digits), "\n", sep = "")
  cat("Standard errors from ", x$nboots, " panel bootstrap replicates ",
      "(cross-sectional units\n  resampled, not rows); cdf = \"", x$cdf,
      "\", ties = \"", x$ties, "\"; optimiser ", x$optim.method,
      if (isTRUE(x$fallback)) " (BFGS fell back)" else "", ".\n", sep = "")
  cat("Pr(>|z|) in both tables: Wald test using the normal approximation with\n",
      "  the bootstrap standard error.\n", sep = "")
  
  d <- x$diagnostics
  cat("\n--- Identification diagnostics ------------------------------------\n")
  cat("\nNon-normality of the transformed endogenous regressors\n",
      "(small p = non-normal, which is what identifies the model):\n", sep = "")
  print(format(d$nonnormality, digits = digits))
  cat("\nCollinearity of the copula data (omega near 0 = weakly identified):\n")
  print(format(d$collinearity, digits = digits))
  cat("\n")
  invisible(x)
}


### REFERENCES ----------------------------------------------------------------
##
## Arellano, M. (1993). On the testing of correlated effects with panel data.
##   Journal of Econometrics 59, 87-97.
##
## Goncalves, S. and L. Kilian (2004). Bootstrapping autoregressions with
##   conditional heteroskedasticity of unknown form. Journal of Econometrics
##   123, 89-120.
##
## Haschka, R. E. (2022). Handling endogenous regressors using copulas: A
##   generalization to linear panel models with fixed effects and correlated
##   regressors. Journal of Marketing Research 59(4), 861-880.
##
## Ngom, P., B. Moussa, and J. De Dieu Nkurunziza (2018). A new nonparametric
##   density estimator.
##
## Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##   estimation using copulas. Marketing Science 31(4), 567-586.
##
## Polansky, A. M. and E. R. Baker (2000). Multistage plug-in bandwidth
##   selection for kernel distribution function estimates. Journal of
##   Statistical Computation and Simulation 65, 63-80.


## =============================================================================
##  6.  Validity check
## =============================================================================
##
##  Three of the checkpoints used for the cross-sectional estimators carry over,
##  and one does not.
##
##  [1] Nonnormality of the endogenous regressors, read off the transformed
##      variables because that is what the copula sees.  Becker, Proksch &
##      Ringle (2022) study fixed-effects panel models in their Study 3 and
##      conclude that "in respect of multilevel (or panel) models, the same
##      recommendations apply regarding a sufficient sample size and
##      nonnormality as do for the simpler, cross-sectional regression model",
##      once the total number of observations is the reference.  Their
##      thresholds are therefore applied to n, not to the number of panels.
##
##  [2] The error term.  Here this is not a matter of interpretation: the
##      likelihood in Equation 12 models the error of the transformed model as
##      normal, so the residuals of that model are the direct check on a stated
##      assumption.
##
##  [3] Standard error inflation relative to the uncorrected within estimator,
##      computed on the same bootstrap resamples.  Haschka (2022) gives the
##      threshold himself: standard deviations "about five to ten times higher
##      compared with those of OLS clearly indicate a lack of identification
##      and may serve as a warning in empirical applications".
##
##  What does not carry over is the assumption that the copula data are
##  uncorrelated with the exogenous regressors.  Modelling that correlation
##  rather than assuming it away is the point of this estimator, so there is
##  nothing to test.

validity.copregpanel <- function(object, level = 0.05, power = 0.8, ...) {
  
  Xv <- object$X.transformed
  nm <- object$endogenous
  n  <- object$n
  
  nn <- .validity_nonnormality(Xv, nm, n, power)
  
  xi <- object$residuals.transformed
  a4 <- .ad_test(xi)
  err <- list(skewness = .skewness(xi), ex.kurtosis = .ex_kurtosis(xi),
              AD = unname(a4["A"]), AD.p = unname(a4["p"]),
              KS.p = unname(.ks_normal(xi)["p"]))
  
  infl <- data.frame(
    `SE (copula MLE)` = object$std.error,
    `SE (within)`     = object$within.std.error,
    ratio             = object$se.ratio,
    row.names = names(object$std.error), check.names = FALSE)
  
  structure(list(method = object$method, n = n, N = object$N,
                 Trange = object$Trange, level = level, power = power,
                 endogenous = nm, step1 = nn$table,
                 thresholds = nn$thresholds, error = err,
                 inflation = infl,
                 ratio.max = max(object$se.ratio, na.rm = TRUE)),
            class = "copregpanel.validity")
}


print.copregpanel.validity <- function(x, digits = 4, ...) {
  
  cat("\nValidity check for ", x$method, "\n", sep = "")
  cat("n = ", x$n, " observations in ", x$N, " panels, target power ",
      100 * x$power, "%\n", sep = "")
  cat("Sources: Becker, Proksch & Ringle (2022); Haschka (2022)\n")
  
  cat("\n[1] Nonnormality of the endogenous regressors, after the ",
      "transformation\n", sep = "")
  print(format(x$step1, digits = digits))
  th <- x$thresholds
  cat("    Becker et al. at n = ", x$n, ": |skewness| >= ",
      if (is.infinite(th["skewness"])) "not attainable" else
        formatC(th["skewness"], digits = 4),
      ", or AD > ", formatC(th["AD"], digits = 5),
      ", or CvM > ", formatC(th["CvM"], digits = 4), "\n", sep = "")
  cat("    Their Study 3 covers fixed-effects panels and finds the ",
      "cross-sectional\n    thresholds carry over once total n is the ",
      "reference.\n", sep = "")
  
  e <- x$error
  cat("\n[2] Error of the transformed model, which Equation 12 assumes normal\n")
  cat("    skewness = ", formatC(e$skewness, digits = digits),
      ", excess kurtosis = ", formatC(e$ex.kurtosis, digits = digits),
      ", AD = ", formatC(e$AD, digits = digits),
      " (p = ", format.pval(e$AD.p, digits = 3), ")\n", sep = "")
  cat("    Unlike in the cross-sectional estimators this is a stated ",
      "assumption of the\n    likelihood rather than a matter of ",
      "interpretation, so a small p is a real\n    warning sign.\n", sep = "")
  
  cat("\n[3] Standard errors against the uncorrected within estimator, on the\n",
      "    same bootstrap resamples\n", sep = "")
  print(format(x$inflation, digits = digits))
  cat("    => largest ratio = ", formatC(x$ratio.max, digits = digits),
      if (x$ratio.max > 5)
        ".\n       Haschka (2022) reads a factor of five to ten as a clear sign that the\n       model is not identified."
      else ".\n       Haschka (2022) reads a factor of five to ten as a sign that the model\n       is not identified; this is below that.",
      "\n\n", sep = "")
  
  invisible(x)
}
