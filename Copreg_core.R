## =============================================================================
##  copreg-core.R
##
##  Shared infrastructure for copula-based endogeneity corrections.
##  Source this file before any of the estimator files (CopRegPG.R, ...).
##
##  Depends on: Formula (CRAN).  Optionally np, for cdf = "kde.cv".
##  Everything else is base R.
## =============================================================================

if (!requireNamespace("Formula", quietly = TRUE)) {
  stop("Package 'Formula' is required. Install it with install.packages('Formula').")
}


## -----------------------------------------------------------------------------
##  Copula-based endogeneity corrections in R
##  https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  Copyright (C) 2026 Rouven E. Haschka
##  ORCID: https://orcid.org/0000-0002-2916-9745
##
##  This program is free software: you can redistribute it and/or modify it
##  under the terms of the GNU General Public License as published by the Free
##  Software Foundation, either version 3 of the License, or (at your option)
##  any later version.
##
##  Additional term under section 7(b) of that licence: you must preserve the
##  author attribution above in any material you convey, and in the notices
##  displayed by works containing it. This applies to modified versions and to
##  larger works that incorporate this material, other R packages included.
##
##  This program is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
##  more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program.  If not, see <https://www.gnu.org/licenses/>.
##
##  Citation, a request rather than a licence condition. If this code
##  contributes to work you publish, please cite the software
##
##    Haschka, R. E. (2026). Copula-based endogeneity corrections in R.
##    https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  together with the papers named in the estimator files, for whichever
##  estimator you used.
## -----------------------------------------------------------------------------


## =============================================================================
##  1.  CDF estimators
## =============================================================================
##
##  Every variant below is taken literally from the paper that proposes it.
##  See ?copreg.cdf (bottom of this file) for the mapping to the literature.

.copreg_cdf_choices  <- c("kde.silverman", "kde.cv", "kde.plugin",
                          "ecdf.fixed", "ecdf.adj",
                          "rank.n", "rank.n1")

.copreg_ties_choices <- c("max", "average")


## Row-wise key for a set of columns, used wherever observations have to be
## grouped into the cells of several discrete variables.  as.data.frame() would
## be the obvious way to feed paste(), but it runs make.names() and
## make.unique() on every call, and this runs once per bootstrap replicate.
.copreg_key <- function(M) {
  if (is.data.frame(M)) M <- as.matrix(M)
  if (ncol(M) == 1L) return(as.character(M[, 1L]))
  do.call(paste, c(lapply(seq_len(ncol(M)), function(j) M[, j]),
                   list(sep = "\r")))
}


## Seconds as a readable duration, for the runtime estimates that the slower
## estimators print before they start.
.copreg_hms <- function(s) {
  if (!is.finite(s))  return("an unknown time")
  if (s < 90)         return(sprintf("%.0f seconds", s))
  if (s < 5400)       return(sprintf("%.0f minutes", s / 60))
  sprintf("%.1f hours", s / 3600)
}


## Counting function  sum_i I(x_i <= x_j).
##
## ties = "max"      reproduces the counting function literally, i.e. the way
##                   every paper writes F(x) = (1/n) sum I(P_i <= x).
## ties = "average"  uses midranks, the convention of the wider copula
##                   literature (Genest, Ghoudi & Rivest 1995) and of
##                   copula::pobs().  This is a deviation from the papers and
##                   only matters for variables with tied values.
.counts_le <- function(x, ties = "max") {
  rank(x, ties.method = if (identical(ties, "max")) "max" else "average")
}


## CDF of the Epanechnikov kernel K(t) = 0.75 (1 - t^2) I(|t| <= 1).
.pepan <- function(u) {
  u <- pmin(pmax(u, -1), 1)
  0.75 * (u - u^3 / 3) + 0.5
}


## Silverman's rule as used by Park & Gupta (2012, p. 571):
##   b = 0.9 * n^(-1/5) * min(s, IQR / 1.34)
.bw_silverman <- function(x) {
  n  <- length(x)
  s  <- stats::sd(x)
  iq <- stats::IQR(x) / 1.34
  sc <- min(s, iq)
  if (!is.finite(sc) || sc <= 0) sc <- if (is.finite(s) && s > 0) s else 1
  0.9 * n^(-1 / 5) * sc
}


## Kernel CDF  F(x_j) = (1/n) sum_i G((x_j - x_i)/b)  with G = .pepan.
##
## Park & Gupta integrate the density numerically; because the Epanechnikov
## kernel integrates in closed form, we use the exact antiderivative instead.
## This is the same estimator, without quadrature error.
##
## Implementation is O(n log n): G is a cubic on its support, so the window
## sums are obtained from prefix sums of x, x^2, x^3.  x is centred and scaled
## first to keep the cubic prefix sums well conditioned.
.cdf_kde_epan <- function(x, b) {
  n <- length(x)
  if (!is.finite(b) || b <= 0)
    stop("Non-positive bandwidth in kernel CDF estimation.", call. = FALSE)
  
  ## affine equivariance: F_x(x_j; b) = F_z(z_j; b/s) for z = (x - m)/s
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  z  <- (x - m) / s
  bb <- b / s
  
  o  <- order(z)
  zs <- z[o]
  
  ## window (lo+1):hi contains all points within +/- bb of zs[j]
  lo <- findInterval(zs - bb, zs)
  hi <- findInterval(zs + bb, zs)
  
  S1 <- c(0, cumsum(zs))
  S2 <- c(0, cumsum(zs^2))
  S3 <- c(0, cumsum(zs^3))
  
  m0 <- hi - lo
  m1 <- S1[hi + 1L] - S1[lo + 1L]
  m2 <- S2[hi + 1L] - S2[lo + 1L]
  m3 <- S3[hi + 1L] - S3[lo + 1L]
  
  ## sum over window of d and d^3, where d = (zs[j] - zs[i]) / bb
  sd1 <- (m0 * zs - m1) / bb
  sd3 <- (m0 * zs^3 - 3 * zs^2 * m1 + 3 * zs * m2 - m3) / bb^3
  
  out <- (lo + 0.75 * sd1 - 0.25 * sd3 + 0.5 * m0) / n
  
  res <- numeric(n)
  res[o] <- out
  pmin(pmax(res, 0), 1)
}


## Plug-in bandwidth for the distribution function (Polansky & Baker 2000).
## This is what ks::hpi.kcde computes, reimplemented here so that the estimator
## carries no extra dependency; it reproduces ks::hpi.kcde(x, binned = FALSE) to
## machine precision.
##
##   psi_r      = int f^(r) f, estimated by n^-2 sum_i sum_j phi^(r)(X_i - X_j; g)
##   two stages: normal scale psi_6 -> g_4 -> psi_4 -> g_2 -> psi_2 -> h

.He <- function(z, r) switch(as.character(r),
                             "2" = z^2 - 1,
                             "4" = z^4 - 6 * z^2 + 3,
                             "6" = z^6 - 15 * z^4 + 45 * z^2 - 15,
                             stop("Hermite polynomial of order ", r, " not needed here.", call. = FALSE))

## Weighted bin counts on 1..M.  tapply() would do the same but builds a factor
## every time, and this runs inside every bootstrap replicate for every
## regressor, so it is worth the few extra lines.
.wbin <- function(idx, w, M) {
  o  <- order(idx)
  i  <- idx[o]
  cs <- cumsum(w[o])
  last <- c(i[-1L] != i[-length(i)], TRUE)
  out <- numeric(M)
  out[i[last]] <- diff(c(0, cs[last]))
  out
}

## r-th derivative of the N(0, sigma^2) density
.dnorm_deriv <- function(x, sigma, r)
  (-1)^r * .He(x / sigma, r) * stats::dnorm(x / sigma) / sigma^(r + 1)

## normal-scale reference value of psi_r
.psins <- function(r, sigma)
  (-1)^(r / 2) * factorial(r) /
  ((2 * sigma)^(r + 1) * factorial(r / 2) * sqrt(pi))

## kernel functional estimate.  Exact for moderate n; above that the double sum
## is replaced by a linear binning plus an FFT autocorrelation, which is what
## ks does internally and costs O(M log M) instead of O(n^2).  With M = 4096 the
## two agree to about 1e-4 relative, far below anything that matters for a
## bandwidth.
.kfe <- function(x, g, r, exact_max = 1500L, M = 4096L) {
  
  n <- length(x)
  
  if (n <= exact_max) {
    s <- 0
    for (a in seq(1L, n, by = 1000L)) {
      b <- min(a + 999L, n)
      s <- s + sum(.dnorm_deriv(outer(x[a:b], x, "-"), g, r))
    }
    return(s / n^2)
  }
  
  lo0 <- min(x); hi0 <- max(x)
  if (hi0 <= lo0) return(.dnorm_deriv(0, g, r))
  delta <- (hi0 - lo0) / (M - 1L)
  
  idx <- (x - lo0) / delta
  lo  <- pmin(pmax(floor(idx) + 1L, 1L), M)
  w   <- idx - (lo - 1L)
  hi  <- pmin(lo + 1L, M)
  cnt <- .wbin(c(lo, hi), c(1 - w, w), M)
  
  L  <- 2^ceiling(log2(2 * M))
  f  <- stats::fft(c(cnt, numeric(L - M)))
  ac <- Re(stats::fft(f * Conj(f), inverse = TRUE)) / L
  sl <- ac[1:M]
  
  ## phi^(r) is an even function for even r, so lags +l and -l contribute alike
  k <- .dnorm_deriv((0:(M - 1L)) * delta, g, r)
  (k[1] * sl[1] + 2 * sum(k[-1] * sl[-1])) / n^2
}

.bw_hpi <- function(x) {
  n  <- length(x)
  K2 <- .dnorm_deriv(0, 1, 2)
  K4 <- .dnorm_deriv(0, 1, 4)
  m1 <- (4 * pi)^(-1 / 2)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx <= 0)
    stop("Cannot choose a bandwidth for a constant variable.", call. = FALSE)
  psi6 <- .psins(6, sx)
  g4   <- (2 * K4 / (-psi6 * n))^(1 / 7)
  psi4 <- .kfe(x, g4, 4)
  g2   <- (2 * K2 / (-psi4 * n))^(1 / 5)
  psi2 <- .kfe(x, g2, 2)
  (2 * m1 / (-psi2 * n))^(1 / 3)
}

## Gaussian kernel CDF, F(x_j) = n^-1 sum_i Phi((x_j - X_i) / h).  Exact for
## moderate n, otherwise binned on a grid and interpolated back, the way
## ks::kcde and its predict method work.
.cdf_kde_gauss <- function(x, h, exact_max = 1500L, M = 4096L) {
  
  n <- length(x)
  if (!is.finite(h) || h <= 0)
    stop("Non-positive bandwidth in kernel CDF estimation.", call. = FALSE)
  
  if (n <= exact_max) {
    out <- numeric(n)
    for (a in seq(1L, n, by = 1000L)) {
      b <- min(a + 999L, n)
      out[a:b] <- rowMeans(stats::pnorm(outer(x[a:b], x, "-") / h))
    }
    return(out)
  }
  
  lo0 <- min(x) - 8 * h; hi0 <- max(x) + 8 * h
  gr  <- seq(lo0, hi0, length.out = M)
  delta <- gr[2] - gr[1]
  
  idx <- (x - lo0) / delta
  lo  <- pmin(pmax(floor(idx) + 1L, 1L), M)
  w   <- idx - (lo - 1L)
  hi  <- pmin(lo + 1L, M)
  cnt <- .wbin(c(lo, hi), c(1 - w, w), M)
  
  ## linear convolution of the counts with Phi(. / h) via zero padding
  L   <- 2^ceiling(log2(2 * M))
  ker <- stats::pnorm((seq_len(L) - 1L - (L / 2)) * delta / h)
  cv  <- Re(stats::fft(stats::fft(c(cnt, numeric(L - M))) *
                         stats::fft(ker), inverse = TRUE)) / L
  Fg  <- cv[(L / 2) + seq_len(M)] / n
  
  Fg <- cummax(pmin(pmax(Fg, 0), 1))
  stats::approx(gr, Fg, xout = x, rule = 2)$y
}


## Cross-validated bandwidth for the distribution function (Li, Li & Racine
## 2017), as referenced by Liengaard et al. (2025, Table 1, estimator F2).
## Delegated to the np package, which is the implementation they name.
.cdf_kde_cv <- function(x) {
  if (!requireNamespace("np", quietly = TRUE))
    stop("cdf = \"kde.cv\" requires the 'np' package. ",
         "Install it with install.packages('np').", call. = FALSE)
  bw <- np::npudistbw(dat = x, ckertype = "epanechnikov")
  as.numeric(np::npudist(bws = bw, edat = x)$dist)
}


## Dispatcher.  Returns u in (0, 1).
.cdf_estimate <- function(x, cdf, ties) {
  n <- length(x)
  switch(cdf,
         
         ## Yang, Qian & Xie (2025); Qian, Koschmann & Xie (2025), Eq. 9
         "rank.n" = {
           u <- .counts_le(x, ties) / n
           u[x == max(x)] <- n / (n + 1)
           u
         },
         
         ## Breitung, Mayer & Wied (2024), Eq. 2.3
         "rank.n1" = .counts_le(x, ties) / (n + 1),
         
         ## Becker, Proksch & Ringle (2022) -- estimator F3 in Liengaard et al.
         "ecdf.fixed" = {
           u <- .counts_le(x, ties) / n
           u[u <= 0] <- 1e-7
           u[u >= 1] <- 1 - 1e-7
           u
         },
         
         ## Liengaard et al. (2025), Eq. 9 -- estimator F4
         "ecdf.adj" = 1 / (2 * n) + (n - 1) / n^2 * .counts_le(x, ties),
         
         ## Park & Gupta (2012), Eq. 3 -- estimator F1
         "kde.silverman" = .cdf_kde_epan(x, .bw_silverman(x)),
         
         ## Li, Li & Racine (2017) -- estimator F2
         "kde.cv" = .cdf_kde_cv(x),
         
         ## Gaussian kernel with the Polansky & Baker (2000) plug-in bandwidth.
         ## Equivalent to ks::kcde, the estimator used in the original panel code.
         "kde.plugin" = .cdf_kde_gauss(x, .bw_hpi(x)),
         
         stop("Unknown cdf: ", cdf, call. = FALSE)
  )
}


## Copula transformation  C(x) = Phi^{-1}(Fhat(x)).
.copula_transform <- function(x, cdf, ties, varname = "regressor") {
  u <- .cdf_estimate(x, cdf, ties)
  cc <- stats::qnorm(u)
  if (anyNA(cc) || any(!is.finite(cc))) {
    stop("Copula transformation of '", varname, "' produced non-finite values ",
         "(the estimated CDF hit 0 or 1). Try a different 'cdf'.", call. = FALSE)
  }
  cc
}


## Column-wise version.
.copula_transform_matrix <- function(M, cdf, ties) {
  out <- matrix(NA_real_, nrow(M), ncol(M), dimnames = dimnames(M))
  for (j in seq_len(ncol(M)))
    out[, j] <- .copula_transform(M[, j], cdf, ties, colnames(M)[j])
  out
}


## =============================================================================
##  2.  Formula handling
## =============================================================================
##
##  Syntax:  y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##
##  Anything lm() accepts is accepted here: transformations (log(p), I(p^2)),
##  interactions (p*m), factors (as.factor(g), g), "-1" to drop the intercept.
##  The only restriction is that endogenous regressors must be numeric.
##
##  Which side of the "|" a term sits on is the whole decision.  Every term in
##  the endogenous part gets a copula term, interactions and polynomials
##  included; every term in the exogenous part gets none.  So
##
##      y ~ p + p:m | m + w      p:m is endogenous and gets a copula term
##      y ~ p | p:m + m + w      p:m is exogenous and gets none
##
##  Qian, Koschmann & Xie (2025) report that leaving higher-order terms out of
##  the copula performs better, because such a term is strongly collinear with
##  the regressor it is built from and the added variance outweighs anything
##  gained.  That is a recommendation, not a rule, so the estimator follows the
##  formula and warns when two endogenous terms share a variable.

.copreg_model <- function(formula, data, subset = NULL, contrasts = NULL) {
  
  F <- Formula::as.Formula(formula)
  nparts <- length(F)
  
  if (nparts[1L] != 1L)
    stop("The formula must have exactly one dependent variable.", call. = FALSE)
  if (nparts[2L] > 3L)
    stop("The formula must have at most three parts: ",
         "endogenous | exogenous | discrete.", call. = FALSE)
  
  has_exog <- nparts[2L] >= 2L
  rhs_all  <- seq_len(nparts[2L])
  
  if (!is.null(subset)) data <- data[subset, , drop = FALSE]
  if (!is.data.frame(data)) data <- as.data.frame(data)
  
  ## --- combined model frame -------------------------------------------------
  tt <- stats::terms(F, rhs = rhs_all, data = data)
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.omit)
  mt <- attr(mf, "terms")
  
  y <- stats::model.response(mf)
  if (!is.numeric(y))
    stop("The dependent variable must be numeric.", call. = FALSE)
  
  X <- stats::model.matrix(mt, mf, contrasts.arg = contrasts)
  
  n_dropped <- nrow(data) - nrow(mf)
  if (n_dropped > 0L)
    message(n_dropped, " observation(s) removed because of missing values.")
  
  if (nrow(X) <= ncol(X))
    stop("Not enough complete observations: ", nrow(X), " rows for ",
         ncol(X), " regressors.", call. = FALSE)
  
  ## --- which columns are endogenous first-order terms? ----------------------
  t_endo    <- stats::terms(F, rhs = 1, lhs = 0, data = data)
  endo_lab  <- attr(t_endo, "term.labels")
  endo_ord  <- attr(t_endo, "order")
  endo_main <- endo_lab
  
  if (length(endo_main) == 0L)
    stop("No endogenous regressor found in the first part of the formula.",
         call. = FALSE)
  
  asg <- attr(X, "assign")
  lab <- character(length(asg))
  lab[asg > 0L] <- attr(mt, "term.labels")[asg[asg > 0L]]
  
  ## An endogenous term must map to exactly one model matrix column carrying
  ## its own name.  Factors expand into contrast columns and fail this test.
  bad <- character(0)
  endo_cols <- integer(0)
  for (l in endo_main) {
    j <- which(lab == l)
    if (length(j) != 1L || !identical(colnames(X)[j], l)) {
      bad <- c(bad, l)
    } else {
      endo_cols <- c(endo_cols, j)
    }
  }
  if (length(bad) > 0L)
    stop("Every term before the \"|\" needs its own copula term, so it has to ",
         "map to a single\n  numeric column. These do not: ",
         paste(bad, collapse = ", "),
         ".\n  A factor, or an interaction involving one, expands into several ",
         "columns and has\n  no single copula term. Move it behind the \"|\" to ",
         "treat it as exogenous.", call. = FALSE)
  
  exo_cols <- setdiff(seq_len(ncol(X)), endo_cols)
  
  for (j in endo_cols) {
    if (length(unique(X[, j])) < 2L)
      stop("Endogenous regressor '", colnames(X)[j], "' is constant.",
           call. = FALSE)
  }
  
  ## --- which formula part does each column come from? ----------------------
  ## part 1 endogenous, part 2 exogenous, part 3 discrete (used by JAMS).
  part <- integer(ncol(X))
  for (q in rhs_all) {
    lq <- attr(stats::terms(F, rhs = q, lhs = 0, data = data), "term.labels")
    part[lab %in% lq] <- q
  }
  
  ## interaction order of the term behind each column, 0 for the intercept.
  ## Only first-order terms are regressors in the sense of P and W in the
  ## papers; interactions belong to g() but must not enter the copula
  ## construction, since an interaction with an endogenous regressor is not
  ## exogenous information.
  ord <- integer(ncol(X))
  ord[asg > 0L] <- attr(mt, "order")[asg[asg > 0L]]
  
  ## columns produced by a factor or character variable, for estimators that
  ## have to tell discrete regressors from continuous ones
  ## A term counts as factor derived if its own label is a factor column of the
  ## model frame -- which is how as.factor(x) shows up, the model frame naming
  ## that column "as.factor(x)" rather than "x" -- or if any variable inside it
  ## is one, which covers interactions between factors.
  fvars <- names(mf)[vapply(mf, function(v) is.factor(v) || is.character(v),
                            logical(1))]
  alab <- attr(mt, "term.labels")
  flab <- alab[vapply(alab, function(l)
    l %in% fvars || any(all.vars(stats::reformulate(l)) %in% fvars),
    logical(1))]
  factor_cols <- which(lab %in% flab)
  
  ## --- endogenous terms built from the same variable ------------------------
  ## A copula term for p and one for p:m or I(p^2) are strongly collinear with
  ## each other and with the regressors they come from.  The formula decides,
  ## but the cost is worth stating.
  evars <- lapply(endo_main, function(l) all.vars(stats::reformulate(l)))
  shared <- vapply(seq_along(endo_main), function(i)
    any(vapply(seq_along(endo_main)[-i], function(j)
      length(intersect(evars[[i]], evars[[j]])) > 0L, logical(1))), logical(1))
  
  if (any(shared))
    warning("These endogenous terms are built from the same variable(s) and ",
            "each gets its own\n  copula term: ",
            paste(endo_main[shared], collapse = ", "),
            ".\n  Such terms are strongly collinear with one another, which ",
            "Qian, Koschmann &\n  Xie (2025) show inflates the variance of the ",
            "structural estimates. Moving the\n  higher-order term behind the ",
            "\"|\" treats it as exogenous and leaves it without\n  a copula term.",
            call. = FALSE)
  
  ## --- ties in endogenous regressors ---------------------------------------
  ## What breaks the copula transformation is not the number of distinct values
  ## as such but a plateau in the CDF: a point mass makes the inverse mapping
  ## non-unique (Qian & Xie 2024).  Flag a regressor if it has few distinct
  ## values in absolute terms, or if any single value carries more than 5% of
  ## the sample.
  tie_msg <- character(0)
  for (j in endo_cols) {
    v  <- X[, j]
    tb <- tabulate(match(v, unique(v)))
    nu <- length(tb)
    mx <- max(tb) / length(v)
    if (nu <= 20L || mx > 0.05)
      tie_msg <- c(tie_msg, sprintf("%s (%d distinct values, largest tie group %.1f%%)",
                                    colnames(X)[j], nu, 100 * mx))
  }
  if (length(tie_msg) > 0L)
    warning("Endogenous regressor(s) with substantial ties: ",
            paste(tie_msg, collapse = "; "),
            ".\n  Copula control functions invert an estimated CDF; a point ",
            "mass creates a plateau whose\n  inverse is not unique (Qian & Xie ",
            "2024).", call. = FALSE)
  
  
  list(y = y, X = X, mf = mf, terms = mt, Formula = F,
       endo_cols = endo_cols, exo_cols = exo_cols,
       part = part, ord = ord, nparts = nparts[2L], factor_cols = factor_cols,
       endo_names = colnames(X)[endo_cols],
       has_intercept = attr(mt, "intercept") == 1L,
       has_exog = has_exog,
       xlevels = stats::.getXlevels(mt, mf),
       contrasts = attr(X, "contrasts"),
       na.action = attr(mf, "na.action"))
}


## =============================================================================
##  3.  Bootstrap
## =============================================================================
##
##  Plain pairs bootstrap: draw row indices with replacement, recompute the
##  copula terms on the resampled data, refit.  This is the procedure used in
##  every one of the underlying papers (e.g. Hu, Qian & Xie 2025, Sec. 2.5).
##
##  No seed is set anywhere.  If you need reproducibility, call set.seed()
##  yourself before calling the estimator.
##
##  Resamples in which a factor level disappears (or the design otherwise
##  loses rank) are rejected and redrawn.

## Which columns can realistically degenerate in a resample?  Only dummies and
## variables with few distinct values; a continuous column drawn n times from n
## distinct values is never constant.  Checking those explicitly is what made
## the bootstrap slow, so we classify once, up front.
.copreg_check_levels <- function(X, has_intercept) {
  
  cols <- setdiff(seq_len(ncol(X)), if (has_intercept) 1L else integer(0))
  n <- nrow(X)
  
  dummy <- integer(0); discrete <- integer(0); thin <- character(0)
  for (j in cols) {
    v <- X[, j]
    m <- NA_integer_
    if (all(v == 0 | v == 1)) {
      dummy <- c(dummy, j)
      m <- min(sum(v), n - sum(v))
    } else if (length(unique(v)) <= 20L) {
      discrete <- c(discrete, j)
      m <- min(table(v))
    }
    ## A resample draws n rows with replacement, so a group of size m survives
    ## with probability 1 - (1 - m/n)^n.  For m = 1 that is only 63%, for m = 3
    ## it is 95%, and every draw that loses a level has to be thrown away.
    if (!is.na(m)) {
      pl <- (1 - m / n)^n
      if (pl > 0.01)
        thin <- c(thin, sprintf("%s: smallest group %d, lost in %.0f%% of draws",
                                colnames(X)[j], m, 100 * pl))
    }
  }
  
  if (length(thin) > 0L)
    warning("Thinly populated column(s). Resamples that lose a level cannot be ",
            "used and are\n  redrawn, so the bootstrap gets slow and can fail ",
            "outright:\n  ", paste(thin, collapse = "\n  "), call. = FALSE)
  
  list(dummy = dummy, discrete = discrete)
}

## What a PSOCK worker needs in its own global environment.
##
## A hand-kept list is the wrong tool here: forgetting one name does not fail
## loudly, it makes every single bootstrap replicate throw "could not find
## function", which tryCatch swallows, and the user is left with the useless
## message that no replicate could be estimated.  So the helpers are collected
## by name instead.  When these files are sourced they live in the global
## environment; inside a package they live in the namespace, and
## environment(.copreg_boot) resolves to the right one either way.
.copreg_export_names <- function(envir) {
  nm <- ls(envir, all.names = TRUE)
  nm <- nm[startsWith(nm, ".")]
  nm[vapply(nm, function(n) {
    o <- get(n, envir = envir)
    is.function(o) || is.character(o)
  }, logical(1))]
}

.copreg_boot <- function(y, X, info, ctor, cdf, ties, nboots,
                         verbose = interactive(), parallel = FALSE,
                         ncores = NULL, max_tries = 200L) {
  
  n <- nrow(X)
  chk <- .copreg_check_levels(X, info$has_intercept)
  dcols <- chk$dummy; scols <- chk$discrete
  
  ## one replicate; returns c(coefficients, rho, ols, tries) or NULL if no
  ## usable draw was found.  The attempt count travels with the result so that
  ## the caller can report how much of the resampling was thrown away.
  one <- function(...) {
    for (a in seq_len(max_tries)) {
      idx <- sample.int(n, n, replace = TRUE)
      Xb  <- X[idx, , drop = FALSE]
      ## Which rows were drawn.  Subsetting a model matrix drops every
      ## attribute but dim and dimnames, so nothing is overwritten here.  A
      ## constructor that needs the model frame (2sCOPE-np conditions on
      ## factors, which the design matrix has already expanded into dummies)
      ## reads this; the others never look at it.
      attr(Xb, "idx") <- idx
      
      ok <- TRUE
      for (j in dcols) {
        s <- sum(Xb[, j])
        if (s == 0 || s == n) { ok <- FALSE; break }
      }
      if (ok) for (j in scols) {
        if (length(unique(Xb[, j])) < 2L) { ok <- FALSE; break }
      }
      if (!ok) next
      
      cb <- tryCatch(ctor(Xb, info, cdf, ties), error = function(e) NULL)
      if (is.null(cb)) next
      
      Ab <- cbind(Xb, cb$C)
      fb <- tryCatch(.lm.fit(Ab, y[idx]), error = function(e) NULL)
      if (is.null(fb) || fb$rank < ncol(Ab)) next
      
      ## the same resample without the copula terms, for the ICON statistic
      fo <- tryCatch(.lm.fit(Xb, y[idx]), error = function(e) NULL)
      if (is.null(fo) || fo$rank < ncol(Xb)) next
      
      cf   <- fb$coefficients
      beta <- cf[seq_len(ncol(Xb))]
      xi   <- y[idx] - as.vector(Xb %*% beta)
      rho  <- vapply(seq_along(info$endo_cols),
                     function(k) stats::cor(xi, cb$Cstar[, k]), numeric(1))
      return(c(cf, rho, fo$coefficients, a))
    }
    NULL
  }
  
  ## --- run ------------------------------------------------------------------
  ## parallel: FALSE | TRUE | "multicore" | "snow".  TRUE picks forking on Unix
  ## and macOS and a PSOCK cluster on Windows; the strings force one of the two.
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
      ## fork: workers inherit everything, no export needed.  mc.set.seed
      ## defaults to TRUE, which is what gives each worker its own RNG stream.
      reps <- parallel::mclapply(seq_len(nboots), one, mc.cores = ncores)
    } else {
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      env <- environment(.copreg_boot)
      parallel::clusterExport(cl, .copreg_export_names(env), envir = env)
      ## distinct L'Ecuyer streams per worker.  iseed stays NULL, so no seed is
      ## set: the streams continue from whatever RNG state the caller is in, and
      ## results are reproducible only if the caller set a seed.
      parallel::clusterSetRNGStream(cl)
      reps <- parallel::parLapply(cl, seq_len(nboots), one)
    }
    
    bad  <- vapply(reps, function(r) is.null(r) || inherits(r, "try-error"),
                   logical(1))
    fail <- sum(bad)
    good <- reps[!bad]
    out  <- if (length(good) == 0L) NULL else do.call(rbind, good)
    
  } else {
    
    if (verbose) {
      message("Bootstrapping standard errors (", nboots, " replicates) ...")
      pb <- utils::txtProgressBar(min = 0, max = nboots, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    
    out  <- NULL
    fail <- 0L
    for (b in seq_len(nboots)) {
      r <- one()
      if (is.null(r)) {
        fail <- fail + 1L
      } else {
        if (is.null(out)) out <- matrix(NA_real_, nboots, length(r))
        out[b, ] <- r
      }
      if (verbose) utils::setTxtProgressBar(pb, b)
    }
    if (!is.null(out)) out <- out[stats::complete.cases(out), , drop = FALSE]
  }
  
  if (is.null(out) || nrow(out) == 0L)
    stop("Every bootstrap resample was degenerate after ", max_tries,
         " attempts. This usually means a factor level or a cell of the ",
         "discrete regressors is too rare. Collapse sparse levels or drop the ",
         "variable.", call. = FALSE)
  if (fail > 0L)
    warning(fail, " of ", nboots, " bootstrap replicates could not be ",
            "computed and were dropped.", call. = FALSE)
  
  ## how much resampling was discarded?
  tries <- out[, ncol(out)]
  out   <- out[, -ncol(out), drop = FALSE]
  if (mean(tries) > 1.05)
    warning(sprintf(paste0("%.1f draws were needed per bootstrap replicate ",
                           "(%.0f%% of resamples discarded because a level or ",
                           "cell was missing).\n  The standard errors are ",
                           "conditional on the draws that survived."),
                    mean(tries), 100 * (1 - 1 / mean(tries))), call. = FALSE)
  
  out
}


## =============================================================================
##  4.  Diagnostics
## =============================================================================

## Anderson-Darling test for composite normality (D'Agostino & Stephens 1986).
.ad_test <- function(x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 8L) return(c(A = NA_real_, p = NA_real_))
  p <- stats::pnorm((x - mean(x)) / stats::sd(x))
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  h <- (2 * seq_len(n) - 1) * (log(p) + log(1 - rev(p)))
  A <- -n - mean(h)
  AA <- A * (1 + 0.75 / n + 2.25 / n^2)
  pv <- if (AA < 0.2)   1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  else if (AA < 0.34)  1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  else if (AA < 0.6)   exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  else if (AA < 10)    exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  else                 3.7e-24
  c(A = A, p = pv)
}

## Kolmogorov-Smirnov test against a fitted normal.  Used by both validity
## methods and by the nonnormality table, so it lives here rather than being
## redefined locally each time.
.ks_normal <- function(v) {
  k <- suppressWarnings(stats::ks.test(v, "pnorm", mean(v), stats::sd(v)))
  c(D = unname(k$statistic), p = k$p.value)
}

## Fisher z test for a single correlation.
.fisher_z_p <- function(r, n) {
  if (!is.finite(r) || n < 4L) return(NA_real_)
  r <- min(max(r, -1 + 1e-12), 1 - 1e-12)
  z <- atanh(r) * sqrt(n - 3)
  2 * stats::pnorm(-abs(z))
}

.copreg_diagnostics <- function(y, X, C, Cstar, info, xi, groups = NULL) {
  
  n  <- nrow(X)
  ec <- info$endo_cols
  K  <- length(ec)
  nm <- info$endo_names
  
  ## (a) non-normality of the endogenous regressors --------------------------
  nonnorm <- data.frame(
    AD          = numeric(K),
    `AD p`      = numeric(K),
    `KS p`      = numeric(K),
    row.names   = nm, check.names = FALSE)
  for (k in seq_len(K)) {
    v  <- X[, ec[k]]
    ad <- .ad_test(v)
    nonnorm$AD[k]     <- ad["A"]
    nonnorm$`AD p`[k] <- ad["p"]
    nonnorm$`KS p`[k] <- suppressWarnings(
      stats::ks.test(v, "pnorm", mean(v), stats::sd(v))$p.value)
  }
  
  ## (b) correlation of exogenous regressors with the copula data ------------
  ##     Assumption 5 of Park & Gupta requires corr(W, P*) = 0.  Two things are
  ##     reported: the largest single correlation with a Holm-adjusted p value
  ##     (it is the maximum over many regressors, so the raw p value would be
  ##     far too small), and a joint test over the whole set, which is what the
  ##     assumption is actually about.  The full correlation matrix is kept in
  ##     $exog.correlation.matrix.
  wcols <- setdiff(info$exo_cols, if (info$has_intercept) 1L else integer(0))
  if (length(wcols) > 0L)
    wcols <- wcols[vapply(wcols, function(j) length(unique(X[, j])) > 1L,
                          logical(1))]
  
  corrW <- NULL; corrM <- NULL
  if (length(wcols) > 0L) {
    
    W <- X[, wcols, drop = FALSE]
    corrM <- matrix(NA_real_, K, ncol(W), dimnames = list(nm, colnames(W)))
    for (k in seq_len(K))
      corrM[k, ] <- suppressWarnings(as.vector(stats::cor(W, Cstar[, k])))
    
    Wd <- cbind(1, W)
    corrW <- data.frame(
      `max |corr|` = numeric(K),  `with`     = character(K),
      `p (Holm)`   = numeric(K),  `joint R2` = numeric(K),
      `joint p`    = numeric(K),
      row.names = nm, check.names = FALSE, stringsAsFactors = FALSE)
    
    for (k in seq_len(K)) {
      r  <- corrM[k, ]
      pa <- stats::p.adjust(vapply(r, .fisher_z_p, numeric(1), n = n),
                            method = "holm")
      i  <- which.max(abs(r))
      corrW$`max |corr|`[k] <- r[i]
      corrW$`with`[k]       <- colnames(W)[i]
      corrW$`p (Holm)`[k]   <- pa[i]
      
      ## joint test: regress the copula data on all exogenous regressors
      fw  <- .lm.fit(Wd, Cstar[, k])
      df1 <- fw$rank - 1L
      df2 <- n - fw$rank
      if (df1 >= 1L && df2 >= 1L) {
        r2 <- 1 - sum(fw$residuals^2) /
          sum((Cstar[, k] - mean(Cstar[, k]))^2)
        Fs <- (r2 / df1) / ((1 - r2) / df2)
        corrW$`joint R2`[k] <- r2
        corrW$`joint p`[k]  <- stats::pf(Fs, df1, df2, lower.tail = FALSE)
      } else {
        corrW$`joint R2`[k] <- NA_real_
        corrW$`joint p`[k]  <- NA_real_
      }
    }
  }
  
  ## (c) collinearity between copula terms and the remaining regressors ------
  ##     omega = 1 - R^2 of the copula term on all other regressors.
  ##     For the bivariate model this reduces to 1 - corr(P, C(P))^2, the bias
  ##     inflation factor of Liengaard et al. (2025, Eq. 4 / Eq. 22).
  ## One row per copula column.  There can be more of them than endogenous
  ## regressors: JAMS produces one per endogenous regressor and category of the
  ## discrete covariates, and those columns are zero outside their own cell, so
  ## the correlation is taken over the rows where the column is supported.
  A  <- cbind(X, C)
  kC <- ncol(X) + seq_len(ncol(C))
  src <- if (!is.null(groups)) match(groups$endogenous, nm) else seq_len(ncol(C))
  
  omega <- numeric(ncol(C)); rPC <- numeric(ncol(C))
  for (k in seq_along(kC)) {
    rest <- A[, -kC[k], drop = FALSE]
    fit  <- .lm.fit(rest, A[, kC[k]])
    r2   <- 1 - sum(fit$residuals^2) /
      sum((A[, kC[k]] - mean(A[, kC[k]]))^2)
    omega[k] <- 1 - r2
    sup <- if (ncol(C) > K) which(C[, k] != 0) else seq_len(nrow(C))
    rPC[k] <- if (length(sup) > 2L)
      suppressWarnings(stats::cor(X[sup, ec[src[k]]], C[sup, k])) else NA_real_
  }
  coll <- data.frame(`corr(P, C)` = rPC, omega = omega,
                     row.names = colnames(C), check.names = FALSE)
  
  list(nonnormality = nonnorm, exog.correlation = corrW,
       exog.correlation.matrix = corrM, collinearity = coll)
}


## Columns of the exogenous formula part(s) usable as W in a first stage.
##
## A term qualifies when none of the variables in it is endogenous.  That keeps
## an interaction with an endogenous regressor out -- it is not exogenous
## information, whichever side of the "|" it was written on -- while letting
## genuine transformations of the exogenous regressors in, including w1:w2 and
## I(w^2).  Breitung, Mayer & Wied (2024) recommend exactly that under their
## Assumption A4: their first stage is linear in x, and the way to loosen it is
## to add nonlinear transformations of the exogenous regressors to x.
.copreg_exog_cols <- function(info, X = info$X) {
  
  asg <- attr(info$X, "assign")
  lab <- character(ncol(info$X))
  lab[asg > 0L] <- attr(info$terms, "term.labels")[asg[asg > 0L]]
  
  ## A variable counts as endogenous when it turns up in an endogenous term and
  ## is not itself a first-order term of the exogenous part.  Writing
  ## y ~ p + p:m | m + w therefore makes p endogenous while leaving m exogenous,
  ## so the first stage uses m and w but not p:m.
  raw <- unique(unlist(lapply(info$endo_names, function(l)
    all.vars(stats::reformulate(lab[match(l, colnames(info$X))])))))
  exo_main <- unique(lab[info$part >= 2L & info$ord == 1L])
  endo_vars <- setdiff(raw, exo_main)
  
  ok <- info$part >= 2L &
    vapply(seq_along(lab), function(j) nzchar(lab[j]) &&
             !any(all.vars(stats::reformulate(lab[j])) %in% endo_vars),
           logical(1))
  wc <- which(ok)
  if (length(wc) > 0L)
    wc <- wc[vapply(wc, function(j) length(unique(X[, j])) > 1L, logical(1))]
  wc
}


## =============================================================================
##  4b.  Two-stage copula control functions
## =============================================================================
##
##  Shared by 2sCOPE (Yang, Qian & Xie 2025) and IMA (Haschka 2025).  Both
##  transform every regressor, endogenous and exogenous, to normal scores and
##  then residualise each endogenous P*_k on the exogenous W* only, never on
##  the other endogenous regressors:
##
##      C_k = P*_k - delta_k' W*.
##
##  The two differ in exactly one place.  Equation 11 of Yang et al. carries no
##  intercept in the first stage, and IMA implements it that way, because the
##  slope then equals the Pearson correlation between P* and W*, which is what
##  Haschka uses to recover rho.  2sCOPE adds an intercept: note 8 of Qian,
##  Koschmann & Xie (2025) reports that the implementation includes it because
##  it is more general and performs well in their simulations.
##
##  .ctor_twostage(intercept) returns the constructor for either variant.

## Which of the exogenous columns actually carry independent information after
## the copula transformation?
##
## The rank-based CDFs are invariant to strictly monotone transformations, so
## Phi^-1(Fhat(w)) and Phi^-1(Fhat(I(w^2))) are the same column whenever w is
## positive, and the first stage would be singular.  The kernel-based CDFs make
## them near rather than exactly collinear, which is no better.  Redundant
## columns are therefore dropped from the first stage -- they stay in the
## structural model -- and the choice is made once and reused, so that the
## design does not change between bootstrap replicates.
.copreg_firststage_cols <- function(X, info, cdf, ties) {
  
  wc <- .copreg_exog_cols(info, X)
  if (length(wc) < 2L) return(wc)
  
  Ws <- .copula_transform_matrix(X[, wc, drop = FALSE], cdf, ties)
  Wc <- sweep(Ws, 2L, colMeans(Ws), "-")
  qz <- qr(Wc, tol = 1e-7)
  if (qz$rank == length(wc)) return(wc)
  
  sel  <- sort(qz$pivot[seq_len(qz$rank)])
  drop <- setdiff(seq_along(wc), sel)
  message("Left out of the first stage because they carry no independent ",
          "information after\n  the copula transformation, which is invariant ",
          "to monotone transformations: ",
          paste(colnames(X)[wc[drop]], collapse = ", "),
          ".\n  They remain regressors in the structural model.")
  wc[sel]
}


.ctor_twostage <- function(intercept = TRUE) {
  
  cols <- NULL
  
  function(X, info, cdf, ties) {
    
    if (is.null(cols)) cols <<- .copreg_firststage_cols(X, info, cdf, ties)
    wc <- cols
    
    P     <- X[, info$endo_cols, drop = FALSE]
    Cstar <- .copula_transform_matrix(P, cdf, ties)
    colnames(Cstar) <- info$endo_names
    
    ## no exogenous regressors: the first stage has nothing to project on and
    ## the estimator collapses to Park & Gupta.  The wrappers redirect before
    ## reaching this point; the branch only guards the bootstrap.
    if (length(wc) == 0L) {
      C <- Cstar
      colnames(C) <- paste0(info$endo_names, "_cop")
      return(list(C = C, Cstar = Cstar, Wstar = NULL))
    }
    
    Wstar <- .copula_transform_matrix(X[, wc, drop = FALSE], cdf, ties)
    D <- if (intercept) cbind(`(Intercept)` = 1, Wstar) else Wstar
    
    C <- matrix(NA_real_, nrow(X), ncol(Cstar))
    for (k in seq_len(ncol(Cstar))) {
      f <- .lm.fit(D, Cstar[, k])
      if (f$rank < ncol(D))
        stop("The first stage is rank deficient.", call. = FALSE)
      C[, k] <- f$residuals
    }
    colnames(C) <- paste0(info$endo_names, "_cop")
    
    list(C = C, Cstar = Cstar, Wstar = Wstar)
  }
}


## =============================================================================
##  5.  Fitting engine
## =============================================================================
##
##  ctor(X, info, cdf, ties) is the estimator-specific constructor.  It returns
##    C     n x K matrix of copula control functions entering the regression
##    Cstar n x K matrix of the plain copula data Phi^{-1}(Fhat(P)), used for
##          the endogeneity measure rho = corr(xi, P*) and for Assumption 5(b)
##    Wstar optional, the copula transformation of the exogenous regressors;
##          kept from the full-sample call so that validity() can redo the
##          first stage.  The bootstrap ignores it.
##  For Park & Gupta C and Cstar coincide; for the two-stage estimators not.

##  id.condition names the identification condition the estimator's own paper
##  states, which decides what summary() and validity() report:
##    "nonnormality"  the endogenous regressor (or, for BMW, the first-stage
##                    error) must be nonnormal, and rho = corr(xi*, P*) is a
##                    meaningful summary of the regressor-error dependence
##    "assumption3"   the copula term must not be a linear combination of
##                    (1, P, W), which is Assumption 3 of Hu, Qian & Xie
##                    (2025).  There, rho is not separately identified -- only
##                    the product sigma_o * rho, which is the coefficient of
##                    the copula term itself -- so no rho is reported.
##
##  cdf and ties may be NA for estimators that use no marginal CDF at all.

.copreg_fit <- function(formula, data, ctor, method, cdf, ties, nboots,
                        subset = NULL, contrasts = NULL,
                        parallel = FALSE, ncores = NULL, assumption5 = TRUE,
                        dhw = FALSE, id.condition = "nonnormality",
                        verbose = interactive(), call = NULL) {
  
  cdf  <- if (is.na(cdf))  NA_character_ else match.arg(cdf,  .copreg_cdf_choices)
  ties <- if (is.na(ties)) NA_character_ else match.arg(ties, .copreg_ties_choices)
  id.condition <- match.arg(id.condition, c("nonnormality", "assumption3"))
  if (!is.numeric(nboots) || length(nboots) != 1L || nboots < 2)
    stop("'nboots' must be a single number >= 2.", call. = FALSE)
  nboots <- as.integer(nboots)
  
  info <- .copreg_model(formula, data, subset, contrasts)
  y <- info$y; X <- info$X
  
  cc <- ctor(X, info, cdf, ties)
  C  <- cc$C; Cstar <- cc$Cstar
  
  ## Two quite different failures land here, so say which one it is.
  A <- cbind(X, C)
  if (qr(A)$rank < ncol(A)) {
    if (qr(X)$rank < ncol(X)) {
      dup <- colnames(X)[-sort(qr(X, tol = 1e-7)$pivot[seq_len(qr(X)$rank)])]
      stop("The design matrix is rank deficient before any copula term is ",
           "added: the regressors\n  themselves are collinear. Redundant ",
           "column(s): ", paste(dup, collapse = ", "),
           ".\n  A constant column next to an intercept, or a duplicated ",
           "regressor, is the usual cause.", call. = FALSE)
    }
    stop("The copula terms are perfectly collinear with the regressors, so the ",
         "model is not\n  identified. This is what happens when an endogenous ",
         "regressor is normally\n  distributed: its copula transform is then a ",
         "linear function of itself.", call. = FALSE)
  }
  
  fit  <- .lm.fit(A, y)
  cf   <- fit$coefficients
  names(cf) <- colnames(A)
  
  ## classical OLS covariance of the augmented regression.  These standard
  ## errors are wrong for the structural coefficients -- the copula terms are
  ## generated regressors -- but Breitung, Mayer & Wied (2024, Corollary 3.2)
  ## show that the textbook t statistic remains valid for testing rho = 0, so
  ## the estimators that rely on that result need them.
  Vc <- tryCatch({
    s2 <- sum(fit$residuals^2) / (nrow(A) - ncol(A))
    V0 <- s2 * chol2inv(qr.R(qr(A)))
    dimnames(V0) <- list(names(cf), names(cf))
    V0
  }, error = function(e) NULL)
  
  beta  <- cf[seq_len(ncol(X))]
  gamma <- cf[ncol(X) + seq_len(ncol(C))]
  
  ## Names follow lm(): the row names of the model frame, so that observations
  ## stay identifiable after na.omit has dropped some.
  nmr <- rownames(info$mf)
  y   <- stats::setNames(as.vector(y), nmr)
  fit_s   <- stats::setNames(as.vector(X %*% beta), nmr)   # mu + P alpha + W beta
  resid_s <- y - fit_s                                      # xi
  fit_a   <- stats::setNames(as.vector(A %*% cf), nmr)      # ... + C gamma
  resid_a <- y - fit_a                                      # u
  
  rho <- vapply(seq_along(info$endo_cols),
                function(k) stats::cor(resid_s, Cstar[, k]), numeric(1))
  names(rho) <- info$endo_names
  
  ## --- bootstrap ------------------------------------------------------------
  B  <- .copreg_boot(y, X, info, ctor, cdf, ties, nboots,
                     verbose = verbose, parallel = parallel, ncores = ncores)
  kp <- length(cf)
  kx <- ncol(X)
  colnames(B) <- c(names(cf), paste0("rho.", info$endo_names),
                   paste0("ols.", colnames(X)))
  
  Bc <- B[, seq_len(kp), drop = FALSE]
  Br <- B[, kp + seq_along(rho), drop = FALSE]
  Bo <- B[, kp + length(rho) + seq_len(kx), drop = FALSE]
  
  V  <- stats::cov(Bc)
  dimnames(V) <- list(names(cf), names(cf))
  se     <- sqrt(diag(V))
  rho.se <- apply(Br, 2, stats::sd)
  
  ## uncorrected OLS on the same resamples, and the resulting ICON statistic
  ## (Qian, Koschmann & Xie 2025, Boundary Condition 1): the inflation of the
  ## standard errors caused by adding the copula terms.  Values above 6 flag
  ## weak identification or a misspecified dependence model.
  ols.cf <- .lm.fit(X, y)$coefficients
  names(ols.cf) <- colnames(X)
  ols.se <- apply(Bo, 2, stats::sd)
  names(ols.se) <- colnames(X)
  icon <- se[colnames(X)] / ols.se
  names(icon) <- colnames(X)
  
  ## --- assemble -------------------------------------------------------------
  structure(list(
    call            = call,
    method          = method,
    cdf             = cdf,
    ties            = ties,
    nboots          = nrow(B),
    parallel        = if (isFALSE(parallel)) "no" else parallel,
    coefficients    = cf,
    std.error       = se,
    vcov            = V,
    boot            = Bc,
    boot.rho        = Br,
    boot.ols        = Bo,
    ols.coefficients = ols.cf,
    ols.std.error   = ols.se,
    icon            = icon,
    rho             = rho,
    rho.se          = rho.se,
    fitted.values   = fit_s,
    residuals       = resid_s,
    fitted.augmented   = fit_a,
    residuals.augmented = resid_a,
    copula.terms    = C,
    copula.data     = Cstar,
    copula.groups   = cc$groups,
    Wstar           = cc$Wstar,
    first.stage.residuals = cc$resid1,
    ols.vcov        = Vc,
    y               = y,
    X               = X,
    endogenous      = info$endo_names,
    exogenous       = setdiff(colnames(X), info$endo_names),
    n               = nrow(X),
    df.residual     = nrow(X) - ncol(A),
    terms           = info$terms,
    formula         = info$Formula,
    xlevels         = info$xlevels,
    contrasts       = info$contrasts,
    na.action       = info$na.action,
    has.intercept   = info$has_intercept,
    assumption5     = assumption5,
    dhw             = dhw,
    id.condition    = id.condition,
    bandwidths      = cc$bandwidths,
    diagnostics     = .copreg_diagnostics(y, X, C, Cstar, info, resid_s,
                                          cc$groups)
  ), class = "copreg")
}


## =============================================================================
##  6.  S3 methods
## =============================================================================

## Bootstrap Wald test of R gamma = 0, using the bootstrap covariance matrix.
## Liengaard et al. (2025) propose this for the copula terms' joint
## significance: with several copula terms, testing them one at a time is a
## multiple testing problem.
## Copula terms are frequently close to collinear, so the covariance matrix of
## the contrasts can be near singular.  Inverting it with solve() would then
## produce an unstable statistic on a degrees-of-freedom count that overstates
## what the data actually identify.  We use the Moore-Penrose inverse instead
## and report its effective rank as the degrees of freedom.
.wald <- function(cf, V, idx, R = NULL, tol = 1e-8) {
  g <- cf[idx]
  if (is.null(R)) R <- diag(length(g))
  Rg <- as.vector(R %*% g)
  RV <- R %*% V[idx, idx, drop = FALSE] %*% base::t(R)
  
  sv  <- svd(RV)
  pos <- sv$d > tol * max(sv$d)
  q   <- sum(pos)
  if (q == 0L) return(c(chisq = NA_real_, df = 0, p = NA_real_))
  Wi <- sv$v[, pos, drop = FALSE] %*% ((1 / sv$d[pos]) *
                                         t(sv$u[, pos, drop = FALSE]))
  st <- as.vector(crossprod(Rg, Wi %*% Rg))
  c(chisq = st, df = q, p = stats::pchisq(st, q, lower.tail = FALSE))
}


coef.copreg  <- function(object, ...) object$coefficients
vcov.copreg  <- function(object, ...) object$vcov
nobs.copreg  <- function(object, ...) object$n

## Notation follows Qian, Koschmann & Xie (2025, Eq. 4):
##   y = mu + P alpha + W beta + C gamma + u,
## with P the endogenous regressors, W the exogenous ones and C the copula
## control functions.  The structural error is xi = y - mu - P alpha - W beta,
## so that xi still contains the endogenous part that C is meant to absorb.
##
## residuals(type = "structural")  xi = y - mu - P alpha - W beta
## residuals(type = "augmented")   u  = xi - C gamma
residuals.copreg <- function(object, type = c("structural", "augmented"), ...) {
  type <- match.arg(type)
  if (type == "structural") object$residuals else object$residuals.augmented
}

fitted.copreg <- function(object, type = c("structural", "augmented"), ...) {
  type <- match.arg(type)
  if (type == "structural") object$fitted.values else object$fitted.augmented
}

## Predictions use the structural model only:  y = mu + P alpha + W beta.
## Copula terms are endogeneity controls, not part of the causal model.
predict.copreg <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  Terms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(Terms, newdata, na.action = stats::na.pass,
                           xlev = object$xlevels)
  X  <- stats::model.matrix(Terms, mf, contrasts.arg = object$contrasts)
  b  <- object$coefficients[colnames(X)]
  if (anyNA(b))
    stop("newdata produced columns that are not in the fitted model.",
         call. = FALSE)
  drop(X %*% b)
}

## update() on a fitted model.  Multi-part formulas need Formula's update
## method; stats::update.formula does not understand the "|" separator, so
## update(fit, . ~ . | . + z) would otherwise fail.
update.copreg <- function(object, formula., ..., evaluate = TRUE) {
  cl <- object$call
  if (is.null(cl))
    stop("The fitted object carries no call to update.", call. = FALSE)
  ## inline the stored Formula rather than the symbol the user happened to
  ## pass, so that update() does not depend on that variable still existing
  cl$formula <- if (missing(formula.)) object$formula
  else stats::update(object$formula, formula.)
  extras <- match.call(expand.dots = FALSE)$`...`
  
  ## Switching estimator: the wrappers take no 'method' argument, so route the
  ## refit through copreg().  Estimator specific arguments have to be left
  ## behind -- bwmethod means nothing to CopRegPG(), and an estimator with a
  ## "..." would silently pass them on to something else -- so only the
  ## arguments the target function actually declares survive the switch.  A cdf
  ## the user set explicitly is among them wherever the target has one; if none
  ## was set, copreg() resolves to the new estimator's own default.
  if ("method" %in% names(extras)) {
    m   <- tryCatch(eval(extras$method, parent.frame()), error = function(e) NULL)
    reg <- if (is.character(m) && length(m) == 1L) .copreg_registry[[m]] else NULL
    fun <- if (!is.null(reg) && exists(reg$fun, mode = "function"))
      get(reg$fun, mode = "function") else NULL
    if (!is.null(fun)) {
      keep <- c(setdiff(names(formals(fun)), "..."), "method")
      nm   <- names(cl)
      drop <- which(nzchar(nm) & !(nm %in% keep))
      if (length(drop) > 0L) {
        message("Dropped when switching estimator, because ", reg$fun,
                "() has no such argument: ",
                paste(nm[drop], collapse = ", "), ".")
        cl[drop] <- NULL
      }
    }
    if (!identical(cl[[1L]], quote(copreg))) cl[[1L]] <- quote(copreg)
  }
  
  if (length(extras) > 0L) {
    exist <- !is.na(match(names(extras), names(cl)))
    for (a in names(extras)[exist]) cl[[a]] <- extras[[a]]
    if (any(!exist)) {
      cl <- c(as.list(cl), extras[!exist])
      cl <- as.call(cl)
    }
  }
  if (evaluate) eval(cl, parent.frame()) else cl
}

formula.copreg <- function(x, ...) x$formula

confint.copreg <- function(object, parm, level = 0.95,
                           type = c("normal", "percentile"), ...) {
  type <- match.arg(type)
  cf <- object$coefficients
  if (missing(parm)) parm <- names(cf)
  a  <- (1 - level) / 2
  if (type == "normal") {
    q  <- stats::qnorm(1 - a)
    ci <- cbind(cf - q * object$std.error, cf + q * object$std.error)
  } else {
    ci <- base::t(apply(object$boot, 2, stats::quantile,
                        probs = c(a, 1 - a), names = FALSE))
  }
  dimnames(ci) <- list(names(cf), paste0(round(100 * c(a, 1 - a), 1), " %"))
  ci[parm, , drop = FALSE]
}

print.copreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCopula endogeneity correction:", x$method, "\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n"); print(x$call)
  }
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

summary.copreg <- function(object, ...) {
  cf <- object$coefficients
  se <- object$std.error
  z  <- cf / se
  p  <- 2 * stats::pnorm(-abs(z))
  tab <- cbind(Estimate = cf, `Std. Error` = se,
               `z value` = z, `Pr(>|z|)` = p)
  
  ## Durbin-Hausman-Wu test of rho = 0.  Breitung, Mayer & Wied (2024,
  ## Corollary 3.2) show that under the null the textbook t statistic built
  ## from classical OLS standard errors keeps a standard normal limit, even
  ## though those standard errors are inconsistent for the structural
  ## coefficients.  Only the estimators whose theory covers this report it.
  dhw <- NULL
  if (isTRUE(object$dhw) && !is.null(object$ols.vcov)) {
    ci <- match(colnames(object$copula.terms), names(cf))
    g  <- cf[ci]
    se0 <- sqrt(diag(object$ols.vcov))[ci]
    dhw <- cbind(Estimate = g, `Std. Error` = se0, `t value` = g / se0,
                 `Pr(>|t|)` = 2 * stats::pnorm(-abs(g / se0)))
    rownames(dhw) <- names(cf)[ci]
  }
  
  ## rho is reported only where the estimator identifies it.  Under the
  ## conditional copula model of Hu, Qian & Xie (2025) the control function is
  ## E(o|p, w) = sigma_o * rho * C_p, so only the product enters the
  ## regression: rho and sigma_o are not separately identified, and the paper
  ## accordingly reports the coefficient of the copula term instead.
  rtab <- NULL
  if (identical(object$id.condition, "nonnormality")) {
    rtab <- cbind(Estimate = object$rho, `Std. Error` = object$rho.se,
                  `z value` = object$rho / object$rho.se,
                  `Pr(>|z|)` = 2 * stats::pnorm(-abs(object$rho / object$rho.se)))
    rownames(rtab) <- paste0("rho(", names(object$rho), "*, xi*)")
  }
  
  ## --- Wald tests on the copula terms ---------------------------------------
  ci   <- match(colnames(object$copula.terms), names(cf))
  grp  <- object$copula.groups
  wald <- NULL
  if (length(ci) >= 2L)
    wald <- .wald(cf, object$vcov, ci)
  
  ## with several categories of the discrete covariates, test whether the
  ## copula structure differs across them: gamma^z_k equal in z, for each k
  ##
  ## Two layouts occur.  JAMS carries one coefficient per cell, so equality
  ## across cells is a set of pairwise contrasts.  2sCOPE-np follows Hu, Qian
  ## & Xie (2025) and parameterises the groups against a reference group, so
  ## the deviation coefficients already are the differences and the test is
  ## simply that all of them vanish.
  wald.z <- NULL
  if (!is.null(grp) && length(unique(grp$cell)) > 1L) {
    R <- NULL
    if (isTRUE(attr(grp, "reference"))) {
      j <- which(nzchar(grp$cell))
      if (length(j) > 0L) {
        R <- matrix(0, length(j), length(ci))
        R[cbind(seq_along(j), j)] <- 1
      }
    } else {
      for (v in unique(grp$endogenous)) {
        j <- which(grp$endogenous == v)
        for (m in seq_along(j)[-1L]) {
          r <- numeric(length(ci)); r[j[1L]] <- 1; r[j[m]] <- -1
          R <- rbind(R, r)
        }
      }
    }
    if (!is.null(R)) wald.z <- .wald(cf, object$vcov, ci, R)
  }
  
  ra  <- object$residuals.augmented          # u  = xi - C gamma
  rs  <- object$residuals                    # xi = y - mu - P alpha - W beta
  df  <- object$df.residual
  y   <- object$fitted.augmented + ra
  ss  <- if (object$has.intercept) sum((y - mean(y))^2) else sum(y^2)
  nn  <- object$n - as.integer(object$has.intercept)
  
  ## Both sets use the same residual degrees of freedom: all coefficients,
  ## including the copula terms, are estimated jointly in the augmented
  ## regression, so that is the number of parameters spent either way.
  fitstat <- function(r) {
    r2 <- 1 - sum(r^2) / ss
    c(`Residual standard error` = sqrt(sum(r^2) / df),
      `R-squared` = r2,
      `Adjusted R-squared` = 1 - (1 - r2) * nn / df)
  }
  fit.tab <- cbind(augmented = fitstat(ra), structural = fitstat(rs))
  
  structure(list(call = object$call, method = object$method,
                 cdf = object$cdf, ties = object$ties, nboots = object$nboots,
                 coefficients = tab, rho = rtab, dhw = dhw,
                 wald = wald, wald.z = wald.z,
                 wald.z.reference = isTRUE(attr(grp, "reference")),
                 id.condition = object$id.condition,
                 bandwidths = object$bandwidths,
                 n.cells = if (is.null(grp)) 1L else length(unique(grp$cell)),
                 residuals = ra, residuals.structural = rs,
                 fit = fit.tab, sigma = fit.tab["Residual standard error", 2],
                 r.squared = fit.tab["R-squared", 1],
                 adj.r.squared = fit.tab["Adjusted R-squared", 1],
                 df.residual = df,
                 n = object$n, endogenous = object$endogenous,
                 assumption5 = object$assumption5,
                 copula.names = colnames(object$copula.terms),
                 diagnostics = object$diagnostics),
            class = "summary.copreg")
}

print.summary.copreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  
  cat("\nCopula endogeneity correction:", x$method, "\n")
  if (!is.null(x$call)) { cat("\nCall:\n"); print(x$call) }
  
  qn <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("\nResiduals of the augmented regression (u = xi - C gamma):\n")
  print(structure(stats::quantile(x$residuals), names = qn), digits = digits)
  cat("\nResiduals of the structural model (xi = y - mu - P alpha - W beta):\n")
  print(structure(stats::quantile(x$residuals.structural), names = qn),
        digits = digits)
  
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, has.Pvalue = TRUE,
                      P.values = TRUE, na.print = "NA")
  
  if (!is.null(x$rho)) {
    cat("\nEndogeneity: rho(P*, xi*) is the correlation between the normal score",
        "\n  of an endogenous regressor and that of the structural error,",
        "xi* = xi / sigma.\n")
    stats::printCoefmat(x$rho, digits = digits, signif.stars = signif.stars,
                        has.Pvalue = TRUE, P.values = TRUE)
  } else {
    cat("\nEndogeneity is read off the copula terms above. The control function is\n",
        "  E(o | p, w) = sigma_o * rho * C_p, so only the product sigma_o * rho\n",
        "  enters the regression and rho is not separately identified; a copula\n",
        "  coefficient of zero means no endogeneity.\n", sep = "")
  }
  
  if (!is.null(x$dhw)) {
    cat("\nDurbin-Hausman-Wu test of no endogeneity, rho = 0. Under that null the\n",
        "  textbook t statistic stays valid even though the copula terms are\n",
        "  generated regressors (Breitung, Mayer & Wied 2024, Corollary 3.2), so\n",
        "  these standard errors are classical rather than bootstrapped:\n", sep = "")
    stats::printCoefmat(x$dhw, digits = digits, signif.stars = signif.stars,
                        has.Pvalue = TRUE, P.values = TRUE)
  }
  
  if (!is.null(x$wald)) {
    cat("\nJoint significance of the copula terms (bootstrap Wald test):\n")
    cat("  chi-squared = ", formatC(x$wald["chisq"], digits = digits),
        " on ", x$wald["df"], " df, p = ",
        format.pval(x$wald["p"], digits = digits), "\n", sep = "")
  }
  if (!is.null(x$wald.z)) {
    cat(if (isTRUE(x$wald.z.reference))
      paste0("Endogeneity constant across the ", x$n.cells,
             " groups, i.e. all deviations from the\n  reference group are zero")
      else
        paste0("Copula structure constant across the ", x$n.cells,
               " categories of the discrete\n  regressors"),
      " (bootstrap Wald test): ",
      "chi-squared = ", formatC(x$wald.z["chisq"], digits = digits),
      " on ", x$wald.z["df"], " df,\n  p = ",
      format.pval(x$wald.z["p"], digits = digits),
      ". A small p favours the grouped estimation.\n", sep = "")
  }
  
  cat("\nFit, on", x$df.residual, "residual degrees of freedom:\n")
  print(format(x$fit, digits = digits), quote = FALSE)
  if (!is.null(x$rho))
    cat("  sigma above is the standard error of the structural model, the one",
        "\n  entering xi* = xi / sigma.\n")
  cat("Standard errors from ", x$nboots, " bootstrap replicates",
      if (is.na(x$cdf)) ".\n" else
        paste0("; cdf = \"", x$cdf, "\", ties = \"", x$ties, "\".\n"), sep = "")
  if (!is.null(x$bandwidths)) {
    cat("Bandwidths of the conditional CDF estimation:\n")
    print(format(x$bandwidths, digits = digits), quote = FALSE)
  }
  cat("Pr(>|z|)", if (is.null(x$rho)) "" else " in both tables",
      ": Wald test using the normal approximation,\n", sep = "")
  cat("  z = Estimate / Std. Error, with the bootstrap standard error.\n",
      "  See confint(object, type = \"percentile\") for bootstrap percentile ",
      "intervals.\n", sep = "")
  
  d <- x$diagnostics
  cat("\n--- Identification diagnostics ------------------------------------\n")
  if (identical(x$id.condition, "nonnormality")) {
    cat("\nNon-normality of the endogenous regressors (small p = non-normal, good):\n")
    print(format(d$nonnormality, digits = digits))
  }
  if (!is.null(d$exog.correlation)) {
    cat("\nCorrelation of the copula terms with the exogenous regressors\n",
        if (isTRUE(x$assumption5))
          "(Park & Gupta assume this is zero; 'joint' tests all of them at once,\n the Holm p value refers to the single largest correlation):\n"
        else if (identical(x$id.condition, "assumption3"))
          "(the conditional CDF conditions on them, so a well estimated one leaves\n these correlations near zero; large values point at a poorly estimated\n conditional CDF rather than at a violated assumption):\n"
        else
          "(this estimator does not assume it is zero; it projects the correlation\n out in the first stage):\n",
        sep = "")
    print(format(d$exog.correlation, digits = digits))
    cat(" full matrix in summary(object)$diagnostics$exog.correlation.matrix\n")
  }
  cat("\nCollinearity of the copula terms (omega near 0 = weakly identified):\n")
  print(format(d$collinearity, digits = digits))
  cat("\n")
  invisible(x)
}


## =============================================================================
##  Documentation of the cdf argument
## =============================================================================
##
##  cdf = "kde.silverman"  F1.  Epanechnikov kernel, Silverman bandwidth
##                         b = 0.9 n^(-1/5) min(s, IQR/1.34), integrated in
##                         closed form.
##                         Park & Gupta (2012), Eq. 3.
##                         Note: Liengaard et al. (2025) point out that this
##                         bandwidth is MSE optimal only under normality, which
##                         is exactly what the copula approach rules out.
##
##  cdf = "kde.cv"         F2.  Epanechnikov kernel, cross-validated bandwidth.
##                         Li, Li & Racine (2017); requires the np package.
##
##  cdf = "kde.plugin"     Gaussian kernel, Polansky & Baker (2000) plug-in
##                         bandwidth for the distribution function.  Same
##                         estimator as ks::kcde, reimplemented without the
##                         dependency.  Default for the panel estimator, where
##                         Haschka (2022) uses a kernel density estimate.
##
##  cdf = "ecdf.fixed"     F3.  ECDF with the endpoints replaced by 1e-7 and
##                         1 - 1e-7.
##                         Becker, Proksch & Ringle (2022); R package REndo.
##
##  cdf = "ecdf.adj"       F4.  1/(2n) + (n-1)/n^2 * sum I(P_i <= x).
##                         Liengaard et al. (2025), Eq. 9.
##
##  cdf = "rank.n"         Rank(P)/n, with the maximum mapped to n/(n+1).
##                         Yang, Qian & Xie (2025);
##                         Qian, Koschmann & Xie (2025), Eq. 9.
##
##  cdf = "rank.n1"        Rank(P)/(n+1).
##                         Breitung, Mayer & Wied (2024), Eq. 2.3.
##
##  ties = "max"           sum I(P_i <= x), i.e. the counting function exactly
##                         as written in all of the papers above.  Default.
##  ties = "average"       midranks, the convention of the wider copula
##                         literature and of copula::pobs().  Only relevant for
##                         variables with tied values.  Note that with midranks
##                         the endpoint corrections of "ecdf.fixed" and "rank.n"
##                         never trigger for a tied maximum.
##
##  For continuous regressors without ties the two settings are identical.  For
##  binary dummies they differ only by an affine transformation and leave the
##  estimates unchanged.  They start to matter for discrete regressors with
##  three or more levels.


## =============================================================================
##  7.  Normality diagnostics used by validity()
## =============================================================================

## Cramer-von Mises test for composite normality (Stephens 1986).
.cvm_test <- function(x) {
  x <- sort(x[is.finite(x)]); n <- length(x)
  if (n < 8L) return(c(W = NA_real_, p = NA_real_))
  p <- stats::pnorm((x - mean(x)) / stats::sd(x))
  W <- 1 / (12 * n) + sum((p - (2 * seq_len(n) - 1) / (2 * n))^2)
  WW <- W * (1 + 0.5 / n)
  pv <- if (WW < 0.0275)     1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
  else if (WW < 0.051)      1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)
  else if (WW < 0.092)      exp(0.886 - 31.62 * WW + 10.897 * WW^2)
  else if (WW < 1.1)        exp(1.111 - 34.242 * WW + 12.832 * WW^2)
  else                      7.37e-10
  c(W = W, p = pv)
}

.skewness <- function(x) {
  x <- x[is.finite(x)]; n <- length(x); m <- x - mean(x)
  (sum(m^3) / n) / (sum(m^2) / n)^1.5
}

.ex_kurtosis <- function(x) {
  x <- x[is.finite(x)]; n <- length(x); m <- x - mean(x)
  (sum(m^4) / n) / (sum(m^2) / n)^2 - 3
}


## =============================================================================
##  8.  Validity check
## =============================================================================
##
##  Brings together the checkpoints of three sources.  They do not fully agree,
##  and the output keeps them apart rather than merging them into one verdict.
##
##  Becker, Proksch & Ringle (2022, Figure 8)
##    Boundary conditions for the Park & Gupta method in models WITH an
##    intercept.  Their C5.0 decision tree on Study 4 gives, for a target power
##    of 80% for the copula term:
##      n <= 200                    no distribution in their design suffices
##      400 <= n <= 1000            |skewness| >= 1.932
##      n > 1000                    |skewness| > 0.774
##      n > 2000                    any skewness
##      or                          AD > 18.964, or CvM > 3.488
##    For 90% power: n > 600 at |skewness| > 1.974, n > 2000 at > 0.998;
##    AD > 67.875 for n <= 2000 and AD > 46.832 for n > 2000.
##    They also require normality of the error term: with a nonnormal error the
##    method is no longer consistent in models with an intercept.
##
##  Yang, Qian & Xie (2025, Figure 2)
##    Step 1  The uncorrelatedness assumption: is the copula transformation
##            term correlated with any exogenous regressor?  Printed only for
##            estimators that rely on it; the two-stage estimators project the
##            correlation out in the first stage.  The consequence of a
##            violation was derived by Haschka (2025).  The assumption carries
##            a different number in each paper, so the output just calls it
##            "Assumption".
##    Step 2  KS test of nonnormality of P, threshold p < .05.
##    Step 3  Where Step 2 fails: a continuous exogenous regressor that is
##            sufficiently nonnormal (KS p < .001) and relevant (F > 10).
##
##  Qian, Koschmann & Xie (2025, Table 4)
##    ICON, the inflation of the standard errors relative to uncorrected OLS.
##    Values above 6 flag weak identification or a misspecified dependence
##    model.
##
##  Two disagreements are worth knowing about and are printed with the output.
##  First, on which normality test to use: Becker et al. report a correlation
##  of .66 between the AD statistic and the copula term's t statistic against
##  -.03 for KS, and therefore recommend AD or CvM; Yang et al. prefer KS
##  precisely because it is conservative.  Second, on the error term: Becker
##  et al. find nonnormal errors to break consistency, while Yang et al. and
##  Qian et al. permit them under the decomposition xi = U + V into a normal
##  endogenous part U that carries the whole regressor dependence and an
##  independent nonnormal V.  Since the residuals show xi rather than U, their
##  shape settles neither question.

## Becker, Proksch & Ringle (2022) thresholds for a target power of the copula
## term, as a function of the sample size.  Factored out because their Study 3
## covers fixed-effects panel models too and concludes that the same
## recommendations apply once the total sample size is the reference, so
## validity() uses them for both model classes.
.becker_thresholds <- function(n, power) {
  if (!power %in% c(0.8, 0.9))
    stop("'power' must be 0.8 or 0.9; Becker et al. report thresholds for ",
         "these two levels only.", call. = FALSE)
  if (power == 0.8)
    c(skewness = if (n <= 200) Inf else if (n <= 1000) 1.932
      else if (n <= 2000) 0.774 else 0,
      AD = 18.964, CvM = 3.488)
  else
    c(skewness = if (n <= 600) Inf else if (n <= 2000) 1.974 else 0.998,
      AD  = if (n <= 2000) 67.875 else 46.832,
      CvM = if (n <= 2000) 12.246 else 7.994)
}

## Non-normality table for a set of columns, shared by both validity methods.
.validity_nonnormality <- function(M, nm, n, power) {
  th <- .becker_thresholds(n, power)
  sk <- vapply(nm, function(v) .skewness(M[, v]),    numeric(1))
  ad <- vapply(nm, function(v) .ad_test(M[, v])["A"],  numeric(1))
  cv <- vapply(nm, function(v) .cvm_test(M[, v])["W"], numeric(1))
  kp <- vapply(nm, function(v) .ks_normal(M[, v])["p"], numeric(1))
  list(table = data.frame(
    skewness      = round(sk, 3),
    `ex.kurtosis` = round(vapply(nm, function(v) .ex_kurtosis(M[, v]),
                                 numeric(1)), 3),
    AD            = round(ad, 3),
    CvM           = round(cv, 3),
    `KS p`        = format.pval(kp, digits = 3, eps = 1e-16),
    `Yang ok`     = kp < 0.05,
    `Becker ok`   = (abs(sk) >= th["skewness"]) | (ad > th["AD"]) |
      (cv > th["CvM"]),
    row.names = nm, check.names = FALSE, stringsAsFactors = FALSE),
    thresholds = th)
}


validity <- function(object, ...) UseMethod("validity")

validity.copreg <- function(object, level = 0.05, power = 0.8, ...) {
  
  X  <- object$X
  y  <- object$y
  Cs <- object$copula.data
  n  <- nrow(X)
  nm <- object$endogenous
  K  <- length(nm)
  
  ic <- if (object$has.intercept) 1L else integer(0)
  wc <- setdiff(seq_len(ncol(X)), c(ic, match(nm, colnames(X))))
  wc <- wc[vapply(wc, function(j) length(unique(X[, j])) > 1L, logical(1))]
  W  <- X[, wc, drop = FALSE]
  
  ## --- [1] nonnormality -----------------------------------------------------
  ## Which variable has to be non-normal depends on the estimator.  For most of
  ## them it is the endogenous regressor.  For BMW it is the first-stage error:
  ## their Theorem 2.1 identifies the model if and only if the distribution of
  ## e is not normal, and the paper says how to check it -- "in practice one
  ## could use the first-stage residuals ehat to test whether e is Gaussian,
  ## which would lead to nonidentification".
  ## Which condition identifies the model is the estimator's own; for
  ## 2sCOPE-np it is Assumption 3 rather than nonnormality of P, because the
  ## conditional CDF, not the marginal one, generates the copula term.
  a3 <- identical(object$id.condition, "assumption3")
  on_resid <- !a3 && !is.null(object$first.stage.residuals)
  step1 <- NULL; nn1 <- list(thresholds = NULL)
  if (!a3) {
    M1    <- if (on_resid) object$first.stage.residuals else X
    nn1   <- .validity_nonnormality(M1, nm, n, power)
    step1 <- nn1$table
  }
  
  ## --- [2] the uncorrelatedness assumption ----------------------------------
  ## Only meaningful for estimators that rely on it.  2sCOPE and the other
  ## two-stage estimators project the correlation out in the first stage, so
  ## the check would just restate what they already handle.
  step2 <- NULL
  if (ncol(W) > 0L && isTRUE(object$assumption5)) {
    fpg <- .lm.fit(cbind(X, Cs), y)
    gpg <- fpg$coefficients[ncol(X) + seq_len(K)]
    CTT <- as.vector(Cs %*% gpg)
    
    r  <- suppressWarnings(as.vector(stats::cor(W, CTT)))
    pa <- stats::p.adjust(vapply(r, .fisher_z_p, numeric(1), n = n),
                          method = "holm")
    fw  <- .lm.fit(cbind(1, W), CTT)
    df1 <- fw$rank - 1L; df2 <- n - fw$rank
    r2  <- 1 - sum(fw$residuals^2) / sum((CTT - mean(CTT))^2)
    Fs  <- (r2 / df1) / ((1 - r2) / df2)
    pj  <- stats::pf(Fs, df1, df2, lower.tail = FALSE)
    
    step2 <- list(
      table = data.frame(`corr(W, CTT)` = round(r, 4),
                         `p (Holm)` = format.pval(pa, digits = 3, eps = 1e-16),
                         row.names = colnames(W), check.names = FALSE,
                         stringsAsFactors = FALSE),
      joint.R2 = r2, joint.F = Fs, joint.p = pj,
      violated = any(pa < level, na.rm = TRUE) || (pj < level))
  }
  
  ## --- [3] exogenous regressors as identifying variation --------------------
  step3 <- NULL
  weak <- if (is.null(step1)) character(0) else nm[!step1$`Yang ok`]
  if (length(weak) > 0L && ncol(W) > 0L && !is.null(object$Wstar)) {
    
    Ws   <- object$Wstar
    cont <- vapply(seq_len(ncol(W)), function(j) length(unique(W[, j])) > 20L,
                   logical(1))
    kw   <- vapply(seq_len(ncol(W)), function(j) .ks_normal(W[, j])["p"],
                   numeric(1))
    
    Fm <- matrix(NA_real_, ncol(Ws), length(weak),
                 dimnames = list(colnames(Ws), paste0("F: ", weak)))
    D  <- cbind(1, Ws)
    XtXi <- tryCatch(chol2inv(chol(crossprod(D))), error = function(e) NULL)
    if (!is.null(XtXi))
      for (k in seq_along(weak)) {
        fk  <- .lm.fit(D, Cs[, match(weak[k], nm)])
        se  <- sqrt(sum(fk$residuals^2) / (n - fk$rank) * diag(XtXi))
        Fm[, k] <- (fk$coefficients / se)[-1]^2
      }
    
    step3 <- data.frame(continuous = cont,
                        `KS p` = format.pval(kw, digits = 3, eps = 1e-16),
                        row.names = colnames(W), check.names = FALSE,
                        stringsAsFactors = FALSE)
    step3 <- cbind(step3, as.data.frame(round(Fm, 3), check.names = FALSE))
    step3$qualifies <- cont & kw < 0.001 &
      apply(Fm, 1, function(z) any(z > 10, na.rm = TRUE))
  }
  
  ## --- [4] error term -------------------------------------------------------
  xi <- object$residuals
  a4 <- .ad_test(xi)
  step4 <- list(skewness = .skewness(xi), ex.kurtosis = .ex_kurtosis(xi),
                AD = unname(a4["A"]), AD.p = unname(a4["p"]),
                KS.p = unname(.ks_normal(xi)["p"]))
  
  ## --- [5] ICON -------------------------------------------------------------
  icon <- data.frame(
    `SE (corrected)`   = object$std.error[names(object$ols.std.error)],
    `SE (uncorrected)` = object$ols.std.error,
    ICON               = object$icon,
    row.names = names(object$ols.std.error), check.names = FALSE)
  
  structure(list(method = object$method, n = n, level = level, power = power,
                 intercept = object$has.intercept, endogenous = nm,
                 has.first.stage = !is.null(object$Wstar),
                 thresholds = nn1$thresholds, nonnormality.of = if (on_resid)
                   "first-stage residuals" else "endogenous regressors",
                 id.condition = object$id.condition,
                 assumption3 = if (a3) object$diagnostics$collinearity else NULL,
                 step1 = step1, step2 = step2, step3 = step3, step4 = step4,
                 icon = icon, icon.max = max(object$icon, na.rm = TRUE),
                 has.exog = ncol(W) > 0L),
            class = "copreg.validity")
}


print.copreg.validity <- function(x, digits = 4, ...) {
  
  cat("\nValidity check for ", x$method, "\n", sep = "")
  cat("n = ", x$n, ", intercept: ", if (x$intercept) "yes" else "no",
      ", target power ", 100 * x$power, "%\n", sep = "")
  cat("Sources: ", if (identical(x$id.condition, "assumption3"))
    "Hu, Qian & Xie (2025); Becker, Proksch & Ringle (2022)\n" else
      paste0("Becker, Proksch & Ringle (2022); Yang, Qian & Xie (2025);\n",
             "         Qian, Koschmann & Xie (2025)\n"), sep = "")
  
  b <- 0L
  nxt <- function() { b <<- b + 1L; paste0("[", b, "] ") }
  
  if (!is.null(x$assumption3)) {
    cat("\n", nxt(), "Assumption 3: the copula term is not a linear combination ",
        "of (1, P, W)\n", sep = "")
    cat("    omega = 1 - R-squared of a copula term on all other regressors; a\n",
        "    value near zero means the augmented regression is close to singular.\n",
        sep = "")
    print(format(x$assumption3, digits = digits))
    cat("    Hu, Qian & Xie (2025, Sec. 3.2.4) make the standard error inflation\n",
        "    the operational test of this assumption, because it is severe\n",
        "    collinearity that inflates it; see the last block.\n", sep = "")
  } else {
    
    cat("\n", nxt(), "Nonnormality of the ", x$nonnormality.of, "\n", sep = "")
    if (identical(x$nonnormality.of, "first-stage residuals"))
      cat("    Theorem 2.1 of Breitung, Mayer & Wied (2024) identifies the model if",
          "and only\n    if the first-stage error is not normal, so that is what is",
          "tested here.\n")
    print(format(x$step1, digits = digits))
    th <- x$thresholds
    cat("    Becker et al. at n = ", x$n, ": |skewness| >= ",
        if (is.infinite(th["skewness"])) "not attainable" else
          formatC(th["skewness"], digits = 4),
        ", or AD > ", formatC(th["AD"], digits = 5),
        ", or CvM > ", formatC(th["CvM"], digits = 4), "\n", sep = "")
    cat("    Yang et al.: KS p < .05\n")
    if (any(!x$step1$`Becker ok`) && any(x$step1$`Yang ok`))
      cat("    Note: the two criteria disagree here. Becker et al. report a\n",
          "    correlation of .66 between the AD statistic and the copula term's\n",
          "    t statistic against -.03 for KS, and recommend AD or CvM; Yang et\n",
          "    al. prefer KS because it is conservative.\n", sep = "")
  }
  
  ## the uncorrelatedness assumption, only where the estimator relies on it
  if (!is.null(x$step2)) {
    cat("\n", nxt(), "Assumption: correlation of the copula transformation ",
        "term\n    with the exogenous regressors\n", sep = "")
    print(format(x$step2$table, digits = digits))
    cat("    Joint test: R2 = ", formatC(x$step2$joint.R2, digits = digits),
        ", F = ", formatC(x$step2$joint.F, digits = digits),
        ", p = ", format.pval(x$step2$joint.p, digits = digits), "\n", sep = "")
    cat("    => ", if (x$step2$violated)
      "violated. Park & Gupta is inconsistent here (Haschka, 2025)." else
        "satisfied. Park & Gupta is consistent here and more efficient.",
      "\n", sep = "")
  }
  
  if (!is.null(x$step1) && any(!x$step1$`Yang ok`)) {
    cat("\n", nxt(), "Exogenous regressors as identifying variation\n",
        "    (continuous, KS p < .001, first-stage F > 10)\n", sep = "")
    if (is.null(x$step3)) {
      cat("    Not available: this estimator has no first stage.\n")
    } else {
      print(format(x$step3, digits = digits))
      q <- x$step3$qualifies
      cat("    => ", if (any(q, na.rm = TRUE))
        paste0("qualifying regressor(s): ",
               paste(rownames(x$step3)[which(q)], collapse = ", "), ".")
        else paste0("none qualifies. The conditions are conservative and not\n",
                    "       necessary; Yang et al. (2025) propose a bootstrap ",
                    "procedure to gauge\n       the finite-sample bias in this ",
                    "situation."), "\n", sep = "")
    }
  }
  
  s4 <- x$step4
  cat("\n", nxt(), "Error term: structural residuals xi\n", sep = "")
  cat("    skewness = ", formatC(s4$skewness, digits = digits),
      ", excess kurtosis = ", formatC(s4$ex.kurtosis, digits = digits),
      ", AD = ", formatC(s4$AD, digits = digits),
      " (p = ", format.pval(s4$AD.p, digits = 3), ")\n", sep = "")
  if (identical(x$id.condition, "assumption3")) {
    cat("    Assumption 1 of Hu, Qian & Xie (2025) requires the error E OR its\n",
        "    endogenous part O to be normal, not both. Their double robustness\n",
        "    property (Web Appendix C.3) places no assumption on the exogenous\n",
        "    part of the error, so E = O + epsilon may be nonnormal and\n",
        "    heteroscedastic; the residuals show E, not O, and their shape\n",
        "    therefore neither establishes nor rules out a violation.\n", sep = "")
  } else
    cat("    Becker et al. (2022) find that with a nonnormal error the approach\n",
        "    is no longer consistent in models with an intercept. Yang et al.\n",
        "    (2025) and Qian et al. (2025) do permit a nonnormal error, but only\n",
        "    under the decomposition xi = U + V into a normally distributed\n",
        "    endogenous part U, which carries the entire dependence with the\n",
        "    regressors, and an independent nonnormal V. The residuals show xi,\n",
        "    not U, so their skewness neither establishes nor rules out a\n",
        "    violation; whether that decomposition holds has to be argued from\n",
        "    the suspected sources of endogeneity.\n", sep = "")
  
  cat("\n", nxt(), "ICON: standard error inflation relative to uncorrected OLS\n",
      sep = "")
  print(format(x$icon, digits = digits))
  cat("    => largest ICON = ", formatC(x$icon.max, digits = digits),
      if (x$icon.max > 6)
        ", above the threshold of 6: weak\n       identification or a misspecified dependence model."
      else ", below the threshold of 6.", "\n\n", sep = "")
  
  invisible(x)
}
