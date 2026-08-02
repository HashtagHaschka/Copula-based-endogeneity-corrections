## =============================================================================
##  CopRegBMW.R
##
##  Nonlinear endogeneity correction of Breitung, Mayer & Wied (2024),
##  The Econometrics Journal 27(3), 362-383.
##
##  Requires copreg-core.R to be sourced first.
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
##    Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of
##      endogeneity corrections using nonlinear transformations. The
##      Econometrics Journal 27(3), 362-383.
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
##  Copula term constructor
## -----------------------------------------------------------------------------
##  Model, Equations 2.1 and 2.2:
##
##      y = x'beta + z gamma + u,     z = x'delta + e,     u = rho f(e) + eps,
##
##  with x the exogenous regressors including a constant, z endogenous,
##  e independent of x and E[eps | x, z] = 0.  The function f is strictly
##  monotone with f(e) ~ N(0, 1), which gives f(e) = Phi^-1(F_e(e)).
##
##  Theorem 2.1: beta, gamma and rho are identified if and only if the
##  distribution of e is not normal.  Note that this is a condition on the
##  first-stage error, not on z itself.
##
##  Estimator, Equation 2.3:
##
##    1. regress z on x by OLS and keep the residuals ehat = z - x'deltahat
##    2. etahat = Phi^-1( Fhat(ehat) ),  Fhat(ehat_i) = rank(ehat_i) / (n + 1)
##    3. regress y on x, z and etahat by OLS
##
##  The order of the two operations is what separates this from 2sCOPE, and the
##  authors put the difference plainly in their abstract: they "remove the
##  (linear) dependence among the exogenous and endogenous regressors by
##  computing the residuals from an auxiliary regression" and only then apply
##  the normal quantile transform to the ranks of those residuals, whereas
##  2sCOPE "calculate residuals from a regression of the transformed
##  regressors, thereby assuming a particular nonlinear dependence among the
##  exogenous and endogenous regressors".  Same two ingredients, opposite order,
##  and the two are not reparameterisations of one another.
##
##  Several endogenous regressors, Remark 2.1: each z_j is regressed on x on its
##  own, never on the other endogenous regressors, and every residual gets its
##  own transform and its own rho_j.
##
##  The first stage always carries an intercept, whether or not the structural
##  model does, because Assumption A4 requires E[e] = 0.

.ctor_bmw <- function() {
  
  cols <- NULL
  
  function(X, info, cdf, ties) {
    
    if (is.null(cols)) cols <<- .copreg_exog_cols(info, X)
    wc <- cols
    
    P  <- X[, info$endo_cols, drop = FALSE]
    dP <- ncol(P)
    
    ## Step 1: first-stage residuals on the untransformed variables
    if (length(wc) == 0L) {
      ## Without exogenous regressors delta is empty and the residual is just z
      ## centred, so the estimator coincides with Park & Gupta computed on
      ## rank/(n+1).  CopRegBMW() redirects before reaching this point; the
      ## branch only guards the bootstrap.
      E <- sweep(P, 2L, colMeans(P), "-")
    } else {
      D <- cbind(`(Intercept)` = 1, X[, wc, drop = FALSE])
      qd <- qr(D)
      if (qd$rank < ncol(D))
        stop("The first stage of BMW is rank deficient.", call. = FALSE)
      E <- matrix(NA_real_, nrow(X), dP)
      for (k in seq_len(dP)) E[, k] <- .lm.fit(D, P[, k])$residuals
    }
    colnames(E) <- info$endo_names
    
    ## Step 2: normal scores of the first-stage residuals
    C <- .copula_transform_matrix(E, cdf, ties)
    colnames(C) <- paste0(info$endo_names, "_cop")
    
    list(C = C, Cstar = C, Wstar = if (length(wc) == 0L) NULL
         else X[, wc, drop = FALSE],
         resid1 = E)
  }
}


## -----------------------------------------------------------------------------
##  CopRegBMW
## -----------------------------------------------------------------------------
##
##  formula   y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##            Everything lm() accepts is allowed.  Endogenous regressors must be
##            numeric.  Without an exogenous part the first stage is empty and
##            the call is redirected to CopRegPG() with a warning.
##
##            Assumption A4 has x entering z linearly.  The authors name the
##            remedy themselves: add nonlinear transformations of the exogenous
##            regressors, so w1:w2 or I(w^2) in the exogenous part go into the
##            first stage as they stand.
##
##  cdf       "rank.n1", Rank/(n+1), the transformation of Equation 2.3.  Unlike
##            the other estimators this one is not a free choice: the asymptotic
##            theory of Proposition 3.1 is derived for this normalisation, so
##            any other value produces a warning.
##
##  ties      "max" (default) or "average", see copreg-core.R.
##
##  Two results of the paper show up in summary().  Corollary 3.2 states that
##  for testing rho = 0 the textbook t statistic with classical OLS standard
##  errors keeps a standard normal limit, even though those standard errors are
##  wrong for the structural coefficients; that Durbin-Hausman-Wu test is
##  reported next to the bootstrap one.  Remark 2.2 rewrites the estimator as a
##  just-identified instrumental variable estimator whose instrument is the
##  residual of z on etahat, which is what the omega diagnostic measures.
##
##  Use validity() for the identification checks.  For this estimator the
##  non-normality requirement applies to the first-stage residuals rather than
##  to the endogenous regressor, and validity() tests them accordingly.

CopRegBMW <- function(formula, data,
                      cdf = "rank.n1",
                      ties = "max",
                      nboots = 199,
                      subset = NULL,
                      contrasts = NULL,
                      parallel = FALSE,
                      ncores = NULL,
                      verbose = interactive()) {
  
  cl <- match.call()
  
  F <- Formula::as.Formula(formula)
  if (length(F)[2L] < 2L) {
    pgcdf <- if (missing(cdf)) formals(CopRegPG)$cdf else cdf
    warning("BMW without exogenous regressors is identical to Park & Gupta ",
            "(the first stage has nothing to regress on). Redirecting to ",
            "CopRegPG() with cdf = \"", pgcdf, "\".", call. = FALSE)
    out <- CopRegPG(formula = formula, data = data, cdf = pgcdf, ties = ties,
                    nboots = nboots, subset = subset, contrasts = contrasts,
                    parallel = parallel, ncores = ncores, verbose = verbose)
    out$call <- cl
    return(out)
  }
  
  if (!missing(cdf) && !identical(cdf, "rank.n1"))
    warning("Proposition 3.1 of Breitung, Mayer & Wied (2024) is derived for ",
            "the rank transformation\n  of their Equation 2.3, cdf = ",
            "\"rank.n1\". With cdf = \"", cdf, "\" the point estimates remain ",
            "sensible\n  but the reported standard errors are no longer ",
            "covered by the published theory.", call. = FALSE)
  
  .copreg_fit(formula = formula, data = data,
              ctor = .ctor_bmw(),
              method = "BMW (Breitung, Mayer & Wied 2024)",
              cdf = cdf, ties = ties, nboots = nboots,
              subset = subset, contrasts = contrasts,
              parallel = parallel, ncores = ncores, assumption5 = FALSE,
              dhw = TRUE, verbose = verbose, call = cl)
}


## register for the generic entry point
.copreg_registry[["bmw"]] <- list(
  fun   = "CopRegBMW",
  label = "BMW (Breitung, Mayer & Wied 2024)",
  cdf   = "rank.n1")


### REFERENCES ----------------------------------------------------------------
##
## Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of
##   endogeneity corrections using nonlinear transformations. The Econometrics
##   Journal 27(3), 362-383.
##
## Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##   estimation using copulas. Marketing Science 31(4), 567-586.
##
## Qian, Y., A. Koschmann, and H. Xie (2025). A practical guide to endogeneity
##   correction using copulas. Journal of Marketing.
##
## Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
##   two-stage copula generated regressor approach. Journal of Marketing
##   Research 62(4), 601-623.
##
## Zhao, Y., I. Gijbels, and I. Van Keilegom (2020). Inference for
##   semiparametric Gaussian copula model adjusted for linear regression using
##   residual ranks. Bernoulli 26(4), 2815-2846.
