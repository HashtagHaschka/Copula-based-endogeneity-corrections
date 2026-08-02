## =============================================================================
##  CopRegIMA.R
##
##  Copula correction exploiting between-regressor correlation, Haschka (2025),
##  IMA Journal of Management Mathematics 36(1), 161-180.
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
##    Haschka, R. E. (2025). Robustness of copula-correction models in causal
##      analysis: Exploiting between-regressor correlation. IMA Journal of
##      Management Mathematics 36(1), 161-180.
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
##  Procedure of Section 3.2:
##
##    1. transform every explanatory variable, endogenous and exogenous, to
##       normal scores: W* = Phi^{-1}(Fhat_W(W)) for W in {X_1..X_K, P_1..P_L}
##    2. regress each endogenous P*_l on X*_1, ..., X*_K and keep the residual
##    3. add those residuals to the structural regression
##
##  This is the same procedure as 2sCOPE with one difference: the auxiliary
##  regression in step 2 carries no intercept.  Haschka motivates it that way
##  because the slope of P* on X* is then the Pearson correlation r between the
##  two, which is what makes the decomposition
##
##      (P* - r X*) / sqrt(1 - r^2)
##
##  in his derivation work and allows rho to be recovered from it.  2sCOPE
##  instead includes the intercept, following note 8 of Qian, Koschmann & Xie
##  (2025).  Since the normal scores are close to mean zero, the two usually
##  differ only slightly; the estimators are otherwise identical.
##
##  The scaling by 1/sqrt(1 - r^2) is not applied.  It rescales the generated
##  regressor and therefore only rescales gamma, leaving alpha and beta
##  untouched.  The endogeneity measure rho of Section 3.1.3 is reported
##  separately in summary() as corr(xi*, P*).

## A fresh closure per call: .ctor_twostage() caches which exogenous columns
## survive the first stage, and that choice belongs to one model, not to the
## session.
.ctor_ima <- function() .ctor_twostage(intercept = FALSE)


## -----------------------------------------------------------------------------
##  CopRegIMA
## -----------------------------------------------------------------------------
##
##  formula   y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##            Everything lm() accepts is allowed.  Endogenous regressors must
##            be numeric.  Without an exogenous part the first stage is empty
##            and the call is redirected to CopRegPG() with a warning.
##
##  cdf       Haschka (2025) leaves the CDF estimator open ("Fhat_W is an
##            estimate of the cdf of W"), so the default follows 2sCOPE:
##            "rank.n", the algorithm of Equation 9 in Qian, Koschmann & Xie
##            (2025).  Any other choice from copreg-core.R is accepted.
##
##  ties      "max" (default) or "average".  Relevant here because the
##            exogenous regressors are transformed as well.
##
##  Use validity() on the result for the identification checks.

CopRegIMA <- function(formula, data,
                      cdf = "rank.n",
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
    warning("IMA without exogenous regressors is identical to Park & Gupta ",
            "(the first stage has nothing to project on). Redirecting to ",
            "CopRegPG() with cdf = \"", pgcdf, "\".", call. = FALSE)
    out <- CopRegPG(formula = formula, data = data, cdf = pgcdf, ties = ties,
                    nboots = nboots, subset = subset, contrasts = contrasts,
                    parallel = parallel, ncores = ncores, verbose = verbose)
    out$call <- cl
    return(out)
  }
  
  .copreg_fit(formula = formula, data = data,
              ctor = .ctor_ima(),
              method = "IMA (Haschka 2025)",
              cdf = cdf, ties = ties, nboots = nboots,
              subset = subset, contrasts = contrasts,
              parallel = parallel, ncores = ncores, assumption5 = FALSE,
              verbose = verbose, call = cl)
}


## register for the generic entry point
.copreg_registry[["ima"]] <- list(
  fun   = "CopRegIMA",
  label = "IMA (Haschka 2025)",
  cdf   = "rank.n")


### REFERENCES ----------------------------------------------------------------
##
## Haschka, R. E. (2025). Robustness of copula-correction models in causal
##   analysis: Exploiting between-regressor correlation. IMA Journal of
##   Management Mathematics 36(1), 161-180.
##
## Qian, Y., A. Koschmann, and H. Xie (2025). A practical guide to endogeneity
##   correction using copulas. Journal of Marketing.
##
## Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
##   two-stage copula generated regressor approach. Journal of Marketing
##   Research 62(4), 601-623.
