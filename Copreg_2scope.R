## =============================================================================
##  CopReg2sCOPE.R
##
##  Two-stage copula generated regressor approach (2sCOPE) of
##  Yang, Qian & Xie (2025), Journal of Marketing Research 62(4), 601-623.
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
##  and the paper this estimator implements
##
##    Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
##      two-stage copula generated regressor approach. Journal of Marketing
##      Research 62(4), 601-623.
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##  Copula term constructor
## -----------------------------------------------------------------------------
##  Estimation procedure, Table 4 of Yang, Qian & Xie (2025):
##
##  Stage 1
##    - obtain the empirical CDFs Hhat(P) and Lhat(W) of every regressor
##    - compute P* = Phi^{-1}(Hhat(P))  and  W* = Phi^{-1}(Lhat(W))
##    - regress each endogenous P*_k separately on W* and keep the residual
##
##  Stage 2
##    - add those residuals to the structural regression as generated regressors
##
##      y_i = mu + sum_k alpha_k P_ik + beta' W_i + sum_k gamma_k C_ik + eps_i,
##      C_ik = P*_ik - delta_k' W*_i.
##
##  Three points where the implementation follows the paper rather than the
##  derivation:
##
##  1. The first stage carries an intercept.  Equation 11 has none, but note 8
##     of Qian, Koschmann & Xie (2025) states that the implementation of 2sCOPE
##     includes it because that is more general and performs well in their
##     simulations.  This is the one place where 2sCOPE and the IMA estimator
##     of Haschka (2025) differ.
##
##  2. Each P*_k is regressed on W* only, separately, never on the other
##     endogenous regressors.
##
##  3. All exogenous regressors are copula transformed, continuous or not
##     (Web Appendices E.2 and E.3 of Yang, Qian & Xie 2025).  For a binary
##     dummy the transformation is affine, so the column space of the first
##     stage is unchanged and the residual is the same as with the raw dummy.

## A fresh closure per call: .ctor_twostage() caches which exogenous columns
## survive the first stage, and that choice belongs to one model, not to the
## session.
.ctor_2scope <- function() .ctor_twostage(intercept = TRUE)


## -----------------------------------------------------------------------------
##  CopReg2sCOPE
## -----------------------------------------------------------------------------
##
##  formula   y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##            Everything lm() accepts is allowed.  Endogenous regressors must
##            be numeric.  Without an exogenous part 2sCOPE is identical to
##            Park & Gupta and the call is redirected there with a warning.
##
##  cdf       default "rank.n", the algorithm of Equation 9 in
##            Qian, Koschmann & Xie (2025), which is what 2sCOPE uses.
##
##  ties      "max" (default) or "average", see copreg-core.R.  This matters
##            here in a way it does not for Park & Gupta, because the
##            exogenous regressors are transformed as well and dummies and
##            other discrete controls carry many ties.
##
##  Use validity() on the result for the decision tree of Yang, Qian & Xie
##  (2025, Figure 2) and the ICON statistic.

CopReg2sCOPE <- function(formula, data,
                         cdf = "rank.n",
                         ties = "max",
                         nboots = 199,
                         subset = NULL,
                         contrasts = NULL,
                         parallel = FALSE,
                         ncores = NULL,
                         verbose = interactive()) {
  
  cl <- match.call()
  
  ## without exogenous regressors the first stage is empty and 2sCOPE reduces
  ## to Park & Gupta; redirect rather than silently fitting a relabelled model
  F <- Formula::as.Formula(formula)
  if (length(F)[2L] < 2L) {
    pgcdf <- if (missing(cdf)) formals(CopRegPG)$cdf else cdf
    warning("2sCOPE without exogenous regressors is identical to Park & Gupta ",
            "(the first stage has nothing to project on). Redirecting to ",
            "CopRegPG() with cdf = \"", pgcdf, "\".", call. = FALSE)
    out <- CopRegPG(formula = formula, data = data, cdf = pgcdf, ties = ties,
                    nboots = nboots, subset = subset, contrasts = contrasts,
                    parallel = parallel, ncores = ncores, verbose = verbose)
    out$call <- cl
    return(out)
  }
  
  .copreg_fit(formula = formula, data = data,
              ctor = .ctor_2scope(),
              method = "2sCOPE (Yang, Qian & Xie 2025)",
              cdf = cdf, ties = ties, nboots = nboots,
              subset = subset, contrasts = contrasts,
              parallel = parallel, ncores = ncores, assumption5 = FALSE,
              verbose = verbose, call = cl)
}


## register for the generic entry point
.copreg_registry[["2scope"]] <- list(
  fun   = "CopReg2sCOPE",
  label = "2sCOPE (Yang, Qian & Xie 2025)",
  cdf   = "rank.n")


### REFERENCES ----------------------------------------------------------------
##
## Hu, Y., Y. Qian, and H. Xie (2025). Two-stage nonparametric copula control
##   function approach.
##
## Qian, Y., A. Koschmann, and H. Xie (2025). A practical guide to endogeneity
##   correction using copulas. Journal of Marketing.
##
## Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
##   two-stage copula generated regressor approach. Journal of Marketing
##   Research 62(4), 601-623.
