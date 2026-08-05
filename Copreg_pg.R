## =============================================================================
##  CopRegPG.R
##
##  Park & Gupta (2012) copula endogeneity correction, least-squares control
##  function representation.
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
##    Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##      estimation using copulas. Marketing Science 31(4), 567-586.
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##  Copula term constructor
## -----------------------------------------------------------------------------
##  Park & Gupta add the copula transformation of each endogenous regressor
##  directly to the regression:
##
##      y = mu + P alpha + W beta + sum_k gamma_k * Phi^{-1}(Fhat(P_k)) + u
##
##  Exogenous regressors are not transformed and do not enter the copula term.
##  The control function and the plain copula data therefore coincide.

.ctor_pg <- function(X, info, cdf, ties) {
  P <- X[, info$endo_cols, drop = FALSE]
  C <- .copula_transform_matrix(P, cdf, ties)
  colnames(C) <- paste0(info$endo_names, "_cop")
  list(C = C, Cstar = C)
}


## -----------------------------------------------------------------------------
##  CopRegPG
## -----------------------------------------------------------------------------
##
##  formula   y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##            Everything lm() accepts is allowed: log(p), I(p^2), p*m,
##            as.factor(g), "-1" to drop the intercept.  The exogenous part may
##            be omitted entirely.  Endogenous regressors must be numeric.
##
##  data      a data.frame
##
##  cdf       CDF estimator for the copula transformation, see the table at the
##            end of copreg-core.R.  Default "kde.silverman", the estimator
##            Park & Gupta use in the original article.
##
##  ties      "max" (default) or "average", see copreg-core.R.
##
##  nboots    number of bootstrap replicates for the standard errors.
##
##  parallel  FALSE (default) or TRUE.  With TRUE the bootstrap is spread over
##            ncores workers: forked processes on Unix and macOS, a PSOCK
##            cluster on Windows.  Each worker gets its own RNG stream, so the
##            replicates stay independent.  No seed is set either way; set one
##            yourself beforehand if you want reproducible standard errors.
##
##  ncores    number of workers; defaults to detectCores(logical = FALSE) - 1.
##
##  Returns an object of class "copreg" with methods coef, vcov, confint,
##  residuals, fitted, predict, print and summary.
##
##  Note that no seed is set internally.  Call set.seed() beforehand if you
##  need reproducible standard errors.

CopRegPG <- function(formula, data,
                     cdf = "kde.silverman",
                     ties = "max",
                     nboots = 199,
                     subset = NULL,
                     contrasts = NULL,
                     parallel = FALSE,
                     ncores = NULL,
                     verbose = interactive()) {
  
  cl <- match.call()
  
  .copreg_fit(formula = formula, data = data,
              ctor = .ctor_pg, method = "PG (Park & Gupta 2012)",
              cdf = cdf, ties = ties, nboots = nboots,
              subset = subset, contrasts = contrasts,
              parallel = parallel, ncores = ncores,
              verbose = verbose, call = cl)
}


## -----------------------------------------------------------------------------
##  Generic entry point
## -----------------------------------------------------------------------------
##
##  copreg() only rewrites the call and hands it to the estimator's own
##  function.  That matters: each wrapper carries rules of its own -- BMW warns
##  when the CDF departs from the one its asymptotics are derived for, the
##  two-stage estimators redirect to Park & Gupta when the formula has no
##  exogenous part, and each tells summary() which assumptions it makes.  An
##  entry point that called the fitting engine directly would have to repeat
##  all of that and would drift out of step the moment one of them changes.
##
##  cdf = NULL, the default, is dropped from the call so that the estimator's
##  own default applies.

.copreg_registry <- list(
  pg = list(fun = "CopRegPG", label = "PG (Park & Gupta 2012)",
            cdf = "kde.silverman")
  ## the other estimators register themselves in their own files
)

copreg <- function(formula, data,
                   method = c("pg", "2scope", "ima", "bmw", "jams", "np"),
                   cdf = NULL, ...) {
  
  method <- match.arg(method)
  reg <- .copreg_registry[[method]]
  if (is.null(reg))
    stop("Estimator \"", method, "\" is not implemented yet.", call. = FALSE)
  
  cl <- match.call()
  cl[[1L]] <- as.name(reg$fun)
  cl$method <- NULL
  if (is.null(cdf)) cl$cdf <- NULL
  
  out <- eval(cl, parent.frame())
  out$call <- match.call()
  out
}


### REFERENCES ----------------------------------------------------------------
##
## Becker, J.-M., D. Proksch, and C. M. Ringle (2022). Revisiting Gaussian
##   copulas to handle endogenous regressors. Journal of the Academy of
##   Marketing Science 50, 46-66.
##
## Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of
##   endogeneity corrections using nonlinear transformations. The Econometrics
##   Journal 27(3), 362-383.
##
## Li, Q., J. Li, and J. S. Racine (2017). Cross-validated mixed-datatype
##   bandwidth selection for nonparametric cumulative distribution/survivor
##   functions. Econometric Reviews 36, 970-987.
##
## Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and
##   C. M. Ringle (2025). Dealing with regression models' endogeneity by means
##   of an adjusted estimator for the Gaussian copula approach. Journal of the
##   Academy of Marketing Science 53, 279-299.
##
## Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##   estimation using copulas. Marketing Science 31(4), 567-586.
##
## Qian, Y., A. Koschmann, and H. Xie (2025). A practical guide to endogeneity
##   correction using copulas. Journal of Marketing.
##
## Qian, Y. and H. Xie (2024). Simplifying bias correction for selective
##   sampling: A unified distribution-free approach.
##
## Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
##   two-stage copula generated regressor approach. Journal of Marketing
##   Research 62(4), 601-623.
