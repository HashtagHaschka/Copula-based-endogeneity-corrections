## =============================================================================
##  Copreg_2scope_np.R
##
##  Nonparametric two-stage copula control function approach (2sCOPE-np),
##  Hu, Qian & Xie (2025), NBER Working Paper 33607.
##
##  Requires copreg-core.R to be sourced first, and the np package.
## =============================================================================

if (!exists(".copreg_fit"))
  stop("Source copreg-core.R before this file.")


## -----------------------------------------------------------------------------
##  Copula term constructor
## -----------------------------------------------------------------------------
##  Table 1 of the paper, for each endogenous regressor P_k separately:
##
##    Stage 1  estimate the CONDITIONAL CDF of P_k given the exogenous
##             regressors W nonparametrically and set
##
##                 C_pk = Phi^{-1}( Fhat(P_k | W) ).
##
##    Stage 2  add C_p1, ..., C_pK to the structural regression as generated
##             regressors (Equation 16).
##
##  This is the one estimator in the toolbox that transforms nothing.  There is
##  no marginal CDF, no first stage in the regression sense and therefore no
##  first-stage residual (the paper says so explicitly in its note 31); the
##  whole correction sits in the conditional CDF.  Everything the other
##  estimators do -- Park & Gupta's marginal transform, the residual of P* on
##  W* in 2sCOPE and IMA, the rank transform of the first-stage residual in BMW
##  -- is what Fhat(P|W) collapses to once a particular auxiliary model for P
##  is imposed (Section 2.8).  Two consequences follow.
##
##  First, cdf and ties are meaningless here and the arguments are absent.
##
##  Second, rho is not identified.  Proposition 1 gives the control function as
##  E(o|p, w) = sigma_o * rho * C_p, so only the product enters the regression;
##  the paper reports the coefficient of the copula term and tests it for zero,
##  which is what summary() does.
##
##  The conditional CDF is estimated with the Li & Racine (2008) kernel
##  estimator of Equation 43, which smooths P as well as W, as implemented in
##  np::npcdist().  That smoothing is what makes the estimator work for
##  noncontinuous endogenous regressors: the conditional CDF takes continuous
##  values with no ties even when P is a count, as long as W is non-trivially
##  related to P.  Bandwidths come from the least-squares cross-validation of
##  Li, Lin & Racine (2013), np's own default.


## -----------------------------------------------------------------------------
##  Which variables to condition on
## -----------------------------------------------------------------------------
##  Taken from the model frame rather than the design matrix, so that factors
##  reach np as factors and are smoothed with the Aitchison & Aitken kernel
##  instead of arriving as numeric dummies.
##
##  Section 3.1.4: polynomial and interaction terms in W "should be removed
##  from the conditional CDF estimation because the nonparametric estimator
##  automatically accommodates such nonlinearities".  A model frame holds
##  variables rather than terms, so an interaction contributes nothing of its
##  own and drops out by itself; a polynomial such as I(w^2) does appear as a
##  column and is dropped here whenever the variable behind it is already
##  represented by a plain column.  Written differently -- log(price) with no
##  bare price anywhere in the model -- the transformed column is the only
##  representation of that variable and is kept.
##
##  Section 3.1.3 recommends dropping covariates that are irrelevant for P;
##  the 'condition' argument does that explicitly.  Everything stays in the
##  structural model either way.

.np_condition <- function(info, condition = NULL) {
  
  mf   <- info$mf
  resp <- attr(info$terms, "response")
  cols <- setdiff(seq_along(mf), resp)
  labs <- names(mf)[cols]
  vars <- lapply(labs, function(l) all.vars(stats::reformulate(l)))
  
  ## A variable counts as endogenous when it appears in an endogenous term and
  ## is not itself a first-order term of the exogenous part; the same rule as
  ## in .copreg_exog_cols().
  asg  <- attr(info$X, "assign")
  labX <- character(ncol(info$X))
  labX[asg > 0L] <- attr(info$terms, "term.labels")[asg[asg > 0L]]
  raw  <- unique(unlist(lapply(info$endo_names, function(l)
    all.vars(stats::reformulate(labX[match(l, colnames(info$X))])))))
  endo_vars <- setdiff(raw, unique(labX[info$part >= 2L & info$ord == 1L]))
  
  keep <- !vapply(vars, function(v) any(v %in% endo_vars), logical(1))
  
  ## a plain column is one whose label is the variable itself
  plain <- vapply(seq_along(labs), function(j)
    length(vars[[j]]) == 1L && identical(labs[j], vars[[j]]), logical(1))
  covered <- unlist(vars[keep & plain])
  drop <- keep & !plain &
    vapply(vars, function(v) length(v) > 0L && all(v %in% covered), logical(1))
  if (any(drop))
    message("Left out of the conditional CDF estimation because the variable ",
            "behind them is already\n  conditioned on and the kernel estimator ",
            "accommodates nonlinearity by itself\n  (Hu, Qian & Xie 2025, Sec. ",
            "3.1.4): ", paste(labs[drop], collapse = ", "),
            ".\n  They remain regressors in the structural model.")
  keep <- keep & !drop
  
  ## user-supplied restriction
  if (!is.null(condition)) {
    if (!is.character(condition))
      stop("'condition' must be a character vector of variable names.",
           call. = FALSE)
    hit  <- vapply(vars, function(v) any(v %in% condition), logical(1))
    miss <- setdiff(condition, unlist(vars[keep]))
    if (length(miss) > 0L)
      stop("Variable(s) named in 'condition' are not exogenous regressors of ",
           "the model: ", paste(miss, collapse = ", "), call. = FALSE)
    keep <- keep & hit
  }
  
  if (!any(keep))
    stop("No exogenous variable left to condition on. Without W the ",
         "conditional CDF is the\n  marginal one and the estimator is Park & ",
         "Gupta; use CopRegPG().", call. = FALSE)
  
  cols[keep]
}


## Columns of the model frame naming a grouping, checked against the model.
.np_group_cols <- function(info, vars, what) {
  if (is.null(vars)) return(integer(0))
  if (!is.character(vars))
    stop("'", what, "' must be a character vector of variable names.",
         call. = FALSE)
  mf   <- info$mf
  resp <- attr(info$terms, "response")
  cols <- setdiff(seq_along(mf), resp)
  hit  <- cols[vapply(cols, function(j)
    any(all.vars(stats::reformulate(names(mf)[j])) %in% vars), logical(1))]
  miss <- setdiff(vars, unlist(lapply(hit, function(j)
    all.vars(stats::reformulate(names(mf)[j])))))
  if (length(miss) > 0L)
    stop("Variable(s) named in '", what, "' are not regressors of the model: ",
         paste(miss, collapse = ", "), call. = FALSE)
  hit
}


## -----------------------------------------------------------------------------
##  The constructor
## -----------------------------------------------------------------------------
##  The layout -- which columns are conditioned on, which define cells, which
##  define the fixed effects -- is derived once from the original model frame
##  and reused for every resample, so that the copula columns keep their
##  meaning across bootstrap replicates.

.ctor_2scope_np <- function(condition = NULL, groups = NULL,
                            heterogeneous = NULL, demean = NULL,
                            bwmethod = "cv.ls", bw.refit = TRUE,
                            nboots = 0L, np.args = list(),
                            verbose = FALSE) {
  
  layout <- NULL      # conditioning columns, cells, fixed effects
  bw0    <- NULL      # full-sample bandwidths, reused when bw.refit = FALSE
  bwtab  <- NULL      # the same, formatted for summary()
  
  function(X, info, cdf, ties) {
    
    if (is.null(layout)) {
      
      cc <- .np_condition(info, condition)
      gc <- .np_group_cols(info, groups,        "groups")
      hc <- .np_group_cols(info, heterogeneous, "heterogeneous")
      dc <- .np_group_cols(info, demean,        "demean")
      
      ## Proposition 3 estimates the conditional CDF within each subsample, so
      ## naming heterogeneous groups implies estimation by those groups unless
      ## a coarser grouping was named explicitly.
      if (length(gc) == 0L) gc <- hc
      
      ## Variables that define cells or fixed effects are constant within their
      ## own cell and drop out of the conditioning set there.
      cc <- setdiff(cc, union(gc, dc))
      if (length(cc) == 0L)
        stop("No exogenous variable left to condition on after removing the ",
             "grouping variables.\n  Name fewer variables in 'groups', ",
             "'heterogeneous' or 'demean'.", call. = FALSE)
      
      mf <- info$mf
      
      ## Cells are numbered once, on the original data, and the resamples only
      ## index into that vector; nothing is re-matched per replicate, so the
      ## copula columns cannot be permuted between them.  interaction() is what
      ## fixes the order: the reference group of Proposition 3 is then the first
      ## factor level rather than whatever sorts first alphabetically.
      mkcell <- function(k) {
        if (length(k) == 0L) return(list(idx = rep(1L, nrow(X)), lab = NULL))
        f <- do.call(interaction,
                     c(unname(lapply(k, function(j) mf[[j]])),
                       list(drop = TRUE, lex.order = TRUE)))
        i <- as.integer(f)
        list(idx = i, lab = vapply(seq_len(nlevels(f)), function(j) {
          r <- which(i == j)[1L]
          paste(vapply(k, function(j2)
            paste0(names(mf)[j2], as.character(mf[[j2]][r])), character(1)),
            collapse = ":")
        }, character(1)))
      }
      cg   <- mkcell(gc)
      ch   <- mkcell(hc)
      cell <- cg$idx
      hcel <- ch$idx
      hlab <- ch$lab
      
      ## Proposition 3 estimates the conditional CDF within each heterogeneous
      ## subsample, so a coarser 'groups' would mix subsamples inside one CDF.
      if (length(gc) > 0L && length(hc) > 0L &&
          any(tapply(hcel, cell, function(z) length(unique(z))) > 1L))
        warning("The cells of 'groups' cut across those of 'heterogeneous', so ",
                "the conditional CDF\n  of a subsample is estimated together ",
                "with parts of other subsamples. Proposition 3\n  of Hu, Qian ",
                "& Xie (2025) estimates it within each subsample.",
                call. = FALSE)
      
      ## Are the cells large enough?  A conditional CDF needs more rows than
      ## conditioning variables, and the bootstrap has to keep them usable.
      need <- length(cc) + 2L
      size <- tabulate(cell, nbins = max(cell))
      if (any(size < need))
        stop("Too few observations to estimate the conditional CDF separately ",
             "in cell(s) of size ", paste(sort(size[size < need]),
                                          collapse = ", "),
             "; ", need, " are needed for ", length(cc),
             " conditioning variables.\n  Name fewer variables in 'groups' or ",
             "'heterogeneous', or merge categories.", call. = FALSE)
      
      layout <<- list(cond = cc, cell = cell, ncell = max(cell),
                      hcel = hcel, hlab = hlab, need = need,
                      dm = if (length(dc) == 0L) NULL else
                        .copreg_key(mf[, dc, drop = FALSE]),
                      cnames = names(mf)[cc])
    }
    
    idx <- attr(X, "idx")
    if (is.null(idx)) idx <- seq_len(nrow(X))
    
    n  <- nrow(X)
    K  <- length(info$endo_cols)
    D  <- info$mf[idx, layout$cond, drop = FALSE]
    P  <- X[, info$endo_cols, drop = FALSE]
    cl <- layout$cell[idx]
    
    ## Section 3.1.5, third route: sweep out group fixed effects from P and the
    ## continuous conditioning variables before estimating the conditional CDF.
    if (!is.null(layout$dm)) {
      g   <- layout$dm[idx]
      num <- vapply(D, is.numeric, logical(1))
      for (j in which(num))
        D[[j]] <- D[[j]] - stats::ave(D[[j]], g, FUN = mean)
      for (k in seq_len(K))
        P[, k] <- P[, k] - stats::ave(P[, k], g, FUN = mean)
    }
    
    Cp <- matrix(NA_real_, n, K, dimnames = list(NULL, info$endo_names))
    if (is.null(bw0)) bw0 <<- vector("list", layout$ncell * K)
    first <- is.null(attr(X, "idx"))
    
    op <- options(np.messages = FALSE); on.exit(options(op), add = TRUE)
    t0 <- proc.time()[["elapsed"]]
    
    for (j in seq_len(layout$ncell)) {
      r <- which(cl == j)
      if (length(r) < layout$need)
        stop("Cell ", j, " holds only ", length(r), " observations in this ",
             "resample, ", layout$need, " needed.", call. = FALSE)
      Dj <- D[r, , drop = FALSE]
      for (k in seq_len(K)) {
        Yj  <- data.frame(P[r, k]); names(Yj) <- info$endo_names[k]
        pos <- (j - 1L) * K + k
        bw  <- if (!bw.refit && !first && !is.null(bw0[[pos]]))
          bw0[[pos]] else
            do.call(np::npcdistbw,
                    c(list(xdat = Dj, ydat = Yj, bwmethod = bwmethod),
                      np.args))
        if (first) bw0[[pos]] <<- bw
        Cp[r, k] <- np::npcdist(bws = bw, txdat = Dj, tydat = Yj)$condist
      }
    }
    
    Cp <- stats::qnorm(Cp)
    if (anyNA(Cp) || any(!is.finite(Cp)))
      stop("The conditional CDF estimate hit 0 or 1, so the copula term is not ",
           "finite. This\n  usually means a cell of the conditioning variables ",
           "is too small for kernel\n  estimation.", call. = FALSE)
    
    if (first) {
      bwtab <<- do.call(rbind, lapply(seq_along(bw0), function(pos) {
        b <- bw0[[pos]]; if (is.null(b)) return(NULL)
        k <- (pos - 1L) %% K + 1L; j <- (pos - 1L) %/% K + 1L
        data.frame(regressor = info$endo_names[k],
                   cell = if (layout$ncell == 1L) "" else as.character(j),
                   variable = c(info$endo_names[k], layout$cnames),
                   bandwidth = c(b$ybw, b$xbw),
                   row.names = NULL, stringsAsFactors = FALSE)
      }))
      ## A condbandwidth object carries a copy of the data it was selected on.
      ## With bw.refit = TRUE nothing reads them again, and dropping them keeps
      ## the closure small when it is shipped to PSOCK workers.
      if (bw.refit) bw0 <<- list()
      if (verbose && nboots > 1L) {
        el <- proc.time()[["elapsed"]] - t0
        message(sprintf(paste0("Stage 1 took %.1f s. The bootstrap repeats it ",
                               "%d times, so expect roughly %s."),
                        el, nboots, .copreg_hms(el * nboots)))
      }
    }
    
    ## Homogeneous: one copula term per endogenous regressor.
    ## Proposition 3: the copula structure may differ across subsamples.  The
    ## paper's own parameterisation is used, group 1 as the reference, so that
    ## the remaining coefficients are the group differences and test them.
    if (is.null(layout$hlab) || length(layout$hlab) < 2L) {
      C <- Cp
      colnames(C) <- paste0(info$endo_names, "_cop")
      grp <- NULL
    } else {
      hcl <- layout$hcel[idx]
      J   <- length(layout$hlab)
      C   <- matrix(0, n, K * J)
      nms <- character(K * J); gen <- character(K * J); gce <- character(K * J)
      for (k in seq_len(K)) {
        C[, k]   <- Cp[, k]
        nms[k]   <- paste0(info$endo_names[k], "_cop")
        gen[k]   <- info$endo_names[k]
        gce[k]   <- ""
        for (j in 2L:J) {
          col <- K * (j - 1L) + k
          C[hcl == j, col] <- Cp[hcl == j, k]
          nms[col] <- paste0(info$endo_names[k], "_cop.", layout$hlab[j])
          gen[col] <- info$endo_names[k]
          gce[col] <- layout$hlab[j]
        }
      }
      colnames(C) <- nms
      grp <- data.frame(endogenous = gen, cell = gce, stringsAsFactors = FALSE)
      attr(grp, "reference") <- TRUE
    }
    
    list(C = C, Cstar = Cp, Wstar = NULL, groups = grp,
         bandwidths = if (first) bwtab else NULL)
  }
}


## -----------------------------------------------------------------------------
##  CopReg2sCOPEnp
## -----------------------------------------------------------------------------
##
##  formula        y ~ endog_1 + endog_2 + ... | exog_1 + exog_2 + ...
##                 Everything lm() accepts is allowed.  Endogenous regressors
##                 must be numeric, but unlike in every other copula method they
##                 may be discrete: the conditional CDF of Equation 43 smooths
##                 P, so counts and binary treatments are admissible (Section
##                 4.5 and Web Appendix C.6).  Without an exogenous part the
##                 conditional CDF is the marginal one and the estimator is
##                 Park & Gupta; the call is redirected there with a warning.
##
##  condition      character vector, which exogenous variables the conditional
##                 CDF conditions on.  The default uses all of them, minus the
##                 polynomial and interaction terms that Section 3.1.4 says to
##                 leave out.  Naming a subset implements the recommendation of
##                 Section 3.1.3 to drop covariates that are irrelevant for P.
##
##  groups         character vector.  Instead of smoothing the categorical
##                 variables with their own kernel, estimate the conditional
##                 CDF separately within each of their joint categories
##                 (Section 3.1.5, second route).  One copula coefficient as
##                 before.
##
##  heterogeneous  character vector, Proposition 3: the copula structure, and
##                 with it the endogeneity, may differ across the joint
##                 categories of these variables.  The conditional CDF is
##                 estimated within each subsample and the copula terms are
##                 parameterised against the first group, so that the further
##                 coefficients are deviations from it and test whether the
##                 endogeneity differs -- the parameterisation the paper uses
##                 in its return-to-education example.
##
##  demean         character vector, Section 3.1.5, third route: treat the
##                 categories of these variables as fixed effects and remove
##                 them by within-group demeaning of P and of the continuous
##                 conditioning variables before Stage 1.  For clustered or
##                 panel-like data whose categorical variables are too many to
##                 smooth and too fine to group by.
##
##  bwmethod       "cv.ls", the least-squares cross-validation of Li, Lin &
##                 Racine (2013) and np's default, which is what the paper
##                 uses; "normal-reference" is the cheap alternative.
##
##  bw.refit       TRUE, the procedure of Section 2.5: every bootstrap replicate
##                 runs the whole two-step procedure, bandwidth selection
##                 included.  This is expensive -- cross-validation is O(n^2)
##                 per evaluation -- and it is what makes the standard errors
##                 account for the variability of the first stage.  FALSE keeps
##                 the full-sample bandwidths and re-estimates only the
##                 conditional CDF, which is a deviation from the paper.
##
##  ...            passed to np::npcdistbw(), e.g. ckertype, ukertype, okertype,
##                 nmulti.  The defaults are np's own and hence the paper's:
##                 second-order Gaussian kernels for continuous variables and
##                 Aitchison & Aitken kernels for unordered factors.
##
##  Use validity() for the identification checks.  For this estimator that is
##  Assumption 3 rather than nonnormality of P, and the standard error inflation
##  relative to uncorrected OLS, which the paper reports for both of its
##  empirical applications.

CopReg2sCOPEnp <- function(formula, data,
                           condition = NULL,
                           groups = NULL,
                           heterogeneous = NULL,
                           demean = NULL,
                           bwmethod = c("cv.ls", "normal-reference"),
                           bw.refit = TRUE,
                           nboots = 199,
                           subset = NULL,
                           contrasts = NULL,
                           parallel = FALSE,
                           ncores = NULL,
                           verbose = interactive(), ...) {
  
  cl <- match.call()
  
  ## cdf and ties are deliberately not formal arguments: there is no marginal
  ## CDF to choose here.  Catching them in "..." keeps the error explicit while
  ## letting update(fit, method = ) drop them when switching away from an
  ## estimator that does have them.
  np.args <- list(...)
  if (any(c("cdf", "ties") %in% names(np.args)))
    stop("2sCOPE-np uses no marginal CDF: the copula term is built from the ",
         "conditional CDF\n  of P given W, estimated nonparametrically. The ",
         "'cdf' and 'ties' arguments\n  therefore do not apply here.",
         call. = FALSE)
  
  if (!requireNamespace("np", quietly = TRUE))
    stop("2sCOPE-np requires the 'np' package for the nonparametric ",
         "conditional CDF.\n  Install it with install.packages('np').",
         call. = FALSE)
  
  bwmethod <- match.arg(bwmethod)
  
  ## Without exogenous regressors Fhat(P|W) is Fhat(P) and the copula term is
  ## the one of Park & Gupta (Section 2.8).
  F <- Formula::as.Formula(formula)
  if (length(F)[2L] < 2L) {
    pgcdf <- formals(CopRegPG)$cdf
    warning("2sCOPE-np without exogenous regressors is identical to Park & ",
            "Gupta: with nothing to\n  condition on, the conditional CDF is ",
            "the marginal one. Redirecting to CopRegPG()\n  with cdf = \"",
            pgcdf, "\".", call. = FALSE)
    out <- CopRegPG(formula = formula, data = data, cdf = pgcdf,
                    nboots = nboots, subset = subset, contrasts = contrasts,
                    parallel = parallel, ncores = ncores, verbose = verbose)
    out$call <- cl
    return(out)
  }
  
  ## The tie warning of .copreg_model() describes the marginal transformation
  ## and does not apply here: handling noncontinuous endogenous regressors is
  ## what this estimator was built for.
  out <- withCallingHandlers(
    .copreg_fit(formula = formula, data = data,
                ctor = .ctor_2scope_np(condition = condition, groups = groups,
                                       heterogeneous = heterogeneous,
                                       demean = demean, bwmethod = bwmethod,
                                       bw.refit = bw.refit, nboots = nboots,
                                       np.args = np.args, verbose = verbose),
                method = "2sCOPE-np (Hu, Qian & Xie 2025)",
                cdf = NA_character_, ties = NA_character_, nboots = nboots,
                subset = subset, contrasts = contrasts,
                parallel = parallel, ncores = ncores, assumption5 = FALSE,
                id.condition = "assumption3", verbose = verbose, call = cl),
    warning = function(w) {
      if (grepl("substantial ties", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  out
}


## register for the generic entry point
.copreg_registry[["np"]] <- list(
  fun   = "CopReg2sCOPEnp",
  label = "2sCOPE-np (Hu, Qian & Xie 2025)",
  cdf   = NA_character_)


### REFERENCES ----------------------------------------------------------------
##
## Hayfield, T. and J. S. Racine (2008). Nonparametric econometrics: The np
##   package. Journal of Statistical Software 27(5), 1-32.
##
## Hu, X., Y. Qian, and H. Xie (2025). Correcting endogeneity via nonparametric
##   copula control functions. NBER Working Paper 33607.
##
## Li, Q. and J. S. Racine (2008). Nonparametric estimation of conditional CDF
##   and quantile functions with mixed categorical and continuous data. Journal
##   of Business & Economic Statistics 26(4), 423-434.
##
## Li, Q., J. Lin, and J. S. Racine (2013). Optimal bandwidth selection for
##   nonparametric conditional distribution and quantile functions. Journal of
##   Business & Economic Statistics 31(1), 57-65.
##
## Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
##   estimation using copulas. Marketing Science 31(4), 567-586.
##
## Trivedi, P. K. and D. M. Zimmer (2017). A note on identification of bivariate
##   copulas for discrete count data. Econometrics 5(1), 10.