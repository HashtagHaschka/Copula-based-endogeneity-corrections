## =============================================================================
##  CopRegJAMS.R
##
##  Generalised Gaussian copula framework of Liengaard, Becker, Bennedsen,
##  Heiler, Taylor & Ringle (2025), Journal of the Academy of Marketing Science
##  53, 279-299.
##
##  Requires copreg-core.R to be sourced first.
## =============================================================================

if (!exists(".copreg_fit"))
  stop("Source copreg-core.R before this file.")


## -----------------------------------------------------------------------------
##  Copula term constructor
## -----------------------------------------------------------------------------
##  Equation 17.  With d_P endogenous regressors P and d_W continuous exogenous
##  regressors W, all transformed to normal scores C(.) = Phi^{-1}(Fhat(.)),
##
##      C(P_i, W_i) = ( C(P_i)'  C(W_i)' ) Sigma^{-1} [ I_dP ; 0_(dW x dP) ],
##
##  where Sigma is the variance-covariance matrix of (C(P), C(W)).  Note that
##  this is the covariance matrix, not the correlation matrix: the two give
##  different linear combinations and therefore different estimates.  The result
##  is a d_P-vector of copula terms, one per endogenous regressor.
##
##  For a single P and a single W the paper writes it out:
##
##      C(P, W) = { C(P) s2_W - C(W) s_PW } / { s2_P s2_W - s2_PW },
##
##  which the unit test in the examples file checks against.
##
##  Equation 19 adds these to the structural model:
##
##      y_i = g(P_i, W_i) + sum_k gamma_k Chat_k(P_i, W_i) + u_i.
##
##  Equations 20 and 21 extend this to discrete exogenous regressors Z.  Z is
##  the categorical variable formed by ALL joint combinations of the discrete
##  regressors: two binary variables give four cells, not two.  Within each
##  cell the CDFs of P and W and the covariance matrix Sigma are estimated from
##  that cell's observations only, and the regression becomes
##
##      y_i = g(P_i, W_i, Z_i)
##            + sum_z sum_k gamma^z_k Chat^z_k(P_i, W_i) I(Z_i = z) + u_i,
##
##  so each copula column is zero outside its own cell.  If the copula
##  structure can be assumed constant across Z, conditional = FALSE collapses
##  this back to Equation 18: Z stays in g as a regressor but does not enter
##  the copula terms.
##
##  Which exogenous regressors count as discrete Z is controlled by the
##  conditional argument, see below.

## The layout -- which columns are discrete, which cells they form, and which
## continuous regressors are usable inside each cell -- is derived once from the
## original design matrix and then reused for every bootstrap resample.  Doing
## it per resample would renumber the cells by order of first appearance and
## silently permute the copula columns between replicates.

.jams_layout <- function(info, conditional) {
  
  X   <- info$X
  exo <- which(info$part >= 2L & info$ord == 1L)
  
  asg <- attr(X, "assign")
  lab <- character(ncol(X))
  lab[asg > 0L] <- attr(info$terms, "term.labels")[asg[asg > 0L]]
  
  ## --- three roles an exogenous regressor can play --------------------------
  ##  cells  the copula structure is estimated separately per category
  ##  disc   discrete, hence never copula transformed: W in Equation 17 is the
  ##         continuous exogenous vector
  ##  wall   continuous, copula transformed and entering Equation 17
  ##
  ##  Every exogenous regressor stays in the structural model g() under all
  ##  three settings of 'conditional'.  What changes is only whether it splits
  ##  the copula structure.
  
  fac <- intersect(exo, info$factor_cols)
  
  if (is.character(conditional)) {
    sel <- vapply(lab, function(l) nzchar(l) &&
                    any(all.vars(stats::reformulate(l)) %in% conditional),
                  logical(1))
    cells <- intersect(exo, which(sel))
    hit   <- unique(unlist(lapply(lab[cells],
                                  function(l) all.vars(stats::reformulate(l)))))
    miss  <- setdiff(conditional, hit)
    if (length(miss) > 0L)
      stop("Variable(s) named in 'conditional' are not exogenous first-order ",
           "regressors: ", paste(miss, collapse = ", "), call. = FALSE)
    ## a factor that is not conditioned on stays a control in g() but is still
    ## not continuous, so it stays out of the copula terms as well
    if (length(setdiff(fac, cells)) > 0L)
      message("Discrete regressor(s) not conditioned on and therefore left ",
              "out of the copula\n  terms as well, while staying in the model: ",
              paste(unique(lab[setdiff(fac, cells)]), collapse = ", "), ".")
  } else if (isTRUE(conditional)) {
    cells <- fac
    ## a numerically coded categorical regressor would otherwise pass as
    ## continuous, so say so rather than guess
    few <- setdiff(exo, fac)[vapply(setdiff(exo, fac), function(j)
      length(unique(X[, j])) <= 10L, logical(1))]
    if (length(few) > 0L)
      message("Treated as continuous although they take few distinct values: ",
              paste(colnames(X)[few], collapse = ", "),
              ".\n  Wrap them in as.factor() or name them in 'conditional' to ",
              "let the copula\n  structure vary over their categories.")
  } else {
    cells <- integer(0)
  }
  
  ## anything used for conditioning is discrete by construction, even if it was
  ## coded numerically
  disc <- union(fac, cells)
  
  ## --- cells ---------------------------------------------------------------
  if (length(cells) == 0L) {
    lev  <- NULL
    cell <- rep(1L, nrow(X))
  } else {
    lev  <- sort(unique(.copreg_key(X[, cells, drop = FALSE])))
    cell <- match(.copreg_key(X[, cells, drop = FALSE]), lev)
  }
  J <- if (is.null(lev)) 1L else length(lev)
  
  ## --- readable cell labels -------------------------------------------------
  ## Cells are found from the dummy pattern, because that is all a bootstrap
  ## resample carries, but the labels are read off the model frame of the
  ## original data, which still holds the factor levels.  Naming follows
  ## model.matrix: term label pasted to level, several terms joined by ":".
  clab <- if (J == 1L) "" else {
    dterms <- unique(lab[cells])
    ok <- all(vapply(dterms, function(tm) !is.null(info$mf[[tm]]), logical(1)))
    if (!ok) paste0("z", seq_len(J)) else
      vapply(seq_len(J), function(jj) {
        r <- which(cell == jj)[1L]
        paste(vapply(dterms, function(tm)
          paste0(tm, as.character(info$mf[[tm]][r])), character(1)),
          collapse = ":")
      }, character(1))
  }
  
  ## --- continuous exogenous regressors usable inside each cell -------------
  ## A regressor can be constant within a cell even though it varies overall,
  ## and its copula transform would then have zero variance and make Sigma
  ## singular.  Such columns are dropped from that cell only.
  wall <- setdiff(exo, disc)
  wc   <- vector("list", J)
  drop <- character(0)
  for (jj in seq_len(J)) {
    r <- which(cell == jj)
    keep <- wall[vapply(wall, function(k)
      length(unique(X[r, k])) > 1L, logical(1))]
    wc[[jj]] <- keep
    if (length(keep) < length(wall))
      drop <- c(drop, paste0(paste(colnames(X)[setdiff(wall, keep)],
                                   collapse = ", "),
                             if (J > 1L) paste0(" in ", clab[jj])))
  }
  if (length(drop) > 0L)
    message("Constant within a category and therefore left out of the copula ",
            "terms there,\n  while staying in the model: ",
            paste(drop, collapse = "; "), ".")
  
  ## --- is every cell large enough, now and under resampling? ---------------
  ## Each cell has to carry dP + dW marginal CDFs and their covariance matrix,
  ## so it needs more observations than that.  Checking it here rather than
  ## inside the constructor means the message names the category and appears
  ## before the bootstrap rather than after several hundred failed draws.
  nn   <- nrow(X)
  dP   <- length(info$endo_cols)
  size <- vapply(seq_len(J), function(jj) sum(cell == jj), integer(1))
  need <- dP + vapply(wc, length, integer(1)) + 2L
  ctag <- if (J == 1L) "the sample" else clab
  
  bad <- which(size < need)
  if (length(bad) > 0L)
    stop("Too few observations to estimate the copula structure separately ",
         "for:\n  ",
         paste(sprintf("%s: %d observations, %d needed for %d marginal CDFs ",
                       ctag[bad], size[bad], need[bad], need[bad] - 2L),
               collapse = "\n  "),
         "\n  Name fewer variables in 'conditional', merge categories, or use ",
         "conditional = FALSE.", call. = FALSE)
  
  if (J > 1L) {
    ## probability that a resample leaves a cell with fewer rows than it needs
    pf <- stats::pbinom(need - 1L, nn, size / nn)
    ok <- prod(1 - pf)
    if (ok < 0.95)
      warning(sprintf(paste0("Small categories: a bootstrap resample keeps ",
                             "every cell usable only %.0f%% of the time, so ",
                             "roughly\n  %.1f draws are needed per replicate. ",
                             "Tightest: %s with %d observations.\n  Name fewer ",
                             "variables in 'conditional' or merge categories."),
                      100 * ok, 1 / max(ok, 1e-6), ctag[which.max(pf)],
                      size[which.max(pf)]), call. = FALSE)
  }
  
  list(cells = cells, lev = lev, J = J, wc = wc, clab = clab, need = need)
}


.ctor_jams <- function(conditional = TRUE) {
  
  layout <- NULL
  
  function(X, info, cdf, ties) {
    
    if (is.null(layout)) layout <<- .jams_layout(info, conditional)
    n  <- nrow(X)
    P  <- X[, info$endo_cols, drop = FALSE]
    dP <- ncol(P)
    J  <- layout$J
    
    cell <- if (is.null(layout$lev)) rep(1L, n)
    else match(.copreg_key(X[, layout$cells, drop = FALSE]), layout$lev)
    
    Cstar <- matrix(NA_real_, n, dP, dimnames = list(NULL, info$endo_names))
    C     <- matrix(0, n, dP * J)
    Wstar <- NULL
    col   <- 0L
    nms   <- character(dP * J)
    gen   <- character(dP * J)
    gce   <- character(dP * J)
    for (jj in seq_len(J)) {
      
      r  <- which(cell == jj)
      wj <- layout$wc[[jj]]
      ## In a resample a category can end up with too few rows even though it
      ## was large enough in the data.  The draw is then unusable and the
      ## bootstrap redraws; .jams_layout has already reported how often that is
      ## to be expected.
      if (length(r) < layout$need[jj])
        stop("Category ", if (J == 1L) "" else layout$clab[jj],
             " holds only ", length(r), " observations in this resample, ",
             layout$need[jj], " needed.", call. = FALSE)
      
      ## Step 1: CDFs within the cell
      cp <- .copula_transform_matrix(P[r, , drop = FALSE], cdf, ties)
      Cstar[r, ] <- cp
      if (length(wj) > 0L) {
        cw <- .copula_transform_matrix(X[r, wj, drop = FALSE], cdf, ties)
        ## only meaningful with a single cell, see the note at the end
        if (J == 1L) Wstar <- cw
        M <- cbind(cp, cw)
      } else {
        M <- cp
      }
      
      ## Step 2: Equation 21 with the empirical variance-covariance matrix
      S  <- stats::cov(M)
      Si <- tryCatch(solve(S), error = function(e) NULL)
      if (is.null(Si))
        stop("The covariance matrix of the copula data is singular in category ",
             if (J == 1L) "the sample" else layout$clab[jj],
             ": two of the regressors are collinear there.", call. = FALSE)
      Ck <- M %*% Si[, seq_len(dP), drop = FALSE]
      
      for (k in seq_len(dP)) {
        col <- col + 1L
        C[r, col] <- Ck[, k]
        nms[col] <- if (J == 1L) paste0(info$endo_names[k], "_cop")
        else paste0(info$endo_names[k], "_cop.", layout$clab[jj])
        gen[col] <- info$endo_names[k]
        gce[col] <- if (J == 1L) "" else layout$clab[jj]
      }
    }
    colnames(C) <- nms
    grp <- data.frame(endogenous = gen, cell = gce, stringsAsFactors = FALSE)
    
    ## Wstar is only a well-defined single first-stage design when there is one
    ## cell; with several, each cell has its own transformation and its own set
    ## of usable regressors, so validity() reports no first stage rather than
    ## splicing incomparable pieces together.
    list(C = C, Cstar = Cstar, Wstar = if (J == 1L) Wstar else NULL,
         groups = grp)
  }
}


## -----------------------------------------------------------------------------
##  CopRegJAMS
## -----------------------------------------------------------------------------
##
##  formula      y ~ endog | exog
##
##               Everything lm() accepts is allowed, including the interactions
##               and non-linear transformations that Equations 13 to 16 cover.
##               Copula terms are generated for the first-order endogenous
##               terms only, matching d_P in Equation 17.
##
##  conditional  TRUE (default)     the copula structure is estimated
##                                   separately per joint category of the
##                                   factor and character regressors in the
##                                   exogenous part, Equations 20 and 21
##               FALSE               one common structure, Equation 18
##               character vector    names the variables whose categories the
##                                   structure may vary over; useful when the
##                                   joint number of cells would be too large
##                                   for the sample
##
##               Under all three settings the discrete regressors stay in the
##               model as ordinary controls.  They are never copula
##               transformed: W in Equation 17 is the continuous exogenous
##               vector.
##
##  cdf          default "ecdf.adj", the adjusted ECDF this paper proposes.
##
##  summary() reports the bootstrap Wald test of the copula terms' joint
##  significance and, with several cells, the test of whether the copula
##  structure varies across them.

CopRegJAMS <- function(formula, data,
                       conditional = TRUE,
                       cdf = "ecdf.adj",
                       ties = "max",
                       nboots = 199,
                       subset = NULL,
                       contrasts = NULL,
                       parallel = FALSE,
                       ncores = NULL,
                       verbose = interactive()) {
  
  cl <- match.call()
  
  .copreg_fit(formula = formula, data = data,
              ctor = .ctor_jams(conditional = conditional),
              method = paste0("JAMS (Liengaard et al. 2025)",
                              if (!isTRUE(conditional)) ", unconditional"),
              cdf = cdf, ties = ties, nboots = nboots,
              subset = subset, contrasts = contrasts,
              parallel = parallel, ncores = ncores, assumption5 = FALSE,
              verbose = verbose, call = cl)
}


## register for the generic entry point
.copreg_registry[["jams"]] <- list(
  fun   = "CopRegJAMS",
  label = "JAMS (Liengaard et al. 2025)",
  cdf   = "ecdf.adj")


### REFERENCES ----------------------------------------------------------------
##
## Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and
##   C. M. Ringle (2025). Dealing with regression models' endogeneity by means
##   of an adjusted estimator for the Gaussian copula approach. Journal of the
##   Academy of Marketing Science 53, 279-299.
##
## Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
##   two-stage copula generated regressor approach. Journal of Marketing
##   Research 62(4), 601-623.