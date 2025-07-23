# Load required packages
pacman::p_load(
  nlme,
  dplyr,
  Matrix,
  copula,
  ks,
  pbapply
)


# function for adj.ecdf
pobs1 <- function (x, na.last = "keep", ties.method = eval(formals(rank)$ties.method),
                   lower.tail = TRUE) {
  ties.method <- match.arg(ties.method)
  U <- if (!is.null(dim(x)))
    apply(x, 2, rank, na.last = na.last, ties.method = ties.method)*((nrow(x) - 1)/(nrow(x)^2)) + 1/(2*nrow(x))
  else rank(x, na.last = na.last, ties.method = ties.method)*((length(x) - 1)/(length(x)^2)) + 1/(2*length(x))
  if (inherits(x, "zoo"))
    attributes(U) <- attributes(x)
  if (lower.tail)
    U
  else 1 - U
}

# Bootstrap functions
boots1 <- function(formula, data_cleaned, dependent_var, independent_vars,
                   has_intercept, X, cdf) {

  data <- data_cleaned

  repeat {

    data_cleaned <- data %>%
      sample_n(size = nrow(data), replace = TRUE)
    design_matrix <- model.matrix(formula, data = data_cleaned)
    YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
    rank_A <- Matrix::rankMatrix(YP)[1]

    if (rank_A == ncol(design_matrix)) { break }

  }

  if (!has_intercept) { full_matrix <- model.matrix(~ . - 1, data = data_cleaned) } else {
    full_matrix <- model.matrix(~ ., data = data_cleaned) }

  # endogenous regressor(s)
  if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
  P1 <- as.matrix(P1)
  P_star1 <- matrix(NA, nrow = nrow(P1), ncol = ncol(P1))

  if (cdf == "kde") {

    for (i in 1:ncol(P1)) {
      Fhat <- ks::kcde(P1[, i])
      P_star1[, i] <- predict(Fhat, x = P1[, i])
    }

  } else if (cdf == "resc.ecdf") {

    P_star1 <- apply(P1, 2, copula::pobs)

  } else if (cdf == "adj.ecdf") {

    P_star1 <- apply(P1, 2, pobs1)

  } else if (cdf == "ecdf") {

    ecdf0 <- apply(P1, 2, ecdf)
    for (i in 1:ncol(P1)) {
      Fhat <- ecdf0[[i]]
      P_star1[, i] <- Fhat(P1[, i])
      P_star1[P_star1[, i] == min(P_star1[, i]), i] <- 10e-7
      P_star1[P_star1[, i] == max(P_star1[, i]), i] <- 1-10e-7
    }

  }

  # rename columns of P_star matrix
  colnames(P_star1) <- paste0(colnames(P1), "_cop")

  # merge columns
  est_matrix <- cbind(full_matrix, P_star1)

  # control function approach
  mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
             data = as.data.frame(est_matrix))

  Estimates <- mod1$coefficients
  return(Estimates)

}
boots_JAMS1 <- function(formula, data_cleaned, dependent_var, independent_P_vars,
                        independent_X_vars, has_intercept, X, cdf,
                        design_matrix1, full_formula) {

  data <- data_cleaned

  repeat {

    data_cleaned <- data %>%
      sample_n(size = nrow(data), replace = TRUE)
    design_matrix <- model.matrix(full_formula, data = data_cleaned)
    YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
    rank_A <- Matrix::rankMatrix(YP)[1]

    if (rank_A == ncol(design_matrix) &
        all(colnames(design_matrix1) %in% colnames(design_matrix))) { break }

  }

  full_matrix <- cbind(data_cleaned[[all.vars(full_formula)[1]]], design_matrix)
  colnames(full_matrix)[1] <- all.vars(full_formula)[1]

  if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
  P_star1 <- matrix(NA, nrow = nrow(P1), ncol = ncol(P1))

  if (cdf == "kde") {

    for (i in 1:ncol(P1)) {
      Fhat <- ks::kcde(P1[, i])
      P_star1[, i] <- predict(Fhat, x = P1[, i])
    }

  } else if (cdf == "resc.ecdf") {

    P_star1 <- apply(P1, 2, copula::pobs)

  } else if (cdf == "adj.ecdf") {

    P_star1 <- apply(P1, 2, pobs1)

  } else if (cdf == "ecdf") {

    ecdf0 <- apply(P1, 2, ecdf)
    for (i in 1:ncol(P1)) {
      Fhat <- ecdf0[[i]]
      P_star1[, i] <- Fhat(P1[, i])
      P_star1[P_star1[, i] == min(P_star1[, i]), i] <- 10e-7
      P_star1[P_star1[, i] == max(P_star1[, i]), i] <- 1-10e-7
    }

  }

  # rename columns of P_star matrix
  colnames(P_star1) <- paste0(colnames(P1))

  # generate copula-correction terms
  P_star2 <- P_star1%*%solve(cor(apply(P_star1, 2, qnorm)))
  P_star2 <- P_star2[, c(1:length(independent_P_vars))]
  P_star2 <- as.matrix(P_star2)

  # rename columns of P_star matrix
  colnames(P_star2) <- paste0(independent_P_vars, "_cop")

  # merge columns
  est_matrix <- cbind(full_matrix, P_star2)

  # control function approach
  mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
             data = as.data.frame(est_matrix))

  Estimates <- mod1$coefficients
  return(Estimates)

}
boots_JAMS2 <- function(formula, data_cleaned, dependent_var, independent_P_vars,
                        independent_X_vars, has_intercept, X, cdf, factor_vars,
                        design_matrix1, full_formula, Ests) {
  
  repeat {
    
    data <- data_cleaned
    
    repeat {
      
      # Loop or logic starts here
      data_cleaned <- data %>%
        sample_n(size = nrow(data), replace = TRUE)
      design_matrix <- model.matrix(full_formula, data = data_cleaned)
      YP <- t(as.matrix(design_matrix))%*%as.matrix(design_matrix)
      rank_A <- Matrix::rankMatrix(YP)[1]
      
      # Check if for each factor variable, there are at least 3 observations for each level
      for (var in factor_vars) {
        levels_var <- levels(data_cleaned[[var]])
        for (lvl in levels_var) {
          num_observations <- sum(data_cleaned[[var]] == lvl)
          
          if (num_observations < 3) { break }
        }
        
      }
      
      # If the rank condition is met, or if the loop was broken due to insufficient observations, we proceed
      if (rank_A == ncol(design_matrix) &
          all(colnames(design_matrix1) %in% colnames(design_matrix))) { break }
      
    }
    
    # all regressor(s)
    data_cleaned_with_P_star2 <- data_cleaned
    full_matrix_with_P_Pstar2 <- data_cleaned_with_P_star2
    
    # Loop over each factor variable
    for (var in factor_vars) {
      
      levels_var <- levels(data_cleaned[[var]])  # Get factor levels
      
      for (lvl in levels_var) {
        
        subset_name <- paste(var, lvl, sep = "_")
        subdat1 <- subset(data_cleaned, data_cleaned[[var]] == lvl)
        
        if (nrow(subdat1) > 3 && 
            any(sapply(subdat1[independent_P_vars], function(col) length(unique(col)) > 1))) {
          
          subdat2 <- subdat1[, setdiff(c(independent_P_vars, independent_X_vars), factor_vars)]
          subdat2 <- as.matrix(subdat2)
          subdat2 <- subdat2[, apply(subdat2, 2, function(x) length(unique(x)) > 1), drop = FALSE]
          
          P_star1 <- matrix(NA, nrow = nrow(subdat2), ncol = ncol(subdat2))
          
          if (cdf == "kde") {
            
            for (i in 1:ncol(subdat2)) {
              Fhat <- ks::kcde(subdat2[, i])
              P_star1[, i] <- predict(Fhat, x = subdat2[, i])
            }
            
          } else if (cdf == "resc.ecdf") {
            
            P_star1 <- apply(subdat2, 2, copula::pobs)
            
          } else if (cdf == "adj.ecdf") {
            
            P_star1 <- apply(subdat2, 2, pobs1)
            
          } else if (cdf == "ecdf") {
            
            ecdf0 <- apply(subdat2, 2, ecdf)
            for (i in 1:ncol(subdat2)) {
              Fhat <- ecdf0[[i]]
              P_star1[, i] <- Fhat(subdat2[, i])
              P_star1[P_star1[, i] == min(P_star1[, i]), i] <- 10e-7
              P_star1[P_star1[, i] == max(P_star1[, i]), i] <- 1-10e-7
            }
            
          }
          
          # generate copula-correction terms
          is_issue <- tryCatch(
            { solve(cor(apply(P_star1, 2, qnorm))) ; FALSE },
            error = function(e) TRUE,
            warning = function(w) TRUE
          )
          
          if (is_issue) {
            
            P_star2 <- matrix(data = 0, nrow = nrow(P_star1), ncol = ncol(P_star1))
            P_star2 <- as.matrix((P_star2))
            
          } else {
            
            P_star2 <- P_star1%*%solve(cor(apply(P_star1, 2, qnorm)))
            
            if (length(independent_P_vars) > ncol(P_star2)) {
              P_star2 <- as.matrix((P_star2))
            } else {
              P_star2 <- P_star2[, c(1:length(independent_P_vars))]
              P_star2 <- as.matrix((P_star2))
            }
            
          }
          
          # rename columns of P_star matrix
          colnames(P_star2) <- paste0(independent_P_vars[c(1:ncol(P_star2))],
                                      "_", var, "_", lvl, "_cop")
          
          # merge data
          P_star2_full <- matrix(0, nrow = nrow(data_cleaned_with_P_star2), ncol = ncol(P_star2))
          P_star2_full[data_cleaned_with_P_star2[[var]] == lvl, ] <- P_star2
          P_star2_df <- as.data.frame(P_star2_full)
          colnames(P_star2_df) <- colnames(P_star2)
          
          full_matrix_with_P_Pstar2 <- cbind(full_matrix_with_P_Pstar2, P_star2_df)
          
        }
        
      }
      
    }
    
    # control function approach
    mod1 <- lm(as.formula(paste(dependent_var, "~ . -",
                                dependent_var)),
               data = full_matrix_with_P_Pstar2)
    
    Estimates <- mod1$coefficients
    Estimates <- Estimates[!is.na(Estimates)]
    
    # Check if the same variables are contained
    trapped_names <- names(Estimates)
    estimates_names <- names(Ests)
    
    # Remove "as.factor()" from the names if it appears
    trapped_names <- gsub("as.factor\\((.*?)\\)", "\\1", trapped_names)
    estimates_names <- gsub("as.factor\\((.*?)\\)", "\\1", estimates_names)
    estimates_names <- gsub("`", "", estimates_names)
    
    if (identical(trapped_names, estimates_names)) { break }
    
  }
  
  return(Estimates)

}

# Liengaard et al. (2024) estimator
CopRegJAMS <- function(formula, data, cdf, nboots = 199) {

  ################################################################################

  # check arguments
  if (rlang::is_formula(formula) == FALSE) { stop("Argument formula is not a formula object", call. = FALSE) }
  if (is.data.frame(data) == FALSE) { stop("Argument data is not a data.frame object", call. = FALSE) }
  if (inherits(data, "data.table")) { stop("Argument data should not be a data.table.", call = FALSE) }
  if (is.numeric(nboots) == FALSE) { stop("Argument nboots is not numeric", call. = FALSE) }
  if (!(cdf %in% c("kde", "resc.ecdf", "adj.ecdf", "ecdf"))) { stop("Invalid cdf Choose from: kde, ecdf, resc.ecdf, adj.ecdf", call. = FALSE) }

  # seperate endogenous and exogenous regressor(s)
  f1 <- nlme::splitFormula(formula, sep = "|")

  # check if intercept is removed
  has_intercept <- attr(terms(f1[[1]]), "intercept") == 1

  ##############################################################################

  # no exogenous regressors
  if (length(f1) == 1) {

    # no exogenous regressors
    f1P <- f1[[1]]

    # Check if all variables exist in the data
    variables <- all.vars(f1P)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {

      stop(paste("The following variables are missing in the data:",
                 paste(missing_vars, collapse=", ")))

    }

    # Check if all variables are numeric and non-constant
    numeric_vars <- sapply(data[variables], is.numeric)
    constant_vars <- sapply(data[variables], function(x) length(unique(x)) == 1)

    if (!all(numeric_vars)) {
      stop("The following are dummy variables and cannot be endogenous: ", paste(variables[!numeric_vars], collapse = ", "))
    }

    if (any(constant_vars)) {
      stop("The following variables are constant: ", paste(variables[constant_vars], collapse = ", "))
    }

    dependent_var <- all.vars(formula)[1]
    independent_vars <- all.vars(f1P)

    data_cleaned <- data %>%
      dplyr::select(all_of(c(dependent_var, independent_vars))) %>%
      na.omit() %>%
      as.data.frame()

    # Check if design matrix has full column rank
    if (!has_intercept) {

      # Use model.matrix to handle categorical variables
      full_matrix <- model.matrix(~ . - 1, data = data_cleaned)
      design_matrix <- full_matrix[, colnames(full_matrix) != dependent_var, drop = FALSE]

      # Check rank of the full and the design matrix
      rank_A <- rankMatrix(full_matrix)[1]
      rank_B <- rankMatrix(design_matrix)[1]

      if (rank_A != ncol(full_matrix)) {
        stop("Perfect fit")
      }
      if (rank_B != ncol(design_matrix)) {
        stop("Design matrix is rank deficient")
      }

    } else {

      # Add an intercept column to the design matrix
      full_matrix <- model.matrix(~ ., data = data_cleaned)
      design_matrix <- full_matrix[, colnames(full_matrix) != dependent_var, drop = FALSE]

      # Check rank of the full and the design matrix
      rank_A <- rankMatrix(full_matrix)[1]
      rank_B <- rankMatrix(design_matrix)[1]

      if (rank_A != ncol(full_matrix)) {
        stop("Perfect fit")
      }
      if (rank_B != ncol(design_matrix)) {
        stop("Design matrix is rank deficient")
      }

    }

    # endogenous regressor(s)
    if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
    P1 <- as.matrix(P1)
    P_star1 <- matrix(NA, nrow = nrow(P1), ncol = ncol(P1))

    if (cdf == "kde") {

      for (i in 1:ncol(P1)) {
        Fhat <- ks::kcde(P1[, i])
        P_star1[, i] <- predict(Fhat, x = P1[, i])
      }

    } else if (cdf == "resc.ecdf") {

      P_star1 <- apply(P1, 2, copula::pobs)

    } else if (cdf == "adj.ecdf") {

      P_star1 <- apply(P1, 2, pobs1)

    } else if (cdf == "ecdf") {

      ecdf0 <- apply(P1, 2, ecdf)
      for (i in 1:ncol(P1)) {
        Fhat <- ecdf0[[i]]
        P_star1[, i] <- Fhat(P1[, i])
        P_star1[P_star1[, i] == min(P_star1[, i]), i] <- 10e-7
        P_star1[P_star1[, i] == max(P_star1[, i]), i] <- 1-10e-7
      }

    }

    # rename columns of P_star matrix
    colnames(P_star1) <- paste0(colnames(P1), "_cop")

    # merge columns
    est_matrix <- cbind(full_matrix, P_star1)

    # control function approach
    mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
               data = as.data.frame(est_matrix))

    Estimates <- mod1$coefficients

    # obtain residuals
    beta <- Estimates[!grepl("_cop", names(Estimates))]
    fitted_values <- design_matrix %*% beta
    residuals_manual <- est_matrix[, dependent_var] - fitted_values

    # Bootstrapping
    print("Estimation done. Calculating bootstrap standard errors")
    trapped <- pbsapply(1:nboots,
                        function(i) boots1(formula = formula,
                                           data_cleaned = data_cleaned,
                                           dependent_var = dependent_var,
                                           independent_vars = independent_vars,
                                           has_intercept = has_intercept,
                                           cdf = cdf, X = i))
    ses <- apply(trapped, 1, sd)

    Estimates1 <- cbind(Estimates, ses)
    colnames(Estimates1) <- c("Estimate", "Std. Error")

  } else {

    ############################################################################

    f1P <- f1[[1]]
    f1X <- f1[[2]]

    variables <- all.vars(f1P)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {

      stop(paste("The following endogenous variables are missing in the data:",
                 paste(missing_vars, collapse=", ")))

    }

    variables <- all.vars(f1X)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {

      stop(paste("The following exogenous variables are missing in the data:",
                 paste(missing_vars, collapse=", ")))

    }

    dependent_var <- all.vars(formula)[1]
    independent_P_vars <- all.vars(f1P)
    independent_X_vars <- all.vars(f1X)

    # Check if all variables are numeric and non-constant
    numeric_vars <- sapply(data[independent_P_vars], is.numeric)
    constant_vars <- sapply(data[variables], function(x) length(unique(x)) == 1)

    if (!all(numeric_vars)) {
      stop("Only continuous variables can be endogenous. The following variables are not numeric: ", paste(variables[!numeric_vars], collapse = ", "))
    }

    if (any(constant_vars)) {
      stop("The following variables are constant: ",
           paste(variables[constant_vars], collapse = ", "))
    }

    # all variables
    data_cleaned <- data %>%
      dplyr::select(all_of(c(dependent_var, independent_P_vars,
                             independent_X_vars))) %>% na.omit() %>%
      as.data.frame()

    # Check if design matrix has full column rank
    full_formula <- as.formula(paste(all.vars(formula)[1], "~",
                                     gsub("\\|", "+", as.character(formula)[3])))
    full_matrix <- model.matrix(full_formula, data = data_cleaned)
    full_matrix <- cbind(DependentVar = data_cleaned[[dependent_var]], full_matrix)

    colnames(full_matrix)[1] <- dependent_var
    design_matrix <- full_matrix[, colnames(full_matrix) != dependent_var,
                                 drop = FALSE]

    rank_A <- Matrix::rankMatrix(full_matrix)[1]
    rank_B <- Matrix::rankMatrix(design_matrix)[1]

    if (rank_A != ncol(full_matrix)) {
      stop("Perfect fit")
    }
    if (rank_B != ncol(design_matrix)) {
      stop("Design matrix is rank deficient")
    }

    ############################################################################

    if (sum(sapply(data_cleaned, is.factor)) == 0 &
        !grepl("as.factor\\(", paste(deparse(full_formula), collapse = " "))) {

      # all regressor(s)
      if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
      P_star1 <- matrix(NA, nrow = nrow(P1), ncol = ncol(P1))

      if (cdf == "kde") {

        for (i in 1:ncol(P1)) {
          Fhat <- ks::kcde(P1[, i])
          P_star1[, i] <- predict(Fhat, x = P1[, i])
        }

      } else if (cdf == "resc.ecdf") {

        P_star1 <- apply(P1, 2, copula::pobs)

      } else if (cdf == "adj.ecdf") {

        P_star1 <- apply(P1, 2, pobs1)

      } else if (cdf == "ecdf") {

        ecdf0 <- apply(P1, 2, ecdf)
        for (i in 1:ncol(P1)) {
          Fhat <- ecdf0[[i]]
          P_star1[, i] <- Fhat(P1[, i])
          P_star1[P_star1[, i] == min(P_star1[, i]), i] <- 10e-7
          P_star1[P_star1[, i] == max(P_star1[, i]), i] <- 1-10e-7
        }

      }

      # rename columns of P_star matrix
      colnames(P_star1) <- paste0(colnames(P1))

      # generate copula-correction terms
      P_star2 <- P_star1%*%solve(cor(apply(P_star1, 2, qnorm)))
      P_star2 <- P_star2[, c(1:length(independent_P_vars))]
      P_star2 <- as.matrix(P_star2)

      # rename columns of P_star matrix
      colnames(P_star2) <- paste0(independent_P_vars, "_cop")

      # merge columns
      est_matrix <- cbind(full_matrix, P_star2)

      # control function approach
      mod1 <- lm(as.formula(paste(dependent_var, "~ . -",
                                  dependent_var, "- 1")),
                 data = as.data.frame(est_matrix))

      Estimates <- mod1$coefficients

      # obtain residuals
      beta <- Estimates[!grepl("_cop", names(Estimates))]
      fitted_values <- design_matrix %*% beta
      residuals_manual <- est_matrix[, dependent_var] - fitted_values

      # Bootstrapping
      print("Estimation done. Calculating bootstrap standard errors")
      trapped <- pbsapply(1:nboots,
                          function(i) boots_JAMS1(formula = formula,
                                                  data_cleaned = data_cleaned,
                                                  dependent_var = dependent_var,
                                                  independent_P_vars = independent_P_vars,
                                                  independent_X_vars = independent_X_vars,
                                                  has_intercept = has_intercept,
                                                  cdf = cdf, X = i,
                                                  full_formula = full_formula,
                                                  design_matrix1 = design_matrix))
      ses <- apply(trapped, 1, sd)

      Estimates1 <- cbind(Estimates, ses)
      colnames(Estimates1) <- c("Estimate", "Std. Error")

    } else {

      ##########################################################################

      # extract factor variables
      formula_str <- paste(deparse(full_formula), collapse = " ")
      factor_vars_formula <- regmatches(formula_str, gregexpr("as.factor\\((.*?)\\)", formula_str))[[1]]
      factor_vars_formula <- gsub("as.factor\\(|\\)", "", factor_vars_formula)

      factor_vars_data <- names(data_cleaned)[sapply(data_cleaned, is.factor)]
      factor_vars <- unique(c(factor_vars_formula, factor_vars_data))

      # all regressor(s)
      data_cleaned_with_P_star2 <- data_cleaned
      full_matrix_with_P_Pstar2 <- full_matrix
      data_cleaned[factor_vars] <- lapply(data_cleaned[factor_vars], as.factor)

      # Loop over each factor variable
      for (var in factor_vars) {

        levels_var <- levels(data_cleaned[[var]])  # Get factor levels

        for (lvl in levels_var) {

          subset_name <- paste(var, lvl, sep = "_")
          subdat1 <- subset(data_cleaned, data_cleaned[[var]] == lvl)

          if (nrow(subdat1) > 3 && 
              any(sapply(subdat1[independent_P_vars], function(col) length(unique(col)) > 1))) {
            
            subdat2 <- subdat1[, setdiff(c(independent_P_vars, independent_X_vars), factor_vars)]
            subdat2 <- as.matrix(subdat2)
            subdat2 <- subdat2[, apply(subdat2, 2, function(x) length(unique(x)) > 1), drop = FALSE]
            
            P_star1 <- matrix(NA, nrow = nrow(subdat2), ncol = ncol(subdat2))

            if (cdf == "kde") {

              for (i in 1:ncol(subdat2)) {
                Fhat <- ks::kcde(subdat2[, i])
                P_star1[, i] <- predict(Fhat, x = subdat2[, i])
              }

            } else if (cdf == "resc.ecdf") {

              P_star1 <- apply(subdat2, 2, copula::pobs)

            } else if (cdf == "adj.ecdf") {

              P_star1 <- apply(subdat2, 2, pobs1)

            } else if (cdf == "ecdf") {

              ecdf0 <- apply(subdat2, 2, ecdf)
              for (i in 1:ncol(subdat2)) {
                Fhat <- ecdf0[[i]]
                P_star1[, i] <- Fhat(subdat2[, i])
                P_star1[P_star1[, i] == min(P_star1[, i]), i] <- 10e-7
                P_star1[P_star1[, i] == max(P_star1[, i]), i] <- 1-10e-7
              }

            }

            # generate copula-correction terms
            is_issue <- tryCatch(
              { solve(cor(apply(P_star1, 2, qnorm))) ; FALSE },
              error = function(e) TRUE,
              warning = function(w) TRUE
            )

            if (is_issue) {

              P_star2 <- matrix(data = 0, nrow = nrow(P_star1), ncol = ncol(P_star1))
              P_star2 <- as.matrix((P_star2))

            } else {

              P_star2 <- P_star1%*%solve(cor(apply(P_star1, 2, qnorm)))

              if (length(independent_P_vars) > ncol(P_star2)) {
                P_star2 <- as.matrix((P_star2))
              } else {
                P_star2 <- P_star2[, c(1:length(independent_P_vars))]
                P_star2 <- as.matrix((P_star2))
              }

            }

            # rename columns of P_star matrix
            colnames(P_star2) <- paste0(independent_P_vars[c(1:ncol(P_star2))],
                                        "_", var, "_", lvl, "_cop")

            # merge data
            P_star2_full <- matrix(0, nrow = nrow(data_cleaned_with_P_star2), ncol = ncol(P_star2))
            P_star2_full[data_cleaned_with_P_star2[[var]] == lvl, ] <- P_star2
            P_star2_df <- as.data.frame(P_star2_full)
            colnames(P_star2_df) <- colnames(P_star2)

            full_matrix_with_P_Pstar2 <- cbind(full_matrix_with_P_Pstar2, P_star2_df)

          }

        }

      }

      # control function approach
      mod1 <- lm(as.formula(paste(dependent_var, "~ . -",
                                  dependent_var)),
                 data = full_matrix_with_P_Pstar2)

      Estimates <- mod1$coefficients
      Estimates <- Estimates[!is.na(Estimates)]

      # obtain residuals
      beta <- Estimates[!grepl("_cop", names(Estimates))]
      fitted_values <- design_matrix %*% beta
      residuals_manual <- data_cleaned_with_P_star2[, dependent_var] - fitted_values

      # Bootstrapping
      print("Estimation done. Calculating bootstrap standard errors")
      trapped <- pbsapply(1:nboots,
                          function(i) boots_JAMS2(formula = formula,
                                                  data_cleaned = data_cleaned,
                                                  dependent_var = dependent_var,
                                                  independent_P_vars = independent_P_vars,
                                                  independent_X_vars = independent_X_vars,
                                                  has_intercept = has_intercept,
                                                  cdf = cdf, X = i,
                                                  factor_vars = factor_vars,
                                                  full_formula = full_formula,
                                                  design_matrix1 = design_matrix,
                                                  Ests = Estimates))
      
      if (is.list(trapped)) {
        common_names <- Reduce(intersect, lapply(trapped, names))
        trapped_df <- do.call(rbind, lapply(trapped, function(x) x[common_names]))
        ses <- apply(trapped_df, 2, sd)
      } else { ses <- apply(trapped, 1, sd) }

      Estimates1 <- cbind(Estimates, ses)
      colnames(Estimates1) <- c("Estimate", "Std. Error")

    }

  }

  return(list(Estimates1, residuals_manual))

}


# This function implements the copula-based endogeneity correction by
# Liengaard et al. (2024) using the least-squares-based correction function
# approach.
#
# formula = depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + as.factor(exog_var2) + ...
#
# data = as.data.frame(datset)
#
# cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf")
# kde is the integral of a density estimator used in Park & Gupta (2012) and Haschka (2022)
# ecdf is the empirical cumulative distribution function (ecdf) with replaced boundary proposed by Becker et al. (2022)
# resc.ecdf is a rescaled ecdf proposed by Qian et al. (2024)
# adj.ecdf is an adjusted ecdf proposed by Liengaard (2024)
#
# CopRegJAMS returns a list of legth 2. First entry are estimates with standard
# errors. Second entry are residuals.


### REFERENCES

# Becker, J.-M., D. Proksch, and C. M. Ringle (2021). Revisiting Gaussian
# copulas to handle endogenous regressors. Journal of the Academy of Marketing
# Science 50, 46–66.
#
# Haschka, R. E (2022). Handling endogenous regressors using copulas: A
# generalisation to linear panel models with fixed effects and correlated
# regressors. Journal of Marketing Research 59(4), 860–881.
#
# Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and
# C. M. Ringle (2024). Dealing with regression models’ endogeneity by means of
# an adjusted estimator for the Gaussian copula approach. Journal of the
# Academy of Marketing Science, 1–21.
#
# Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
# estimation using copulas. Marketing Science 31 (4), 567–586.
#
# Qian, Y., A. Koschmann, and H. Xie (2024). A practical guide to endogeneity
# correction using copulas. NBER Working Paper.
# https://www.nber.org/system/files/workingpapers/w32231/w32231.pdf.





# 1. Example

library(bayesm)

data("orangeJuice")
dat1 <- orangeJuice[[1]]

dat1_Tropicana <- dat1 %>% filter(brand == 4)
dat1_Tropicana <- dat1_Tropicana %>%
  mutate(across(starts_with("price"), log))


CopRegJAMS1 <- CopRegJAMS(formula = logmove ~ price4 | price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10,
                           data = dat1_Tropicana, cdf = "adj.ecdf")
CopRegJAMS1[[1]]
hist(CopRegJAMS1[[2]])


CopRegJAMS2 <- CopRegJAMS(formula = logmove ~ price4 + price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10,
                           data = dat1_Tropicana, cdf = "adj.ecdf")
CopRegJAMS2[[1]]
hist(CopRegJAMS2[[2]])


CopRegJAMS3 <- CopRegJAMS(formula = logmove ~ price4 + price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10 | feat + as.factor(deal),
                           data = dat1_Tropicana, cdf = "adj.ecdf")
CopRegJAMS3[[1]]
hist(CopRegJAMS3[[2]])


CopRegJAMS4 <- CopRegJAMS(formula = logmove ~ price4 + price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10 | feat + as.factor(deal) + as.factor(store),
                           data = dat1_Tropicana, cdf = "adj.ecdf")
CopRegJAMS4[[1]]
hist(CopRegJAMS4[[2]])



# 2. Example

library(AER)

data1 <- data("CPS1988")
data1 <- CPS1988
data1$lwage <- log(data1$wage)
data1$experience_sq <- data1$experience^2

modJAMS1 <- CopRegJAMS(formula = lwage ~ education + experience | as.factor(experience_sq) + as.factor(parttime) + smsa + ethnicity,
                           data = data1, cdf = "ecdf")
modJAMS1[[1]]


modJAMS2 <- CopRegJAMS(formula = lwage ~ education + experience | experience_sq + parttime + smsa + ethnicity,
                           data = data1, cdf = "ecdf")
modJAMS2[[1]]



# 3. Example

library(ISLR)

data("Carseats")
dat1 <- Carseats

dat1$lsales <- log(dat1$Sales)
dat1$lprice <- log(dat1$Price)
dat1$lcompprice <- log(dat1$CompPrice)
dat1 <- subset(dat1, is.finite(log(Sales)))

modJAMS1 <- CopRegJAMS(formula = lsales ~ lprice + lcompprice | ShelveLoc + Income + Advertising + Population + Age + Education + Urban + US,
                           data = dat1, cdf = "ecdf")
modJAMS1[[1]]


dat1$ShelveLoc_num <- as.numeric(dat1$ShelveLoc)

modJAMS2 <- CopRegJAMS(formula = lsales ~ lprice + lcompprice | Income + Advertising + Population + Age + Education + as.factor(ShelveLoc_num) + as.factor(Urban) + US,
                           data = dat1, cdf = "ecdf")
modJAMS2[[1]]



