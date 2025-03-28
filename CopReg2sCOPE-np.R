# Load required packages
pacman::p_load(
  nlme,
  dplyr,
  Matrix,
  ks,
  pbapply,
  np
)


# functions for bootstrapping
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
boots_np <- function(formula, data_cleaned, dependent_var, independent_P_vars,
                     independent_X_vars, X, cdf, design_matrix1, regressors, 
                     full_formula, has_intercept) {
  
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
  
  # conditional cdf estimation
  P_star <- matrix(NA, nrow = nrow(data_cleaned), 
                   ncol = length(independent_P_vars))
  colnames(P_star) <- independent_P_vars
  
  if (has_intercept) { design_matrix <- design_matrix[, -1] } else { design_matrix <- design_matrix }
  regressors <- setdiff(colnames(design_matrix), independent_P_vars)
  regressors <- sapply(regressors, function(x) {
    if (grepl("\\)", x)) paste0("`", x, "`") else x
  })
  
  for (i in seq_along(independent_P_vars)) {
    p_var <- independent_P_vars[i]
    bw <- npcdistbw(ydat = as.data.frame(full_matrix)[p_var], 
                    xdat = as.data.frame(full_matrix)[names(regressors)])
    cdf1 <- npcdist(bws = bw)
    P_star[, i] <- cdf1$condist
  }
  
  # rename columns of P_star matrix
  colnames(P_star) <- paste0(independent_P_vars, "_cop")
  P_star <- apply(P_star, 2, qnorm)
  
  # merge columns
  est_matrix <- cbind(full_matrix, P_star)
  
  # control function approach
  mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
             data = as.data.frame(est_matrix))
  
  Estimates <- mod1$coefficients
  return(Estimates)
  
}

# Hu et. al. (2025) estimator
CopReg2sCOPEnp <- function(formula, data, nboots = 3) {
  
  ################################################################################
  
  # check arguments
  if (rlang::is_formula(formula) == FALSE) { stop("Argument formula is not a formula object", call. = FALSE) }
  if (is.data.frame(data) == FALSE) { stop("Argument data is not a data.frame object", call. = FALSE) }
  if (inherits(data, "data.table")) { stop("Argument data should not be a data.table.", call = FALSE) }
  if (is.numeric(nboots) == FALSE) { stop("Argument nboots is not numeric", call. = FALSE) }

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
      rank_A <- Matrix::rankMatrix(full_matrix)[1]
      rank_B <- Matrix::rankMatrix(design_matrix)[1]
      
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
      rank_A <- Matrix::rankMatrix(full_matrix)[1]
      rank_B <- Matrix::rankMatrix(design_matrix)[1]
      
      if (rank_A != ncol(full_matrix)) {
        stop("Perfect fit")
      }
      if (rank_B != ncol(design_matrix)) {
        stop("Design matrix is rank deficient")
      }
      
    }
    
    # endogenous regressor(s)
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
      stop("Only continuous variables can be endogenous. The following variables are not numeric: ", 
           paste(variables[!numeric_vars], collapse = ", "))
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
    
    P_star <- matrix(NA, nrow = nrow(data_cleaned), 
                     ncol = length(independent_P_vars))
    colnames(P_star) <- independent_P_vars
    
    if (has_intercept) { design_matrix <- design_matrix[, -1] } else { design_matrix <- design_matrix }
    regressors <- setdiff(colnames(design_matrix), independent_P_vars)
    regressors <- sapply(regressors, function(x) {
      if (grepl("\\)", x)) paste0("`", x, "`") else x
    })
    
    for (i in seq_along(independent_P_vars)) {
      p_var <- independent_P_vars[i]
      
      print(paste("Computing cdf of endogenous regressor", paste(p_var), 
                  "conditional on exogenous regressors", paste(names(regressors), collapse=", ")))
      
      bw <- npcdistbw(ydat = as.data.frame(full_matrix)[p_var], 
                      xdat = as.data.frame(full_matrix)[names(regressors)])
      cdf1 <- npcdist(bws = bw)
      P_star[, i] <- cdf1$condist
    }
    
    # rename columns of P_star matrix
    colnames(P_star) <- paste0(independent_P_vars, "_cop")
    P_star <- apply(P_star, 2, qnorm)
    
    # merge columns
    est_matrix <- cbind(full_matrix, P_star)
    
    # control function approach
    mod1 <- lm(as.formula(paste(dependent_var, "~ . -", 
                                dependent_var, "- 1")),
               data = as.data.frame(est_matrix))
    
    Estimates <- mod1$coefficients
    
    # obtain residuals
    beta <- Estimates[!grepl("_cop", names(Estimates))]
    if (has_intercept) { fitted_values <- cbind(rep(1, nrow(design_matrix)), 
                                                design_matrix) %*% beta } else {
                                                  fitted_values <- design_matrix %*% beta }
    residuals_manual <- est_matrix[, dependent_var] - fitted_values
    
    # Bootstrapping
    print("Estimation done. Calculating bootstrap standard errors")
    trapped <- pbsapply(1:nboots, 
                        function(i) boots_np(formula = formula,
                                             data_cleaned = data_cleaned,
                                             dependent_var = dependent_var,
                                             independent_P_vars = independent_P_vars,
                                             independent_X_vars = independent_X_vars,
                                             cdf = cdf, design_matrix1 = design_matrix,
                                             regressors = regressors, X = i,
                                             full_formula = full_formula,
                                             has_intercept))
    ses <- apply(trapped, 1, sd)
    
    Estimates1 <- cbind(Estimates, ses)
    colnames(Estimates1) <- c("Estimate", "Std. Error")
    
  }
  
  return(list(Estimates1, residuals_manual))
  
}


# This function implements the copula-based endogeneity correction by 
# Hu et. al. (2025) using the least-squares-based correction function 
# approach.
#
# formula = depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + exog_var2 + ...
#
# data = as.data.frame(datset)
#
# CopReg2sCOPEnp returns a list of legth 2. First entry are estimates with 
# standard errors. Second entry are residuals.


### REFERENCES

# Hu, X., Qian, Y., A., and H. Xie (2025). Correcting endogeneity via 
# instrument-free two-stage nonparametric copula control functions
# http://www.nber.org/papers/w33607





# 1. Example

library(bayesm)

data("orangeJuice")
dat1 <- orangeJuice[[1]]

dat1_Tropicana <- dat1 %>% filter(brand == 4)
dat1_Tropicana <- dat1_Tropicana %>%
  mutate(across(starts_with("price"), log))

mod2sCOPEnp1 <- CopReg2sCOPEnp(formula = logmove ~ price4 | price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10,
                           data = dat1_Tropicana)
mod2sCOPEnp1[[1]]
hist(mod2sCOPEnp1[[2]])


mod2sCOPEnp2 <- CopReg2sCOPEnp(formula = logmove ~ price4 + price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10,
                               data = dat1_Tropicana)
mod2sCOPEnp2[[1]]
hist(mod2sCOPEnp2[[2]])


mod2sCOPEnp3 <- CopReg2sCOPEnp(formula = logmove ~ price4 + price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10 | feat + as.factor(deal),
                               data = dat1_Tropicana)
mod2sCOPEnp3[[1]]
hist(mod2sCOPEnp3[[2]])


mod2sCOPEnp4 <- CopReg2sCOPEnp(formula = logmove ~ price4 + price1 + price2 + price3 + price5 + price6 + price7 + price8 + price9 + price10 | feat + as.factor(deal) + as.factor(store),
                               data = dat1_Tropicana)
mod2sCOPEnp4[[1]]
hist(mod2sCOPEnp4[[2]])



# 2. Example

library(AER)

data1 <- data("CPS1988")
data1 <- CPS1988
data1$lwage <- log(data1$wage)
data1$experience_sq <- data1$experience^2

mod2sCOPEnp1 <- CopReg2sCOPEnp(formula = lwage ~ education + experience | as.factor(experience_sq) + as.factor(parttime) + smsa + ethnicity,
                     data = data1)
mod2sCOPEnp1[[1]]

mod2sCOPEnp2 <- CopReg2sCOPEnp(formula = lwage ~ education + experience | experience_sq + as.factor(parttime) + smsa + ethnicity,
                               data = data1)
mod2sCOPEnp2[[1]]



# 3. Example

library(ISLR)

data("Carseats")
dat1 <- Carseats

dat1$lsales <- log(dat1$Sales)
dat1$lprice <- log(dat1$Price)
dat1$lcompprice <- log(dat1$CompPrice)
dat1 <- subset(dat1, is.finite(log(Sales)))

mod2sCOPEnp1 <- CopReg2sCOPEnp(formula = lsales ~ lprice + lcompprice | ShelveLoc + Income + Advertising + Population + Age + Education + Urban + US,
                     data = dat1)
mod2sCOPEnp1[[1]]


dat1$ShelveLoc_num <- as.numeric(dat1$ShelveLoc)

mod2sCOPEnp1 <- CopReg2sCOPEnp(formula = lsales ~ lprice + lcompprice | Income + Advertising + Population + Age + Education + as.factor(ShelveLoc_num) + as.factor(Urban) + US,
                     data = dat1, cdf = "ecdf")
mod2sCOPE2[[1]]



