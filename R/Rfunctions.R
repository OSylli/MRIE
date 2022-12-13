#' Mendelian Randomization (MR) Analysis with Incomplete Measurements on the Exposure
#' 
#' This function is used for MR analysis with a continuous exposure and a continuous 
#' outcome, where the exposure variable is potentially unmeasured or non-detectable. 
#' The estimations are carried out based on the maximum likelihood estimation method, 
#' and the EM algorithm is applied to cope with the incomplete measurements on the exposure.
#' 
#' @param IV.dat A data frame containing data on the instrumental variables (IVs) for all 
#' the subjects, with a header line indicating the names of the IVs. Each row corresponds 
#' to one subject, and each column corresponds to one IV. Values can be SNP dosages, 
#' polygenic risk scores, etc. Missing value or \code{NA} is not allowed.
#' @param pheno.dat A data frame containing data on the outcome and the exposure, with a 
#' header line indicating their names. Missing value or \code{NA} is not allowed. Each row 
#' corresponds to one subject, and the subjects should be listed in the same order as 
#' in \code{IV.dat}. The dataset\code{pheno.dat} should have 3 columns. The 1st column 
#' contains the floating-point values for the continuous outcome. The 2nd column (integer 
#' valued) and the 3rd column (floating-point valued) contain information on the exposure 
#' variable. See the "Details" section for instructions on how the 2nd and the 3rd column 
#' should be specified.
#' @param covar.dat An optional data frame containing data on the measured covariates 
#' (e.g., age, gender, and race), with a header line indicating the names of these variables. 
#' Each row corresponds to one subject, and the subjects should be listed in the same order 
#' as in \code{IV.dat}. Each column corresponds to one covariate. Missing value or \code{NA} 
#' is not allowed.
#' @param epsilon The convergence threshold of the EM algorithm. The iterations will be 
#' terminated if the Euclidean distance between the parameter values at two successive 
#' iterations is less than this value. Default as \eqn{10^{-5}}.
#' 
#' @details 
#' Let \eqn{Y} be a continuous outcome, \eqn{S} be a continuous exposure which is 
#' potentially unmeasured or subject to detection limits, \eqn{G} be a 
#' vector of instrumental variables for \eqn{S}, and \eqn{Z} be a vector of measured covariates. 
#' Assume that the first component of \eqn{Z} is 1, and let \eqn{X = (G^T, Z^T)^T}. The 
#' following model is considered:
#' \deqn{S = \alpha'X + \epsilonS,}
#' \deqn{Y = \gammaS + \beta'Z + \epsilonY,}
#' where \eqn{\alpha} and \eqn{\beta} are vectors of regression parameters, \eqn{\gamma} 
#' represents the causal effect of the exposure on the outcome, and 
#' \eqn{(\epsilonS, \epsilonY)^T} is a bivariate normal random vector with mean zero and 
#' an unstructured covariance matrix. We allow a nonzero correlation between \eqn{\epsilonS} 
#' and \eqn{\epsilonY} to account for the potential unmeasured confounders of the 
#' exposure-outcome relationship.
#' 
#' The exposure variable \eqn{S} may be unmeasured or subject to detection limits. This 
#' function employs the EM algorithm to cope with these incomplete measurements and derive 
#' the parameter estimates, and the covariance matrix of the estimated parameters is derived 
#' from the Louis formula (Little and Rubin, 2019).
#' 
#' For \code{pheno.dat}, the i-th element of the 2nd column indicates whether the 
#' exposure of i-th subject is observed, beyond detection limits, or unmeasured, and the 
#' corresponding value in the 3rd column should be the exact measurement on the exposure, the 
#' detection limit, and a dummy value of -999, respectively. Below is a table showing how 
#' these two columns should be specified when the exposure variable is observed, beyond 
#' detection limits, or unmeasured.
#' 
#' |						                                   | The 2nd column | The 3rd column		        |
#' |:--------------------------------------------- |:--------------:|:-------------------------:|
#' | Measured and Detectable			                 |       0        | the exact measurement	    |
#' | Measured but below the lower detection limit	 |       1        | the lower detection limit	|
#' | Measured but above the upper detection limit	 |       2        | the upper detection limit	|
#' | Unmeasured					                           |       3        | -999				              |
#' 
#' @return 
#' A list will be returned with the following elements:
#' \itemize{
#' \item{\code{results_reg} contains the parameter estimates, the standard error estimates, and the p-values for \eqn{\gamma} and the components in \eqn{\alpha} and \eqn{\beta}.}
#' \item{\code{results_var} contains the estimates for the variance components (i.e., the variances of \eqn{\epsilonS} and \eqn{\epsilonY}, and the correlation between \eqn{\epsilonS} and \eqn{\epsilonY}).}
#' \item{\code{full_cov_mat} is the estimated covariance matrix derived from the Louis formula.}
#' }
#' 
#' @importFrom stats pchisq
#' 
#' @references 
#' Little, R. J., & Rubin, D. B. (2019). Statistical analysis with missing data (3rd Edition). John Wiley & Sons.
#' 
#' @useDynLib MRIE
#' @export
MRIE_EM <- function(IV.dat, pheno.dat, covar.dat = NULL, epsilon = 1e-5){
  ## input column names
  outcome.name <- colnames(pheno.dat)[1]
  exposure.name <- colnames(pheno.dat)[3]
  
  ## checking the input datasets - # observations
  n <- nrow(IV.dat)
  if(nrow(pheno.dat) != n){
    stop("The instrumental variable dataset and the phenotype dataset have different numbers of observations.")
  }
  
  ## checking the input datasets - NA's
  if(sum(is.na(IV.dat)) > 0){
    stop("NA is not allowed in the instrumental variable dataset.")
  }
  if(sum(is.na(pheno.dat)) > 0){
    stop("NA is not allowed in the phenotype dataset.")
  }
  if(!is.null(covar.dat)){
    if(sum(is.na(covar.dat)) > 0){
      stop("NA is not allowed in the instrumental variable dataset.")
    }
  }
  
  ## Preparations for the arguments of the CPP function
  #### Outcome (Y)
  Y <- as.vector(pheno.dat[, 1])
  
  #### Exposure (S)
  S <- as.vector(pheno.dat[, 3])
  R <- as.vector(pheno.dat[, 2])
  Check_R <- (R != 0) & (R != 1) & (R != 2) & (R != 3)
  if(sum(Check_R) > 0){
    stop("Values in the 2nd column of the phenotype dataset should only be 0, 1, 2, or 3.")
  }
  
  #### intercept & the independent variables for the model equations
  Z <- data.frame(intercept = rep(1, n))
  if(!is.null(covar.dat)){
    if(nrow(covar.dat) != n){
      stop("The instrumental variable dataset and the covariate dataset have different numbers of observations.")
    }
    Z <- cbind(Z, covar.dat)
  }
  X <- cbind(IV.dat, Z)
  
  colnames_Z <- colnames(Z)
  colnames_X <- colnames(X)
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  
  ## Preparation for Imp-Mid
  R_sub <- R[R != 3]
  Y_sub <- Y[R != 3]
  S_sub <- S[R != 3]
  S_sub[which(R_sub == 1)] <- S_sub[which(R_sub == 1)] - log(2)
  Z_sub <- Z[R != 3, ]
  X_sub <- X[R != 3, ]
  
  Z_sub <- as.matrix(Z_sub)
  X_sub <- as.matrix(X_sub)
  
  ## Analysis - Imputation at mid-point
  gamma_index <- ncol(X_sub) + ncol(Z_sub) + 1
  init_ImpMid <- make_initial(Y_sub, S_sub, Z_sub, X_sub)
  fit <- MR_estimation_NoMissing_CORR(X_sub, Y_sub, Z_sub, S_sub,
                                      init_ImpMid[[1]], init_ImpMid[[2]], init_ImpMid[[3]],
                                      init_ImpMid[[4]], init_ImpMid[[5]], init_ImpMid[[6]], epsilon)
  
  ## Analysis - Proposed method using the EM algorithm
  EM_fit <- MR_estimation_rho(Y, S, R, Z, X, fit[[1]], fit[[2]], fit[[3]], fit[[4]], fit[[5]], fit[[6]], epsilon)
  
  ## output row & column names
  output.colnames <- c("Estimate", "Std. Error", "P-value")
  rowname1 <- paste(colnames_X, " on ", exposure.name, sep = "", collapse = NULL)
  rowname2 <- paste(colnames_Z, " on ", outcome.name, sep = "", collapse = NULL)
  rowname3 <- paste(exposure.name, " on ", outcome.name, sep = "", collapse = NULL)
  output.rownames <- c(rowname1, rowname2, rowname3)
  
  ## result data.frame & List
  results_reg <- matrix(NA, nrow = length(output.rownames), ncol = length(output.colnames), dimnames = list(output.rownames, output.colnames))
  results_var <- matrix(NA, nrow = 3, ncol = 1, dimnames = list(c("variance-S", "variance-Y", "correlation"), c("Estimate")))
  
  #### estimates
  results_reg[, 1] <- c(EM_fit[[1]], EM_fit[[2]], EM_fit[[3]])
  
  #### standard error estimates
  results_reg[, 2] <- sqrt(diag(EM_fit[[7]])[1 : length(output.rownames)])
  
  #### two-sided p-value
  z_score <- results_reg[, 1] / results_reg[, 2]
  results_reg[ ,3] <- 1 - pchisq(z_score^2, 1)
  
  #### variance components
  results_var[1, 1] <- EM_fit[[4]]
  results_var[2, 1] <- EM_fit[[5]]
  results_var[3, 1] <- EM_fit[[6]]
  
  #### full covariance matrix
  full_para_name <- c(output.rownames, "variance-S", "variance-Y", "correlation")
  full_cov_mat <- EM_fit[[7]]
  dimnames(full_cov_mat) <- list(full_para_name, full_para_name)
  
  results_list <- list(results_reg, results_var, full_cov_mat)
  return(results_list)
}
