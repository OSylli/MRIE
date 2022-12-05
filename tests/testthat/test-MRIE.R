test_that("NA values in the input dataset", {
  # model construction
  n <- 5000
  prob <- 0.3

  # independent variables and error terms
  epsilon_S <- rnorm(n)
  epsilon_Y <- rnorm(n)
  G <- rbinom(n, 2, prob)

  # parameters
  alpha <- 0.25
  gamma <- 0.20
  intercept <- 0.5

  # dependent variables
  S <- alpha * G + intercept + epsilon_S
  R <- rep(0, n)
  Y <- gamma * S + intercept + epsilon_Y
  Y[1 : (0.01*n)] <- NA

  IV.dat <- data.frame(IV1 = G)
  pheno.dat <- data.frame(outcome = Y, type = R, exposure = S)
  expect_error(MRIE_EM(IV.dat, pheno.dat), "NA is not allowed in the phenotype dataset.")
})


test_that("different number of observations", {
  # model construction
  n <- 5000
  prob <- 0.3

  # independent variables and error terms
  epsilon_S <- rnorm(n)
  epsilon_Y <- rnorm(n)
  G <- rbinom(n, 2, prob)

  # parameters
  alpha <- 0.25
  gamma <- 0.20
  intercept <- 0.5

  # dependent variables
  S <- alpha * G + intercept + epsilon_S
  R <- rep(0, n)
  Y <- gamma * S + intercept + epsilon_Y

  IV.dat <- data.frame(IV1 = G[1 : (0.5 * n)])
  pheno.dat <- data.frame(outcome = Y, type = R, exposure = S)
  expect_error(MRIE_EM(IV.dat, pheno.dat), "The instrumental variable dataset and the phenotype dataset have different numbers of observations.")
})


test_that("Wrong specifications of the measurement type", {
  # model construction
  n <- 5000
  prob <- 0.3

  # independent variables and error terms
  epsilon_S <- rnorm(n)
  epsilon_Y <- rnorm(n)
  G <- rbinom(n, 2, prob)

  # parameters
  alpha <- 0.25
  gamma <- 0.20
  intercept <- 0.5

  # dependent variables
  S <- alpha * G + intercept + epsilon_S
  R <- rep(0, n)
  Y <- gamma * S + intercept + epsilon_Y

  S[1 : (0.5 * n)] <- -999
  R[1 : (0.5 * n)] <- 4

  detection_limit <- 0.2
  for(i in ((0.5 * n + 1):n)){
    if(S[i] < detection_limit){
      S[i] <- detection_limit
      R[i] <- 1
    }
  }

  IV.dat <- data.frame(IV1 = G)
  pheno.dat <- data.frame(outcome = Y, type = R, exposure = S)
  expect_error(MRIE_EM(IV.dat, pheno.dat), "Values in the 2nd column of the phenotype dataset should only be 0, 1, 2, or 3.")
})


