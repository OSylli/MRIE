# MRIE (Mendelian Randomization with Incomplete Measurements on the Exposure)

## Description

This package is used for Mendelian Randomization (MR) analysis with a continuous outcome $Y$ and a continuous exposure $S$, which is potentially unmeasured or subject to detection limits. The estimations are carried out based on the maximum likelihood estimation (MLE) method, and the expectation-maximization (EM) algorithm is used for the computation. The covariance matrix of the estimated parameters is derived from the Louis formula described in Section 8.4 in [Little and Rubin (2019)](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119482260).

To be more specific, let $G$ be a vector of instrumental variables (IVs) for $S$ and ${Z}$ be a vector of measured covariates. Assume that the first component of ${Z}$ is the constant 1, and let the vector ${X} = ({G}^T, {Z}^T)^T$ . We consider the following models:
$$S = {\alpha}^T{X} + \epsilon_S,$$
$$Y = \gamma S + {\beta}^T{Z} + \epsilon_Y,$$
where ${\alpha}$ and ${\beta}$ are vectors of regression parameters, $\gamma$ represents the causal effect of $S$ on $Y$, and $(\epsilon_S, \epsilon_Y)^T$ is a bivariate normal random vector with mean zero and an unstructured covariance matrix. To accommodate the unmeasured confounders that affect both $S$ and $Y$, we allow the error terms $\epsilon_S$ and $\epsilon_Y$ to be correlated. 



## Installation

The package `MRIE` requires the following R packages: `Rcpp` (>= 0.11.0), `RcppArmadillo`, `data.table`, and `stats`. You can install them by
```{r}
install.packages(c("Rcpp", "RcppArmadillo", "data.table", "stats"), dependencies=TRUE)
```

To install `MRIE`, you can first install and load the package `devtools` and then run the following code 
```{r}
install_github("OSylli/MRIE")
```
for the latest version. You can also clone the latest repository here and install it by yourself using `devtools`.



## Tutorials

### Function arguments

The primary function of the package is `MRIE_EM`, which has four arguments:

* The first argument `IV.dat` is a required data frame containing data on the IVs for all the subjects, with a header line indicating the names of the IVs. Each row corresponds to one subject, and each column corresponds to one IV, such as SNP dosages and polygenic risk scores. Missing value or `NA` is not allowed.
* The second argument `pheno.dat` is a required data frame containing data on the outcome and the exposure, with a header line indicating their names. Missing value or `NA` is not allowed. Each row corresponds to one subject, and the subjects should be listed in the same order as in `IV.dat`. The dataset `pheno.dat` should have three columns. 
  + The first column contains the floating-point values for the continuous outcome. 
  + The second column (integer-valued) and the third column (floating-point valued) contain information on the exposure variable. Specifically, the $i$th element of the second column indicates whether the exposure of the $i$th subject is observed, beyond detection limits, or unmeasured, and the corresponding value in the third column should be the exact measurement on the exposure, the detection limit, and a dummy value of -999, respectively. Below is a table showing how the second and the third column should be specified when the exposure variable is observed, beyond detection limits, or unmeasured:

Measurement type | The 2nd column | The 3rd column 
:-----------------:|:----------------:|:---------------:|
Measured and detectable | 0 | the exact measurement  
Measured but below the lower detection limit | 1 | the lower detection limit  
Measured but above the upper detection limit | 2 | the upper detection limit  
Unmeasured | 3 | -999  

* The third argument `covar.dat` is an optional data frame containing data on the measured covariates (e.g., age, gender, and race), with a header line indicating the names of these variables. Each row corresponds to one subject, and the subjects should be listed in the same order as in `IV.dat`. Each column corresponds to one covariate. Missing value or `NA` is not allowed.
* The last argument `epsilon` is the convergence threshold of the EM algorithm. The iterations will be terminated if the Euclidean distance between the parameter values at two successive iterations is less than this value. Default as $10^{-5}$.

### Output values

When the function `MRIE_EM` is run correctly, a list will be returned with the following elements:

* `results_reg` contains the parameter estimates for $\gamma$ and the components in $\alpha$ and $\beta$. The standard error estimates of the estimated parameters and the corresponding $p$-values are also provided. The inference on the causal effect of interest is presented in the last row.
* `results_var` contains the estimates for the variance components (i.e., the variances of $\epsilon_S$ and $\epsilon_Y$ and the correlation between $\epsilon_S$ and $\epsilon_Y$).
* `full_cov_mat` is the estimated covariance matrix of the parameter estimates derived from the Louis formula.

### Example

We show the implementation of the proposed function with the following example. Firstly, we load the example dataset.
```{r}
# set the working directory as the path where the example data were stored
library(MRIE)
setwd(system.file("extdata", package="MRIE"))

data <- read.table("example_data.txt", header = TRUE, sep = "\t")
dim(data)
head(data)
tail(data)
unique(data$Type)
```

The example dataset has seven columns:

* The first column contains the measurements on the continuous outcome of interest.
* The second column shows the exact measurement values of the exposure if the exposure is measured and detectable, the lower detection limit if the exposure is measured but below the lower detection limit, and a dummy value of `-999` if it is unmeasured.
* The third column indicates the measurement type of the exposure, where the value `0` suggests that the exposure is measured and detectable, the value `1` corresponds to the exposure that is measured but below the lower detection limit, and the value `3` means that the exposure is unmeasured, as shown in the "Function argument" section above.
* The fourth column contains the data on the IV.
* The last three columns contain measurements on the measured covariates of age, gender, and the first principal component for ancestry.

Then, we can prepare for the function arguments:
```{r}
# The argument `IV.dat` is a data frame containing measurements on the IVs
IV.dat <- data.frame(IV = data$IV)

# The argument `pheno.dat` is a data frame containing information on the outcome and the exposure
pheno.dat <- data.frame(Outcome = data$Y, Type = data$Type, Exposure = data$S)

# The argument `covar.dat` is a data frame containing measurements on the measured covariates
covar.dat <- data.frame(AGE = data$AGE, GENDER = data$GENDER, PC = data$PC)

# Analyze the example data with the function `MRIE_EM` in the package
fit <- MRIE_EM(IV.dat, pheno.dat, covar.dat)
```

We can then check the results as following
```{r}
# The inference on the regression parameters
fit$results_reg

# The estimates of the variance components
fit$results_var

# The full estimated covariance matrix derived from the Louis formula
fit$full_cov_mat
```

Below is the output (of `results_reg` and `results_var`) for the analysis using the example data. The inference on the causal effect of interest is presented in the last row of `results_reg`.

```
> fit$results_reg
                       Estimate Std. Error      p-value
IV on Exposure        0.2559015 0.02219588 0.000000e+00
intercept on Exposure 0.4667074 0.03609327 0.000000e+00
AGE on Exposure       0.5386735 0.05003400 0.000000e+00
GENDER on Exposure    0.5059728 0.02884068 0.000000e+00
PC on Exposure        0.4859232 0.01475425 0.000000e+00
intercept on Outcome  0.5181836 0.04560314 0.000000e+00
AGE on Outcome        0.4645173 0.04889943 0.000000e+00
GENDER on Outcome     0.4578556 0.03803557 0.000000e+00
PC on Outcome         0.4852580 0.03227270 0.000000e+00
Exposure on Outcome   0.2645298 0.06298764 2.672571e-05

> fit$results_var
                                   Estimate
variance of error term (exposure) 1.0101424
variance of error term (outcome)  0.9975058
correlation                       0.1571034
```



## Support
Documentation can also be obtained by running `help` from `R` as following:
```{r}
help(MRIE_EM)
```

If there are any further questions or problems with installing or running `MRIE`, please email me at <yilunli1997@gmail.com>.
