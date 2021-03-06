% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regModel.R
\name{regModel}
\alias{regModel}
\title{regModel}
\usage{
regModel(
  mxModelObject,
  alpha = 1,
  gamma = 0,
  regOn,
  regIndicators,
  regValues = 0
)
}
\arguments{
\item{mxModelObject}{an already run mxModel}

\item{alpha}{alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization}

\item{gamma}{gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.}

\item{regOn}{string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")}

\item{regIndicators}{list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter}

\item{regValues}{numeric value depicting the penalty size}
}
\description{
Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
regModel creates a regularized model from an mxModel.
}
\examples{
# The following example is adapted from the regsem help to demonstrate the equivalence of both methods:

library(lavaan)
library(OpenMx)
# put variables on same scale for regsem
HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))

# define variables:
latent = c("f1")
manifest = c("x1","x2","x3","x4","x5", "x6", "x7", "x8", "x9")

# define paths:
loadings <- mxPath(from = latent, to = manifest, free = c(F,T,T,T,T,T,T,T,T), values = 1)
lcov <- mxPath(from = latent, arrows = 2, free = T, values = 1)
lmanif <- mxPath(from = manifest, arrows =2 , free =T, values = 1)

# define model:
myModel <- mxModel(name = "myModel", latentVars = latent, manifestVars = manifest, type = "RAM",
                   mxData(observed = HS, type = "raw"), loadings, lcov, lmanif,
                   mxPath(from = "one", to = manifest, free = T)
)

fit_myModel <- mxRun(myModel)

# Show the names of the matrices in the model:
names(fit_myModel$matrices)
# Show the values of the directional paths:
round(fit_myModel$A$values,5)

# Penalize specific parameters from the A matrix (directional paths):
regOn <- c("A")

selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
selectedA[c(2,3,7,8,9),10] <-1 # parameters that should be regularized have to be marked with 1
regIndicators <- list("A" = selectedA) # save in a list. Note the naming of the list element

# size of the penalty:
regValues = .2

# implement lasso regularization:
reg_model <- regModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = c("A"), regIndicators = regIndicators, regValues = regValues)
fit_reg_model <- mxRun(reg_model)

# extract the A matrix
round(fit_reg_model$Submodel$A$values,5) # Note: the values are stored in the Submodel
# Compare to unregularized parameter values:
round(fit_myModel$A$values,5)

}
\author{
Jannik Orzek
}
