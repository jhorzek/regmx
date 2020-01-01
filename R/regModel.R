#' regModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' regModel creates a regularized model from an mxModel.
#'
#' @param mxModelObject an already run mxModel
#' @param alpha alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization
#' @param gamma gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.
#' @param regValues numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#'
#' @examples
#' # The following example is adapted from the regsem help to demonstrate the equivalence of both methods:
#'
#' library(lavaan)
#' library(OpenMx)
#' # put variables on same scale for regsem
#' HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
#'
#' # define variables:
#' latent = c("f1")
#' manifest = c("x1","x2","x3","x4","x5", "x6", "x7", "x8", "x9")
#'
#' # define paths:
#' loadings <- mxPath(from = latent, to = manifest, free = c(F,T,T,T,T,T,T,T,T), values = 1)
#' lcov <- mxPath(from = latent, arrows = 2, free = T, values = 1)
#' lmanif <- mxPath(from = manifest, arrows =2 , free =T, values = 1)
#'
#' # define model:
#' myModel <- mxModel(name = "myModel", latentVars = latent, manifestVars = manifest, type = "RAM",
#'                    mxData(observed = HS, type = "raw"), loadings, lcov, lmanif,
#'                    mxPath(from = "one", to = manifest, free = T)
#' )
#'
#' fit_myModel <- mxRun(myModel)
#'
#' # Show the names of the matrices in the model:
#' names(fit_myModel$matrices)
#' # Show the values of the directional paths:
#' round(fit_myModel$A$values,5)
#'
#' # Penalize specific parameters from the A matrix (directional paths):
#' regOn <- c("A")
#'
#' selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
#' selectedA[c(2,3,7,8,9),10] <-1 # parameters that should be regularized have to be marked with 1
#' regIndicators <- list("A" = selectedA) # save in a list. Note the naming of the list element
#'
#' # size of the penalty:
#' regValues = .2
#'
#' # implement lasso regularization:
#' reg_model <- regModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = c("A"), regIndicators = regIndicators, regValues = regValues)
#' fit_reg_model <- mxRun(reg_model)
#'
#' # extract the A matrix
#' round(fit_reg_model$Submodel$A$values,5) # Note: the values are stored in the Submodel
#' # Compare to unregularized parameter values:
#' round(fit_myModel$A$values,5)
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

regModel <- function(mxModelObject, alpha = 1, gamma = 0, regOn, regIndicators, regValues = 0){

  ## Checks

  ### regOn
  checkRegularizedMatrixExistance(regOn = regOn, mxModelObject = mxModelObject, regIndicators = regIndicators)

  ### matrix dimensions
  for(matrix in regOn){
    if(nrow(mxModelObject[[matrix]]$values) == nrow(regIndicators[[matrix]]) &
       ncol(mxModelObject[[matrix]]$values) == ncol(regIndicators[[matrix]])){}else{
         stop("Dimensions of Matrix ", regOn[[matrix]], " and provided regIndicator with index ", matrix, " do not match.", sep = "")
       }
  }

  ## get number of observations:
  if(mxModelObject$data$type == "raw"){
    numObs <- nrow(mxModelObject$data$observed)
  }else if (mxModelObject$data$type == "cov"){
    numObs <- mxModelObject$data$numObs
  }else{
    stop("Could not extract the number of observations from the mxModelObject provided")
  }

  mxNumObs <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = numObs,name = "numObs") # define numObs as mxMatrix

  # Define provided mxModelObject as submodel
  Submodel <- mxModel(mxModelObject, name = "Submodel") # has all the parameters and the base fit function (FML or FIML)

  outModel <- mxModel(model= "regmxModel", # model that is returned
                      Submodel, # BaseModel is a submodel of outModel. The elements of BaseModel can be accessed by the outModel
                      mxNumObs)

  # Define the new fitting function:

  # Basis: unregularized fitting function from the provided mxModel

  fitfun_string <- "Submodel.fitfunction"

  # list of mxMatrices:
  mxRegIndicators <- vector("list", length = length(regOn))
  names(mxRegIndicators) <- paste("selected",regOn, "Values", sep ="")
  mxRegFunctions <- vector("list", length = length(regOn))
  names(mxRegFunctions) <- paste("penaltyOn",regOn, sep ="")
  mxRegValues <- vector("list", length = length(regOn))
  names(mxRegValues) <- paste("regValues",regOn, sep ="")
  MLEEstimates <- vector("list", length = length(regOn))
  names(MLEEstimates) <- paste("MLE",regOn, "Estimate", sep ="")

  # iterate through the matrices that should be regularized:
  for (matrix in regOn){
    # save MLE Estimates for adaptive LASSO
    MLEEstimates[[paste("MLE",matrix, "Estimate", sep ="")]] <- mxMatrix(type = "Full", values = mxModelObject[[matrix]]$values, free = F, name =names(MLEEstimates[paste("MLE",matrix, "Estimate", sep ="")]))

    mxRegIndicators[[paste("selected",matrix, "Values", sep ="")]] <- mxMatrix(type = "Full", values = regIndicators[[matrix]], free = F, name =names(mxRegIndicators[paste("selected",matrix, "Values", sep ="")]))

    # create mxMatrix from regValues:
    if(is.list(regValues)){
      mxRegValues[[paste("regValues",matrix, sep ="")]] <- mxMatrix(type = "Full", values = regValues[[matrix]], free = F, nrow = 1, ncol = 1, name = names(mxRegValues[paste("regValues",matrix, sep ="")]))
    }else{
      mxRegValues[[paste("regValues",matrix, sep ="")]] <- mxMatrix(type = "Full", values = regValues, free = F,nrow = 1, ncol = 1, name = names(mxRegValues[paste("regValues",matrix, sep ="")]))
    }

    # create mxAlgebra:
    mxRegIndicators[[paste("selected",matrix, "Values", sep ="")]] <- mxMatrix(type = "Full", values = regIndicators[[matrix]], free = F, name = names(mxRegIndicators[paste("selected",matrix, "Values", sep ="")]))

    regularizationString <- paste("numObs*(regValues",matrix,"*((1-",alpha,")*sum(omxSelectRows(cvectorize((Submodel.",
                                  matrix,")^2), cvectorize(selected",
                                  matrix,"Values)))+",alpha,"*(sum(omxSelectRows(cvectorize(abs(MLE",
                                  matrix,"Estimate^(-",gamma,"))), cvectorize(selected",
                                  matrix,"Values)) * omxSelectRows(cvectorize(abs(Submodel.",
                                  matrix,")), cvectorize(selected",matrix,"Values))))))", sep = "")

    mxRegFunctions[[paste("penaltyOn",matrix, sep ="")]] <- mxAlgebraFromString(algString = regularizationString, name = paste("penaltyOn",matrix, sep =""))

    # Add mxRegIndicator and mxRegFunction to the model:
    outModel <- mxModel(outModel,
                        mxRegIndicators[[paste("selected",matrix, "Values", sep ="")]],
                        mxRegFunctions[[paste("penaltyOn",matrix, sep ="")]],
                        MLEEstimates[[paste("MLE",matrix, "Estimate", sep ="")]],
                        mxRegValues[[paste("regValues",matrix, sep ="")]]
    )

    # expand the fitting function:
    fitfun_string <- paste(fitfun_string,names(mxRegFunctions[paste("penaltyOn",matrix, sep ="")]), sep = " + ")

  }

  # define fitfunction:
  regFitAlgebra <- mxAlgebraFromString(fitfun_string, name = "regFitAlgebra")
  regFitFunction <- mxFitFunctionAlgebra("regFitAlgebra")

  # complete model
  outModel <- mxModel(outModel,
                      regFitAlgebra,
                      regFitFunction
  )
  outModel <- mxOption(outModel, "Calculate Hessian", "No")
  outModel <- mxOption(outModel, "Standard Errors", "No")

  # return model
  return(outModel)

}
