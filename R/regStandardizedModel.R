#' regStandardizedModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' regStandardizedModel creates a regularized model from an mxModel in RAM notation. WARNING: Can only handle models, where
#' the "A" matrix contains directeded paths, the "S" matrix undirected paths and "F" is the filter matrix. Only the " matrix can be regularized.
#' The penalty will be introduced to standardized paths. Standardized paths are calculated unsing the formulas provided in Bollen (1989, p.349)
#'
#' @param mxModelObject an already run mxModel
#' @param regType so far only "lasso" and "ridge" implemented
#' @param regValues numeric value depicting the penalty size
#' @param regIndicator matrix indicating which parameters to regularize in "A"
#'
#' @examples
#' # The following example is adapted from the regsem help to demonstrate the equivalence of both methods:
#'
#' library(lavaan)
#' library(OpenMx)
#' library(regmx)
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
#' mxStandardizeRAMpaths(fit_myModel)
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
#' reg_model <- regStandardizedModel(mxModelObject = fit_myModel, regType = "lasso",regIndicator = selectedA, regValues = regValues)
#' fit_reg_model <- mxRun(reg_model)
#' fit_reg_model$standardizedPaths
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

regStandardizedModel <- function(mxModelObject, regType = "lasso", regIndicator, regValues = 0){
  warning("Work in Progress.")
  ## Checks
  ### Matrices
  if(is.null(mxModelObject$S) | is.null(mxModelObject$A) | is.null(mxModelObject$F)){
    stop("The mxModel object has to contain the A, S, and F matrices.")
  }
  if(!mxModelObject$S$name == "S" && mxModelObject$A$name == "A" && mxModelObject$F$name == "F"){
    stop("The A, S, and F matrices have to be named 'A', 'S', 'F'.")
  }

  ### fitfunction:

  ### model Type
  #if(!class(mxModelObject)[1] == "MxRAMModel"){
  #  stop("Provided mxModelObject is not of type MxRAMModel")
  #}


  ### matrix dimensions
  if(nrow(mxModelObject$A$values) == nrow(regIndicator) &&
     ncol(mxModelObject$A$values) == ncol(regIndicator)){}else{
       stop("Dimensions of Matrix A and provided regIndicator with index do not match.", sep = "")
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

  # create mxMatrix from regValues:
  mxRegValue <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = regValues,name = "regValues") # define peanlty value

  # Define provided mxModelObject as submodel
  Submodel <- mxModel(mxModelObject, name = "Submodel") # has all the parameters and the base fit function (FML or FIML)

  outModel <- mxModel(model= "regmxModel", # model that is returned
                      Submodel, # Submodel is a submodel of outModel. The elements of Submodel can be accessed by the outModel
                      mxNumObs,
                      mxRegValue)

  # Define the new fitting function:

  # Basis: unregularized fitting function from the provided mxModel

  # Parameters need to be standardized:
  ## Compute covariance matrix
  I <- diag(length(c(mxModelObject$manifestVars,mxModelObject$latentVars)))
  Ident <- mxMatrix(type = "Full", values = I, name = "Ident")
  # Covariance-Matrix including latent variances and covariances:
  covMat <- mxAlgebra(solve(Ident-Submodel.A)%*%Submodel.S%*%t(solve(Ident-Submodel.A)), name = "covMat", dimnames = dimnames(Submodel$A))
  standardizer <- mxAlgebra(diag2vec(covMat)^(-1/2)%*%t(diag2vec(covMat)^(1/2)), name = "standardizer", dimnames = dimnames(Submodel$A)) # Die einzelnen Einträge in der A-Matrix müssen zur Standardisierung durch diese Einträge geteilt werden
  standardizedPaths <- mxAlgebra(Submodel.A*standardizer, "standardizedPaths", dimnames = dimnames(Submodel$A))


  if(regType == "lasso"){
    regularizationString <- paste("numObs*(regValues*(t(abs(cvectorize(standardizedPaths))) %*% cvectorize(regIndicator)))", sep = "")
   }else if(
    regType == "ridge"
  ){
     regularizationString <- paste("numObs*(regValues*(t((cvectorize(Submodel.",matrix,"^2))) %*% cvectorize(selected",matrix,"Values)))", sep = "")
  }

  regIndicator <- mxMatrix(type = "Full", values = regIndicator, name = "regIndicator")

  # define the regularized fitfunction:
  fitfun_string <- paste("Submodel.fitfunction",regularizationString, sep = " + ")

  # define fitfunction:
  regFitAlgebra <- mxAlgebraFromString(fitfun_string, name = "regFitAlgebra")
  regFitFunction <- mxFitFunctionAlgebra("regFitAlgebra")

  # complete model
  outModel <- mxModel(outModel,
                      Ident,
                      covMat,
                      standardizer,
                      standardizedPaths,
                      regIndicator,
                      regFitAlgebra,
                      regFitFunction
  )
  outModel <- mxOption(outModel, "Calculate Hessian", "No")
  outModel <- mxOption(outModel, "Standard Errors", "No")

  # return model
  return(outModel)

}
