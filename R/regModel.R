#' regModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' regModel creates a regularized model from an mxModel.
#'
#' @param mxModelObject run mxModel
#' @param regType so far only "lasso" and "ridge" implemented
#' @param regValue numeric value of penalty size
#' @param regOn string vector with matrices that should be regularized.
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. Has to be in the same order as regOn. 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#'
#' @examples
#' # The following example is taken from the regsem help to demonstrate the equivalence of both methods:
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
#' round(fit_myModel$A$values,5)
#'
#' # create regularized model:
#'
#' selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
#' selectedA[c(2,3,7,8,9),10] <-1
#' regIndicators <- list("selectedA" = selectedA)
#'
#'
#' reg_model <- regModel(mxModelObject = fit_myModel, regType = "lasso", regOn  = "A", regValue = .05
#' )
#' fit_reg_model <- mxRun(reg_model)
#'
#' round(fit_reg_model$BaseModel$A$values,5)
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

regModel <- function(mxModelObject, regType = "lasso", regOn, regIndicators, regValue = 0){

  ## Checks
  ### fitfunction:

  ### model Type
  if(!class(mxModelObject)[1] == "MxRAMModel"){
    Stop("Provided mxModelObject is not of type MxRAMModel")
  }
  ### regOn
  for(matrixName in regOn) {
    if(!matrixName %in% names(mxModelObject)){
      stop(paste("Matrix ", matrixName, " provided in regOn was not found in the provided regModel", sep = ""))
    }
  }

  ### matrix dimensions
  for(matrix in 1:length(regOn)){
    if(nrow(mxModelObject[[regOn[matrix]]]$values) == nrow(regIndicators[[matrix]]) &
       ncol(mxModelObject[[regOn[matrix]]]$values) == ncol(regIndicators[[matrix]])){}else{
         stop("Dimensions of Matrix ", regOn[matrix], " and provided regIndicator with index ", matrix, " do not match.", sep = "")
       }
  }

  ## get number of observations:
    if(mxModelObject$data$type == "raw")
      numObs <- nrow(mxModelObject$data$observed)
    else if (mxModelObject$data$type == "cov"){
      numObs <- mxModelObject$data$numObs
    }else{
      stop("Could not extract the number of observations from the mxModelObject provided")
    }

    mxNumObs <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = numObs,name = "numObs") # define numObs as mxMatrix

    # create mxMatrix from regValue:
    mxRegValue <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = regValue,name = "regValue") # define peanlty value

    # Define provided mxModelObject as submodel
    Submodel <- mxModel(mxModelObject, name = "Submodel") # has all the parameters and the base fit function (FML or FIML)

    outModel <- mxModel(model= "regmxModel", # model that is returned
                        Submodel, # BaseModel is a submodel of outModel. The elements of BaseModel can be accessed by the outModel
                        mxNumObs,
                        mxRegValue)

    # Define the new fitting function:

    # Basis: unregularized fitting function from the provided mxModel

    fitfun_string <- "Submodel.fitfunction"

    # list of mxMatrices:
    mxRegIndicators <- vector("list", length = length(regOn))
    names(mxRegIndicators) <- paste("selected",regOn, "Values", sep ="")
    mxRegFunctions <- vector("list", length = length(regOn))
    names(mxRegFunctions) <- paste("penaltyOn",regOn, sep ="")

    # iterate through the matrices that should be regularized:
    for (matrix in 1:length(regOn)){

      # create mxAlgebra:
      mxRegIndicators[[matrix]] <- mxMatrix(type = "Full", values = regIndicators[[matrix]], free = F, name =names(mxRegIndicators[matrix]))

      if(regType == "lasso"){
      regularizationString <- paste("numObs*(regValue*(t(abs(cvectorize(Submodel.",regOn[matrix],"))) %*% cvectorize(",names(mxRegIndicators[matrix]),")))", sep = "")}else if(
        regType == "ridge"
      ){
        regularizationString <- paste("numObs*(regValue*(t((cvectorize(Submodel.",regOn[matrix],"^2))) %*% cvectorize(",names(mxRegIndicators[matrix]),")))", sep = "")
      }
      mxRegFunctions[[matrix]] <- mxAlgebraFromString(algString = regularizationString, name = names(mxRegFunctions[matrix]))

      # Add mxRegIndicator and mxRegFunction to the model:
      outModel <- mxModel(outModel,
                          mxRegIndicators[[matrix]],
                          mxRegFunctions[[matrix]]
      )

      # expand the fitting function:
      fitfun_string <- paste(fitfun_string,names(mxRegFunctions[matrix]), sep = " + ")

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
