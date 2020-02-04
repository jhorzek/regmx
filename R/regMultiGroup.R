#' regMultiGroup
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' regMultiGroup a mutliple Group Model with a penalty on the difference in parameter estimates between the groups. Its is not tested and should not be used!
#'
#' @param mxModelObjects a list with mxModels for the different groups
#' @param regType so far only "lasso" implemented
#' @param regValue numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObjects (e.g., "A").
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#'
#' @examples
#' # The following example is adapted from the OpenMx documentation
#' library(OpenMx)
#'
#' # Simulate some data
#'
#' # Group 1
#' N1 = 100
#' x = rnorm(N1, mean= 0, sd= 1)
#' y = 0.5*x + rnorm(N1, mean= 0, sd= 1)
#' ds1 <- data.frame(x, y)
#' dsNames <- names(ds1)
#'
#' # Group 2: higher regression weight
#' N2= 1000
#' x= rnorm(N2, mean= 0, sd= 1)
#' y= 0.6*x + rnorm(N2, mean= 0, sd= 1)
#' ds2 <- data.frame(x, y)
#'
#'
#' # Define the matrices (A matrix implementation of 2 RAM models)
#'
#' # I is identical for both groups:
#' I <- mxMatrix(name="I", type="Iden", nrow=2, ncol=2)
#' # means are identical for both groups:
#' M <- mxMatrix(name = "M", type = "Full", nrow = 1, ncol = 2, values=0,
#'               free=TRUE, labels=c("Mean_x", "Mean_y"))
#' # The A matrix can differ:
#' A1 <- mxMatrix(name = "A", type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0),
#'                free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b1", NA, NA))
#' A2 <- mxMatrix(name = "A", type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0),
#'                free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b2", NA, NA))
#' # S is identical:
#' S <- mxMatrix(name = "S", type = "Diag", nrow = 2, ncol = 2, values=1,
#'               free=TRUE, labels=c("Var_x", "Resid"))
#'
#'
#' # Define the expectation
#' expect <- mxExpectationRAM('A', 'S', 'I', 'M', dimnames= dsNames)
#'
#' # Choose a fit function
#' fitFunction <- mxFitFunctionML(rowDiagnostics=TRUE)
#' # Also return row likelihoods (the fit function value is still 1x1)
#'
#' # Multiple-group fit function sums the model likelihoods
#' # from its component models
#' mgFitFun <- mxFitFunctionMultigroup(c('g1model', 'g2model'))
#'
#'
#' # Define model 1 and model 2
#' m1 = mxModel(model="g1model",
#'              M, S, A1, I, expect, fitFunction,
#'              mxData(cov(ds1), type="cov", numObs=N1, means=colMeans(ds1))
#' )
#' m2 = mxModel(model="g2model",
#'              M, S, A2, I, expect, fitFunction,
#'              mxData(cov(ds2), type="cov", numObs=N2, means=colMeans(ds2))
#' )
#'
#' mg <- mxModel(model='multipleGroup', m1, m2, mgFitFun)
#' # All paths except for b are constrained to equality
#'
#' # Fit the model and print a summary
#' mg <- mxRun(mg)
#' summary(mg)
#'
#' ############### Regularize A matrix:
#'
#' regIndicators <- list("A" = matrix(c(0,1,1,0), nrow = 2))
#'
#' regValue <- 1
#'
#' regOn <- "A"
#'
#' mxModelObjects <- list("model1" = m1, "model2" = m2)
#' regModel <- regMultiGroup(mxModelObjects = mxModelObjects, regOn = regOn, regIndicators = regIndicators, regType = "lasso",regValue = regValue )
#' fit_regModel <- mxRun(regModel)
#' summary(fit_regModel)
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export
regMultiGroup <- function(mxModelObjects, regOn, regType = "lasso", regIndicators, regValue = 0){
  warning("This function has not been tested! It is experimental and I recommend against its use")
  ### regOn
  for(matrixName in regOn) {
    for (model in mxModelObjects){
      if(!matrixName %in% names(model)){
        stop(paste("Matrix ", matrixName, " provided in regOn was not found in the provided Model", sep = ""))
      }
    }
    if(!matrixName %in% names(regIndicators)){
      stop(paste("Matrix ", matrixName, " provided in regOn was not found in the regIndicator list.", sep = ""))
    }
  }

  ### matrix dimensions
  for(matrix in regOn){
    for (model in mxModelObjects){
      if(nrow(model[[matrix]]$values) == nrow(regIndicators[[matrix]]) &
         ncol(model[[matrix]]$values) == ncol(regIndicators[[matrix]])){}else{
           stop("Dimensions of Matrix ", regOn[[matrix]], " and provided regIndicator with index ", matrix, " do not match.", sep = "")
         }
    }
  }

  ## get number of observations:
  numObsTotal <- 0
  for (model in mxModelObjects){
    if(model$data$type == "raw"){
      numObsTotal <- numObsTotal + nrow(model$data$observed)
    }else if (model$data$type == "cov"){
      numObsTotal <- numObsTotal + model$data$numObs
    }else{
      stop("Could not extract the number of observations from the mxModelObject provided")
    }
  }

  mxNumObsTotal <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = numObsTotal,name = "mxNumObsTotal") # define numObs as mxMatrix

  # create mxMatrix from regValue:
  mxRegValue <- mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = regValue,name = "regValue") # define peanlty value

  # create regMultiGroupModel
  regMultiGroupModel <- mxModel(model= "regMultiGroupModel", # model that is returned
                                mxNumObsTotal,
                                mxRegValue)

  # add the mxModelObects
  for (model in mxModelObjects){
    regMultiGroupModel <- mxModel(
      regMultiGroupModel,
      model # the models are now submodels of reMultiGroupModel
    )
  }

  # Define the new fitting function:

  # Basis: unregularized fitting function from the provided mxModels
  # get modelnames
  modelNames <- c()

  for (model in mxModelObjects){
    modelNames <- c(modelNames, model$name)
  }

  fitfun_string <- paste(paste(modelNames, ".fitfunction", sep =""), collapse = " + ")

  # list of mxMatrices:
  mxRegIndicators <- vector("list", length = length(regOn))
  names(mxRegIndicators) <- paste("selected",regOn, "Values", sep ="")
  mxRegFunctions <- vector("list", length = length(regOn))
  names(mxRegFunctions) <- paste("penaltyOn",regOn, sep ="")

  # iterate through the matrices that should be regularized:
  for (matrix in regOn){

    # create mxAlgebra:
    mxRegIndicators[[paste("selected",matrix, "Values", sep ="")]] <- mxMatrix(type = "Full", values = regIndicators[[matrix]], free = F, name =names(mxRegIndicators[paste("selected",matrix, "Values", sep ="")]))

    if(regType == "lasso"){
      regularizationString_1 <- c("mxNumObsTotal*regValue*(t(abs(cvectorize(")
      regularizationString_2 <- c()
      for (model in mxModelObjects){
        regularizationString_2 <- c(regularizationString_2, paste(model$name, matrix, sep = "."))
      }
      regularizationString_2 <- paste(regularizationString_2, collapse = " - ")
      regularizationString_3 <- paste(")))%*%cvectorize(selected",matrix, "Values))", sep ="")

      regularizationString <- paste(regularizationString_1, regularizationString_2, regularizationString_3, collapse = "")

    }else if(
      regType == "ridge"
    ){
      stop("not yet implemented")
    }
    mxRegFunctions[[paste("penaltyOn",matrix, sep ="")]] <- mxAlgebraFromString(algString = regularizationString, name = paste("penaltyOn",matrix, sep =""))

    # Add mxRegIndicator and mxRegFunction to the model:
    regMultiGroupModel <- mxModel(regMultiGroupModel,
                                  mxRegIndicators[[paste("selected",matrix, "Values", sep ="")]],
                                  mxRegFunctions[[paste("penaltyOn",matrix, sep ="")]]
    )

    # expand the fitting function:
    fitfun_string <- paste(fitfun_string,names(mxRegFunctions[paste("penaltyOn",matrix, sep ="")]), sep = " + ")

  }

  # define fitfunction:
  regFitAlgebra <- mxAlgebraFromString(fitfun_string, name = "regFitAlgebra")
  regFitFunction <- mxFitFunctionAlgebra("regFitAlgebra")

  # complete model
  regMultiGroupModel <- mxModel(regMultiGroupModel,
                                regFitAlgebra,
                                regFitFunction
  )
  regMultiGroupModel <- mxOption(regMultiGroupModel, "Calculate Hessian", "No")
  regMultiGroupModel <- mxOption(regMultiGroupModel, "Standard Errors", "No")

  # return model
  return(regMultiGroupModel)




}
