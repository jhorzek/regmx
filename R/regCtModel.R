regCtModel <- function(ctsemModelObject, regType = "lasso", regOn, regIndicators, regValue = 0, link = list("exp"), dt = list(c(1))){

  call <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Checks
  ### fitfunction:

  ### model Type
  if(!class(ctsemModelObject)[1] == "ctsemFit"){
    Stop("Provided ctsemModelObject is not of type ctsemFit")
  }

  ### Extract mxObject
  mxModelObject <- ctsemModelObject$mxobj
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
    stop("Could not extract the number of observations from the ctsemModelObject provided")
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

    if(link[[matrix]] == "exp"){
      dt_count = 1
      matrix_regstring <- rep(NA, length(dt[[matrix]]))
      for (dt in dt[[matrix]]){
        if(regType == "lasso"){
          matrix_regstring[dt_count] <- paste("numObs*(regValue*(t(abs(cvectorize(expm(Submodel.",regOn[matrix],"*",dt[matrix],")))) %*% cvectorize(",names(mxRegIndicators[matrix]),")))", sep = "")
          }else if(regType == "ridge"){
            matrix_regstring[dt_count] <- paste("numObs*(regValue*(t(cvectorize(expm((Submodel.",regOn[matrix],"*",dt[matrix],")^2))) %*% cvectorize(",names(mxRegIndicators[matrix]),")))", sep = "")
          }
        dt_count <- dt_count+1
      }
    }else if(link[[matrix]] == "QdeltaT"){
      stop("Regularization using QdeltaT not yet implemented.")

    }else{
      # link provided by user as string
      userlink <- link[[matrix]]

      dt_count = 1
      matrix_regstring <- rep(NA, length(dt[[matrix]]))
      for (dt in dt[[matrix]]){
        if(regType == "lasso"){
          matrix_regstring[dt_count] <- paste("numObs*(regValue*(t(abs(cvectorize(",userlink,"(Submodel.",regOn[matrix],"*",dt[matrix],")))) %*% cvectorize(",names(mxRegIndicators[matrix]),")))", sep = "")
          }else if(regType == "ridge"){
            matrix_regstring[dt_count] <- paste("numObs*(regValue*(t(cvectorize(",userlink,"((Submodel.",regOn[matrix],"*",dt[matrix],")^2))) %*% cvectorize(",names(mxRegIndicators[matrix]),")))", sep = "")
      }
      dt_count <- dt_count+1
      }
    }

      comb_matrix_regstring <- paste(matrix_regstring, collapse = "+")

      mxRegFunctions[[matrix]] <- mxAlgebraFromString(algString = comb_matrix_regstring, name = names(mxRegFunctions[matrix]))

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
