#' regCtModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' regModel creates a regularized model from a ctsem.
#'
#' @param ctsemModelObject an already run ctsem object
#' @param alpha alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization
#' @param gamma gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.
#' @param regValues numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#' @param link list with functions that will be applied to the regularized matrices. For example, if the discrete time autoregressive and cross-lagged parameters should be regularized, use list("DRIFT" = "expm"). If link = "ident" is provided, the ct parameters will be regularized directly.
#' @param dt list with vectors depicting the time points for which the discrete time parameters should be regularized (e.g., list("DRIFT" = c(1)))
#'
#'
#' @examples
#' library(ctsem)
#'
#' set.seed(175446)
#'
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                          MANIFESTVAR=diag(0,2),
#'                          LAMBDA=diag(1,2),
#'                          DRIFT=ct_drift,
#'                          DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                          CINT=matrix(c(0,0),nrow=2),
#'                          T0MEANS=matrix(0,ncol=1,nrow=2),
#'                          T0VAR=diag(1,2))
#'
#' # simulate a training data set
#' traindata <- ctGenerate(generatingModel,n.subjects = 100)
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#'
#' myModel <- ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                    LAMBDA=diag(1,2),
#'                    MANIFESTVAR=diag(0,2),
#'                    CINT=matrix(c(0,0),nrow=2),
#'                    DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                    T0MEANS=matrix(0,ncol=1,nrow=2),
#'                    T0VAR=diag(1,2))
#'
#' # fit the model using ctsem:
#' fit_myModel <- ctFit(traindata, myModel)
#' # discrete time parameters for deltaT = 1
#' expm(fit_myModel$mxobj$DRIFT$values)
#'
#' # select DRIFT values:
#' regOn = "DRIFT"
#' regIndicators = list("DRIFT" = matrix(c(0,1,1,0), byrow = T, ncol = 2))
#'
#' # regularize discrete time parameters for deltaT = 1
#' link = list("DRIFT" = "expm")
#' dt = list("DRIFT"= 1)
#'
#' # set regValues
#' regValues = .2
#'
#' # build regularized model
#' myRegCtModel <- regCtModel(ctsemModelObject = fit_myModel, alpha = 1, gamma = 0, regOn = regOn, regIndicators = regIndicators, regValues = regValues, link = link, dt = dt)
#' fit_myRegCtModel <- mxRun(myRegCtModel)
#' # extract DRIFT parameters for deltaT = 1
#' expm(fit_myRegCtModel$Submodel$DRIFT$values)
#' @author Jannik Orzek
#' @import OpenMx ctsem
#' @export

regCtModel <- function(ctsemModelObject, alpha = 1, gamma = 0, regOn, regIndicators, regValues = 0, link, dt = NULL, scaleFactors = NULL){

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
    if(!matrixName %in% names(regIndicators)){
      stop(paste("Matrix ", matrixName, " provided in regOn was not found in the regIndicator list.", sep = ""))
    }
  }

  ### matrix dimensions
  for(matrix in regOn){
    if(!is.null(OpenMx::imxExtractSlot(mxModelObject[[matrix]], "values"))){
      tempMat <- imxExtractSlot(mxModelObject[[matrix]], "values")
    }else if(!is.null(OpenMx::imxExtractSlot(mxModelObject[[matrix]], "result"))){
      tempMat <- imxExtractSlot(mxModelObject[[matrix]], "result")
    }else{next}
    if(nrow(tempMat) == nrow(regIndicators[[matrix]]) &
       ncol(tempMat) == ncol(regIndicators[[matrix]])){}else{
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
    if(!is.null(OpenMx::imxExtractSlot(mxModelObject[[matrix]], "values"))){
      MLEEstimates[[paste("MLE",matrix, "Estimate", sep ="")]] <- mxMatrix(type = "Full", values = imxExtractSlot(mxModelObject[[matrix]], "values"), free = F, name =names(MLEEstimates[paste("MLE",matrix, "Estimate", sep ="")]))
      }else if(!is.null(OpenMx::imxExtractSlot(mxModelObject[[matrix]], "result"))){
        MLEEstimates[[paste("MLE",matrix, "Estimate", sep ="")]] <- mxMatrix(type = "Full", values = imxExtractSlot(mxModelObject[[matrix]], "result"), free = F, name =names(MLEEstimates[paste("MLE",matrix, "Estimate", sep ="")]))
      }else{stop("Could not extract values from matrix ", regOn[[matrix]], " is matrix ", matrix, " an mxMatrix or mxAlgebra?", sep = "")}

    # create mxMatrix from regValues:
    if(is.list(regValues)){
      mxRegValues[[paste("regValues",matrix, sep ="")]] <- mxMatrix(type = "Full", values = regValues[[matrix]], free = F, nrow = 1, ncol = 1, name = names(mxRegValues[paste("regValues",matrix, sep ="")]))
    }else{
      mxRegValues[[paste("regValues",matrix, sep ="")]] <- mxMatrix(type = "Full", values = regValues, free = F,nrow = 1, ncol = 1, name = names(mxRegValues[paste("regValues",matrix, sep ="")]))
    }

    # create mxAlgebra:

    mxRegIndicators[[paste("selected",matrix, "Values", sep ="")]] <- mxMatrix(type = "Full", values = regIndicators[[matrix]], free = F, name = names(mxRegIndicators[paste("selected",matrix, "Values", sep ="")]))

    if(link[[matrix]] == "expm"){
      dt_count = 1
      matrix_regstring <- rep(NA, length(t))
      for (t in dt){
          matrix_regstring[dt_count] <- paste("numObs*(regValues",matrix,"*((1-",alpha,")*sum(omxSelectRows(cvectorize(expm(Submodel.",
                                              matrix,"*",t,")^2), cvectorize(selected",
                                              matrix,"Values)))+",alpha,"*(sum(omxSelectRows(cvectorize(abs(expm(MLE",
                                              matrix,"Estimate *",t,")^(-",gamma,"))), cvectorize(selected",
                                              matrix,"Values)) * omxSelectRows(cvectorize(abs(expm(Submodel.",
                                              matrix,"*",t,"))), cvectorize(selected",matrix,"Values))))))", sep = "")
        dt_count <- dt_count+1
      }
    }else if(link[[matrix]] == "QdeltaT"){
      stop("Regularization using QdeltaT not yet implemented.")

    }else if(link[[matrix]] == "ident" && is.null(scaleFactors)){
      # regularize ct parameter
          matrix_regstring <- paste("numObs*(regValues",matrix,"*((1-",alpha,")*sum(omxSelectRows(cvectorize((Submodel.",
                                    matrix,")^2), cvectorize(selected",
                                    matrix,"Values)))+",alpha,"*(sum(omxSelectRows(cvectorize(abs(MLE",
                                    matrix,"Estimate^(-",gamma,"))), cvectorize(selected",
                                    matrix,"Values)) * omxSelectRows(cvectorize(abs(Submodel.",
                                    matrix,")), cvectorize(selected",matrix,"Values))))))", sep = "")
    }else if(link[[matrix]] == "ident" && !is.null(scaleFactors)){
      warning("Experimental! Only implemented for penalty on DRIFT! See Harrell 2015, p. 209. Function ignores gamma, i.e. adaptive lasso is not yet implemented!")
      if(!(matrix == "DRIFT")){
        stop("only implemented for drift matrix")
      }
      if(scaleFactors == "autoScaleDRIFT"){
        fullScaleFactorMatrixString <- "(diag2vec(Submodel.T0VAR)^(-1))%*%t(diag2vec(Submodel.T0VAR))"
        selectedScaleFactorString <- paste("omxSelectRows(cvectorize(",fullScaleFactorMatrixString, "), cvectorize(selected",
                                           matrix,"Values))", sep = "")
      }
      # regularize ct parameter
      matrix_regstring <- paste("numObs*(regValues",matrix,"*((1-",alpha,")*t(",selectedScaleFactorString,")%*%(omxSelectRows(cvectorize((Submodel.",
                                matrix,")^2), cvectorize(selected",
                                matrix,"Values)))+",
                                alpha,"*(t(",selectedScaleFactorString,")%*% omxSelectRows(cvectorize(abs(Submodel.",
                                matrix,")), cvectorize(selected",matrix,"Values)))))", sep = "")
    }

      comb_matrix_regstring <- paste(matrix_regstring, collapse = "+")

      mxRegFunctions[[paste("penaltyOn",matrix, sep ="")]] <- mxAlgebraFromString(algString = comb_matrix_regstring, name = paste("penaltyOn",matrix, sep =""))

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
