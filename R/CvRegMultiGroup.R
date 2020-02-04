#' CvRegMultiGroup
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#'
#' computes cross-validated multi-group model
#' @author Jannik Orzek
#' @import OpenMx
#' @examples
#'
#' @export
#'


CvRegMultiGroup <- function(mxModelObjects, regOn, regType = "lasso", regIndicators, regValue_start = 0, regValue_end = 1,
                            regValue_stepsize = .01,autoCV = FALSE, k = 5, manualCV = NULL, zeroThresh = .001, scaleCV = TRUE){

  stop("not finished")

  # save call:
  call <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Checks
  if(!autoCV & is.null(manualCV)){
    stop("Currently only working with cross-validation. Either provide a validation set using manualCV or set autoCV to TRUE")
  }

  modelNames <- c()

  for (model in mxModelObjects){
    modelNames <- c(modelNames, model$name)
  }

  if(!is.null(manualCV)){
    for(modeName in modelNames){
      if(!modelName %in% names(manualCV)){
        stop(paste("ModelName ", modelName, " provided in mxModelObjects was not found in the provided manualCV data set list", sep = ""))
      }
    }
  }
  # create regValues:
  regValues = seq(from = regValue_start, to = regValue_end, by = regValue_stepsize)

  if(!autoCV){

    results <- matrix(NA, nrow = 4, ncol = length(regValues))
    rownames(results) <- c("penalty", "CV.m2LL", "negative variances","convergence")

    counter <- 1

    # create progress bar:
    pb <- txtProgressBar(min = regValue_start, max = regValue_end, style = 3)

    # iterate over regValues:
    for(regValue in regValues){
      results["penalty",counter] <- regValue
      reg_Model <- regMultiGroup(mxModelObjects = mxModelObjects, regType = regType, regOn = regOn,
                                 regIndicators = regIndicators, regValue = regValue)

      reg_Model <- mxOption(reg_Model, "Calculate Hessian", "No") # might cause errors; check
      reg_Model <- mxOption(reg_Model, "Standard Errors", "No") # might cause errors; check
      fit_reg_Model <- mxRun(reg_Model, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

      results["convergence",counter] <- fit_reg_Model$output$status$code# check convergence

      if("S" %in% names(fit_reg_Model$Submodel)){
        variances = diag(nrow(fit_reg_Model$Submodel$S$values))==1

        if(any(fit_reg_Model$Submodel$S$values[variances] <0)){
          results["negative variances",counter] <- 1 # check negative variances
        }else(
          results["negative variances",counter] <- 0
        )}

      # Compute CV fit:

      sumCvM2LL <- 0
      for (submodel in fit_reg_Model$submodels){
        ExpCov <- mxGetExpected(submodel, component = "covariance")
        ExpMean <- mxGetExpected(submodel, component = "means")

        ObsData <- manualCV[[submodel$name]]$observed

        CvM2LL <- getM2LogL_FIML(ExpCov = ExpCov, ExpMean = ExpMean, ObsData = ObsData)
        sumCvM2LL <- sumCvM2LL + CvM2LL
      }

      results["CV.m2LL", counter] <- sumCvM2LL

    }

    # find minimal CV m2LL:

    minimum_CV.m2LL <- convergeSubset[1,which(convergeSubset["CV.m2LL",]==min(convergeSubset["CV.m2LL",]))]

    best_penalty = minimum_CV.m2LL
    reg_Model_CVm2LL <- regMultiGroup(mxModelObjects = mxModelObjects, regType = regType, regOn = regOn,
                                      regIndicators = regIndicators, regValue = best_penalty)
    reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Calculate Hessian", "No") # might cause errors; check
    reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Standard Errors", "No") # might cause errors; check

    fit_reg_Model_CVm2LL <- mxRun(reg_Model_CVm2LL, silent = T)
    out <- list("best penalty" = reg_Model_CVm2LL, "bestmodel" = fit_reg_Model_CVm2LL, "fit measures" = t(results), "call" = call)

    class(out) <- "CvRegMultiGroupObject"

    return(out)



  }else if(autoCV){
    # create List of folds

    FoldsList <- vector("list", length = length(modelNames))
    names(FoldsList) <- modelNames

    for (model in mxModelObjects){
      for(modelName in modelNames){
        subjects <- 1:nrow(model$data$observed)
        Folds <- split(sample(subjects, length(subjects),replace=FALSE), f = as.factor(1:k))
        FoldsList[[modelName]] <- Folds
      }
    }

    rawDataList <- vector("list", length = length(modelNames))
    names(rawDataList) <- modelNames

    for (model in mxModelObjects){
      for(modelName in modelNames){
        subjects <- 1:nrow(model$data$observed)
        Folds <- split(sample(subjects, length(subjects),replace=FALSE), f = as.factor(1:k))
        FoldsList[[modelName]] <- Folds
      }
    }
    ################################ LAST UPDATE

    full_raw_data <- mxModelObject$data$observed

    Res <- matrix(NA, nrow = length(seq(from = regValue_start, to = regValue_end, by = regValue_stepsize)), ncol = k+4)
    colnames(Res) <- c("penalty", "mean CV/Boot_m2LL",paste("fold", 1:k), "negative variances", "convergence problems")
    Res[,"penalty"] <- seq(from = regValue_start, to = regValue_end, by = regValue_stepsize)
    Res[,"negative variances"] <- 0
    Res[,"convergence problems"] <- 0


  }

}
