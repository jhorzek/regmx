#' SingleCoreOptimRegModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' SingleCoreOptimRegModel creates a range of regularized models from an mxModel. It automatically tests different penalty values and returns the best model. It only uses a single core.
#'
#' @param mxModelObject an already run mxModel
#' @param alpha alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization
#' @param gamma gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.
#' @param regValue numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#' @param regValue_start initial penalty value (recommended: 0)
#' @param regValue_end highest penalty value tested
#' @param regValue_stepsize increase in penValue between iterations
#' @param criterion Criterion for chosing the best final model. Possible are: AIC, BIC, CV.m2LL (only if manualCV is provided)
#' @param autoCV logical indicating if cross-validation should be performed automatically
#' @param k number of splits performed when autoCV = TRUE
#' @param Boot logical indicating if Bootstrap should be performed. Not yet implemented!
#' @param manualCV if cross-validation should be performed manually, provide a cross-validation sample (has to be of the same class as the data in the mxModelObject; e.g., mxData)
#' @param zeroThresh Threshold at which parameters are evaluated as being zero (necessary for AIC and BIC)
#' @param scaleCV indicate if the CV samples should be scaled automatically
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
#'
#' # Penalize specific parameters from the A matrix (directional paths):
#' regOn <- c("A")
#'
#' selectedA <- matrix(0, ncol = ncol(fit_myModel$A$values), nrow = nrow(fit_myModel$A$values))
#' selectedA[c(2,3,7,8,9),10] <-1 # parameters that should be regularized have to be marked with 1
#' regIndicators <- list("A" = selectedA) # save in a list. Note the naming of the list element
#'
#' # Run the models:
#'
#' reg_model <- optimRegModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = regOn, regIndicators = regIndicators)
#'
#' reg_model$`fit measures`
#'
#' reg_model$`best penalty`
#'
#' # Run the same model with 5-fold cross-validation
#'
#' CV_reg_model <- optimRegModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = regOn, regIndicators = regIndicators,
#'                               autoCV = T, k = 5)
#' CV_reg_model$`CV results`
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

SingleCoreOptimRegModel <- function(mxModelObject, alpha = 1, gamma = 0, regOn, regIndicators,
                          regValue_start = 0, regValue_end = 1, regValue_stepsize = .01,
                          criterion = "BIC", autoCV = FALSE, k = 5, Boot = FALSE, manualCV = NULL, zeroThresh = .001, scaleCV = TRUE, cores = 1){

  # save call:
  call <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Checks
  if(Boot){
    stop("Bootstrap not yet implemented")
  }
  ### fitfunction:

  ### model Type
  #if(!class(mxModelObject)[1] == "MxRAMModel"){
  #  stop("Provided mxModelObject is not of type MxRAMModel")
  #}
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
    if(nrow(mxModelObject[[matrix]]$values) == nrow(regIndicators[[matrix]]) &
       ncol(mxModelObject[[matrix]]$values) == ncol(regIndicators[[matrix]])){}else{
         stop("Dimensions of Matrix ", regOn[[matrix]], " and provided regIndicator with index ", matrix, " do not match.", sep = "")
       }
  }



  # create regValues:
  regValues = seq(from = regValue_start, to = regValue_end, by = regValue_stepsize)

  if (!autoCV & !Boot){
    results <- matrix(NA, nrow = 8, ncol = length(regValues))
    rownames(results) <- c("penalty", "estimated Parameters", "m2LL","AIC",
                           "BIC", "CV.m2LL", "negative variances","convergence")
    counter <- 1

    # create progress bar:
    pb <- txtProgressBar(min = regValue_start, max = regValue_end, style = 3)

    # iterate over regValues:
    for(regValue in regValues){
      results["penalty",counter] <- regValue
      reg_Model <- regModel(mxModelObject = mxModelObject,
                            alpha = alpha, gamma = gamma, regOn = regOn,
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

      ### compute AIC and BIC:
      FitM <- getFitMeasures(regModel = fit_reg_Model, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, cvSample = manualCV, zeroThresh = zeroThresh)

      results["estimated Parameters",counter] <- FitM$estimated_params # estimated parameters
      results["m2LL",counter] <- FitM$m2LL # -2LogL
      results["AIC",counter] <- FitM$AIC # AIC
      results["BIC",counter] <- FitM$BIC # BIC
      results["CV.m2LL",counter] <- FitM$CV.m2LL # CV.m2LL

      setTxtProgressBar(pb, regValue)
      counter <- counter+1
    }

    # Find Minima / best penalty value
    convergeSubset <- (results["convergence",] == 0) == (results["negative variances",]==0)
    convergeSubset <- results[,convergeSubset]

    minimum_m2LL <- convergeSubset[1,which(convergeSubset["m2LL",]==min(as.numeric(convergeSubset["m2LL",])))]

    minimum_AIC <- convergeSubset[1,which(convergeSubset["AIC",]==min(as.numeric(convergeSubset["AIC",])))]

    minimum_BIC <- convergeSubset[1,which(convergeSubset["BIC",]==min(convergeSubset["BIC",]))]

    minimum_CV.m2LL <- convergeSubset[1,which(convergeSubset["CV.m2LL",]==min(convergeSubset["CV.m2LL",]))]


    # getting parameters:
    if(criterion == "m2LL"){
      best_penalty = minimum_m2LL
      reg_Model_m2LL <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_m2LL <- mxRun(reg_Model_m2LL, silent = T)
      out <- list("best penalty" = minimum_m2LL, "bestmodel" = fit_reg_Model_m2LL, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "AIC"){
      best_penalty = minimum_AIC
      reg_Model_AIC <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_AIC <- mxRun(reg_Model_AIC, silent = T)
      out <- list("best penalty" = minimum_AIC, "bestmodel" = fit_reg_Model_AIC, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "BIC"){
      best_penalty = minimum_BIC
      reg_Model_BIC <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
      out <- list("best penalty" = minimum_BIC, "bestmodel" = fit_reg_Model_BIC, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "CV.m2LL"){
      best_penalty = minimum_CV.m2LL
      reg_Model_CVm2LL <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_CVm2LL <- mxRun(reg_Model_CVm2LL, silent = T)
      out <- list("best penalty" = minimum_CV.m2LL, "bestmodel" = fit_reg_Model_CVm2LL, "fit measures" = t(results), "call" = call)
    }
    class(out) <- "OptimRegModelObject"

    return(out)



  }else if(autoCV | Boot){
    if(autoCV & Boot){
      stop("autoCV and Boot can not be combined.")
    }
    if(!mxModelObject$data$type == "raw"){
      stop("Cross-Validation and Bootstrap only available for raw data")
    }

    if(autoCV){
      subjects <- 1:nrow(mxModelObject$data$observed)
      Folds <- split(sample(subjects, length(subjects),replace=FALSE), f = as.factor(1:k))
    }
    if(Boot){
      stop("Bootstrap not yet implemented")
      subjects <- 1:nrow(mxModelObject$data$observed)
      Folds <- vector("list", length = k)
      for(Fold in 1:k){
        Folds[[Fold]] <- sample(subjects, length(subjects),replace=TRUE)
      }
    }

    full_raw_data <- mxModelObject$data$observed

    Res <- matrix(NA, nrow = length(seq(from = regValue_start, to = regValue_end, by = regValue_stepsize)), ncol = k+4)
    colnames(Res) <- c("penalty", "mean CV/Boot_m2LL",paste("fold", 1:k), "negative variances", "convergence problems")
    Res[,"penalty"] <- seq(from = regValue_start, to = regValue_end, by = regValue_stepsize)
    Res[,"negative variances"] <- 0
    Res[,"convergence problems"] <- 0

    fold <- 1
    for(Fold in Folds){
      Test_Sample <- full_raw_data[Fold,]
      Train_Sample <- full_raw_data[-Fold,]

      if(scaleCV){
        Test_Sample <- scale(Test_Sample)
        Train_Sample <- scale(Train_Sample)
      }
      trainModel <- mxModelObject
      trainModel$data <- mxData(observed = Train_Sample, type = "raw")
      Test_Sample <- mxData(observed = Test_Sample, type = "raw")

      fit_trainModel <- optimRegModel(mxModelObject = trainModel, alpha = alpha, gamma = gamma, regOn = regOn,
                                      regIndicators = regIndicators,
                                      regValue_start = regValue_start,
                                      regValue_end = regValue_end,
                                      regValue_stepsize = regValue_stepsize,
                                      criterion = "CV.m2LL", autoCV = FALSE,
                                      Boot = FALSE, manualCV = Test_Sample,
                                      k = 5, zeroThresh = zeroThresh, scaleCV = scaleCV)

      Res[,paste("fold", fold)] <- fit_trainModel$`fit measures`[,'CV.m2LL']
      Res[,"negative variances"] <- Res[,"negative variances"] + fit_trainModel$`fit measures`[,'negative variances']
      Res[,"convergence problems"] <- Res[,"convergence problems"] + fit_trainModel$`fit measures`[,'convergence']

      cat(paste("\n Completed CV for fold", fold, "of", k, "\n"))

      fold <- fold + 1
    }

    # mean the m2LLs:
    Res[,"mean CV/Boot_m2LL"] <- apply(Res[, paste("fold", 1:k)], 1, mean)

    # only use runs without problems:
    Res_final <- Res[Res[,"negative variances"] == 0 && Res[,"convergence problems"] == 0,]

    # find best penalty value:
    best_penalty <- Res_final[which(Res_final[,"mean CV/Boot_m2LL"] == min(Res_final[,"mean CV/Boot_m2LL"])), "penalty"]

    # fit best penalty model with full data set:
    if(scaleCV){
      scale_full_raw_data <- scale(full_raw_data)
      mxModelObject$data <- mxData(observed = scale_full_raw_data, type = "raw")
    }
    finalModel <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma,
                           regOn = regOn, regIndicators = regIndicators,
                           regValue = best_penalty)

    ffinalModel <- mxRun(finalModel, silent = T)

    ret <- list("CV results" = Res, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k, "call" = call)
    class(ret) <- "CvOptimRegModelObject"

    return(ret)

  }



}
