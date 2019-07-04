optimRegModel <- function(mxModelObject, regType = "lasso", regOn, regIndicators,
                          regValue_start = 0, regValue_end = 1, regValue_stepsize = .01,
                          criterion = "BIC", autoCV = FALSE, Boot = FALSE, manualCV = NULL, k = 5, zeroThresh = .001, scaleCV = TRUE){

  # save call:
  call <- mget(names(formals()),sys.frame(sys.nframe()))

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
                            regType = regType, regOn = regOn,
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
      FitM <- getFitMeasures(regModel = fit_reg_Model, regType = regType, regOn = regOn, regIndicators = regIndicators, cvSample = manualCV, zeroThresh = zeroThresh)

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
      reg_Model_m2LL <- regModel(mxModelObject = mxModelObject, regType = regType, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_m2LL <- mxOption(reg_Model_m2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_m2LL <- mxRun(reg_Model_m2LL, silent = T)
      out <- list("best penalty" = minimum_m2LL, "bestmodel" = fit_reg_Model_m2LL, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "AIC"){
      best_penalty = minimum_AIC
      reg_Model_AIC <- regModel(mxModelObject = mxModelObject, regType = regType, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_AIC <- mxOption(reg_Model_AIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_AIC <- mxRun(reg_Model_AIC, silent = T)
      out <- list("best penalty" = minimum_AIC, "bestmodel" = fit_reg_Model_AIC, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "BIC"){
      best_penalty = minimum_BIC
      reg_Model_BIC <- regModel(mxModelObject = mxModelObject, regType = regType, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_BIC <- mxOption(reg_Model_BIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_BIC <- mxRun(reg_Model_BIC, silent = T)
      out <- list("best penalty" = minimum_BIC, "bestmodel" = fit_reg_Model_BIC, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "CV.m2LL"){
      best_penalty = minimum_BIC
      reg_Model_CVm2LL <- regModel(mxModelObject = mxModelObject, regType = regType, regOn = regOn, regIndicators = regIndicators, regValue = best_penalty)
      reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_Model_CVm2LL <- mxOption(reg_Model_CVm2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_Model_CVm2LL <- mxRun(reg_Model_CVm2LL, silent = T)
      out <- list("best penalty" = reg_Model_CVm2LL, "bestmodel" = fit_reg_Model_CVm2LL, "fit measures" = t(results), "call" = call)
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

      trainModel <- optimRegModel(mxModelObject, regType = regType, regOn = regOn,
                                  regIndicators = regIndicators,
                                  regValue_start = regValue_start,
                                  regValue_end = regValue_end,
                                  regValue_stepsize = regValue_stepsize,
                                  criterion = criterion, autoCV = FALSE,
                                  Boot = FALSE, manualCV = Test_Sample,
                                  k = 5, zeroThresh = zeroThresh, scaleCV = scaleCV)

      Res[,paste("fold", fold)] <- trainModel$`fit measures`[,'CV.m2LL']
      Res[,"negative variances"] <- Res[,"negative variances"] + trainModel$`fit measures`[,'negative variances']
      Res[,"convergence problems"] <- Res[,"convergence problems"] + trainModel$`fit measures`[,'convergence']

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
    finalModel <- regModel(mxModelObject = mxModelObject, regType = regType,
                           regOn = regOn, regIndicators = regIndicators,
                           regValue = best_penalty)

    ffinalModel <- mxRun(finalModel, silent = T)

    ret <- list("CV results" = Res, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k)
    class(ret) <- "CVlaremm"

    return(ret)

  }



}
