#' MultiCoreOptimRegCtModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' MultiCoreOptimRegCtModel creates a range of regularized models from a ctsem. It automatically tests different penalty values and returns the best model. it uses multiple cores.
#'
#' @param ctsemModelObject an already run ctsem object
#' @param regType so far only "lasso" and "ridge" implemented
#' @param regValue numeric value depicting the penalty size
#' @param regOn string vector with matrices that should be regularized. The matrices must have the same name as the ones provided in the mxModelObject (e.g., "A")
#' @param regIndicators list of matrices indicating which parameters to regularize in the matrices provided in regOn. The matrices in regIndicators must to have the same names as the matrices they correspond to (e.g., regIndicators = list("A" = diag(10))). 1 Indicates a parameter that will be regularized, 0 an unregularized parameter
#' @param link list with functions that will be applied to the regularized matrices. For example, if the discrete time autoregressive and cross-lagged parameters should be regularized, use list("DRIFT" = "expm"). You can also define your own link function. DESCRIPTION NOT YET IMPLEMENTED
#' @param dt list with vectors depicting the time points for which the discrete time parameters should be regularized (e.g., list("DRIFT" = c(1)))
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
#' @param cores number of cores
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
#' regIndicators = list("DRIFT" = matrix(c(0,1,0,0), byrow = T, ncol = 2))
#'
#' # regularize discrete time parameters for deltaT = 1
#' link = list("DRIFT" = "expm")
#' dt = list("DRIFT"= 1)
#'
#' # optimize regularized ctsem
#'
#' myRegCtModel <- optimRegCtModel(ctsemModelObject = fit_myModel, regType = "lasso",
#'                                 regOn = regOn, regIndicators = regIndicators,
#'                                 link = link, dt = dt, autoCV = F, criterion = "BIC", cores = 2)
#'
#' # extract discrete time parameters for deltaT = 1
#' expm(myRegCtModel$bestmodel$Submodel$DRIFT$values)
#'
#' # optimize regularized ctsem with cross-validation
#'
#' CV_myRegCtModel <- optimRegCtModel(ctsemModelObject = fit_myModel, regType = "lasso",
#'                                    regOn = regOn, regIndicators = regIndicators,
#'                                    link = link, dt = dt, autoCV = T, k = 5, regValue_start = 0, regValue_end = .5, regValue_stepsize = .01, cores = 2)
#'
#' # extract discrete time parameters for deltaT = 1
#' expm(CV_myRegCtModel$`final Model`$Submodel$DRIFT$values)
#'
#' @author Jannik Orzek
#' @import OpenMx ctsem doParallel foreach iterators parallel
#' @export

MultiCoreOptimRegCtModel <- function(ctsemModelObject, alpha = 1, gamma = 0, regOn, regIndicators,
                                     link = list("exp"), dt,
                                     regValues,
                                     criterion = "BIC", autoCV = FALSE, Boot = FALSE, manualCV = NULL, k = 5, zeroThresh = .001, scaleCV = TRUE, scaleFactors = NULL, cores){

  # save call:
  call <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Checks

  ### Number of Cores:

  numCores <- detectCores()
  if(cores > numCores){
    stop(c("The specified number of cores higher than the number of cores detected by doParallel: Your PC seems to only have ", numCores, " cores"))
  }else{
    if (cores == numCores){
      warning("Are you sure you want to use all available cores for this R session? This might result in errors.")
      continue <- readline(prompt="Please enter 'Y' if you want to continue with all specified cores: ")

      if(!continue == "Y"){
        stop("Programm stopped by user.")
      }
    }
  }

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  if(Boot){
    stop("Bootstrap not yet implemented")
  }
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


  if (!autoCV & !Boot){
    # iterate over regValues:
    RegValuesGrid <- createRegValueGrid(regValues = regValues, regOn = regOn)
    # create list to store results

    results <- list("penalty"= RegValuesGrid, "estimated Parameters"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "m2LL"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),"AIC"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),
                    "BIC"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "CV.m2LL"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "negative variances"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),"convergence"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1))


    # iterate over regValues:
    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }
    ParRes <- foreach(row = 1:nrow(RegValuesGrid), .combine='comb', .multicombine=TRUE,.packages = c("OpenMx","regmx"),
                      .inorder = FALSE,
                      .errorhandling = "stop",
                      .verbose = TRUE) %dopar% {
                        regValue <- RegValuesGrid[row,]
                        reg_ctModel <- regCtModel(ctsemModelObject = ctsemModelObject, link = link, dt = dt,alpha = alpha, gamma = gamma,regValues = regValue,
                                                  regOn = regOn,
                                                  regIndicators = regIndicators,
                                                  scaleFactors = scaleFactors)

                        reg_ctModel <- mxOption(reg_ctModel, "Calculate Hessian", "No") # might cause errors; check
                        reg_ctModel <- mxOption(reg_ctModel, "Standard Errors", "No") # might cause errors; check
                        fit_reg_ctModel <- mxRun(reg_ctModel, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

                        convergence <- fit_reg_ctModel$output$status$code# check convergence

                        if("S" %in% names(fit_reg_ctModel$Submodel)){
                          variances = diag(nrow(fit_reg_ctModel$Submodel$S$values))==1

                          if(any(fit_reg_ctModel$Submodel$S$values[variances] <0)){
                            negative.variances <- 1 # check negative variances
                          }else(
                            negative.variances <- 0
                          )}

                        ### compute AIC and BIC:
                        FitM <- getCtFitMeasures(regCtModel = fit_reg_ctModel, alpha = alpha,
                                                 gamma = gamma, regOn = regOn,
                                                 regIndicators = regIndicators, cvSample = manualCV,
                                                 zeroThresh = zeroThresh)

                        estimated.Parameters <- FitM$estimated_params # estimated parameters
                        m2LL <- FitM$m2LL # -2LogL
                        AIC <- FitM$AIC # AIC
                        BIC <- FitM$BIC # BIC
                        CV.m2LL <- FitM$CV.m2LL # CV.m2LL

                        list("penalty"= regValue, "estimated Parameters"=estimated.Parameters, "m2LL"=m2LL,"AIC"=AIC,
                             "BIC"=BIC, "CV.m2LL"=CV.m2LL, "negative variances"=negative.variances,"convergence"=convergence)
                      }

    results <- ParRes

    # Find Minima / best penalty value
    convergedSubset <- (results[["convergence"]] == 0) & (results[["negative variances"]]==0)

    rowIndicators <- which(results[["m2LL"]][convergedSubset] == min(results[["m2LL"]][convergedSubset]))
    minimum_m2LL <- matrix(results[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
    colnames(minimum_m2LL) = regOn

    rowIndicators <- which(results[["AIC"]][convergedSubset] == min(results[["AIC"]][convergedSubset]))
    minimum_AIC <- matrix(results[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
    colnames(minimum_AIC) = regOn

    rowIndicators <- which(results[["BIC"]][convergedSubset] == min(results[["BIC"]][convergedSubset]))
    minimum_BIC <- matrix(results[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
    colnames(minimum_BIC) = regOn

    rowIndicators <- which(results[["CV.m2LL"]][convergedSubset] == min(results[["CV.m2LL"]][convergedSubset]))
    minimum_CV.m2LL <- matrix(results[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
    colnames(minimum_CV.m2LL) = regOn


    # getting parameters:
    if(criterion == "m2LL"){
      if(nrow(minimum_m2LL)>1){
        warning(paste("Multiple minima for m2LL found. Will use the following penalty values:", paste(as.matrix(minimum_m2LL)[nrow(minimum_m2LL),],collapse = ","), sep = ":"))
        best_penalty = t(as.matrix(as.matrix(minimum_m2LL)[nrow(minimum_m2LL),]))
      }else{
        best_penalty = minimum_m2LL
      }

      if(is.matrix(best_penalty)){
        regValue <- vector("list", length = ncol(best_penalty))
        names(regValue) <- colnames(best_penalty)
        regValue[colnames(best_penalty)] <- best_penalty
      }else{
        regValue <- best_penalty
      }

      reg_CtModel_m2LL <- regCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma,
                                     regOn = regOn, regIndicators = regIndicators,
                                     regValues = regValue, link = link, dt = dt,
                                     scaleFactors = scaleFactors)
      reg_CtModel_m2LL <- mxOption(reg_CtModel_m2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_m2LL <- mxOption(reg_CtModel_m2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_m2LL <- mxRun(reg_CtModel_m2LL, silent = T)
      out <- list("best penalty" = minimum_m2LL, "bestmodel" = fit_reg_CtModel_m2LL, "fit measures" = results, "call" = call)
    }

    if(criterion == "AIC"){
      if(nrow(minimum_AIC)>1){
        warning(paste("Multiple minima for AIC found. Will use the following penalty values:", paste(as.matrix(minimum_AIC)[nrow(minimum_AIC),],collapse = ","), sep = ":"))
        best_penalty = t(as.matrix(as.matrix(minimum_AIC)[nrow(minimum_AIC),]))
      }else{
        best_penalty = minimum_AIC
      }

      if(is.matrix(best_penalty)){
        regValue <- vector("list", length = ncol(best_penalty))
        names(regValue) <- colnames(best_penalty)
        regValue[colnames(best_penalty)] <- best_penalty
      }else{
        regValue <- best_penalty
      }
      reg_CtModel_AIC <- regCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma,
                                    regOn = regOn, regIndicators = regIndicators,
                                    regValues = regValue, link = link, dt = dt,
                                    scaleFactors = scaleFactors)
      reg_CtModel_AIC <- mxOption(reg_CtModel_AIC, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_AIC <- mxOption(reg_CtModel_AIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_AIC <- mxRun(reg_CtModel_AIC, silent = T)
      out <- list("best penalty" = minimum_AIC, "bestmodel" = fit_reg_CtModel_AIC, "fit measures" = results, "call" = call)
    }

    if(criterion == "BIC"){
      if(nrow(minimum_BIC)>1){
        warning(paste("Multiple minima for BIC found. Will use the following penalty values:", paste(as.matrix(minimum_BIC)[nrow(minimum_BIC),],collapse = ","), sep = ":"))
        best_penalty = t(as.matrix(as.matrix(minimum_BIC)[nrow(minimum_BIC),]))
      }else{
        best_penalty = minimum_BIC
      }

      if(is.matrix(best_penalty)){
        regValue <- vector("list", length = ncol(best_penalty))
        names(regValue) <- colnames(best_penalty)
        regValue[colnames(best_penalty)] <- best_penalty
      }else{
        regValue <- best_penalty
      }
      reg_CtModel_BIC <- regCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma,
                                    regOn = regOn, regIndicators = regIndicators,
                                    regValues = regValue, link = link, dt = dt,
                                    scaleFactors = scaleFactors)
      reg_CtModel_BIC <- mxOption(reg_CtModel_BIC, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_BIC <- mxOption(reg_CtModel_BIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_BIC <- mxRun(reg_CtModel_BIC, silent = T)
      out <- list("best penalty" = minimum_BIC, "bestmodel" = fit_reg_CtModel_BIC, "fit measures" = results, "call" = call)
    }

    if(criterion == "CV.m2LL"){
      if(nrow(minimum_CV.m2LL)>1){
        warning(paste("Multiple minima for CV.m2LL found. Will use the following penalty values:", paste(as.matrix(minimum_CV.m2LL)[nrow(minimum_CV.m2LL),],collapse = ","), sep = ":"))
        best_penalty = t(as.matrix(as.matrix(minimum_CV.m2LL)[nrow(minimum_CV.m2LL),]))
      }else{
        best_penalty = minimum_CV.m2LL
      }

      if(is.matrix(best_penalty)){
        regValue <- vector("list", length = ncol(best_penalty))
        names(regValue) <- colnames(best_penalty)
        regValue[colnames(best_penalty)] <- best_penalty
      }else{
        regValue <- best_penalty
      }
      reg_CtModel_CVm2LL <- regCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma,
                                       regOn = regOn, regIndicators = regIndicators,
                                       regValues = regValue, link = link, dt = dt,
                                       scaleFactors = scaleFactors)
      reg_CtModel_CVm2LL <- mxOption(reg_CtModel_CVm2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_CVm2LL <- mxOption(reg_CtModel_CVm2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_CVm2LL <- mxRun(reg_CtModel_CVm2LL, silent = T)
      out <- list("best penalty" =minimum_CV.m2LL, "bestmodel" = fit_reg_CtModel_CVm2LL, "fit measures" = results, "call" = call)
    }

    stopCluster(cl)

    class(out) <- "OptimRegCtModelObject"

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
    # separate variables from dT and intervalID
    colsWithData <- grep("Y", colnames(full_raw_data))

    RegValuesGrid <- createRegValueGrid(regValues = regValues, regOn = regOn)

    Res <- vector("list", length = k+4)
    names(Res) <- c("penalty", "mean CV/Boot_m2LL", paste("fold", 1:k),
                    "negative variances","convergence")
    Res[["penalty"]] <- RegValuesGrid
    Res[["mean CV/Boot_m2LL"]] <- matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1)
    Res[["negative variances"]] <- matrix(0, nrow = nrow(RegValuesGrid), ncol = 1)
    Res[["convergence"]] <- matrix(0, nrow = nrow(RegValuesGrid), ncol = 1)
    for(i in 1:k){
      Res[[paste("fold", i)]] <- matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1)
    }

    fold <- 1
    for(Fold in Folds){
      Test_Sample <- full_raw_data[Fold,]
      Train_Sample <- full_raw_data[-Fold,]

      if(scaleCV){
        Test_Sample[,colsWithData] <- scale(Test_Sample[,colsWithData])
        Train_Sample[,colsWithData] <- scale(Train_Sample[,colsWithData])
      }
      #trainModel <- mxModelObject
      #trainModel$data <- mxData(observed = Train_Sample, type = "raw")
      Test_Sample <- mxData(observed = Test_Sample, type = "raw")
      trainModel <- ctsemModelObject
      trainModel$mxobject$data <- mxData(observed = Train_Sample, type = "raw")

      fit_trainModel <- optimRegCtModel(ctsemModelObject = trainModel, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators,
                                        link = link, dt = dt,
                                        regValues = regValues,
                                        criterion = "CV.m2LL", autoCV = FALSE, Boot = FALSE, manualCV = Test_Sample, k = k, zeroThresh = zeroThresh, scaleCV = scaleCV, scaleFactors = scaleFactors, cores = cores)

      Res[[paste("fold", fold)]] <- fit_trainModel$`fit measures`[['CV.m2LL']]
      Res[["negative variances"]] <- Res[["negative variances"]] + fit_trainModel$`fit measures`[['negative variances']]
      Res[["convergence problems"]] <- Res[["convergence problems"]] + fit_trainModel$`fit measures`[['convergence']]

      cat(paste("\n Completed CV for fold", fold, "of", k, "\n"))

      fold <- fold + 1
    }

    # mean the m2LLs:
    m2LLs <- Res[[paste("fold", 1)]]
    for(i in 2:k){
      m2LLs <- cbind(m2LLs, Res[[paste("fold", i)]])
    }

    Res[["mean CV/Boot_m2LL"]] <- matrix(apply(m2LLs, 1, mean), ncol = 1)

    # only use runs without problems:
    convergedSubset <- (Res[["convergence"]] == 0) & (Res[["negative variances"]]==0)

    # find best penalty value:
    # find best penalty value:

    rowIndicators <- which(Res[["mean CV/Boot_m2LL"]][convergedSubset] == min(Res[["mean CV/Boot_m2LL"]][convergedSubset]))
    best_penalty <- matrix(Res[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
    colnames(best_penalty) = regOn

    # fit best penalty model with full data set:
    if(scaleCV){
      scale_full_raw_data <- full_raw_data
      scale_full_raw_data[,colsWithData] <- scale(full_raw_data[,colsWithData])
      ctsemModelObject$mxobj$data <- mxData(observed = scale_full_raw_data, type = "raw")
    }
    if(nrow(best_penalty)>1){
      warning(paste("Multiple minima for CV.m2LL found. Will use the following penalty values:", paste(as.matrix(best_penalty)[nrow(best_penalty),],collapse = ","), sep = ":"))
      best_penalty = t(as.matrix(as.matrix(best_penalty)[nrow(best_penalty),]))
    }else{
      best_penalty = best_penalty
    }

    if(is.matrix(best_penalty)){
      regValue <- vector("list", length = ncol(best_penalty))
      names(regValue) <- colnames(best_penalty)
      regValue[colnames(best_penalty)] <- best_penalty
    }else{
      regValue <- best_penalty
    }
    finalModel <- regCtModel(ctsemModelObject = ctsemModelObject, alpha = alpha, gamma = gamma,
                             regOn = regOn, regIndicators = regIndicators,
                             regValue = regValue, link = link, dt = dt, scaleFactors = scaleFactors)
    finalModel <- mxOption(finalModel, "Calculate Hessian", "No") # might cause errors; check
    finalModel <- mxOption(finalModel, "Standard Errors", "No") # might cause errors; check

    ffinalModel <- mxRun(finalModel, silent = T)

    ret <- list("CV results" = Res, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k, "call" = call)
    class(ret) <- "CvOptimRegCtModelObject"

    return(ret)

  }



}
