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
#' @import OpenMx ctsem doParallel
#' @export

MultiCoreOptimRegCtModel <- function(ctsemModelObject, regType = "lasso", regOn, regIndicators,
                                     link = list("exp"), dt,
                                     regValue_start = 0, regValue_end = 1, regValue_stepsize = .01,
                                     criterion = "BIC", autoCV = FALSE, Boot = FALSE, manualCV = NULL, k = 5, zeroThresh = .001, scaleCV = TRUE, cores){

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

    # iterate over regValues:
    ParRes <- foreach(regValue = regValues, .combine = "cbind",.packages = c("OpenMx","regmx"),
                      .inorder = FALSE,
                      .errorhandling = "stop",
                      .verbose = TRUE) %dopar% {
                        results <- matrix(NA, nrow = 8, ncol = 1)
                        rownames(results) <- c("penalty", "estimated Parameters", "m2LL","AIC",
                                               "BIC", "CV.m2LL", "negative variances","convergence")
                        results["penalty",1] <- regValue
                        reg_ctModel <- regCtModel(ctsemModelObject = ctsemModelObject, link = link, dt = dt,
                                                  regType = regType, regOn = regOn,
                                                  regIndicators = regIndicators, regValue = regValue)

                        reg_ctModel <- mxOption(reg_ctModel, "Calculate Hessian", "No") # might cause errors; check
                        reg_ctModel <- mxOption(reg_ctModel, "Standard Errors", "No") # might cause errors; check
                        fit_reg_ctModel <- mxRun(reg_ctModel, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

                        results["convergence",1] <- fit_reg_ctModel$output$status$code# check convergence

                        if("S" %in% names(fit_reg_ctModel$Submodel)){
                          variances = diag(nrow(fit_reg_ctModel$Submodel$S$values))==1

                          if(any(fit_reg_ctModel$Submodel$S$values[variances] <0)){
                            results["negative variances",1] <- 1 # check negative variances
                          }else(
                            results["negative variances",1] <- 0
                          )}

                        ### compute AIC and BIC:
                        FitM <- getCtFitMeasures(regCtModel = fit_reg_ctModel, regType = regType, regOn = regOn, regIndicators = regIndicators, cvSample = manualCV, zeroThresh = zeroThresh)

                        results["estimated Parameters",1] <- FitM$estimated_params # estimated parameters
                        results["m2LL",1] <- FitM$m2LL # -2LogL
                        results["AIC",1] <- FitM$AIC # AIC
                        results["BIC",1] <- FitM$BIC # BIC
                        results["CV.m2LL",1] <- FitM$CV.m2LL # CV.m2LL

                        # foreach expects a vector:
                        results <- as.vector(results)
                      }

    rownames(ParRes) <- c("penalty", "estimated Parameters", "m2LL","AIC",
                          "BIC", "CV.m2LL", "negative variances","convergence")

    results <- ParRes

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
      reg_CtModel_m2LL <- regCtModel(ctsemModelObject = ctsemModelObject, regType = regType,
                                     regOn = regOn, regIndicators = regIndicators,
                                     regValue = best_penalty, link = link, dt = dt)
      reg_CtModel_m2LL <- mxOption(reg_CtModel_m2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_m2LL <- mxOption(reg_CtModel_m2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_m2LL <- mxRun(reg_CtModel_m2LL, silent = T)
      out <- list("best penalty" = minimum_m2LL, "bestmodel" = fit_reg_CtModel_m2LL, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "AIC"){
      best_penalty = minimum_AIC
      reg_CtModel_AIC <- regCtModel(ctsemModelObject = ctsemModelObject, regType = regType,
                                    regOn = regOn, regIndicators = regIndicators,
                                    regValue = best_penalty, link = link, dt = dt)
      reg_CtModel_AIC <- mxOption(reg_CtModel_AIC, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_AIC <- mxOption(reg_CtModel_AIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_AIC <- mxRun(reg_CtModel_AIC, silent = T)
      out <- list("best penalty" = minimum_AIC, "bestmodel" = fit_reg_CtModel_AIC, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "BIC"){
      best_penalty = minimum_BIC
      reg_CtModel_BIC <- regCtModel(ctsemModelObject = ctsemModelObject, regType = regType,
                                    regOn = regOn, regIndicators = regIndicators,
                                    regValue = best_penalty, link = link, dt = dt)
      reg_CtModel_BIC <- mxOption(reg_CtModel_BIC, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_BIC <- mxOption(reg_CtModel_BIC, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_BIC <- mxRun(reg_CtModel_BIC, silent = T)
      out <- list("best penalty" = minimum_BIC, "bestmodel" = fit_reg_CtModel_BIC, "fit measures" = t(results), "call" = call)
    }

    if(criterion == "CV.m2LL"){
      best_penalty = minimum_CV.m2LL
      reg_CtModel_CVm2LL <- regCtModel(ctsemModelObject = ctsemModelObject, regType = regType,
                                       regOn = regOn, regIndicators = regIndicators,
                                       regValue = best_penalty, link = link, dt = dt)
      reg_CtModel_CVm2LL <- mxOption(reg_CtModel_CVm2LL, "Calculate Hessian", "No") # might cause errors; check
      reg_CtModel_CVm2LL <- mxOption(reg_CtModel_CVm2LL, "Standard Errors", "No") # might cause errors; check

      fit_reg_CtModel_CVm2LL <- mxRun(reg_CtModel_CVm2LL, silent = T)
      out <- list("best penalty" =minimum_CV.m2LL, "bestmodel" = fit_reg_CtModel_CVm2LL, "fit measures" = t(results), "call" = call)
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
        Test_Sample[,colsWithData] <- scale(Test_Sample[,colsWithData])
        Train_Sample[,colsWithData] <- scale(Train_Sample[,colsWithData])
      }
      #trainModel <- mxModelObject
      #trainModel$data <- mxData(observed = Train_Sample, type = "raw")
      Test_Sample <- mxData(observed = Test_Sample, type = "raw")
      trainModel <- ctsemModelObject
      trainModel$mxobject$data <- mxData(observed = Train_Sample, type = "raw")

      fit_trainModel <- optimRegCtModel(ctsemModelObject = trainModel, regType = regType, regOn = regOn, regIndicators = regIndicators,
                                        link = link, dt = dt,
                                        regValue_start = regValue_start, regValue_end = regValue_end, regValue_stepsize = regValue_stepsize,
                                        criterion = "CV.m2LL", autoCV = FALSE, Boot = FALSE, manualCV = Test_Sample, k = 5, zeroThresh = zeroThresh, scaleCV = scaleCV, cores = cores)

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
      scale_full_raw_data <- full_raw_data
      scale_full_raw_data[,colsWithData] <- scale(full_raw_data[,colsWithData])
      mxModelObject$data <- mxData(observed = scale_full_raw_data, type = "raw")
    }
    finalModel <- regCtModel(ctsemModelObject = ctsemModelObject, regType = regType,
                             regOn = regOn, regIndicators = regIndicators,
                             regValue = best_penalty, link = link, dt = dt)

    ffinalModel <- mxRun(finalModel, silent = T)

    ret <- list("CV results" = Res, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k, "call" = call)
    class(ret) <- "CvOptimRegCtModelObject"

    return(ret)

  }



}
