#' SingleCoreOptimRegCtModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' SingleCoreOptimRegCtModel creates a range of regularized models from a ctsem. It automatically tests different penalty values and returns the best model. It uses a single core.
#'
#' @param ctsemModelObject an already run ctsem object
#' @param alpha alpha controls the type of penalty. For lasso regularization, set alpha = 1, for ridge alpha = 0. Values between 0 and 1 implement elastic net regularization
#' @param gamma gamma sets the power in the denominator of parameter specific weights when using adaptive lasso regularization. Make sure to set alpha to 1 when using a gamma other than 0.
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
#'                                 link = link, dt = dt, autoCV = F, criterion = "BIC")
#'
#' # extract discrete time parameters for deltaT = 1
#' expm(myRegCtModel$bestmodel$Submodel$DRIFT$values)
#'
#' # optimize regularized ctsem with cross-validation
#'
#' CV_myRegCtModel <- optimRegCtModel(ctsemModelObject = fit_myModel, regType = "lasso",
#'                                    regOn = regOn, regIndicators = regIndicators,
#'                                    link = link, dt = dt, autoCV = T, k = 5, regValue_start = 0, regValue_end = .5, regValue_stepsize = .01)
#'
#' # extract discrete time parameters for deltaT = 1
#' expm(CV_myRegCtModel$`final Model`$Submodel$DRIFT$values)
#'
#' @author Jannik Orzek
#' @import OpenMx ctsem
#' @export

SingleCoreOptimRegCtModel <- function(ctsemModelObject, alpha = 1, gamma = 0, regOn, regIndicators,
                                      link = list("exp"), dt,
                                      regValues,
                                      criterion = "BIC", autoCV = FALSE, Boot = FALSE, manualCV = NULL, k = 5, zeroThresh = .001, scaleCV = TRUE, cores = 1){

  # save call:
  call <- mget(names(formals()),sys.frame(sys.nframe()))

  ## Checks
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
  checkRegularizedMatrixExistance(regOn = regOn, mxModelObject = mxModelObject, regIndicators = regIndicators)

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

    counter <- 1

    # iterate over regValues:
    # create a grid with all possible combinations of regValues:

    RegValuesGrid <- createRegValueGrid(regValues = regValues, regOn = regOn)

    # create progress bar:
    pb <- txtProgressBar(min = 1, max = nrow(RegValuesGrid), style = 3)

    # create list to store results

    results <- list("penalty"= RegValuesGrid, "estimated Parameters"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "m2LL"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),"AIC"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),
                    "BIC"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "CV.m2LL"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "negative variances"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),"convergence"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1))

    for(row in 1:nrow(RegValuesGrid)){
      if(is.matrix(RegValuesGrid)){
        regValue <- vector("list", length = ncol(RegValuesGrid))
        names(regValue) <- colnames(RegValuesGrid)
        regValue[colnames(RegValuesGrid)] <- RegValuesGrid[row,colnames(RegValuesGrid)]
      }else{
        regValue <- RegValuesGrid[row]
      }

      reg_ctModel <- regCtModel(ctsemModelObject = ctsemModelObject, link = link, dt = dt,
                                alpha = alpha, gamma = gamma, regOn = regOn,
                                regIndicators = regIndicators, regValue = regValue)

      reg_ctModel <- mxOption(reg_ctModel, "Calculate Hessian", "No") # might cause errors; check
      reg_ctModel <- mxOption(reg_ctModel, "Standard Errors", "No") # might cause errors; check
      fit_reg_ctModel <- mxRun(reg_ctModel, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

      results[["convergence"]][counter] <- fit_reg_ctModel$output$status$code# check convergence

      if("S" %in% names(fit_reg_ctModel$Submodel)){
        variances = diag(nrow(fit_reg_ctModel$Submodel$S$values))==1

        if(any(fit_reg_ctModel$Submodel$S$values[variances] <0)){
          results[["negative variances"]][counter] <- 1 # check negative variances
        }else(
          results[["negative variances"]][counter] <- 0
        )}

      ### compute AIC and BIC:
      FitM <- getCtFitMeasures(regCtModel = fit_reg_ctModel, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, cvSample = manualCV, zeroThresh = zeroThresh)

      results[["estimated Parameters"]][counter] <- FitM$estimated_params # estimated parameters
      results[["m2LL"]][counter] <- FitM$m2LL # -2LogL
      results[["AIC"]][counter] <- FitM$AIC # AIC
      results[["BIC"]][counter] <- FitM$BIC # BIC
      results[["CV.m2LL"]][counter] <- FitM$CV.m2LL # CV.m2LL

      setTxtProgressBar(pb, counter)
      counter <- counter+1
    }

    # Find Minima / best penalty value

    convergedSubset <- (results[["convergence"]] == 0) & (results[["negative variances"]]==0)

    minimum_m2LL <- findMinumumCriterion(results = results, criterion = "m2LL", convergedSubset = convergedSubset, regOn = regOn)

    minimum_AIC <- findMinumumCriterion(results = results, criterion = "AIC", convergedSubset = convergedSubset, regOn = regOn)

    minimum_BIC <- findMinumumCriterion(results = results, criterion = "BIC", convergedSubset = convergedSubset, regOn = regOn)

    minimum_CV.m2LL <- findMinumumCriterion(results = results, criterion = "CV.m2LL", convergedSubset = convergedSubset, regOn = regOn)


    # getting parameters:
    if(criterion == "m2LL"){
      out <- computeFinalctParameters(minimum_criterion = minimum_m2LL, ctsemModelObject = ctsemModelObject,
                                    alpha = alpha, gamma = gamma, regOn = regOn,
                                    regIndicators = regIndicators, link = link, dt = dt, results = results)
    }

    if(criterion == "AIC"){
      out <- computeFinalctParameters(minimum_criterion = minimum_AIC, ctsemModelObject = ctsemModelObject,
                                      alpha = alpha, gamma = gamma, regOn = regOn,
                                      regIndicators = regIndicators, link = link, dt = dt, results = results)
    }

    if(criterion == "BIC"){
      out <- computeFinalctParameters(minimum_criterion = minimum_BIC, ctsemModelObject = ctsemModelObject,
                                      alpha = alpha, gamma = gamma, regOn = regOn,
                                      regIndicators = regIndicators, link = link, dt = dt, results = results)
    }

    if(criterion == "CV.m2LL"){
      out <- computeFinalctParameters(minimum_criterion = minimum_CV.m2LL, ctsemModelObject = ctsemModelObject,
                                      alpha = alpha, gamma = gamma, regOn = regOn,
                                      regIndicators = regIndicators, link = link, dt = dt, results = results)
    }
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

    # create a grid with all possible combinations of regValues:

    RegValuesGrid <- createRegValueGrid(regValues = regValues, regOn = regOn)

    # create list to store results

    results <- createResultsList(RegValuesGrid = RegValuesGrid, autoCV = autoCV, k = k)

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
                                        criterion = "CV.m2LL", autoCV = FALSE, Boot = FALSE, manualCV = Test_Sample, k = 5, zeroThresh = zeroThresh, scaleCV = scaleCV)

      results[[paste("fold", fold)]] <- fit_trainModel$`fit measures`[['CV.m2LL']]
      results[["negative variances"]] <- results[["negative variances"]] + fit_trainModel$`fit measures`[['negative variances']]
      results[["convergence problems"]] <- results[["convergence problems"]] + fit_trainModel$`fit measures`[['convergence']]

      cat(paste("\n Completed CV for fold", fold, "of", k, "\n"))

      fold <- fold + 1
    }

    # mean of the m2LLs:
    m2LLs <- results[[paste("fold", 1)]]
    for(i in 2:k){
      m2LLs <- cbind(m2LLs, results[[paste("fold", i)]])
    }

    results[["mean CV/Boot_m2LL"]] <- matrix(apply(m2LLs, 1, mean), ncol = 1)

    # only use runs without problems:
    convergeSubset <- (results[["convergence"]] == 0) & (results[["negative variances"]]==0)

    # find best penalty value:

    rowIndicators <- which(results[["mean CV/Boot_m2LL"]][convergeSubset] == min(results[["mean CV/Boot_m2LL"]][convergeSubset]))
    best_penalty <- matrix(results[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
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
                             regValue = regValue, link = link, dt = dt)
    finalModel <- mxOption(finalModel, "Calculate Hessian", "No") # might cause errors; check
    finalModel <- mxOption(finalModel, "Standard Errors", "No") # might cause errors; check

    ffinalModel <- mxRun(finalModel, silent = T)

    ret <- list("CV results" = results, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k, "call" = call)
    class(ret) <- "CvOptimRegCtModelObject"

    return(ret)

  }



}
