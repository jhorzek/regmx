#' MultiCoreOptimRegModel
#'
#' Note: regmx is based on the R package \pkg{regsem}. Because of the early status of regmx, it is recommended to use regsem instead!
#' MultiCoreOptimRegModel creates a range of regularized models from an mxModel. It automatically tests different penalty values and returns the best model. It only uses multiple cores.
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
#' reg_model <- optimRegModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0, regOn  = regOn, regIndicators = regIndicators, cores = 2)
#'
#' reg_model$`fit measures`
#'
#' reg_model$`best penalty`
#'
#' # Run the same model with 5-fold cross-validation
#'
#' CV_reg_model <- optimRegModel(mxModelObject = fit_myModel, alpha = 1, gamma = 0 regOn  = regOn, regIndicators = regIndicators,
#'                               autoCV = T, k = 5, cores = 2)
#' CV_reg_model$`CV results`
#'
#' @author Jannik Orzek
#' @import OpenMx doParallel
#' @export

MultiCoreOptimRegModel <- function(mxModelObject, alpha = 1, gamma = 0, regOn, regIndicators,
                                   regValues,
                                   criterion = "BIC", autoCV = FALSE, k = 5, Boot = FALSE, manualCV = NULL, zeroThresh = .001, scaleCV = TRUE, cores = 1){

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
  #if(!class(mxModelObject)[1] == "MxRAMModel"){
  #  stop("Provided mxModelObject is not of type MxRAMModel")
  #}
  ### regOn
  checkRegularizedMatrixExistance(regOn = regOn, mxModelObject = mxModelObject, regIndicators = regIndicators)

  ### matrix dimensions
  for(matrix in regOn){
    if(nrow(mxModelObject[[matrix]]$values) == nrow(regIndicators[[matrix]]) &
       ncol(mxModelObject[[matrix]]$values) == ncol(regIndicators[[matrix]])){}else{
         stop("Dimensions of Matrix ", regOn[[matrix]], " and provided regIndicator with index ", matrix, " do not match.", sep = "")
       }
  }


  if (!autoCV & !Boot){

    # create a grid with all possible combinations of regValues:

    RegValuesGrid <- createRegValueGrid(regValues = regValues, regOn = regOn)

    # create list to store results

    results <- createResultsList(RegValuesGrid = RegValuesGrid, autoCV = autoCV, k = k)

    # iterate over regValues:
    comb <- function(...) {
      mapply('rbind', ..., SIMPLIFY=FALSE)
    }

    ParRes <- foreach(row = 1:nrow(RegValuesGrid), .combine='comb', .multicombine=TRUE,.packages = c("OpenMx","regmx"),
                       .inorder = FALSE,
                       .errorhandling = "stop",
                       .verbose = TRUE) %dopar% {
                        regValue <- RegValuesGrid[row,]
                        reg_Model <- regModel(mxModelObject = mxModelObject,
                                              alpha = alpha, gamma = gamma, regOn = regOn,
                                              regIndicators = regIndicators, regValue = regValue)

                        reg_Model <- mxOption(reg_Model, "Calculate Hessian", "No") # might cause errors; check
                        reg_Model <- mxOption(reg_Model, "Standard Errors", "No") # might cause errors; check
                        fit_reg_Model <- mxRun(reg_Model, silent = T) # run Model; starting values can be very critical as the model tends to get stuck in local minima either close to the model parameters without penalty or all parameters set to 0

                        convergence <- fit_reg_Model$output$status$code# check convergence

                        if("S" %in% names(fit_reg_Model$Submodel)){
                          variances = diag(nrow(fit_reg_Model$Submodel$S$values))==1

                          if(any(fit_reg_Model$Submodel$S$values[variances] <0)){
                            negative.variances <- 1 # check negative variances
                          }else(
                            negative.variances <- 0
                          )}

                        ### compute AIC and BIC:
                        FitM <- getFitMeasures(regModel = fit_reg_Model, alpha = alpha, gamma = gamma, regOn = regOn, regIndicators = regIndicators, cvSample = manualCV, zeroThresh = zeroThresh)

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

    minimum_m2LL <- findMinumumCriterion(results = results, criterion = "m2LL", convergedSubset = convergedSubset, regOn = regOn)

    minimum_AIC <- findMinumumCriterion(results = results, criterion = "AIC", convergedSubset = convergedSubset, regOn = regOn)

    minimum_BIC <- findMinumumCriterion(results = results, criterion = "BIC", convergedSubset = convergedSubset, regOn = regOn)

    minimum_CV.m2LL <- findMinumumCriterion(results = results, criterion = "CV.m2LL", convergedSubset = convergedSubset, regOn = regOn)

    # getting parameters:
    if(criterion == "m2LL"){
      out <- computeFinalParameters(minimum_criterion = minimum_m2LL, mxModelObject = mxModelObject,
                                    alpha = alpha, gamma = gamma, regOn = regOn,
                                    regIndicators = regIndicators, results = results)
    }

    if(criterion == "AIC"){
      out <- computeFinalParameters(minimum_criterion = minimum_AIC, mxModelObject = mxModelObject,
                                    alpha = alpha, gamma = gamma, regOn = regOn,
                                    regIndicators = regIndicators, results = results)
    }

    if(criterion == "BIC"){
      out <- computeFinalParameters(minimum_criterion = minimum_BIC, mxModelObject = mxModelObject,
                                    alpha = alpha, gamma = gamma, regOn = regOn,
                                    regIndicators = regIndicators, results = results)
    }

    if(criterion == "CV.m2LL"){
      out <- computeFinalParameters(minimum_criterion = minimum_CV.m2LL, mxModelObject = mxModelObject,
                                    alpha = alpha, gamma = gamma, regOn = regOn,
                                    regIndicators = regIndicators, results = results)
    }
    class(out) <- "OptimRegModelObject"
    stopCluster(cl)
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

    # create a grid with all possible combinations of regValues:

    RegValuesGrid <- createRegValueGrid(regValues = regValues, regOn = regOn)

    # create list to store results

    results <- createResultsList(RegValuesGrid = RegValuesGrid, autoCV = autoCV, k = k)

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
                                      regValues = regValues,
                                      criterion = "CV.m2LL", autoCV = FALSE,
                                      Boot = FALSE, manualCV = Test_Sample,
                                      k = 5, zeroThresh = zeroThresh, scaleCV = scaleCV, cores = cores)

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
      scale_full_raw_data <- scale(full_raw_data)
      mxModelObject$data <- mxData(observed = scale_full_raw_data, type = "raw")
    }
    finalModel <- regModel(mxModelObject = mxModelObject, alpha = alpha, gamma = gamma,
                           regOn = regOn, regIndicators = regIndicators,
                           regValue = best_penalty)

    ffinalModel <- mxRun(finalModel, silent = T)

    ret <- list("CV results" = results, "final Model" = ffinalModel, "best penalty" = best_penalty, "k" = k, "call" = call)
    class(ret) <- "CvOptimRegModelObject"

    return(ret)


  }



}
