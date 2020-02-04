#' createResultsList
#'
#'
#' creates a list to store the results
#' @author Jannik Orzek
#'
#' @export
createResultsList <- function(RegValuesGrid, autoCV, k){
  if(autoCV){
    results <- vector("list", length = k+4)
    names(results) <- c("penalty", "mean CV/Boot_m2LL", paste("fold", 1:k),
                    "negative variances","convergence")
    results[["penalty"]] <- RegValuesGrid
    results[["mean CV/Boot_m2LL"]] <- matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1)
    results[["negative variances"]] <- matrix(0, nrow = nrow(RegValuesGrid), ncol = 1)
    results[["convergence"]] <- matrix(0, nrow = nrow(RegValuesGrid), ncol = 1)
    for(i in 1:k){
      results[[paste("fold", i)]] <- matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1)
    }
  }else{
  results <- list("penalty"= RegValuesGrid, "estimated Parameters"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),
                  "m2LL"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),"AIC"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),
                  "BIC"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1), "CV.m2LL"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),
                  "negative variances"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1),"convergence"=matrix(NA, nrow = nrow(RegValuesGrid), ncol = 1))
  }
  return(results)
}
