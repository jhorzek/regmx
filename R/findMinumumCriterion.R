#' findMinumumCriterion
#'
#'
#' returns the best tuning parameter value for a specified criterion
#' @author Jannik Orzek
#' @examples
#'
#' @export
findMinumumCriterion <- function(results, criterion, convergedSubset, regOn){
  rowIndicators <- which(results[[criterion]][convergedSubset] == min(results[[criterion]][convergedSubset]))
  minimum_criterion <- matrix(results[["penalty"]][rowIndicators,], ncol = length(regOn), byrow = T)
  colnames(minimum_criterion) = regOn
  return(minimum_criterion)
}
