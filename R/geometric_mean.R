#'  Geometric mean of a vector ignoring NAs.
#' 
#' @param x A vector.
#' @return Geometric mean of  \code{x} ignoring NAs.
#'
#' @examples
#' geometric.mean(1:1000)
#' @export

geometric_mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}


