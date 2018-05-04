#'  Vertical to horizontal data reformatting and aggregation for whitefly counts
#' just a wrapper for \code{reshape2::dcast} from \code{reshape2}
#' @param formula formula object specifiying the \code{reshape2::dcast} reformatting
#' @param wf A whitefly count dataframe.
#' @param fun function for aggregation deafult \code{geometric.mean}
#' @return A wide formated dataframe
#'
#' @examples
#' # clone.means <- wf.wide(nymphs ~ clone,fun = mean)

wf_wide <- function(formula, wf, fun = geometric_mean){
  reshape2::dcast(wf,formula,
        value.var = "nymphs",
        fun.aggregate = fun)
}

# plot_wide.cor<-function(wf.wide, fun = "Geometric Means"){
#   pairs.panels(log10(wf.wide[,-1]), smooth = FALSE, density=FALSE,
#                ellipses=FALSE, rug = FALSE, pch = 1, lm= TRUE, stars = TRUE,
#                main = paste(cross,fun, "Correlation"))
# }

