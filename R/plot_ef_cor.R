#'  Plot wide formated dataframe correlations per clone
#' in arc sine scale
#' @param ef.par wide formatted data by clone
#' @param main plot title  
#' @return \code{pairs} plot of correlation
#'
#' @examples
#' # plot_wide.cor(efficacy_list(ef)$mean)
#' @export

plot_ef_cor <- function(ef_par, main = "Correlation", ...){
  pairs(asin(ef_par[,-1]),
        lower.panel = panel_points,
        diag.panel = panel_hist,
        upper.panel = panel_cor,
        main = main)
}

  
