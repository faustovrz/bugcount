#'  Plot wide formated dataframe correlations per clone
#' in log scale
#' @param wf_wide wide formatted data by clone
#' @return \code{pairs} plot of correlation
#'
#' @examples
#' # plot_wide.cor(efficacy.list(ef)$mean)
#' @export
plot_wide_cor <- function(wf_wide, main = "Correlation", ...) {
  pairs(log10(wf_wide[,-1]),
        lower.panel = panel_points,
        diag.panel = panel_hist,
        upper.panel = panel_cor, 
        main = main)
}
