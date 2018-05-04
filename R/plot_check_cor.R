#'  Plot wide formated dataframe correlations per clone
#' with point labels, asumes plotted points to be ~ 10
#' more points would be extremely unreadable
#' @param check_wide wide formatted data by clone for checks only
#' @param main plot title 
#' @return \code{pairs} plot of correlation
#'
#' @examples
#' # plot_check.cor(eplot_check.cor)
#' @export

plot_check_cor <- function(check_wide,
                          main = "Correlation",
                          ...){
  pairs(check_wide[,-1],
        lower.panel = function(x,y){
          text(x,y, labels = check_wide$clone, cex = 1)
          abline(stats::lm(y ~ x), col = "red")
        },
        labels = gsub("_","\n", colnames(check_wide[,-1])),
        upper.panel = panel_cor, 
        gap = 0,
        main = main
        )
}
