#'  Plot Probability Density over histogram to appreciate fit
#' @param x A \code{count.fit()} list
#' @return plot_
#'
#' @examples
#' # plot_count.fit(count.fit(wf.per.leaf))
#' @export

plot_count_fit <- function(x,...){
  plot(x$hist, xlab = "Nypmh Count",...)
  lines(x$hist$mids,max(x$hist$counts) * x$nb_d / max(x$nb_d),
         col = gg_color_hue(4)[1], lwd = 2, lty = 1)
  lines(x$hist$mids,max(x$hist$counts) * x$lnorm_d / max(x$lnorm_d),
        col = gg_color_hue(4)[2], lwd = 2, lty = 2)
  lines(x$hist$mids,max(x$hist$counts) * x$pois_d / max(x$pois_d),
        col = gg_color_hue(4)[3], lwd = 2, lty = 3)
  lines(x$hist$mids,max(x$hist$counts) * x$norm_d / max(x$norm_d),
        col = gg_color_hue(4)[4], lwd = 2, lty = 4)
  
  # logspline aproximation, not good.
  # lines(wf.hist$mids,max(wf.hist$counts)*logspline.d/max(logspline.d))
  
  legend("topright",legend=factor(c("Negative Binomial", "Lognormal",
                                    "Poisson","Normal"),
                                  levels = c("Negative Binomial","Lognormal",
                                             "Poisson","Normal")),
         title = "Model Aproximation",
         col= gg_color_hue(4), lwd = 2, lty = 1:4, bty = "n")
}
