#'  Plot nb GLM per cross
#' @param posthoc \code{multcomp} common letter display \code{cld} object from
#'        posthoc comparison in \code{plot_fit.nb.glm()}
#' @param exp_allowed list of vectors containing indices of experiments
#'        to merge in order to make nb GLM analysis
#' @return Nothing
#'
#' @examples
#' # nb.fit(posthoc, exp_allowed)
#' @export
#

plot_cross_nb_fit <- function(posthoc, exp_allowed, breaks = NULL){
  
  ### Order the levels for printing
  posthoc$clone <- factor(posthoc$clone, levels = posthoc$clone)

  ###  Remove spaces in .group  
  
  posthoc$.group = gsub(" ", "", posthoc$.group)
  
  posthoc_str <- paste("posthoc groups :", length(unique(posthoc$.group)))
  
  ggplot2::ggplot(posthoc,
            ggplot2::aes(x     = clone,
                         y     = response,
                         color = .group,
                         label = .group)) +
          ggplot2::geom_point(shape  = 15,
                              size = ggplot2::rel(0.5)) +
          ggplot2::geom_errorbar(
          ggplot2::aes(ymin  =  asymp.LCL, 
                        ymax  =  asymp.UCL),
                        width = ggplot2::rel(0.5)) +
          ggplot2::theme(legend.position = "none",
            axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major =  ggplot2::element_line(colour = "grey"),
            panel.background = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()) +
          ggplot2::scale_y_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
            limits = c(0.01,100000)) +
          ggplot2::scale_x_discrete(breaks = breaks) +
          ggplot2::ggtitle(paste(c(exp_allowed, posthoc_str),
                                 collapse = ", "))
}
