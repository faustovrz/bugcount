#'  Plot ANOVA per cross
#' @param postcoc \code{multcomp} common letter display \code{cld} object from
#'        posthoc comparison in \code{plot_fit.nb.glm()}
#' @param exp_allowed list of vectors containing indices of experiments
#'        to merge in order to make ANOVA analysis
#' @return Nothing
#'
#' @examples
#' # plot_cross.anova(posthoc, exp_allowed)
#' @export

plot_cross_anova <- function(posthoc, exp_allowed){
  
  ### Order the levels for printing
  
  clone_order <- posthoc$clone
  posthoc$clone <- factor(posthoc$clone,
                          levels = clone_order)
  ###  Remove spaces in .group  
  
  posthoc$.group = gsub(" ", "", posthoc$.group)
  posthoc_str <- paste("posthoc groups :", length(unique(posthoc$.group)))
  
  print(ggplot2::ggplot(posthoc,
        ggplot2::aes(x     = clone,
                     y     = lsmean,
                     color = .group,
                     label = .group),
                     log10 = "y") +
        ggplot2::geom_point(shape  = 15, size = 1) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
            width =  0.2,
            size  =  0.7) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       legend.position = "none") +
          # TODO: Calculate interval (min, max)
          # scale_y_continuous(limits = c(-1000,6000)) +  
        ggplot2::ggtitle(
          paste(c(exp_allowed, posthoc_str, 
                  "ANOVA, Mean + 95% Confidence Interval"),
                  collapse = "\n")))
}

