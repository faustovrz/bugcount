#'  Emulate ggplot color palette.
#' @param n number of colors to generate
#' @return vector of hexadecimal color strings
#'
#' # @examples
#' gg_color_hue(length(levels(wf$clone)))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# Correlation Analysis, reproducibility ########################################
