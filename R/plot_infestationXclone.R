#'  Plot infestation X clone interaction as boxplot
#' @param wf A whitefly count dataframe.
#' @return Nothing
#'
#' @examples
#' # plot_infestationXclone(wf_per_leaf)
#' @export

plot_infestationXclone <- function(wf) {
  # TODO:  validate that infestation is a column of wf
  wf_aggregate(wf, by_stat = "gmean")
  clone_order <- wf_aggregate(wf, by_stat = "gmean")$clone
  n_clones <- length(clone_order)
  
  wf$clone <- factor(wf$clone, levels = clone_order)
  
  par(mar = c(6,4,4,2)) 
  boxplot(nymphs+1 ~  infestation + clone, data = wf, frame = FALSE,
          border = "white", log = "y", xaxt = "n", yaxt = "n")
  
  rect(seq(0.5, 2 * n_clones, by = 4),
       rep(1,n_clones),
       seq(2.5, 2 * n_clones + 0.5, by = 4),
       rep(100000,n_clones),
       col = "grey90",lty = 0)
  
  boxplot(nymphs+1 ~  infestation + clone, las = 2, log = "y",
          main = "Interaction: Infestation Regime x Clone",
          ylab = "Nymphs per Unit",xaxt = "n", xlab = "",
          col = c("white", "grey"),
          data = wf, add = TRUE)
  
  axis(1, at = 2 * (1:n_clones) - 0.5, labels = clone_order, las = 2)
  
  mtext("Clone (sorted by geometric mean of nymphs)", side = 1, line = 5)
  
  legend(2 * n_clones - 3, 6, legend = c("High", "Low"),
         title = "Infestation",
         bg = "white",
         fill = c("grey", "white"))
}
