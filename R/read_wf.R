#'  Read whitefly count table
#' 
#' @param file A filename for the whitefly counts per picture.
#' @return Dataframe containing whitefly counts
#' @details 
#' Whitefly counts per picture. Each picture must have a unique combination of
#' experiment + replicate + plant + pic
#' experiment description consists of:
#' exp.year + exp.cross + exp.propagation + exp.substrate
#' @examples
#' # my.wf <- read.wf("whitefly_counts.tab")
#' @export

read_wf <- function(file,...){
  wf <- read.table(file, sep="\t",header=TRUE)
  wf <- wf[!(is.na(wf$nymphs)),]
  # Add experiment as factor = exp.propagation + exp.substrate + exp.year
  wf$experiment <- apply( wf[, exp_descriptor(wf) ] ,
                          1 , paste , collapse = "_" )
  wf
  # TODO: order columns with expriment at the beginning. new function?
}
