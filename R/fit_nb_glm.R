#'  Fit MASS Negative binomial General Linear Model for whitefly counts
#' @param formula formula object specifiying the GLM
#' @param wf A whitefly count dataframe. 
#' @return List with model fit, lsmeans and posthoc groups
#'
#' @examples
#' # fit.nb.glm(nymphs ~ year + leaf, wf.per.leaf)
#' @export

fit_nb_glm <-function(formula, wf){
  # TODO: validate formula for one predictor only
  response <- as.character(terms(formula)[[2]])
  predictor <- as.character(terms(formula)[[3]])
  fit <- MASS::glm.nb(formula, data = wf)

  lsm_fit <- lsmeans::lsmeans(
    fit, predictor,
    adjust = "tuckey",
    data = wf)
  ls_posthoc <- lsmeans::cld(lsm_fit, type = "response")
  #ls.posthoc$experiment <- gsub(" ","",ls.posthoc$.group)
  list( formula = formula,
        response = response,
        predictor = predictor,
        glm = fit, lsm = lsm_fit,
        posthoc = ls_posthoc,
        order = ls_posthoc[,predictor],
        group = as.factor(ls_posthoc$.group)
        )
}
