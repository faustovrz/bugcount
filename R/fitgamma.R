#'  Equivalent to `MASS::fitdistr(x, densfun = "gamma")`, where x are first
#' rescaled to the appropriate scale for a gamma distribution.
#' Useful for fitting the gamma distribution to 
#' data which, when multiplied by a constant, follows this distribution
#' @param x A vector.
#' @return MASS::fitdistr with scaling multiplier
#'
#' @examples
#' fitgamma(1:1000)
#' @export

fitgamma <- function(x) {
  if (!requireNamespace("MASS")) stop("Requires MASS package.")
  fit <- glm(formula = x ~ 1, family = Gamma)
  out <- MASS::fitdistr(x * coef(fit), "gamma")
  out$scaling_multiplier <- unname(coef(fit))
  out
}
