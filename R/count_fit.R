#'  Fit different distributions to count data
#' @param wf A whitefly count dataframe.
#' @return List of distribution parameters for diferent fits.
#'         and histogram data for plotting.
#'
#' @examples
#' # fit.list <- count.fit(plant.count)
#' @export
count_fit <- function(wf){
  poisson <- MASS::fitdistr(wf$nymphs, "Poisson")
  nbinom  <- MASS::fitdistr(wf$nymphs, "Negative Binomial")
  nymph_hist <-  hist(wf$nymphs, breaks = 50, plot = FALSE)
  mids <- nymph_hist$mids
  
  list( hist = nymph_hist,
        nb_d = dnbinom(mids,
                       mu = nbinom$estimate[[2]],
                       size = nbinom$estimate[[1]]),
        
        lnorm_d = dlnorm(mids,
                         meanlog = mean(log(wf$nymphs + 1)),
                         sdlog = sd(log((wf$nymphs + 1)))),
        
        norm_d = dnorm(mids,
                        mean = mean(wf$nymphs),
                        sd = sd(wf$nymphs)),
        
        pois_d = dpois(mids,
                        lambda = poisson$estimate[[1]])
  )
}
