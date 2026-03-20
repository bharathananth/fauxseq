#' Estimate the mean-dispersion of count dataset
#'
#' @param counts A integer matrix.
#' @param parallel A boolean, whether samples must be parallelly processed,
#'   default TRUE
#' @param ncores An integer, the number of cores to use.
#'
#' @returns A data.frame of mean and size of negative binomial fits.
#' @export
#'
estimate_disp_dist <- function(counts,
                               parallel = TRUE,
                               ncores = NULL) {
  counts <- downsample(counts, parallel = parallel, ncores = ncores)
  ests <- apply(counts, 1, nb_fit)
  final_ests <- t(ests) |>
    stats::na.omit() |>
    data.frame()
  return(final_ests)
}

#' Fit negative binomial distribution using MLE
#'
#' @param x A numeric vector of counts.
#'
#' @returns A numieric vector containing the mean and size of fit.
#' @export
#'
nb_fit <- function(x) {
  if (stats::var(x) < mean(x)) {
    return(c(size = 1e6, mu = mean(x)))
  }
  tryCatch({
    m <- mean(x)
    v <- stats::var(x)
    suppressWarnings(
      fit <- MASS::fitdistr(
        x,
        "negative binomial",
        method = "Brent",
        start = list(size = m / (v - m)),
        mu = m,
        lower = 0,
        upper = 1e6
      )
    )
    return(c(size = fit$estimate, mu = m))
  }, error = function(e) {
    return(c(size = NA, mu = NA))
  })
}
