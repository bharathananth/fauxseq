#' Generate rhythmic count RNA-seq data
#'
#' Generate rhythmic count RNA-seq data with same properties of a reference
#' dataset.
#'
#' @param t A numeric vector.
#' @param reps NULL or numeric scalar or numeric vector with the same length as
#'   `t`
#' @param period A numeric scalar
#' @param n_genes An integer.
#' @param rhy_frac A fraction between 0 and 1.
#' @param min_A_effect A positive numeric scalar.
#' @param A_spread A positive numeric scalar
#' @param emp_dist A data.frame with two columns containing the empirical mu
#'   (means) and empirical size (reciprocal of dispersion)
#' @param depth A positive integer.
#' @param lib_size_var A numeric vector of two positive elements.
#' @param seed A positive integer.
#'
#' @returns A list containing a integer matrix, a data.frame of feature
#'   properties and a numeric vector of sampling times.
#' @export
#'
generate_rhy_rnaseq <- function(t = seq(0, 20, 4),
                                reps = NULL,
                                period = 24,
                                n_genes = 10000,
                                rhy_frac = 0.1,
                                min_A_effect = 0.0,
                                A_spread = 1,
                                emp_dist = NULL,
                                depth = 1e7,
                                lib_size_var = c(0.8, 1.2),
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(reps)) {
    reps <- rep_len(1, length(t))
  } else {
    if (length(reps) != length(t)) {
      stop("The length of timepoints and replicates do not match.")
    }
  }
  t <- rep(t, times = reps)

  N <- length(t)
  G_rhy <- (stats::runif(n_genes) <= rhy_frac)
  A <- matrix((1 + min_A_effect + stats::rexp(n_genes, A_spread)),
                nrow = n_genes)
  A[!G_rhy] <- 1.0
  phi <- matrix(stats::runif(n_genes) * 2 * pi * G_rhy, nrow = n_genes)
  params <- data.frame(rhy_genes = which(G_rhy), A = A[G_rhy], phi = phi[G_rhy])

  t_pattern <- cos(phi) %*% cos(2 * pi * t / period) +
                  sin(phi) %*% sin(2 * pi * t / period)

  lib_size_fct <- stats::runif(N, min = lib_size_var[1], max = lib_size_var[2])
  lib_size <- lib_size_fct * depth

  draw <- sample(nrow(emp_dist), n_genes, replace = TRUE)
  lambda <- matrix(emp_dist$mu[draw], nrow = n_genes)

  size <- matrix(emp_dist$size[draw], nrow = n_genes)
  size <- size[, rep(1, length(t))]

  lambda <- lambda[, rep(1, N)] * A[, rep(1, N)]^t_pattern

  mu <- sweep(lambda, 2, colSums(lambda) / lib_size, FUN = "/")

  counts <- matrix(mapply(\(x, y) stats::rnbinom(1, size = x, mu = y), size, mu),
                   nrow = n_genes, byrow = FALSE)

  return(list(counts = counts, params = params, t = t))
}
