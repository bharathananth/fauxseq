#' Generate synthetic rhythmic RNA-seq count data in one condition.
#'
#' This function generates artificial timeseries transcriptomic data under one
#' condition, with rhythmically expressed genes and empirically estimated
#' variability between replicates/samples.
#'
#' @param t A numeric vector, time points at which samples are generated, by
#'   default seq(0, 20, 4)
#' @param reps NULL or numeric scalar or numeric vector with the same length as
#'   `t`, number of replicates at each time point, by default 1
#' @param period A numeric scalar, period of the rhythmic genes, by default 24
#' @param n_genes An integer, number of genes in the dataset, by default 10000
#' @param rhy_frac A fraction between 0 and 1, fraction of genes that are
#'   rhythmic, by default 0.1
#' @param min_A_effect A positive numeric scalar, minimum amplitude (in log2
#'   fold) of rhythmic genes, by default 0.5
#' @param A_spread A positive numeric scalar, mean of the exponential
#'   distribution of rhythmic gene amplitudes (in log2 fold), by default 1
#' @param emp_dist A data.frame with two columns containing the empirical mu
#'   (means) and empirical size (reciprocal of dispersion), by default NULL
#'   (read from one of five mouse tissue datasets at random)
#' @param depth A positive integer, average sequencing depth of different
#'   samples, by default 4e7
#' @param lib_size_var A numeric vector of two positive elements, window of
#'   variability of sequencing depth of samples about 'depth', by default c(0.8,
#'   1.2)
#' @param seed A positive integer, seed to ensure reproducible datasets, by
#'   default NULL
#'
#' @returns A list containing an integer count matrix, a data.frame of feature
#'   properties and a data.frame of the sample data.
#' @export
#'
generate_rhythmic_rnaseq <- function(t = seq(0, 20, 4),
                                     reps = NULL,
                                     period = 24,
                                     n_genes = 10000,
                                     rhy_frac = 0.1,
                                     min_A_effect = 0.5,
                                     A_spread = 1,
                                     emp_dist = NULL,
                                     depth = 4e7,
                                     lib_size_var = c(0.8, 1.2),
                                     seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)

  if (is.null(reps)) {
    reps <- rep_len(1, length(t))
  } else {
    stopifnot((length(reps) == 1) || (length(reps) == length(t)))
  }

  if (is.null(emp_dist)) {
    emp_dist <- fauxseq::Mm_multitissue_LD_ALF
  }

  if (is.list(emp_dist)) {
    pick_emp_dist <- sample(length(emp_dist), 1)
    emp_dist <- emp_dist[[pick_emp_dist]]
  }

  t <- rep(t, times = reps)
  N <- length(t)

  G_rhy <- (stats::runif(n_genes) <= rhy_frac)

  A <- matrix((1 + min_A_effect + stats::rexp(n_genes, A_spread)),
              nrow = n_genes)
  A[!G_rhy] <- 1.0
  phi <- matrix(stats::runif(n_genes) * 2 * pi * G_rhy, nrow = n_genes)
  phi[!G_rhy] <- 0
  params <- data.frame(id = paste0("g", seq(n_genes)),
                       A = log2(A),
                       phi = phi)

  t_pattern <- cos(phi) %*% cos(2 * pi * t / period) +
    sin(phi) %*% sin(2 * pi * t / period)

  lib_size_fct <- stats::runif(N, min = lib_size_var[1], max = lib_size_var[2])
  lib_size <- lib_size_fct * depth

  draw <- sample(nrow(emp_dist), n_genes, replace = TRUE)
  lambda <- matrix(emp_dist$mu[draw], nrow = n_genes)
  params$mu <- emp_dist$mu[draw]
  size <- matrix(emp_dist$size[draw], nrow = n_genes)
  params$dispersion <- 1 / emp_dist$size[draw]
  size <- size[, rep(1, length(t))]

  lambda <- lambda[, rep(1, N)] * A[, rep(1, N)]^t_pattern

  mu <- sweep(lambda, 2, colSums(lambda) / lib_size, FUN = "/")

  counts <- matrix(
    mapply(\(x, y) stats::rnbinom(1, size = x, mu = y), size, mu),
    nrow = n_genes,
    byrow = FALSE,
    dimnames = list(
      id = paste0("g", seq(n_genes)),
      samples = make.names(replicate(N, paste0(
        sample(letters, 5), collapse = ""
      )), unique = TRUE)
    )
  )

  exp_design <- data.frame(time = t, row.names = colnames(counts))

  list(counts = counts,
       params = params,
       exp_design = exp_design)
}

#' Generate rhythmic RNA-seq count data under two conditions
#'
#' This function generates artificial timeseries transcriptomic data under two
#' conditions, which includes both differentially rhythmic and differential
#' expressed genes, as well as empirically estimated variability between
#' replicates/samples.
#'
#' @param t A numeric vector, time points at which samples are generated, by
#'   default seq(0, 20, 4)
#' @param reps NULL or numeric scalar or numeric vector with the same length as
#'   `t`. number of replicates at each time point, by default 1
#' @param period A numeric scalar, period of the rhythmic genes, by default 24
#' @param n_genes An integer, number of genes in the dataset, by default 10000
#' @param rhy_frac A fraction between 0 and 1, fraction of genes that are
#'   rhythmic, by default 0.1
#' @param DR_probs A vector of length 4. Determines the relative average number
#'   of  ("same", "gain", "loss", "change") group elements. The larger the
#'   numbers more concentrated are th actual number around the averages, by
#'   default rep(1.25, 4)
#' @param min_A_effect A positive numeric scalar, minimum amplitude (in log2
#'   fold) of rhythmic genes, by default 0.5
#' @param A_spread A positive numeric scalar, mean of the exponential
#'   distribution of rhythmic gene amplitudes (in log2 fold), by default 1
#' @param DE_frac A fraction between 0 and 1, fraction of differentially
#'   expressed (DE) genes, by default 0.1
#' @param min_DE_effect A positive numeric scalar, minimum log2 fold change in
#'   expression of DE genes, by default 0.5
#' @param DE_spread A positive numeric scalar, mean of the exponential
#'   distribution of DE fold changes, by default 1.0
#' @param groups A character vector of length 2, labels for the two
#'   groups/conditions, by default c("ctrl", "expt")
#' @param emp_dist A data.frame with two columns containing the empirical mu
#'   (means) and empirical size (reciprocal of dispersion), by default read from
#'   one of 5 mouse tissues at random.
#' @param depth A positive integer, average sequencing depth of different
#'   samples, by default 4e7
#' @param lib_size_var A numeric vector of two positive elements, window of
#'   variability of sequencing depth of samples about 'depth', by default c(0.8,
#'   1.2)
#' @param seed A positive integer, seed to ensure reproducible datasets, by
#'   default NULL
#'
#' @returns A list containing a integer count matrix, a data.frame of feature
#'   properties and a data.frame of the sample metadata.
#' @export
#'
generate_diffrhythmic_rnaseq <- function(t = seq(0, 20, 4),
                                         reps = NULL,
                                         period = 24,
                                         n_genes = 10000,
                                         rhy_frac = 0.1,
                                         DR_probs = rep(1.25, 4),
                                         min_A_effect = 0.5,
                                         A_spread = 1,
                                         DE_frac = 0.1,
                                         min_DE_effect = 0.5,
                                         DE_spread = 1,
                                         groups = c("ctrl", "expt"),
                                         emp_dist = NULL,
                                         depth = 4e7,
                                         lib_size_var = c(0.8, 1.2),
                                         seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)

  if (is.null(reps)) {
    reps <- rep_len(1, length(t))
  } else {
    stopifnot((length(reps) == 1) || (length(reps) == length(t)))
  }

  if (is.null(emp_dist)) {
    emp_dist <- fauxseq::Mm_multitissue_LD_ALF
  }

  if (is.list(emp_dist)) {
    pick_emp_dist <- sample(length(emp_dist), 1)
    emp_dist <- emp_dist[[pick_emp_dist]]
  }

  t <- rep(t, times = reps)

  N <- length(t)

  G_rhy <- (stats::runif(n_genes) <= rhy_frac)

  if (sum(G_rhy) > 0) {
    prior_p <- stats::rgamma(4, shape = DR_probs, rate = 1)
    prior_p <- prior_p / sum(prior_p)
    DR_counts <- stats::rmultinom(1, sum(G_rhy), prior_p)
    DR_groups <- c("same", "gain", "loss", "change")[rep.int(seq(4), DR_counts)]

    A <- matrix((1 + min_A_effect + stats::rexp(n_genes * 2, A_spread)),
                nrow = n_genes)
    A[!G_rhy, ] <- 1.0
    phi <- matrix(stats::runif(n_genes * 2) * 2 * pi * G_rhy, nrow = n_genes)
    phi[!G_rhy, ] <- 0

    G_rhy_index <- which(G_rhy)
    for (i in seq_along(DR_groups)) {
      if (DR_groups[i] == "gain") {
        A[G_rhy_index[i], 1] <- 1.0
        phi[G_rhy_index[i], 1] <- 0.0
      } else if (DR_groups[i] == "loss") {
        A[G_rhy_index[i], 2] <- 1.0
        phi[G_rhy_index[i], 2] <- 0.0
      } else if (DR_groups[i] == "same") {
        A[G_rhy_index[i], 1] <- A[G_rhy_index[i], 2]
        phi[G_rhy_index[i], 1] <- phi[G_rhy_index[i], 2]
      }
    }
    params <- data.frame(
      id = paste0("g", seq(n_genes)),
      A_1 = log2(A[, 1]),
      A_2 = log2(A[, 2]),
      phi_1 = phi[, 1],
      phi_2 = phi[, 2]
    )
    colnames(params) <- c("id",
                          paste0(rep(c("A", "phi"), each = 2),
                                 "_",
                                 rep(groups, 2)))

    params <- merge(
      params,
      data.frame(id = paste0("g", which(G_rhy)), category = DR_groups),
      all.x = TRUE,
      by = "id"
    )

    t_pattern <- cbind(
      cos(phi[, 1, drop = FALSE]) %*% cos(2 * pi * t / period) +
        sin(phi[, 1, drop = FALSE]) %*% sin(2 * pi * t / period),
      cos(phi[, 2, drop = FALSE]) %*% cos(2 * pi * t / period) +
        sin(phi[, 2, drop = FALSE]) %*% sin(2 * pi * t / period)
    )
  } else {
    params <- data.frame(matrix(0, nrow = n_genes, ncol = 6))
    colnames(params) <- c("id", paste0(rep(c("A", "phi"), each = 2),
                                       "_", rep(groups, 2)),
                          "category")
    params$id <- paste0("g", seq(n_genes))
    params$category <- NA
  }

  G_DE <- (stats::runif(n_genes) <= DE_frac)
  DE_effects <- matrix((1 + min_DE_effect + stats::rexp(n_genes, DE_spread)),
                       nrow = n_genes)^sign(2 * stats::runif(n_genes) - 1)
  DE_effects[!G_DE] <- 1.0
  params_de <- data.frame(id = paste0("g", seq(n_genes)),
                          DE_effect = ifelse(G_DE, log2(DE_effects[G_DE]), NA))

  params <- merge(params, params_de, by = "id", all = TRUE)
  params <- params[order(as.numeric(gsub("g", "", params$id))), ]

  lib_size_fct <- stats::runif(2 * N, min = lib_size_var[1],
                               max = lib_size_var[2])
  lib_size <- lib_size_fct * depth

  draw <- sample(nrow(emp_dist), n_genes, replace = TRUE)

  size <- matrix(emp_dist$size[draw], nrow = n_genes)
  params$dispersion <- 1 / emp_dist$size[draw]
  size <- size[, rep(1, 2 * N)]

  lambda <- matrix(emp_dist$mu[draw], nrow = n_genes)
  params$mu <- emp_dist$mu[draw]
  DE_effects <- cbind(rep(1, n_genes), DE_effects)

  if (sum(G_rhy) > 0) {
    if (sum(G_DE) > 0) {
      lambda <- (lambda[, rep(1, 2 * N)] *
                   DE_effects[, rep(c(1, 2), each = N)] *
                   A[, rep(c(1, 2), each = N)]^t_pattern)
    } else {
      lambda <- lambda[, rep(1, 2 * N)] * A[, rep(c(1, 2), each = N)]^t_pattern
    }
  } else {
    if (sum(G_DE) > 0) {
      lambda <- lambda[, rep(1, 2 * N)] * DE_effects[, rep(c(1, 2), each = N)]
    } else {
      lambda <- lambda[, rep(1, 2 * N)]
    }
  }

  mu <- sweep(lambda, 2, colSums(lambda) / lib_size, FUN = "/")

  counts <- matrix(
    mapply(\(x, y) stats::rnbinom(1, size = x, mu = y), size, mu),
    nrow = n_genes,
    byrow = FALSE,
    dimnames = list(
      id = paste0("g", seq(n_genes)),
      samples = make.names(replicate(2 * N, paste0(
        sample(letters, 5), collapse = ""
      )), unique = TRUE)
    )
  )

  exp_design <- data.frame(
    time = rep(t, 2),
    group = rep(groups, each = N),
    row.names = colnames(counts)
  )

  list(counts = counts,
       params = params,
       exp_design = exp_design)
}
