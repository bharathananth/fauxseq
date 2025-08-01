#' Downsample all samples to have same number of reads as the smallest sample
#'
#' @param counts An integer matrix.
#' @param parallel A boolean.
#' @param ncores An integer.
#'
#' @returns An integer matrix
#' @export
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_apply
downsample <- function(counts,
                       parallel = TRUE,
                       ncores = NULL) {
  min_lib_size <- min(colSums(counts))
  old_plan <- future::plan()

  # Setup parallelism safely
  if (parallel) {
    if (!is.null(ncores)) {
      future::plan(future::multisession, workers = ncores)
    } else {
      future::plan(future::multisession)
    }
  } else {
    future::plan(future::sequential)
  }

  # Perform parallel apply
  downsampled_counts <- future.apply::future_apply(
    counts,
    MARGIN = 2,
    FUN = drop_counts,
    N = min_lib_size,
    future.packages = "fauxseq",
    future.seed = TRUE
  )

  # Reset plan to original (to not break parent/global state)
  future::plan(old_plan)

  return(downsampled_counts)
}

#' Reduce the total column of a sample by downsampling
#'
#' @param y An integer vector to downsample.
#' @param N An integer scalar.
#'
#' @returns An integer vector.
#'
drop_counts <- function(y, N) {
  total <- sum(y)

  if (total < N) {
    stop("N is larger than total number of objects.")
  }

  if (total == N) {
    return(y)
  } else {
    all_reads <- rep(seq_along(y), times = y)
    removed <- sample(all_reads, N)

    return(tabulate(removed, nbins = length(y)))
  }
}
