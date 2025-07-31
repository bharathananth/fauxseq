#' Downsample all samples to have same number of reads as the smallest sample
#'
#' @param counts An integer matrix.
#'
#' @returns An integer matrix
#' @export
#'
downsample <- function(counts) {
  min_lib_size <- min(colSums(counts))
  downsampled_counts <- apply(counts, 2, drop_counts, min_lib_size)
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
