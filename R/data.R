#' Empirical mean-dispersion data from 5 different tissues
#'
#' A subset of tissues from Deota et al. data (adlib feeding) were selected
#' based on clustering of their empirical mean-dispersion estimates. All the
#' original count data were from paired end libraries run on Novaseq 6000.
#'
#' @format ## `Mm_multitissue_LD_ALF` A list of 5 data frames with 10,000 rows
#'   and 2 columns:
#' \describe{
#'   \item{mu}{Mean of the negative binomial distributed read counts}
#'   \item{size}{Reciprocal of the dispersion parameter}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190389>
"Mm_multitissue_LD_ALF"
