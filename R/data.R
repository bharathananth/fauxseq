#' Empirical mean-dispersion data from 5 different tissues
#'
#' A subset of tissues from Deota et al. data (adlib feeding) were selected
#' based on clustering of their empirical mean-dispersion estimates. All the
#' original count data were generated from paired end libraries run on Novaseq
#' 6000.
#'
#' @format ## `Mm_multitissue_LD_ALF` A list of 5 data frames with 10,000 rows
#'   and 2 columns:
#' \describe{
#'   \item{mu}{Mean of the negative binomial distributed read counts}
#'   \item{size}{Reciprocal of the dispersion parameter}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190389>
"Mm_multitissue_LD_ALF"

#' Empirical mean-dispersion data from mouse liver
#'
#' The control data from Quagliarini et al. data (normal chow) used to study the
#' effect of high-fat diet on the liver transcriptome. The original count data
#' were from paired end libraries run on HiSeq 2500.
#'
#' @format ## `Mm_liver_LD_NC` A data frame with 10,000 rows and 2 columns:
#' \describe{
#'   \item{mu}{Mean of the negative binomial distributed read counts}
#'   \item{size}{Reciprocal of the dispersion parameter}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108688>
"Mm_liver_LD_NC"
