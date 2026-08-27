## code to prepare `Mm_multitissue` dataset goes here
library(GEOquery)
library(tidyverse)
library(fauxseq)

geo <- getGEOSuppFiles("GSE190389")

data <- read_csv(rownames(geo)[1], locale = locale(encoding = "ISO-8859-1"))

data <- data[, grepl("ALF", colnames(data))]

tissues <- str_extract(colnames(data), "^[A-Za-z]+") |> unique()

tissue_data <- map(tissues, \(tis) {
  data[, grepl(tis, names(data))] |>
    filter(if_any(everything(), \(x) x > 0))
})

tissue_stats <- map(tissue_data, \(x) {
  data.matrix(x) |>
    estimate_disp_dist() |>
    slice_sample(n = 10000)
})


names(tissue_stats) <- tissues

# Mm_multitissue <- list_rbind(tissue_stats, names_to = "tissue")  # Complete data used for exploring mean-dispersion relationship

Mm_multitissue_LD_ALF <- tissue_stats[c("SCN", "Pancreas", "Lung", "Liver", "Muscle")]

names(Mm_multitissue_LD_ALF) <- c("SCN", "Pancreas", "Lung", "Liver", "Muscle")

usethis::use_data(Mm_multitissue_LD_ALF, overwrite = TRUE)
