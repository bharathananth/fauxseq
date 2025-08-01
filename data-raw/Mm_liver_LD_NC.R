library(compareRhythms)

data("high_fat_diet_rnaseq")

nc_data <- high_fat_diet_rnaseq[, grepl("^NC", colnames(high_fat_diet_rnaseq))]
nc_data <- nc_data[rowMeans(nc_data) >= 10, ]

Mm_liver_LD_NC <- estimate_disp_dist(data.matrix(nc_data))

usethis::use_data(Mm_liver_LD_NC, overwrite = TRUE)
