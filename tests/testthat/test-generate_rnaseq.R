library(quickcheck)

tests <- 20L

is.wholenumber <- function(x, tol = 2 * .Machine$double.eps)  abs(x - round(x)) < tol

test_that("counts have n_genes rows, n_genes cols, and are positive integers", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    reps = integer_bounded(1, 4, len = 1),
    t = numeric_bounded(0, 100, len = c(4, 8)),
    property = function(n_genes, reps, t) {
      sim_data = generate_rhythmic_rnaseq(t = t, reps = reps, n_genes=n_genes)
      expect_equal(n_genes, nrow(sim_data$counts))
      expect_equal(length(t)*reps, ncol(sim_data$counts))
      expect_true(all(sim_data$count>=0))
      expect_true(all(is.wholenumber(sim_data$count)))
    },
    tests = tests
  )
})
