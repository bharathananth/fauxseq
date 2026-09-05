library(quickcheck)

# Tests for generate_rhythmic_rnaseq

numeric_df <- data_frame_(
  mu = double_bounded(1, 1e4),
  size = double_bounded(1, 100),
  rows = c(100, 1000)
)

test_that("counts have n_genes rows, n_genes cols, and are positive integers", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    reps = integer_bounded(1, 4, len = 1),
    t = numeric_bounded(0, 100, len = c(4, 8)),
    emp_dist = one_of(
      constant("liver"), constant("multitissue"),
      numeric_df, list_of(numeric_df)
    ),
    property = function(n_genes, reps, t, emp_dist) {
      sim_data <- generate_rhythmic_rnaseq(
        t = t,
        reps = reps,
        n_genes = n_genes,
        emp_dist = emp_dist
      )
      expect_shape(sim_data$counts,
        dim = c(n_genes, length(t) * reps)
      )
      expect_true(all(sim_data$count >= 0))
      expect_true(all(is.wholenumber(sim_data$count)))
    }
  )
})

test_that("counts have unique colnames and rownames of params match colnames of
          count", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    reps = integer_bounded(1, 4, len = 1),
    t = numeric_bounded(0, 100, len = c(4, 8)),
    property = function(n_genes, reps, t) {
      sim_data <- generate_rhythmic_rnaseq(
        t = t,
        reps = reps,
        n_genes = n_genes
      )
      expect_equal(nrow(sim_data$params), nrow(sim_data$counts))
      expect_equal(anyDuplicated(colnames(sim_data$counts)), 0)
    }
  )
})

test_that("parameters satisfy conditions", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    rhy_frac = double_bounded(0, 1.0, len = 1),
    min_A_effect = double_bounded(0, 2.0, len = 1),
    mean_disp_trend = logical_(len = 1),
    property = function(n_genes, rhy_frac, min_A_effect, mean_disp_trend) {
      sim_data <- generate_rhythmic_rnaseq(
        n_genes = n_genes,
        rhy_frac = rhy_frac,
        min_A_effect = min_A_effect,
        mean_disp_trend = mean_disp_trend
      )
      params <- sim_data$params
      expect_shape(params, dim = c(n_genes, 5))
      expect_identical(rownames(sim_data$counts), params$id)
      expect_all_true(ifelse(params$A > 0,
        params$A >= min_A_effect, TRUE
      ))
      expect_all_true(params$phi >= 0 & params$phi <= 2 * pi)
    }
  )
  test_that("seeding produces identical objects", {
    for_all(
      n_genes = integer_bounded(2, 20000, len = 1),
      rhy_frac = double_bounded(0, 1.0, len = 1),
      min_A_effect = double_bounded(0, 2.0, len = 1),
      reps = integer_bounded(1, 8, len = 1),
      seed = integer_bounded(0, 2^16 - 1, len = 1),
      property = function(n_genes, rhy_frac, min_A_effect, reps, seed) {
        data_1 <- generate_rhythmic_rnaseq(
          n_genes = n_genes,
          rhy_frac = rhy_frac,
          min_A_effect = min_A_effect,
          reps = reps,
          seed = seed
        )
        data_2 <- generate_rhythmic_rnaseq(
          n_genes = n_genes,
          rhy_frac = rhy_frac,
          min_A_effect = min_A_effect,
          reps = reps,
          seed = seed
        )
        expect_identical(data_1, data_2)
      }
    )
  })
})

# Tests for generate_diffrhythmic_rnaseq

test_that("counts have n_genes rows, n_genes cols, and are positive integers", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    reps = integer_bounded(1, 4, len = 1),
    t = numeric_bounded(0, 100, len = c(4, 8)),
    emp_dist = one_of(
      constant("liver"), constant("multitissue"),
      numeric_df, list_of(numeric_df)
    ),
    property = function(n_genes, reps, t, emp_dist) {
      sim_data <- generate_diffrhythmic_rnaseq(
        t = t,
        reps = reps,
        n_genes = n_genes,
        emp_dist = emp_dist
      )
      expect_shape(sim_data$counts,
        dim = c(n_genes, 2 * length(t) * reps)
      )
      expect_true(all(sim_data$count >= 0))
      expect_true(all(is.wholenumber(sim_data$count)))
    }
  )
})

test_that("counts have unique colnames and rownames of params match colnames of
          count", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    reps = integer_bounded(1, 4, len = 1),
    t = numeric_bounded(0, 100, len = c(4, 8)),
    property = function(n_genes, reps, t) {
      sim_data <- generate_rhythmic_rnaseq(
        t = t,
        reps = reps,
        n_genes = n_genes
      )
      expect_equal(nrow(sim_data$params), nrow(sim_data$counts))
      expect_equal(anyDuplicated(colnames(sim_data$counts)), 0)
    }
  )
})

test_that("parameters satisfy conditions", {
  for_all(
    n_genes = integer_bounded(2, 20000, len = 1),
    rhy_frac = double_bounded(0, 1.0, len = 1),
    DE_frac = double_bounded(0, 1.0, len = 1),
    min_A_effect = double_bounded(0, 2.0, len = 1),
    min_DE_effect = double_bounded(0, 2.0, len = 1),
    mean_disp_trend = logical_(len = 1),
    property = function(n_genes, rhy_frac, DE_frac,
                        min_A_effect, min_DE_effect,
                        mean_disp_trend) {
      sim_data <- generate_diffrhythmic_rnaseq(
        n_genes = n_genes,
        rhy_frac = rhy_frac,
        DE_frac = DE_frac,
        min_A_effect = min_A_effect,
        min_DE_effect = min_DE_effect,
        mean_disp_trend = mean_disp_trend
      )
      params <- sim_data$params
      expect_shape(params, dim = c(n_genes, 9))
      expect_identical(rownames(sim_data$counts), params$id)
      expect_in(params$category, c("same", "gain", "loss", "change", NA))
      expect_all_true(ifelse(!is.na(params$DE_effect), abs(params$DE_effect) >=
        min_DE_effect, TRUE))
      expect_all_true(ifelse(params$category %in% c("same", "loss", "change"),
        params$A_ctrl >= min_A_effect, TRUE
      ))
      expect_all_true(ifelse(params$category %in% c("same", "gain", "change"),
        params$A_expt >= min_A_effect, TRUE
      ))
      expect_all_true(ifelse(!is.na(params$category), (params$phi_ctrl >= 0) &
        (params$phi_ctrl <= 2 * pi), TRUE))
      expect_all_true(ifelse(!is.na(params$category), (params$phi_expt >= 0) &
        (params$phi_expt <= 2 * pi), TRUE))
      expect_all_true(xor(
        (params$A_ctrl == 0) & (params$A_expt == 0),
        !is.na(params$category)
      ))
    }
  )
  test_that("seeding produces identical objects", {
    for_all(
      n_genes = integer_bounded(2, 20000, len = 1),
      rhy_frac = double_bounded(0, 1.0, len = 1),
      DE_frac = double_bounded(0, 1.0, len = 1),
      min_A_effect = double_bounded(0, 2.0, len = 1),
      min_DE_effect = double_bounded(0, 2.0, len = 1),
      reps = integer_bounded(1, 8, len = 1),
      seed = integer_bounded(0, 2^16 - 1, len = 1),
      property = function(n_genes, rhy_frac, DE_frac, min_A_effect,
                          min_DE_effect, reps, seed) {
        data_1 <- generate_diffrhythmic_rnaseq(
          n_genes = n_genes,
          rhy_frac = rhy_frac,
          DE_frac = DE_frac,
          min_A_effect = min_A_effect,
          min_DE_effect = min_DE_effect,
          reps = reps,
          seed = seed
        )
        data_2 <- generate_diffrhythmic_rnaseq(
          n_genes = n_genes,
          rhy_frac = rhy_frac,
          DE_frac = DE_frac,
          min_A_effect = min_A_effect,
          min_DE_effect = min_DE_effect,
          reps = reps,
          seed = seed
        )
        expect_identical(data_1, data_2)
      }
    )
  })
})
