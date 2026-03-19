test_that("new_saturation_analyzer validates inputs", {
  expect_error(modsite:::new_saturation_analyzer(""), regexp = "non-empty")
  expect_error(modsite:::new_saturation_analyzer(NA_character_), regexp = "non-empty")
  expect_error(modsite:::new_saturation_analyzer("a.bam", min_mod_rate = 1.5), regexp = "0, 1")
})


test_that("simulate_downsampling works with in-memory site_stats fixture", {
  a <- modsite:::new_saturation_analyzer("fake.bam", min_depth = 10L, min_mod_rate = 0.05)
  a$site_stats <- data.frame(
    total_depth = c(100L, 100L, 5L),
    mod_count = c(20L, 5L, 2L),
    stringsAsFactors = FALSE
  )

  sim <- modsite:::simulate_downsampling(a, ratios = c(0.5, 1.0), random_seed = 1L)
  expect_true(all(c("ratio", "valid_sites") %in% colnames(sim)))
  expect_equal(nrow(sim), 2L)
  expect_equal(sim$ratio[2], 1.0)
  # At 100%: site1 rate=0.2 ok; site2 rate=0.05 ok; site3 depth=5 < 10
  expect_equal(sim$valid_sites[2], 2L)
})


test_that("plot_saturation writes a PNG file", {
  a <- modsite:::new_saturation_analyzer("fake.bam")
  sim <- data.frame(ratio = c(0.1, 1.0), valid_sites = c(1L, 10L))
  out <- tempfile(fileext = ".png")
  expect_silent(modsite:::plot_saturation(a, sim, output_file = out))
  expect_true(file.exists(out))
})

