test_that("chorus_example_mqtls_path resolves to a real, well-formed CSV", {
  path <- chorus_example_mqtls_path()
  expect_true(file.exists(path))

  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_setequal(
    names(df),
    c("mqtl_id", "chrom", "snp_pos", "ref", "alt", "cpg_pos", "beta")
  )
  expect_equal(nrow(df), 4)
  expect_true(all(grepl("^chr", df$chrom)))
  expect_true(all(nchar(df$ref) == 1))
  expect_true(all(nchar(df$alt) == 1))
  expect_true(all(df$ref != df$alt))
  expect_true(is.numeric(df$snp_pos))
  expect_true(is.numeric(df$cpg_pos))
  expect_true(is.numeric(df$beta))

  # Documented in chorus_example_mqtls_path(): two rows are inside
  # Cherimoya's 2,114 bp window, two are outside it.
  dist <- abs(df$cpg_pos - df$snp_pos)
  expect_equal(sum(dist < 1000), 2)
  expect_equal(sum(dist >= 1000), 2)
})
