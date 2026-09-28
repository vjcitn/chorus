test_that("run_cherimoya_mqtl errors clearly on a missing mqtl_csv", {
  expect_error(
    run_cherimoya_mqtl("does-not-exist.csv", repo_dir = "."),
    "does not exist"
  )
})

test_that("run_cherimoya_mqtl validates the assays argument", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)
  expect_error(
    run_cherimoya_mqtl(tmp_csv, assays = "CHIP", repo_dir = "."),
    "DNASE.*ATAC"
  )
})

test_that("run_cherimoya_mqtl builds the expected script call and reads results back", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)

  tmp_repo <- tempfile()
  dir.create(tmp_repo)
  dir.create(file.path(tmp_repo, "chorus"))
  file.create(file.path(tmp_repo, "score_mqtls_cherimoya.py"))

  tmp_out <- tempfile(fileext = ".csv")
  tmp_bg <- tempfile()
  fake_results <- data.frame(mqtl_id = "mqtl_0", track_id = "ATAC:ENCSR000EOT",
                              effect_percentile = 0.8)

  captured <- NULL
  local_mocked_bindings(
    .run_chorus_script = function(script, script_args,
                                   mamba_env = "chorus",
                                   mamba_bin = "mamba") {
      captured <<- list(script = script, args = script_args)
      utils::write.csv(fake_results, tmp_out, row.names = FALSE)
      list(status = 0L, stdout = "ok", stderr = character(0),
           command = mamba_bin, args = script_args, log_file = tempfile())
    }
  )

  res <- run_cherimoya_mqtl(
    tmp_csv,
    output = tmp_out,
    biosamples = c("B cell", "CD14-positive monocyte"),
    assays = "ATAC",
    save_bedgraph = tmp_bg,
    repo_dir = tmp_repo
  )

  expect_equal(res$status, 0L)
  expect_true(is.data.frame(res$results))
  expect_equal(res$results$track_id, "ATAC:ENCSR000EOT")

  expect_true(grepl("score_mqtls_cherimoya\\.py$", captured$script))
  expect_true("--biosamples" %in% captured$args)
  expect_true(all(c("B cell", "CD14-positive monocyte") %in% captured$args))
  expect_true("--assays" %in% captured$args)
  expect_equal(captured$args[which(captured$args == "--assays") + 1], "ATAC")
  expect_true("--save-bedgraph" %in% captured$args)
  expect_equal(captured$args[which(captured$args == "--save-bedgraph") + 1], tmp_bg)
  # device is NULL by default for Cherimoya (no bf16 constraint), so
  # --device should not appear unless the caller asks for it.
  expect_false("--device" %in% captured$args)
  expect_false("--resume" %in% captured$args)
  expect_false(is.null(res$log_file))
})

test_that("run_cherimoya_mqtl passes --resume only when requested", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)

  tmp_repo <- tempfile()
  dir.create(tmp_repo)
  dir.create(file.path(tmp_repo, "chorus"))
  file.create(file.path(tmp_repo, "score_mqtls_cherimoya.py"))

  captured <- NULL
  local_mocked_bindings(
    .run_chorus_script = function(script, script_args, ...) {
      captured <<- script_args
      list(status = 0L, stdout = "ok", stderr = character(0),
           command = "mamba", args = script_args, log_file = tempfile())
    }
  )

  suppressWarnings(
    run_cherimoya_mqtl(tmp_csv, resume = TRUE, repo_dir = tmp_repo)
  )
  expect_true("--resume" %in% captured)
})

test_that("run_cherimoya_mqtl passes --device only when supplied", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)

  tmp_repo <- tempfile()
  dir.create(tmp_repo)
  dir.create(file.path(tmp_repo, "chorus"))
  file.create(file.path(tmp_repo, "score_mqtls_cherimoya.py"))

  captured <- NULL
  local_mocked_bindings(
    .run_chorus_script = function(script, script_args, ...) {
      captured <<- script_args
      list(status = 0L, stdout = "ok", stderr = character(0),
           command = "mamba", args = script_args)
    }
  )

  suppressWarnings(
    run_cherimoya_mqtl(tmp_csv, device = "cuda", repo_dir = tmp_repo)
  )
  expect_true("--device" %in% captured)
  expect_equal(captured[which(captured == "--device") + 1], "cuda")
})
