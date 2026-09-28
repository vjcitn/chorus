test_that("run_alphagenome_mqtl errors clearly on a missing mqtl_csv", {
  expect_error(
    run_alphagenome_mqtl("does-not-exist.csv", repo_dir = "."),
    "does not exist"
  )
})

test_that("run_alphagenome_mqtl builds the expected script call and reads results back", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)

  tmp_repo <- tempfile()
  dir.create(tmp_repo)
  dir.create(file.path(tmp_repo, "chorus"))
  file.create(file.path(tmp_repo, "score_mqtls.py"))

  tmp_out <- tempfile(fileext = ".csv")
  fake_results <- data.frame(mqtl_id = "mqtl_0", effect_percentile = 0.9)

  captured <- NULL
  local_mocked_bindings(
    .run_chorus_script = function(script, script_args,
                                   mamba_env = "chorus",
                                   mamba_bin = "mamba") {
      captured <<- list(script = script, args = script_args,
                         mamba_env = mamba_env)
      utils::write.csv(fake_results, tmp_out, row.names = FALSE)
      list(status = 0L, stdout = "ok", stderr = character(0),
           command = mamba_bin, args = script_args, log_file = tempfile())
    }
  )

  res <- run_alphagenome_mqtl(
    tmp_csv,
    output = tmp_out,
    flip_beta = TRUE,
    device = "cpu",
    cell_types = c("B cell", "neutrophil"),
    repo_dir = tmp_repo
  )

  expect_equal(res$status, 0L)
  expect_true(is.data.frame(res$results))
  expect_equal(res$results$mqtl_id, "mqtl_0")
  expect_false(is.null(res$log_file))

  expect_true(grepl("score_mqtls\\.py$", captured$script))
  expect_true("--flip-beta" %in% captured$args)
  expect_true("--device" %in% captured$args)
  expect_equal(captured$args[which(captured$args == "--device") + 1], "cpu")
  expect_true("--cell-types" %in% captured$args)
  expect_true(all(c("B cell", "neutrophil") %in% captured$args))
  expect_true("-o" %in% captured$args)
  expect_false("--resume" %in% captured$args)
})

test_that("run_alphagenome_mqtl passes --resume only when requested", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)

  tmp_repo <- tempfile()
  dir.create(tmp_repo)
  dir.create(file.path(tmp_repo, "chorus"))
  file.create(file.path(tmp_repo, "score_mqtls.py"))

  captured <- NULL
  local_mocked_bindings(
    .run_chorus_script = function(script, script_args, ...) {
      captured <<- script_args
      list(status = 0L, stdout = "ok", stderr = character(0),
           command = "mamba", args = script_args, log_file = tempfile())
    }
  )

  suppressWarnings(
    run_alphagenome_mqtl(tmp_csv, resume = TRUE, repo_dir = tmp_repo)
  )
  expect_true("--resume" %in% captured)
})

test_that("run_alphagenome_mqtl leaves results NULL on a nonzero exit status", {
  tmp_csv <- tempfile(fileext = ".csv")
  writeLines(c("chrom,snp_pos,ref,alt,cpg_pos,beta",
               "chr1,100,A,G,105,0.3"), tmp_csv)

  tmp_repo <- tempfile()
  dir.create(tmp_repo)
  dir.create(file.path(tmp_repo, "chorus"))
  file.create(file.path(tmp_repo, "score_mqtls.py"))

  local_mocked_bindings(
    .run_chorus_script = function(...) {
      list(status = 1L, stdout = "boom", stderr = character(0),
           command = "mamba", args = character(0))
    }
  )

  res <- run_alphagenome_mqtl(tmp_csv, repo_dir = tmp_repo)
  expect_equal(res$status, 1L)
  expect_null(res$results)
  expect_true(any(grepl("boom", res$log)))
})
