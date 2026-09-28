test_that(".build_mamba_run_call assembles the mamba/conda invocation", {
  call <- .build_mamba_run_call("score_mqtls.py", c("in.csv", "-o", "out.csv"))
  expect_equal(call$command, "mamba")
  expect_equal(
    call$args,
    c("run", "-n", "chorus", "python", "score_mqtls.py",
      "in.csv", "-o", "out.csv")
  )
})

test_that(".build_mamba_run_call respects mamba_env and mamba_bin", {
  call <- .build_mamba_run_call("s.py", "x", mamba_env = "myenv",
                                 mamba_bin = "conda")
  expect_equal(call$command, "conda")
  expect_equal(call$args[1:3], c("run", "-n", "myenv"))
})

test_that(".build_mamba_run_call validates its inputs", {
  expect_error(.build_mamba_run_call(c("a", "b"), "x"))
  expect_error(.build_mamba_run_call("s.py", "x", mamba_env = 1))
})

test_that("chorus_repo_dir resolves from the chorusR.repo_dir option", {
  tmp <- tempfile()
  dir.create(tmp)
  dir.create(file.path(tmp, "chorus"))
  file.create(file.path(tmp, "score_mqtls.py"))

  old <- getOption("chorusR.repo_dir")
  options(chorusR.repo_dir = tmp)
  on.exit(options(chorusR.repo_dir = old), add = TRUE)

  expect_equal(chorus_repo_dir(), normalizePath(tmp))
})

test_that("chorus_repo_dir resolves from the CHORUS_REPO_DIR env var", {
  tmp <- tempfile()
  dir.create(tmp)
  dir.create(file.path(tmp, "chorus"))
  file.create(file.path(tmp, "score_mqtls.py"))

  old_opt <- getOption("chorusR.repo_dir")
  options(chorusR.repo_dir = NULL)
  old_env <- Sys.getenv("CHORUS_REPO_DIR", unset = NA)
  Sys.setenv(CHORUS_REPO_DIR = tmp)
  on.exit({
    options(chorusR.repo_dir = old_opt)
    if (is.na(old_env)) {
      Sys.unsetenv("CHORUS_REPO_DIR")
    } else {
      Sys.setenv(CHORUS_REPO_DIR = old_env)
    }
  }, add = TRUE)

  expect_equal(chorus_repo_dir(), normalizePath(tmp))
})

test_that("chorus_repo_dir errors clearly when nothing resolves", {
  old_opt <- getOption("chorusR.repo_dir")
  options(chorusR.repo_dir = NULL)
  old_env <- Sys.getenv("CHORUS_REPO_DIR", unset = NA)
  Sys.unsetenv("CHORUS_REPO_DIR")
  old_wd <- getwd()
  tmp <- tempfile()
  dir.create(tmp)
  setwd(tmp)
  on.exit({
    setwd(old_wd)
    options(chorusR.repo_dir = old_opt)
    if (!is.na(old_env)) Sys.setenv(CHORUS_REPO_DIR = old_env)
  }, add = TRUE)

  expect_error(chorus_repo_dir(), "Could not locate")
})

test_that(".find_chorus_script errors clearly when the script is missing", {
  tmp <- tempfile()
  dir.create(tmp)
  expect_error(
    .find_chorus_script("score_mqtls.py", repo_dir = tmp),
    "is not there"
  )
})

test_that(".find_chorus_script finds a script that is present", {
  tmp <- tempfile()
  dir.create(tmp)
  file.create(file.path(tmp, "score_mqtls.py"))
  expect_equal(
    .find_chorus_script("score_mqtls.py", repo_dir = tmp),
    normalizePath(file.path(tmp, "score_mqtls.py"))
  )
})

test_that(".check_mamba_env errors clearly when mamba_bin is not on PATH", {
  local_mocked_bindings(Sys.which = function(...) c(nosuchbin = ""), .package = "base")
  expect_error(
    .check_mamba_env(mamba_bin = "nosuchbin"),
    "Could not find 'nosuchbin' on PATH"
  )
})

test_that(".check_mamba_env errors clearly when the env list command fails", {
  local_mocked_bindings(Sys.which = function(...) c(mamba = "/usr/bin/mamba"), .package = "base")
  local_mocked_bindings(
    system2 = function(...) {
      out <- character(0)
      attr(out, "status") <- 1L
      out
    },
    .package = "base"
  )
  expect_error(
    .check_mamba_env(mamba_bin = "mamba"),
    "env list' failed"
  )
})

test_that(".check_mamba_env errors clearly when mamba_env does not exist", {
  local_mocked_bindings(Sys.which = function(...) c(mamba = "/usr/bin/mamba"), .package = "base")
  local_mocked_bindings(
    system2 = function(...) {
      c(
        "# conda environments:",
        "#",
        "base                     /opt/conda",
        "chorus                *  /opt/conda/envs/chorus",
        "chorus-alphagenome       /opt/conda/envs/chorus-alphagenome"
      )
    },
    .package = "base"
  )
  expect_error(
    .check_mamba_env(mamba_env = "chorus-nope", mamba_bin = "mamba"),
    "does not exist"
  )
})

test_that(".check_mamba_env succeeds when mamba_bin and mamba_env are both fine", {
  local_mocked_bindings(Sys.which = function(...) c(mamba = "/usr/bin/mamba"), .package = "base")
  local_mocked_bindings(
    system2 = function(...) {
      c(
        "# conda environments:",
        "#",
        "base                     /opt/conda",
        "chorus                *  /opt/conda/envs/chorus"
      )
    },
    .package = "base"
  )
  expect_true(.check_mamba_env(mamba_env = "chorus", mamba_bin = "mamba"))
})

test_that(".run_chorus_script runs the preflight check by default", {
  local_mocked_bindings(
    .check_mamba_env = function(...) stop("preflight check ran")
  )
  expect_error(
    .run_chorus_script("s.py", "x"),
    "preflight check ran"
  )
})

test_that(".run_chorus_script skips the preflight check when check_env = FALSE", {
  local_mocked_bindings(
    .check_mamba_env = function(...) stop("preflight check should not run")
  )
  local_mocked_bindings(
    system2 = function(...) {
      out <- "ok"
      attr(out, "status") <- 0L
      out
    },
    .package = "base"
  )
  run <- .run_chorus_script("s.py", "x", check_env = FALSE)
  expect_equal(run$status, 0L)
})

test_that(".run_chorus_script warns (but does not error) on a nonzero exit status", {
  local_mocked_bindings(.check_mamba_env = function(...) TRUE)
  local_mocked_bindings(
    system2 = function(...) {
      out <- c("Traceback (most recent call last):", "ModuleNotFoundError")
      attr(out, "status") <- 1L
      out
    },
    .package = "base"
  )
  expect_warning(
    run <- .run_chorus_script("s.py", "x"),
    "chorus script failed"
  )
  expect_equal(run$status, 1L)
  expect_true(any(grepl("ModuleNotFoundError", run$stdout)))
})
