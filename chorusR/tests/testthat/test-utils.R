test_that(".build_python_call assembles the direct-python invocation", {
  call <- .build_python_call("/opt/conda/envs/chorus/bin/python",
                              "score_mqtls.py", c("in.csv", "-o", "out.csv"))
  expect_equal(call$command, "/opt/conda/envs/chorus/bin/python")
  expect_equal(call$args, c("-u", "score_mqtls.py", "in.csv", "-o", "out.csv"))
})

test_that(".build_python_call validates its inputs", {
  expect_error(.build_python_call(c("a", "b"), "s.py", "x"))
  expect_error(.build_python_call("py", c("a", "b"), "x"))
})

test_that(".list_mamba_envs parses names and paths", {
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
  envs <- .list_mamba_envs("mamba")
  expect_equal(envs$name, c("base", "chorus", "chorus-alphagenome"))
  expect_equal(
    envs$path,
    c("/opt/conda", "/opt/conda/envs/chorus", "/opt/conda/envs/chorus-alphagenome")
  )
})

test_that(".list_mamba_envs errors clearly when the env list command fails", {
  local_mocked_bindings(
    system2 = function(...) {
      out <- character(0)
      attr(out, "status") <- 1L
      out
    },
    .package = "base"
  )
  expect_error(.list_mamba_envs("mamba"), "env list' failed")
})

test_that(".resolve_env_python finds the interpreter for an existing environment", {
  tmp <- tempfile()
  dir.create(file.path(tmp, "bin"), recursive = TRUE)
  file.create(file.path(tmp, "bin", "python"))
  local_mocked_bindings(
    .list_mamba_envs = function(...) {
      data.frame(name = "chorus", path = tmp, stringsAsFactors = FALSE)
    }
  )
  expect_equal(
    .resolve_env_python("chorus", "mamba"),
    file.path(tmp, "bin", "python")
  )
})

test_that(".resolve_env_python errors clearly when the environment does not exist", {
  local_mocked_bindings(
    .list_mamba_envs = function(...) {
      data.frame(name = "base", path = "/opt/conda", stringsAsFactors = FALSE)
    }
  )
  expect_error(.resolve_env_python("chorus-nope", "mamba"), "does not exist")
})

test_that(".resolve_env_python errors clearly when the python binary is missing", {
  tmp <- tempfile()
  dir.create(tmp)
  local_mocked_bindings(
    .list_mamba_envs = function(...) {
      data.frame(name = "chorus", path = tmp, stringsAsFactors = FALSE)
    }
  )
  expect_error(.resolve_env_python("chorus", "mamba"), "is not there")
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
    suppressMessages(.run_chorus_script("s.py", "x")),
    "preflight check ran"
  )
})

test_that(".run_chorus_script skips the preflight check when check_env = FALSE", {
  local_mocked_bindings(
    .check_mamba_env = function(...) stop("preflight check should not run")
  )
  local_mocked_bindings(
    .resolve_env_python = function(...) "/opt/conda/envs/chorus/bin/python"
  )
  local_mocked_bindings(
    system2 = function(command, args, stdout, stderr, ...) {
      writeLines("ok", stdout)
      0L
    },
    .package = "base"
  )
  run <- suppressMessages(
    .run_chorus_script("s.py", "x", check_env = FALSE)
  )
  expect_equal(run$status, 0L)
  expect_equal(run$stdout, "ok")
  expect_equal(run$command, "/opt/conda/envs/chorus/bin/python")
})

test_that(".run_chorus_script runs the environment's python directly, not mamba run", {
  local_mocked_bindings(.check_mamba_env = function(...) TRUE)
  local_mocked_bindings(
    .resolve_env_python = function(...) "/opt/conda/envs/chorus/bin/python"
  )
  local_mocked_bindings(
    system2 = function(command, args, stdout, stderr, ...) {
      writeLines("ok", stdout)
      0L
    },
    .package = "base"
  )
  run <- suppressMessages(.run_chorus_script("s.py", "x"))
  expect_equal(run$command, "/opt/conda/envs/chorus/bin/python")
  expect_equal(run$args, c("-u", "s.py", "x"))
})

test_that(".run_chorus_script writes a live log file the caller can tail", {
  local_mocked_bindings(.check_mamba_env = function(...) TRUE)
  local_mocked_bindings(
    .resolve_env_python = function(...) "/opt/conda/envs/chorus/bin/python"
  )
  local_mocked_bindings(
    system2 = function(command, args, stdout, stderr, ...) {
      writeLines("ok", stdout)
      0L
    },
    .package = "base"
  )
  log_file <- tempfile(fileext = ".log")
  run <- suppressMessages(
    .run_chorus_script("s.py", "x", log_file = log_file)
  )
  expect_equal(run$log_file, log_file)
  expect_true(file.exists(log_file))
  expect_equal(readLines(log_file), "ok")
})

test_that(".read_mqtl_results returns NULL without a warning on nonzero status", {
  run <- list(status = 1L, log_file = tempfile())
  expect_null(.read_mqtl_results(tempfile(), run))
})

test_that(".read_mqtl_results warns and returns NULL on status 0 but a missing file", {
  run <- list(status = 0L, log_file = tempfile())
  expect_warning(
    res <- .read_mqtl_results(tempfile(fileext = ".csv"), run),
    "missing or empty"
  )
  expect_null(res)
})

test_that(".read_mqtl_results warns and returns NULL on status 0 but an empty file", {
  out <- tempfile(fileext = ".csv")
  file.create(out)
  run <- list(status = 0L, log_file = tempfile())
  expect_warning(
    res <- .read_mqtl_results(out, run),
    "missing or empty"
  )
  expect_null(res)
})

test_that(".read_mqtl_results reads a well-formed CSV without warning", {
  out <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(mqtl_id = "a", x = 1), out, row.names = FALSE)
  run <- list(status = 0L, log_file = tempfile())
  res <- .read_mqtl_results(out, run)
  expect_true(is.data.frame(res))
  expect_equal(res$mqtl_id, "a")
})

test_that(".run_chorus_script warns (but does not error) on a nonzero exit status", {
  local_mocked_bindings(.check_mamba_env = function(...) TRUE)
  local_mocked_bindings(
    .resolve_env_python = function(...) "/opt/conda/envs/chorus/bin/python"
  )
  local_mocked_bindings(
    system2 = function(command, args, stdout, stderr, ...) {
      writeLines(c("Traceback (most recent call last):", "ModuleNotFoundError"), stdout)
      1L
    },
    .package = "base"
  )
  expect_warning(
    run <- suppressMessages(.run_chorus_script("s.py", "x")),
    "chorus script failed"
  )
  expect_equal(run$status, 1L)
  expect_true(any(grepl("ModuleNotFoundError", run$stdout)))
})
