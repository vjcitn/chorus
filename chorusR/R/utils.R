#' Locate the chorus Python repository
#'
#' chorusR does not vendor the `chorus` Python package or its scripts; it
#' shells out to scripts that live in a checkout of the `chorus` repository
#' (<https://github.com/pinellolab/chorus>). This function resolves that
#' checkout's path, in order:
#'
#' 1. the `chorusR.repo_dir` R option (set with
#'    `options(chorusR.repo_dir = "/path/to/chorus")`),
#' 2. the `CHORUS_REPO_DIR` environment variable,
#' 3. the current working directory, if it looks like a chorus checkout
#'    (i.e. it contains a `chorus/` subdirectory and a `score_mqtls.py`
#'    file).
#'
#' @return A normalized path to the chorus repository checkout.
#' @export
#' @examples
#' # Requires the CHORUS_REPO_DIR environment variable to point at a
#' # chorus repository checkout.
#' chorus_repo_dir()
chorus_repo_dir <- function() {
  opt <- getOption("chorusR.repo_dir")
  if (!is.null(opt)) {
    return(normalizePath(opt, mustWork = TRUE))
  }

  env <- Sys.getenv("CHORUS_REPO_DIR", unset = NA)
  if (!is.na(env) && nzchar(env)) {
    return(normalizePath(env, mustWork = TRUE))
  }

  candidate <- getwd()
  if (dir.exists(file.path(candidate, "chorus")) &&
      file.exists(file.path(candidate, "score_mqtls.py"))) {
    return(normalizePath(candidate, mustWork = TRUE))
  }

  stop(
    "Could not locate the chorus repository checkout. Set it with ",
    "options(chorusR.repo_dir = \"/path/to/chorus\") or the ",
    "CHORUS_REPO_DIR environment variable.",
    call. = FALSE
  )
}

#' Find a chorus script inside the repository checkout
#'
#' @param script_name File name of the script, e.g. `"score_mqtls.py"`.
#' @param repo_dir Chorus repository checkout; defaults to
#'   [chorus_repo_dir()].
#' @return Normalized path to the script.
#' @keywords internal
.find_chorus_script <- function(script_name, repo_dir = chorus_repo_dir()) {
  path <- file.path(repo_dir, script_name)
  if (!file.exists(path)) {
    stop(
      "Expected to find '", script_name, "' in the chorus repository at '",
      repo_dir, "', but it is not there. Is chorusR.repo_dir / ",
      "CHORUS_REPO_DIR pointing at the right checkout?",
      call. = FALSE
    )
  }
  normalizePath(path, mustWork = TRUE)
}

#' List conda/mamba environments, with their names and install paths
#'
#' Shared parsing for [.check_mamba_env()] and [.resolve_env_python()], both
#' of which need `mamba env list`'s output -- one just the names, the other
#' also the install path.
#'
#' @param mamba_bin Name or path of the mamba/conda executable.
#' @return A data frame with columns `name` and `path`.
#' @keywords internal
.list_mamba_envs <- function(mamba_bin = "mamba") {
  envs <- suppressWarnings(
    system2(mamba_bin, c("env", "list"), stdout = TRUE, stderr = TRUE)
  )
  status <- attr(envs, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      "'", mamba_bin, " env list' failed (exit status ", status, "):\n",
      paste(envs, collapse = "\n"),
      call. = FALSE
    )
  }

  # Each non-comment, non-blank line is "<name>  [*]  <path>"; strip the
  # active-env marker before taking the first token as the name and the
  # last as the path.
  lines <- envs[!grepl("^\\s*#", envs) & nzchar(trimws(envs))]
  tokens <- lapply(strsplit(trimws(lines), "\\s+"), function(x) x[nzchar(x)])
  tokens <- tokens[lengths(tokens) > 0]
  data.frame(
    name = vapply(tokens, `[`, character(1), 1),
    path = vapply(tokens, function(x) x[length(x)], character(1)),
    stringsAsFactors = FALSE
  )
}

#' Resolve the python interpreter inside a conda/mamba environment
#'
#' chorusR used to invoke scripts via `mamba run -n <env> ...`, but that
#' wrapper's shell-script generation has, on different mamba/conda
#' versions and installs, produced a broken `exec -- ...` line -- bash's
#' `exec` builtin does not treat `--` as an end-of-options marker the way
#' most other commands do, so it fails with `exec: --: invalid option`
#' regardless of where `--no-capture-output` is placed on the `mamba run`
#' command line. Resolving and invoking the environment's own `python`
#' binary directly sidesteps `mamba run`'s wrapper entirely: the
#' interpreter's own path is enough to make it use that environment's
#' site-packages, without needing full shell activation.
#'
#' @param mamba_env Name of the conda/mamba environment.
#' @param mamba_bin Name or path of the mamba/conda executable.
#' @return Path to that environment's `python` executable.
#' @keywords internal
.resolve_env_python <- function(mamba_env = "chorus", mamba_bin = "mamba") {
  envs <- .list_mamba_envs(mamba_bin)
  match <- envs[envs$name == mamba_env, , drop = FALSE]
  if (!nrow(match)) {
    stop(
      "conda/mamba environment '", mamba_env, "' does not exist. ",
      "Available environments: ", paste(envs$name, collapse = ", "), ". ",
      "Set up the chorus environments first (see the chorus repository's ",
      "README/CLAUDE.md), or pass mamba_env= to point at the right one.",
      call. = FALSE
    )
  }
  python_bin <- file.path(match$path[1], "bin", "python")
  if (!file.exists(python_bin)) {
    stop(
      "Expected a python executable at '", python_bin, "' for conda/mamba ",
      "environment '", mamba_env, "', but it is not there.",
      call. = FALSE
    )
  }
  python_bin
}

#' Build the command-line invocation to run a chorus script directly
#'
#' Pure argument assembly, factored out of [.run_chorus_script()] so it can
#' be unit-tested without a mamba/conda installation or a Python
#' interpreter present. `-u` (unbuffered) defeats Python's own stdout
#' buffering when it isn't attached to a terminal (as it isn't here,
#' invoked via `system2()`), so progress prints reach `log_file` as they
#' happen rather than only once the script exits.
#'
#' @param python_bin Path to the python interpreter inside the target
#'   conda/mamba environment; see [.resolve_env_python()].
#' @param script Path to the Python script to run.
#' @param script_args Character vector of arguments to pass to the script.
#'
#' @return A list with `command` (the executable to run) and `args`
#'   (character vector of arguments), suitable for passing to
#'   [base::system2()] as `command` and `args`.
#' @keywords internal
.build_python_call <- function(python_bin, script, script_args) {
  stopifnot(
    is.character(python_bin), length(python_bin) == 1,
    is.character(script), length(script) == 1,
    is.character(script_args)
  )
  list(
    command = python_bin,
    args = c("-u", script, script_args)
  )
}

#' Check that mamba/conda and the target environment are actually usable
#'
#' Run before every script invocation so a bad Python configuration fails
#' immediately with a specific, actionable message, instead of surfacing
#' later as a status-127 `system2()` result or a buried Python traceback.
#'
#' @param mamba_env Name of the conda/mamba environment the caller wants
#'   to run in.
#' @param mamba_bin Name or path of the mamba/conda executable.
#' @return `TRUE`, invisibly, if `mamba_bin` resolves on `PATH` and
#'   `mamba_env` is one of its environments.
#' @keywords internal
.check_mamba_env <- function(mamba_env = "chorus", mamba_bin = "mamba") {
  if (!nzchar(Sys.which(mamba_bin))) {
    stop(
      "Could not find '", mamba_bin, "' on PATH. chorusR shells out to a ",
      "conda/mamba environment to run chorus -- install mamba (or conda) ",
      "and make sure it's on this R session's PATH, or pass a full path ",
      "via mamba_bin=.",
      call. = FALSE
    )
  }

  envs <- .list_mamba_envs(mamba_bin)

  if (!mamba_env %in% envs$name) {
    stop(
      "conda/mamba environment '", mamba_env, "' does not exist. ",
      "Available environments: ",
      paste(envs$name, collapse = ", "), ". ",
      "Set up the chorus environments first (see the chorus repository's ",
      "README/CLAUDE.md), or pass mamba_env= to point at the right one.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Run a chorus Python script via system2() and capture its output
#'
#' Redirects the child's combined stdout+stderr to `log_file` on disk
#' rather than capturing it only in memory, and only after the script
#' exits. Long chorus scoring runs (each row can be a multi-minute model
#' forward pass) print progress as they go -- see [.build_python_call()]
#' for how that reaches the log file unbuffered -- and writing straight
#' to disk means that progress is visible (`tail -f log_file` from
#' another terminal) and preserved even if this R session is killed
#' while `system2()` is still blocked on the child.
#'
#' Runs the target environment's `python` interpreter directly (see
#' [.resolve_env_python()]) rather than through `mamba run`/`conda run`,
#' whose wrapper script generation has proven unreliable across
#' mamba/conda versions.
#'
#' @param script Path to the Python script to run.
#' @param script_args Character vector of arguments to pass to the script.
#' @param mamba_env Name of the conda/mamba environment to run the script
#'   in. Defaults to `"chorus"` (the base environment); see
#'   [.resolve_env_python()].
#' @param mamba_bin Name or path of the mamba/conda executable.
#' @param log_file Path to write the child's combined stdout+stderr to,
#'   live, as it runs. Defaults to a fresh temp file; pass a stable path
#'   to keep it around after the call returns.
#' @param check_env Logical; run [.check_mamba_env()] first. Default
#'   `TRUE`; set `FALSE` only if the caller already checked (or is a test
#'   that wants to bypass it).
#' @return A list with `status` (integer exit code), `stdout` (the
#'   contents of `log_file`, read back in as a character vector once the
#'   script exits), `stderr` (always empty; stderr is merged into
#'   `log_file`/`stdout`), `command`, `args` (the invocation actually
#'   run, useful for debugging a failure), and `log_file` (the path
#'   itself).
#' @keywords internal
.run_chorus_script <- function(script,
                                script_args,
                                mamba_env = "chorus",
                                mamba_bin = "mamba",
                                log_file = tempfile(fileext = ".log"),
                                check_env = TRUE) {
  if (isTRUE(check_env)) {
    .check_mamba_env(mamba_env, mamba_bin)
  }

  python_bin <- .resolve_env_python(mamba_env, mamba_bin)
  call <- .build_python_call(python_bin, script, script_args)

  # chorusR shells out to a script living in a *separate* chorus repository
  # checkout rather than vendoring it, so "which version of chorusR" and
  # "which chorus commit is score_mqtls*.py actually running from" can
  # silently diverge -- e.g. this checkout ahead of or behind the one a fix
  # was tested against. Passed through as an env var (not a script_args
  # entry) so it reaches the log via the script's own startup banner
  # without changing either script's argument parsing.
  chorusr_version <- tryCatch(
    as.character(utils::packageVersion("chorusR")),
    error = function(e) "unknown (not installed as a package)"
  )
  message(
    "chorusR version: ", chorusr_version,
    " | python: ", call$command,
    " | script: ", script
  )
  message(
    "chorus script log: ", log_file,
    " (tail -f it from another terminal to watch progress)"
  )

  status <- system2(
    call$command, args = call$args,
    stdout = log_file, stderr = log_file,
    env = paste0("CHORUSR_VERSION=", chorusr_version)
  )
  if (is.null(status)) status <- 0L

  log <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)

  if (!identical(as.integer(status), 0L)) {
    tail_n <- min(length(log), 20L)
    warning(
      "chorus script failed (exit status ", status, "): ",
      call$command, " ", paste(call$args, collapse = " "),
      "\nLast ", tail_n, " line(s) of output (full log: ", log_file, "):\n",
      paste(utils::tail(log, tail_n), collapse = "\n"),
      call. = FALSE
    )
  }

  list(
    status = as.integer(status),
    stdout = log,
    stderr = character(0),
    command = call$command,
    args = call$args,
    log_file = log_file
  )
}

#' Read back a chorus script's results CSV, without crashing on a bad file
#'
#' A `run$status` of `0` from [.run_chorus_script()] is not by itself proof
#' that `output` holds a valid, non-empty results table: `mamba run`/`conda
#' run` are known to occasionally mis-report a wrapped process's exit code
#' (see this package's parent repository's CLAUDE.md), and a long scoring
#' run that is OOM-killed or otherwise dies mid-write can leave a
#' zero-byte or header-only CSV behind despite that. Reading such a file
#' with [utils::read.csv()] directly raises a generic, unhelpful
#' `"no lines available in input"` error. This wraps that read so a bad
#' file instead produces `NULL` results and a warning that points at
#' `run$log_file` for diagnosis.
#'
#' @param output Path to the results CSV written by the chorus script.
#' @param run The list returned by [.run_chorus_script()].
#' @return A data frame of results, or `NULL` if `output` does not exist,
#'   is empty, or could not be parsed as CSV.
#' @keywords internal
.read_mqtl_results <- function(output, run) {
  if (!identical(run$status, 0L)) {
    return(NULL)
  }
  if (!file.exists(output) || file.info(output)$size == 0) {
    warning(
      "chorus script exited with status 0 but '", output, "' is missing ",
      "or empty -- the underlying run likely died partway through (e.g. ",
      "killed for memory) without that being reflected in the exit code. ",
      "Check the log for what actually happened: ", run$log_file,
      call. = FALSE
    )
    return(NULL)
  }
  results <- tryCatch(
    utils::read.csv(output, stringsAsFactors = FALSE),
    error = function(e) {
      warning(
        "chorus script exited with status 0, but '", output, "' could not ",
        "be read as CSV (", conditionMessage(e), "). Check the log: ",
        run$log_file,
        call. = FALSE
      )
      NULL
    }
  )
  results
}
