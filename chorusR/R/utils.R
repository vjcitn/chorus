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
#' \dontrun{
#' options(chorusR.repo_dir = "/path/to/your/chorus/checkout")
#' chorus_repo_dir()
#' }
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

#' Build the command-line invocation for a chorus mamba-run script call
#'
#' Pure argument assembly, factored out of [.run_chorus_script()] so it can
#' be unit-tested without a mamba/conda installation or a Python
#' interpreter present.
#'
#' @param script Path to the Python script to run.
#' @param script_args Character vector of arguments to pass to the script.
#' @param mamba_env Name of the conda/mamba environment to run the script
#'   in (chorus isolates each oracle's dependencies per environment; the
#'   mQTL-scoring scripts themselves only need the base environment
#'   because they call oracles with `use_environment=True`, which spawns
#'   the oracle-specific environment as a subprocess internally).
#' @param mamba_bin Name or path of the mamba/conda executable.
#' @return A list with `command` (the executable to run) and `args`
#'   (character vector of arguments), suitable for passing to
#'   [base::system2()] as `command` and `args`.
#' @keywords internal
.build_mamba_run_call <- function(script,
                                   script_args,
                                   mamba_env = "chorus",
                                   mamba_bin = "mamba") {
  stopifnot(
    is.character(script), length(script) == 1,
    is.character(script_args),
    is.character(mamba_env), length(mamba_env) == 1,
    is.character(mamba_bin), length(mamba_bin) == 1
  )
  list(
    command = mamba_bin,
    args = c("run", "-n", mamba_env, "python", script, script_args)
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
  # active-env marker before taking the first token as the env name.
  lines <- envs[!grepl("^\\s*#", envs) & nzchar(trimws(envs))]
  names <- vapply(strsplit(trimws(lines), "\\s+"), `[`, character(1), 1)
  names <- names[nzchar(names)]

  if (!mamba_env %in% names) {
    stop(
      "conda/mamba environment '", mamba_env, "' does not exist. ",
      "Available environments: ",
      paste(names, collapse = ", "), ". ",
      "Set up the chorus environments first (see the chorus repository's ",
      "README/CLAUDE.md), or pass mamba_env= to point at the right one.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Run a chorus Python script via system2() and capture its output
#'
#' @param script Path to the Python script to run.
#' @param script_args Character vector of arguments to pass to the script.
#' @param mamba_env Name of the conda/mamba environment to run the script
#'   in. Defaults to `"chorus"` (the base environment); see
#'   [.build_mamba_run_call()].
#' @param mamba_bin Name or path of the mamba/conda executable.
#' @param stdout,stderr Passed through to [base::system2()]; default
#'   `TRUE` for both, which per `system2()`'s own semantics merges stderr
#'   into the captured stdout character vector (so `stderr` in the
#'   returned value is always empty in that default case).
#' @param check_env Logical; run [.check_mamba_env()] first. Default
#'   `TRUE`; set `FALSE` only if the caller already checked (or is a test
#'   that wants to bypass it).
#' @return A list with `status` (integer exit code), `stdout` (character
#'   vector, combined stdout+stderr by default), `stderr` (character
#'   vector, only populated if `stderr` was given as a file path rather
#'   than `TRUE`), `command`, and `args` (the invocation actually run,
#'   useful for debugging a failure).
#' @keywords internal
.run_chorus_script <- function(script,
                                script_args,
                                mamba_env = "chorus",
                                mamba_bin = "mamba",
                                stdout = TRUE,
                                stderr = TRUE,
                                check_env = TRUE) {
  if (isTRUE(check_env)) {
    .check_mamba_env(mamba_env, mamba_bin)
  }

  call <- .build_mamba_run_call(script, script_args, mamba_env, mamba_bin)
  out <- system2(
    call$command, args = call$args,
    stdout = stdout, stderr = stderr
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  log <- if (is.character(out)) out else character(0)

  if (status != 0L) {
    tail_n <- min(length(log), 20L)
    warning(
      "chorus script failed (exit status ", status, "): ",
      call$command, " ", paste(call$args, collapse = " "),
      "\nLast ", tail_n, " line(s) of output:\n",
      paste(utils::tail(log, tail_n), collapse = "\n"),
      call. = FALSE
    )
  }

  list(
    status = status,
    stdout = log,
    stderr = character(0),
    command = call$command,
    args = call$args
  )
}
