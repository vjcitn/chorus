#' Score whole-blood mQTLs against AlphaGenome chromatin tracks
#'
#' Wraps `score_mqtls.py` (in the chorus repository root): for each
#' SNP-CpG pair in `mqtl_csv`, predicts the variant's effect on AlphaGenome
#' DNase/ATAC/CAGE/histone-ChIP tracks in matched blood cell types, and
#' reports whether the predicted chromatin-state change is concordant with
#' the mQTL's reported methylation direction. See
#' `alphagenome4mqtl.md` in the chorus repository for the full column
#' spec of `mqtl_csv` and the scoring methodology.
#'
#' AlphaGenome's JAX backend cannot run on pre-Ampere GPUs (Tesla T4,
#' V100, P100, RTX 20-series) -- see the chorus repo's
#' `CHORUS_RESOURCES_SUMMARY.md` for the underlying `ALG_DOT_BF16_BF16_F32`
#' error. `device` defaults to `"cpu"` for that reason.
#'
#' @param mqtl_csv Path to a CSV of SNP-CpG pairs. Required columns:
#'   `chrom`, `snp_pos`, `ref`, `alt`, `cpg_pos`, `beta`. Optional:
#'   `mqtl_id`, `gene_symbol`. See `alphagenome4mqtl.md` for the exact
#'   conventions (1-based positions, hg38, beta sign convention).
#' @param output Path to write the long-format results CSV (one row per
#'   mqtl x track). Defaults to a temp file; read back into R either way
#'   via the returned `results` data frame.
#' @param flip_beta Logical; pass `TRUE` if your table's beta sign is the
#'   opposite convention (negative = hypermethylation).
#' @param device Device for AlphaGenome: `"cpu"` (default, safe
#'   everywhere), `"gpu"`, or `"cuda"`.
#' @param cell_types Character vector of AlphaGenome cell-type labels to
#'   score against; `NULL` (default) uses the script's built-in default
#'   panel (monocyte, CD4+/CD8+ T cell, B cell, neutrophil, PBMC).
#' @param repo_dir Chorus repository checkout; defaults to
#'   [chorus_repo_dir()].
#' @param mamba_env Conda/mamba environment to run the script in; default
#'   `"chorus"` (the base environment -- AlphaGenome itself runs in its
#'   own environment as a subprocess via `use_environment=True`).
#' @param mamba_bin Name or path of the mamba/conda executable.
#'
#' @return A list with `status` (integer exit code from the script),
#'   `log` (character vector of the script's combined stdout/stderr, for
#'   troubleshooting a nonzero `status`), `output` (path to the results
#'   CSV), and `results` (that CSV read in as a data frame, or `NULL` if
#'   the script did not exit successfully).
#' @export
#' @examples
#' # Requires the CHORUS_REPO_DIR environment variable to point at a
#' # chorus repository checkout, and mamba/conda with the
#' # chorus-alphagenome environment set up (see chorus_repo_dir()).
#' res <- run_alphagenome_mqtl(chorus_example_mqtls_path())
#' head(res$results)
run_alphagenome_mqtl <- function(mqtl_csv,
                                  output = tempfile(fileext = ".csv"),
                                  flip_beta = FALSE,
                                  device = "cpu",
                                  cell_types = NULL,
                                  repo_dir = chorus_repo_dir(),
                                  mamba_env = "chorus",
                                  mamba_bin = "mamba") {
  if (!file.exists(mqtl_csv)) {
    stop("mqtl_csv '", mqtl_csv, "' does not exist.", call. = FALSE)
  }

  script <- .find_chorus_script("score_mqtls.py", repo_dir)

  args <- c(normalizePath(mqtl_csv, mustWork = TRUE), "-o", output,
            "--device", device)
  if (isTRUE(flip_beta)) args <- c(args, "--flip-beta")
  if (!is.null(cell_types)) args <- c(args, "--cell-types", cell_types)

  run <- .run_chorus_script(script, args, mamba_env = mamba_env,
                             mamba_bin = mamba_bin)

  results <- NULL
  if (identical(run$status, 0L) && file.exists(output)) {
    results <- utils::read.csv(output, stringsAsFactors = FALSE)
  }

  list(
    status = run$status,
    log = run$stdout,
    output = output,
    results = results
  )
}
