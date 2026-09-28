#' Score whole-blood mQTLs against Cherimoya (CATv1) chromatin tracks
#'
#' Wraps `score_mqtls_cherimoya.py` (in the chorus repository root): for
#' each SNP-CpG pair in `mqtl_csv`, predicts the variant's effect on
#' Cherimoya DNase/ATAC accessibility tracks in matched blood cell types,
#' and reports whether the predicted accessibility change is concordant
#' with the mQTL's reported methylation direction. Uses the same
#' `mqtl_csv` format as [run_alphagenome_mqtl()] -- see
#' `alphagenome4mqtl.md` in the chorus repository for the column spec.
#'
#' Two structural differences from AlphaGenome shape this wrapper:
#' Cherimoya's atlas covers only DNase/ATAC (no CAGE or histone marks, so
#' only the accessibility side of the chromatin/methylation concordance
#' check is available), and its input window is 2,114 bp rather than
#' AlphaGenome's 1,048,576 bp, so SNP-CpG pairs more than ~1,000 bp apart
#' are dropped before scoring (reported, not silently skipped) -- see the
#' script's own `WARNING:` lines in `log` if that happens to your table.
#' Unlike AlphaGenome, Cherimoya has no bf16 tensor-core dependency, so
#' `device` defaults to `NULL` (chorus auto-detects; a Tesla T4 is fine).
#'
#' @param mqtl_csv Path to a CSV of SNP-CpG pairs. Required columns:
#'   `chrom`, `snp_pos`, `ref`, `alt`, `cpg_pos`, `beta`. Optional:
#'   `mqtl_id`, `gene_symbol`.
#' @param output Path to write the long-format results CSV (one row per
#'   mqtl x track). Defaults to a temp file; read back into R either way
#'   via the returned `results` data frame.
#' @param flip_beta Logical; pass `TRUE` if your table's beta sign is the
#'   opposite convention (negative = hypermethylation).
#' @param device Device for Cherimoya; `NULL` (default) lets chorus
#'   auto-detect (safe on a Tesla T4, unlike AlphaGenome).
#' @param biosamples Character vector of CATv1 biosample term names to
#'   score against; `NULL` (default) uses the script's built-in default
#'   panel (CD14-positive monocyte, CD4+/CD8+ alpha-beta T cell, B cell).
#'   Note CATv1 has no neutrophil, PBMC, or generic "blood" biosample,
#'   unlike AlphaGenome.
#' @param assays Character vector, subset of `c("DNASE", "ATAC")`; `NULL`
#'   (default) uses both.
#' @param save_bedgraph Optional directory; if given, also saves each
#'   mqtl's full-resolution ref/alt predicted signal as BedGraph files
#'   there (one pair per mqtl x track) via the script's
#'   `--save-bedgraph` option. Off by default.
#' @param resume Logical; if `TRUE` and `output` already has rows (e.g.
#'   from a run that was interrupted), skip (mqtl_id, track_id) pairs
#'   already scored there and append only the rest, instead of
#'   overwriting `output` from scratch. Only useful with a stable
#'   (non-default, non-tempfile) `output` path, since a resumed run needs
#'   to find the same file the earlier run wrote to.
#' @param repo_dir Chorus repository checkout; defaults to
#'   [chorus_repo_dir()].
#' @param mamba_env Conda/mamba environment to run the script in; default
#'   `"chorus"` (the base environment -- Cherimoya itself runs in its own
#'   environment as a subprocess via `use_environment=True`).
#' @param mamba_bin Name or path of the mamba/conda executable.
#'
#' @return A list with `status` (integer exit code from the script),
#'   `log` (character vector of the script's combined stdout/stderr,
#'   including any window-filter or track-resolution warnings), `log_file`
#'   (path to that same log on disk, written live as the script runs --
#'   Cherimoya loads one track/biosample at a time and scores every mQTL
#'   against it before moving to the next, so `tail -f log_file` from
#'   another terminal shows progress in real time, and the file survives
#'   even if this R session is interrupted), `output` (path to the
#'   results CSV, itself written incrementally one (mqtl, track) pair at
#'   a time rather than only at the end, so a crash or interrupt only
#'   loses the row in progress -- rerun with `resume = TRUE` to pick up
#'   where it left off), and `results` (that CSV read in as a data frame,
#'   or `NULL` if the script did not exit successfully).
#' @export
#' @examples
#' # Requires the CHORUS_REPO_DIR environment variable to point at a
#' # chorus repository checkout, and mamba/conda with the
#' # chorus-cherimoya environment set up (see chorus_repo_dir()). Two of
#' # the four bundled example rows are farther apart than Cherimoya's
#' # window and will be dropped with a WARNING: see
#' # chorus_example_mqtls_path().
#' res <- run_cherimoya_mqtl(chorus_example_mqtls_path())
#' head(res$results)
run_cherimoya_mqtl <- function(mqtl_csv,
                                output = tempfile(fileext = ".csv"),
                                flip_beta = FALSE,
                                device = NULL,
                                biosamples = NULL,
                                assays = NULL,
                                save_bedgraph = NULL,
                                resume = FALSE,
                                repo_dir = chorus_repo_dir(),
                                mamba_env = "chorus",
                                mamba_bin = "mamba") {
  if (!file.exists(mqtl_csv)) {
    stop("mqtl_csv '", mqtl_csv, "' does not exist.", call. = FALSE)
  }
  if (!is.null(assays) && !all(assays %in% c("DNASE", "ATAC"))) {
    stop("assays must be a subset of c(\"DNASE\", \"ATAC\").", call. = FALSE)
  }

  script <- .find_chorus_script("score_mqtls_cherimoya.py", repo_dir)

  args <- c(normalizePath(mqtl_csv, mustWork = TRUE), "-o", output)
  if (isTRUE(flip_beta)) args <- c(args, "--flip-beta")
  if (!is.null(device)) args <- c(args, "--device", device)
  if (!is.null(biosamples)) args <- c(args, "--biosamples", biosamples)
  if (!is.null(assays)) args <- c(args, "--assays", assays)
  if (!is.null(save_bedgraph)) args <- c(args, "--save-bedgraph", save_bedgraph)
  if (isTRUE(resume)) args <- c(args, "--resume")

  run <- .run_chorus_script(script, args, mamba_env = mamba_env,
                             mamba_bin = mamba_bin)

  results <- .read_mqtl_results(output, run)

  list(
    status = run$status,
    log = run$stdout,
    log_file = run$log_file,
    output = output,
    results = results
  )
}
