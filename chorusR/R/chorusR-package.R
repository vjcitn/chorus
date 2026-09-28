#' chorusR: R Wrappers Around the chorus Genomic Oracle Toolkit
#'
#' @description
#' `chorusR` drives the Python `chorus` package's per-oracle scripts via
#' [base::system2()], rather than an in-process bridge such as
#' 'reticulate'. chorus already isolates each oracle's ML stack in its own
#' conda/mamba environment (JAX for AlphaGenome, PyTorch+Triton for
#' Cherimoya, etc.); shelling out to a script that itself calls
#' `mamba run -n <oracle-env> ...` internally is a natural fit for that
#' design, and keeps the R session decoupled from any one Python
#' environment.
#'
#' Currently two oracles are wrapped, both via mQTL chromatin-concordance
#' scoring scripts shipped in the chorus repository:
#'
#' - [run_alphagenome_mqtl()] -- AlphaGenome (DNase/ATAC/CAGE/histone
#'   ChIP, 1,048,576 bp window).
#' - [run_cherimoya_mqtl()] -- Cherimoya / CATv1 (DNase/ATAC only,
#'   2,114 bp window).
#'
#' Both expect a chorus repository checkout on disk; point at it with
#' [chorus_repo_dir()] (via `options(chorusR.repo_dir = ...)` or the
#' `CHORUS_REPO_DIR` environment variable) before calling either wrapper.
#'
#' @keywords internal
"_PACKAGE"
