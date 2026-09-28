#' Path to a small bundled example mQTL table
#'
#' Four real, replicated whole-blood cis-mQTLs, for exercising
#' [run_alphagenome_mqtl()] and [run_cherimoya_mqtl()] without needing your
#' own data first.
#'
#' **Provenance.** All four SNP-CpG pairs and their effect sizes (`beta`,
#' the `beta_a1` column) come from the GoDMC whole-blood mQTL meta-analysis
#' (<http://www.godmc.org.uk>, `assoc_meta_all.csv.gz`), filtered to rows
#' with `cistrans == TRUE` and `clumped == TRUE` (i.e. a replicated,
#' index/lead cis-mQTL, not merely a candidate). GoDMC reports SNP
#' positions as `hg19`/GRCh37 and does not include the CpG probe's own
#' coordinate, so:
#'
#' - each CpG's genomic position was looked up from the UCSC
#'   `snpArrayIllumina450k` track (`api.genome.ucsc.edu/search`);
#' - both the SNP and the CpG position were lifted from hg19 to hg38 with
#'   `pyliftover` (UCSC `hg19ToHg38` chain), since chorus's oracles are
#'   hg38-only;
#' - `ref`/`alt` were assigned by querying the hg38 reference base at the
#'   lifted SNP position (UCSC `api.genome.ucsc.edu/getData/sequence`) and
#'   matching it against GoDMC's `allele1`/`allele2` -- in all four rows
#'   the hg38 reference base matched GoDMC's non-effect allele, so
#'   `beta` (GoDMC's `beta_a1`, the effect-allele coefficient) already
#'   follows this package's convention: positive = alt allele associated
#'   with increased methylation.
#'
#' The four rows deliberately span a range of SNP-CpG distances (326 bp,
#' 800 bp, 2,180 bp, and 37,326 bp): the first two fit inside Cherimoya's
#' 2,114 bp window and the last two don't, so running this table through
#' [run_cherimoya_mqtl()] exercises its window pre-filter (with a
#' `WARNING:` for the two dropped rows) as well as its scoring path,
#' while all four run through [run_alphagenome_mqtl()] (1,048,576 bp
#' window) unfiltered.
#'
#' This is a demonstration table, not a curated set of the "best" or most
#' biologically important mQTLs -- it's four rows from partway through a
#' 5.9 GB file, chosen only for being real, replicated, and covering both
#' short and long SNP-CpG distances.
#'
#' @return A file path to a bundled CSV, in the format documented for
#'   [run_alphagenome_mqtl()] / [run_cherimoya_mqtl()].
#' @export
#' @examples
#' chorus_example_mqtls_path()
#' read.csv(chorus_example_mqtls_path())
chorus_example_mqtls_path <- function() {
  system.file("extdata", "example_mqtls.csv", package = "chorusR",
              mustWork = TRUE)
}
