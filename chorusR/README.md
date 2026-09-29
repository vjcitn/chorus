# chorusR

Thin R wrappers around the Python [`chorus`](https://github.com/pinellolab/chorus)
genomic oracle toolkit. Rather than bridging into Python in-process (e.g. via
`reticulate`), chorusR shells out with `system2()` to scripts that run inside
chorus's own per-oracle conda/mamba environments — the same way chorus runs
each oracle's forward pass in isolation already. This avoids tying the R
session's process to a specific Python environment, and each oracle keeps its
own (sometimes conflicting) dependency set.

Currently wraps two oracles via the mQTL chromatin-concordance scoring
scripts shipped in the chorus repository root:

- **AlphaGenome** — `run_alphagenome_mqtl()`, wraps `score_mqtls.py`
- **Cherimoya (CATv1)** — `run_cherimoya_mqtl()`, wraps `score_mqtls_cherimoya.py`

Both take a CSV of SNP-CpG pairs, predict each variant's effect on chromatin
tracks in matched cell types, and score whether the predicted direction of
change is concordant with the mQTL's reported methylation direction. See
`alphagenome4mqtl.md` in the chorus repository for the full column spec and
scoring methodology.

## Requirements

- A checkout of the [chorus](https://github.com/pinellolab/chorus) repository
  (point `CHORUS_REPO_DIR` at it, or pass `repo_dir=` explicitly).
- `mamba`/`conda` with the relevant per-oracle environments set up (see the
  chorus repository's own README/CLAUDE.md) — at minimum `chorus` (base) and
  the oracle-specific env (`chorus-alphagenome`, `chorus-cherimoya`, ...).
  chorusR calls scripts via each environment's `python` binary directly
  (resolved from `mamba env list`), not via `mamba run`/`conda run`.

## Basic usage

```r
library(chorusR)

res <- run_cherimoya_mqtl(chorus_example_mqtls_path())
head(res$results)
```

Each wrapper returns a list with `status` (script exit code), `log` /
`log_file` (combined stdout/stderr, written live so you can `tail -f` it from
another terminal during a long run), `output` (path to the results CSV,
written incrementally one row at a time), and `results` (that CSV read back
in, or `NULL` if the run failed or produced nothing).

## Device / GPU notes

Both wrappers accept a `device` argument (`"cpu"`, `"cuda"`/`"gpu"`, or
`NULL`/omitted to let chorus auto-detect). A few hard-won gotchas, from
debugging real runs on a cloud GPU host:

- **AlphaGenome needs Ampere or newer.** Its JAX backend requires bf16
  tensor-core matmul, which pre-Ampere GPUs (Tesla T4, V100, P100, RTX
  20-series) don't support — auto-detect picks the GPU and fails with
  `UNIMPLEMENTED: ... ALG_DOT_BF16_BF16_F32`. Pass `device = "cpu"` on such
  hardware (this is `run_alphagenome_mqtl()`'s default for that reason).
- **Virtualized/cloud GPUs can fail a different way.** Some GPU-passthrough
  virtualization layers don't forward CUDA's Virtual Memory Management (VMM)
  API to the guest even though the underlying hardware supports it. JAX then
  fails to initialize the CUDA backend at all (`Device 0 does not support
  CUDA Virtual Memory Management (VMM)`), and PyTorch's
  `torch.cuda.is_available()` reports `False` on the same host — so this
  affects Cherimoya too, and isn't fixable by picking a different jaxlib
  version or config flag. `device = "cpu"` is the only option on such a host,
  for both oracles.
- **Cherimoya itself has no bf16 dependency**, so on real (non-virtualized)
  hardware a T4 is fine for it even where AlphaGenome isn't — the CPU-only
  fallback above is specifically about GPU-passthrough hosts, not a general
  Cherimoya limitation.
- There is no `CHORUS_DEVICE` environment variable set by chorusR or by
  chorus's own environments by default; if a run logs a device you didn't
  ask for, check whether something in the shell (e.g. left over from an
  earlier experiment) is exporting it, or whether `device=` was passed
  explicitly somewhere upstream in the call.

## Known past issues (fixed)

These bit real runs during development and are recorded here so a similar
symptom doesn't get re-diagnosed from scratch:

- **Cherimoya: "No reference genome provided"` on every row.** Unlike
  AlphaGenome, whose oracle fetches sequence itself, Cherimoya needs a local
  reference FASTA to pull raw sequence around each variant.
  `score_mqtls_cherimoya.py` now resolves this once via
  `chorus.utils.get_genome("hg38")` and passes `reference_fasta=` to
  `create_oracle()`, mirroring how `score_mqtls.py` already did it for
  AlphaGenome.
- **Every row scored `FAILED: '<allele letter>'`.** `score_variant_effect()`
  keys its return dict by `"alt_1"`, `"alt_2"`, ... — never by the literal
  allele base. Both scoring scripts were looking up `scores[row["alt"]]`
  (e.g. `scores["A"]`), which is always a `KeyError`, silently caught by the
  per-row try/except so the script still exited 0 having scored nothing.
  Fixed to key by `"alt_1"` (always correct, since exactly one alt allele is
  ever scored per call). This was masked for a long time because earlier
  runs died even earlier (the reference-genome and GPU issues above) before
  ever reaching this line.
- **`mamba run -n <env> ...` failing with `exec: --: invalid option`.** A
  `mamba run` wrapper-script quirk (bash's `exec` builtin doesn't treat `--`
  as an end-of-options marker), not fixable by reordering flags. Resolved by
  having chorusR resolve and invoke each environment's `python` binary
  directly instead of going through `mamba run`/`conda run` at all.

## Development

```r
devtools::test()      # run the test suite
devtools::document()  # regenerate man/ and NAMESPACE from roxygen comments
```

Bump `Version:` in `DESCRIPTION` with every commit that touches chorusR code,
tests, or docs (this does not apply to commits that touch only the upstream
`chorus` Python scripts without touching chorusR itself).
