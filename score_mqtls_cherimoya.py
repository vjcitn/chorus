"""Score whole-blood mQTL SNPs with Cherimoya (CATv1) for chromatin-based support.

Cherimoya-specific variant of score_mqtls.py (AlphaGenome). Same input CSV
format and column conventions -- see alphagenome4mqtl.md -- but adapted for
two structural differences between the oracles:

1. Track resolution. Cherimoya has no JSON track catalog; tracks are ENCODE
   experiments (assay + biosample -> `ASSAY:ENCSR...`), resolved through
   CATv1Metadata().resolve(). Coverage is DNase/ATAC only -- no CAGE, no
   histone ChIP -- so this script only checks the "active" (accessibility)
   side of the chromatin/methylation concordance that score_mqtls.py checks
   for both active and repressive marks.

2. Window size. Cherimoya's input is 2,114 bp, vs AlphaGenome's 1,048,576 bp.
   The predicted region is centered on the SNP, so a CpG more than ~1,000 bp
   from its SNP falls outside the window entirely. Rows like that are
   dropped up front (reported, not silently skipped).

3. One-model-per-instance. Cherimoya's load_pretrained_model() loads exactly
   one ENCODE experiment at a time (like ChromBPNet), unlike AlphaGenome
   where one oracle instance serves all 5,168 tracks. So the loop here is
   inverted relative to score_mqtls.py: outer loop over resolved
   tracks/biosamples (load once each), inner loop over the qualifying
   mQTLs -- instead of AlphaGenome's outer-loop-over-mQTLs shape.

chorus does not predict methylation itself -- this script uses chromatin
accessibility as a mechanistic, orthogonal check on whether an mQTL SNP
plausibly disrupts a regulatory element, not as a replacement for the
methylation measurement.
"""
import argparse
import os
import sys
import traceback

import numpy as np
import pandas as pd

import chorus
from chorus.utils import get_genome
from chorus.core.result import score_variant_effect
from chorus.analysis.normalization import get_pertrack_normalizer
from chorus.oracles.cherimoya_source.catv1_globals import CATV1_INPUT_LENGTH
from chorus.oracles.cherimoya_source.catv1_metadata import get_metadata

from score_mqtls import load_mqtl_table, load_existing_keys, append_checkpoint

# Biosamples confirmed present in CATv1-metadata.tsv. Cherimoya's atlas has
# no neutrophil, PBMC, or generic "blood" biosample term, unlike
# AlphaGenome -- those are simply not covered here.
DEFAULT_BIOSAMPLES = [
    "CD14-positive monocyte",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "B cell",
]
DEFAULT_ASSAYS = ["DNASE", "ATAC"]

# Half the input window, minus a margin, as the max |cpg_pos - snp_pos| a
# variant-centered 2,114 bp region can still cover.
MAX_SNP_CPG_DISTANCE = CATV1_INPUT_LENGTH // 2 - 50


def resolve_tracks(biosamples, assays):
    """Return a list of (assay, encode_id, biosample) for every
    (assay, biosample) pair that resolves, skipping pairs CATv1 doesn't
    cover.
    """
    meta = get_metadata()
    resolved = []
    for biosample in biosamples:
        for assay in assays:
            try:
                resolved_assay, encode_id = meta.resolve(
                    assay=assay, cell_type=biosample
                )
                resolved.append((resolved_assay, encode_id, biosample))
            except KeyError as exc:
                print(
                    f"  skipping {assay} / {biosample!r}: {exc}",
                    file=sys.stderr,
                )
    return resolved


def filter_by_window(mqtls):
    dist = (mqtls["cpg_pos"] - mqtls["snp_pos"]).abs()
    too_far = mqtls[dist > MAX_SNP_CPG_DISTANCE]
    if len(too_far):
        print(
            f"WARNING: dropping {len(too_far)} row(s) where the CpG is "
            f"more than {MAX_SNP_CPG_DISTANCE} bp from the SNP -- outside "
            f"Cherimoya's {CATV1_INPUT_LENGTH} bp window even centered on "
            f"the variant: {too_far['mqtl_id'].tolist()}",
            file=sys.stderr,
        )
    return mqtls[dist <= MAX_SNP_CPG_DISTANCE].reset_index(drop=True)


def score_one_mqtl(oracle, row, norm, window_bins=2, bedgraph_dir=None):
    chrom = row["chrom"]
    snp_pos = int(row["snp_pos"])
    cpg_pos = int(row["cpg_pos"])
    half = CATV1_INPUT_LENGTH // 2
    region = f"{chrom}:{snp_pos - half}-{snp_pos + half}"
    variant_pos = f"{chrom}:{snp_pos}"
    track_id = oracle.track_id

    variant_effects = oracle.predict_variant_effect(
        region, variant_pos, [row["ref"], row["alt"]], [track_id]
    )

    if bedgraph_dir is not None:
        # Full-resolution ref and alt predicted signal over the window, for
        # inspection in IGV -- separate from the single scalar effect score
        # computed below. track_id is included since the outer loop calls
        # this once per (mqtl, track) pair.
        clean_track = track_id.replace(":", "_")
        prefix = f"{row['mqtl_id']}_{clean_track}"
        variant_effects["predictions"]["reference"].save_predictions_as_bedgraph(
            output_dir=bedgraph_dir, prefix=f"{prefix}_ref"
        )
        variant_effects["predictions"]["alt_1"].save_predictions_as_bedgraph(
            output_dir=bedgraph_dir, prefix=f"{prefix}_alt"
        )

    if abs(cpg_pos - snp_pos) <= 5:
        scores = score_variant_effect(
            variant_effects, at_variant=True, window_bins=window_bins
        )
    else:
        scores = score_variant_effect(
            variant_effects,
            chrom=chrom,
            start=cpg_pos - 1,
            end=cpg_pos + 1,
            scoring_strategy="mean",
        )

    # score_variant_effect() keys its return dict by "alt_1", "alt_2", ...
    # (chorus/core/base.py: `intervals[f'alt_{i+1}'] = ...`), never by the
    # literal allele letter. Exactly one alt allele is ever passed above, so
    # the key is always "alt_1".
    effect = scores["alt_1"][track_id]["effect"]
    eff_pct = norm.effect_percentile("cherimoya", track_id, abs(effect), signed=False)
    # Accessibility is an "active" mark: more signal -> less methylation.
    concordant = np.sign(effect) == np.sign(-row["beta"])

    return {
        "mqtl_id": row["mqtl_id"],
        "chrom": chrom,
        "snp_pos": snp_pos,
        "cpg_pos": cpg_pos,
        "ref": row["ref"],
        "alt": row["alt"],
        "beta": row["beta"],
        "track_id": track_id,
        "log2fc_effect": effect,
        "effect_percentile": eff_pct,
        "concordant_with_methylation": bool(concordant),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mqtl_csv", help="CSV of mQTL SNP-CpG pairs")
    ap.add_argument("-o", "--output", default="mqtl_cherimoya_scores.csv")
    ap.add_argument("--flip-beta", action="store_true",
                     help="Pass if your table's beta sign is the opposite "
                          "convention (negative = hypermethylation)")
    ap.add_argument("--device", default=None,
                     help="Cherimoya has no bf16 dependency, so a T4 (or "
                          "any CUDA GPU) is fine; default lets chorus "
                          "auto-detect.")
    ap.add_argument("--biosamples", nargs="+", default=DEFAULT_BIOSAMPLES)
    ap.add_argument("--assays", nargs="+", default=DEFAULT_ASSAYS,
                     choices=["DNASE", "ATAC"])
    ap.add_argument("--save-bedgraph", metavar="DIR", default=None,
                     help="Also save each mQTL's full-resolution ref/alt "
                          "predicted signal (over the Cherimoya window) as "
                          "BedGraph files under this directory, one pair "
                          "per (mqtl, track). Off by default -- this is a "
                          "lot of files for a large mQTL table.")
    ap.add_argument("--resume", action="store_true",
                     help="If --output already has rows (e.g. from a run "
                          "that was interrupted), skip (mqtl_id, track_id) "
                          "pairs already scored there and append only the "
                          "rest. Without this flag, --output is "
                          "overwritten from scratch, as before.")
    args = ap.parse_args()

    mqtls = load_mqtl_table(args.mqtl_csv, flip_beta=args.flip_beta)
    print(f"Loaded {len(mqtls)} SNP-CpG pairs", flush=True)
    mqtls = filter_by_window(mqtls)
    print(f"{len(mqtls)} pair(s) fit inside Cherimoya's "
          f"{CATV1_INPUT_LENGTH} bp window", flush=True)

    tracks = resolve_tracks(args.biosamples, args.assays)
    print(f"Resolved {len(tracks)} track(s) across "
          f"{len(args.biosamples)} biosample(s): "
          f"{[(a, e) for a, e, _ in tracks]}", flush=True)

    norm = get_pertrack_normalizer("cherimoya")
    genome_path = get_genome("hg38")

    if args.resume:
        done_keys = load_existing_keys(args.output, ["mqtl_id", "track_id"])
        header_written = bool(done_keys)
        if done_keys:
            print(f"--resume: {len(done_keys)} (mqtl_id, track_id) pair(s) "
                  f"already in {args.output}, skipping those", flush=True)
    else:
        done_keys = set()
        header_written = False
        if os.path.exists(args.output):
            os.remove(args.output)

    # Each (mqtl, track) row is appended to --output as soon as it's
    # scored (see append_checkpoint), so a crash or kill part-way through
    # only costs the row in progress, not the whole table -- rerun with
    # --resume to pick up where it left off.
    for assay, encode_id, biosample in tracks:
        print(f"\n=== {assay} {encode_id} ({biosample}) ===", flush=True)
        oracle = chorus.create_oracle(
            "cherimoya",
            use_environment=True,
            reference_fasta=str(genome_path),
            device=args.device,
        )
        oracle.load_pretrained_model(assay=assay, encode_id=encode_id)
        track_id = oracle.track_id

        for i, row in mqtls.iterrows():
            if (row["mqtl_id"], track_id) in done_keys:
                continue
            print(f"  [{i + 1}/{len(mqtls)}] {row['mqtl_id']} "
                  f"{row['chrom']}:{row['snp_pos']} {row['ref']}>{row['alt']}",
                  flush=True)
            try:
                result_row = score_one_mqtl(
                    oracle, row, norm, bedgraph_dir=args.save_bedgraph
                )
                result_row["biosample"] = biosample
                header_written = append_checkpoint(
                    args.output, [result_row], header_written
                )
            except Exception as exc:
                # A bare str(exc) reads as a data problem for a KeyError
                # (e.g. "FAILED: 'A'") and gives no way to tell which line
                # raised it. The full traceback costs nothing.
                print(f"    FAILED: {exc}", file=sys.stderr, flush=True)
                traceback.print_exc(file=sys.stderr)

    result = (
        pd.read_csv(args.output) if os.path.exists(args.output) else pd.DataFrame()
    )
    print(f"\nWrote {len(result)} rows to {args.output}", flush=True)

    if len(result):
        summary = (
            result.groupby("mqtl_id")
            .agg(
                max_effect_percentile=("effect_percentile", "max"),
                frac_concordant=("concordant_with_methylation", "mean"),
            )
            .sort_values("max_effect_percentile", ascending=False)
        )
        print("\nTop candidates by max effect percentile across tracks:")
        print(summary.head(20).to_string())


if __name__ == "__main__":
    main()
