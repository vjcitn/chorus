"""Score whole-blood mQTL SNPs with AlphaGenome for chromatin-based support.

Input: a CSV with one row per SNP-CpG pair. Required columns:
    chrom, snp_pos, ref, alt, cpg_pos, beta
Optional columns:
    mqtl_id, gene_symbol

See the accompanying README note / conversation for the exact column
conventions (1-based positions, hg38, beta sign = alt-allele direction of
methylation change).

Output: one long-format CSV, one row per (mqtl, track), with the predicted
effect, its background percentile, and whether the chromatin direction is
concordant with the reported methylation direction.

chorus does not predict methylation itself -- this script uses chromatin
accessibility and histone marks as a mechanistic, orthogonal check on
whether an mQTL SNP plausibly disrupts a regulatory element, not as a
replacement for the methylation measurement.
"""
import argparse
import sys

import numpy as np
import pandas as pd

import chorus
from chorus.utils import get_genome
from chorus.core.result import score_variant_effect
from chorus.analysis.normalization import get_pertrack_normalizer

# Marks where MORE signal is expected to go with LESS methylation
# (open chromatin / active enhancer or promoter marks).
ACTIVE_ASSAY_TYPES = {"DNASE", "ATAC", "CAGE"}
ACTIVE_HISTONE_KEYWORDS = ("H3K4me1", "H3K4me3", "H3K27ac", "H3K9ac")

# Marks where MORE signal is expected to go with MORE methylation
# (repressive/heterochromatic marks).
REPRESSIVE_HISTONE_KEYWORDS = ("H3K27me3", "H3K9me3")

DEFAULT_CELL_TYPES = [
    "CD14-positive monocyte",
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "B cell",
    "neutrophil",
    "peripheral blood mononuclear cell",
]


def select_tracks(oracle, cell_types):
    """Return {assay_id: expected_direction} for DNase/ATAC/CAGE/histone
    tracks matching the given cell types. expected_direction is +1
    (more signal -> less methylation) or -1 (more signal -> more
    methylation).
    """
    wanted = {}
    for assay_type in ["DNASE", "ATAC", "CAGE", "CHIP_HISTONE"]:
        info = oracle.get_track_info(assay_type)
        if info is None or len(info) == 0:
            continue
        hit = info[info["cell_type"].isin(cell_types)]
        for _, row in hit.iterrows():
            assay_id = row["identifier"]
            name = row.get("name", "")
            if assay_type in ACTIVE_ASSAY_TYPES:
                wanted[assay_id] = +1
            elif any(k in name for k in ACTIVE_HISTONE_KEYWORDS):
                wanted[assay_id] = +1
            elif any(k in name for k in REPRESSIVE_HISTONE_KEYWORDS):
                wanted[assay_id] = -1
            # else: a histone mark we don't have a directional prior for
            # (e.g. H3K36me3) -- skip rather than guess.
    return wanted


def load_mqtl_table(path, flip_beta=False):
    df = pd.read_csv(path)
    required = ["chrom", "snp_pos", "ref", "alt", "cpg_pos", "beta"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"mQTL table is missing required columns: {missing}")

    df = df.copy()
    if "mqtl_id" not in df.columns:
        df["mqtl_id"] = [f"mqtl_{i}" for i in range(len(df))]
    if "gene_symbol" not in df.columns:
        df["gene_symbol"] = None

    df["chrom"] = df["chrom"].astype(str).apply(
        lambda c: c if c.startswith("chr") else f"chr{c}"
    )
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()

    bad_alleles = df[(df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)]
    if len(bad_alleles):
        print(
            f"WARNING: dropping {len(bad_alleles)} row(s) with non-SNV "
            f"ref/alt (this script scores single-base variants only): "
            f"{bad_alleles['mqtl_id'].tolist()}",
            file=sys.stderr,
        )
        df = df.drop(bad_alleles.index)

    if flip_beta:
        df["beta"] = -df["beta"]

    return df.reset_index(drop=True)


def score_one_mqtl(oracle, row, assay_ids, norm, window_bins=2):
    chrom = row["chrom"]
    snp_pos = int(row["snp_pos"])
    cpg_pos = int(row["cpg_pos"])
    region = f"{chrom}:{snp_pos - 500}-{snp_pos + 500}"  # auto-widened to full window, centered on the SNP
    variant_pos = f"{chrom}:{snp_pos}"

    variant_effects = oracle.predict_variant_effect(
        region, variant_pos, [row["ref"], row["alt"]], assay_ids
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

    out_rows = []
    alt_allele_name = row["alt"]
    for assay_id, per_track in scores[alt_allele_name].items():
        effect = per_track["effect"]
        expected_dir = assay_ids[assay_id]  # +1 or -1
        eff_pct = norm.effect_percentile(
            "alphagenome", assay_id, abs(effect), signed=False
        )
        # Concordant: predicted chromatin change direction matches what
        # the mQTL's beta implies, given whether this mark is active or
        # repressive.
        # expected_dir=+1 (active mark): effect and -beta should agree in sign.
        # expected_dir=-1 (repressive mark): effect and beta should agree in sign.
        concordant = np.sign(effect) == np.sign(-row["beta"] * expected_dir)

        out_rows.append(
            {
                "mqtl_id": row["mqtl_id"],
                "chrom": chrom,
                "snp_pos": snp_pos,
                "cpg_pos": cpg_pos,
                "ref": row["ref"],
                "alt": row["alt"],
                "beta": row["beta"],
                "assay_id": assay_id,
                "expected_direction": "active" if expected_dir > 0 else "repressive",
                "log2fc_effect": effect,
                "effect_percentile": eff_pct,
                "concordant_with_methylation": bool(concordant),
            }
        )
    return out_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mqtl_csv", help="CSV of mQTL SNP-CpG pairs")
    ap.add_argument("-o", "--output", default="mqtl_chromatin_scores.csv")
    ap.add_argument("--flip-beta", action="store_true",
                     help="Pass if your table's beta sign is the opposite "
                          "convention (negative = hypermethylation)")
    ap.add_argument("--device", default="cpu",
                     help="'cpu', 'gpu', or 'cuda'; default cpu (safe on "
                          "pre-Ampere GPUs, e.g. Tesla T4)")
    ap.add_argument("--cell-types", nargs="+", default=DEFAULT_CELL_TYPES)
    args = ap.parse_args()

    mqtls = load_mqtl_table(args.mqtl_csv, flip_beta=args.flip_beta)
    print(f"Loaded {len(mqtls)} SNP-CpG pairs")

    genome_path = get_genome("hg38")
    oracle = chorus.create_oracle(
        "alphagenome",
        use_environment=True,
        reference_fasta=str(genome_path),
        device=args.device,
    )
    oracle.load_pretrained_model()

    assay_ids = select_tracks(oracle, args.cell_types)
    print(f"Scoring against {len(assay_ids)} tracks across "
          f"{len(args.cell_types)} cell types: {sorted(assay_ids)}")

    norm = get_pertrack_normalizer("alphagenome")

    all_rows = []
    # predict_session() keeps one subprocess + loaded model alive across
    # the whole loop instead of reloading AlphaGenome (and refetching its
    # ~330 MB GCS reference tables) once per mQTL.
    with oracle.predict_session():
        for i, row in mqtls.iterrows():
            print(f"[{i + 1}/{len(mqtls)}] {row['mqtl_id']} "
                  f"{row['chrom']}:{row['snp_pos']} {row['ref']}>{row['alt']}")
            try:
                all_rows.extend(
                    score_one_mqtl(oracle, row, assay_ids, norm)
                )
            except Exception as exc:
                print(f"  FAILED: {exc}", file=sys.stderr)

    result = pd.DataFrame(all_rows)
    result.to_csv(args.output, index=False)
    print(f"\nWrote {len(result)} rows to {args.output}")

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
