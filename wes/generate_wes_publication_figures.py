#!/usr/bin/env python3
"""Generate the WES comparison tables and manuscript panels used for this case.

The analysis is deliberately restricted to the operations represented in the
manuscript: exact PASS-variant membership, AD-based VAF calculation, comparison
of primary and recurrent tumors, and highlighting of exact-shared variants with
SnpEff MODERATE/HIGH predicted impact.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

VariantKey = Tuple[str, int, str, str]


@dataclass(frozen=True)
class Call:
    key: VariantKey
    ref_ad: int
    alt_ad: int
    ad_depth: int
    vaf: float
    format_af: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wes-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--no-strict", action="store_true")
    return parser.parse_args()


def natural_chrom_key(chrom: str):
    value = chrom.removeprefix("chr")
    if value.isdigit():
        return int(value), ""
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return special.get(value, 100), value


def variant_sort_key(key: VariantKey):
    rank, suffix = natural_chrom_key(key[0])
    return rank, suffix, key[1], key[2], key[3]


def read_pass_vcf(path: Path, tumor_sample: str) -> Dict[VariantKey, Call]:
    calls: Dict[VariantKey, Call] = {}
    sample_names = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                sample_names = line.rstrip().split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, pos, _id, ref, alt, _qual, filt, _info, fmt = fields[:9]
            if filt != "PASS":
                continue
            if "," in alt:
                raise ValueError(f"Multiallelic record not decomposed: {chrom}:{pos} {ref}>{alt}")
            idx = sample_names.index(tumor_sample)
            values = dict(zip(fmt.split(":"), fields[9 + idx].split(":")))
            ad = [int(x) for x in values["AD"].split(",")]
            if len(ad) != 2 or sum(ad) <= 0:
                raise ValueError(f"Invalid biallelic AD at {chrom}:{pos}")
            key = (chrom, int(pos), ref, alt)
            calls[key] = Call(key, ad[0], ad[1], sum(ad), ad[1] / sum(ad), values.get("AF", ""))
    return calls


def read_annotation_table(path: Path) -> Dict[VariantKey, Dict[str, str]]:
    out = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])
            out[key] = row
    return out


def select_annotation(key, primary_annotations, recurrent_annotations):
    return primary_annotations.get(key) or recurrent_annotations.get(key) or {}


def compute_results(wes_root: Path):
    variant_dir = wes_root / "results" / "variants"
    annotation_dir = wes_root / "results" / "annotated"

    primary = read_pass_vcf(
        variant_dir / "Primary_CIC_vs_Blood.filtered.PASS.vcf.gz", "Primary_CIC"
    )
    recurrent = read_pass_vcf(
        variant_dir / "Recurrent_CIC_vs_Blood.filtered.PASS.vcf.gz", "Recurrent_CIC"
    )
    primary_annotations = read_annotation_table(
        annotation_dir / "Primary_CIC_vs_Blood.annotated_variants.tsv"
    )
    recurrent_annotations = read_annotation_table(
        annotation_dir / "Recurrent_CIC_vs_Blood.annotated_variants.tsv"
    )

    pkeys, rkeys = set(primary), set(recurrent)
    shared = pkeys & rkeys
    primary_specific = pkeys - rkeys
    recurrent_specific = rkeys - pkeys

    highlighted = {
        key for key in shared
        if select_annotation(key, primary_annotations, recurrent_annotations).get("IMPACT")
        in {"MODERATE", "HIGH"}
    }

    ordered_shared = sorted(shared, key=variant_sort_key)
    x = np.array([primary[k].vaf for k in ordered_shared])
    y = np.array([recurrent[k].vaf for k in ordered_shared])
    pearson_r = float(np.corrcoef(x, y)[0, 1]) if len(ordered_shared) > 1 else np.nan

    summary = {
        "primary_pass": len(primary),
        "recurrent_pass": len(recurrent),
        "primary_specific": len(primary_specific),
        "exact_shared": len(shared),
        "recurrent_specific": len(recurrent_specific),
        "union": len(pkeys | rkeys),
        "shared_moderate_high": len(highlighted),
        "pearson_r": pearson_r,
        "primary_specific_keys": primary_specific,
        "shared_keys": shared,
        "recurrent_specific_keys": recurrent_specific,
        "highlighted_keys": highlighted,
    }
    return primary, recurrent, primary_annotations, recurrent_annotations, summary


def format_float(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}f}"


def validate_current_package(summary, primary_annotations, recurrent_annotations):
    expected = {
        "primary_pass": 505,
        "recurrent_pass": 453,
        "primary_specific": 391,
        "exact_shared": 114,
        "recurrent_specific": 339,
        "union": 844,
        "shared_moderate_high": 8,
    }
    observed = {k: int(summary[k]) for k in expected}
    if observed != expected:
        raise ValueError(f"Reviewed publication counts changed: expected {expected}; observed {observed}")


def write_variant_tables(output_dir, primary, recurrent, pa, ra, summary):
    output_dir.mkdir(parents=True, exist_ok=True)

    membership = output_dir / "CIC_WES_exact_variant_membership.tsv"
    with membership.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "MEMBERSHIP"])
        for key in sorted(summary["primary_specific_keys"], key=variant_sort_key):
            writer.writerow([*key, "Primary-specific"])
        for key in sorted(summary["shared_keys"], key=variant_sort_key):
            writer.writerow([*key, "Exact shared"])
        for key in sorted(summary["recurrent_specific_keys"], key=variant_sort_key):
            writer.writerow([*key, "Recurrent-specific"])

    shared_table = output_dir / "CIC_WES_exact_shared_variant_VAF.tsv"
    with shared_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "CHROM", "POS", "REF", "ALT", "GENE", "EFFECT", "IMPACT",
            "PRIMARY_REF_AD", "PRIMARY_ALT_AD", "PRIMARY_VAF_AD",
            "RECURRENT_REF_AD", "RECURRENT_ALT_AD", "RECURRENT_VAF_AD",
        ])
        for key in sorted(summary["shared_keys"], key=variant_sort_key):
            ann = select_annotation(key, pa, ra)
            writer.writerow([
                *key, ann.get("GENE", ""), ann.get("EFFECT", ""), ann.get("IMPACT", ""),
                primary[key].ref_ad, primary[key].alt_ad, format_float(primary[key].vaf),
                recurrent[key].ref_ad, recurrent[key].alt_ad, format_float(recurrent[key].vaf),
            ])

    with (output_dir / "CIC_WES_figure_statistics.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key in [
            "primary_pass", "recurrent_pass", "primary_specific", "exact_shared",
            "recurrent_specific", "union", "shared_moderate_high", "pearson_r"
        ]:
            writer.writerow([key, summary[key]])


def configure_matplotlib():
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def draw_membership_panel(summary, output_dir):
    counts = [summary["primary_specific"], summary["exact_shared"], summary["recurrent_specific"]]
    labels = ["Primary-specific", "Exact shared", "Recurrent-specific"]
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    bars = ax.bar(labels, counts)
    ax.set_ylabel("PASS variants")
    ax.set_title("Exact somatic variant membership")
    ax.tick_params(axis="x", rotation=20)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), str(count), ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(output_dir / "CIC_WES_panel_B_exact_variant_overview.png", dpi=600)
    fig.savefig(output_dir / "CIC_WES_panel_B_exact_variant_overview.pdf")
    plt.close(fig)


def draw_vaf_panel(primary, recurrent, pa, ra, summary, output_dir):
    shared = sorted(summary["shared_keys"], key=variant_sort_key)
    highlighted = set(summary["highlighted_keys"])
    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot([0, 1], [0, 1], "--", color="0.5", linewidth=1)
    other = [k for k in shared if k not in highlighted]
    ax.scatter([primary[k].vaf for k in other], [recurrent[k].vaf for k in other],
               s=24, alpha=0.6, label=f"Other exact shared (n={len(other)})")
    for impact, marker in [("MODERATE", "o"), ("HIGH", "D")]:
        keys = [k for k in highlighted if select_annotation(k, pa, ra).get("IMPACT") == impact]
        ax.scatter([primary[k].vaf for k in keys], [recurrent[k].vaf for k in keys],
                   s=55, marker=marker, label=f"{impact.title()} impact (n={len(keys)})")
        for k in keys:
            gene = select_annotation(k, pa, ra).get("GENE", "")
            ax.annotate(gene, (primary[k].vaf, recurrent[k].vaf), xytext=(5, 4),
                        textcoords="offset points", fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Primary tumor VAF (AD-based)")
    ax.set_ylabel("Recurrent tumor VAF (AD-based)")
    ax.set_title("Primary vs recurrent VAF of exact shared PASS variants")
    ax.text(0.04, 0.96, f"n = {summary['exact_shared']}\nPearson r = {summary['pearson_r']:.2f}",
            transform=ax.transAxes, va="top")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_dir / "CIC_WES_panel_C_exact_shared_VAF_scatter.png", dpi=600)
    fig.savefig(output_dir / "CIC_WES_panel_C_exact_shared_VAF_scatter.pdf")
    plt.close(fig)


def write_checksums(output_dir: Path) -> None:
    import hashlib
    path = output_dir / "SHA256SUMS.txt"
    with path.open("w") as handle:
        for file in sorted(p for p in output_dir.iterdir() if p.is_file() and p.name != path.name):
            handle.write(f"{hashlib.sha256(file.read_bytes()).hexdigest()}  {file.name}\n")


def main():
    args = parse_args()
    output_dir = args.output_dir.resolve()
    primary, recurrent, pa, ra, summary = compute_results(args.wes_root.resolve())
    if not args.no_strict:
        validate_current_package(summary, pa, ra)
    configure_matplotlib()
    write_variant_tables(output_dir, primary, recurrent, pa, ra, summary)
    draw_membership_panel(summary, output_dir)
    draw_vaf_panel(primary, recurrent, pa, ra, summary, output_dir)
    write_checksums(output_dir)
    print(summary)


if __name__ == "__main__":
    main()
