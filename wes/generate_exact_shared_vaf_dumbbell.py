#!/usr/bin/env python3
"""Draw the eight highlighted exact-shared WES variants as a VAF dumbbell plot.

The values are recalculated from tumor AD fields in the two PASS VCFs. Variant
membership uses the exact key CHROM + POS + REF + ALT, and annotations come from
the reviewed snpEff tables used by generate_wes_publication_figures.py.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np

from generate_wes_publication_figures import (
    Call,
    VariantKey,
    compute_results,
    configure_matplotlib,
    format_float,
    select_annotation,
    validate_current_package,
    write_checksums,
)


STEM = "CIC_WES_8_shared_MODERATE_HIGH_VAF_dumbbell"
TAB10 = (
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--wes-root",
        type=Path,
        required=True,
        help="WES_GATK analysis root containing results/variants and results/annotated",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir,
        help="Output directory (default: directory containing this script)",
    )
    parser.add_argument(
        "--no-strict",
        action="store_true",
        help="Allow inputs that no longer match the reviewed counts/gene set",
    )
    parser.add_argument(
        "--panel-label",
        default="B",
        help="Panel label; use an empty string to omit it (default: B)",
    )
    return parser.parse_args()


def compact_effect(effect: str) -> str:
    """Return the short publication label used in the supplied example."""
    if "splice_acceptor_variant" in effect:
        return "splice acceptor"
    if "frameshift_variant" in effect:
        return "frameshift"
    if "stop_gained" in effect:
        return "stop-gain"
    if "missense_variant" in effect:
        return "missense"
    return effect.replace("_variant", "").replace("_", " ")


def highlighted_rows(
    primary: Mapping[VariantKey, Call],
    recurrent: Mapping[VariantKey, Call],
    primary_annotations: Mapping[VariantKey, Mapping[str, str]],
    recurrent_annotations: Mapping[VariantKey, Mapping[str, str]],
    summary: Mapping[str, object],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key in summary["highlighted_keys"]:
        annotation = select_annotation(key, primary_annotations, recurrent_annotations)
        rows.append(
            {
                "key": key,
                "gene": annotation.get("GENE", "") or f"{key[0]}:{key[1]}",
                "effect": annotation.get("EFFECT", ""),
                "display_effect": compact_effect(annotation.get("EFFECT", "")),
                "impact": annotation.get("IMPACT", ""),
                "hgvs_c": annotation.get("HGVS_C", ""),
                "hgvs_p": annotation.get("HGVS_P", ""),
                "primary": primary[key],
                "recurrent": recurrent[key],
            }
        )
    rows.sort(key=lambda row: (-row["primary"].vaf, row["gene"]))
    return rows


def write_plot_data(rows: list[dict[str, object]], output_dir: Path) -> None:
    path = output_dir / f"{STEM}_data.tsv"
    fieldnames = [
        "DISPLAY_ORDER", "CHROM", "POS", "REF", "ALT", "GENE",
        "DISPLAY_EFFECT", "SNPEFF_EFFECT", "SNPEFF_IMPACT", "HGVS_C", "HGVS_P",
        "PRIMARY_REF_AD", "PRIMARY_ALT_AD", "PRIMARY_AD_DEPTH", "PRIMARY_VAF_AD",
        "RECURRENT_REF_AD", "RECURRENT_ALT_AD", "RECURRENT_AD_DEPTH", "RECURRENT_VAF_AD",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for order, row in enumerate(rows, start=1):
            key = row["key"]
            primary = row["primary"]
            recurrent = row["recurrent"]
            writer.writerow({
                "DISPLAY_ORDER": order,
                "CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3],
                "GENE": row["gene"], "DISPLAY_EFFECT": row["display_effect"],
                "SNPEFF_EFFECT": row["effect"], "SNPEFF_IMPACT": row["impact"],
                "HGVS_C": row["hgvs_c"], "HGVS_P": row["hgvs_p"],
                "PRIMARY_REF_AD": primary.ref_ad, "PRIMARY_ALT_AD": primary.alt_ad,
                "PRIMARY_AD_DEPTH": primary.ad_depth, "PRIMARY_VAF_AD": format_float(primary.vaf),
                "RECURRENT_REF_AD": recurrent.ref_ad, "RECURRENT_ALT_AD": recurrent.alt_ad,
                "RECURRENT_AD_DEPTH": recurrent.ad_depth, "RECURRENT_VAF_AD": format_float(recurrent.vaf),
            })


def marker_area(vaf: float) -> float:
    return 600.0 + 1000.0 * vaf


def draw_figure(rows: list[dict[str, object]], panel_label: str) -> plt.Figure:
    configure_matplotlib()
    fig, ax = plt.subplots(figsize=(12.0, 5.55))
    y_positions = np.arange(len(rows) - 1, -1, -1, dtype=float)

    for y, row, color in zip(y_positions, rows, TAB10):
        primary = row["primary"]
        recurrent = row["recurrent"]
        ax.plot([0, 1], [y, y], color=color, alpha=0.28, linewidth=1.0,
                solid_capstyle="round", zorder=1)
        ax.scatter([0, 1], [y, y],
                   s=[marker_area(primary.vaf), marker_area(recurrent.vaf)],
                   color=color, edgecolor=color, linewidth=1.0, alpha=0.88, zorder=2)
        for x, call in ((0, primary), (1, recurrent)):
            ax.text(x, y, f"{call.vaf:.1%}", ha="center", va="center",
                    fontsize=8.0, fontweight="bold", color="white", zorder=3)

    y_labels = []
    for row in rows:
        star = "*" if row["impact"] == "HIGH" else ""
        y_labels.append(f"{row['gene']}{star}  ({row['display_effect']})")
    ax.set_yticks(y_positions, labels=y_labels)
    ax.set_xticks([0, 1], labels=["Primary tumor", "Recurrent tumor"])
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(-0.76, len(rows) - 0.24)
    ax.tick_params(axis="y", length=0, pad=17, labelsize=10.2)
    ax.tick_params(axis="x", direction="out", length=4.2, width=0.8, pad=7, labelsize=10.3)
    for label in ax.get_yticklabels():
        label.set_horizontalalignment("right")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.spines["bottom"].set_color("#222222")
    ax.set_title("Eight exact shared MODERATE/HIGH PASS variants and their allele fractions",
                 fontsize=13.0, fontweight="semibold", pad=16)
    if panel_label:
        ax.text(-0.19, 1.055, panel_label, transform=ax.transAxes,
                ha="left", va="top", fontsize=18, fontweight="bold", color="#111111")

    fig.subplots_adjust(left=0.30, right=0.95, top=0.86, bottom=0.20)
    fig.text(
        0.63, 0.045,
        "* snpEff HIGH impact; this denotes predicted functional consequence, not an established cancer driver.\n"
        "VAF is AD-based and is not corrected for tumor purity or local copy-number state.",
        ha="center", va="bottom", fontsize=8.0, fontstyle="italic", color="#252525",
    )
    return fig


def save_figure(fig: plt.Figure, output_dir: Path) -> None:
    for extension in ("png", "pdf", "svg"):
        kwargs: dict[str, object] = {"bbox_inches": "tight", "facecolor": "white"}
        if extension == "png":
            kwargs["dpi"] = 600
        fig.savefig(output_dir / f"{STEM}.{extension}", **kwargs)


def main() -> None:
    args = parse_args()
    wes_root = args.wes_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    primary, recurrent, primary_annotations, recurrent_annotations, summary = compute_results(wes_root)
    if not args.no_strict:
        validate_current_package(summary, primary_annotations, recurrent_annotations)
    rows = highlighted_rows(primary, recurrent, primary_annotations, recurrent_annotations, summary)
    if len(rows) != 8:
        raise ValueError(f"Expected eight highlighted exact-shared variants; observed {len(rows)}")

    write_plot_data(rows, output_dir)
    fig = draw_figure(rows, args.panel_label)
    save_figure(fig, output_dir)
    plt.close(fig)
    write_checksums(output_dir)


if __name__ == "__main__":
    main()
