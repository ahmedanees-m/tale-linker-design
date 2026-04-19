"""
Publication-quality figure generation for the C5 manuscript.

Implements Step 9 of the C5 Execution Plan.
All figures use matplotlib + seaborn. Plotly 3D for interactive HTML supplement.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Colorblind-friendly palette (Wong 2011)
CLASS_COLORS = {
    "F": "#0072B2",   # blue
    "H": "#D55E00",   # vermillion
    "M": "#009E73",   # green
    "N": "#CC79A7",   # pink
    "P": "#E69F00",   # orange
}
CLASS_LABELS = {
    "F": "Flexible (GGS)",
    "H": "Helical (EAAAK)",
    "M": "Mixed (GGS-EAAAK-GGS)",
    "N": "Natural (TALEN C+63)",
    "P": "Proline-rich (PG)",
}


# ─── Figure 1: TALE-DNA geometry and canonical reference frame ───────────────

def figure1_reference_frame(
    scissile_table: dict,
    output_path: str | Path,
) -> plt.Figure:
    """
    Figure 1: Schematic of TALE-DNA geometry with canonical reference frame axes
    and scissile phosphate positions annotated.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Left panel: 2D projection showing TALE superhelix and DNA axis
    ax = axes[0]
    ax.set_title("A   TALE-DNA Canonical Reference Frame")

    # Draw DNA double helix (simplified as two sinusoidal traces)
    z = np.linspace(0, 50, 300)
    twist_rad = np.radians(36 / 3.4) * z   # helical twist
    r_helix = 10.0
    x_top = r_helix * np.cos(twist_rad)
    x_bot = r_helix * np.cos(twist_rad + np.pi)

    ax.plot(x_top, z, color="steelblue", lw=1.5, label="Top strand (sense)")
    ax.plot(x_bot, z, color="tomato", lw=1.5, label="Bottom strand (antisense)")

    # TALE body (simplified as a rectangle below origin)
    tale_patch = mpatches.FancyBboxPatch((-32, -45), 64, 45,
                                          boxstyle="round,pad=2", linewidth=1.5,
                                          edgecolor="#444", facecolor="#D8EAF5", alpha=0.8)
    ax.add_patch(tale_patch)
    ax.text(0, -22, "TALE\nArray", ha="center", va="center",
            fontsize=10, fontweight="bold", color="#1a4b8c")

    # Reference frame axes
    ax.annotate("", xy=(0, 28), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="black", lw=2))
    ax.text(1, 22, "Z (3′ direction)", fontsize=8, color="black")

    ax.annotate("", xy=(22, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=2))
    ax.text(14, -3, "X (major groove)", fontsize=8, color="#2ca02c")

    # Origin
    ax.scatter([0], [0], s=80, c="black", zorder=5)
    ax.text(1.5, -3, "Origin\n(TALE C-term Cα)", fontsize=7.5)

    # Scissile phosphate positions
    for (strand, bp), entry in scissile_table.items():
        if bp > 6:
            continue
        coords = entry["coords"]
        color = CLASS_COLORS.get("F" if strand == "top" else "H", "grey")
        ax.scatter([coords[0]], [coords[2]], s=40, c=color, zorder=4, alpha=0.9)
        if strand == "top":
            ax.text(coords[0] + 1.5, coords[2], f"+{bp}", fontsize=7)

    ax.scatter([], [], s=40, c=CLASS_COLORS["F"], label="Top strand phosphates")
    ax.scatter([], [], s=40, c=CLASS_COLORS["H"], label="Bottom strand phosphates")

    ax.set_xlim(-36, 36)
    ax.set_ylim(-50, 55)
    ax.set_xlabel("X — major groove axis (Å)")
    ax.set_ylabel("Z — helical axis (Å)")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.set_aspect("equal")

    # Right panel: distance from C-terminus to each phosphate
    ax2 = axes[1]
    ax2.set_title("B   Distance from TALE C-terminus to Scissile Phosphates")

    for strand, color, label in [("top", CLASS_COLORS["F"], "Top strand"),
                                   ("bottom", CLASS_COLORS["H"], "Bottom strand")]:
        bps = sorted([bp for (s, bp) in scissile_table if s == strand])
        dists = [scissile_table[(strand, bp)]["distance_from_origin"] for bp in bps]
        ax2.plot(bps, dists, "o-", color=color, label=label, lw=2, ms=6)

    # GENESIS target region
    ax2.axvspan(3, 5, alpha=0.15, color="#FF8C00", label="GENESIS target region")
    ax2.axvline(4, ls="--", color="#FF8C00", lw=1.5)
    ax2.text(4.05, 2, "Primary\ntarget\n(bp +4)", fontsize=7.5, color="#FF8C00")

    ax2.set_xlabel("Base pair offset from TALE footprint")
    ax2.set_ylabel("Distance to TALE C-terminus (Å)")
    ax2.legend()
    ax2.set_xticks(range(0, 11))

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    return fig


# ─── Figure 2: Linker library characterization ───────────────────────────────

def figure2_linker_library(
    linker_library: dict,
    output_path: str | Path,
) -> plt.Figure:
    """
    Figure 2: End-to-end distance distributions for each linker class and length.
    """
    classes = ["F", "H", "M", "N", "P"]
    lengths_to_plot = [5, 10, 15, 20, 30]

    fig, axes = plt.subplots(1, len(classes), figsize=(14, 4), sharey=False)
    r_vals = np.linspace(0, 120, 500)

    for ax, cls_str in zip(axes, classes):
        from .linkers import LinkerClass
        cls = LinkerClass(cls_str)
        for n in lengths_to_plot:
            ens = linker_library.get((cls, n))
            if ens is None:
                continue
            pdf = ens.end_to_end_pdf(r_vals)
            pdf_norm = pdf / (pdf.max() + 1e-30)
            ax.plot(r_vals, pdf_norm, color=CLASS_COLORS[cls_str],
                    alpha=0.3 + 0.14 * lengths_to_plot.index(n),
                    label=f"n={n}", lw=1.5)
        ax.set_title(f"{CLASS_LABELS[cls_str]}", fontsize=9)
        ax.set_xlabel("End-to-end distance (Å)")
        if ax == axes[0]:
            ax.set_ylabel("Normalized probability density")
        ax.legend(fontsize=7, title="Length (res)", title_fontsize=7)

    axes[2].set_title("C   Mixed (GGS–EAAAK–GGS)", fontsize=9)
    panel_labels = "ABCDE"
    for i, ax in enumerate(axes):
        ax.set_title(f"{panel_labels[i]}   {CLASS_LABELS[classes[i]]}", fontsize=9)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    return fig


# ─── Figure 3: Reachability maps ─────────────────────────────────────────────

def figure3_reachability_maps(
    reachability_maps: dict,
    scissile_table: dict,
    output_path: str | Path,
) -> plt.Figure:
    """
    Figure 3: 2D projections (XZ and YZ) of reachability maps for representative linkers.
    """
    representatives = [
        ("F", 15), ("H", 15), ("M", 15), ("N", 20), ("F", 25)
    ]

    ncols = len(representatives)
    fig, axes = plt.subplots(2, ncols, figsize=(14, 7))

    for col, (cls_str, n) in enumerate(representatives):
        rm = reachability_maps.get((cls_str, n))
        label = f"{CLASS_LABELS.get(cls_str, cls_str)}\nn={n} res"

        for row, (xi, yi, xlabel, ylabel) in enumerate([
            (0, 2, "X (Å)", "Z (Å)"),   # XZ projection
            (1, 2, "Y (Å)", "Z (Å)"),   # YZ projection
        ]):
            ax = axes[row, col]
            if rm is None or len(rm.samples) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9)
                continue

            samples = rm.samples
            x_data = samples[:, xi]
            y_data = samples[:, yi]

            ax.hexbin(x_data, y_data, gridsize=40, cmap="Blues",
                      mincnt=1, linewidths=0.2)

            # Plot target positions
            for (strand, bp), entry in scissile_table.items():
                if bp > 6:
                    continue
                coords = entry["coords"]
                marker = "*" if (strand == "top" and bp == 4) else "o"
                ms = 120 if marker == "*" else 40
                ax.scatter([coords[xi]], [coords[yi]], marker=marker,
                           s=ms, c="#FF4500", zorder=5, edgecolors="k", lw=0.5)

            if row == 0:
                ax.set_title(label, fontsize=8.5)
            ax.set_xlabel(xlabel)
            if col == 0:
                ax.set_ylabel(ylabel)
            ax.set_xlim(-60, 60)
            ax.set_ylim(-20, 90)

    # Legend: star = primary target bp+4 top
    star = mpatches.Patch(color="#FF4500", label="Scissile phosphates\n(★ = bp+4 primary target)")
    fig.legend(handles=[star], loc="lower center", ncol=1, fontsize=8)

    fig.suptitle("Figure 3. Reachability maps — catalytic domain attachment point density",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    return fig


# ─── Figure 4: Published fusion benchmark ────────────────────────────────────

def figure4_published_benchmark(
    published_data: list[dict],
    reachability_maps: dict,
    scissile_table: dict,
    output_path: str | Path,
) -> plt.Figure:
    """
    Figure 4: Predicted geometric distance vs. reported efficiency for published fusions.
    """
    from scipy.stats import spearmanr

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    names, distances, efficiencies, classes = [], [], [], []
    for entry in published_data:
        cls_str = entry["linker_class"]
        n = entry["n_residues"]
        strand = entry.get("strand", "top")
        bp = entry["cut_site_bp_offset"]
        eff = entry["reported_efficiency_pct"]

        rm = reachability_maps.get((cls_str, n))
        sc_key = (strand, min(bp, 10))
        if rm is None or sc_key not in scissile_table:
            continue

        target_coords = scissile_table[sc_key]["coords"]
        dist = rm.nearest_sample_distance(target_coords)
        if not np.isfinite(dist):
            continue

        names.append(entry.get("name", "?"))
        distances.append(dist)
        efficiencies.append(eff)
        classes.append(cls_str)

    distances = np.array(distances)
    efficiencies = np.array(efficiencies)

    # Left: scatter
    ax = axes[0]
    for cls_str in set(classes):
        mask = [c == cls_str for c in classes]
        ax.scatter(
            distances[mask], efficiencies[mask],
            c=CLASS_COLORS.get(cls_str, "grey"),
            label=CLASS_LABELS.get(cls_str, cls_str),
            s=60, alpha=0.85, edgecolors="k", lw=0.5,
        )

    if len(distances) >= 3:
        rho, pval = spearmanr(distances, efficiencies)
        ax.text(0.98, 0.95, f"ρ = {rho:.2f}\np = {pval:.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    ax.set_xlabel("Nearest attachment-point distance to cut site (Å)")
    ax.set_ylabel("Reported editing efficiency (%)")
    ax.set_title("A   Geometric distance vs. editing efficiency")
    ax.legend(fontsize=7)

    # Right: bar chart of distance by system
    ax2 = axes[1]
    order = np.argsort(distances)
    bar_names = [names[i][:22] for i in order]
    bar_dists = distances[order]
    bar_colors = [CLASS_COLORS.get(classes[i], "grey") for i in order]

    bars = ax2.barh(range(len(bar_dists)), bar_dists, color=bar_colors, alpha=0.85)
    ax2.set_yticks(range(len(bar_dists)))
    ax2.set_yticklabels(bar_names, fontsize=8)
    ax2.set_xlabel("Nearest attachment-point distance to cut site (Å)")
    ax2.set_title("B   Ranked geometric reach by published system")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    return fig


# ─── Figure 5: GENESIS-specific linker recommendations ───────────────────────

def figure5_genesis_recommendations(
    designs: list,
    reachability_maps: dict,
    scissile_table: dict,
    output_path: str | Path,
) -> plt.Figure:
    """
    Figure 5: Top-3 GENESIS linker designs — reachability, orientation, and geometry.
    """
    ncols = 3
    fig, axes = plt.subplots(1, ncols, figsize=(12, 4.5))

    for col, design in enumerate(designs[:ncols]):
        ax = axes[col]
        rm = reachability_maps.get((design.linker_class.value, design.n_residues))

        if rm is not None and len(rm.samples) > 0:
            samples = rm.samples
            ax.hexbin(samples[:, 0], samples[:, 2], gridsize=35, cmap="Blues",
                      mincnt=1, linewidths=0.2, alpha=0.9)

        # Target position bp+4 top strand
        primary = scissile_table.get(("top", 4))
        if primary is not None:
            ax.scatter([primary["coords"][0]], [primary["coords"][2]],
                       marker="*", s=250, c="#FF4500", zorder=6,
                       edgecolors="k", lw=0.7, label="bp+4 top (primary)")

        ax.set_title(
            f"Design {col+1}  ({design.linker_class.value}, n={design.n_residues})\n"
            f"P(reach) = {design.p_reach_primary_pct:.1f}%   "
            f"Nearest = {design.nearest_A:.1f} Å",
            fontsize=8.5
        )
        ax.set_xlabel("X (Å)")
        if col == 0:
            ax.set_ylabel("Z (Å)")
        ax.set_xlim(-55, 55)
        ax.set_ylim(-10, 80)
        ax.legend(fontsize=7)

    fig.suptitle("Figure 5. GENESIS-specific linker design recommendations for bp+4 target",
                 fontsize=10)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    return fig


# ─── Figure 6: Decision flowchart ────────────────────────────────────────────

def figure6_decision_flowchart(output_path: str | Path) -> plt.Figure:
    """
    Figure 6: TALE-fusion linker design decision tree (matplotlib-rendered).
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Figure 6. Decision framework for TALE-fusion linker design", fontsize=11)

    boxes = [
        (5.0, 9.2, "Define catalytic domain position\nrelative to TALE footprint", "#AED6F1"),
        (5.0, 7.8, "Compute scissile phosphate\ncoordinates (canonical frame)", "#AED6F1"),
        (5.0, 6.4, "Query reachability maps:\nWhich (class, length) pairs\nhave P(reach) > 10%?", "#D5E8D4"),
        (5.0, 4.8, "Required directionality?", "#FFF2CC"),
        (2.2, 3.2, "YES → Helical (H)\nor Mixed (M)", "#F8CECC"),
        (7.8, 3.2, "NO → Flexible (F)\nor Natural (N)", "#DAE8FC"),
        (5.0, 1.8, "Select shortest linker\nmeeting P(reach) threshold", "#E1D5E7"),
        (5.0, 0.5, "Freeze linker specification\nfor C1/C2 design", "#D5E8D4"),
    ]

    for (x, y, text, color) in boxes:
        w, h = 3.2, 0.9
        rect = mpatches.FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.1", linewidth=1.2,
            edgecolor="#555", facecolor=color
        )
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=8.5, wrap=True)

    # Arrows
    arrow_pairs = [
        (5, 8.75, 5, 8.25),
        (5, 7.35, 5, 6.85),
        (5, 5.95, 5, 5.3),
        (3.8, 4.55, 2.2, 3.65),
        (6.2, 4.55, 7.8, 3.65),
        (2.2, 2.75, 5, 2.25),
        (7.8, 2.75, 5, 2.25),
        (5, 1.35, 5, 0.95),
    ]
    for x1, y1, x2, y2 in arrow_pairs:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="#333", lw=1.2))

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight")
    return fig


def plot_reachability(
    tale_structure,
    linker_class: str,
    length: int,
    reachability_maps: dict,
    scissile_table: dict,
) -> plt.Figure:
    """Quick-look reachability plot for a single linker class/length (API helper)."""
    rm = reachability_maps.get((linker_class, length))
    fig, ax = plt.subplots(figsize=(5, 5))
    if rm is not None and len(rm.samples):
        ax.hexbin(rm.samples[:, 0], rm.samples[:, 2], gridsize=30, cmap="Blues", mincnt=1)
    primary = scissile_table.get(("top", 4))
    if primary is not None:
        ax.scatter([primary["coords"][0]], [primary["coords"][2]], marker="*",
                   s=200, c="red", zorder=5, label="bp+4 target")
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Z (Å)")
    ax.set_title(f"Reachability: class={linker_class}, n={length} res")
    ax.legend()
    return fig


def plot_linker_distributions(linker_library: dict) -> plt.Figure:
    """Quick-look figure of end-to-end distance distributions."""
    fig, ax = plt.subplots(figsize=(7, 4))
    r_vals = np.linspace(0, 120, 400)
    from .linkers import LinkerClass
    for cls_str in ["F", "H", "M"]:
        cls = LinkerClass(cls_str)
        for n in [10, 20, 30]:
            ens = linker_library.get((cls, n))
            if ens:
                pdf = ens.end_to_end_pdf(r_vals)
                pdf /= (pdf.max() + 1e-30)
                ax.plot(r_vals, pdf, color=CLASS_COLORS[cls_str],
                        alpha=0.4 + 0.2 * [10, 20, 30].index(n),
                        label=f"{cls_str} n={n}")
    ax.set_xlabel("End-to-end distance (Å)")
    ax.set_ylabel("Normalized P(r)")
    ax.legend(fontsize=7)
    ax.set_title("Linker end-to-end distance distributions")
    return fig
