"""
PNG figures per sample (150 dpi): histograms of the QC metrics with the
thresholds, and spatial scatters of log total counts and of the removed
observations, on the lowres H&E when the sample has one. One colour system:
a single-hue blue ramp for magnitude, fixed categorical slots per removal
reason, the reserved critical red only for thresholds and removed marks,
muted inks for chrome.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from processor.qc import REASONS

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
CRITICAL = "#d03b3b"
KEPT = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
REASON_COLOUR = dict(zip(REASONS, SERIES))   # colour follows the reason, never its rank
REASON_LABEL = {"low_counts": "low counts", "low_genes": "few genes", "high_mt": "high % mt", "high_control": "control fraction"}
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": INK, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.spines.top": False, "axes.spines.right": False,
    })
    return plt


def _hist(ax, values, removed, bins, label, threshold=None, direction="below", log_x=False):
    plt = _plt()
    kept = values[~removed]
    gone = values[removed]
    ax.hist([kept, gone], bins=bins, stacked=True, color=[SERIES[0], CRITICAL], linewidth=0,
            label=[f"kept ({kept.size:,})", f"removed ({gone.size:,})"])
    if threshold is not None and np.isfinite(threshold):
        ax.axvline(threshold, color=CRITICAL, linewidth=1.2, linestyle=(0, (4, 2)))
        text = f"{threshold:,.0f}" if threshold >= 100 else f"{threshold:.3g}"
        ax.text(threshold, ax.get_ylim()[1] * 0.98, f" {'<' if direction == 'below' else '>'} {text}",
                ha="left", va="top", fontsize=8, color=CRITICAL)
    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel(label)
    ax.set_ylabel("observations")
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def qc_histograms(sample_id: str, kind: str, obs, removed: np.ndarray, thr, n_before: int, n_after: int, path: Path) -> Path:
    plt = _plt()
    has_control = "control_fraction" in obs
    n_panels = 4 if has_control else 3
    fig, axes = plt.subplots(1, n_panels, figsize=(3.4 * n_panels + 0.6, 3.3))
    counts = obs["total_counts"].to_numpy(dtype=float)
    genes = obs["n_genes_by_counts"].to_numpy(dtype=float)
    mt = obs["pct_counts_mt"].to_numpy(dtype=float)
    pos = counts[counts > 0]
    lo = max(1.0, pos.min()) if pos.size else 1.0
    hi = max(pos.max(), lo * 10) if pos.size else 10.0
    _hist(axes[0], np.maximum(counts, 1.0), removed, np.geomspace(lo, hi, 50), "total counts (log scale)",
          thr.min_counts, "below", log_x=True)
    _hist(axes[1], genes, removed, 50, "genes detected", thr.min_genes, "below")
    _hist(axes[2], mt, removed, 50, "% mitochondrial counts", thr.max_pct_mt, "above")
    if mt.size and np.nanmax(mt) == 0:
        axes[2].text(0.5, 0.55, "no mitochondrial counts\n(none in the panel?)", transform=axes[2].transAxes,
                     ha="center", va="center", fontsize=8.5, color=INK_2)
    if has_control:
        ctrl = obs["control_fraction"].to_numpy(dtype=float)
        _hist(axes[3], ctrl, removed, np.linspace(0, max(ctrl.max(), 0.2), 50), "control fraction",
              thr.max_control_fraction, "above")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, loc="upper right", ncol=2)
    pct = 100.0 * (n_before - n_after) / n_before if n_before else 0.0
    fig.suptitle(f"{sample_id} · QC ({kind}s): {n_before:,} → {n_after:,} observations, {pct:.1f} % removed · "
                 f"FILTER_MODE={thr.mode}", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _image_and_coords(adata, sample_id: str):
    """(image or None, coordinates to draw, scale factor): lowres H&E when the sample has one."""
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    try:
        entry = adata.uns["spatial"][sample_id]
        images, scales = entry["images"], entry["scalefactors"]
        for key, scale_key in (("lowres", "tissue_lowres_scalef"), ("hires", "tissue_hires_scalef")):
            if key in images and scale_key in scales:
                img = np.asarray(images[key])
                s = float(scales[scale_key])
                return img, xy * s, s
    except (KeyError, TypeError, AttributeError):
        pass
    return None, xy, 1.0


def spatial_counts(sample_id: str, adata_in, removed: np.ndarray, reasons, path: Path) -> Path:
    """Left: log1p total counts of every input observation; right: kept (grey) and removed by reason."""
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap

    img, xy, _ = _image_and_coords(adata_in, sample_id)
    n = xy.shape[0]
    size = float(np.clip(40000.0 / max(n, 1), 0.15, 14.0))
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    values = np.log1p(adata_in.obs["total_counts"].to_numpy(dtype=float))
    cmap = LinearSegmentedColormap.from_list("counts", SEQUENTIAL)
    for ax in axes:
        if img is not None:
            ax.imshow(img, interpolation="nearest")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    sca = axes[0].scatter(xy[:, 0], xy[:, 1], c=values, cmap=cmap, s=size, linewidths=0, rasterized=True)
    cb = fig.colorbar(sca, ax=axes[0], fraction=0.04, pad=0.02)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, colors=MUTED, length=0)
    cb.set_label("log1p(total counts)", fontsize=8, color=INK_2)
    axes[0].set_title(f"total counts ({n:,} observations)", loc="left")
    kept = ~removed
    axes[1].scatter(xy[kept, 0], xy[kept, 1], c=KEPT, s=size, linewidths=0, rasterized=True, label=f"kept ({kept.sum():,})")
    primary = reasons.str.split(";").str[0].to_numpy()
    for r in REASONS:   # fixed order, fixed colour per reason
        m = removed & (primary == r)
        if m.any():
            axes[1].scatter(xy[m, 0], xy[m, 1], c=REASON_COLOUR[r], s=max(size, 2.0), linewidths=0, rasterized=True,
                            label=f"{REASON_LABEL[r]} ({m.sum():,})")
    axes[1].set_title(f"removed ({removed.sum():,}), by first reason", loc="left")
    axes[1].legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0, -0.01), ncol=3,
                   markerscale=max(1.0, 6.0 / max(size, 1.0)))
    if img is None:
        for ax in axes:
            ax.invert_yaxis()   # image convention: y down
    fig.suptitle(f"{sample_id} · spatial QC", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
