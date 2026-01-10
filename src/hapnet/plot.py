from __future__ import annotations

import math
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from .haplotypes import Haplotype


def _make_population_colors(populations: List[str]) -> Dict[str, Tuple[float, float, float, float]]:
    pops = sorted(populations)
    cmap = plt.get_cmap("tab20")
    return {p: cmap(i % cmap.N) for i, p in enumerate(pops)}


def _normalize_and_spread(pos: Dict[int, np.ndarray], n: int, spread: float) -> Dict[int, np.ndarray]:
    coords = np.array([pos[i] for i in range(n)], dtype=float)
    coords = coords - coords.mean(axis=0, keepdims=True)

    norms = np.linalg.norm(coords, axis=1)
    span = float(np.max(norms)) if len(norms) else 1.0
    if span == 0.0:
        span = 1.0

    coords = (coords / span) * spread
    return {i: coords[i] for i in range(n)}


def _resolve_overlaps(
    pos: Dict[int, np.ndarray],
    radii: List[float],
    *,
    pad: float = 0.15,
    steps: int = 250,
    step_size: float = 0.25,
    seed: int = 0,
) -> Dict[int, np.ndarray]:
    """
    Simple collision-avoidance: iteratively pushes overlapping circles apart.
    Deterministic given seed.
    """
    rng = np.random.default_rng(seed)
    n = len(radii)
    P = np.array([pos[i] for i in range(n)], dtype=float)

    for _ in range(steps):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = P[j] - P[i]
                dist = float(np.linalg.norm(dx))
                min_dist = radii[i] + radii[j] + pad

                if dist == 0.0:
                    dx = rng.normal(size=2)
                    dist = float(np.linalg.norm(dx)) or 1.0

                if dist < min_dist:
                    direction = dx / dist
                    overlap = (min_dist - dist)

                    shift = direction * overlap * 0.5 * step_size
                    P[i] -= shift
                    P[j] += shift
                    moved = True

        if not moved:
            break

    return {i: P[i] for i in range(n)}


def _draw_edge_with_ticks(
    ax,
    p1: np.ndarray,
    p2: np.ndarray,
    *,
    weight: int,
    lw: float = 3.0,
    edge_color: str = "0.4",
    draw_ticks: bool = True,
    tick_length: float = 0.38,
    tick_gap: float = 0.55,
    tick_color: str = "0.15",
    tick_lw: float = 2.0,
    zorder: int = 1,
):
    """
    Draw an edge line and (optionally) mutation tick marks.

    weight: integer mutational distance. We draw (weight - 1) ticks so a
            1-step edge has no ticks (reduces clutter).
    """
    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])

    # Edge line
    ax.plot([x1, x2], [y1, y2], lw=lw, color=edge_color, zorder=zorder)

    if not draw_ticks:
        return

    n_ticks = max(0, int(weight) - 1)
    if n_ticks == 0:
        return

    v = np.array([x2 - x1, y2 - y1], dtype=float)
    L = float(np.linalg.norm(v))
    if L == 0.0:
        return

    u = v / L                      # unit along-edge
    perp = np.array([-u[1], u[0]]) # unit perpendicular

    # Center ticks at the midpoint. If the edge is short, compress spacing so ticks still fit.
    total_span = (n_ticks - 1) * tick_gap
    max_span = max(0.0, L - 2.0 * tick_gap)
    if total_span > max_span and n_ticks > 1:
        tick_gap_eff = max_span / (n_ticks - 1) if (n_ticks - 1) else tick_gap
    else:
        tick_gap_eff = tick_gap

    mid = np.array([(x1 + x2) / 2.0, (y1 + y2) / 2.0], dtype=float)
    start_offset = -0.5 * (n_ticks - 1) * tick_gap_eff

    for k in range(n_ticks):
        center = mid + (start_offset + k * tick_gap_eff) * u
        a = center - 0.5 * tick_length * perp
        b = center + 0.5 * tick_length * perp
        ax.plot([a[0], b[0]], [a[1], b[1]], lw=tick_lw, color=tick_color, zorder=zorder + 0.1)


def plot_network(
    G: nx.Graph,
    haplotypes: List[Haplotype],
    out: str = "hapnet.png",
    *,
    # Node sizing controls
    min_radius: float = 0.35,
    max_radius: float = 2.10,
    alpha: float = 0.35,            # <0.5 shrinks rare haplotypes strongly
    singleton_shrink: float = 0.60, # extra shrink for n=1
    # Layout / overlap controls
    spread: float | None = None,
    avoid_overlap: bool = True,
    overlap_pad: float = 0.25,
    overlap_steps: int = 400,
    overlap_step_size: float = 0.35,
    # Tick controls
    draw_ticks: bool = True,
    tick_length: float = 0.38,
    tick_gap: float = 0.55,
    # Aesthetics
    figsize: Tuple[float, float] = (10, 7),
    dpi: int = 300,
    show_counts_in_label: bool = False,
):
    if not haplotypes:
        raise ValueError("No haplotypes to plot.")

    n = len(haplotypes)

    # 1) Initial layout
    pos0 = nx.kamada_kawai_layout(G, weight="weight")

    # 2) Radii: compressed scaling for rare/private haplotypes
    max_n = max(h.n_total for h in haplotypes)

    def radius_for(count: int) -> float:
        if max_n <= 1:
            r = (min_radius + max_radius) / 2.0
        else:
            a = max(0.05, min(alpha, 1.0))
            t = (count / max_n) ** a
            r = min_radius + t * (max_radius - min_radius)

        if count == 1:
            r *= singleton_shrink

        return max(r, 0.18)

    radii = [radius_for(h.n_total) for h in haplotypes]
    rmax = max(radii)

    # 3) Normalize and spread
    if spread is None:
        spread = max(14.0, 4.0 * rmax * math.sqrt(n))

    pos = _normalize_and_spread({i: np.array(pos0[i]) for i in range(n)}, n, spread)

    # 4) Overlap resolution
    if avoid_overlap and n >= 3:
        pos = _resolve_overlaps(
            pos,
            radii,
            pad=overlap_pad,
            steps=overlap_steps,
            step_size=overlap_step_size,
            seed=0,
        )

    # 5) Colors
    all_pops = sorted({p for h in haplotypes for p in h.counts_by_pop.keys()})
    pop_colors = _make_population_colors(all_pops)

    # 6) Figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_axis_off()
    ax.set_aspect("equal")

    # 7) Edges + ticks
    for u, v, data in G.edges(data=True):
        w = data.get("weight", 1)
        try:
            w_int = int(round(float(w)))
        except Exception:
            w_int = 1

        _draw_edge_with_ticks(
            ax,
            np.array(pos[u]),
            np.array(pos[v]),
            weight=w_int,
            lw=3.2,
            edge_color="0.4",
            draw_ticks=draw_ticks,
            tick_length=tick_length,
            tick_gap=tick_gap,
            tick_color="0.15",
            tick_lw=2.0,
            zorder=1,
        )

    # 8) Nodes (pies)
    for i, h in enumerate(haplotypes):
        x, y = pos[i]
        r = radii[i]

        pops = sorted(h.counts_by_pop.keys())
        counts = [h.counts_by_pop[p] for p in pops]
        total = sum(counts)

        start = 0.0
        for pop, c in zip(pops, counts):
            frac = c / total
            wedge = plt.matplotlib.patches.Wedge(
                (float(x), float(y)),
                r,
                start * 360.0,
                (start + frac) * 360.0,
                facecolor=pop_colors[pop],
                edgecolor="black",
                linewidth=1.0,
                zorder=2,
            )
            ax.add_patch(wedge)
            start += frac

        outline = plt.matplotlib.patches.Circle(
            (float(x), float(y)),
            r,
            facecolor="none",
            edgecolor="black",
            linewidth=1.4,
            zorder=3,
        )
        ax.add_patch(outline)

        label = f"{h.hap_id}\n(n={h.n_total})" if show_counts_in_label else h.hap_id
        fs = 11 if not show_counts_in_label else 9
        ax.text(float(x), float(y), label, ha="center", va="center", fontsize=fs, zorder=4)

    # 9) Limits
    xs = [pos[i][0] for i in range(n)]
    ys = [pos[i][1] for i in range(n)]
    pad = 2.8 * rmax + 0.8
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)

    # 10) Legend
    handles = [
        plt.matplotlib.patches.Patch(facecolor=pop_colors[p], edgecolor="black", label=p)
        for p in all_pops
    ]
    ax.legend(
        handles=handles,
        title="Population",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=True,
    )

    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
