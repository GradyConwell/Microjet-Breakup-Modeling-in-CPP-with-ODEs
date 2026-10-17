#!/usr/bin/env python3
"""
Microjet Breakup Visualization Toolchain
========================================
Standalone, publication-ready figure generator for 1D lubrication-approximation
simulations of Rayleigh microjets (Van Hoeve et al., 2010).

Design Principles:
- Standalone figures: Every chart is saved as its own individual image file
  inside a dedicated per-run directory (`plots/<csv_stem>_<timestamp>/`) ready
  for inclusion in a paper or technical report.
- Typographic hierarchy:
  - Titles / section accents: Editorial serif (`Fraunces`, `Martina Plantijn`,
    `New York`, or `Charter`), normal or light weight, sentence case.
  - Axis labels, legends & prose annotations: Grotesque sans-serif (`Inter` or
    `Helvetica Neue`), sentence case (never all-caps).
  - Numerals, ticks & dimension values: Technical monospace (`ABC Monument Grotesk Mono`,
    `Monument Grotesk Mono`, or `SF Mono`).
- Generous whitespace, crisp ISO-weight contours, zero decorative fluff.

Programmatic Usage Example:
---------------------------
    import plotFrames as pf

    cfg = pf.PlotConfig(
        time_step_us=2.0,
        steady_span_us=100.0,
        color_by_velocity=False,
        show_titles=True,
    )
    ds = pf.MicrojetDataset.from_csv("microjet10-15-2026.csv", cfg)
    outputs = pf.generate_all_figures(ds, cfg)
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
from dataclasses import dataclass, field
from datetime import datetime
import glob
import math
import os
import shutil
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.gridspec as gridspec


# =============================================================================
# 0. Font Resolution Helpers
# =============================================================================

def _register_system_extra_fonts():
    """Registers newly installed user and system fonts with Matplotlib's fontManager."""
    user_font_dir = os.path.expanduser("~/Library/Fonts")
    patterns = [
        os.path.join(user_font_dir, "Splentinex*.otf"),
        os.path.join(user_font_dir, "IBMPlexMono*.ttf"),
        "/Library/Fonts/Artifakt Element*.ttf",
    ]
    for pat in patterns:
        for font_path in glob.glob(pat):
            try:
                fm.fontManager.addfont(font_path)
            except Exception:
                pass


_register_system_extra_fonts()


def _resolve_font(candidates: Sequence[str], fallback: str = "sans-serif") -> str:
    """Returns the first installed font family matching `candidates`."""
    available = {f.name for f in fm.fontManager.ttflist}
    for c in candidates:
        if c in available:
            return c
    return fallback


DEFAULT_TITLE_FONT = _resolve_font(
    ("Splentinex", "New York", "Charter", "Georgia"),
    fallback="serif",
)
DEFAULT_ACCENT_FONT = _resolve_font(
    ("Splentinex", "New York", "Charter"),
    fallback="serif",
)
DEFAULT_BODY_FONT = _resolve_font(
    ("Artifakt Element", "Helvetica Neue", "Helvetica", "Arial"),
    fallback="sans-serif",
)
DEFAULT_BODY_LIGHT_FONT = _resolve_font(
    ("Artifakt Element", "Helvetica Neue", "Helvetica"),
    fallback="sans-serif",
)
DEFAULT_MONO_FONT = _resolve_font(
    ("IBM Plex Mono", "PT Mono", "Menlo"),
    fallback="monospace",
)


# =============================================================================
# 1. Configuration Dataclasses (Programmatic Control API)
# =============================================================================

@dataclass
class StyleConfig:
    """
    Visual and typographic parameters controlling figure presentation.
    - `title_font`: Editorial serif (`Splentinex`) for titles
    - `accent_font`: Serif (`Splentinex`) for subtle accents
    - `body_font`: Technical sans-serif (`Artifakt Element`) for labels & prose
    - `mono_font`: Technical monospace (`IBM Plex Mono`) strictly for numerals & ticks
    """
    name: str = "cad_light"
    title_font: str = DEFAULT_TITLE_FONT
    accent_font: str = DEFAULT_ACCENT_FONT
    body_font: str = DEFAULT_BODY_FONT
    body_light_font: str = DEFAULT_BODY_LIGHT_FONT
    mono_font: str = DEFAULT_MONO_FONT

    bg_color: str = "#FFFFFF"
    panel_bg_color: str = "#FFFFFF"
    grid_color: str = "#EEF2F6"
    axis_color: str = "#111827"
    text_color: str = "#111827"
    muted_text_color: str = "#52525B"

    # Jet & droplet geometry colors
    jet_edge_color: str = "#111827"
    jet_fill_color: str = "#D8E0E9"
    primary_edge_color: str = "#111827"
    primary_fill_color: str = "#1E3A8A"
    satellite_edge_color: str = "#991B1B"
    satellite_fill_color: str = "#DC2626"

    # Callouts, datum lines & trajectories
    datum_color: str = "#0369A1"
    locus_color: str = "#B91C1C"
    traj_primary_color: str = "#1D4ED8"
    traj_satellite_color: str = "#B91C1C"
    centerline_color: str = "#9CA3AF"
    dimension_color: str = "#111827"

    # Line weights (points)
    lw_border: float = 0.65
    lw_contour: float = 0.50
    lw_datum: float = 0.60
    lw_centerline: float = 0.35
    lw_grid: float = 0.40
    lw_dimension: float = 0.65

    # Velocity colormap stops
    velocity_cmap_colors: Tuple[str, ...] = (
        "#0284C7",  # Recoil / negative velocity
        "#E2E8F0",  # Stagnation
        "#334155",  # Nominal jet velocity v0
        "#DC2626",  # Capillary pinch-off acceleration peak
    )

    @classmethod
    def cad_light(cls) -> "StyleConfig":
        """Crisp white technical drawing style with hard ink lines."""
        return cls(name="cad_light")

    @classmethod
    def cad_blueprint(cls) -> "StyleConfig":
        """Dark slate technical drafting theme."""
        return cls(
            name="cad_blueprint",
            bg_color="#0B1120",
            panel_bg_color="#0F172A",
            grid_color="#1E293B",
            axis_color="#E2E8F0",
            text_color="#F8FAFC",
            muted_text_color="#94A3B8",
            jet_edge_color="#38BDF8",
            jet_fill_color="#1E3A5F",
            primary_edge_color="#38BDF8",
            primary_fill_color="#0284C7",
            satellite_edge_color="#FBBF24",
            satellite_fill_color="#F59E0B",
            datum_color="#38BDF8",
            locus_color="#F59E0B",
            traj_primary_color="#38BDF8",
            traj_satellite_color="#F59E0B",
            centerline_color="#475569",
            dimension_color="#E2E8F0",
            velocity_cmap_colors=("#F59E0B", "#1E293B", "#0284C7", "#38BDF8"),
        )

    @classmethod
    def paper_classic(cls) -> "StyleConfig":
        """Minimal academic journal style matching Van Hoeve et al. (2010)."""
        return cls(
            name="paper_classic",
            bg_color="#FFFFFF",
            panel_bg_color="#FFFFFF",
            grid_color="#F4F4F5",
            axis_color="#000000",
            text_color="#000000",
            muted_text_color="#3F3F46",
            jet_edge_color="#1E293B",
            jet_fill_color="#334B8E",
            primary_edge_color="#1E293B",
            primary_fill_color="#334B8E",
            satellite_edge_color="#1E293B",
            satellite_fill_color="#334B8E",
            datum_color="#3F3F46",
            locus_color="#991B1B",
            traj_primary_color="#334B8E",
            traj_satellite_color="#991B1B",
            centerline_color="#A1A1AA",
            dimension_color="#000000",
        )


@dataclass
class PlotConfig:
    """
    Master configuration for physics reconstruction, sampling, and figure layout.
    """
    # Physical parameters (Van Hoeve et al., 2010, Sec. III.A default)
    rho: float = 1098.0               # kg/m^3
    gamma: float = 0.0679             # N/m
    mu: float = 0.00365               # Pa*s
    h0_m: float = 18.5e-6             # m (nozzle radius)
    v0_ms: float = 5.4                # m/s (jet velocity)
    amp: float = 0.05                 # dimensionless perturbation amplitude
    domain_length_m: float = 700e-6   # m (axial domain length)
    pinch_threshold_m: float = 1.5e-6 # m (neck pinch-off threshold)

    # Droplet tracking & coalescence rules (Van Hoeve et al., 2010, Sec. II & Fig. 5)
    coalesce_droplets: bool = True
    coalesce_primary_satellite: bool = False  # Paper Fig. 5: primary-satellite merging excluded
    absorb_recoiling_at_tip: bool = True      # Re-absorb microdroplets that recoil into jet tip
    crit_distance_h0_frac: float = 0.005      # Critical overlap distance 0.005 * h0
    primary_radius_cutoff_um: float = 14.0    # Threshold distinguishing primary vs satellite (µm)
    min_render_radius_um: float = 0.5         # Minimum droplet radius to render (µm)

    # Stroboscopic montage sampling, windows & coordinate scaling
    time_step_us: float = 2.0
    t_min_us: Optional[float] = None
    t_max_us: Optional[float] = None
    steady_start_us: Optional[float] = None   # Defaults to clean 10 µs multiple after first exit (~150 µs)
    steady_span_us: float = 80.0              # Focused 80 µs window (~2.6 Rayleigh periods)
    z_min_um: float = -28.0                   # Bottom breathing room below nozzle (z = 0)
    z_max_um: float = 840.0                   # Top breathing room above domain exit (z = 700 µm)
    aspect_um_per_us: float = 44.0            # Balanced spatial-to-temporal aspect ratio

    # Layer visibility toggles
    show_titles: bool = True                  # Quiet editorial serif title above plot
    color_by_velocity: bool = False           # Shade stroboscopic jet interior by v(z, t)
    show_grid: bool = True                    # Show subtle grid on line plots
    show_montage_grid: bool = False           # Keep stroboscopic montage background pure white by default
    show_centerlines: bool = False
    show_trajectories: bool = False           # Optional velocity slope triangle callout
    show_breakup_locus: bool = True           # Show mean breakup length datum line
    show_dimension_callouts: bool = True
    show_scale_bar: bool = True
    show_z_axis: bool = True
    highlight_satellites: bool = True

    # Output settings
    dpi: int = 350
    pad_inches: float = 0.26                  # Generous outer margin around standalone figures
    output_root: str = "plots"
    also_save_legacy_names: bool = True       # Keep microjet_before_leave.png & microjet_after_leave.png updated
    style: StyleConfig = field(default_factory=StyleConfig.cad_light)

    @property
    def scale_x(self) -> float:
        return 1.0 / self.aspect_um_per_us

    @property
    def lambda_opt_um(self) -> float:
        return 2.0 * math.sqrt(2.0) * math.pi * (self.h0_m * 1e6)

    @property
    def freq_khz(self) -> float:
        return (self.v0_ms / (2.0 * math.sqrt(2.0) * math.pi * self.h0_m)) * 1e-3

    @property
    def tau_c_us(self) -> float:
        return math.sqrt(self.rho * (self.h0_m ** 3) / self.gamma) * 1e6

    @property
    def weber(self) -> float:
        return self.rho * (self.v0_ms ** 2) * self.h0_m / self.gamma

    @property
    def ohnesorge(self) -> float:
        return self.mu / math.sqrt(self.rho * self.gamma * self.h0_m)

    @property
    def reynolds(self) -> float:
        return self.rho * self.v0_ms * self.h0_m / self.mu

    @property
    def rayleigh_primary_radius_um(self) -> float:
        vol_um3 = math.pi * ((self.h0_m * 1e6) ** 2) * self.lambda_opt_um
        return ((3.0 * vol_um3) / (4.0 * math.pi)) ** (1.0 / 3.0)


# =============================================================================
# 2. Data Structures & Physics-Accurate Simulation Ingestion
# =============================================================================

@dataclass
class DropletState:
    """Instantaneous state of a detached droplet at a single time slice."""
    droplet_id: int
    z_m: float
    r_m: float
    v_ms: float
    vol_m3: float
    t_birth_s: float
    is_primary: bool
    coalescence_count: int = 0

    @property
    def z_um(self) -> float:
        return self.z_m * 1e6

    @property
    def r_um(self) -> float:
        return self.r_m * 1e6


@dataclass
class PinchoffEvent:
    """Metadata recorded whenever the jet neck pinches off."""
    event_id: int
    row_index: int
    t_s: float
    neck_node: int
    neck_lo: int
    drop_end: int
    z_pinch_m: float
    z_centroid_m: float
    r_m: float
    v_ms: float
    vol_m3: float
    is_primary: bool

    @property
    def t_us(self) -> float:
        return self.t_s * 1e6

    @property
    def z_pinch_um(self) -> float:
        return self.z_pinch_m * 1e6

    @property
    def r_um(self) -> float:
        return self.r_m * 1e6


class MicrojetDataset:
    """
    Holds the full spatiotemporal state of a microjet simulation, including:
    - Raw and truncated jet profiles h(z, t) and velocity fields v(z, t)
    - Exact pinch-off events matching the C++ solver's `performPinchoff`
    - Time-resolved droplet trajectories with coalescence and tip-recoil absorption
    """

    def __init__(
        self,
        csv_path: str,
        times_s: np.ndarray,
        h_grid_m: np.ndarray,
        v_grid_ms: np.ndarray,
        last_nodes: np.ndarray,
        pinchoff_events: List[PinchoffEvent],
        droplets_by_row: List[List[DropletState]],
        config: PlotConfig,
    ):
        self.csv_path = csv_path
        self.times_s = times_s
        self.times_us = times_s * 1e6
        self.h_grid_m = h_grid_m
        self.v_grid_ms = v_grid_ms
        self.last_nodes = last_nodes
        self.pinchoff_events = pinchoff_events
        self.droplets_by_row = droplets_by_row
        self.config = config

        self.num_rows, self.num_nodes = h_grid_m.shape
        self.dz_m = config.domain_length_m / (self.num_nodes - 1)
        self.z_grid_m = np.arange(self.num_nodes, dtype=float) * self.dz_m
        self.z_grid_um = self.z_grid_m * 1e6

    @classmethod
    def from_csv(cls, csv_path: Optional[str] = None, config: Optional[PlotConfig] = None) -> "MicrojetDataset":
        if config is None:
            config = PlotConfig()
        if csv_path is None:
            csv_path = find_latest_csv()

        with open(csv_path, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            raw_rows = [list(map(float, r)) for r in reader if r]

        if not raw_rows:
            raise ValueError(f"CSV file '{csv_path}' contains no data rows.")

        data = np.array(raw_rows, dtype=float)
        times_s = data[:, 0]
        num_nodes = (len(header) - 1) // 2
        h_grid_m = data[:, 1::2][:, :num_nodes].copy()
        v_grid_ms = data[:, 2::2][:, :num_nodes].copy()

        dz_m = config.domain_length_m / (num_nodes - 1)

        bad_mask = ~np.isfinite(h_grid_m) | ~np.isfinite(v_grid_ms) | (h_grid_m < 0.0) | (h_grid_m > 100e-6)
        if np.any(bad_mask):
            h_grid_m[bad_mask] = 0.0
            v_grid_ms[bad_mask] = 0.0

        active_mask = h_grid_m > 1e-12
        last_nodes = np.zeros(len(times_s), dtype=int)
        for r_idx in range(len(times_s)):
            nz = np.flatnonzero(active_mask[r_idx])
            last_nodes[r_idx] = int(nz[-1]) if len(nz) > 0 else 0

        pinchoff_events: List[PinchoffEvent] = []
        droplets_by_row: List[List[DropletState]] = [[] for _ in range(len(times_s))]
        live_droplets: List[DropletState] = []
        next_drop_id = 1
        crit_dist_m = config.crit_distance_h0_frac * config.h0_m
        prim_cutoff_m = config.primary_radius_cutoff_um * 1e-6

        for r_idx in range(len(times_s)):
            t_cur = times_s[r_idx]
            dt = (t_cur - times_s[r_idx - 1]) if r_idx > 0 else 0.0

            if dt > 0.0:
                for d in live_droplets:
                    d.z_m += d.v_ms * dt

            h_row = h_grid_m[r_idx]
            v_row = v_grid_ms[r_idx]
            last_n = last_nodes[r_idx]

            pinched_here = False
            neck_lo = -1
            drop_end = last_n
            k_min = 1

            if last_n >= 3:
                interior_h = h_row[1:last_n]
                k_min = 1 + int(np.argmin(interior_h))
                if h_row[k_min] <= config.pinch_threshold_m:
                    pinched_here = True
                    neck_lo = k_min
                    while neck_lo > 1 and h_row[neck_lo - 1] <= config.pinch_threshold_m:
                        neck_lo -= 1

            if not pinched_here and r_idx > 0 and last_n < last_nodes[r_idx - 1]:
                prev_last = last_nodes[r_idx - 1]
                h_prev = h_grid_m[r_idx - 1]
                v_prev = v_grid_ms[r_idx - 1]
                neck_lo = last_n + 1
                drop_end = prev_last
                k_min = neck_lo
                if neck_lo <= drop_end:
                    dV = math.pi * (h_prev[neck_lo:drop_end + 1] ** 2) * dz_m
                    vol = float(np.sum(dV))
                    if vol > 0.0:
                        v_d = float(np.sum(dV * v_prev[neck_lo:drop_end + 1]) / vol)
                        z_d = float(np.sum(dV * (np.arange(neck_lo, drop_end + 1) * dz_m)) / vol)
                        r_d = (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0)
                        is_prim = r_d >= prim_cutoff_m
                        ev = PinchoffEvent(
                            event_id=len(pinchoff_events) + 1,
                            row_index=r_idx,
                            t_s=t_cur,
                            neck_node=k_min,
                            neck_lo=neck_lo,
                            drop_end=drop_end,
                            z_pinch_m=k_min * dz_m,
                            z_centroid_m=z_d,
                            r_m=r_d,
                            v_ms=v_d,
                            vol_m3=vol,
                            is_primary=is_prim,
                        )
                        pinchoff_events.append(ev)
                        live_droplets.append(
                            DropletState(
                                droplet_id=next_drop_id,
                                z_m=z_d,
                                r_m=r_d,
                                v_ms=v_d,
                                vol_m3=vol,
                                t_birth_s=t_cur,
                                is_primary=is_prim,
                            )
                        )
                        next_drop_id += 1

            if pinched_here and neck_lo >= 3:
                dV = math.pi * (h_row[neck_lo:drop_end + 1] ** 2) * dz_m
                vol = float(np.sum(dV))
                if vol > 0.0:
                    v_d = float(np.sum(dV * v_row[neck_lo:drop_end + 1]) / vol)
                    z_d = float(np.sum(dV * (np.arange(neck_lo, drop_end + 1) * dz_m)) / vol)
                    r_d = (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0)
                    is_prim = r_d >= prim_cutoff_m
                    ev = PinchoffEvent(
                        event_id=len(pinchoff_events) + 1,
                        row_index=r_idx,
                        t_s=t_cur,
                        neck_node=k_min,
                        neck_lo=neck_lo,
                        drop_end=drop_end,
                        z_pinch_m=k_min * dz_m,
                        z_centroid_m=z_d,
                        r_m=r_d,
                        v_ms=v_d,
                        vol_m3=vol,
                        is_primary=is_prim,
                    )
                    pinchoff_events.append(ev)
                    live_droplets.append(
                        DropletState(
                            droplet_id=next_drop_id,
                            z_m=z_d,
                            r_m=r_d,
                            v_ms=v_d,
                            vol_m3=vol,
                            t_birth_s=t_cur,
                            is_primary=is_prim,
                        )
                    )
                    next_drop_id += 1
                    last_nodes[r_idx] = neck_lo - 1

            z_tip_m = last_nodes[r_idx] * dz_m
            if config.absorb_recoiling_at_tip and z_tip_m > 0.0:
                live_droplets = [
                    d for d in live_droplets
                    if not (d.v_ms < 0.0 and (d.z_m - d.r_m) <= (z_tip_m + crit_dist_m))
                ]

            if config.coalesce_droplets and len(live_droplets) > 1:
                merged = True
                while merged:
                    merged = False
                    for i in range(len(live_droplets)):
                        for j in range(i + 1, len(live_droplets)):
                            d1, d2 = live_droplets[i], live_droplets[j]
                            if not config.coalesce_primary_satellite and (d1.is_primary != d2.is_primary):
                                continue
                            if abs(d1.z_m - d2.z_m) <= (d1.r_m + d2.r_m + crit_dist_m):
                                V = d1.vol_m3 + d2.vol_m3
                                v_new = (d1.vol_m3 * d1.v_ms + d2.vol_m3 * d2.v_ms) / V
                                z_new = (d1.vol_m3 * d1.z_m + d2.vol_m3 * d2.z_m) / V
                                r_new = (3.0 * V / (4.0 * math.pi)) ** (1.0 / 3.0)
                                d1.vol_m3 = V
                                d1.v_ms = v_new
                                d1.z_m = z_new
                                d1.r_m = r_new
                                d1.is_primary = r_new >= prim_cutoff_m
                                d1.t_birth_s = min(d1.t_birth_s, d2.t_birth_s)
                                d1.coalescence_count += d2.coalescence_count + 1
                                live_droplets.pop(j)
                                merged = True
                                break
                        if merged:
                            break

            droplets_by_row[r_idx] = [dataclasses.replace(d) for d in live_droplets]

        return cls(
            csv_path=csv_path,
            times_s=times_s,
            h_grid_m=h_grid_m,
            v_grid_ms=v_grid_ms,
            last_nodes=last_nodes,
            pinchoff_events=pinchoff_events,
            droplets_by_row=droplets_by_row,
            config=config,
        )

    @property
    def primary_events(self) -> List[PinchoffEvent]:
        return [e for e in self.pinchoff_events if e.is_primary]

    @property
    def satellite_events(self) -> List[PinchoffEvent]:
        return [e for e in self.pinchoff_events if not e.is_primary]

    def first_droplet_leave_time_us(self) -> float:
        prim = self.primary_events
        if not prim:
            return float(self.times_us[-1] * 0.35)
        first = prim[0]
        return (first.t_s + (self.config.domain_length_m - first.z_centroid_m) / first.v_ms) * 1e6

    def closest_row_index(self, t_us: float) -> int:
        return int(np.argmin(np.abs(self.times_us - t_us)))

    def unique_coalesced_droplets(self) -> Tuple[List[Tuple[float, float, float]], List[Tuple[float, float, float]]]:
        final_by_id: Dict[int, DropletState] = {}
        for row_drops in self.droplets_by_row:
            for d in row_drops:
                if d.r_um >= self.config.min_render_radius_um and d.v_ms > 0.05:
                    final_by_id[d.droplet_id] = d

        prim_list = []
        sat_list = []
        for d in sorted(final_by_id.values(), key=lambda x: x.t_birth_s):
            item = (d.t_birth_s * 1e6, d.v_ms, d.r_um)
            if d.is_primary:
                prim_list.append(item)
            else:
                sat_list.append(item)
        return prim_list, sat_list

    def summary_metrics(self) -> Dict[str, float]:
        prim = self.primary_events
        _, sat_list = self.unique_coalesced_droplets()

        prim_v = [e.v_ms for e in prim]
        prim_r = [e.r_um for e in prim]
        prim_lb = [e.z_pinch_um for e in prim]
        periods = np.diff([e.t_us for e in prim]) if len(prim) > 1 else [0.0]

        sat_velocities = [x[1] for x in sat_list]
        sat_radii_um = [x[2] for x in sat_list]

        return {
            "n_primary": len(prim),
            "n_satellite": len(sat_list),
            "v_prim_mean": float(np.mean(prim_v)) if prim_v else 0.0,
            "v_prim_std": float(np.std(prim_v)) if prim_v else 0.0,
            "r_prim_mean_um": float(np.mean(prim_r)) if prim_r else 0.0,
            "lb_mean_um": float(np.mean(prim_lb)) if prim_lb else 0.0,
            "period_mean_us": float(np.mean(periods)) if len(periods) > 0 else 0.0,
            "v_sat_mean": float(np.mean(sat_velocities)) if sat_velocities else 0.0,
            "v_sat_std": float(np.std(sat_velocities)) if sat_velocities else 0.0,
            "r_sat_mean_um": float(np.mean(sat_radii_um)) if sat_radii_um else 0.0,
            "t_leave_us": self.first_droplet_leave_time_us(),
        }


def find_latest_csv(prefix: str = "microjet") -> str:
    files = sorted(glob.glob(f"{prefix}*.csv"), key=os.path.getmtime)
    valid_files = [
        f for f in files
        if not os.path.basename(f).startswith("scratch_")
        and not os.path.basename(f).endswith("_5th_frame.csv")
    ]
    if not valid_files:
        raise FileNotFoundError(f"No {prefix}*.csv files found in current directory.")
    return valid_files[-1]


# =============================================================================
# 3. Typographic & Technical Styling Helpers
# =============================================================================

def _apply_cad_axes_style(ax: plt.Axes, style: StyleConfig, show_grid: bool = True):
    """
    Configures axes with crisp borders, inward ticks, and numerical tick labels
    rendered exclusively in `style.mono_font`.
    """
    ax.set_facecolor(style.panel_bg_color)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(style.axis_color)
        spine.set_linewidth(style.lw_border)

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        color=style.axis_color,
        labelcolor=style.text_color,
        width=style.lw_border * 0.85,
        length=3.8,
        labelsize=8.2,
        pad=4.5,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        color=style.axis_color,
        width=style.lw_border * 0.6,
        length=2.0,
    )
    if show_grid:
        ax.grid(True, which="major", color=style.grid_color, linestyle="-", linewidth=style.lw_grid, zorder=0)
    else:
        ax.grid(False)

    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontfamily(style.mono_font)


def _set_figure_title(ax: plt.Axes, title_text: str, config: PlotConfig):
    """
    Renders an understated, sentence-case heading in `style.title_font`
    (editorial serif, normal weight) when `config.show_titles` is enabled.
    """
    if not config.show_titles or not title_text:
        return
    style = config.style
    ax.set_title(
        title_text,
        loc="left",
        pad=10,
        fontsize=11.5,
        fontfamily=style.title_font,
        fontweight="normal",
        color=style.text_color,
    )


def _draw_cad_dimension_bar(
    ax: plt.Axes,
    x_start: float,
    x_end: float,
    y_pos: float,
    label: str,
    style: StyleConfig,
    tick_height: float = 18.0,
):
    """Draws a crisp dimension bar |--- 200 µm ---| with monospace numerical value."""
    c = style.dimension_color
    lw = style.lw_dimension
    ax.plot([x_start, x_end], [y_pos, y_pos], color=c, linewidth=lw, zorder=10, solid_capstyle="butt")
    ax.plot([x_start, x_start], [y_pos - tick_height * 0.5, y_pos + tick_height * 0.5], color=c, linewidth=lw, zorder=10)
    ax.plot([x_end, x_end], [y_pos - tick_height * 0.5, y_pos + tick_height * 0.5], color=c, linewidth=lw, zorder=10)
    ax.text(
        0.5 * (x_start + x_end),
        y_pos + tick_height * 0.65,
        label,
        color=style.text_color,
        fontsize=7.4,
        fontfamily=style.mono_font,
        ha="center",
        va="bottom",
        bbox=dict(boxstyle="square,pad=0.10", facecolor=style.panel_bg_color, edgecolor="none", alpha=0.90),
        zorder=11,
    )


def _save_standalone(fig: plt.Figure, output_path: str, config: PlotConfig):
    """Saves a standalone figure with generous outer whitespace."""
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(
        output_path,
        dpi=config.dpi,
        facecolor=config.style.bg_color,
        bbox_inches="tight",
        pad_inches=config.pad_inches,
    )
    plt.close(fig)


# =============================================================================
# 4. Stroboscopic Montage Figures (01 Transient & 02 Steady-State)
# =============================================================================

def render_stroboscopic_panel(
    ax: plt.Axes,
    dataset: MicrojetDataset,
    config: PlotConfig,
    t_start_us: float,
    t_end_us: float,
    title_text: str = "",
    tick_step_us: float = 20.0,
    velocity_norm: Optional[Normalize] = None,
    velocity_cmap: Optional[matplotlib.colors.Colormap] = None,
    show_right_axis: bool = True,
):
    """
    Renders a stroboscopic time-series montage of the microjet and detached droplets
    onto `ax` over `[t_start_us, t_end_us]` with generous vertical & horizontal whitespace.
    """
    style = config.style
    scale_x = config.scale_x
    domain_um = config.domain_length_m * 1e6

    sample_times: List[float] = []
    t_val = math.ceil(t_start_us / config.time_step_us) * config.time_step_us
    while t_val <= t_end_us + 1e-6:
        if t_val >= dataset.times_us[0] - 1e-6 and t_val <= dataset.times_us[-1] + 1e-6:
            sample_times.append(t_val)
        t_val += config.time_step_us

    jet_polys = []
    jet_edge_polys = []
    vel_quads = []
    vel_values = []

    primary_circles = []
    satellite_circles = []
    centerline_segments = []

    num_dome_pts = 14
    dome_angles = np.linspace(0.0, math.pi, num_dome_pts)
    num_circ_pts = 36
    circ_angles = np.linspace(0.0, 2.0 * math.pi, num_circ_pts, endpoint=False)
    cos_c = np.cos(circ_angles)
    sin_c = np.sin(circ_angles)

    for t_sample in sample_times:
        r_idx = dataset.closest_row_index(t_sample)
        last_i = dataset.last_nodes[r_idx]
        if last_i < 1:
            continue

        h_um = dataset.h_grid_m[r_idx, :last_i + 1] * 1e6
        v_ms = dataset.v_grid_ms[r_idx, :last_i + 1]
        z_um = dataset.z_grid_um[:last_i + 1]

        x_right = t_sample + h_um * scale_x
        y_right = z_um

        h_tip = h_um[-1]
        z_tip = z_um[-1]
        x_dome = t_sample + np.cos(dome_angles[1:-1]) * h_tip * scale_x
        y_dome = z_tip + np.sin(dome_angles[1:-1]) * h_tip

        x_left = t_sample - h_um[::-1] * scale_x
        y_left = z_um[::-1]

        x_poly = np.concatenate([x_right, x_dome, x_left, [x_right[0]]])
        y_poly = np.concatenate([y_right, y_dome, y_left, [y_right[0]]])
        poly_verts = np.column_stack([x_poly, y_poly])
        jet_edge_polys.append(poly_verts)

        if config.color_by_velocity:
            for i in range(last_i):
                quad = [
                    (x_left[last_i - i], z_um[i]),
                    (x_right[i], z_um[i]),
                    (x_right[i + 1], z_um[i + 1]),
                    (x_left[last_i - 1 - i], z_um[i + 1]),
                ]
                vel_quads.append(quad)
                vel_values.append(0.5 * (v_ms[i] + v_ms[i + 1]))
            dome_verts = np.column_stack([
                np.concatenate([[x_right[-1]], x_dome, [x_left[0]]]),
                np.concatenate([[y_right[-1]], y_dome, [y_left[0]]]),
            ])
            vel_quads.append(dome_verts)
            vel_values.append(v_ms[-1])
        else:
            jet_polys.append(poly_verts)

        if config.show_centerlines:
            centerline_segments.append([(t_sample, 0.0), (t_sample, z_tip + h_tip)])

        for d in dataset.droplets_by_row[r_idx]:
            if d.r_um < config.min_render_radius_um:
                continue
            if d.z_um - d.r_um <= domain_um:
                cx = t_sample + (d.r_um * scale_x) * cos_c
                cy = d.z_um + d.r_um * sin_c
                c_verts = np.column_stack([cx, cy])
                if d.is_primary or not config.highlight_satellites:
                    primary_circles.append(c_verts)
                else:
                    satellite_circles.append(c_verts)

    # Subtle nozzle datum line at z = 0
    ax.axhline(0.0, color=style.grid_color, linestyle="-", linewidth=style.lw_grid, zorder=1)

    if config.color_by_velocity and vel_quads:
        if velocity_norm is None:
            velocity_norm = Normalize(vmin=-1.5, vmax=9.5)
        if velocity_cmap is None:
            velocity_cmap = LinearSegmentedColormap.from_list("cad_vel", style.velocity_cmap_colors, N=256)
        vel_coll = PolyCollection(
            vel_quads,
            array=np.array(vel_values),
            cmap=velocity_cmap,
            norm=velocity_norm,
            edgecolors="face",
            linewidths=0.15,
            zorder=2,
        )
        ax.add_collection(vel_coll)
        contour_coll = PolyCollection(
            jet_edge_polys,
            facecolors="none",
            edgecolors=style.jet_edge_color,
            linewidths=style.lw_contour,
            zorder=3,
        )
        ax.add_collection(contour_coll)
    elif jet_polys:
        jet_coll = PolyCollection(
            jet_polys,
            facecolors=style.jet_fill_color,
            edgecolors=style.jet_edge_color,
            linewidths=style.lw_contour,
            zorder=2,
        )
        ax.add_collection(jet_coll)

    if config.show_centerlines and centerline_segments:
        cl_coll = LineCollection(
            centerline_segments,
            colors=style.centerline_color,
            linestyles="-.",
            linewidths=style.lw_centerline,
            zorder=3.5,
        )
        ax.add_collection(cl_coll)

    if primary_circles:
        prim_coll = PolyCollection(
            primary_circles,
            facecolors=style.primary_fill_color,
            edgecolors=style.primary_edge_color,
            linewidths=style.lw_contour,
            zorder=4,
        )
        ax.add_collection(prim_coll)

    if satellite_circles:
        sat_coll = PolyCollection(
            satellite_circles,
            facecolors=style.satellite_fill_color,
            edgecolors=style.satellite_edge_color,
            linewidths=style.lw_contour * 1.15,
            zorder=5,
        )
        ax.add_collection(sat_coll)

    # Mean breakup length datum line
    if config.show_breakup_locus and dataset.primary_events:
        mean_lb = float(np.mean([e.z_pinch_um for e in dataset.primary_events]))
        ax.axhline(
            mean_lb,
            color=style.locus_color,
            linestyle="--",
            linewidth=style.lw_datum,
            alpha=0.65,
            zorder=3.2,
        )
        if config.show_dimension_callouts:
            ax.text(
                t_end_us + 2.3,
                mean_lb + 12.0,
                f"{mean_lb:.0f} µm",
                color=style.locus_color,
                fontsize=7.2,
                fontfamily=style.mono_font,
                ha="left",
                va="bottom",
                bbox=dict(boxstyle="square,pad=0.08", facecolor=style.panel_bg_color, edgecolor="none", alpha=0.90),
                zorder=10,
            )

    # Optional velocity slope trajectory callout
    if config.show_trajectories:
        prim_in_win = [e for e in dataset.primary_events if t_start_us + 2.0 <= e.t_us <= t_end_us - 25.0]
        if prim_in_win:
            ref_ev = prim_in_win[len(prim_in_win) // 2]
            t0 = ref_ev.t_us
            z0 = ref_ev.z_centroid_m * 1e6
            v_p = ref_ev.v_ms
            dt_line = min(40.0, (domain_um - z0) / v_p, t_end_us - t0 - 2.0)
            if dt_line > 10.0:
                t1 = t0 + dt_line
                z1 = z0 + v_p * dt_line
                ax.plot(
                    [t0, t1], [z0, z1],
                    color=style.traj_primary_color,
                    linestyle="--",
                    linewidth=style.lw_datum * 1.1,
                    zorder=7,
                )

    x_pad_left = 3.5
    x_pad_right = 10.2

    if config.show_scale_bar:
        bar_len_um = 200.0
        bar_len_x = bar_len_um * scale_x
        bar_xstart = t_end_us + 2.3
        bar_xend = bar_xstart + bar_len_x
        _draw_cad_dimension_bar(ax, bar_xstart, bar_xend, 60.0, "200 µm", style, tick_height=18.0)

    _set_figure_title(ax, title_text, config)

    # Generous horizontal and vertical breathing room
    ax.set_xlim(t_start_us - x_pad_left, t_end_us + x_pad_right)
    ax.set_ylim(config.z_min_um, config.z_max_um)
    ax.set_aspect(scale_x, adjustable="box")

    first_tick = math.ceil(t_start_us / tick_step_us) * tick_step_us
    xticks = np.arange(first_tick, t_end_us + 1e-6, tick_step_us)
    ax.set_xticks(xticks)
    ax.set_xlabel("Time, t (µs)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)

    if config.show_z_axis:
        ax.set_yticks([0, 200, 400, 600, 800])
        ax.set_ylabel("Axial position, z (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
        if show_right_axis:
            h0_um = config.h0_m * 1e6
            sec_ax = ax.secondary_yaxis("right", functions=(lambda z: z / h0_um, lambda zh: zh * h0_um))
            sec_ax.set_yticks([0, 10, 20, 30, 40])
            sec_ax.set_ylabel("z / r₀", fontsize=8.8, fontfamily=style.body_light_font, color=style.muted_text_color, labelpad=6)
            sec_ax.tick_params(direction="in", labelsize=8.0, colors=style.muted_text_color, length=3.2, pad=4.0)
            sec_ax.spines["right"].set_color(style.axis_color)
            sec_ax.spines["right"].set_linewidth(style.lw_border)
            for lbl in sec_ax.get_yticklabels():
                lbl.set_fontfamily(style.mono_font)
    else:
        ax.set_yticks([])

    _apply_cad_axes_style(ax, style, show_grid=config.show_montage_grid)


def plot_transient_montage(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_path: str = "01_transient_montage.png",
) -> str:
    """
    Standalone Figure 1: Transient startup regime (0 <= t <= t_leave).
    """
    if config is None:
        config = dataset.config
    style = config.style

    t_leave_us = dataset.first_droplet_leave_time_us()
    t_min_us = config.t_min_us if config.t_min_us is not None else float(dataset.times_us[0])
    t_end_before = math.ceil(t_leave_us / 10.0) * 10.0

    span_us = max(60.0, t_end_before - t_min_us)
    fig_w = max(8.8, (span_us / 150.0) * 10.2)
    fig, ax = plt.subplots(figsize=(fig_w, 3.8), dpi=config.dpi, facecolor=style.bg_color)

    render_stroboscopic_panel(
        ax, dataset, config,
        t_start_us=t_min_us,
        t_end_us=t_end_before,
        title_text="Transient jet emergence and initial breakup",
        tick_step_us=20.0,
    )
    _save_standalone(fig, output_path, config)
    return output_path


def plot_steady_montage(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_path: str = "02_steady_montage.png",
) -> str:
    """
    Standalone Figure 2: Steady-state periodic Rayleigh breakup over a focused,
    uncrowded time window (`config.steady_span_us`, default 80 µs).
    """
    if config is None:
        config = dataset.config
    style = config.style

    t_leave_us = dataset.first_droplet_leave_time_us()
    t_max_sim = config.t_max_us if config.t_max_us is not None else float(dataset.times_us[-1])

    if config.steady_start_us is not None:
        t_start = config.steady_start_us
    else:
        t_start = math.ceil(t_leave_us / 10.0) * 10.0

    t_end = min(t_max_sim, t_start + config.steady_span_us)
    if t_end <= t_start:
        t_start = max(0.0, t_max_sim - config.steady_span_us)
        t_end = t_max_sim

    span_us = max(50.0, t_end - t_start)
    fig_w = max(6.8, (span_us / 80.0) * 7.2)
    fig, ax = plt.subplots(figsize=(fig_w, 3.8), dpi=config.dpi, facecolor=style.bg_color)

    render_stroboscopic_panel(
        ax, dataset, config,
        t_start_us=t_start,
        t_end_us=t_end,
        title_text="Steady-state periodic Rayleigh breakup",
        tick_step_us=10.0 if span_us <= 120.0 else 20.0,
    )
    _save_standalone(fig, output_path, config)
    return output_path


# =============================================================================
# 5. Standalone Quantitative Figure: Neck Thinning Cycles & Breakup Length
# =============================================================================

def plot_neck_cycles_and_breakup_length(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_path: str = "03_neck_cycles_and_breakup_length.png",
) -> str:
    """
    Standalone Figure 3: Minimum neck radius across periodic pinch-off
    cycles and primary breakup length.
    """
    if config is None:
        config = dataset.config
    style = config.style

    fig, ax_evol = plt.subplots(figsize=(7.8, 4.2), dpi=config.dpi, facecolor=style.bg_color)
    _set_figure_title(ax_evol, "Periodic neck thinning cycles and breakup length", config)

    h_min_series_um = []
    for r_i in range(dataset.num_rows):
        ln = dataset.last_nodes[r_i]
        if ln >= 2:
            h_min_series_um.append(float(np.min(dataset.h_grid_m[r_i, 1:ln + 1])) * 1e6)
        else:
            h_min_series_um.append(config.h0_m * 1e6)

    ax_evol.plot(
        dataset.times_us,
        h_min_series_um,
        color=style.jet_edge_color,
        linewidth=0.95,
        label="Minimum neck radius",
        zorder=3,
    )
    ax_evol.axhline(
        config.pinch_threshold_m * 1e6,
        color=style.locus_color,
        linestyle=":",
        linewidth=0.90,
        label=f"Pinch-off threshold ({config.pinch_threshold_m*1e6:.1f} µm)",
        zorder=2,
    )

    ax_lb = ax_evol.twinx()
    prim_ev = dataset.primary_events
    if prim_ev:
        ax_lb.plot(
            [e.t_us for e in prim_ev],
            [e.z_pinch_um for e in prim_ev],
            color=style.datum_color,
            marker="o",
            markersize=4.0,
            linestyle="--",
            linewidth=1.0,
            label="Primary breakup length",
            zorder=5,
        )
    ax_lb.set_ylabel("Breakup length (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.datum_color, labelpad=6)
    ax_lb.tick_params(direction="in", labelsize=8.2, colors=style.datum_color, length=3.8, pad=4.5)
    ax_lb.set_ylim(240.0, 480.0)
    for lbl in ax_lb.get_yticklabels():
        lbl.set_fontfamily(style.mono_font)

    h0_um = config.h0_m * 1e6
    ax_evol.set_xlim(float(dataset.times_us[0]), float(dataset.times_us[-1]))
    ax_evol.set_ylim(0.0, h0_um * 1.28)
    ax_evol.set_xlabel("Time, t (µs)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
    ax_evol.set_ylabel("Minimum neck radius (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
    _apply_cad_axes_style(ax_evol, style, show_grid=config.show_grid)

    l1, lb1 = ax_evol.get_legend_handles_labels()
    l2, lb2 = ax_lb.get_legend_handles_labels()
    leg = ax_evol.legend(
        l1 + l2, lb1 + lb2,
        loc="upper right",
        ncol=1,
        frameon=True,
        facecolor=style.panel_bg_color,
        edgecolor=style.axis_color,
        prop={"family": style.body_font, "size": 7.6},
    )
    leg.get_frame().set_linewidth(0.5)

    _save_standalone(fig, output_path, config)
    return output_path


# =============================================================================
# 6. Standalone Pinch-Off Detail Figures (04 Sequence, 05 Neck Stages, 06 Profile)
# =============================================================================

def plot_pinchoff_sequence(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_path: str = "04_pinchoff_sequence.png",
    event_index: int = 2,
    offsets_us: Sequence[float] = (-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0),
) -> str:
    """
    Standalone Figure 4: Synchronized cross-sectional pinch-off sequence
    (Δt = -4 to +10 µs around primary breakup, cf. Fig. 5 of Van Hoeve et al., 2010).
    """
    if config is None:
        config = dataset.config
    style = config.style

    prim_events = dataset.primary_events
    if not prim_events:
        raise ValueError("No primary pinch-off events found in dataset.")
    ref_event = prim_events[min(event_index, len(prim_events) - 1)]
    t_ref_us = ref_event.t_us
    z_pinch_um = ref_event.z_pinch_um

    n_frames = len(offsets_us)
    fig = plt.figure(figsize=(8.0, 8.2), dpi=config.dpi, facecolor=style.bg_color)
    gs = gridspec.GridSpec(n_frames, 1, figure=fig, hspace=0.22)

    vel_norm = Normalize(vmin=-1.5, vmax=9.5)
    vel_cmap = LinearSegmentedColormap.from_list("cad_vel", style.velocity_cmap_colors, N=256)

    z_win_min = max(0.0, z_pinch_um - 180.0)
    z_win_max = min(config.domain_length_m * 1e6, z_pinch_um + 220.0)
    r_lim = 58.0

    circ_angles = np.linspace(0.0, 2.0 * math.pi, 48, endpoint=False)
    cos_c = np.cos(circ_angles)
    sin_c = np.sin(circ_angles)
    dome_angles = np.linspace(-0.5 * math.pi, 0.5 * math.pi, 16)

    for idx_f, dt_us in enumerate(offsets_us):
        ax = fig.add_subplot(gs[idx_f, 0])

        if idx_f == 0:
            _set_figure_title(ax, "Synchronized pinch-off sequence", config)

        t_target = t_ref_us + dt_us
        r_idx = dataset.closest_row_index(t_target)
        last_i = dataset.last_nodes[r_idx]

        ax.axhline(0.0, color=style.centerline_color, linestyle="-.", linewidth=style.lw_centerline, zorder=1)
        ax.axvline(z_pinch_um, color=style.locus_color, linestyle=":", linewidth=style.lw_datum * 0.85, alpha=0.65, zorder=1)

        if last_i >= 1:
            z_um = dataset.z_grid_um[:last_i + 1]
            h_um = dataset.h_grid_m[r_idx, :last_i + 1] * 1e6
            v_ms = dataset.v_grid_ms[r_idx, :last_i + 1]

            quads = []
            q_vals = []
            for i in range(last_i):
                if z_um[i + 1] < z_win_min or z_um[i] > z_win_max:
                    continue
                quads.append([
                    (z_um[i], -h_um[i]),
                    (z_um[i + 1], -h_um[i + 1]),
                    (z_um[i + 1], h_um[i + 1]),
                    (z_um[i], h_um[i]),
                ])
                q_vals.append(0.5 * (v_ms[i] + v_ms[i + 1]))

            h_tip = h_um[-1]
            z_tip = z_um[-1]
            z_dome = z_tip + np.cos(dome_angles) * h_tip
            r_dome = np.sin(dome_angles) * h_tip
            dome_verts = np.column_stack([z_dome, r_dome])
            quads.append(dome_verts)
            q_vals.append(v_ms[-1])

            if quads:
                vc = PolyCollection(quads, array=np.array(q_vals), cmap=vel_cmap, norm=vel_norm, edgecolors="face", linewidths=0.1, zorder=2)
                ax.add_collection(vc)

            z_contour = np.concatenate([z_um, z_dome[1:-1], z_um[::-1]])
            r_contour = np.concatenate([h_um, r_dome[1:-1], -h_um[::-1]])
            ax.plot(z_contour, r_contour, color=style.jet_edge_color, linewidth=style.lw_contour * 1.3, zorder=4)

        for d in dataset.droplets_by_row[r_idx]:
            if d.r_um < config.min_render_radius_um:
                continue
            if z_win_min - 40.0 <= d.z_um <= z_win_max + 40.0:
                dz_c = d.z_um + d.r_um * cos_c
                dr_c = d.r_um * sin_c
                fill_c = vel_cmap(vel_norm(d.v_ms)) if d.is_primary else style.satellite_fill_color
                edge_c = style.satellite_edge_color if not d.is_primary else style.jet_edge_color
                ax.fill(dz_c, dr_c, facecolor=fill_c, edgecolor=edge_c, linewidth=style.lw_contour * 1.3, zorder=5)

        sign_str = "+" if dt_us > 0 else ""
        badge_txt = f"t = {sign_str}{dt_us:.0f} µs"
        ax.text(
            z_win_min + 6.0, r_lim - 6.5, badge_txt,
            fontsize=7.6, fontfamily=[style.mono_font, style.body_font],
            color=style.text_color, ha="left", va="top",
            bbox=dict(boxstyle="square,pad=0.10", facecolor=style.panel_bg_color, edgecolor="none", alpha=0.88),
            zorder=10,
        )

        ax.set_xlim(z_win_min, z_win_max)
        ax.set_ylim(-r_lim, r_lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_yticks([-30, 0, 30])
        if idx_f == n_frames // 2:
            ax.set_ylabel("Radius, r (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
        if idx_f < n_frames - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Axial position, z (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
            _draw_cad_dimension_bar(ax, z_win_min + 12.0, z_win_min + 112.0, -45.0, "100 µm", style, tick_height=10.0)

        _apply_cad_axes_style(ax, style, show_grid=config.show_grid)

    _save_standalone(fig, output_path, config)
    return output_path


def plot_neck_collapse_stages(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_path: str = "05_neck_collapse_stages.png",
    event_index: int = 2,
) -> str:
    """
    Standalone Figure 5: Three-stage capillary neck collapse approaching breakup
    (2, 1, 0 µs prior to pinch-off, cf. Fig. 7 of Van Hoeve et al., 2010).
    """
    if config is None:
        config = dataset.config
    style = config.style

    prim_events = dataset.primary_events
    if not prim_events:
        raise ValueError("No primary pinch-off events found in dataset.")
    ref_event = prim_events[min(event_index, len(prim_events) - 1)]
    t_ref_us = ref_event.t_us
    z_pinch_um = ref_event.z_pinch_um

    fig, ax_neck = plt.subplots(figsize=(7.4, 4.8), dpi=config.dpi, facecolor=style.bg_color)
    _set_figure_title(ax_neck, "Capillary neck collapse stages approaching breakup", config)

    stage_offsets = [(-2.0, "t = −2 µs"), (-1.0, "t = −1 µs"), (0.0, "t = 0 µs")]
    y_shifts = [88.0, 0.0, -88.0]
    neck_z_min = z_pinch_um - 95.0
    neck_z_max = z_pinch_um + 80.0

    for (dt_s, stage_lbl), y0 in zip(stage_offsets, y_shifts):
        r_idx = dataset.closest_row_index(t_ref_us + dt_s)
        h_um = dataset.h_grid_m[r_idx] * 1e6
        z_um = dataset.z_grid_um
        mask = (z_um >= neck_z_min) & (z_um <= neck_z_max) & (h_um > 1e-3)
        z_seg = z_um[mask]
        h_seg = h_um[mask]

        ax_neck.axhline(y0, color=style.centerline_color, linestyle="-.", linewidth=style.lw_centerline, zorder=1)
        if len(z_seg) > 2:
            ax_neck.fill_between(
                z_seg, y0 - h_seg, y0 + h_seg,
                facecolor=style.jet_fill_color, edgecolor=style.jet_edge_color,
                linewidth=style.lw_contour * 1.25, zorder=3,
            )
            k_loc = int(np.argmin(h_seg))
            z_min_loc = z_seg[k_loc]
            h_min_loc = h_seg[k_loc]
            ax_neck.plot(
                [z_min_loc, z_min_loc],
                [y0 - h_min_loc - 10.0, y0 + h_min_loc + 10.0],
                color=style.locus_color,
                linewidth=0.8,
                zorder=6,
            )
            ax_neck.text(
                z_min_loc, y0 + h_min_loc + 12.0,
                f"2h = {2.0*h_min_loc:.2f} µm",
                fontsize=7.4, fontfamily=[style.mono_font, style.body_font],
                color=style.locus_color, ha="center", va="bottom",
                bbox=dict(boxstyle="square,pad=0.10", facecolor=style.panel_bg_color, edgecolor="none", alpha=0.90),
                zorder=8,
            )

        ax_neck.text(
            neck_z_min + 4.5, y0 + 38.0, stage_lbl,
            fontsize=7.6, fontfamily=[style.mono_font, style.body_font],
            color=style.text_color, ha="left", va="top",
            bbox=dict(boxstyle="square,pad=0.12", facecolor=style.panel_bg_color, edgecolor="none", alpha=0.90),
            zorder=9,
        )

    ax_neck.set_xlim(neck_z_min, neck_z_max)
    ax_neck.set_ylim(-142.0, 142.0)
    ax_neck.set_yticks(y_shifts)
    ax_neck.set_yticklabels(["i", "ii", "iii"], fontsize=8.5, fontfamily=style.mono_font)
    ax_neck.set_ylabel("Collapse stage", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
    ax_neck.set_xlabel("Axial position, z (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
    _apply_cad_axes_style(ax_neck, style, show_grid=config.show_grid)

    _save_standalone(fig, output_path, config)
    return output_path


def plot_pinchoff_contour_and_velocity(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_path: str = "06_pinchoff_contour_and_velocity.png",
    event_index: int = 2,
) -> str:
    """
    Standalone Figure 6: Instantaneous pinch-off contour ±h(z), equivalent
    detached droplet sphere, and axial velocity spike v(z) (cf. Fig. 11).
    """
    if config is None:
        config = dataset.config
    style = config.style

    prim_events = dataset.primary_events
    if not prim_events:
        raise ValueError("No primary pinch-off events found in dataset.")
    ref_event = prim_events[min(event_index, len(prim_events) - 1)]
    z_pinch_um = ref_event.z_pinch_um

    fig, ax_prof = plt.subplots(figsize=(7.6, 4.4), dpi=config.dpi, facecolor=style.bg_color)
    _set_figure_title(ax_prof, "Pinch-off contour and axial velocity profile", config)

    r_idx_pinch = ref_event.row_index
    h_pinch_um = dataset.h_grid_m[r_idx_pinch] * 1e6
    v_pinch_ms = dataset.v_grid_ms[r_idx_pinch]
    z_rel_um = dataset.z_grid_um - z_pinch_um
    active_p = h_pinch_um > 1e-3

    z_p = z_rel_um[active_p]
    h_p = h_pinch_um[active_p]
    v_p = v_pinch_ms[active_p]

    circ_angles = np.linspace(0.0, 2.0 * math.pi, 64, endpoint=False)
    cos_c = np.cos(circ_angles)
    sin_c = np.sin(circ_angles)

    ax_prof.fill_between(z_p, -h_p, h_p, facecolor=style.jet_fill_color, alpha=0.55, zorder=2)
    ax_prof.plot(z_p, h_p, color=style.jet_edge_color, linewidth=1.15, label="Lubrication contour, ±h(z)", zorder=4)
    ax_prof.plot(z_p, -h_p, color=style.jet_edge_color, linewidth=1.15, zorder=4)
    ax_prof.axhline(0.0, color=style.centerline_color, linestyle="-.", linewidth=style.lw_centerline, zorder=1)
    ax_prof.axvline(0.0, color=style.locus_color, linestyle=":", linewidth=style.lw_datum, label="Pinch-off plane", zorder=3)

    z_drop_rel = (ref_event.z_centroid_m * 1e6) - z_pinch_um
    r_drop_um = ref_event.r_um
    ax_prof.plot(
        z_drop_rel + r_drop_um * cos_c,
        r_drop_um * sin_c,
        color=style.datum_color,
        linestyle="--",
        linewidth=1.05,
        label=f"Equivalent sphere (R = {r_drop_um:.1f} µm)",
        zorder=5,
    )

    ax_vel = ax_prof.twinx()
    ax_vel.plot(z_p, v_p, color=style.locus_color, linewidth=1.05, linestyle="-", label="Axial velocity, v(z)", zorder=6)
    ax_vel.set_ylabel("Axial velocity, v (m/s)", fontsize=9.2, fontfamily=style.body_font, color=style.locus_color, labelpad=6)
    ax_vel.tick_params(direction="in", labelsize=8.2, colors=style.locus_color, length=3.8, pad=4.5)
    ax_vel.set_ylim(-2.5, max(12.5, float(np.max(v_p)) * 1.18))
    for lbl in ax_vel.get_yticklabels():
        lbl.set_fontfamily(style.mono_font)

    ax_prof.set_xlim(-195.0, 95.0)
    ax_prof.set_ylim(-56.0, 56.0)
    ax_prof.set_xlabel("Axial distance from pinch-off (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
    ax_prof.set_ylabel("Radius, r (µm)", fontsize=9.2, fontfamily=style.body_font, color=style.text_color, labelpad=6)
    _apply_cad_axes_style(ax_prof, style, show_grid=config.show_grid)

    lines1, labels1 = ax_prof.get_legend_handles_labels()
    lines2, labels2 = ax_vel.get_legend_handles_labels()
    leg = ax_prof.legend(
        lines1 + lines2, labels1 + labels2,
        loc="lower left",
        frameon=True,
        facecolor=style.panel_bg_color,
        edgecolor=style.axis_color,
        prop={"family": style.body_font, "size": 7.5},
    )
    leg.get_frame().set_linewidth(0.5)

    _save_standalone(fig, output_path, config)
    return output_path


# =============================================================================
# 7. Per-Run Output Directory & Pipeline Entrypoint
# =============================================================================

def create_run_output_dir(csv_path: str, output_root: str = "plots", custom_dir: Optional[str] = None) -> str:
    """
    Creates a dedicated output folder for this run of the script, e.g.:
    `plots/microjet10-15-2026_20261015_215500/` and updates `plots/latest/`.
    """
    if custom_dir:
        os.makedirs(custom_dir, exist_ok=True)
        return custom_dir

    stem = os.path.splitext(os.path.basename(csv_path))[0]
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(output_root, f"{stem}_{stamp}")
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def generate_all_figures(
    dataset: MicrojetDataset,
    config: Optional[PlotConfig] = None,
    output_dir: Optional[str] = None,
    ext: str = "png",
) -> Dict[str, str]:
    """
    Generates all 6 standalone publication figures inside a dedicated per-run
    output directory and returns a dictionary mapping figure keys to file paths.
    """
    if config is None:
        config = dataset.config

    run_dir = create_run_output_dir(dataset.csv_path, config.output_root, output_dir)

    paths = {
        "transient_montage": os.path.join(run_dir, f"01_transient_montage.{ext}"),
        "steady_montage": os.path.join(run_dir, f"02_steady_montage.{ext}"),
        "neck_cycles": os.path.join(run_dir, f"03_neck_cycles_and_breakup_length.{ext}"),
        "pinchoff_sequence": os.path.join(run_dir, f"04_pinchoff_sequence.{ext}"),
        "neck_collapse_stages": os.path.join(run_dir, f"05_neck_collapse_stages.{ext}"),
        "pinchoff_contour_velocity": os.path.join(run_dir, f"06_pinchoff_contour_and_velocity.{ext}"),
    }

    plot_transient_montage(dataset, config, output_path=paths["transient_montage"])
    plot_steady_montage(dataset, config, output_path=paths["steady_montage"])
    plot_neck_cycles_and_breakup_length(dataset, config, output_path=paths["neck_cycles"])
    plot_pinchoff_sequence(dataset, config, output_path=paths["pinchoff_sequence"])
    plot_neck_collapse_stages(dataset, config, output_path=paths["neck_collapse_stages"])
    plot_pinchoff_contour_and_velocity(dataset, config, output_path=paths["pinchoff_contour_velocity"])

    # Maintain `plots/latest/` mirror for quick access
    latest_dir = os.path.join(config.output_root, "latest")
    if os.path.abspath(run_dir) != os.path.abspath(latest_dir):
        if os.path.islink(latest_dir) or os.path.isfile(latest_dir):
            os.unlink(latest_dir)
        elif os.path.isdir(latest_dir):
            shutil.rmtree(latest_dir)
        shutil.copytree(run_dir, latest_dir)

    # Also update top-level convenience files if enabled
    if config.also_save_legacy_names and ext == "png":
        shutil.copy2(paths["transient_montage"], "microjet_before_leave.png")
        shutil.copy2(paths["steady_montage"], "microjet_after_leave.png")

    return paths


def ingest_and_plot(
    csv_file: Optional[str] = None,
    output_dir: Optional[str] = None,
    time_step_us: float = 2.0,
    t_max_us: Optional[float] = None,
    steady_span_us: float = 80.0,
    config: Optional[PlotConfig] = None,
    ext: str = "png",
) -> Tuple[MicrojetDataset, Dict[str, str]]:
    """
    Main pipeline entry point: loads the simulation CSV, reconstructs droplet & coalescence
    physics, and writes all standalone figures into a per-run output directory.
    """
    if config is None:
        config = PlotConfig(time_step_us=time_step_us, t_max_us=t_max_us, steady_span_us=steady_span_us)
    else:
        config.time_step_us = time_step_us
        config.steady_span_us = steady_span_us
        if t_max_us is not None:
            config.t_max_us = t_max_us

    dataset = MicrojetDataset.from_csv(csv_file, config)
    m = dataset.summary_metrics()

    print(f"Ingested: {dataset.csv_path}")
    print(f"  Fonts: title='{config.style.title_font}', body='{config.style.body_font}', mono='{config.style.mono_font}'")
    print(f"  Grid: {dataset.num_nodes} nodes | dz = {dataset.dz_m*1e6:.3f} µm | {dataset.num_rows} slices (t_max = {dataset.times_us[-1]:.1f} µs)")
    print(f"  Dimensionless: We = {config.weber:.2f}, Oh = {config.ohnesorge:.3f}, Re = {config.reynolds:.1f}, τc = {config.tau_c_us:.2f} µs")
    print(f"  Primary droplets: {m['n_primary']} formed | <R> = {m['r_prim_mean_um']:.2f} µm (Rayleigh: {config.rayleigh_primary_radius_um:.2f} µm) | <v> = {m['v_prim_mean']:.2f} +/- {m['v_prim_std']:.2f} m/s")
    print(f"  Satellite droplets: {m['n_satellite']} formed | <R> = {m['r_sat_mean_um']:.2f} µm | <v> = {m['v_sat_mean']:.2f} +/- {m['v_sat_std']:.2f} m/s")
    print(f"  First primary droplet exits domain (z = {config.domain_length_m*1e6:.0f} µm) at t = {m['t_leave_us']:.2f} µs")

    paths = generate_all_figures(dataset, config, output_dir=output_dir, ext=ext)
    run_folder = os.path.dirname(next(iter(paths.values())))
    print(f"  -> Saved {len(paths)} standalone figures in: {run_folder}/")
    for k, p in paths.items():
        print(f"     - {os.path.basename(p)}")

    return dataset, paths


def _build_cli_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Standalone Scientific Visualization Toolchain for Rayleigh Microjet Breakup."
    )
    p.add_argument("csv_file", nargs="?", default=None, help="Path to simulation CSV (defaults to latest microjet*.csv)")
    p.add_argument("--out-dir", type=str, default=None, help="Explicit output directory (defaults to plots/<csv_stem>_<timestamp>/)")
    p.add_argument("--theme", choices=["cad_light", "cad_blueprint", "paper_classic"], default="cad_light", help="Visual theme preset")
    p.add_argument("--dt", type=float, default=2.0, help="Stroboscopic sample spacing in µs (default: 2.0)")
    p.add_argument("--t-min", type=float, default=None, help="Minimum time (µs) for transient montage")
    p.add_argument("--t-max", type=float, default=None, help="Maximum time (µs) for simulation dataset")
    p.add_argument("--steady-start", type=float, default=None, help="Start time (µs) for steady-state montage (default: ~150 µs)")
    p.add_argument("--steady-span", type=float, default=80.0, help="Duration (µs) of steady-state montage window (default: 80.0)")
    p.add_argument("--no-titles", action="store_true", help="Omit figure titles for direct LaTeX/report embedding")
    p.add_argument("--velocity-fill", action="store_true", help="Shade jet interior by axial velocity field v(z, t)")
    p.add_argument("--montage-grid", action="store_true", help="Show background grid on stroboscopic montages")
    p.add_argument("--centerlines", action="store_true", help="Draw axis-of-symmetry centerlines on stroboscopic montage")
    p.add_argument("--trajectories", action="store_true", help="Overlay droplet velocity trajectory slope line")
    p.add_argument("--no-locus", action="store_true", help="Hide breakup length datum line <Lb>")
    p.add_argument("--no-coalesce", action="store_true", help="Disable droplet-droplet coalescence")
    p.add_argument("--merge-primary-satellite", action="store_true", help="Allow primary droplets to absorb satellite droplets")
    p.add_argument("--dpi", type=int, default=350, help="Output image resolution in DPI (default: 350)")
    p.add_argument("--format", choices=["png", "pdf", "svg"], default="png", help="Output file format (default: png)")
    return p


if __name__ == "__main__":
    parser = _build_cli_parser()
    args = parser.parse_args()

    theme_map = {
        "cad_light": StyleConfig.cad_light,
        "cad_blueprint": StyleConfig.cad_blueprint,
        "paper_classic": StyleConfig.paper_classic,
    }
    style_cfg = theme_map[args.theme]()

    cfg = PlotConfig(
        time_step_us=args.dt,
        t_min_us=args.t_min,
        t_max_us=args.t_max,
        steady_start_us=args.steady_start,
        steady_span_us=args.steady_span,
        show_titles=not args.no_titles,
        color_by_velocity=args.velocity_fill,
        show_montage_grid=args.montage_grid,
        show_centerlines=args.centerlines,
        show_trajectories=args.trajectories,
        show_breakup_locus=not args.no_locus,
        coalesce_droplets=not args.no_coalesce,
        coalesce_primary_satellite=args.merge_primary_satellite,
        dpi=args.dpi,
        style=style_cfg,
    )

    ingest_and_plot(
        csv_file=args.csv_file,
        output_dir=args.out_dir,
        time_step_us=args.dt,
        t_max_us=args.t_max,
        steady_span_us=args.steady_span,
        config=cfg,
        ext=args.format,
    )