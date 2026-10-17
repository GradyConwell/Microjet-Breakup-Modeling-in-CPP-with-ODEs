# Authored by Opus 5.5, Gemini Flash 3.8, Gemini 4 Argon, with direction and material provided by Grady Conwell

# Figures for a microjet run, laid out like Van Hoeve et al. (2010) Figs. 3-6, in the Grady Conwell style.
# Usage: .venv/bin/python plotFrames.py [results/<run>]   (defaults to the latest run)
import glob
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Ellipse

H0, TC = 18.5, 10.119  # jet radius (µm), capillary time (µs)
PAPER = {"primary": (4.68, 4.9, 0.2), "satellite": (1.6, 1.9, 0.6)}  # model, experiment ± (m/s)

# Style tokens (Grady Conwell design system)
GROUND = "#f7f6f3"     # stalactite
INK = "#38373d"        # charcoal-art: text, axes, rules
MUTED = "#554340"      # french-roast: secondary text, reference lines
PRIMARY = "#28587a"    # blue-sapphire: jet and primary drops
SATELLITE = "#902c32"  # scarlet-sage: satellite drops (the one highlighted detail)
for f in glob.glob("fonts/*.[ot]tf"):  # put Fraunces here if it is not installed system-wide
    font_manager.fontManager.addfont(f)
installed = {f.name for f in font_manager.fontManager.ttflist}
DISPLAY = [f for f in ["Fraunces", "Georgia", "DejaVu Serif"] if f in installed]  # titles; details in Light
BODY = [f for f in ["Artifakt Element", "Founders Grotesk", "Inter", "DejaVu Sans"] if f in installed]
if "Fraunces" not in installed:
    print(f"Fraunces not found: using {DISPLAY[0]} for titles and details")
plt.rcParams.update({
    "figure.facecolor": GROUND, "axes.facecolor": GROUND, "savefig.facecolor": GROUND,
    "text.color": INK, "axes.labelcolor": MUTED, "axes.edgecolor": INK, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "font.family": BODY, "axes.titlesize": 15, "axes.titlelocation": "left", "axes.titlepad": 14,
    "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.frameon": False,
})
TITLE = {"family": DISPLAY, "weight": "normal", "color": INK}
DETAIL = {"family": DISPLAY, "weight": "light" if "Fraunces" in installed else "normal"}  # Fraunces Light
MONOGRAM = plt.imread("assets/monogram-black.png")[::40, ::40]

run = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob("results/*/"))[-1]
jet = np.loadtxt(f"{run}/jet.csv", delimiter=",", skiprows=1, ndmin=2)
t, L, H = jet[:, 0], jet[:, 1], jet[:, 2:]
drops = np.loadtxt(f"{run}/drops.csv", delimiter=",", skiprows=1, ndmin=2).reshape(-1, 4)
ev = np.genfromtxt(f"{run}/events.csv", delimiter=",", names=True, dtype=None, encoding=None)
neck = np.loadtxt(f"{run}/neck.csv", delimiter=",", skiprows=1, ndmin=2)
pinch = ev[ev["type"] == "pinch"]
steady = t[-1] - 100.0  # analyse the last 100 µs, after the start-up transient


def frame(i):
    """Jet profile (z, h) of frame i, closed at the tip, and the drops alive in that frame."""
    n = H.shape[1]
    z = np.append((np.arange(n) + 0.5) * L[i] / n, L[i])
    return z, np.append(H[i], 0.0), drops[np.isclose(drops[:, 0], t[i])]


def dropColor(R):
    return PRIMARY if R > H0 else SATELLITE


def style(ax, title="", xlabel="", ylabel=""):
    ax.set_title(title, fontdict=TITLE)
    ax.set_xlabel(xlabel, fontdict=DETAIL), ax.set_ylabel(ylabel, fontdict=DETAIL)
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontfamily(DISPLAY), lab.set_fontweight(DETAIL["weight"])


def save(fig, name):
    w, h = fig.get_size_inches()  # monogram watermark: bottom-right, 0.3 in tall, 10% opacity
    mark = fig.add_axes([1 - 0.45 / w, 0.1 / h, 0.3 / w, 0.3 / h], zorder=10)
    mark.imshow(MONOGRAM, alpha=0.1), mark.axis("off")
    fig.savefig(f"{run}/{name}.png", dpi=200, bbox_inches="tight", pad_inches=0.4)
    print(f"saved {run}/{name}.png")


# Fig. 3: time series of the jet, one vertical silhouette every 2 µs
fig, ax = plt.subplots(figsize=(12, 4))
s = 2 / (2 * 2 * H0)  # µs of x-axis per µm of radius
for i in np.flatnonzero((t >= steady) & (np.round(t) % 2 == 0)):
    z, h, d = frame(i)
    ax.fill_betweenx(z, t[i] - s * h, t[i] + s * h, color=PRIMARY, lw=0)
    for _, zc, R, _ in d:
        ax.add_patch(Ellipse((t[i], zc), 2 * R * s, 2 * R, color=dropColor(R), lw=0))
ax.set_ylim(0, min(drops[:, 1].max(), 1.2 * L.max()))  # nozzle at the bottom, flow upward
style(ax, "Jet time series (cf. paper Fig. 3b)", "time (µs)", "z (µm)")
save(fig, "fig3_timeseries")

# Fig. 4: drop velocities
fig, ax = plt.subplots(figsize=(7, 4))
late = drops[(drops[:, 0] >= steady) & (drops[:, 2] < 2 * H0)]  # skip the oversized start-up drops
print(f"\n{'':10s} {'R (µm)':>14s} {'v (m/s)':>14s}   paper v: model / experiment")
for kind, sel, c in [("primary", late[:, 2] > H0, PRIMARY), ("satellite", late[:, 2] <= H0, SATELLITE)]:
    model, exp, err = PAPER[kind]
    ax.axhspan(exp - err, exp + err, color=MUTED, alpha=0.08, lw=0)
    ax.axhline(model, color=MUTED, ls="--", lw=0.8)
    ax.plot(late[sel, 0], late[sel, 3], ".", color=c, ms=3, label=f"{kind} drops")
    born = pinch[(pinch["t_us"] >= steady) & ((pinch["R_um"] > H0) == (kind == "primary"))]
    if len(born):  # size and velocity of each drop at its pinch-off
        R, v = born["R_um"], born["v_ms"]
        print(f"{kind:10s} {R.mean():7.2f} ± {R.std():4.2f} {v.mean():7.2f} ± {v.std():4.2f}   {model} / {exp} ± {err}")
ax.set_ylim(0, 6)
ax.legend(loc="center right", labelcolor=INK)
style(ax, "Drop velocity (cf. paper Fig. 4)", "time (µs)", "velocity (m/s)")
ax.text(0, -0.22, "Dashed: paper model. Band: paper experiment.", transform=ax.transAxes, color=MUTED,
        fontsize=9, **DETAIL)
save(fig, "fig4_velocity")

# Fig. 5: pinch-off sequence after a primary pinch-off, aligned on the pinch location
prim = pinch[(pinch["R_um"] > H0) & (pinch["t_us"] >= steady)]
if len(prim):
    t0, zp = prim["t_us"][0], prim["z_um"][0] - prim["R_um"][0]
    print(f"\nprimary pinch-offs every {np.diff(prim['t_us']).mean():.1f} µs (forcing period 30.4 µs), "
          f"breakup length {(prim['z_um'] - prim['R_um']).mean():.0f} µm")
    print("pinch-off times after primary at", t0, "µs:", np.round(pinch["t_us"][(pinch["t_us"] >= t0)][:3] - t0, 2))
    fig, axs = plt.subplots(4, 3, figsize=(12, 7), sharex=True, sharey=True)
    for k, ax in enumerate(axs.T.flat):
        i = np.searchsorted(t, t0 + 2 * k - 1e-6)
        if i >= len(t):
            break
        z, h, d = frame(i)
        ax.fill_between(z - zp, -h, h, color=PRIMARY, lw=0)
        for _, zc, R, _ in d:
            ax.add_patch(plt.Circle((zc - zp, 0), R, color=dropColor(R), lw=0))
        ax.set_aspect("equal")
        ax.text(0.02, 0.82, f"{t[i] - t0:.0f} µs", transform=ax.transAxes, color=MUTED, fontsize=9, **DETAIL)
        style(ax)
    axs[0, 0].set(xlim=(-200, 60), ylim=(-45, 45))
    fig.supxlabel("distance from pinch-off (µm)", color=MUTED, **DETAIL)
    fig.supylabel("radius (µm)", color=MUTED, **DETAIL)
    fig.suptitle("Primary and satellite pinch-off (cf. paper Fig. 5)", x=0.06, ha="left", fontsize=15, **TITLE)
    save(fig, "fig5_pinchoff")

# Fig. 6: minimum neck radius vs time to pinch-off
fig, ax = plt.subplots(figsize=(5, 4))
for k, tp in enumerate(pinch["t_us"]):
    rows = neck[(neck[:, 0] == k) & (neck[:, 1] < tp)]
    if tp >= steady and len(rows):
        ax.plot(np.log10((tp - rows[:, 1]) / TC), np.log10(rows[:, 2] / H0), lw=1.2,
                color=dropColor(pinch["R_um"][k]))
x = np.array([-2.5, 0])
ax.plot(x, 2 / 3 * x - 0.2, color=MUTED, ls="--", lw=0.8, label="slope 2/3")
ax.set(xlim=(-2.5, 0.5), ylim=(-2.5, 0.2))
ax.legend(labelcolor=INK)
style(ax, "Neck collapse (cf. paper Fig. 6)", r"log($\tau/\tau_c$)", r"log($h_{min}/h_0$)")
save(fig, "fig6_neck")
