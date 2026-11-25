import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(REPO, "images")
os.makedirs(OUT, exist_ok=True)
NL = "\n"

C_CRUST, C_KEEL, C_ASTH, C_WATER, C_SKY, C_LAYER = "#b08a5a", "#5b7aa6", "#e8a35a", "#7fb8e6", "#f4f8fc", "#c0392b"

# schematic units: 1.2 units = 35 km crust. Layer and subsidence drawn 7x too large.
CR = 1.2
L_PRE = 1.44      # 7.0 km drawn 6x
L_POST = 1.25     # 6.1 km drawn 6x
DROP = L_PRE - L_POST   # 0.19 -> 0.91 km drawn 6x
KT = CR + L_PRE   # top of keel (does not move)
KB = 4.9          # keel base (drawn compressed)
BOT = 8.2

fig, axes = plt.subplots(1, 3, figsize=(18, 10.5), dpi=150)
fig.suptitle("A thick Moho layer collapses to a thinner, denser layer, the crust drops, and the ocean covers the land",
             fontsize=15.5, fontweight="bold", y=0.985)
titles = ["1. BEFORE" + NL + "Thick Moho layer, 7.0 km (4.3 mi)",
          "2. DURING" + NL + "Collapse spreads as a cascade",
          "3. AFTER" + NL + "Layer 6.1 km (3.8 mi), crust down 0.91 km"]

x = np.linspace(0.5, 9.5, 300)


def smooth(t):
    t = np.clip(t, 0, 1)
    return t * t * (3 - 2 * t)


for i, ax in enumerate(axes):
    ax.set_xlim(0, 10)
    ax.set_ylim(BOT + 0.6, -2.2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(titles[i], fontsize=12, fontweight="bold", pad=6)
    ax.add_patch(Rectangle((0, -2.2), 10, 2.2, color=C_SKY, zorder=0))
    ax.add_patch(Rectangle((0.5, CR), 9, BOT + 0.6 - CR, fc=C_ASTH, ec="none", zorder=1))

    if i == 0:
        surf = np.zeros_like(x)
    elif i == 1:
        surf = DROP * (1 - smooth((x - 3.5) / 3.0))       # collapsed on the left, still thick on the right
    else:
        surf = np.full_like(x, DROP)

    crust_bot = CR + surf
    ax.add_patch(Polygon(np.vstack([np.column_stack([x, surf]), np.column_stack([x[::-1], crust_bot[::-1]])]),
                         fc=C_CRUST, ec="k", lw=0.8, zorder=3))
    ax.text(5, 0.55 + surf.mean(), "CRUST (35 km, 21.7 mi)", ha="center", va="center", fontsize=10.5,
            fontweight="bold", color="white", zorder=8)

    # weak layer between crust and keel
    ax.add_patch(Polygon(np.vstack([np.column_stack([x, crust_bot]), np.column_stack([x[::-1], np.full_like(x, KT)[::-1]])]),
                         fc=C_LAYER, ec="k", lw=0.8, zorder=3, hatch="////" if i >= 1 else None))
    # keel (never changes)
    ax.add_patch(Polygon([[1.2, KT], [8.8, KT], [8.2, KB], [1.8, KB]], fc=C_KEEL, ec="k", lw=0.8, zorder=2))
    ax.text(5, (KT + KB) / 2, "KEEL: UNCHANGED" + NL + "150-250 km (93-155 mi)", ha="center", va="center",
            fontsize=10, fontweight="bold", color="white", zorder=8)
    ax.text(0.9, 7.8, "ASTHENOSPHERE" + NL + "(hot, weak)", ha="left", va="center", fontsize=10, fontweight="bold",
            color="#6b3a00", zorder=8)

    # in-layer labels
    if i == 0:
        ax.text(5, CR + L_PRE / 2, "THICK MOHO LAYER 7.0 km (4.3 mi)", ha="center", va="center", fontsize=9.5,
                fontweight="bold", color="white", zorder=9)
    elif i == 1:
        ax.text(2.4, CR + (L_POST) / 2 + DROP, "collapsed", ha="center", va="center", fontsize=9, fontweight="bold", color="white", zorder=9)
        ax.text(7.9, CR + L_PRE / 2, "still thick", ha="center", va="center", fontsize=9, fontweight="bold", color="white", zorder=9)
        ax.annotate("", xy=(2.9, 5.9), xytext=(7.4, 5.9), arrowprops=dict(arrowstyle="-|>", color="k", lw=2.5), zorder=9)
        ax.text(5.1, 6.4, "collapse front moves this way", ha="center", va="center", fontsize=9, zorder=9)
    else:
        ax.text(5, CR + DROP + L_POST / 2, "COLLAPSED LAYER 6.1 km (3.8 mi)", ha="center", va="center", fontsize=9.5,
                fontweight="bold", color="white", zorder=9)

    # water
    if i == 0:
        ax.plot([0, 10], [0, 0], color="#1f6fb2", lw=2, zorder=6)
        ax.text(5, -1.0, "Land: mean 0.2 km (0.12 mi) above sea level," + NL + "highest peaks under 0.61 km (0.38 mi)",
                ha="center", va="center", fontsize=10, zorder=7)
    elif i == 1:
        ax.plot([0, 10], [0, 0], color="#1f6fb2", lw=2, zorder=6)
        ax.add_patch(Polygon(np.vstack([np.column_stack([x, np.zeros_like(x)]), np.column_stack([x[::-1], surf[::-1]])]),
                             fc=C_WATER, alpha=0.9, ec="none", zorder=4))
        ax.text(5, -1.0, "Land sags where the layer has collapsed." + NL + "Water spreads over the low ground.",
                ha="center", va="center", fontsize=10, zorder=7)
    else:
        ax.add_patch(Rectangle((0.5, 0), 9, DROP, color=C_WATER, alpha=0.9, zorder=4))
        ax.plot([0, 10], [0, 0], color="#1f6fb2", lw=2, zorder=6)
        ax.annotate("OCEAN", xy=(5, 0.1), xytext=(5, -0.5), ha="center", va="center", fontsize=13, fontweight="bold",
                    color="#0b3d6e", arrowprops=dict(arrowstyle="-|>", color="#0b3d6e", lw=1.8), zorder=9)
        ax.text(5, -1.4, "Land 0.91 km (0.57 mi) lower; ocean covers a 0.61 km peak.", ha="center", va="center",
                fontsize=10, zorder=7)

# left depth scale on panel 1
ax = axes[0]
for y, t in [(CR, "Crust base" + NL + "depth 35 km (21.7 mi)"),
             (KT, "Keel top (does not move)" + NL + "depth 42 km (26.1 mi)"),
             (KB, "Keel base" + NL + "depth about 200 km (124 mi)")]:
    ax.plot([0.0, 1.1], [y, y], color="k", lw=0.8, ls="--", zorder=9)
    ax.text(-0.15, y, t, ha="right", va="center", fontsize=8.6, clip_on=False)

fig.text(0.5, 0.135,
         "Key numbers (15% density increase, eclogite, middle case):  subsidence = layer thickness lost = 7.0 km - 6.1 km = 0.91 km (0.57 mi), a 13% collapse   |   range: 5.5-10 km before, 9-17% collapse" + NL +
         "With 1.3x10^8 km3 of shallow water released only 0.48 km (0.30 mi) is needed: 3.7 km before, 3.2 km after." + NL +
         "Low end (garnet granulite, 10% denser, more likely at craton pressures): 10.0 km before for 0.91 km of drop, 5.3 km before for 0.48 km." + NL +
         "UNTESTED: needs fluid-driven gabbro-to-eclogite change; speed at craton temperature is unmeasured (scaled: days to years).",
         ha="center", va="center", fontsize=9.4, fontweight="bold", bbox=dict(fc="#f1f1f1", ec="#999", boxstyle="round,pad=0.5"))
fig.text(0.5, 0.032,
         "Schematic, NOT to scale: the layer and the drop are drawn about 6 times too thick compared with the crust, and the keel is compressed. The keel is untouched, so the idea satisfies premise P9." + NL +
         "The layer thicknesses are computed from the density change (Section 15.2), not measured. Reported high-velocity lower-crust bodies in the Fennoscandian Shield are 6-8 km thick." + NL +
         "Likely acts only where the lower crust is cold (cratons); the warm margins are uncovered (O16)." + NL +
         "Model crust is 35 km; the craton median Moho is 39 km at about 1.07 GPa, below the expected eclogite field (O18, Section 19).",
         ha="center", va="center", fontsize=9.5)

fig.subplots_adjust(left=0.15, right=0.97, top=0.89, bottom=0.20, wspace=0.35)
fig.savefig(os.path.join(OUT, "01_collapse_overview.png"))
fig.savefig(os.path.join(OUT, "01_collapse_overview.svg"))
print("saved")
