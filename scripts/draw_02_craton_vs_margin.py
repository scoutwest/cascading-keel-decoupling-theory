import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(REPO, "images")
NL = "\n"
C_CRUST, C_KEEL, C_ASTH, C_WATER, C_SKY, C_LAYER = "#b08a5a", "#5b7aa6", "#e8a35a", "#7fb8e6", "#f4f8fc", "#c0392b"
C_MARGIN = "#c9a27a"

# schematic units (not to scale). 1.2 units = 35 km craton crust. Layer 7.0 km drawn 6x thick.
CR = 1.2
L_PRE, L_POST = 1.44, 1.25
DROP = L_PRE - L_POST
KB = 4.9
BOT = 8.0
XC0, XC1 = 0.5, 6.2     # craton
XM0, XM1 = 6.2, 9.0     # margin
XO1 = 10.0              # ocean to the right
MC = 0.85               # margin crust thickness (thinned), schematic

fig, axes = plt.subplots(1, 2, figsize=(18, 10.5), dpi=150)
fig.suptitle("Where the collapse can act: cold craton interior and warm margin",
             fontsize=16, fontweight="bold", y=0.985)
titles = ["BEFORE" + NL + "Thick layer under the craton, unreacted", "AFTER" + NL + "Craton layer collapsed and flooded; margin still dry"]

for i, ax in enumerate(axes):
    ax.set_xlim(0, XO1)
    ax.set_ylim(BOT + 0.6, -1.9)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(titles[i], fontsize=12.5, fontweight="bold", pad=6)
    ax.add_patch(Rectangle((0, -1.9), XO1, 1.9, color=C_SKY, zorder=0))
    ax.add_patch(Rectangle((XC0, 0), XO1 - XC0, BOT + 0.6, fc=C_ASTH, ec="none", zorder=1))

    top_c = 0.0 if i == 0 else DROP
    sea = 0.15 if i == 0 else 0.10
    lay = L_PRE if i == 0 else L_POST

    # craton crust and layer and keel
    ax.add_patch(Rectangle((XC0, top_c), XC1 - XC0, CR, fc=C_CRUST, ec="k", lw=0.8, zorder=3))
    ax.add_patch(Rectangle((XC0, CR + top_c), XC1 - XC0, lay, fc=C_LAYER, ec="k", lw=0.8, zorder=3, hatch="////" if i else None))
    kt = CR + L_PRE
    ax.add_patch(Polygon([[1.0, kt], [XC1 - 0.5, kt], [XC1 - 1.0, KB], [1.5, KB]], fc=C_KEEL, ec="k", lw=0.8, zorder=2))
    ax.text((XC0 + XC1) / 2, 0.55 + top_c, "CRATON CRUST", ha="center", va="center", fontsize=10.5, fontweight="bold", color="white", zorder=8)
    ax.text((XC0 + XC1) / 2, CR + top_c + lay / 2, ("THICK LAYER 7.0 km (4.3 mi)" if i == 0 else "COLLAPSED LAYER 6.1 km (3.8 mi)"),
            ha="center", va="center", fontsize=9.5, fontweight="bold", color="white", zorder=9)
    ax.text((1.5 + XC1 - 1.0) / 2 + 0.15, (kt + KB) / 2, "KEEL" + NL + "UNCHANGED", ha="center", va="center", fontsize=10, fontweight="bold", color="white", zorder=8)
    # margin: thinner crust, no layer, warm
    ax.add_patch(Rectangle((XM0, 0.0), XM1 - XM0, MC, fc=C_MARGIN, ec="k", lw=0.8, zorder=3))
    ax.text((XM0 + XM1) / 2, MC / 2, "MARGIN CRUST", ha="center", va="center", fontsize=9, fontweight="bold", color="white", zorder=8)
    # ocean on the right
    ax.add_patch(Rectangle((XM1, sea), XO1 - XM1, 1.6 - sea, fc=C_WATER, alpha=0.9, zorder=2))
    ax.add_patch(Rectangle((XM1, 1.6), XO1 - XM1, 0.5, fc="#8a7a68", ec="k", lw=0.6, zorder=3))
    ax.text((XM1 + XO1) / 2, 0.9, "OCEAN", ha="center", va="center", fontsize=11, fontweight="bold", color="#0b3d6e", zorder=9)
    ax.plot([0, XO1], [sea, sea], color="#1f6fb2", lw=2, zorder=6)
    ax.text(3.3, -1.25, "COLD, DRY lower crust: eclogite stability" + NL + "uncertain at this depth (O18)", ha="center", va="center", fontsize=9.5, fontweight="bold", color="#6b3a00", zorder=9)
    ax.text(7.6, -1.25, "WARM lower crust: granulite stays stable," + NL + "no collapse (O16)", ha="center", va="center", fontsize=9.5, fontweight="bold", color="#6b3a00", zorder=9)

    if i == 1:
        ax.add_patch(Rectangle((XC0, sea), XC1 - XC0, top_c - sea, color=C_WATER, alpha=0.9, zorder=4))
        ax.text(3.3, -0.35, "Craton land is under water", ha="center", va="center", fontsize=10, fontweight="bold", color="#0b3d6e", zorder=9)
        ax.text(7.6, -0.35, "Margin land is above water", ha="center", va="center", fontsize=10, fontweight="bold", color="#6b3a00", zorder=9)
    else:
        ax.text(3.3, -0.35, "Craton: land above sea level", ha="center", va="center", fontsize=10, zorder=9)
        ax.text(7.6, -0.35, "Margin: land above sea level", ha="center", va="center", fontsize=10, zorder=9)

    # fluid path (dashed) from below the margin up to the craton layer
    ax.annotate("", xy=(XC1 - 0.05, CR + top_c + lay * 0.75), xytext=(6.6, 7.4),
                arrowprops=dict(arrowstyle="-|>", color="#0b3d6e", lw=2.2, ls="--", connectionstyle="arc3,rad=0.08"), zorder=10)
    ax.text(7.9, 5.4, "fluid needed at the" + NL + "crust base: source" + NL + "unresolved (O17)", ha="center", va="center", fontsize=9.2,
            fontweight="bold", color="#0b3d6e", zorder=10)
    ax.text(0.9, 7.6, "ASTHENOSPHERE (hot, weak)", ha="left", va="center", fontsize=9.5, fontweight="bold", color="#6b3a00", zorder=8)

fig.text(0.5, 0.127,
         "Craton: crust 35 km (21.7 mi), layer 7.0 km (4.3 mi) before and 6.1 km (3.8 mi) after, crust drops 0.91 km (0.57 mi), keel 150-250 km (93-155 mi) unchanged." + NL +
         "Margin: not covered by the collapse because the lower crust is warm; with stretching removed as a flood source, the margins (18% of Pangea) need another source or stay dry (O16)." + NL +
         "Craton Moho (CRUST1.0): median 39 km (24 mi), 1.07 GPa; only 3% of craton area reaches 1.3 GPa, where eclogite is expected at 500-600 C (Section 19).",
         ha="center", va="center", fontsize=9.8, fontweight="bold", bbox=dict(fc="#f1f1f1", ec="#999", boxstyle="round,pad=0.5"))
fig.text(0.5, 0.04,
         "Schematic, NOT to scale: the layer and the drop are drawn about 6 times too thick compared with the crust, and the keel is compressed. The margin crust is drawn thinner for illustration." + NL +
         "Cold-crust coverage and the fluid path are inferences from the findings in Sections 15-17, not measurements.",
         ha="center", va="center", fontsize=9.5)
fig.subplots_adjust(left=0.03, right=0.97, top=0.89, bottom=0.2, wspace=0.08)
fig.savefig(os.path.join(OUT, "02_craton_vs_margin.png"))
fig.savefig(os.path.join(OUT, "02_craton_vs_margin.svg"))
print("saved")
