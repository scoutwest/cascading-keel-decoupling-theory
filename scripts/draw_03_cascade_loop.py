import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(REPO, "images")
NL = "\n"
GREEN, YELLOW, RED = "#2e8b57", "#d4a017", "#c0392b"

fig, ax = plt.subplots(figsize=(18, 10.5), dpi=150)
ax.set_xlim(0, 18)
ax.set_ylim(0, 10.5)
ax.axis("off")
fig.suptitle("The cascade: how one collapse could trigger the next", fontsize=16, fontweight="bold", y=0.975)

W, H = 5.0, 2.1
nodes = [
    (9.0, 9.3, "1. LOWER-CRUSTAL EARTHQUAKE", "Cold, dry granulite fails (reported to start" + NL + "the process in the Bergen Arcs, Norway)", GREEN),
    (14.6, 6.9, "2. FLUID PATHWAYS OPEN", "Faulting and fine fragments let fluid in;" + NL + "the reaction is limited by fluid supply", GREEN),
    (12.8, 3.9, "3. ECLOGITIZATION BEGINS", "Denser eclogite, or garnet granulite if pressure is low;" + NL + "days to months of fluid flow in one study (820 C)", YELLOW),
    (5.2, 3.9, "4. LAYER SHRINKS 9-17%, PORES OPEN", "Crust drops (my arithmetic); volume loss keeps" + NL + "the fluid pathways open (reported feedback)", YELLOW),
    (3.4, 6.9, "5. STRESS CHANGE LOADS NEIGHBOURS", "The next earthquake follows nearby. Transfer in" + NL + "dry lower crust is my synthesis, unverified", RED),
]
boxes = []
for (x, y, title, sub, col) in nodes:
    ax.add_patch(FancyBboxPatch((x - W / 2, y - H / 2), W, H, boxstyle="round,pad=0.05,rounding_size=0.2", fc="white", ec=col, lw=3.2, zorder=3))
    ax.add_patch(Rectangle((x - W / 2 + 0.02, y - H / 2 + 0.05), 0.22, H - 0.1, fc=col, ec="none", zorder=4))
    ax.text(x + 0.12, y + 0.45, title, ha="center", va="center", fontsize=11, fontweight="bold", zorder=5)
    ax.text(x + 0.12, y - 0.35, sub, ha="center", va="center", fontsize=9.3, zorder=5)
    boxes.append((x, y))


def arrow(p, q, rad=0.25):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=28, lw=3, color="#444", connectionstyle="arc3,rad=%s" % rad, zorder=2))


# clockwise arrows between nodes
arrow((11.6, 9.0), (13.8, 8.0), -0.25)
arrow((14.0, 5.8), (13.4, 5.05), -0.2)
arrow((10.2, 3.9), (7.8, 3.9), 0.0)
arrow((4.8, 4.9), (3.9, 5.8), -0.2)
arrow((3.8, 7.95), (6.4, 9.0), -0.25)
ax.text(2.6, 8.9, "repeats" + NL + "across the craton", ha="center", va="center", fontsize=9.5, style="italic", color="#444")

# centre
ax.add_patch(FancyBboxPatch((6.3, 5.2), 5.4, 2.1, boxstyle="round,pad=0.05,rounding_size=0.25", fc="#fff3e0", ec="#6b3a00", lw=2.5, zorder=3))
ax.text(9.0, 6.95, "RESULT", ha="center", va="center", fontsize=11, fontweight="bold", color="#6b3a00", zorder=5)
ax.text(9.0, 6.3, "Crust drops 0.48-0.91 km (0.30-0.57 mi);" + NL + "the ocean covers the land", ha="center", va="center", fontsize=10.2, zorder=5)
ax.text(9.0, 5.5, "keel untouched (P9)", ha="center", va="center", fontsize=9.3, style="italic", zorder=5)

# needs box
ax.add_patch(FancyBboxPatch((0.4, 0.35), 8.2, 1.15, boxstyle="round,pad=0.05,rounding_size=0.15", fc="#f1f1f1", ec="#999", lw=1.5, zorder=3))
ax.text(4.5, 0.92, "Needs: cold, dry, brittle lower crust (cratons only, margins uncovered: O16);" + NL +
        "fluid at the crust base, 35-42 km deep (source unresolved: O17);" + NL + "and enough eclogite driving force there (O18)", ha="center", va="center", fontsize=9.6, fontweight="bold", zorder=5)

# legend
lx, ly = 9.2, 0.35
ax.add_patch(FancyBboxPatch((lx, ly), 8.4, 1.15, boxstyle="round,pad=0.05,rounding_size=0.15", fc="white", ec="#999", lw=1.5, zorder=3))
for k, (col, txt) in enumerate([(GREEN, "reported in the literature (search summaries, not full papers)"),
                                (YELLOW, "single source, or partly my arithmetic"),
                                (RED, "my synthesis, unverified")]):
    yy = ly + 0.92 - k * 0.33
    ax.add_patch(Rectangle((lx + 0.2, yy - 0.1), 0.35, 0.2, fc=col, ec="none", zorder=5))
    ax.text(lx + 0.7, yy, txt, ha="left", va="center", fontsize=9, zorder=5)

fig.text(0.5, 0.012, "Schematic: timing is not shown. Evidence from Sections 15-17 of the model; the loop as a whole is my synthesis from the cited findings.",
         ha="center", va="center", fontsize=9.3)
fig.savefig(os.path.join(OUT, "03_cascade_loop.png"))
fig.savefig(os.path.join(OUT, "03_cascade_loop.svg"))
print("saved")
