import os
import numpy as np
import shapefile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.geometry import shape, Point
from shapely.prepared import prep
from shapely.strtree import STRtree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, ".."))
# The CRUST1.0 and craton data sit in the sibling project folders (see scripts/README.md).
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(REPO, "..")))
OUT = os.path.join(REPO, "images")
NL = "\n"

bnds = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.bnds")
rho = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.rho")
lat = 89.5 - np.arange(180)
lon = -179.5 + np.arange(360)
LAT, LON = np.meshgrid(lat, lon, indexing="ij")
LAT, LON = LAT.ravel(), LON.ravel()

sf = shapefile.Reader(ROOT + "/Noah's Pangea v2/data/hasterok_2022/cratons.shp")
geoms = [shape(s.__geo_interface__) for s in sf.shapes()]
prepped = [prep(g) for g in geoms]
tree = STRtree(geoms)
inside = np.zeros(len(LAT), dtype=bool)
for k in range(len(LAT)):
    pt = Point(LON[k], LAT[k])
    for idx in tree.query(pt):
        if prepped[int(idx)].contains(pt):
            inside[k] = True
            break

w = np.cos(np.radians(LAT))
g = 9.81
P = np.zeros(len(LAT))
for i in range(8):
    th = np.clip(bnds[:, i] - bnds[:, i + 1], 0, None)
    P += rho[:, i] * 1000.0 * g * th * 1000.0
P /= 1e9
lower = bnds[:, 7] - bnds[:, 8]
m = inside
print("cells", m.sum())

fig, axes = plt.subplots(1, 2, figsize=(18, 9), dpi=150)
fig.suptitle("Craton Moho pressure and lower-crust thickness (CRUST1.0, Hasterok craton polygons)", fontsize=15.5, fontweight="bold", y=0.975)

# left: pressure histogram
ax = axes[0]
bins = np.arange(0.6, 1.7, 0.04)
h, e = np.histogram(P[m], bins=bins, weights=w[m])
h = 100 * h / h.sum()
ax.bar(e[:-1], h, width=0.04, align="edge", color="#5b7aa6", ec="white", zorder=3)
ax.axvspan(1.3, 1.6, color="#c0392b", alpha=0.18, zorder=1)
ax.set_ylim(0, h.max() * 1.3)
ax.text(1.45, h.max() * 1.22, "eclogite expected" + NL + "at 500-600 C" + NL + "(one source:" + NL + "1.3-1.6 GPa)", ha="center", va="top", fontsize=10, fontweight="bold", color="#8b2a1f")
ax.axvline(1.07, color="k", lw=2, ls="--", zorder=4)
ax.text(1.09, h.max() * 1.16, "median 1.07 GPa", ha="left", va="center", fontsize=10.5, fontweight="bold")
ax.axvline(1.0, color="#2e8b57", lw=2, ls=":", zorder=4)
ax.text(0.985, h.max() * 1.1, "about 1 GPa: a source puts eclogite" + NL + "stable here for one rock", ha="right", va="center", fontsize=9, color="#2e8b57", fontweight="bold")
ax.set_xlabel("Pressure at the Moho (GPa)", fontsize=11)
ax.set_ylabel("Share of craton area (%)", fontsize=11)
ax.set_xlim(0.6, 1.7)
ticks = [0.8, 1.0, 1.2, 1.4, 1.6]
ax.set_xticks(ticks)
ax.set_xticklabels(["%.1f GPa" % t + NL + "~%d km" % round(t * 1e9 / (3000 * g) / 1000) + NL + "(%d mi)" % round(t * 1e9 / (3000 * g) / 1609.34) for t in ticks], fontsize=9.5)
ax.set_title("Only 3% of craton area reaches 1.3 GPa", fontsize=12.5, fontweight="bold")
ax.grid(axis="y", alpha=0.3, zorder=0)

# right: lower crust thickness
ax = axes[1]
bins2 = np.arange(4, 20, 0.5)
h2, e2 = np.histogram(lower[m], bins=bins2, weights=w[m])
h2 = 100 * h2 / h2.sum()
ax.bar(e2[:-1], h2, width=0.5, align="edge", color="#b08a5a", ec="white", zorder=3)
ax.axvspan(2.9, 10.0, color="#2e8b57", alpha=0.15, zorder=1)
ax.set_ylim(0, h2.max() * 1.3)
ax.text(4.3, h2.max() * 1.22, "layer needed before collapse:" + NL + "2.9-10 km (1.8-6.2 mi)", ha="left", va="top", fontsize=10.5, fontweight="bold", color="#1f6b3f")
ax.axvline(12.5, color="k", lw=2, ls="--", zorder=4)
ax.text(12.8, h2.max() * 1.16, "median 12.5 km (7.8 mi)", ha="left", va="center", fontsize=10.5, fontweight="bold")
ax.set_xlabel("CRUST1.0 lower-crust layer thickness (km)", fontsize=11)
ax.set_ylabel("Share of craton area (%)", fontsize=11)
ax.set_xlim(4, 20)
t2 = [5, 8, 10, 12.5, 15, 18]
ax.set_xticks(t2)
ax.set_xticklabels(["%g km" % t + NL + "(%.1f mi)" % (t * 0.621371) for t in t2], fontsize=9.5)
ax.set_title("The lower crust is thick enough (97% of area at least 7 km)", fontsize=12.5, fontweight="bold")
ax.grid(axis="y", alpha=0.3, zorder=0)

fig.text(0.5, 0.075,
         "Reading: the layer is thick enough, but the pressure at the craton Moho is mostly below where one source expects eclogite at craton temperature (O18)." + NL +
         "If eclogite is not reached the product is garnet granulite, about 10% denser, which needs a thicker layer for the same drop (Section 19).",
         ha="center", va="center", fontsize=10, fontweight="bold", bbox=dict(fc="#f1f1f1", ec="#999", boxstyle="round,pad=0.5"))
fig.text(0.5, 0.022,
         "Area-weighted over 5,709 one-degree cells inside the craton polygons (8.0% of the sphere). Pressure from CRUST1.0 layer densities; the km scale on the left is approximate (uniform 3,000 kg/m3)." + NL +
         "CRUST1.0 is a seismic model with assumed velocity-density relations, and whether the lower crust is mafic is disputed. Eclogite-field ranges are from cited studies and are uncertain.",
         ha="center", va="center", fontsize=9.3)
fig.subplots_adjust(left=0.07, right=0.97, top=0.90, bottom=0.2, wspace=0.18)
fig.savefig(os.path.join(OUT, "04_craton_moho_pressure.png"))
fig.savefig(os.path.join(OUT, "04_craton_moho_pressure.svg"))
print("saved")
