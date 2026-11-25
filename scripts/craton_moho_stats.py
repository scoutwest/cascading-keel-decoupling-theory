import os
import numpy as np
import shapefile
from shapely.geometry import shape, Point
from shapely.prepared import prep
from shapely.strtree import STRtree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, ".."))
# The CRUST1.0 and craton data sit in the sibling project folders (see scripts/README.md).
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(REPO, "..")))
bnds = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.bnds")
rho = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.rho")
assert bnds.shape == (64800, 9) and rho.shape == (64800, 9)

lat = 89.5 - np.arange(180)
lon = -179.5 + np.arange(360)
LAT, LON = np.meshgrid(lat, lon, indexing="ij")
LAT, LON = LAT.ravel(), LON.ravel()

sf = shapefile.Reader(ROOT + "/Noah's Pangea v2/data/hasterok_2022/cratons.shp")
print("craton fields:", [f[0] for f in sf.fields[1:]])
geoms = [shape(s.__geo_interface__) for s in sf.shapes()]
print("craton polygons:", len(geoms))
prepped = [prep(g) for g in geoms]
tree = STRtree(geoms)

inside = np.zeros(len(LAT), dtype=bool)
for k in range(len(LAT)):
    pt = Point(LON[k], LAT[k])
    for idx in tree.query(pt):
        if prepped[int(idx)].contains(pt):
            inside[k] = True
            break
print("craton cells:", inside.sum())

w = np.cos(np.radians(LAT))
# Moho depth below the surface (km): surface = col 0, Moho = top of mantle = col 8
surf = bnds[:, 0]
moho_depth = surf - bnds[:, 8]
lower_thick = bnds[:, 7] - bnds[:, 8]
lower_rho = rho[:, 7]

# pressure at the Moho (GPa) from the layer densities (water, ice, sediments, crust layers)
g = 9.81
P = np.zeros(len(LAT))
for i in range(8):
    th = bnds[:, i] - bnds[:, i + 1]
    th = np.clip(th, 0, None)
    P += rho[:, i] * 1000.0 * g * th * 1000.0
P /= 1e9

m = inside
def wq(x, q):
    xs = x[m]; ws = w[m]
    o = np.argsort(xs); xs = xs[o]; ws = ws[o]
    c = np.cumsum(ws) / ws.sum()
    return [float(xs[np.searchsorted(c, qq)]) for qq in q]

def wfrac(cond):
    return float((w[m] * cond[m]).sum() / w[m].sum())

print("area-weighted Moho depth percentiles (5,25,50,75,95) km:", [round(v, 1) for v in wq(moho_depth, [.05, .25, .5, .75, .95])])
print("area-weighted Moho pressure percentiles GPa:", [round(v, 2) for v in wq(P, [.05, .25, .5, .75, .95])])
for z in (35, 38, 40, 42, 45, 50):
    print(f"fraction of craton area with Moho deeper than {z} km: {wfrac(moho_depth > z):.2f}")
for pp in (1.0, 1.2, 1.3, 1.4, 1.5, 1.6):
    print(f"fraction with Moho pressure >= {pp} GPa: {wfrac(P >= pp):.2f}")
print("lower-crust layer thickness percentiles km:", [round(v, 1) for v in wq(lower_thick, [.05, .25, .5, .75, .95])])
for t in (3, 5.5, 7, 10):
    print(f"fraction with CRUST1.0 lower-crust layer >= {t} km: {wfrac(lower_thick >= t):.2f}")
print("lower-crust density percentiles g/cm3:", [round(v, 2) for v in wq(lower_rho, [.05, .5, .95])])
print("craton area share of sphere (cos-weighted cells):", float(w[m].sum() / w.sum()))
# crude cross-check: all land
land = surf > 0
print("all-land median Moho depth (area-weighted) km:", end=" ")
xs = moho_depth[land]; ws = w[land]; o = np.argsort(xs); c = np.cumsum(ws[o]) / ws.sum()
print(round(float(xs[o][np.searchsorted(c, .5)]), 1))
