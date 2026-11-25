"""First-pass geometry for Section 22: how far is any point on continental crust from a
province boundary (candidate suture / vent route)? Uses the Hasterok et al. (2022) provinces.
Boundaries counted as 'internal' are those shared between two continental/transitional provinces
(coastlines against oceanic provinces are excluded). Distances are great-circle, via 3-D chord lengths.
Deep-time ages in the data are not used (premise P1)."""
import os
import numpy as np
import shapefile
import shapely
from shapely.geometry import shape
from shapely import STRtree
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(HERE, "..", "..")))
R = 6371.0
import sys
STRICT = len(sys.argv) > 1 and sys.argv[1] == 'strict'  # only boundaries between provinces of different type
print('mode:', 'strict (different province types only)' if STRICT else 'all shared boundaries')
r = shapefile.Reader(ROOT + "/Noah's Pangea v2/data/hasterok_2022/global_gprv.shp")
geoms, recs = [], []
for sr in r.shapeRecords():
    g = shape(sr.shape.__geo_interface__)
    if not g.is_valid:
        g = g.buffer(0)
    geoms.append(g); recs.append(sr.record.as_dict())
land = [i for i, x in enumerate(recs) if x["crust_type"] in ("continental", "transitional")]
lg = [geoms[i] for i in land]
tree = STRtree(lg)

# Shared boundaries between neighbouring land provinces
lines = []
for a, g in enumerate(lg):
    for b in tree.query(g):
        if b <= a:
            continue
        sh = g.boundary.intersection(lg[b].boundary)
        if not sh.is_empty and sh.length > 0:
            same = recs[land[a]]["prov_type"] == recs[land[int(b)]]["prov_type"]
            if STRICT and same:
                continue
            lines.append((a, int(b), sh))

def pts_on(geom, step=0.05):
    out = []
    for part in getattr(geom, "geoms", [geom]):
        if part.geom_type not in ("LineString", "LinearRing"):
            continue
        n = max(2, int(part.length / step))
        out += [part.interpolate(t, normalized=True).coords[0] for t in np.linspace(0, 1, n)]
    return out

def xyz(lon, lat):
    lon, lat = np.radians(lon), np.radians(lat)
    return np.c_[R*np.cos(lat)*np.cos(lon), R*np.cos(lat)*np.sin(lon), R*np.sin(lat)]

bpts = np.array([p for _, _, sh in lines for p in pts_on(sh)])
kd = cKDTree(xyz(bpts[:, 0], bpts[:, 1]))
print("province pairs sharing a boundary:", len(lines), " boundary sample points:", len(bpts))

# Sample points: 1-degree grid cells inside continental provinces (cos-latitude weights)
lons, lats = np.meshgrid(np.arange(-179.5, 180, 1.0), np.arange(-89.5, 90, 1.0))
P = shapely.points(lons.ravel(), lats.ravel())
cont = [i for i in land if recs[i]["crust_type"] == "continental"]
ctree = STRtree([geoms[i] for i in cont])
pi, gi = ctree.query(P, predicate="within")
first = {}
for p, g in zip(pi, gi):
    first.setdefault(int(p), cont[int(g)])
idx = np.array(sorted(first))
ptype = np.array([recs[first[i]]["prov_type"] for i in idx])
w = np.cos(np.radians(lats.ravel()[idx]))
d, _ = kd.query(xyz(lons.ravel()[idx], lats.ravel()[idx]))

def wq(x, w, qs):
    o = np.argsort(x); c = np.cumsum(w[o]) / w.sum()
    return [float(x[o][np.searchsorted(c, q)]) for q in qs]

qs = (0.25, 0.5, 0.75, 0.95)
groups = {"all continental": np.ones(len(idx), bool),
          "craton or shield": np.isin(ptype, ("craton", "shield")),
          "orogenic/accretionary/arc": np.isin(ptype, ("orogenic belt", "accretionary complex", "volcanic arc")),
          "rift/passive margin/basin": np.isin(ptype, ("passive margin", "wide rift", "narrow rift", "basin", "foredeep basin"))}
print(f"{'group':30s} {'cells':>6s} {'area%':>6s}  distance to nearest shared boundary, km (25/50/75/95%)")
for k, m in groups.items():
    print(f"{k:30s} {m.sum():6d} {100*w[m].sum()/w.sum():6.1f}  " + " / ".join(f"{v:6.0f}" for v in wq(d[m], w[m], qs)))
for L in (100, 250, 500, 1000):
    print(f"share of continental area within {L:4d} km of a shared boundary: {100*w[d<=L].sum()/w.sum():5.1f}%")
