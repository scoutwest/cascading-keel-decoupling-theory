"""Section 25: the craton selection problem for a serpentinite layer at the crust base.
A craton can host the layer only if its Moho is cold enough to keep antigorite stable (T <= L, the antigorite limit)
and warm enough that the heat to be added to trigger breakdown (L - T) is small (<= dT).
Moho depth per craton cell: CRUST1.0 over the 114 Hasterok craton polygons (as in craton_moho_stats.py).
Moho temperature: the Section 24 conductive geotherm with the measured heat production (0.3-0.6 uW/m3).
Surface heat flow is the free variable: no global heat-flow grid is on disk. Archean heat flow is taken as
41 +/- 11 mW/m2 (garbled search snippet: unverified). Antigorite limit L at ~1.1 GPa: 500-550 C (long runs,
Gemba and Nakamura 2006 abstract) or 600-650 C (one-week runs, 2005 abstract)."""
import os, itertools
import numpy as np
import shapefile
from shapely.geometry import shape, Point
from shapely.prepared import prep
from shapely.strtree import STRtree
from math import erf, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(HERE, "..", "..")))
bnds = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.bnds")
lat = 89.5 - np.arange(180); lon = -179.5 + np.arange(360)
LAT, LON = np.meshgrid(lat, lon, indexing="ij"); LAT, LON = LAT.ravel(), LON.ravel()
sf = shapefile.Reader(ROOT + "/Noah's Pangea v2/data/hasterok_2022/cratons.shp")
geoms = [shape(s.__geo_interface__) for s in sf.shapes()]
prepped = [prep(g) for g in geoms]; tree = STRtree(geoms)
inside = np.zeros(len(LAT), bool)
for k in range(len(LAT)):
    pt = Point(LON[k], LAT[k])
    for i in tree.query(pt):
        if prepped[int(i)].contains(pt):
            inside[k] = True; break
w = np.cos(np.radians(LAT))[inside]
z = (bnds[:, 0] - bnds[:, 8])[inside]            # Moho depth below surface, km
print("craton cells:", inside.sum(), " Moho depth median %.1f km" % np.median(z))

T0 = 15.0
def moho_T(q0, a0, hr, k, zk):
    q0 = q0 * 1e-3; a0 = a0 * 1e-6; hr = hr * 1e3; zz = zk * 1e3
    return T0 + (q0 - a0 * hr) * zz / k + a0 * hr**2 / k * (1 - np.exp(-zz / hr))

PARAMS = list(itertools.product((0.3, 0.45, 0.6), (8.0, 10.0, 12.0), (2.5, 3.0)))
def cell_T(q0):                                   # temperature per cell, median over the parameter sets
    return np.median([moho_T(q0, a, h, k, z) for a, h, k in PARAMS], axis=0)

def wshare(mask): return float(w[mask].sum() / w.sum())

qs = np.arange(20, 62, 2.0)
Tq = {q: cell_T(q) for q in qs}
print("\nCraton-area-weighted median Moho temperature vs surface heat flow (median over heat production, decay depth, conductivity):")
for q in qs[::2]:
    t = Tq[q]; o = np.argsort(t); c = np.cumsum(w[o]) / w.sum()
    print(f"  q0 {q:4.0f} mW/m2: median {t[o][np.searchsorted(c,0.5)]:4.0f} C  (10-90%: {t[o][np.searchsorted(c,0.1)]:.0f}-{t[o][np.searchsorted(c,0.9)]:.0f})")

def cdf(x, mu=41.0, sd=11.0): return 0.5 * (1 + erf((x - mu) / (sd * sqrt(2))))
print("\nWINDOW: surface heat flow band that puts the median craton cell at L-dT <= T <= L")
print("  (and the share of Archean cratons, q0 ~ N(41, 11), that fall in that band)")
for L in (525, 625):
    for dT in (50, 100, 200):
        ok = [q for q in qs if np.median(Tq[q]) <= L and np.median(Tq[q]) >= L - dT]
        if ok:
            lo, hi = min(ok) - 1, max(ok) + 1
            print(f"  L {L} C, heat addable <= {dT:3d} C: q0 {lo:.0f}-{hi:.0f} mW/m2  -> {100*(cdf(hi)-cdf(lo)):4.1f}% of Archean heat-flow range")
        else:
            print(f"  L {L} C, heat addable <= {dT:3d} C: no band on the grid")
print("\nShare of craton AREA inside the window at fixed heat flow (L = 525 C, dT = 100 C):")
for q in (26, 30, 34, 38, 41, 44, 47):
    q = float(q); t = cell_T(q)
    print(f"  q0 {q:3.0f}: hold {100*wshare(t<=525):5.1f}%   window {100*wshare((t<=525)&(t>=425)):5.1f}%   above limit {100*wshare(t>525):5.1f}%")
