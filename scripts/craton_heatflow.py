"""Section 29 (uses Section 25 setup): the craton selection problem for a serpentinite layer at the crust base.
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


import csv
hf = os.environ.get("HEATFLOW_GRID", os.path.join(HERE, "..", "data", "heatflow", "2019GC008389-sup-0003-Data_Set_SI-S01", "HFgrid14.csv"))
G = np.full((360, 720), np.nan); O = np.full((360, 720), np.nan)
with open(hf) as f:
    r = csv.reader(f, delimiter=";"); next(r)
    for lo, la, p, s, o in r:
        i = int((89.75 - float(la)) / 0.5 + 0.5); j = int((float(lo) + 179.75) / 0.5 + 0.5)
        G[i, j] = float(p)
        if o.strip(): O[i, j] = float(o)
def cell_mean(A):
    B = A.reshape(180, 2, 360, 2)
    return np.nanmean(B, axis=(1, 3)).ravel()
q_pred = cell_mean(G)[inside]; q_obs = cell_mean(O)[inside]
print("craton cells with predicted HF: %d of %d; with any observation in cell: %d" % (np.isfinite(q_pred).sum(), inside.sum(), np.isfinite(q_obs).sum()))
def wq(x, ww, p):
    m = np.isfinite(x); x, ww = x[m], ww[m]; o = np.argsort(x); c = np.cumsum(ww[o]) / ww.sum()
    return x[o][np.searchsorted(c, p)]
print("\nCraton surface heat flow (mW/m2), area-weighted percentiles 10/50/90:")
print("  predicted (Lucazeau 2019 grid): %.0f / %.0f / %.0f" % tuple(wq(q_pred, w, p) for p in (.1, .5, .9)))
print("  observed cells only:            %.0f / %.0f / %.0f" % tuple(wq(q_obs, w, p) for p in (.1, .5, .9)))
ok = np.isfinite(q_pred); wk = w[ok]; qk = q_pred[ok]
Tcell = np.array([np.median([moho_T(qk[n], a, h, k, z[ok][n]) for a, h, k in PARAMS]) for n in range(len(qk))])
print("\nMoho temperature per craton cell from its own heat flow (C), area-weighted 10/50/90: %.0f / %.0f / %.0f" % tuple(wq(Tcell, wk, p) for p in (.1, .5, .9)))
print("\nShare of craton area by antigorite limit L and addable heat dT (window L-dT <= T <= L):")
for L in (525, 595, 625):
    print("  L = %d C: holds (T<=L) %.1f%%;  T>L %.1f%%" % (L, 100*wk[Tcell<=L].sum()/wk.sum(), 100*wk[Tcell>L].sum()/wk.sum()))
    for dT in (50, 100, 200):
        m = (Tcell <= L) & (Tcell >= L - dT)
        print("      dT <= %3d C: %5.1f%% of craton area" % (dT, 100*wk[m].sum()/wk.sum()))
# sensitivity: observed-only cells
ob = np.isfinite(q_obs)
if ob.sum():
    To = np.array([np.median([moho_T(q_obs[ok_i], a, h, k, z[ok_i]) for a, h, k in PARAMS]) for ok_i in np.where(ob)[0]])
    wo = w[ob]
    print("\nObserved-only cells (%d): Moho T 10/50/90 = %.0f / %.0f / %.0f C; window L=525, dT<=100: %.1f%%" % (ob.sum(), *(wq(To, wo, p) for p in (.1,.5,.9)), 100*wo[(To<=525)&(To>=425)].sum()/wo.sum()))
