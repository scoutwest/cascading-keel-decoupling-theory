"""Section 30: resolve the heat-flow disagreement (Section 29) with the NGHF point measurements.
1. Measured heat flow (NGHF, Lucazeau 2019 S02) at points inside the 114 craton polygons, versus the grid prediction there.
2. Grid and point values at the two xenolith sites (Udachnaya, Siberia; Tanzanian craton).
3. Measured heat production and conductivity inside cratons.
4. Moho temperature from measured heat flow."""
import os, csv, itertools
import numpy as np
import shapefile
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
from shapely.prepared import prep
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(HERE, "..", "..")))
D = os.path.join(HERE, "..", "data", "heatflow")
sf = shapefile.Reader(ROOT + "/Noah's Pangea v2/data/hasterok_2022/cratons.shp")
geoms = [shape(s.__geo_interface__) for s in sf.shapes()]; prepped = [prep(g) for g in geoms]; tree = STRtree(geoms)
def craton_of(lon, lat):
    pt = Point(lon, lat)
    for i in tree.query(pt):
        if prepped[int(i)].contains(pt): return int(i)
    return -1
rows = list(csv.DictReader(open(os.path.join(D, "NGHF.csv"), encoding="utf-8", errors="replace")))
keys = list(rows[0].keys())
kq = [k for k in keys if k.startswith("heat-flow")][0]; kh = [k for k in keys if k.startswith("heat production (")][0]
kc = [k for k in keys if k.startswith("conductivity (")][0]; ku = [k for k in keys if k.startswith("uncertainty")][0]
def f(x):
    try: return float(x)
    except: return np.nan
lon = np.array([f(r["longitude"]) for r in rows]); lat = np.array([f(r["latitude"]) for r in rows])
q = np.array([f(r[kq]) for r in rows]); hp = np.array([f(r[kh]) for r in rows]); kk = np.array([f(r[kc]) for r in rows])
ok = np.isfinite(q) & np.isfinite(lon) & np.isfinite(lat)
cid = np.full(len(rows), -2)
for n in np.where(ok)[0]: cid[n] = craton_of(lon[n], lat[n])
inc = cid >= 0
print("NGHF records: %d; with heat flow %d; inside craton polygons %d (%.1f%%)" % (len(rows), ok.sum(), inc.sum(), 100*inc.sum()/ok.sum()))
# grid
G = np.full((360, 720), np.nan)
with open(os.path.join(D, "2019GC008389-sup-0003-Data_Set_SI-S01", "HFgrid14.csv")) as fh:
    r = csv.reader(fh, delimiter=";"); next(r)
    for lo, la, p, s, o in r:
        G[int((89.75-float(la))/0.5+0.5), int((float(lo)+179.75)/0.5+0.5)] = float(p)
def grid_at(lo, la): return G[int((89.75-la)/0.5+0.5), int((lo+179.75)/0.5+0.5)]
def pct(x): x = x[np.isfinite(x)]; return "n=%d  10/50/90 = %.0f / %.0f / %.0f" % (len(x), *np.percentile(x, [10, 50, 90]))
print("\nHeat flow at measurement points inside cratons (mW/m2, unweighted):  " + pct(q[inc]))
gp = np.array([grid_at(lon[n], lat[n]) for n in np.where(inc)[0]])
print("Grid prediction at those same points:                                " + pct(gp))
print("Grid prediction minus measured, median: %.1f mW/m2" % np.nanmedian(gp - q[inc]))
# per-craton medians, only cratons with >= 3 measurements
ids = np.unique(cid[inc]); med = []; ct = []
for i in ids:
    m = cid == i
    if m.sum() >= 3: med.append(np.median(q[m])); ct.append(m.sum())
print("Per-craton median measured heat flow, cratons with >=3 points (%d of %d cratons with any): %s" % (len(med), len(ids), pct(np.array(med))))
print("Cratons with no measurement at all: %d of %d" % (len(geoms) - len(ids), len(geoms)))
print("\nMeasured heat production inside cratons (uW/m3): " + pct(hp[inc]))
print("Measured conductivity inside cratons (W/m/K):   " + pct(kk[inc]))
# the two xenolith sites
sites = {"Udachnaya (Siberia)": (112.3, 66.4), "Tanzanian craton (Lashaine/Labait)": (36.0, -3.0), "Tanzanian craton (west, ~Lake Victoria)": (33.0, -3.5)}
print("\nGrid and nearby measurements at the xenolith sites (within 1.5 degrees):")
for name, (lo, la) in sites.items():
    d = np.hypot((lon - lo) * np.cos(np.radians(la)), lat - la); m = ok & (d < 1.5)
    print("  %-42s grid %.0f mW/m2; %d points, median %s, range %s" % (name, grid_at(lo, la), m.sum(),
          ("%.0f" % np.median(q[m])) if m.sum() else "-", ("%.0f-%.0f" % (q[m].min(), q[m].max())) if m.sum() else "-"))
# Moho T from measured point heat flow, using CRUST1 Moho depth 39 km and the Section 24.7 parameter sweep
T0 = 15.0
def moho_T(q0, a0, hr, k, zk):
    q0 *= 1e-3; a0 *= 1e-6; hr *= 1e3; zz = zk * 1e3
    return T0 + (q0 - a0*hr)*zz/k + a0*hr**2/k*(1 - np.exp(-zz/hr))
P = list(itertools.product((0.3, 0.45, 0.6), (8., 10., 12.), (2.5, 3.0)))
def Tmed(qv): return np.array([np.median([moho_T(x, a, h, k, 39.0) for a, h, k in P]) for x in qv])
print("\nMoho temperature at 39 km from heat flow (assumed heat production 0.3-0.6 uW/m3):")
for lab, qv in (("measured points in cratons", q[inc]), ("per-craton medians", np.array(med)), ("grid at those points", gp)):
    T = Tmed(qv[np.isfinite(qv)]); print("  %-28s T 10/50/90 = %.0f / %.0f / %.0f C;  share <=525 C: %.0f%%; <=625 C: %.0f%%" % (lab, *np.percentile(T, [10, 50, 90]), 100*np.mean(T <= 525), 100*np.mean(T <= 625)))
print("\nHeat flow that gives each xenolith Moho T at 39 km (median over the sweep): 400 C -> %.0f; 500 C -> %.0f; 600 C -> %.0f mW/m2" %
      tuple(np.interp(t, Tmed(np.arange(15, 90, 1.0)), np.arange(15, 90, 1.0)) for t in (400, 500, 600)))

# Sensitivity to crustal heat production (the xenolith study's 0.3-0.6 against the continental average 0.8-1.0, Section 24.6)
print("\nSensitivity of the per-craton-median Moho T to heat production (decay depth and conductivity swept as before):")
medq = np.array(med)
for lab, hps in (("0.3-0.6 (xenolith)", (0.3, 0.45, 0.6)), ("0.6-0.9", (0.6, 0.75, 0.9)), ("0.8-1.0 (continental mean)", (0.8, 0.9, 1.0))):
    PP = list(itertools.product(hps, (8., 10., 12.), (2.5, 3.0)))
    T = np.array([np.median([moho_T(x, a, h, k, 39.0) for a, h, k in PP]) for x in medq])
    print("  heat production %-26s Moho T 10/50/90 = %.0f / %.0f / %.0f C; share <=525 C: %.0f%%; <=625 C: %.0f%%" % (lab, *np.percentile(T, [10, 50, 90]), 100*np.mean(T <= 525), 100*np.mean(T <= 625)))
