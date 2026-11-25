"""Section 31: per-craton table. For each of the 114 Hasterok craton polygons: area, CRUST1.0 median Moho depth,
NGHF measured heat flow (median of points inside), Lucazeau grid heat flow (median of the 0.5 degree cells inside),
and the Moho temperature from the Section 24.7 geotherm (heat production 0.3-0.6; decay depth and conductivity swept;
median). The temperature uses the measured median where there are 3 or more points, else the grid median.
Window: T <= L and T >= L - dT (L = 525 C, dT = 100 C)."""
import os, csv
import itertools
import numpy as np
import shapefile
from shapely.geometry import shape, Point
from shapely.prepared import prep
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(HERE, "..", "..")))
D = os.path.join(HERE, "..", "data", "heatflow")
sf = shapefile.Reader(ROOT + "/Noah's Pangea v2/data/hasterok_2022/cratons.shp", encoding="latin-1")
geoms = [shape(s.__geo_interface__) for s in sf.shapes()]; recs = sf.records()
bnds = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.bnds")
moho = (bnds[:, 0] - bnds[:, 8]).reshape(180, 360)
rows = list(csv.DictReader(open(os.path.join(D, "NGHF.csv"), encoding="utf-8", errors="replace")))
kq = [k for k in rows[0] if k.startswith("heat-flow")][0]
def f(x):
    try: return float(x)
    except Exception: return np.nan
plon = np.array([f(r["longitude"]) for r in rows]); plat = np.array([f(r["latitude"]) for r in rows])
pq = np.array([f(r[kq]) for r in rows])
G = np.full((360, 720), np.nan)
with open(os.path.join(D, "2019GC008389-sup-0003-Data_Set_SI-S01", "HFgrid14.csv")) as fh:
    r = csv.reader(fh, delimiter=";"); next(r)
    for lo, la, p, s, o in r:
        G[int((89.75 - float(la)) / 0.5 + 0.5), int((float(lo) + 179.75) / 0.5 + 0.5)] = float(p)
T0 = 15.0
def moho_T(q0, a0, hr, k, zk):
    q0 *= 1e-3; a0 *= 1e-6; hr *= 1e3; zz = zk * 1e3
    return T0 + (q0 - a0 * hr) * zz / k + a0 * hr ** 2 / k * (1 - np.exp(-zz / hr))
P = list(itertools.product((0.3, 0.45, 0.6), (8., 10., 12.), (2.5, 3.0)))
def Tm(q, z): return float(np.median([moho_T(q, a, h, k, z) for a, h, k in P]))
L, dT = 525.0, 100.0
out = []
for g, rc in zip(geoms, recs):
    pg = prep(g); minx, miny, maxx, maxy = g.bounds
    sel = np.where((plon >= minx) & (plon <= maxx) & (plat >= miny) & (plat <= maxy) & np.isfinite(pq))[0]
    qm = np.array([pq[i] for i in sel if pg.contains(Point(plon[i], plat[i]))])
    gq = []; zz = []
    for la in np.arange(np.floor(miny * 2) / 2 + 0.25, maxy + 0.5, 0.5):
        for lo in np.arange(np.floor(minx * 2) / 2 + 0.25, maxx + 0.5, 0.5):
            if -90 < la < 90 and -180 < lo < 180 and pg.contains(Point(lo, la)):
                gq.append(G[int((89.75 - la) / 0.5 + 0.5), int((lo + 179.75) / 0.5 + 0.5)])
                zz.append(moho[min(179, int(90 - la)), min(359, int(lo + 180))])
    gq = np.array(gq); zz = np.array(zz)
    area = float(rc["area"]) if rc["area"] else np.nan
    if len(gq) == 0:
        out.append((rc["prov_name"], area, np.nan, len(qm), np.nan, np.nan, np.nan, "no cells")); continue
    zm = float(np.median(zz)); gm = float(np.nanmedian(gq))
    if len(qm) >= 3: qu, basis = float(np.median(qm)), "measured"
    else: qu, basis = gm, "grid"
    T = Tm(qu, zm)
    out.append((rc["prov_name"], area, zm, len(qm), (float(np.median(qm)) if len(qm) else np.nan), gm, T, basis))
rows_ok = [o for o in out if o[7] != "no cells"]
print("%d polygons; %d too small for a 0.5 degree cell centre (no grid cells)" % (len(out), len(out) - len(rows_ok)))
print("\n%-34s %9s %6s %4s %9s %6s %6s  %s" % ("craton", "area km2", "Moho", "n", "meas q", "grid q", "T C", "basis / class"))
for o in sorted(rows_ok, key=lambda o: o[6]):
    T = o[6]
    cls = "hold, window" if (T <= L and T >= L - dT) else ("hold, cold" if T < L - dT else "above limit")
    print("%-34s %9.2e %6.1f %4d %9s %6.0f %6.0f  %s / %s" % (str(o[0])[:34], o[1], o[2], o[3], ("%.0f" % o[4]) if o[3] else "-", o[5], T, o[7], cls))
Ts = np.array([o[6] for o in rows_ok]); ar = np.array([o[1] for o in rows_ok])
print("\nBy count: window %d, cold (T < %.0f) %d, above limit %d" % (np.sum((Ts <= L) & (Ts >= L - dT)), L - dT, np.sum(Ts < L - dT), np.sum(Ts > L)))
print("By area:  window %.1f%%, cold %.1f%%, above limit %.1f%%" % tuple(100 * ar[m].sum() / ar.sum() for m in ((Ts <= L) & (Ts >= L - dT), Ts < L - dT, Ts > L)))
