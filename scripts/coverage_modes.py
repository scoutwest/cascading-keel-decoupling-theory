"""Section 34: how much continental crust could a collapse reach?
Two modes, both on CRUST1.0 1-degree cells with the Lucazeau 2019 grid heat flow and the Section 24.7 geotherm:
  S  serpentinite window: Moho T between L - dT and L (L = 525 C, dT = 100 C), layer assumed to exist.
  G  garnet-granulite densification of the lower crust (dry, no water vented): Moho T >= TG and Moho P >= PG
     (assumed boundary, see Section 34.1) and T <= TMAX (conductive geotherm not trusted above).
     Subsidence available = lower-crust thickness x mafic fraction x 10% (garnet granulite is about 10% denser, Section 19.2).
Continental crust = CRUST1.0 cells with crustal thickness >= 25 km (area compared with Pangea's 2.105e8 km2).
Prints area shares and the share where the density mode can supply the 0.48 or 0.91 km the flood needs."""
import os, csv, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("FLOOD_MODELS_ROOT", os.path.normpath(os.path.join(HERE, "..", "..")))
bnds = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.bnds")
rho = np.loadtxt(ROOT + "/Noah's Pangea v1/data/reference/crust1/crust1.rho")
hf = os.environ.get("HEATFLOW_GRID", os.path.join(HERE, "..", "data", "heatflow", "2019GC008389-sup-0003-Data_Set_SI-S01", "HFgrid14.csv"))
lat = 89.5 - np.arange(180)
LAT = np.repeat(lat, 360)
w = np.cos(np.radians(LAT))
AREA_EARTH = 5.10e8
cell_area = w / w.sum() * AREA_EARTH

G_ = np.full((360, 720), np.nan)
with open(hf) as fh:
    r = csv.reader(fh, delimiter=";"); next(r)
    for lo, la, p, s, o in r:
        G_[int((89.75 - float(la)) / 0.5 + 0.5), int((float(lo) + 179.75) / 0.5 + 0.5)] = float(p)
q = np.nanmean(G_.reshape(180, 2, 360, 2), axis=(1, 3)).ravel()

moho = bnds[:, 0] - bnds[:, 8]                       # km below the surface
crust_thick = bnds[:, 5] - bnds[:, 8]                # top of upper crust to Moho (km)
lower = bnds[:, 7] - bnds[:, 8]                      # lower-crust thickness (km)
P = np.zeros(len(LAT))
for i in range(8):
    P += rho[:, i] * 1000 * 9.81 * np.clip(bnds[:, i] - bnds[:, i + 1], 0, None) * 1000
P /= 1e9

cont = (crust_thick >= 25) & np.isfinite(q)
print("continental cells: %d; area %.3e km2 (Pangea 2.105e8)" % (cont.sum(), cell_area[cont].sum()))

T0 = 15.0
def moho_T(q0, a0, hr, k, zk):
    q0 = q0 * 1e-3; a0 = a0 * 1e-6; hr = hr * 1e3; zz = zk * 1e3
    return T0 + (q0 - a0 * hr) * zz / k + a0 * hr ** 2 / k * (1 - np.exp(-zz / hr))
HP = tuple(float(x) for x in os.environ.get('HEAT_PRODUCTION', '0.3,0.45,0.6').split(','))
PAR = list(itertools.product(HP, (8., 10., 12.), (2.5, 3.0)))
print('heat production set (uW/m3):', HP)
T = np.full(len(LAT), np.nan)
ix = np.where(cont)[0]
T[ix] = np.median([moho_T(q[ix], a, h, k, moho[ix]) for a, h, k in PAR], axis=0)

L, dT, TG, PG, TMAX = 525.0, 100.0, 600.0, 1.0, 900.0
A = cell_area[cont].sum()
def share(mask): return 100 * cell_area[mask & cont].sum() / A
S = (T >= L - dT) & (T <= L)
cold = T < L - dT
above = T > L
print("\nMoho temperature over continental crust (area-weighted percentiles 10/50/90): ", end="")
o = np.argsort(T[ix]); c = np.cumsum(cell_area[ix][o]) / A
print(" / ".join("%.0f" % T[ix][o][np.searchsorted(c, p)] for p in (.1, .5, .9)))
print("S window: %.1f%% | cold (< %.0f C): %.1f%% | above the antigorite limit: %.1f%% of continental area" % (share(S), L - dT, share(cold), share(above)))
print("Pangea-equivalent area of S: %.2e km2 (%.1f%% of 2.105e8)" % (cell_area[S & cont].sum(), 100 * cell_area[S & cont].sum() / 2.105e8))

Gm = (T >= TG) & (T <= TMAX) & (P >= PG)
print("\nG candidate (T %.0f-%.0f C and Moho P >= %.2f GPa): %.1f%% of continental area" % (TG, TMAX, PG, share(Gm)))
print("Moho pressure over continental crust: median %.2f GPa; share >= 1.0 GPa %.0f%%, >= 1.2 GPa %.0f%%" %
      (float(np.median(P[ix])), share(P >= 1.0), share(P >= 1.2)))
print("Lower-crust thickness over continental crust: median %.1f km" % float(np.median(lower[ix])))
print("\nSubsidence available from G = lower-crust thickness x mafic fraction x 10%:")
for need, lab in ((0.48, "0.48 km (1.3e8 km3 shallow water)"), (0.91, "0.91 km (no water)")):
    for f in (1.0, 0.5, 0.3):
        ok = Gm & (lower * f * 0.10 >= need)
        print("  need %-36s mafic fraction %.1f: G supplies it over %.1f%% of continental area" % (lab, f, share(ok)))
print("\nUnion for the 0.48 km case (S window, or G with mafic fraction f):")
for f in (1.0, 0.5, 0.3):
    u = S | (Gm & (lower * f * 0.10 >= 0.48))
    print("  f = %.1f: S or G covers %.1f%% of continental area (%.2e km2)" % (f, share(u), cell_area[u & cont].sum()))
print("\nSensitivity of G to its assumed boundary (union, f = 0.5, need 0.48 km):")
for tg, pg in ((550, 0.9), (600, 1.0), (650, 1.0), (700, 1.1), (750, 1.2)):
    g_ = (T >= tg) & (T <= TMAX) & (P >= pg) & (lower * 0.5 * 0.10 >= 0.48)
    print("  T >= %d C, P >= %.1f GPa: G alone %.1f%%; with S %.1f%%" % (tg, pg, share(g_), share(g_ | S)))

# Section 34.3: the collapse need only lower the highest ground (flood_fraction_needed.py), so the needed drop is smaller.
print("\nShare of continental area where S (window) or G supplies a given drop (heat production set above):")
print("  %-34s %-9s %-9s %-9s" % ("drop needed at the highest ground", "f = 1.0", "f = 0.5", "f = 0.3"))
for need, lab in ((0.08, "0.08 km (2.3e8 km3 water)"), (0.28, "0.28 km (1.3e8 km3 water)"), (0.41, "0.41 km (7.3e7 km3 water)"), (0.48, "0.48 km (uniform, 1.3e8 km3)"), (0.52, "0.52 km (net water about 3e7 km3)"), (0.91, "0.91 km (no water)")):
    row = []
    for f in (1.0, 0.5, 0.3):
        u = S | (Gm & (lower * f * 0.10 >= need))
        row.append("%.1f%%" % share(u))
    print("  %-34s %-9s %-9s %-9s" % (lab, *row))

# Is the capable area where the high ground is? Present elevation as an isostatic proxy (physics, not age): capable share by elevation band.
elev = bnds[:, 0] * 1000.0
cap = S | (Gm & (lower * 0.5 * 0.10 >= 0.28))
o = np.argsort(elev[ix])[::-1]; cs = np.cumsum(cell_area[ix][o]) / A
print("\nShare of each present-elevation band of continental crust that S or G (f = 0.5, 0.28 km) can lower:")
edges = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
for lo_, hi_ in zip(edges[:-1], edges[1:]):
    sel = ix[o][(cs > lo_) & (cs <= hi_)]
    print("  elevation rank %3d-%3d%% (highest first): median elevation %5.0f m, capable %.0f%% of its area" % (100 * lo_, 100 * hi_, np.median(elev[sel]), 100 * cell_area[sel][cap[sel]].sum() / cell_area[sel].sum()))
