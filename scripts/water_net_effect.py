"""Section 35: net sea-level effect of the shallow water reservoirs (O9).
The Section 5.1 budget holds the ocean floor fixed. Water that leaves rock under the ocean floor lets the floor settle, and the
basin gains as much volume as the sediment or rock loses. Only the difference raises sea level.
  groundwater (continental)    : net = gross (the land pore space does not become ocean basin).
  sediment pore water          : sediment compacts as the water leaves, basin gain = water lost, net = 0.
  serpentinite under the ocean : solid shrinks by s of the layer volume while the water released has volume w of it; net = (w - s)/w.
Then solve the flood geometry for the net water (as flood_fraction_needed.py) and print the drop and share of Pangea to be lowered."""
import numpy as np

A_E = 5.10e8; f = 2.105e8 / A_E; HMAX = 610.0; CUB = 6.86
GW = 2.26e7                                   # km3, Gleeson et al. 2016 (sourced)
SED_MODEL, SED_LAROWE = 6e7, 8.46e7           # km3: model assumption; LaRowe et al. 2017 (sourced)
SERP = 5e7                                    # km3, model assumption (5 km, 10%, 13 wt%)
RHO_ROCK, RHO_W = 2600.0, 1000.0
print("Serpentinite net fraction (w - s)/w, w = water volume per unit rock volume, s = solid shrink:")
nets = []
for wt in (0.107, 0.13):
    w = wt * RHO_ROCK / RHO_W
    for s in (0.23, 0.25):
        nets.append((w - s) / w)
        print("  water %.1f wt%%: w = %.3f, s = %.2f -> net fraction %.2f" % (100 * wt, w, s, (w - s) / w))
lo_n, hi_n = min(nets), max(nets)
print("  range of net fraction: %.2f-%.2f" % (lo_n, hi_n))

print("\nReservoir table (km3):")
print("  %-34s %-10s %-10s" % ("reservoir", "gross", "net"))
print("  %-34s %-10.2e %-10.2e" % ("continental groundwater (sourced)", GW, GW))
print("  %-34s %-10.2e %-10.2e" % ("sediment pore water (model 6e7)", SED_MODEL, 0))
print("  %-34s %-10.2e %-10.2e" % ("sediment pore water (LaRowe 8.46e7)", SED_LAROWE, 0))
print("  %-34s %-10.2e %.1e-%.1e" % ("oceanic serpentinite (assumed)", SERP, SERP * lo_n, SERP * hi_n))
gross_model = GW + SED_MODEL + SERP
net_lo = GW + SERP * lo_n; net_hi = GW + SERP * hi_n
print("  %-34s %-10.2e %.2e-%.2e" % ("TOTAL (model reservoirs)", gross_model, net_lo, net_hi))

rng = np.random.default_rng(1656)
def solve(h, Vr_km3):
    Vr = Vr_km3 / A_E * 1000.0
    lo, hi = 0.0, HMAX
    for _ in range(80):
        X = 0.5 * (lo + hi)
        if X + CUB > Vr + f * np.minimum(h, X).mean(): hi = X
        else: lo = X
    return 0.5 * (lo + hi)
print("\nFlood geometry (pre-Flood heights Beta(1.5, 3) on 0-610 m, mean about 200 m):")
h = HMAX * rng.beta(1.5, 3.0, 400000)
print("  %-34s %-12s %-24s %-12s" % ("water counted", "V_r (km3)", "share of Pangea lowered", "max drop (m)"))
for lab, Vr in (("gross model total (Section 5.1)", gross_model), ("net, upper (serpentinite 0.17)", net_hi), ("net, lower (serpentinite 0.10)", net_lo),
                ("groundwater only", GW), ("none", 0.0)):
    X = solve(h, Vr); m = h > X
    print("  %-34s %-12.2e %-24.1f %-12.0f" % (lab, Vr, 100 * m.mean(), h.max() - X))
print("\nUniform lowering (Section 16 case A), 334 m per 1e8 km3 from 910 m: gross %.0f m; net %.0f-%.0f m" %
      (910 - 334 * gross_model / 1e8, 910 - 334 * net_hi / 1e8, 910 - 334 * net_lo / 1e8))
