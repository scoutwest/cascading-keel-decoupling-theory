"""Section 34.3: how much of Pangea must the collapse lower, as a function of released water?
Water balance (Section 16.1): final sea level h_s = V_r/A_Earth + f * mean(land height after subsidence); every land cell must end
below h_s - 6.86 m. Land is lowered only as far as needed: a cell of pre-Flood height h ends at min(h, X), where X + 6.86 = h_s.
Solve X for each V_r. Pre-Flood heights (free parameter, Section 16.4): a Beta(a, b) distribution on [0, 610 m] with mean about 200 m.
Prints, per V_r: the cap X (m), the share of Pangea that must be lowered, the mean and maximum drop in that share."""
import numpy as np

A_E = 5.10e8                 # km2
f = 2.105e8 / A_E            # Pangea fraction of Earth
HMAX = 610.0                 # m
CUB = 6.86                   # m
rng = np.random.default_rng(1656)

def solve(h, Vr_km3):
    Vr = Vr_km3 / A_E * 1000.0                      # m of sea level over the whole Earth
    lo, hi = 0.0, HMAX
    for _ in range(80):                              # bisection on X
        X = 0.5 * (lo + hi)
        lhs = X + CUB
        rhs = Vr + f * np.minimum(h, X).mean()
        if lhs > rhs: hi = X
        else: lo = X
    return 0.5 * (lo + hi)

print("Pangea fraction of Earth f = %.4f; sea level per 1e8 km3 = %.0f m" % (f, 1e8 / A_E * 1000))
for a, b in ((1.5, 3.0), (1.0, 2.0), (2.0, 4.0)):
    h = HMAX * rng.beta(a, b, 400000)
    print("\nBeta(%.1f, %.1f): mean land height %.0f m, share of land above 400 m %.0f%%" % (a, b, h.mean(), 100 * (h > 400).mean()))
    print("  %-10s %-10s %-22s %-14s %-10s" % ("V_r km3", "cap X (m)", "share of Pangea lowered", "mean drop (m)", "max drop (m)"))
    for Vr in (0.0, 2.26e7, 7.3e7, 1.3e8, 2.3e8):
        X = solve(h, Vr)
        m = h > X
        print("  %-10.2e %-10.0f %-22.1f %-14.0f %-10.0f" % (Vr, X, 100 * m.mean(), (h[m] - X).mean() if m.any() else 0.0, (h.max() - X) if m.any() else 0.0))
