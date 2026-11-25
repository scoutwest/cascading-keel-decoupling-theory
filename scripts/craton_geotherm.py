"""Steady-state conductive geotherm for a present-day craton crust (Section 24).
T(z) = T0 + q_m z / k + (A0 hr^2 / k) (1 - exp(-z/hr)),  q_m = q0 - A0 hr
(exponentially decreasing crustal heat production, constant conductivity; the Lachenbruch/Pollack-Chapman form).
Surface heat flow q0 (40-47 mW/m2) and crust thickness are from search results (Section 21.7).
Surface temperature, heat production A0, its decay depth hr and conductivity k are assumed, from memory, unsourced.
Sweeping them gives the spread; nothing is fitted."""
import itertools
import numpy as np

T0 = 15.0                                    # surface temperature, C (same start as the project)
Q0 = (40.0, 44.0, 47.0)                      # mW/m2 (search: 40-50; Archean mean 41; Kalahari 47 +/- 2)
A0 = (0.8, 1.2, 1.6)                         # surface heat production, uW/m3 (assumed)
HR = (8.0, 10.0, 12.0)                       # decay depth, km (assumed)
K = (2.5, 3.0)                               # conductivity, W/m/K (assumed)
ZM = (35.0, 39.0, 42.0)                      # Moho depth, km (35-40 from search; 39 and 42 from Section 19)

def moho_T(q0, a0, hr, k, z):
    q0 *= 1e-3; a0 *= 1e-6; hr *= 1e3; z *= 1e3
    qm = q0 - a0 * hr
    if qm < 0.008:                           # reduced heat flow below 8 mW/m2 is implausible; skip
        return None
    return T0 + qm * z / k + a0 * hr**2 / k * (1 - np.exp(-z / hr)), qm * 1e3

rows = []
for q0, a0, hr, k, z in itertools.product(Q0, A0, HR, K, ZM):
    r = moho_T(q0, a0, hr, k, z)
    if r:
        rows.append((q0, a0, hr, k, z, r[0], r[1]))
T = np.array([r[5] for r in rows])
print(f"{len(rows)} valid combinations (of {len(Q0)*len(A0)*len(HR)*len(K)*len(ZM)})")
print(f"Moho temperature, C: min {T.min():.0f}  5% {np.percentile(T,5):.0f}  25% {np.percentile(T,25):.0f}  "
      f"median {np.median(T):.0f}  75% {np.percentile(T,75):.0f}  95% {np.percentile(T,95):.0f}  max {T.max():.0f}")
print("\nBy surface heat flow and Moho depth (all A0, hr, k):")
for q0 in Q0:
    print(f"  q0 {q0:.0f} mW/m2: " + "  ".join(
        f"Moho {z:.0f} km: {np.median([r[5] for r in rows if r[0]==q0 and r[4]==z]):4.0f} C "
        f"({min(r[5] for r in rows if r[0]==q0 and r[4]==z):.0f}-{max(r[5] for r in rows if r[0]==q0 and r[4]==z):.0f})"
        for z in ZM))
print("\nBy conductivity (all else):")
for k in K:
    t = np.array([r[5] for r in rows if r[3]==k]); print(f"  k {k}: median {np.median(t):.0f} C, range {t.min():.0f}-{t.max():.0f}")
print("\nShare of valid cases with Moho temperature at or above:")
for lim in (500, 550, 600, 650, 680):
    print(f"  {lim} C: {100*(T>=lim).mean():.0f}%")
