"""Arithmetic behind Sections 14-18 of the model (draft 6).

Pure arithmetic: no input files, no plotting. Run:  python section_arithmetic.py
Every table printed here is quoted in the model; see the section named in each heading.
Inputs marked 'assumption' are placeholders (see Section 14 of the model).
"""
import math

G = 9.81
R = 8.314

# ---- Section 14 (Phase A) -------------------------------------------------
def mckenzie_initial_subsidence(beta, a=125e3, tc=35e3, rm=3330.0, rc=2800.0, al=3.28e-5, Tm=1333.0, rw=1030.0):
    """McKenzie (1978) initial subsidence (m), water-loaded, textbook form.
    v2's code applied (1 - al*Tm) twice in the denominator (overstated by 4.6%)."""
    t1 = (rm - rc) * (tc / a) * (1 - al * Tm * tc / (2 * a))
    t2 = al * Tm * rm / 2
    return a * (t1 - t2) * (1 - 1 / beta) / (rm * (1 - al * Tm) - rw)


print("== Section 14: McKenzie initial subsidence, 35 km crust, 125 km lithosphere ==")
for b in (1.2, 1.5, 2.0):
    print(f"beta {b}: {mckenzie_initial_subsidence(b):.0f} m")

print("\n== Section 14: ocean heat capacity and budget B (assumptions: rho 1030, cp 3990, 1.335e18 m3) ==")
m_ocean = 1.335e18 * 1030
cap = m_ocean * 3990
print(f"heat capacity {cap:.3e} J/K; warming 15->40 C: {cap * 25:.3e} J")

print("\n== Section 14: sea-level equivalent of released water ==")
A_EARTH = 5.10e8  # km2
for V in (2.26e7, 1.3e8, 2.3e8, 1e8):
    print(f"{V:.2e} km3 -> {V / A_EARTH * 1000:.0f} m of sea level")

# ---- Section 15.2 ----------------------------------------------------------
print("\n== Section 15.2: layer needed before collapse (thickness loss = d/(1+d)) ==")
print("density rise | thickness loss | layer for 0.91 km | layer for 0.48 km")
for d in (0.10, 0.15, 0.20):
    loss = d / (1 + d)
    print(f"{d:>11.0%} | {loss:>13.1%} | {0.91 / loss:>14.1f} km | {0.48 / loss:>14.1f} km")

print("\n== Section 15.2: heat of sinking the crust by 0.91 km (35 km crust, 2800 kg/m3) ==")
for area, name in ((1.73e14, "continental"), (2.1e14, "continental + transitional"), (4.08e13, "cratons")):
    E = 2800 * 35e3 * G * 910 * area
    print(f"{name}: {E:.2e} J = {E / 1.4e26:.2f} x B")

# ---- Section 16 ------------------------------------------------------------
print("\n== Section 16: subsidence needed (m) ==")
f = 2.105e8 / A_EARTH      # Pangea fraction of Earth
w = 3.8e7 / 2.105e8        # margin share of Pangea
P = 610 + 6.86             # peak plus 15 cubits (0.457 m, assumption)
S = 1445.0                 # beta = 1.5 margin subsidence (m)
print(f"f = {f:.4f}, w = {w:.3f}")
print("water (km3)     | water (m) | A all drops | B peaks only | C margins stretched, uniform | D margins stretched, peaks")
for V, name in ((0, "none"), (2.26e7, "groundwater"), (1.3e8, "all shallow"), (2.3e8, "present-depth")):
    va = V / A_EARTH * 1000
    A_ = (P - va - f * 200) / (1 - f)
    B_ = P - (va + f * 200)
    C_ = (P - va - f * (200 - w * S)) / (1 - f * (1 - w))
    D_ = P - (va + f * ((1 - w) * 200 + w * (200 - S)))
    print(f"{name:<15} | {va:>9.0f} | {A_:>11.0f} | {B_:>12.0f} | {C_:>28.0f} | {D_:>26.0f}")

# ---- Section 18 ------------------------------------------------------------
print("\n== Section 18: Arrhenius slowdown relative to 820 C ==")
T0 = 820 + 273.15
print("T (C) | E=20 | E=60 | E=90 | E=100 kJ/mol")
for Tc in (600, 550, 500):
    row = [math.exp(E * 1000 / R * (1 / (Tc + 273.15) - 1 / T0)) for E in (20, 60, 90, 100)]
    print(f"{Tc:>5} | " + " | ".join(f"x{v:.1f}" for v in row))
print("\nScaled time (days) for base times 7, 14, 90 days (window = 150 days)")
for Tc in (600, 550, 500):
    for E in (60, 90):
        f_ = math.exp(E * 1000 / R * (1 / (Tc + 273.15) - 1 / T0))
        print(f"{Tc} C, E={E}: " + ", ".join(f"{b * f_:.0f}" for b in (7, 14, 90)))

print("\n== Sections 15-19: pressure at depth (3000 kg/m3 uniform; approximate) ==")
for z in (35, 38.3, 42, 45, 50):
    print(f"{z} km: {3000 * G * z * 1e3 / 1e9:.2f} GPa")
