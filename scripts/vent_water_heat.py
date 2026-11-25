"""Section 26: heat carried by water vented from the crust base, against the heat budget (H5).
Inputs: water volume 1.3e8 km3 (model's shallow-water target, Section 5.1) or the amount from a 3 km serpentinite layer
(0.83 km3 per km2 over Pangea's 2.105e8 km2; 10.7 wt% water, 2.6 g/cm3); source temperature 500-580 C (craton Moho);
ocean cap 40 C; heat budget B = 1.4e26 J (heating the 1.335e9 km3 ocean from 15 to 40 C; Section 3.2, unverified inputs).
Constant specific heat 4186 J/kg/K (real hot water is higher, so heat is a lower bound). Rock diffusivity 1e-6 m2/s,
rock rho 2800 kg/m3 and cp 1000 J/kg/K, hot-water viscosity 1e-4 Pa s: all assumed from memory."""
import numpy as np
CP, B, TOC = 4186.0, 1.4e26, 40.0
LAND = 2.105e8                         # km2, Pangea
V_TARGET = 1.3e8                       # km3
V_LAYER = 0.83 * LAND                  # km3 from a 3 km serpentinite layer
kg = lambda Vkm3: Vkm3 * 1e9 * 1000.0

print("HEAT CARRIED UP AS HOT WATER")
for V in (V_TARGET, V_LAYER):
    for T in (300, 400, 500, 580):
        Q = kg(V) * CP * (T - TOC)
        print(f"  V {V:.2e} km3, source {T} C: {Q:.2e} J = {Q/B:.2f} x B")
print("\nMAXIMUM SOURCE TEMPERATURE IF ALL OF B IS AVAILABLE TO THE VENTED WATER")
for V in (V_TARGET, V_LAYER):
    print(f"  V {V:.2e} km3: {TOC + B/(kg(V)*CP):.0f} C")
print("\nMAXIMUM WATER VOLUME AT A GIVEN SOURCE TEMPERATURE (km3), and the extra collapse subsidence needed to make up the rest")
for T in (300, 400, 500, 580):
    Vmax = B / (CP * (T - TOC)) / 1e3 / 1e9 * 1e0   # kg of water -> m3 -> km3
    Vmax = B / (CP * (T - TOC)) / 1000.0 / 1e9      # m3 -> km3
    short = max(0.0, V_TARGET - Vmax)
    print(f"  {T} C: {Vmax:.2e} km3; shortfall {short:.2e} km3 -> extra subsidence {330*short/1e8:.0f} m (330 m per 1e8 km3, Section 12)")
print("\nCOOLING BY THE CRACK WALLS (150 days)")
kappa, t = 1e-6, 150 * 86400
d = np.sqrt(kappa * t)
spacing = 2 * 81.0                     # km between cracks (twice the median distance to a province boundary, Section 22)
trace = LAND / spacing * 1e3 * 35e3    # wall area, m2 (35 km tall)
Qwall = trace * 2 * d * 2800 * 1000 * 250
Q = kg(V_TARGET) * CP * (540)
print(f"  conduction length {d:.1f} m; wall heat capacity (250 K) {Qwall:.1e} J = {100*Qwall/Q:.2f}% of the heat to shed")
print("\nDIFFUSE FLOW THROUGH A FRACTURE NETWORK")
# rock that must be cooled by 250 K to absorb the heat of V_TARGET at 540 K excess
rockvol = Q / (2800 * 1000 * 250)      # m3
print(f"  rock cooled by 250 K needed: {rockvol/1e9:.1e} km3 = {rockvol/1e9/LAND:.1f} km thick over Pangea")
print(f"  full heat exchange needs fracture spacing of about {2*d:.0f} m or less")
q = 0.83e3 / (150 * 86400)             # Darcy flux, m/s (0.83 km of water per km2 over 150 days)
mu, grad = 1e-4, (2800 - 1000) * 9.81  # Pa s, Pa/m (buoyancy)
k = q * mu / grad
b = (12 * (2 * d) * k) ** (1 / 3)
print(f"  Darcy flux {q:.1e} m/s; permeability needed {k:.1e} m2 (typical crust 1e-18 to 1e-16)")
print(f"  parallel fractures {b*1e3:.2f} mm wide every {2*d:.0f} m would give that permeability")
