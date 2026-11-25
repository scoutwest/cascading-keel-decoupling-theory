"""Section 27: can a free water layer sit under the crust for the ~1,650 pre-Flood years?
Three first-order checks. ALL inputs are assumed from memory (not sourced): rock fracture toughness 1-3 MPa m^0.5,
hot-water viscosity 1e-4 Pa s, rock and water densities, lower-crust viscosity 1e20-1e23 Pa s, permeability 1e-19 to 1e-15 m2,
connected porosity 1-3%."""
import numpy as np
YR = 3.156e7
g, rho_r, rho_w = 9.81, 2800.0, 1100.0          # water at ~1 GPa, 500 C: ~1.0-1.2 g/cm3 (assumed)
drho = rho_r - rho_w
H = 35e3                                        # crust thickness, m

print("1) BRITTLE REGIME: how tall can a connected water body be before it drives a crack upward by buoyancy?")
for K in (1e6, 2e6, 3e6):
    Lc = (K / (drho * g)) ** (2 / 3)
    print(f"   fracture toughness {K/1e6:.0f} MPa m^0.5: critical height about {Lc:.0f} m  (buoyant pressure {drho*g*Lc/1e6:.1f} MPa)")
print(f"   (a free water layer 480-910 m thick is far taller than this, so in brittle rock it escapes upward at once)")

print("\n2) DUCTILE REGIME: Darcy ascent time through 35 km of crust for buoyant water, versus permeability")
mu = 1e-4
for k in (1e-19, 1e-18, 1e-17, 1e-16, 1e-15):
    for phi in (0.01, 0.03):
        v = k * drho * g / (mu * phi)           # m/s (porosity-wave speed order of magnitude)
        print(f"   k {k:.0e} m2, porosity {phi*100:.0f}%: ascent {v*YR:.2e} m/yr -> {H/v/YR:.2e} yr to cross 35 km")

print("\n3) DENSITY INVERSION: Rayleigh-Taylor sinking time of the crust into a water layer, versus crust viscosity")
for eta in (1e19, 1e20, 1e21, 1e22, 1e23):
    t = eta / (drho * g * H)                    # s, characteristic viscous overturn time
    print(f"   viscosity {eta:.0e} Pa s: about {t/YR:.1e} yr")
