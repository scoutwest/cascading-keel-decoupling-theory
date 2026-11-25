"""Section 27: state of water vented from the crust base, using the IAPWS-95 equation of state (iapws package).
Water at the Moho (P0, T0) rises through a crack with no heat loss to the walls (the crack-wall cooling is 0.08%, Section 26)
and no kinetic energy at the exit, losing only gravitational potential energy g*dz. The enthalpy at the surface is then
h_surf = h0 - g*dz, and the state follows at the surface pressure. IAPWS-95 is valid to 1000 MPa and 1273 K, so the start
pressure is capped at 1000 MPa (1 GPa). Heat to remove to bring it to liquid water at 40 C = (h_surf - h_liq(40 C)).
Budget B = 1.4e26 J (Section 26, unverified)."""
from iapws import IAPWS95
g, dz = 9.81, 35e3
B = 1.4e26
KG = lambda Vkm3: Vkm3 * 1e9 * 1000.0           # km3 liquid -> kg
liq40 = IAPWS95(T=313.15, P=0.1013)             # liquid at 40 C, surface pressure
print(f"liquid water at 40 C, 0.1 MPa: h = {liq40.h:.1f} kJ/kg")
print(f"gravitational potential term g*dz = {g*dz/1e3:.0f} kJ/kg")
print(f"{'start':>14s} {'h0 kJ/kg':>9s} {'h_surf':>8s} {'T_surf C':>9s} {'steam frac':>10s} {'heat to 40C kJ/kg':>18s} {'x budget B (1.3e8 km3)':>24s}")
for P0 in (1000.0, 600.0, 300.0):
    for T0c in (400, 500, 580, 700):
        s0 = IAPWS95(T=T0c + 273.15, P=P0)
        h0 = s0.h
        hs = h0 - g * dz / 1e3
        sf = IAPWS95(P=0.1013, h=hs)
        Q = KG(1.3e8) * (hs - liq40.h) * 1e3
        print(f"{P0:6.0f} MPa {T0c:3d}C {h0:9.0f} {hs:8.0f} {sf.T-273.15:9.0f} {sf.x:10.2f} {hs-liq40.h:18.0f} {Q/B:24.2f}")
# constant-cp comparison (Section 26)
print("\nConstant-cp estimate used in Section 26 for 580 C:", round(4.186*(580-40)), "kJ/kg  ->", round(KG(1.3e8)*4186*(580-40)/B,2), "x B")

# ---- where the heat goes (first order; Earth values from memory) ----
CP_AIR, M_AIR, A_EARTH, F_OUT = 1005.0, 5.15e18, 5.1e14, 240.0
s = IAPWS95(T=580 + 273.15, P=1000.0); hs = s.h - g * dz / 1e3
sf = IAPWS95(P=0.1013, h=hs)
Q = KG(1.3e8) * (hs - liq40.h) * 1e3
latent = sf.x * 2257.0 / (hs - liq40.h)
print(f"\n580 C start, 1.3e8 km3: heat to shed {Q:.2e} J; steam fraction {sf.x:.2f}; latent heat about {100*latent:.0f}% of it")
print(f"atmosphere heat capacity {CP_AIR*M_AIR:.1e} J/K -> the heat would warm it by {Q/(CP_AIR*M_AIR):.0f} K if it stayed there")
print(f"steam mass {sf.x*KG(1.3e8):.1e} kg (present atmospheric water vapour is about 1.3e16 kg)")
print(f"radiation to space at present outgoing flux {F_OUT} W/m2 over the globe: {F_OUT*A_EARTH:.1e} W -> {Q/(F_OUT*A_EARTH)/3.156e7:.0f} years")
print(f"1.3e8 km3 spread as rain over the whole globe: {1.3e8*1e9/A_EARTH:.0f} m")
