"""Volume change on serpentinite breakdown (Section 21). Inputs are assumptions, stated below."""
# Densities in g/cm3. Antigorite ~2.55-2.60 and 13 wt% water are literature-range values
# (water content seen in search results; density from memory, unverified). Peridotite products ~3.2-3.3.
for rho_serp in (2.55, 2.60):
    for wt in (0.12, 0.13):          # water released, mass fraction
        for rho_solid in (3.20, 3.30):
            for rho_w in (1.0, 1.1): # fluid density at ~1 GPa, 600 C (assumed range)
                v0 = 1 / rho_serp
                v_solid = (1 - wt) / rho_solid
                v_fluid = wt / rho_w
                solid_loss = 1 - v_solid / v0
                net_if_fluid_stays = 1 - (v_solid + v_fluid) / v0
                print(f"serp {rho_serp} wt {wt} solid {rho_solid} water {rho_w}: "
                      f"solid skeleton loss {solid_loss*100:5.1f}%  net if fluid stays {net_if_fluid_stays*100:5.1f}%")
# Layer thickness needed for 0.91 km and 0.48 km of subsidence (Section 4.8) at 25-31% loss
for s in (0.91, 0.48):
    print(s, "km needs", [round(s / f, 1) for f in (0.25, 0.28, 0.31)], "km of serpentinite at 25/28/31% loss")
