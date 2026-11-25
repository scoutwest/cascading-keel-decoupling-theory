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

# ---- Stoichiometric version for the low-pressure reaction (added 2026-11-24) ----
# Reaction (1) of Perrillat et al. 2005, as quoted in their introduction (pressure below 2 GPa):
#   antigorite = 18 forsterite + 4 talc + 27 H2O      (antigorite Mg48Si34O85(OH)62)
# Molar masses (g/mol) from atomic masses; densities (g/cm3) are from memory: antigorite 2.55-2.60,
# forsterite 3.22, talc 2.75-2.80, free water 1.0-1.1 at ~1 GPa and ~600 C (assumed).
M_O, M_H, M_Mg, M_Si = 15.999, 1.008, 24.305, 28.086
M_ant = 48*M_Mg + 34*M_Si + 85*M_O + 62*(M_O + M_H)
M_fo = 2*M_Mg + M_Si + 4*M_O
M_tc = 3*M_Mg + 4*M_Si + 10*M_O + 2*(M_O + M_H)
M_w = 2*M_H + M_O
assert abs(M_ant - (18*M_fo + 4*M_tc + 27*M_w)) < 0.5      # the reaction balances
print(f"\nantigorite {M_ant:.1f} g/mol; water released by reaction (1): {27*M_w:.1f} g = {100*27*M_w/M_ant:.1f} wt%"
      f" (all structural water, 31 H2O: {100*31*M_w/M_ant:.1f} wt%)")
for rho_a in (2.55, 2.60):
    for rho_t in (2.75, 2.80):
        for rho_w in (1.0, 1.1):
            v0 = M_ant / rho_a
            vs = 18*M_fo/3.22 + 4*M_tc/rho_t
            vw = 27*M_w/rho_w
            print(f"ant {rho_a} talc {rho_t} water {rho_w}: solid loss {100*(1-vs/v0):4.1f}%  net if water stays {100*(1-(vs+vw)/v0):+5.1f}%")
