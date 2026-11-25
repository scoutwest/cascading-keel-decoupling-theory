"""Upper bound on the volume fraction of serpentinite taken by trapped hydrogen at lower-crust pressure (Section 21).
H2 yield per kg of peridotite: 0.47 mol is the stoichiometric ceiling (10 wt% FeO, all Fe2+ -> magnetite);
modelled/measured values 0.1-0.35 mol (search results, Section 21.6). Molar volume of H2 at ~1 GPa, 500-600 C is NOT
sourced here: 10-20 cm3/mol is an assumed range from memory (ideal gas would give 5 cm3/mol, real hydrogen is much
less compressible at this pressure). Serpentinite density 2.6 g/cm3, 13 wt% water gained (12-13 wt% assumed)."""
rock_kg = 1.13                # 1 kg peridotite + 0.13 kg water
rock_cm3 = rock_kg * 1000 / 2.6
print(f"rock volume per kg peridotite as serpentinite: {rock_cm3:.0f} cm3")
for mol in (0.10, 0.35, 0.47):
    for vm in (10, 15, 20):
        f = mol * vm / rock_cm3
        print(f"H2 {mol:.2f} mol/kg, molar volume {vm:2d} cm3/mol: gas volume {mol*vm:5.1f} cm3 = {100*f:4.2f}% of rock")
fmax = 0.47 * 20 / rock_cm3
fmin = 0.10 * 10 / rock_cm3
print(f"range {100*fmin:.2f}% - {100*fmax:.2f}%")
for s in (0.91, 0.48):
    print(f"subsidence {s} km from gas alone needs a layer of {s/fmax:.0f} km (best case) to {s/fmin:.0f} km (low case)")
