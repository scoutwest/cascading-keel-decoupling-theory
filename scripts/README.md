# Scripts

Scripts that produce the images and the numbers quoted in `Cascading Keel Decoupling Theory (CKDT).md`. They are the working scripts from draft 6, not the Phase C rebuild; Phase C will replace them with a single tested package (see Section 11 of the model).

| Script | Produces | Inputs |
|--------|----------|--------|
| `section_arithmetic.py` | The tables in Sections 14-18 (McKenzie subsidence, water sea-level equivalents, layer thickness, the Section 16 budgets, the Arrhenius scaling, pressure at depth) | none |
| `craton_moho_stats.py` | The Section 19 statistics (craton Moho depth and pressure, lower-crust thickness) | CRUST1.0 and the Hasterok craton polygons |
| `serpentinite_volume.py` | The Section 21 volume arithmetic for serpentinite breakdown | none |
| `gas_volume.py` | The Section 21.3 upper bound on hydrogen volume in serpentinite | none |
| `craton_geotherm.py` | The Section 24 Moho temperature sweep (conductive geotherm; heat production and conductivity assumed) | none |
| `craton_geotherm_measured.py` | The Section 24.7 forward and inverse geotherm with the measured xenolith values | none |
| `craton_selection.py` | The Section 25 craton selection test (Moho temperature against surface heat flow, window against the antigorite limit) | CRUST1.0 and the Hasterok craton polygons |
| `craton_heatflow.py` | The Section 29 craton heat-flow sampling and Moho-temperature window | CRUST1.0, Hasterok craton polygons, and the Lucazeau 2019 grid `HFgrid14.csv` (`data/heatflow/` in this repository's folder; override with `HEATFLOW_GRID`) |
| `craton_heatflow_points.py` | The Section 30 comparison of NGHF point heat flow, the grid and the xenolith sites, and the Moho temperature spread | Hasterok craton polygons; `data/heatflow/NGHF.csv` and `HFgrid14.csv` |
| `craton_table.py` | The Section 31 per-craton Moho temperature table and window classes | Hasterok craton polygons, CRUST1.0, `data/heatflow/` |
| `coverage_modes.py` | The Section 34 coverage of the serpentinite-window and garnet-granulite modes over continental crust (set `HEAT_PRODUCTION=0.8,0.9,1.0` for the second run) | CRUST1.0 and `data/heatflow/` |
| `flood_fraction_needed.py` | The Section 34.3 share of Pangea that must be lowered, by released water | none |
| `vent_water_heat.py` | The Section 26 heat carried by vented water against the heat budget | none |
| `free_water_layer.py` | The Section 27.1 free-water-layer checks (buoyant crack height, ascent time, density-inversion time) | none |
| `vent_steam.py` | The Section 27.2-27.3 steam-state and heat calculation (IAPWS-95) | needs `iapws` (`pip install iapws`) |
| `province_suture_distance.py` | The Section 22 distances from continental crust to province boundaries (`python province_suture_distance.py strict` for different-type boundaries only) | Hasterok provinces (`global_gprv.shp`); needs `scipy` too |
| `draw_01_collapse_overview.py` | `images/01_collapse_overview` (.png, .svg) | none |
| `draw_02_craton_vs_margin.py` | `images/02_craton_vs_margin` | none |
| `draw_03_cascade_loop.py` | `images/03_cascade_loop` | none |
| `draw_04_craton_moho_pressure.py` | `images/04_craton_moho_pressure` | CRUST1.0 and the Hasterok craton polygons |

## Running

```
python section_arithmetic.py
python draw_01_collapse_overview.py        # writes into ../images
python craton_moho_stats.py                # prints the Section 19 statistics
```

Needs Python 3 with `numpy` and `matplotlib`. The two data scripts also need `pyshp` (`pip install pyshp`) and `shapely`.

## Data the data scripts need

The data are not in this repository. The scripts look for them in the folder that holds this repository's folder (override with the `FLOOD_MODELS_ROOT` environment variable):

- `Noah's Pangea v1/data/reference/crust1/crust1.bnds` and `crust1.rho` (CRUST1.0, Laske et al., 1°x1° grid; nine layers per cell).
- `Noah's Pangea v2/data/hasterok_2022/cratons.shp` (and its companion files): the 114 craton polygons of Hasterok et al. (2022).

## Notes

- Moho depth is the surface elevation minus the top-of-mantle boundary (the last CRUST1.0 column); pressure at the Moho is the sum of layer density x g x thickness. Statistics are area-weighted by the cosine of latitude over the 5,709 cells that fall inside the craton polygons (8.0% of the sphere).
- Inputs marked as assumptions in the scripts (densities, cubit length, ocean volume) are placeholders recorded in Section 14 of the model.
- The images regenerate to identical PNGs. The SVGs differ only in embedded metadata (date and identifiers).

`data/` (git-ignored, not in the repository) holds the Lucazeau (2019) supplement files for `craton_heatflow.py`: `2019gc008389-sup-0003-data_set_si-s01.zip` (grid `HFgrid14.csv`), `-sup-0004-...-s02.zip` (point data `NGHF.csv`) and the figure PDF. Source: doi:10.1029/2019GC008389, Supporting Information.
