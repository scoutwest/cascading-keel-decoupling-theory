# Scripts

Scripts that produce the images and the numbers quoted in `Cascading Moho Collapse Theory (CMCT).md`. They are the working scripts from draft 6, not the Phase C rebuild; Phase C will replace them with a single tested package (see Section 11 of the model).

| Script | Produces | Inputs |
|--------|----------|--------|
| `section_arithmetic.py` | The tables in Sections 14-18 (McKenzie subsidence, water sea-level equivalents, layer thickness, the Section 16 budgets, the Arrhenius scaling, pressure at depth) | none |
| `craton_moho_stats.py` | The Section 19 statistics (craton Moho depth and pressure, lower-crust thickness) | CRUST1.0 and the Hasterok craton polygons |
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

## Data the two data scripts need

The data are not in this repository. The scripts look for them in the folder that holds this repository's folder (override with the `FLOOD_MODELS_ROOT` environment variable):

- `Noah's Pangea v1/data/reference/crust1/crust1.bnds` and `crust1.rho` (CRUST1.0, Laske et al., 1°x1° grid; nine layers per cell).
- `Noah's Pangea v2/data/hasterok_2022/cratons.shp` (and its companion files): the 114 craton polygons of Hasterok et al. (2022).

## Notes

- Moho depth is the surface elevation minus the top-of-mantle boundary (the last CRUST1.0 column); pressure at the Moho is the sum of layer density x g x thickness. Statistics are area-weighted by the cosine of latitude over the 5,709 cells that fall inside the craton polygons (8.0% of the sphere).
- Inputs marked as assumptions in the scripts (densities, cubit length, ocean volume) are placeholders recorded in Section 14 of the model.
- The images regenerate to identical PNGs. The SVGs differ only in embedded metadata (date and identifiers).
