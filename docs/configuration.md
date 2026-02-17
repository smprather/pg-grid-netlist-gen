# Configuration Reference

The generator is configured by `grid_specs.yaml`.

## Top-Level Keys

- `random_seed_value`
- `itf`
- `units`
- `beol_thickness`
- `standard_cells`
- `ploc`
- `visualizer` (optional)
- `spice_netlist`
- `standard_cell_placement`
- `pg_nets`
- `grid`

## Key Notes

- `grid.layer_usage.<LAYER>.width` and `pitch` are multipliers:
  - width: `WMIN * width`
  - pitch: `(WMIN + SMIN) * pitch`
- The lowest metal layer is implicit and must not be listed in `grid.layer_usage`.
- Standard-cell PG pins are matched by pin `type`:
  - `power` pins connect to `pg_nets.power.name`
  - `ground` pins connect to `pg_nets.ground.name`
- Optional external SPICE model include:
  - `standard_cells[].spice_netlist_file`
- Optional visualization defaults:
  - `visualizer.initial_visible_objects`
  - Example: `["Cells", "M0", "V0", "M1"]`

## Run

```bash
./pg_grid_netlist_gen grid_specs.yaml
```
