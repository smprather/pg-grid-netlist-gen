# Configuration Reference

The generator is configured by `config.yaml`. This document describes all
configuration parameters.

## Global Settings

| Key | Type | Description |
|-----|------|-------------|
| `random_seed_value` | `number` or `"random"` | Controls determinism. A numeric value produces repeatable results; `"random"` generates a fresh seed per run. |
| `itf.file` | `string` | Path to the ITF (Interconnect Technology Format) file. Resolved relative to the config file. |
| `itf.units` | `string` | Distance unit used in the ITF file (e.g. `um`, `nm`). |
| `units.distance` | `string` | Distance unit for all config values (e.g. `um`). |
| `units.capacitance` | `string` | Capacitance unit (e.g. `ff`, `pf`). |
| `units.resistance` | `string` | Resistance unit (e.g. `ohms`). |
| `units.time` | `string` | Time unit (e.g. `ps`). |
| `feol_thickness` | `number` | Distance from the bottom of the lowest metal layer to the substrate, used for capacitance-to-substrate calculation. |

## Standard Cells (`standard_cells[]`)

| Key | Type | Description |
|-----|------|-------------|
| `name` | `string` | Cell name, used as the SPICE subcircuit name. |
| `height` | `number` | Physical cell height in distance units. |
| `width` | `number` | Physical cell width in distance units. |
| `pins[].name` | `string` | Pin name. |
| `pins[].type` | `string` | Pin type: `signal`, `power`, or `ground`. |
| `pins[].direction` | `string` | Required for signal pins: `input`, `output`, or `inout`. |
| `pins[].location` | `string` | Pin side: `left`, `right`, `top`, or `bottom`. Signal pins are placed at the midpoint of their declared side. |
| `unateness` | `string` | `positive` (non-inverting) or `negative` (inverting). Default: `positive`. Affects IR-drop measurement window edge mapping. |
| `spice_netlist_file` | `string` | Path to a SPICE file containing the `.subckt` definition. Port order is parsed from the `.subckt` line; ports must match declared pins exactly. |

## PLOC Settings (`ploc`)

| Key | Type | Description |
|-----|------|-------------|
| `pitch` | `number` | Spacing between PLOC points in distance units. |
| `staggered` | `bool` | If true, alternating rows of PLOCs are offset by half a pitch. |
| `offset_from_origin.x` | `number` | X offset from the grid origin for the first PLOC. |
| `offset_from_origin.y` | `number` | Y offset from the grid origin for the first PLOC. |
| `visualizer_render_diameter` | `number` | Diameter of PLOC circles in the 2D visualization. Default: `1.0`. |

## Visualizer (`visualizer`)

| Key | Type | Description |
|-----|------|-------------|
| `initial_visible_objects` | `list[string]` | Legend entry names to show on initial load. All other entries start as `legendonly`. Example: `["Cells", "M0", "V0", "M1"]`. |

## SPICE Netlist (`spice_netlist`)

### R/C Scaling

| Key | Type | Description |
|-----|------|-------------|
| `scaling.resistance` | `number` | Multiplier applied to all grid extraction resistance values. Default: `1.0`. |
| `scaling.capacitance` | `number` | Multiplier applied to all grid extraction capacitance values. Default: `1.0`. |

### Cell Chains (`cell_chains`)

| Key | Type | Description |
|-----|------|-------------|
| `cell` | `string` | Name of the standard cell used for chains (must exist in `standard_cells[]` and have signal pins). |
| `interconnect.resistance` | `number` | Total resistance for the pi-model interconnect between consecutive cells. |
| `interconnect.capacitance` | `number` | Total capacitance for the pi-model interconnect between consecutive cells. |
| `interconnect.number_pi_segments` | `int` | Number of pi sections in the interconnect. Must be >= 1. |
| `end_of_chain_load.resistance` | `number` | Resistance of the RC load on the last cell's output. |
| `end_of_chain_load.capacitance` | `number` | Capacitance of the RC load on the last cell's output. |
| `chain_input_stimulus.period` | `number` | PULSE source period in time units. |
| `chain_input_stimulus.transition_time` | `object` | Gaussian variation for PULSE rise/fall time. See below. |
| `chain_input_stimulus.initial_delay` | `object` | Gaussian variation for PULSE initial delay. See below. |
| `max_instance_count_per_chain` | `int` | Maximum number of instances per chain. Must be >= 1. |

#### Gaussian Variation Object

Used by `transition_time` and `initial_delay`. Each chain's PULSE source samples independently from the seeded RNG.

| Key | Type | Description |
|-----|------|-------------|
| `nominal` | `number` | Mean value. |
| `sigma` | `number` | Standard deviation. Must be >= 0. |
| `floor` | `number` | Minimum clamp. Must be <= `nominal`. |
| `ceiling` | `number` | Maximum clamp. Must be >= `nominal`. |

### IR Drop Measurement

| Key | Type | Description |
|-----|------|-------------|
| `ir_drop_measurement.averaging_window.start` | `number` | Start threshold as a percentage (0-100) of VDD on the input transition. |
| `ir_drop_measurement.averaging_window.end` | `number` | End threshold as a percentage (0-100) of VDD on the output transition. Must satisfy `start < end`. |

### User Defined Lines

| Key | Type | Description |
|-----|------|-------------|
| `user_defined_lines` | `list[string]` | Raw SPICE lines written verbatim after the netlist header. Use for `.include`, `.option`, `.param`, etc. Default: empty. |

## Standard Cell Placement (`standard_cell_placement`)

| Key | Type | Description |
|-----|------|-------------|
| `row_height` | `number` | Placement row pitch in distance units. |
| `site_width` | `number` | Placement site pitch (horizontal quantum). Cell left edges must be on integer multiples. |
| `min_space.x` | `number` | Minimum center-to-center spacing in X between chain cells. |
| `min_space.y` | `number` | Minimum center-to-center spacing in Y between chain cells. |
| `stagger_row_start.range` | `number` | Row-start stagger amplitude. Set to `0.0` to disable staggering. |
| `stagger_row_start.random` | `bool` | If true, randomize the row-start offset; if false, alternate between `0.0` and the range value. |
| `dcap_cells.enabled` | `bool` | Enable decoupling capacitance cell insertion. Default: `false`. |
| `dcap_cells.cell` | `string` | Name of the standard cell to use for dcaps (must exist in `standard_cells[]`). |
| `dcap_cells.max_density_pct` | `number` | Maximum dcap area density as a percentage (0-100) of total grid area. |

## PG Nets (`pg_nets`)

| Key | Type | Description |
|-----|------|-------------|
| `power.name` | `string` | Power net name. Connectivity is by pin type, not pin name. |
| `power.voltage` | `number` | Ideal supply voltage applied at power PLOCs. |
| `ground.name` | `string` | Ground net name. Connectivity is by pin type, not pin name. |

## Grid (`grid`)

| Key | Type | Description |
|-----|------|-------------|
| `via_min_space_factor` | `number` | Via min-space = via side length x factor. Used for multi-via packing. Default: `1.5`. |
| `size.rows` | `int` | Number of placement rows. Must be >= 1. |
| `size.sites` | `int` | Number of placement sites per row. Must be >= 1. |

### Layer Usage (`grid.layer_usage.<LAYER>`)

Each entry configures a metal layer by its ITF conductor name.

| Key | Type | Description |
|-----|------|-------------|
| `type` | `string` | `grid` for full-span grid stripe layer, `staple` for staple (local patch) layer. |
| `width` | `number` | Width multiplier: `actual_width = width * WMIN` from ITF. |
| `pitch` | `number` | Pitch multiplier: `actual_pitch = pitch * (WMIN + SMIN)` from ITF. |

**Important notes:**

- The lowest metal layer in the ITF is implicit and must **not** be listed in `layer_usage`.
- For the implicit lowest layer: type is always `grid`, width is `2.0 * WMIN`, pitch equals `row_height`.
- All layer names must exist as ITF `CONDUCTOR` names.
- For grid-to-grid via crossings, the maximum number of minimum-spaced vias are packed
  in the overlap rectangle for visualization. For extraction, the group is modeled as a
  single via with `R = RPV / via_count`.

## Run

```bash
./pg_grid_netlist_gen config.yaml
```
