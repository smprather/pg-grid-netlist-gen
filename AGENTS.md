# AGENTS.md

## Requirement Keywords

- MUST: Required behavior.
- MUST NOT: Strictly forbidden behavior.
- SHOULD: Strong recommendation; deviations must be intentional.
- MAY: Optional behavior.

## Project Overview

This project is a script that generates:

1. The power and ground grid for an ASIC place-and-route block.
1. A 2D visualization of the grid using Plotly and a browser.
1. A SPICE netlist containing:
   1. Resistors and capacitors for the PG grid.
   1. Standard-cell instances.
   1. HSPICE `.MEAS` statements that measure per-instance IR drop, plus the maximum and average of all per-instance IR drops.
1. An ASCII report with detailed generation and extraction metrics.

## Definitions

- Standard cell: A logic-cell abstraction defined by cell shape and pin locations.
  In this project, cells are vertically flipped on alternating rows.
- Instance: A placed instantiation of a standard cell.
- ITF file: Interconnect Technology Format file describing conductor and via layers.
- PLOC: Pad locations, i.e., `(x, y)` ideal supply-connection points.
- PG: Power and ground.
- Track: Router placement/routing grid unit. Cell height is commonly specified in tracks.
- Grid layer: A routing layer containing full-span horizontal or vertical stripes.
- Staple layer: A layer represented by local square metal patches with vias above and below.
- Flight line: A visualization line from one instance output pin to another instance input pin.
- IR drop: `V(power_pin) - V(ground_pin)` averaged over a configured transient window.
- CPP: Contacted Poly Pitch (minimum horizontal spacing between FinFET transistors).
- Site width: Horizontal placement quantum. A cell left edge MUST lie on an integer multiple of site width.
- Unateness: Determines whether a logic cell is inverting from input to output.
  This is needed to know how to construct the measurement of the averaging windows
  for IR drop measurement. Positive means non-inverting. Negative means inverting.
- Dcap cell: Decoupling capacitance cell. These cells have no signal input or output. They are used to reduce
  IR drop locally.

## YAML Input File

The YAML file (currently `grid_specs.yaml`) defines technology, grid, placement, PLOC, and SPICE settings.

### Global Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `random_seed_value` | `number` or `"random"` | Yes | N/A | Global randomness control. |
| `itf.file` | `string` | Yes | N/A | Path to ITF file. |
| `itf.units` | `string` | Yes | N/A | ITF distance unit label. |
| `units.distance` | `string` | Yes | N/A | Distance unit label (e.g. `um`). |
| `units.capacitance` | `string` | Yes | N/A | Capacitance unit label (e.g. `ff`). |
| `units.resistance` | `string` | Yes | N/A | Resistance unit label (e.g. `ohms`). |
| `units.time` | `string` | Yes | N/A | Time unit label (e.g. `ps`). |
| `beol_thickness` | `number` | Yes | distance | BEOL stack thickness used for reporting/visualization context. |

### `standard_cells[]` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `standard_cells[].name` | `string` | Yes | N/A | Cell name/subckt name basis. |
| `standard_cells[].height` | `number` | Yes | distance | Physical cell height. |
| `standard_cells[].width` | `number` | Yes | distance | Physical cell width. |
| `standard_cells[].pins[]` | `list` | Yes | N/A | Pin definitions. |
| `standard_cells[].pins[].name` | `string` | Yes | N/A | Pin name. |
| `standard_cells[].pins[].type` | `string` | Yes | N/A | Allowed: `signal`, `power`, `ground`. |
| `standard_cells[].pins[].direction` | `string` | Conditionally | N/A | Required when `type=signal`; allowed: `input`, `output`, `inout`. |
| `standard_cells[].pins[].location` | `string` | Yes | N/A | Allowed: `left`, `right`, `top`, `bottom`. |
| `standard_cells[].spice_port_order` | `string` | Yes | N/A | Space-delimited ordered port list. |
| `standard_cells[].unateness` | `string` | No | N/A | Allowed: `positive` (non-inverting), `negative` (inverting). Default: `positive`. |
| `standard_cells[].spice_netlist_file` | `string` | No | path | Optional external SPICE file to include with `.include`; file existence is not validated at generation time. |

### `ploc` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `ploc.pitch` | `number` | Yes | distance | PLOC pitch. |
| `ploc.staggered` | `bool` | Yes | N/A | Enables staggered rows/columns of PLOCs. |
| `ploc.offset_from_origin.x` | `number` | Yes | distance | X offset from origin. |
| `ploc.offset_from_origin.y` | `number` | Yes | distance | Y offset from origin. |
| `ploc.visualizer_render_diameter` | `number` | No | distance | Diameter of PLOC circles in 2D visualization. |

### `spice_netlist` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `spice_netlist.cell_chains.chain_cell` | `string` | Yes | N/A | `standard_cells` name for chains. |
| `spice_netlist.cell_chains.cell_output_loads.in_chain.resistance` | `number` | Yes | resistance | Per-link R. |
| `spice_netlist.cell_chains.cell_output_loads.in_chain.capacitance` | `number` | Yes | capacitance | Per-link C. |
| `spice_netlist.cell_chains.cell_output_loads.in_chain.number_pi_segments` | `int` | Yes | N/A | Pi sections. |
| `spice_netlist.cell_chains.cell_output_loads.end_of_chain.*` | same | Yes | same | Last-stage load (also pi model). |
| `spice_netlist.cell_chains.chain_input_stimulus.period` | `number` | Yes | time | PULSE period. |
| `spice_netlist.cell_chains.chain_input_stimulus.transition_time` | `number` | Yes | time | Rise/fall time. |
| `spice_netlist.cell_chains.max_instance_count_per_chain` | `int` | Yes | N/A | Max instances per chain. |
| `spice_netlist.transient_simulation.total_time` | `number` | Yes | time | Transient stop time. |
| `spice_netlist.transient_simulation.time_step` | `number` | Yes | time | Transient step. |
| `spice_netlist.ir_drop_measurement.averaging_window.start` | `number` | Yes | % | Start threshold as % of power voltage on the input-net transition. |
| `spice_netlist.ir_drop_measurement.averaging_window.end` | `number` | Yes | % | End threshold as % of power voltage on the output-net transition. |

### `visualizer` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `visualizer.initial_visible_objects` | `list[string]` | No | N/A | Legend entry names to show on load; all others start as `legendonly`. |

### `standard_cell_placement` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `standard_cell_placement.row_height` | `number` | Yes | distance | Placement row pitch. |
| `standard_cell_placement.site_width` | `number` | Yes | distance | Placement site pitch. |
| `standard_cell_placement.min_space.x` | `number` | Yes | distance | Minimum center-to-center spacing in X. |
| `standard_cell_placement.min_space.y` | `number` | Yes | distance | Minimum center-to-center spacing in Y. |
| `standard_cell_placement.stagger_row_start.range` | `number` | Yes | distance | Row-start stagger amplitude. |
| `standard_cell_placement.stagger_row_start.random` | `bool` | Yes | N/A | Randomize row-start offset if true. |
| `standard_cell_placement.dcap_cells.enabled` | `bool` | No | N/A | Enable dcap insertion. Default: false. |
| `standard_cell_placement.dcap_cells.cell` | `string` | No | N/A | `standard_cells` name for dcaps. |
| `standard_cell_placement.dcap_cells.max_density_pct` | `number` | No | % | Max area density (0-100). |

### `pg_nets` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `pg_nets.power.name` | `string` | Yes | N/A | Power-net name; MUST match cell power pin name. |
| `pg_nets.power.voltage` | `number` | Yes | volts | Ideal supply voltage at power PLOCs. |
| `pg_nets.ground.name` | `string` | Yes | N/A | Ground-net name; MUST match cell ground pin name. |

### `grid` Schema

| Key path | Type | Required | Units | Notes |
|---|---|---|---|---|
| `grid.size.rows` | `int` | Yes | rows | Number of placement rows. |
| `grid.size.sites` | `int` | Yes | sites | Number of placement sites per row. |
| `grid.layer_usage.<LAYER>.type` | `string` | Yes | N/A | Enum: `g` (grid stripe), `s` (staple). |
| `grid.layer_usage.<LAYER>.width` | `number` | Yes | distance | Metal width for configured layer. |
| `grid.layer_usage.<LAYER>.pitch` | `number` | Yes | distance | Stripe pitch for `g`; retained for uniformity on `s`. |

#### Layer-usage enum and implicit defaults

- The layer width and pitch are in terms of WMIN, and SMIN.
  - For example, if width=1.0, then the actual width is 1.0 * WMIN
  - For example, if pitch=1.0, then the actual pitch is 1.0 * (WMIN+SMIN)
- `grid.layer_usage.<LAYER>.type` MUST be one of:
  - `g`: full-span stripe layer.
  - `s`: staple layer.
- The lowest metal layer in ITF is implicit and SHOULD be omitted from `grid.layer_usage`.
- For the implicit lowest layer:
  - `type` MUST be treated as `g`.
  - `width` MUST be the ITF `WMIN` of that lowest conductor.
  - `pitch` MUST equal `standard_cell_placement.row_height`.
- All explicitly configured layer names MUST exist as ITF `CONDUCTOR` names.

## Netlist and Extraction Rules

- The grid MUST be constructed bottom-up.
- A via MUST be instantiated anywhere the metal above and below intersect.
- Via shape MUST be square and satisfy ITF `AREA`.
- PLOC generation MUST start with the power net and alternate horizontally between power and ground.
- If a PLOC does not land on a top-layer stripe of the intended PG net, it MUST snap to the nearest valid stripe centerline.
- Anywhere an instance power or ground pin lands on the lowest layer of metal, that segement of metal must be broken into two
  segments where the pin lands. This is not required if the cell pin lands exactly on the edge of a segment of metal. In this
  case, just tap the cell into the grid at the pre-existing node.
- Add HSPICE .PROBE to probe the voltage of every pin of every instance. `.PROBE V(X*)`

### RC extraction formulas

- Metal segment resistance MUST use `R = RPSQ * (L / W)`.
- Via resistance MUST use ITF `RPV`.
- Segment capacitance MUST include only:
  - Plate-to-substrate capacitance.
  - Fringe-to-substrate capacitance (isolated-wire assumption).
- Capacitance to other conductors MUST NOT be included.
- Via capacitance and staple-shape capacitance MUST NOT be included.

### Staple resistance model

For a staple layer between adjacent routed layers (example `M7-V6-M6-V5-M5`), vertical connection resistance MUST include three series elements:

1. `R(V6)` from ITF `RPV`.
1. `R(M6)` as one square: `RPSQ * 1.0`.
1. `R(V5)` from ITF `RPV`.

## Placement and Connectivity Rules

- The bottom ITF metal layer MUST be treated as horizontal.
- Higher layers MUST alternate horizontal/vertical by level.
- Standard-cell signal pins MUST be placed at the midpoint of their declared side.
- Cell rows MUST alternate vertical flip orientation.
- Cell top/bottom edges MUST lie on PG metal centerlines so power pins land on power metal and ground pins on ground metal.
- Cell placement MUST honor row and site grids.
- Cells SHOULD be distributed according to `standard_cell_placement.min_space`.
- PG net names need not match standard-cell PG pin names.
  - Connectivity is established by matching a power-type PG net to a power-type pin on the cell.
    Same for ground type net to pin.
- MUST NOT allow placement overlap with any other cells, dcap or chain.

### Decoupling capacitance cell insertion

- If enabled, insert decoupling caps (dcaps) until the maximum density is reached. Density is defined as the total area of the
  dcap cells divided by the total area of the grid.
- Draw them in the visualizer using a different legend group than the chain cell instances so that visibility can be independently
  controlled.

### Chain generation

- Start a chain from a random instance and connect output-to-input across randomly chosend instances.
- A cell can only be used in a chain once
- Continue until `max_instance_count_per_chain` is reached per chain.
- Continue creating chains until all instances are assigned.
- Chain-loading values MUST come from `spice_netlist.cell_chains.cell_output_loads.in_chain`, with last-stage values from `end_of_chain`.
- Input stimulus MUST use a SPICE `PULSE` source with configured period and transition time.

## Determinism and Randomness

- If `random_seed_value` is numeric, all randomized operations MUST be deterministic and reproducible for that seed.
- If `random_seed_value` is `"random"`, a fresh seed MAY be generated per run.
- All random subsystems MUST draw from a single seeded RNG stream (placement staggering, chain construction, and any randomized PLOC behaviors).
- The ASCII report MUST include the effective seed used for the run.

## Measurement Timing Details

- Per-instance IR-drop averaging windows MUST be measured using transition crossing events, not fixed absolute times.
- For each instance, two measurement windows MUST be generated:
  - output-rise window (`..._RISE`)
  - output-fall window (`..._FALL`)
- For `unateness=positive` (non-inverting):
  - output rise uses input rise for start trigger
  - output fall uses input fall for start trigger
- For `unateness=negative` (inverting):
  - output rise uses input fall for start trigger
  - output fall uses input rise for start trigger
- Start threshold level MUST use `averaging_window.start` as a percentage of `pg_nets.power.voltage`.
  - For falling-edge start triggers, use the complementary level (`VDD * (1 - start%)`).
- End threshold level MUST use `averaging_window.end` as a percentage of `pg_nets.power.voltage`.
  - For falling-edge end triggers, use the complementary level (`VDD * (1 - end%)`).

## Validation Rules

The generator MUST fail fast with a clear error message when validation fails.

- YAML keys and value types MUST match the schema above.
- `grid.layer_usage` keys MUST be valid ITF conductor names.
- Explicitly configured layer stack MUST be vertically connectable using ITF vias.
- `grid.layer_usage` MUST NOT include the implicit lowest layer.
- `width > 0`, `pitch > 0`, `row_height > 0`, `site_width > 0`.
- `grid.size.rows >= 1` and `grid.size.sites >= 1`.
- `max_instance_count_per_chain >= 1`.
- Averaging-window percentages MUST satisfy `0 <= start < end <= 100`.
- `spice_port_order` pins MUST exactly match declared pin names for each standard cell.

## Visualization Requirements

- Use Plotly legend entries for layer visibility control. A dropdown MUST NOT be used.
- Initial visibility MUST show:
  - Standard-cell instances.
  - Bottom 4 metal layers of the generated grid.
  - Via layers that touch those bottom 4 metal layers.
- All other layers SHOULD initialize as `visible='legendonly'`.
- Flight lines MUST be drawn in their own legend trace group so users can toggle them.
- A cross-section visualization of the ITF stack MUST be placed below the 2D render.

## Output Contract

Unless overridden by explicit CLI options, outputs MUST be written to `output/` with these default filenames:

1. `output/pg_grid_visualization.html`
   - Contains the 2D grid/cell render and ITF cross-section visualization.
1. `output/pg_grid_netlist.sp`
   - Contains PG RC network, instances, sources, stimulus, and `.MEAS` statements.
1. `output/pg_grid_summary.txt`
   - ASCII report with at least the following fields:
     - Effective random seed.
     - Total PG capacitance.
     - Total resistor count.
     - Total capacitor count.
     - Total chain instance count.
     - Total chain count.
     - Total dcap cell count.
     - Total dcap cell density.

## Generation Order

Implementation SHOULD follow this sequence:

1. Parse YAML and ITF.
1. Validate schema and technology/stack consistency.
1. Build the PG grid from bottom to top (including implicit lowest layer behavior).
1. Insert vias and staple structures.
1. Place and orient standard-cell instances.
1. Generate/snap PLOCs.
1. Build instance chains and signal loading.
1. Extract RC network and generate SPICE.
1. Add `.MEAS` statements for per-instance, max, and average IR drop.
1. Generate visualization and ASCII summary outputs.

## Documentation Generation Rules

1. Create and maintain the top-level README.md.

- MUST contain a link to the live-demo at https://smprather.github.io/pg-grid-netlist-gen
- MUST contain a link to the configuration file documentation.

2. Create the configuration file documentation in the docs/ directory.

## Tech Stack

- Python 3.14
- Astral uv for project management
- Rich-Click for CLI
- Plotly for visualization
- Use python's Pathlib whenever possible
