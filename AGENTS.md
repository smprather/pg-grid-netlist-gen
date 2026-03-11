# AGENTS.md

## Requirement Keywords

- MUST: Required behavior.
- MUST NOT: Strictly forbidden behavior.
- SHOULD: Strong recommendation; deviations must be intentional.
- MAY: Optional behavior.

## Formatting Rules

- Markdown tables in this file MUST NOT exceed 110 columns per line.
  - Use relative key paths (omit the section-header prefix).
  - Abbreviate: "Req" for Required; "Y" / "N" / "Cond" for values.
  - Use line-wrapping within a cell if needed to stay below the limit.
  - Keep Notes column concise.

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
- IR drop: `avg(V(power_pin) - V(ground_pin)) / VDD_nominal` — percent drop averaged over a configured transient window.
- CPP: Contacted Poly Pitch (minimum horizontal spacing between FinFET transistors).
- Site width: Horizontal placement quantum. A cell left edge MUST lie on an integer multiple of site width.
- Unateness: Determines whether a logic cell is inverting from input to output.
  This is needed to know how to construct the measurement of the averaging windows
  for IR drop measurement. Positive means non-inverting. Negative means inverting.
- Dcap cell: Decoupling capacitance cell. These cells have no signal input or output. They are used to reduce
  IR drop locally.

## YAML Input File

The YAML file (currently `config.yaml`) defines technology, grid, placement, PLOC, and SPICE settings.

### Global Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `random_seed_value` | `number` or `"random"` | Y | N/A | Global randomness control. |
| `itf.file` | `string` | Y | N/A | Path to ITF file. |
| `itf.units` | `string` | Y | N/A | ITF distance unit label. |
| `units.distance` | `string` | Y | N/A | Distance unit (e.g. `um`). |
| `units.capacitance` | `string` | Y | N/A | Capacitance unit (e.g. `ff`). |
| `units.resistance` | `string` | Y | N/A | Resistance unit (e.g. `ohms`). |
| `units.time` | `string` | Y | N/A | Time unit (e.g. `ps`). |
| `feol_thickness` | `number` | Y | distance | M0-to-substrate distance for cap calc. |

### `standard_cells[]` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `name` | `string` | Y | N/A | Cell name / subckt name. |
| `height` | `number` | Y | distance | Physical cell height. |
| `width` | `number` | Y | distance | Physical cell width. |
| `pins[]` | `list` | Y | N/A | Pin definitions. |
| `pins[].name` | `string` | Y | N/A | Pin name. |
| `pins[].type` | `string` | Y | N/A | `signal`, `power`, or `ground`. |
| `pins[].direction` | `string` | Cond | N/A | For `type=signal`. `input`/`output`/`inout`. |
| `pins[].location` | `string` | Y | N/A | `left`, `right`, `top`, or `bottom`. |
| `unateness` | `string` | N | N/A | `positive` or `negative`. Default: `positive`. |
| `spice_netlist_file` | `string` | Y | path | SPICE file with `.subckt`; port order is parsed from it. |

### `ploc` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `pitch` | `number` | Y | distance | PLOC pitch. |
| `staggered` | `bool` | Y | N/A | Staggered rows/columns of PLOCs. |
| `offset_from_origin.x` | `number` | Y | distance | X offset from origin. |
| `offset_from_origin.y` | `number` | Y | distance | Y offset from origin. |
| `visualizer_render_diameter` | `number` | N | distance | PLOC circle diameter in visualization. |

### `spice_netlist` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `scaling.resistance` | `number` | N | N/A | Grid R multiplier. Default: 1.0. |
| `scaling.capacitance` | `number` | N | N/A | Grid C multiplier. Default: 1.0. |
| `cell_chains.cell` | `string` | Y | N/A | `standard_cells` name for chains. |
| `cell_chains.interconnect.resistance` | `number` | Y | resistance | Total R for pi interconnect. |
| `cell_chains.interconnect.capacitance` | `number` | Y | capacitance | Total C for pi interconnect. |
| `cell_chains.interconnect.number_pi_segments` | `int` | Y | N/A | Pi sections between cells. |
| `cell_chains.end_of_chain_load.resistance` | `number` | Y | resistance | Last-stage RC load R. |
| `cell_chains.end_of_chain_load.capacitance` | `number` | Y | capacitance | Last-stage RC load C. |
| `cell_chains.chain_input_stimulus.period` | `number` | Y | time | PULSE period. |
| `cell_chains.chain_input_stimulus.transition_time` | `object` | Y | time | Gaussian variation: see below. |
| `cell_chains.chain_input_stimulus.initial_delay` | `object` | Y | time | Gaussian variation: see below. |
| `cell_chains.max_instance_count_per_chain` | `int` | Y | N/A | Max instances per chain. |
| `ir_drop_measurement.averaging_window.start` | `number` | Y | % | Start % on input transition. |
| `ir_drop_measurement.averaging_window.end` | `number` | Y | % | End % on output transition. |
| `user_defined_lines` | `list[string]` | N | N/A | Raw SPICE lines emitted after header. |

#### Gaussian variation object schema

Used by `transition_time` and `initial_delay`. Each chain samples independently.

| Key | Type | Req | Notes |
|-----|------|-----|-------|
| `nominal` | `number` | Y | Mean value. |
| `sigma` | `number` | Y | Standard deviation (>= 0). |
| `floor` | `number` | Y | Minimum clamp (<= nominal). |
| `ceiling` | `number` | Y | Maximum clamp (>= nominal). |

### `visualizer` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `initial_visible_objects` | `list[string]` | N | N/A | Legend entries visible on load; rest hidden. |

### `standard_cell_placement` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `row_height` | `number` | Y | distance | Placement row pitch. |
| `site_width` | `number` | Y | distance | Placement site pitch. |
| `min_space.x` | `number` | Y | distance | Min center-to-center X spacing. |
| `min_space.y` | `number` | Y | distance | Min center-to-center Y spacing. |
| `stagger_row_start.range` | `number` | Y | distance | Row-start stagger amplitude. |
| `stagger_row_start.random` | `bool` | Y | N/A | Randomize row-start offset. |
| `dcap_cells.enabled` | `bool` | N | N/A | Enable dcap insertion. Default: false. |
| `dcap_cells.cell` | `string` | N | N/A | `standard_cells` name for dcaps. |
| `dcap_cells.max_density_pct` | `number` | N | % | Max area density (0-100). |

### `pg_nets` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `power.name` | `string` | Y | N/A | Power-net name; MUST match cell power pin. |
| `power.voltage` | `number` | Y | volts | Ideal supply voltage at power PLOCs. |
| `ground.name` | `string` | Y | N/A | Ground-net name; MUST match cell ground pin. |

### `grid` Schema

| Key | Type | Req | Units | Notes |
|-----|------|-----|-------|-------|
| `via_min_space_factor` | `number` | N | N/A | Via min-space = via_side * factor. Default: 1.5. |
| `size.rows` | `int` | Y | rows | Number of placement rows. |
| `size.sites` | `int` | Y | sites | Number of sites per row. |
| `layer_usage.<LAYER>.type` | `string` | Y | N/A | `grid` or `staple`. |
| `layer_usage.<LAYER>.width` | `number` | Y | distance | Metal width for layer. |
| `layer_usage.<LAYER>.pitch` | `number` | Y | distance | Stripe pitch for `grid`. |

#### Layer-usage enum and implicit defaults

- The layer width and pitch are in terms of WMIN and SMIN.
  - For example, if width=1.0, then the actual width is 1.0 * WMIN
  - For example, if pitch=1.0, then the actual pitch is 1.0 * (WMIN+SMIN)
- `grid.layer_usage.<LAYER>.type` MUST be one of:
  - `grid`: A full-span metal stripe layer.
  - `staple`: A staple layer. Just enough metal to get to the above and below VIA layer.
- The lowest metal layer in ITF is implicit and SHOULD be omitted from `grid.layer_usage`.
- For the implicit lowest layer:
  - `type` MUST be treated as `grid`.
  - `width` MUST be `2.0 * WMIN` of that lowest conductor.
  - `pitch` MUST equal `standard_cell_placement.row_height`.
- All explicitly configured layer names MUST exist as ITF `CONDUCTOR` names.

## Netlist and Extraction Rules

- The grid MUST be constructed bottom-up.
- A via MUST be instantiated anywhere the metal above and below intersect.
- Via shape MUST be square and satisfy ITF `AREA`.
- PLOC generation MUST start with the power net and alternate horizontally between power and ground.
- If a PLOC does not land on a top-layer stripe of the intended PG net, it MUST snap to the nearest valid stripe centerline.
- Anywhere an instance power or ground pin lands on the lowest layer of metal, that segment of metal must be broken into two
  segments where the pin lands. This is not required if the cell pin lands exactly on the edge of a segment of metal. In this
  case, just tap the cell into the grid at the pre-existing node.
- Add HSPICE .PROBE to probe the voltage of every pin of every instance. `.PROBE V(X*)`
- Netlist node names should begin with the net name of the segment the R or C is associated with.
- New netlist nodes MUST be created at:
  - The top and bottom of VIAs.
  - Where cell or dcap cell PG pins touch the lowest layer of metal.
  - Where PLOC points touch the layer of metal they connect to.
  - Where any metal stripe touches the edge of the defined grid area.
- For the capacitance of each metal segment, connect half of the cap at each end of the segment (pi model).
- Connect all capacitors to the ground net. Do not attempt to connect to the nearest ground net segment. For example,
  if VSS is the ground net name, then all caps should connect to VSS.

### RC extraction formulas

- Only do capacitance calculation on PG nets of power type.
- For capacitance calculation to substrate, the distance from the bottom of the lowest layer of metal to the substrate is feol_thickness.
- Metal segment resistance MUST use `R = RPSQ * (L / W)`.
- Via resistance MUST use ITF `RPV`.
- Segment capacitance MUST include:
  - Plate-to-substrate capacitance.
  - Fringe-to-substrate capacitance (isolated-wire assumption).
  - Plate-to-plate capacitance from each power net segment to the nearest ground net on same layer.
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
- For each VIA layer, calculate the min-space as the length of a side of a via, times grid.via_min_space_factor.
- For the visualizer, place the maximum number of minimum-spaced VIAs between adjacent-layer metal stripes of
  the same net. However, for the purposes of node creation and resistance extraction, treat the group of vias
  as a single via with R=RPV/number_of_vias_used_in_visualizer. In other words, for extraction, model the group
  of vias as a single r-scaled via located in the middle of the metal stripe.

### Decoupling capacitance cell insertion

- If enabled, insert decoupling caps (dcaps) until the maximum density is reached. Density is defined as the total area of the
  dcap cells divided by the total area of the grid.
- Draw them in the visualizer using a different legend group than the chain cell instances so that visibility can be independently controlled.

### Chain generation

- Start a chain from a random instance and connect output-to-input across randomly chosen instances.
- A cell can only be used in a chain once.
- Continue until `max_instance_count_per_chain` is reached per chain.
- Continue creating chains until all instances are assigned.
- Between each pair of consecutive cells, a multi-segment pi-model interconnect MUST bridge the output of cell N to the input of cell N+1, using `cell_chains.interconnect` (R, C, number_pi_segments).
- The last cell in each chain MUST have a simple RC load on its output from `cell_chains.end_of_chain_load`.
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
- `.subckt` ports in `spice_netlist_file` MUST exactly match declared pin names for each standard cell.

## Visualization Requirements

- Use Plotly legend entries for layer visibility control. A dropdown MUST NOT be used.

- Initial visibility MUST show:

  - Standard-cell instances.
  - Bottom 4 metal layers of the generated grid.
  - Via layers that touch those bottom 4 metal layers.

- All other layers SHOULD initialize as `visible='legendonly'`.

- Flight lines MUST be drawn in their own legend trace group so users can toggle them.

- A cross-section visualization of the ITF stack MUST be placed below the 2D render.

  - Use hover-data to label the layers.
  - Any layers with the same thickness should share the same color.
  - Render FEOL as a layer. The bottom of the FEOL layer is the y=0 point.
  - Render the substrate as a layer 3x the thickness of the FEOL layer.
  - Do not include the "\<via_layer>\_diel" layers. Just render the VIAs.
  - For substrate only, add a label to the layer.
  - Don't use a legend.

- Use include_plotlyjs=True to support offline usage, but only in the first div (see below).

- Use this method to make the 2D layout and ITF cross section two independent plots in
  a single html file.

  ```python
      # Create two independent plots
      fig1 = go.Figure(go.Scatter(x=[1,2], y=[1,2]))
      fig2 = go.Figure(go.Bar(x=[1,2], y=[2,1]))
      # Generate HTML divs
      div1 = pyo.plot(fig1, include_plotlyjs=True, output_type='div')
      div2 = pyo.plot(fig2, include_plotlyjs=False, output_type='div')
      # Combine in HTML
      html_content = f"<html><body>{div1}{div2}</body></html>"
      with open('two_plots.html', 'w') as f:
          f.write(html_content)
  ```

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
     - Total capacitance to substrate.
     - Total capacitance to same-layer.
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

1. Create the configuration file documentation in the docs/ directory.

   - MUST include a description of all configuration parameters.

## Tech Stack

- Python 3.14
- Astral uv for project management
- Rich-Click for CLI
- Plotly for visualization
- Use python's Pathlib whenever possible
