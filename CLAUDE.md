# CLAUDE.md — Agent Instructions

Canonical spec: `AGENTS.md`. Config schema: `docs/configuration.md`.

## Tech Stack

- Python 3.14, Astral uv for project management
- Rich-Click for CLI, Plotly for visualization
- Use Pathlib whenever possible

## Project Structure

```
config.yaml                          # Primary YAML config
freepdk3_rctyp.itf                   # ITF technology file
src/pg_grid_netlist_gen/
  cli.py                             # Rich-Click CLI entry point
  config.py                          # YAML parsing & validation
  itf_parser.py                      # ITF file parser
  grid_builder.py                    # PG grid construction
  geometry.py                        # Geometry utilities
  physics.py                         # RC extraction formulas
  netlist.py                         # SPICE netlist generation
  visualize.py                       # Plotly visualization
  reporter.py                        # ASCII summary report
output/                              # Default output directory
  pg_grid_visualization.html         # 2D grid + ITF cross-section
  pg_grid_netlist.sp                 # SPICE RC + instances + .MEAS
  pg_grid_summary.txt                # ASCII metrics report
```

## Requirement Keywords

- MUST / MUST NOT: Required / forbidden.
- SHOULD: Strong recommendation.
- MAY: Optional.

## Definitions

- **Standard cell**: Logic-cell abstraction; vertically flipped on alternating rows.
- **Instance**: Placed instantiation of a standard cell.
- **ITF file**: Interconnect Technology Format — conductor and via layers.
- **PLOC**: Pad location — ideal supply-connection point.
- **Grid layer** (`grid`): Full-span horizontal or vertical metal stripes.
- **Staple layer** (`staple`): Local square metal patches with vias above and below.
- **Flight line**: Visualization line from instance output pin to another instance input pin.
- **IR drop**: `avg(V(power_pin) - V(ground_pin)) / VDD_nominal` — percent drop averaged over a configured transient window.
- **Site width**: Horizontal placement quantum; cell left edges MUST lie on integer multiples.
- **Unateness**: Inverting (negative) or non-inverting (positive). Affects IR-drop measurement edge mapping.
- **Dcap cell**: Decoupling capacitance cell — no signal pins, reduces IR drop locally.

## Layer-Usage Critical Rules

- Width/pitch are multiples of ITF WMIN/SMIN: `actual_width = width * WMIN`; `actual_pitch = pitch * (WMIN + SMIN)`.
- `grid.layer_usage.<LAYER>.type` MUST be `grid` or `staple`.
- The lowest ITF metal layer is **implicit** — MUST be omitted from `grid.layer_usage`.
- Implicit lowest layer: type=`grid`, width=`2.0 * WMIN`, pitch=`standard_cell_placement.row_height`.
- All explicit layer names MUST exist as ITF `CONDUCTOR` names.

## Netlist and Extraction Rules

- Grid MUST be constructed bottom-up.
- A via MUST be instantiated anywhere metal above and below intersect.
- Via shape MUST be square and satisfy ITF `AREA`.
- PLOC generation MUST start with power net, alternate horizontally between power/ground.
- PLOC not on intended top-layer stripe MUST snap to nearest valid stripe centerline.
- Where instance PG pin lands on lowest metal, break that segment into two at pin location. Exception: pin on segment edge taps pre-existing node.
- Add `.PROBE V(X*)` for all instance pin voltages.
- Netlist node names MUST begin with the associated net name.
- New nodes MUST be created at: VIA top/bottom; cell/dcap PG pin on lowest metal; PLOC connection point; metal stripe at grid edge.
- Capacitance: half at each segment end (pi model). All caps connect to ground net directly (e.g., `VSS`).

### RC Extraction

- Capacitance only on power-type PG nets. Substrate distance = `feol_thickness`.
- Resistance: `R = RPSQ * (L / W)`. Via resistance: ITF `RPV`.
- Segment capacitance MUST include: plate-to-substrate, fringe-to-substrate (isolated-wire), plate-to-plate (power to nearest same-layer ground).
- Via and staple-shape capacitance MUST NOT be included.

### Staple Resistance Model

Between adjacent routed layers (e.g., `M7-V6-M6-V5-M5`), three series elements:
1. `R(V6)` — ITF `RPV`.
2. `R(M6)` — one square: `RPSQ * 1.0`.
3. `R(V5)` — ITF `RPV`.

## Placement and Connectivity Rules

- Bottom ITF metal = horizontal. Higher layers alternate H/V by level.
- Signal pins at midpoint of declared side.
- Cell rows alternate vertical flip. Top/bottom edges on PG metal centerlines.
- Placement MUST honor row/site grids. MUST NOT overlap any cells.
- PG connectivity by pin **type** (power/ground), not pin name.
- Via min-space = via side * `grid.via_min_space_factor`.
- Visualizer: pack max min-spaced VIAs in overlap. Extraction: model as single via with `R = RPV / via_count`.

### Dcap Insertion

- Insert until `max_density_pct` (dcap area / total grid area). Separate legend group from chain cells.

### Chain Generation

- Random instance start, output-to-input across random instances. Each cell in one chain only.
- Up to `max_instance_count_per_chain` per chain; create chains until all instances assigned.
- Pi-model interconnect (`cell_chains.interconnect`) bridges output→input between consecutive cells. Simple RC load (`end_of_chain_load`) on last cell. Input: SPICE `PULSE`.

## Determinism

- Numeric seed: deterministic/reproducible. `"random"`: fresh seed per run.
- Single seeded RNG for all random subsystems. Report MUST include effective seed.

## Measurement Timing

- Per-instance IR-drop windows MUST use transition crossing events, not fixed times.
- Two windows per instance: `..._RISE` and `..._FALL`.
- Positive unateness: rise←input rise, fall←input fall.
- Negative unateness: rise←input fall, fall←input rise.
- Start: `averaging_window.start` % of VDD. Falling-edge: complement `VDD * (1 - start%)`.
- End: `averaging_window.end` % of VDD. Falling-edge: complement `VDD * (1 - end%)`.

## Validation

Generator MUST fail fast with clear error messages.

- YAML keys/types match schema. Layer names valid ITF conductors. Stack vertically connectable.
- `grid.layer_usage` MUST NOT include implicit lowest layer.
- `width > 0`, `pitch > 0`, `row_height > 0`, `site_width > 0`.
- `rows >= 1`, `sites >= 1`, `max_instance_count_per_chain >= 1`.
- Averaging window: `0 <= start < end <= 100`.
- `.subckt` ports in `spice_netlist_file` MUST exactly match declared pin names.

## Visualization

- Plotly legend for layer visibility (NO dropdown).
- Initial: instances + bottom 4 metal layers + their via layers visible. Rest: `legendonly`.
- Flight lines in own legend trace group.
- ITF cross-section below 2D render: hover labels, same-thickness layers share colors, FEOL bottom=y=0, substrate=3x FEOL, omit `_diel` layers, substrate label, no legend.
- `include_plotlyjs=True` only in first div. Two independent plots via `pyo.plot(output_type='div')`.

## Output Contract

1. `output/pg_grid_visualization.html` — 2D grid + ITF cross-section.
2. `output/pg_grid_netlist.sp` — PG RC network, instances, sources, stimulus, `.MEAS`.
3. `output/pg_grid_summary.txt` — seed, total PG cap, R count, C count, cap-to-substrate, cap-to-same-layer, chain instance count, chain count, dcap count, dcap density.

## Generation Order

1. Parse YAML and ITF.
2. Validate schema and technology/stack consistency.
3. Build PG grid bottom-up (including implicit lowest layer).
4. Insert vias and staple structures.
5. Place and orient standard-cell instances.
6. Generate/snap PLOCs.
7. Build instance chains and signal loading.
8. Extract RC network and generate SPICE.
9. Add `.MEAS` for per-instance, max, and average IR drop.
10. Generate visualization and ASCII summary.

## Documentation Rules

- `README.md` MUST link to live-demo at https://smprather.github.io/pg-grid-netlist-gen and to `docs/configuration.md`.
- `docs/configuration.md` MUST describe all configuration parameters.

## Formatting Rules (for AGENTS.md)

- Markdown tables MUST NOT exceed 110 columns per line.
- Use relative key paths, abbreviate: "Req"/"Y"/"N"/"Cond". Keep notes concise.
