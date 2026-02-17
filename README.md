# Power/Ground Grid Netlist Generator

Generates an ASIC PG grid model from YAML + ITF inputs:

1. `output/pg_grid_netlist.sp` (SPICE RC + instances + `.MEAS`)
2. `output/pg_grid_visualization.html` (2D grid + stack cross-section)
3. `output/pg_grid_summary.txt` (ASCII report)

## Links

- Live demo: https://smprather.github.io/pg-grid-netlist-gen
- Configuration documentation: `docs/configuration.md`

## Requirements

- Python `>=3.14`
- `uv` (recommended)

## Setup

```bash
uv pip install -e .
```

## Run

```bash
./pg_grid_netlist_gen grid_specs.yaml
```

Common options:

```bash
./pg_grid_netlist_gen grid_specs.yaml --open-browser
./pg_grid_netlist_gen grid_specs.yaml --output-dir output
```

## Configuration

Primary inputs:

- `grid_specs.yaml`
- ITF file referenced by `itf.file` (for example `freepdk3_rctyp.itf`)

Highlights:

- Grid layer usage supports stripe (`g`) and staple (`s`) layers.
- Metal resistance uses ITF `RPSQ`.
- Via resistance uses ITF `RPV`.
- Standard-cell PG pins tap directly onto lowest metal grid nodes by pin `type` (`power` / `ground`), not by pin name.
- IR-drop measurements are generated for both output rise and output fall per instance, with edge mapping based on `standard_cells[].unateness`.
- Optional external cell model include: `standard_cells[].spice_netlist_file`.
- Netlist includes HSPICE probing for instance pin voltages via `.PROBE V(X*)`.

## Output Regeneration Hook

This repo includes a versioned pre-commit hook at:

- `.githooks/pre-commit`

It regenerates outputs from `grid_specs.yaml` and stages `output/` artifacts before each commit.

Enable it once per clone:

```bash
git config core.hooksPath .githooks
```

## Notes

- Default outputs are written under `output/`.
- The HTML visualization uses legend toggles for layers (no dropdown selector).
- The cross-section plot is included below the 2D grid view.
