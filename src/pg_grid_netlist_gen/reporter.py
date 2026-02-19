"""Functions for generating a detailed ASCII report of the generated grid."""

from __future__ import annotations

from collections import Counter, defaultdict

from pg_grid_netlist_gen.config import Config
from pg_grid_netlist_gen.geometry import Grid


def _count_chain_load_elements(grid: Grid, config: Config) -> tuple[int, int]:
    """Return (resistor_count, capacitor_count) for chain load elements."""
    chain_cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
    output_pin = next((p.name for p in chain_cell_cfg.pins if p.direction == "output"), None)
    input_pin = next(
        (p.name for p in chain_cell_cfg.pins if p.direction == "input" and p.type == "signal"),
        None,
    )
    if not output_pin or not input_pin:
        return 0, 0

    loads_cfg = config.spice_netlist.cell_chains.cell_output_loads
    in_chain_segs = max(1, loads_cfg.in_chain.number_pi_segments)
    end_chain_segs = max(1, loads_cfg.end_of_chain.number_pi_segments)

    res_count = 0
    cap_count = 0
    for cell in grid.cells:
        output_net = cell.pin_connections.get(output_pin)
        if not output_net or not str(output_net).lower().startswith("chain_"):
            continue

        is_last_in_chain = not any(output_net == c.pin_connections.get(input_pin) for c in grid.cells)
        if is_last_in_chain:
            res_count += end_chain_segs
            cap_count += end_chain_segs + 1
        else:
            res_count += in_chain_segs
            cap_count += in_chain_segs + 1
    return res_count, cap_count


def generate_report(grid: Grid, config: Config) -> str:
    """Generates a detailed, multi-line ASCII report of the grid."""
    
    report_lines = [
        "--- Power Grid Generation Report ---",
        "",
        "** Grid Dimensions **",
        f"  - Technology: {config.beol_stack.technology}",
        f"  - Effective Seed: {config.effective_seed if config.effective_seed is not None else config.random_seed_value}",
    ]

    # Grid Size Calculation
    grid_size_y_um = config.grid.size["rows"] * config.standard_cell_placement.row_height
    grid_size_x_um = config.grid.size["sites"] * config.standard_cell_placement.site_width
    report_lines.append(f"  - Size (X x Y): {grid_size_x_um:.2f} um x {grid_size_y_um:.2f} um")
    report_lines.append("")

    # Component Counts
    report_lines.append("** Component Counts **")
    report_lines.append(f"  - Metal Stripes: {len(grid.stripes)}")
    report_lines.append(f"  - Metal Segments (after breaking): {len(grid.segments)}")
    report_lines.append(f"  - Vias: {len(grid.vias)}")
    report_lines.append(f"  - Netlist Nodes: {len(grid.nodes)}")
    report_lines.append(f"  - Placed Chain Cells: {len(grid.cells)}")
    report_lines.append(f"  - Placed Dcap Cells: {len(grid.dcap_cells)}")
    dcap_cfg = config.standard_cell_placement.dcap_cells
    if dcap_cfg is not None and dcap_cfg.enabled and grid.dcap_cells:
        dcap_cell_cfg = config.get_cell_by_name(dcap_cfg.cell)
        dcap_area = config.distance_to_nm(dcap_cell_cfg.width) * config.distance_to_nm(dcap_cell_cfg.height)
        grid_size_x_nm = config.grid.size["sites"] * config.distance_to_nm(config.standard_cell_placement.site_width)
        grid_size_y_nm = config.grid.size["rows"] * config.distance_to_nm(config.standard_cell_placement.row_height)
        grid_area = grid_size_x_nm * grid_size_y_nm
        density = (len(grid.dcap_cells) * dcap_area / grid_area) * 100.0 if grid_area > 0 else 0.0
        report_lines.append(f"  - Dcap Cell Density: {density:.2f}%")
    else:
        report_lines.append(f"  - Dcap Cell Density: 0.00%")
    report_lines.append("")

    # Netlist Details
    grid_resistors = len(grid.segments) + len(grid.vias)
    source_resistors = len(grid.plocs) if grid.plocs else 2

    chain_resistors, chain_caps = _count_chain_load_elements(grid, config)
    num_resistors = grid_resistors + source_resistors + chain_resistors

    node_cap: dict[str, float] = defaultdict(float)
    for seg in grid.segments:
        c_total = seg.cap_plate + seg.cap_fringe
        half_c = c_total / 2.0
        node_cap[seg.node_a.name] += half_c
        node_cap[seg.node_b.name] += half_c
    grid_node_caps = sum(1 for cap in node_cap.values() if cap > 0)
    num_capacitors = grid_node_caps + chain_caps

    total_cap_pf = sum(s.cap_plate + s.cap_fringe for s in grid.segments) * 1e12
    power_net = config.pg_nets.power.name
    ground_net = config.pg_nets.ground.name
    
    report_lines.append("** Netlist Details **")
    report_lines.append(f"  - Resistors: {num_resistors}")
    report_lines.append(f"  - Capacitors: {num_capacitors}")
    report_lines.append(f"  - Total Grid Capacitance (to substrate): {total_cap_pf:.4f} pF")
    report_lines.append(f"  - Power Net: {power_net}")
    report_lines.append(f"  - Ground Net: {ground_net}")
    chain_cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
    input_pin = next(
        (p.name for p in chain_cell_cfg.pins if p.type == "signal" and p.direction == "input"),
        None,
    )
    if input_pin:
        chain_count = len(
            {
                c.pin_connections.get(input_pin)
                for c in grid.cells
                if str(c.pin_connections.get(input_pin, "")).startswith("CHAIN_")
            }
        )
    else:
        chain_count = 0
    report_lines.append(f"  - Total Chain Instance Count: {len(grid.cells)}")
    report_lines.append(f"  - Total Chains: {chain_count}")
    report_lines.append("")

    # Per-Layer Counts
    report_lines.append("** Details Per Layer **")
    stripe_counts = Counter(s.layer for s in grid.stripes)
    via_counts = Counter(v.via_layer for v in grid.vias)
    
    for layer in config.beol_stack.layers:
        if layer.type == 'metal' and layer.name in stripe_counts:
            report_lines.append(f"  - {layer.name}: {stripe_counts[layer.name]} stripes")
        if layer.type == 'via' and layer.name in via_counts:
            report_lines.append(f"  - {layer.name}: {via_counts[layer.name]} vias")

    report_lines.append("\n--- End of Report ---")
    
    return "\n".join(report_lines)
