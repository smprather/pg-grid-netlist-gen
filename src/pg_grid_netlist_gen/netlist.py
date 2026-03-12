"""SPICE netlist writer."""

from __future__ import annotations

import random
from collections import defaultdict
from pathlib import Path
from typing import TextIO

from pg_grid_netlist_gen.config import Config, GaussianVariationConfig
from pg_grid_netlist_gen.geometry import Grid


def _first_pin_name_of_type(cell_cfg, pin_type: str) -> str | None:
    return next((p.name for p in cell_cfg.pins if p.type == pin_type), None)


def _sample_variation(cfg: GaussianVariationConfig, rng: random.Random) -> float:
    """Sample a value from a Gaussian distribution, clamped to [floor, ceiling]."""
    value = rng.gauss(cfg.nominal, cfg.sigma)
    return max(cfg.floor, min(cfg.ceiling, value))


def _write_chain_stimulus(f: TextIO, grid: Grid, config: Config, rng: random.Random):
    f.write("* === Chain Input Stimulus ===\n")
    stimulus_cfg = config.spice_netlist.cell_chains.chain_input_stimulus
    power_voltage = config.pg_nets.power.voltage

    chain_cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.cell)
    input_pin = next(
        (p.name for p in chain_cell_cfg.pins if p.direction == "input" and p.type == "signal"),
        None,
    )
    if input_pin is None:
        return

    # Find all unique chain input nets
    input_nets = {c.pin_connections.get(input_pin) for c in grid.cells if c.pin_connections.get(input_pin, '').startswith("CHAIN_")}

    tu = config.units.time
    for net in sorted(input_nets):
        if net:
            tt = _sample_variation(stimulus_cfg.transition_time, rng)
            delay = _sample_variation(stimulus_cfg.initial_delay, rng)
            half_period = stimulus_cfg.period / 2
            f.write(
                f"V_{net} {net} 0 PULSE(0 {power_voltage} "
                f"{delay}{tu} {tt}{tu} {tt}{tu} "
                f"{half_period}{tu} {stimulus_cfg.period}{tu})\n"
            )
    f.write("\n")


def _write_pi_interconnect(
    f: TextIO, from_net: str, to_net: str, icn_cfg, config: Config,
) -> None:
    """Write a multi-segment pi-model interconnect from from_net to to_net.

    For N segments the topology is:
        from_net ─┬─ C ─┬─ R ─┬─ C ─┬─ R ─┬─ C ─┬─ to_net
                  GND        GND        GND
    with N resistors and N+1 capacitors.
    """
    r = icn_cfg.resistance
    c_f = _cap_value_to_f(icn_cfg.capacitance, config.units.capacitance)
    n_seg = max(1, icn_cfg.number_pi_segments)
    r_per = r / n_seg
    c_per_f = c_f / (n_seg + 1)

    node = from_net
    for seg_idx in range(n_seg):
        f.write(f"C_{from_net}_piC_{seg_idx} {node} 0 {_format_value(c_per_f)}\n")
        if seg_idx < n_seg - 1:
            next_node = f"{from_net}_pi_{seg_idx + 1}"
        else:
            next_node = to_net
        f.write(f"R_{from_net}_piR_{seg_idx} {node} {next_node} {_format_value(r_per)}\n")
        node = next_node
    # Final cap at the to_net end
    f.write(f"C_{from_net}_piC_{n_seg} {node} 0 {_format_value(c_per_f)}\n")


def _write_end_of_chain_load(
    f: TextIO, output_net: str, load_cfg, config: Config,
) -> None:
    """Write a simple RC load to ground on the last cell's output."""
    r = load_cfg.resistance
    c_f = _cap_value_to_f(load_cfg.capacitance, config.units.capacitance)
    load_node = f"{output_net}_eoc"
    f.write(f"R_{output_net}_eoc {output_net} {load_node} {_format_value(r)}\n")
    f.write(f"C_{output_net}_eoc {load_node} 0 {_format_value(c_f)}\n")


def _write_chain_loads(f: TextIO, grid: Grid, config: Config):
    f.write("* === Chain Interconnects and Loads ===\n")
    chains_cfg = config.spice_netlist.cell_chains

    chain_cell_cfg = config.get_cell_by_name(chains_cfg.cell)
    output_pin = next(
        (p.name for p in chain_cell_cfg.pins if p.direction == "output" and p.type == "signal"),
        None,
    )
    input_pin = next(
        (p.name for p in chain_cell_cfg.pins if p.direction == "input" and p.type == "signal"),
        None,
    )
    if output_pin is None or input_pin is None:
        return

    for cell in grid.cells:
        output_net = cell.pin_connections.get(output_pin)
        if not output_net or not output_net.lower().startswith("chain_"):
            continue

        # End-of-chain outputs are named CHAIN_<n>_OUT
        if output_net.endswith("_OUT"):
            _write_end_of_chain_load(f, output_net, chains_cfg.end_of_chain_load, config)
        else:
            # Interconnect links: output is chain_<n>_link_<j>,
            # next cell's input is chain_<n>_link_<j>_in
            to_net = f"{output_net}_in"
            _write_pi_interconnect(f, output_net, to_net, chains_cfg.interconnect, config)
    f.write("\n")


def write_netlist(grid: Grid, config: Config, output_path: str | Path, rng: random.Random) -> None:
    """Write the SPICE netlist to a file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        _write_header(f, config)
        _write_user_defined_lines(f, config)
        _write_sources(f, grid, config)
        _write_chain_stimulus(f, grid, config, rng)
        _write_subckt_stubs(f, config)
        _write_resistors(f, grid)
        _write_via_resistors(f, grid)
        _write_capacitors(f, grid, config)
        _write_cells(f, grid, config)
        _write_chain_loads(f, grid, config)
        _write_analysis(f)
        _write_measurements(f, grid, config)
        f.write("\n.end\n")


def _write_user_defined_lines(f: TextIO, config: Config) -> None:
    lines = config.spice_netlist.user_defined_lines
    if lines:
        f.write("* === User Defined Lines ===\n")
        for line in lines:
            f.write(f"{line}\n")
        f.write("\n")


def _write_header(f: TextIO, config: Config) -> None:
    f.write(f"* Power Grid Netlist - {config.beol_stack.technology}\n")
    f.write("* Generated by pg_grid_netlist_gen\n")

    grid_size_y = config.grid.size["rows"] * config.standard_cell_placement.row_height
    grid_size_x = config.grid.size["sites"] * config.standard_cell_placement.site_width
    dist_unit = config.units.distance

    f.write(
        f"* Grid: {grid_size_x:.2f}{dist_unit} x "
        f"{grid_size_y:.2f}{dist_unit}\n"
    )
    f.write(f"* Layers: {len(config.beol_stack.layers)}\n")
    f.write("\n")

def _write_sources(f: TextIO, grid: Grid, config: Config) -> None:
    f.write("* === Voltage Sources ===\n")
    power_net = config.pg_nets.power.name
    power_voltage = config.pg_nets.power.voltage
    ground_net = config.pg_nets.ground.name
    if grid.plocs:
        for i, ploc in enumerate(grid.plocs):
            if ploc.net == power_net:
                src = f"{power_net}_src_{i}"
                f.write(f"V{power_net}_{i} {src} 0 DC {power_voltage}\n")
                f.write(f"R{power_net}_src_{i} {src} {ploc.node_name} 1e-3\n")
            elif ploc.net == ground_net:
                src = f"{ground_net}_src_{i}"
                f.write(f"V{ground_net}_{i} {src} 0 DC 0\n")
                f.write(f"R{ground_net}_src_{i} {src} {ploc.node_name} 1e-3\n")
    else:
        # Fallback when no PLOCs are generated.
        power_node = next((n.name for n in grid.nodes.values() if n.net == power_net), power_net)
        ground_node = next((n.name for n in grid.nodes.values() if n.net == ground_net), ground_net)
        f.write(f"V{power_net} {power_net}_src 0 DC {power_voltage}\n")
        f.write(f"V{ground_net} {ground_net}_src 0 DC 0\n")
        f.write(f"R{power_net}_src {power_net}_src {power_node} 1e-3\n")
        f.write(f"R{ground_net}_src {ground_net}_src {ground_node} 1e-3\n")
    f.write("\n")


def _write_subckt_stubs(f: TextIO, config: Config) -> None:
    f.write("* === Standard Cell Model Includes ===\n")
    included_files: set[str] = set()
    for cell_cfg in config.standard_cells:
        netlist_file = cell_cfg.spice_netlist_file
        if netlist_file not in included_files:
            f.write(f'.include "{netlist_file}"\n')
            included_files.add(netlist_file)
    f.write("\n")


def _write_resistors(f: TextIO, grid: Grid) -> None:
    f.write("* === Metal Segments (R) ===\n")
    for i, seg in enumerate(grid.segments):
        f.write(
            f"R_{seg.layer}_seg_{i} {seg.node_a.name} {seg.node_b.name} "
            f"{_format_value(seg.resistance)}\n"
        )
    f.write("\n")


def _write_via_resistors(f: TextIO, grid: Grid) -> None:
    f.write("* === Via Resistances ===\n")
    for i, via in enumerate(grid.vias):
        f.write(
            f"R_{via.via_layer}_{i} {via.node_top.name} {via.node_bot.name} "
            f"{_format_value(via.resistance)}\n"
        )
    f.write("\n")


def _write_capacitors(f: TextIO, grid: Grid, config: Config) -> None:
    f.write("* === Node Capacitances ===\n")
    ground_net = config.pg_nets.ground.name
    node_cap: dict[str, float] = defaultdict(float)
    for seg in grid.segments:
        c_total = seg.cap_plate + seg.cap_fringe + seg.cap_coupling
        half_c = c_total / 2.0
        node_cap[seg.node_a.name] += half_c
        node_cap[seg.node_b.name] += half_c

    for node_name, cap in sorted(node_cap.items()):
        if cap > 0:
            f.write(f"C_{node_name} {node_name} {ground_net} {_format_value(cap)}\n")
    f.write("\n")


def _write_cells(f: TextIO, grid: Grid, config: Config) -> None:
    f.write("* === Standard Cell Instances ===\n")

    for cell in grid.cells:
        cell_cfg = next(c for c in config.standard_cells if c.name == cell.cell_name)
        pin_nodes: list[str] = []
        for pin in cell_cfg.spice_port_order:
            metal_node = cell.pin_connections.get(pin, f"NC_{pin}")
            pin_nodes.append(metal_node)
        pins_str = " ".join(pin_nodes)
        f.write(f"{cell.instance_name} {pins_str} {cell.cell_name}\n")

    if grid.dcap_cells:
        f.write("\n* === Dcap Cell Instances ===\n")
        for cell in grid.dcap_cells:
            cell_cfg = next(c for c in config.standard_cells if c.name == cell.cell_name)
            pin_nodes: list[str] = []
            for pin in cell_cfg.spice_port_order:
                metal_node = cell.pin_connections.get(pin, f"NC_{pin}")
                pin_nodes.append(metal_node)
            pins_str = " ".join(pin_nodes)
            f.write(f"{cell.instance_name} {pins_str} {cell.cell_name}\n")
    f.write("\n")

def _write_analysis(f: TextIO) -> None:
    f.write("* === Analysis ===\n")
    f.write(".PROBE V(X*)\n")
    f.write("\n")

def _write_long_measure(f: TextIO, name: str, expr: str, max_line: int = 800) -> None:
    """Write a .measure PARAM statement, splitting with HSPICE '+' continuation lines."""
    prefix = f".measure tran {name} PARAM='"
    suffix = "'\n"
    # If it fits on one line, just write it.
    if len(prefix) + len(expr) + len(suffix) <= max_line:
        f.write(f"{prefix}{expr}{suffix}")
        return
    # Split the expression at ', ' or ' + ' boundaries to stay under max_line.
    f.write(prefix)
    line_len = len(prefix)
    i = 0
    while i < len(expr):
        # Find next split point (after ', ' or ' + ')
        next_comma = expr.find(", ", i + 1)
        next_plus = expr.find(" + ", i + 1)
        candidates = [c for c in (next_comma, next_plus) if c > i]
        if not candidates:
            # No more split points; write the rest.
            f.write(expr[i:])
            break
        split_at = min(candidates)
        # Include the delimiter in the current chunk.
        delim_len = 2 if expr[split_at : split_at + 2] == ", " else 3
        chunk_end = split_at + delim_len
        chunk = expr[i:chunk_end]
        if line_len + len(chunk) > max_line and i > 0:
            f.write("\n+ ")
            line_len = 2
        f.write(chunk)
        line_len += len(chunk)
        i = chunk_end
    f.write(suffix)


def _write_measurements(f: TextIO, grid: Grid, config: Config) -> None:
    f.write("* === IR Drop Measurements ===\n")
    # Only measure chain cells (dcap cells have no signal pins → no measurements)
    cell_cfg_by_name = {c.name: c for c in config.standard_cells}

    window_cfg = config.spice_netlist.ir_drop_measurement["averaging_window"]
    power_voltage = config.pg_nets.power.voltage
    start_pct = float(window_cfg["start"]) / 100.0
    end_pct = float(window_cfg["end"]) / 100.0
    rise_start_level = power_voltage * start_pct
    rise_end_level = power_voltage * end_pct
    fall_start_level = power_voltage * (1.0 - start_pct)
    fall_end_level = power_voltage * (1.0 - end_pct)

    measure_names: list[str] = []
    probe_nodes: set[str] = set()
    probe_diffs: set[tuple[str, str]] = set()

    for cell in grid.cells:
        cell_cfg = cell_cfg_by_name.get(cell.cell_name)
        if cell_cfg is None:
            continue

        power_pin_name = _first_pin_name_of_type(cell_cfg, "power")
        gnd_pin_name = _first_pin_name_of_type(cell_cfg, "ground")
        if not power_pin_name or not gnd_pin_name:
            continue

        vdd_node = cell.pin_connections.get(power_pin_name)
        vss_node = cell.pin_connections.get(gnd_pin_name)
        if not vdd_node or not vss_node:
            continue

        probe_nodes.add(vdd_node)
        probe_nodes.add(vss_node)
        probe_diffs.add((vdd_node, vss_node))

        input_pin_name = next(
            (p.name for p in cell_cfg.pins if p.type == "signal" and p.direction == "input"),
            None,
        )
        output_pin_name = next(
            (p.name for p in cell_cfg.pins if p.type == "signal" and p.direction == "output"),
            None,
        )
        in_node = cell.pin_connections.get(input_pin_name) if input_pin_name else None
        out_node = cell.pin_connections.get(output_pin_name) if output_pin_name else None

        if in_node:
            probe_nodes.add(in_node)
        if out_node:
            probe_nodes.add(out_node)

        for out_edge in ("RISE", "FALL"):
            measure_name = f"IR_DROP_{cell.instance_name}_{out_edge}"
            if in_node and out_node:
                if cell_cfg.unateness == "negative":
                    in_edge = "FALL" if out_edge == "RISE" else "RISE"
                else:
                    in_edge = out_edge

                start_level = rise_start_level if in_edge == "RISE" else fall_start_level
                end_level = rise_end_level if out_edge == "RISE" else fall_end_level

                t_start = f"T_START_{cell.instance_name}_{out_edge}"
                t_end = f"T_END_{cell.instance_name}_{out_edge}"
                f.write(
                    f".measure tran {t_start} WHEN v({in_node})={_format_value(start_level)} {in_edge}=1\n"
                )
                f.write(
                    f".measure tran {t_end} WHEN v({out_node})={_format_value(end_level)} {out_edge}=1\n"
                )
                f.write(
                    f".measure tran {measure_name} "
                    f"AVG '(v({vdd_node}) - v({vss_node}))/{power_voltage}' FROM='{t_start}' TO='{t_end}'\n"
                )
            else:
                # Fallback if signal pins are unavailable for this cell.
                measure_time = config.spice_netlist.cell_chains.chain_input_stimulus.period / 2.0
                f.write(
                    f".measure tran {measure_name} "
                    f"find '(v({vdd_node}) - v({vss_node}))/{power_voltage}' at={measure_time}{config.units.time}\n"
                )
            measure_names.append(measure_name)

    f.write("\n")
    if measure_names:
        if len(measure_names) == 1:
            f.write(f".measure tran MAX_IR_DROP PARAM='{measure_names[0]}'\n")
            f.write(f".measure tran AVG_IR_DROP PARAM='{measure_names[0]}'\n")
        else:
            _write_long_measure(f, "MAX_IR_DROP", f"max({', '.join(measure_names)})")
            avg_expr = f"({' + '.join(measure_names)})/{len(measure_names)}"
            _write_long_measure(f, "AVG_IR_DROP", avg_expr)
    f.write("\n")

    # Single-ended probes for every voltage node used in .measure lines.
    f.write("* === Measurement Probes ===\n")
    for node in sorted(probe_nodes):
        f.write(f".probe v({node})\n")
    f.write("\n")

    # Differential probes mirroring the AVG expressions.
    f.write("* === Differential Probes ===\n")
    for vdd_n, vss_n in sorted(probe_diffs):
        f.write(f".probe v({vdd_n}, {vss_n})\n")
    f.write("\n")


def _format_value(value: float) -> str:
    """Format a value with engineering notation for SPICE."""
    if value == 0:
        return "0"
    abs_val = abs(value)
    if abs_val >= 1e6:
        return f"{value:.4g}MEG"
    if abs_val >= 1e3:
        return f"{value:.4g}k"
    if abs_val >= 1:
        return f"{value:.4g}"
    if abs_val >= 1e-3:
        return f"{value * 1e3:.4g}m"
    if abs_val >= 1e-6:
        return f"{value * 1e6:.4g}u"
    if abs_val >= 1e-9:
        return f"{value * 1e9:.4g}n"
    if abs_val >= 1e-12:
        return f"{value * 1e12:.4g}p"
    if abs_val >= 1e-15:
        return f"{value * 1e15:.4g}f"
    return f"{value:.4e}"


def _cap_value_to_f(value: float, unit: str) -> float:
    unit_map = {
        "f": 1.0,
        "pf": 1e-12,
        "ff": 1e-15,
        "nf": 1e-9,
        "uf": 1e-6,
    }
    key = unit.strip().lower()
    if key not in unit_map:
        raise ValueError(
            f"Unsupported units.capacitance '{unit}'. Supported: f, pf, ff, nf, uf"
        )
    return float(value) * unit_map[key]
