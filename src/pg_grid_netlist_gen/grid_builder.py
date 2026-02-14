"""Core grid generation algorithm."""

from __future__ import annotations

from dataclasses import dataclass

from pg_grid_netlist_gen.config import Config
from pg_grid_netlist_gen.geometry import (
    CellPlacement,
    Grid,
    Node,
    Segment,
    Staple,
    Stripe,
    ViaConnection,
)
from pg_grid_netlist_gen.physics import (
    metal_resistance,
    total_segment_capacitance,
    via_resistance,
)


@dataclass
class ResolvedLayer:
    """BEOL + grid layer properties merged together."""

    name: str
    layer_type: str  # "grid" or "staple"
    direction: str  # "vertical" or "horizontal"
    actual_pitch_nm: float
    width_nm: float
    thickness_nm: float
    z_bottom_nm: float
    z_top_nm: float
    material_name: str
    conductivity_ms_m: float
    beol_type: str  # "metal" or "via"


def _compute_z_coordinates(config: Config) -> dict[str, tuple[float, float]]:
    """Walk the BEOL stack bottom-to-top, computing z_bottom and z_top for each layer.

    The BEOL stack is listed top-to-bottom in the YAML. We reverse it so SUB is at z=0.
    The substrate top surface is our z=0 reference.
    """
    layers_bottom_up = list(reversed(config.beol_stack.layers))
    z_coords: dict[str, tuple[float, float]] = {}
    z = 0.0

    for layer in layers_bottom_up:
        if layer.type == "substrate":
            z_coords[layer.name] = (-layer.thickness, 0.0)
            z = 0.0
        else:
            z_bottom = z
            z_top = z + layer.thickness
            z_coords[layer.name] = (z_bottom, z_top)
            z = z_top

    return z_coords


def _resolve_layers(config: Config) -> list[ResolvedLayer]:
    """Resolve grid layer properties from config + BEOL stack."""
    z_coords = _compute_z_coordinates(config)
    resolved = []

    for layer_name in config.grid_layer_order():
        grid_layer = config.grid.layers[layer_name]
        beol_layer = config.get_beol_layer(layer_name)

        actual_pitch = grid_layer.pitch * beol_layer.pitch
        width = grid_layer.width if grid_layer.width is not None else beol_layer.min_width
        z_bottom, z_top = z_coords[layer_name]

        mat_name = grid_layer.material or "copper"
        mat = config.get_material(mat_name)

        resolved.append(ResolvedLayer(
            name=layer_name,
            layer_type=grid_layer.type,
            direction=grid_layer.direction,
            actual_pitch_nm=actual_pitch,
            width_nm=width,
            thickness_nm=beol_layer.thickness,
            z_bottom_nm=z_bottom,
            z_top_nm=z_top,
            material_name=mat_name,
            conductivity_ms_m=mat.conductivity,
            beol_type=beol_layer.type,
        ))

    return resolved


def _node_name(layer: str, x: float, y: float) -> str:
    return f"{layer}_X_{int(x)}_Y_{int(y)}"


def _get_or_create_node(
    nodes: dict[str, Node],
    layer: str,
    x: float,
    y: float,
    z: float,
    net: str,
) -> Node:
    name = _node_name(layer, x, y)
    if name not in nodes:
        nodes[name] = Node(name=name, x=x, y=y, z=z, layer=layer, net=net)
    return nodes[name]


def _get_resolved(resolved: list[ResolvedLayer], name: str) -> ResolvedLayer:
    for rl in resolved:
        if rl.name == name:
            return rl
    raise ValueError(f"Resolved layer '{name}' not found")


def _find_via_layer_between(config: Config, top_layer: str, bot_layer: str) -> str | None:
    """Find the BEOL via layer between two metal layers."""
    beol_names = [l.name for l in config.beol_stack.layers]
    beol_types = {l.name: l.type for l in config.beol_stack.layers}

    try:
        top_idx = beol_names.index(top_layer)
        bot_idx = beol_names.index(bot_layer)
    except ValueError:
        return None

    if top_idx > bot_idx:
        top_idx, bot_idx = bot_idx, top_idx

    for i in range(top_idx + 1, bot_idx):
        if beol_types[beol_names[i]] == "via":
            return beol_names[i]

    return None


def _generate_stripes(
    resolved: list[ResolvedLayer],
    grid_size_x: float,
    grid_size_y: float,
    cycle: list[str],
) -> list[Stripe]:
    """Generate stripes for all grid-type layers."""
    stripes = []

    for rl in resolved:
        if rl.layer_type != "grid":
            continue

        if rl.direction == "vertical":
            n_stripes = int(grid_size_x / rl.actual_pitch_nm) + 1
            for i in range(n_stripes):
                pos = i * rl.actual_pitch_nm
                if pos > grid_size_x:
                    break
                net = cycle[i % len(cycle)]
                stripes.append(Stripe(
                    layer=rl.name, direction="vertical",
                    position=pos, start=0.0, end=grid_size_y,
                    width=rl.width_nm, net=net,
                ))
        else:
            n_stripes = int(grid_size_y / rl.actual_pitch_nm) + 1
            for i in range(n_stripes):
                pos = i * rl.actual_pitch_nm
                if pos > grid_size_y:
                    break
                net = cycle[i % len(cycle)]
                stripes.append(Stripe(
                    layer=rl.name, direction="horizontal",
                    position=pos, start=0.0, end=grid_size_x,
                    width=rl.width_nm, net=net,
                ))

    return stripes


def _generate_staples(
    resolved: list[ResolvedLayer],
    config: Config,
    grid_size_x: float,
    grid_size_y: float,
    cycle: list[str],
    stripes: list[Stripe],
) -> list[Staple]:
    """Generate staple pads for all staple-type layers."""
    staples = []
    layer_order = config.grid_layer_order()

    # Build indexes: layer -> set of positions, and layer -> position -> net
    stripe_v_positions: dict[str, dict[float, str]] = {}  # layer -> {x_pos -> net} for vertical
    stripe_h_positions: dict[str, dict[float, str]] = {}  # layer -> {y_pos -> net} for horizontal
    for s in stripes:
        if s.direction == "vertical":
            stripe_v_positions.setdefault(s.layer, {})[s.position] = s.net
        else:
            stripe_h_positions.setdefault(s.layer, {})[s.position] = s.net

    for rl in resolved:
        if rl.layer_type != "staple":
            continue

        idx = layer_order.index(rl.name)
        above_name = layer_order[idx - 1] if idx > 0 else None
        below_name = layer_order[idx + 1] if idx < len(layer_order) - 1 else None

        # Staple size = 2 × min_width of the layer above
        if above_name:
            above_rl = _get_resolved(resolved, above_name)
            staple_size = 2.0 * above_rl.width_nm
        else:
            staple_size = 2.0 * rl.width_nm

        # Collect x and y positions from adjacent layers
        x_positions: set[float] = set()
        y_positions: set[float] = set()

        for adj_name in (above_name, below_name):
            if adj_name is None:
                continue
            if adj_name in stripe_v_positions:
                x_positions.update(stripe_v_positions[adj_name].keys())
            if adj_name in stripe_h_positions:
                y_positions.update(stripe_h_positions[adj_name].keys())

        # Fill in from pitch if we don't have positions from stripes
        if not x_positions:
            n = int(grid_size_x / rl.actual_pitch_nm) + 1
            x_positions = {i * rl.actual_pitch_nm for i in range(n) if i * rl.actual_pitch_nm <= grid_size_x}
        if not y_positions:
            n = int(grid_size_y / rl.actual_pitch_nm) + 1
            y_positions = {i * rl.actual_pitch_nm for i in range(n) if i * rl.actual_pitch_nm <= grid_size_y}

        # Build net lookup functions using indexed data
        above_v = stripe_v_positions.get(above_name, {}) if above_name else {}
        above_h = stripe_h_positions.get(above_name, {}) if above_name else {}
        below_v = stripe_v_positions.get(below_name, {}) if below_name else {}
        below_h = stripe_h_positions.get(below_name, {}) if below_name else {}

        for x in sorted(x_positions):
            for y in sorted(y_positions):
                # Determine net from adjacent stripe positions
                net = (
                    above_v.get(x)
                    or above_h.get(y)
                    or below_v.get(x)
                    or below_h.get(y)
                )
                if net is None:
                    pos_idx = int(x / rl.actual_pitch_nm) if rl.direction == "vertical" else int(y / rl.actual_pitch_nm)
                    net = cycle[pos_idx % len(cycle)]

                staples.append(Staple(
                    layer=rl.name, x=x, y=y, size=staple_size, net=net,
                ))

    return staples


def _place_vias(
    config: Config,
    resolved: list[ResolvedLayer],
    stripes: list[Stripe],
    staples: list[Staple],
    nodes: dict[str, Node],
    z_coords: dict[str, tuple[float, float]],
) -> list[ViaConnection]:
    """Place vias between adjacent grid layers using indexed lookups."""
    vias = []
    layer_order = config.grid_layer_order()

    # Build indexes for fast lookup
    # Stripes: layer -> direction -> {position -> net}
    stripe_index: dict[str, dict[str, dict[float, str]]] = {}
    for s in stripes:
        stripe_index.setdefault(s.layer, {}).setdefault(s.direction, {})[s.position] = s.net

    # Staples: layer -> {(x, y) -> net}
    staple_index: dict[str, dict[tuple[float, float], str]] = {}
    # Also: layer -> {x -> set of y} and layer -> {y -> set of x} for fast position lookups
    staple_by_x: dict[str, dict[float, set[float]]] = {}
    staple_by_y: dict[str, dict[float, set[float]]] = {}
    for s in staples:
        staple_index.setdefault(s.layer, {})[(s.x, s.y)] = s.net
        staple_by_x.setdefault(s.layer, {}).setdefault(s.x, set()).add(s.y)
        staple_by_y.setdefault(s.layer, {}).setdefault(s.y, set()).add(s.x)

    for i in range(len(layer_order) - 1):
        top_name = layer_order[i]
        bot_name = layer_order[i + 1]
        top_rl = _get_resolved(resolved, top_name)
        bot_rl = _get_resolved(resolved, bot_name)

        via_layer_name = _find_via_layer_between(config, top_name, bot_name)
        if via_layer_name is None:
            continue

        via_beol = config.get_beol_layer(via_layer_name)
        conductivity = top_rl.conductivity_ms_m

        top_si = stripe_index.get(top_name, {})
        bot_si = stripe_index.get(bot_name, {})
        top_sti = staple_index.get(top_name, {})
        bot_sti = staple_index.get(bot_name, {})

        if top_rl.layer_type == "grid" and bot_rl.layer_type == "grid":
            # Grid-to-grid: find crossings of perpendicular stripes with matching nets
            top_v = top_si.get("vertical", {})
            top_h = top_si.get("horizontal", {})
            bot_v = bot_si.get("vertical", {})
            bot_h = bot_si.get("horizontal", {})

            # Case: top=vertical, bot=horizontal
            for tx, tnet in top_v.items():
                for by, bnet in bot_h.items():
                    if tnet == bnet:
                        _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                                 via_beol, top_rl, bot_rl, tx, by, tnet, conductivity)

            # Case: top=horizontal, bot=vertical
            for ty, tnet in top_h.items():
                for bx, bnet in bot_v.items():
                    if tnet == bnet:
                        _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                                 via_beol, top_rl, bot_rl, bx, ty, tnet, conductivity)

        elif top_rl.layer_type == "grid" and bot_rl.layer_type == "staple":
            # For each top stripe, find staples that lie on it
            for direction, positions in top_si.items():
                for pos, tnet in positions.items():
                    if direction == "vertical":
                        # Stripe at x=pos; find staples with matching x
                        for y in sorted(staple_by_x.get(bot_name, {}).get(pos, set())):
                            snet = bot_sti.get((pos, y))
                            if snet == tnet:
                                _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                                         via_beol, top_rl, bot_rl, pos, y, tnet, conductivity)
                    else:
                        # Stripe at y=pos; find staples with matching y
                        for x in sorted(staple_by_y.get(bot_name, {}).get(pos, set())):
                            snet = bot_sti.get((x, pos))
                            if snet == tnet:
                                _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                                         via_beol, top_rl, bot_rl, x, pos, tnet, conductivity)

        elif top_rl.layer_type == "staple" and bot_rl.layer_type == "grid":
            for direction, positions in bot_si.items():
                for pos, bnet in positions.items():
                    if direction == "vertical":
                        for y in sorted(staple_by_x.get(top_name, {}).get(pos, set())):
                            snet = top_sti.get((pos, y))
                            if snet == bnet:
                                _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                                         via_beol, top_rl, bot_rl, pos, y, bnet, conductivity)
                    else:
                        for x in sorted(staple_by_y.get(top_name, {}).get(pos, set())):
                            snet = top_sti.get((x, pos))
                            if snet == bnet:
                                _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                                         via_beol, top_rl, bot_rl, x, pos, bnet, conductivity)

        elif top_rl.layer_type == "staple" and bot_rl.layer_type == "staple":
            # Staple-to-staple: match by position and net
            # Iterate over the smaller set
            for (x, y), tnet in top_sti.items():
                bnet = bot_sti.get((x, y))
                if bnet == tnet:
                    _add_via(vias, nodes, top_name, bot_name, via_layer_name,
                             via_beol, top_rl, bot_rl, x, y, tnet, conductivity)

    return vias


def _add_via(
    vias: list[ViaConnection],
    nodes: dict[str, Node],
    top_name: str,
    bot_name: str,
    via_layer_name: str,
    via_beol,
    top_rl: ResolvedLayer,
    bot_rl: ResolvedLayer,
    x: float,
    y: float,
    net: str,
    conductivity: float,
) -> None:
    node_top = _get_or_create_node(nodes, top_name, x, y, top_rl.z_bottom_nm, net)
    node_bot = _get_or_create_node(nodes, bot_name, x, y, bot_rl.z_bottom_nm, net)
    r = via_resistance(via_beol.thickness, via_beol.min_width, conductivity)
    vias.append(ViaConnection(
        node_top=node_top, node_bot=node_bot,
        via_layer=via_layer_name,
        width=via_beol.min_width, height=via_beol.thickness,
        net=net, resistance=r,
    ))


def _segment_stripes(
    resolved: list[ResolvedLayer],
    config: Config,
    stripes: list[Stripe],
    nodes: dict[str, Node],
) -> list[Segment]:
    """Break stripes into segments at via connection points."""
    segments = []

    # Build index: (layer, position_along_stripe_axis, net) -> [via_coords_along_length]
    stripe_via_points: dict[tuple[str, float, str], list[float]] = {}
    for stripe in stripes:
        key = (stripe.layer, stripe.position, stripe.net)
        stripe_via_points.setdefault(key, [])

    for node in nodes.values():
        # Check if this node's layer has stripes
        for stripe in stripes:
            if stripe.layer != node.layer or stripe.net != node.net:
                continue
            if stripe.direction == "vertical" and stripe.position == node.x:
                key = (stripe.layer, stripe.position, stripe.net)
                stripe_via_points.setdefault(key, []).append(node.y)
                break  # each node only matches one stripe per layer
            elif stripe.direction == "horizontal" and stripe.position == node.y:
                key = (stripe.layer, stripe.position, stripe.net)
                stripe_via_points.setdefault(key, []).append(node.x)
                break

    oxide_mat = config.get_material("oxide")

    for stripe in stripes:
        rl = _get_resolved(resolved, stripe.layer)
        key = (stripe.layer, stripe.position, stripe.net)
        via_coords = sorted(set(stripe_via_points.get(key, [])))

        all_points = sorted(set([stripe.start] + via_coords + [stripe.end]))

        for j in range(len(all_points) - 1):
            c_start = all_points[j]
            c_end = all_points[j + 1]
            length = c_end - c_start
            if length <= 0:
                continue

            if stripe.direction == "vertical":
                x = stripe.position
                node_a = _get_or_create_node(nodes, stripe.layer, x, c_start, rl.z_bottom_nm, stripe.net)
                node_b = _get_or_create_node(nodes, stripe.layer, x, c_end, rl.z_bottom_nm, stripe.net)
            else:
                y = stripe.position
                node_a = _get_or_create_node(nodes, stripe.layer, c_start, y, rl.z_bottom_nm, stripe.net)
                node_b = _get_or_create_node(nodes, stripe.layer, c_end, y, rl.z_bottom_nm, stripe.net)

            r = metal_resistance(length, rl.width_nm, rl.thickness_nm, rl.conductivity_ms_m)
            cp, cf, _ = total_segment_capacitance(
                rl.width_nm, length, rl.thickness_nm,
                rl.z_bottom_nm,
                oxide_mat.relative_permittivity,
            )

            segments.append(Segment(
                node_a=node_a, node_b=node_b,
                layer=stripe.layer, width=rl.width_nm,
                thickness=rl.thickness_nm, length=length,
                net=stripe.net, resistance=r,
                cap_plate=cp, cap_fringe=cf,
            ))

    return segments


def _place_standard_cells(
    config: Config,
    resolved: list[ResolvedLayer],
    stripes: list[Stripe],
    nodes: dict[str, Node],
) -> list[CellPlacement]:
    """Place standard cells along M1 stripes."""
    cells = []
    layer_order = config.grid_layer_order()

    m1_name = layer_order[-1]
    m1_rl = _get_resolved(resolved, m1_name)
    m1_stripes = sorted(
        [s for s in stripes if s.layer == m1_name],
        key=lambda s: s.position,
    )

    if len(m1_stripes) < 2:
        return cells

    cell_cfg = config.standard_cells[0]
    cell_width = cell_cfg.size["x"]
    distance_apart = cell_cfg.distance_apart_um

    nc_counter = 0

    for row_idx in range(len(m1_stripes) - 1):
        stripe_a = m1_stripes[row_idx]
        stripe_b = m1_stripes[row_idx + 1]

        y_top = stripe_a.position
        y_bot = stripe_b.position
        row_height = y_bot - y_top

        flipped = row_idx % 2 == 1

        if not flipped:
            vdd_stripe = stripe_a
            vss_stripe = stripe_b
        else:
            vss_stripe = stripe_a
            vdd_stripe = stripe_b

        cell_y = y_top + row_height / 2.0
        x = distance_apart / 2.0

        col_idx = 0
        while x + cell_width / 2.0 <= stripe_a.end:
            instance_name = f"XCELL_R{row_idx}_C{col_idx}"
            cell_x = x + cell_width / 2.0

            pin_connections = {}
            for pin in cell_cfg.instance_pins:
                if pin == "VDD":
                    vdd_node = _get_or_create_node(
                        nodes, m1_name, cell_x, vdd_stripe.position,
                        m1_rl.z_bottom_nm, "VDD",
                    )
                    pin_connections["VDD"] = vdd_node.name
                elif pin == "VSS":
                    vss_node = _get_or_create_node(
                        nodes, m1_name, cell_x, vss_stripe.position,
                        m1_rl.z_bottom_nm, "VSS",
                    )
                    pin_connections["VSS"] = vss_node.name
                else:
                    pin_connections[pin] = f"NC_{nc_counter}"
                    nc_counter += 1

            cells.append(CellPlacement(
                instance_name=instance_name,
                cell_name=cell_cfg.cell,
                x=cell_x, y=cell_y,
                row=row_idx, flipped=flipped,
                pin_connections=pin_connections,
            ))

            x += distance_apart
            col_idx += 1

    return cells


def build_grid(config: Config) -> Grid:
    """Build the complete power grid from configuration."""
    resolved = _resolve_layers(config)
    z_coords = _compute_z_coordinates(config)
    nodes: dict[str, Node] = {}

    grid_size_x = config.grid.grid_size.x
    grid_size_y = config.grid.grid_size.y
    cycle = config.grid.cycle

    stripes = _generate_stripes(resolved, grid_size_x, grid_size_y, cycle)
    staples = _generate_staples(resolved, config, grid_size_x, grid_size_y, cycle, stripes)
    vias = _place_vias(config, resolved, stripes, staples, nodes, z_coords)
    segments = _segment_stripes(resolved, config, stripes, nodes)
    cells = _place_standard_cells(config, resolved, stripes, nodes)

    return Grid(
        stripes=stripes, staples=staples,
        segments=segments, vias=vias,
        cells=cells, nodes=nodes,
    )
