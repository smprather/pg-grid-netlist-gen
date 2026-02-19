"""Core grid generation algorithm."""

from __future__ import annotations

from dataclasses import dataclass
import math
import random

from pg_grid_netlist_gen.config import Config
from pg_grid_netlist_gen.geometry import (
    CellPlacement,
    Grid,
    Node,
    PlocPoint,
    Segment,
    Staple,
    Stripe,
    ViaConnection,
)
from pg_grid_netlist_gen.physics import total_segment_capacitance


@dataclass
class ResolvedLayer:
    """BEOL + grid-layer properties merged together."""

    name: str
    layer_type: str  # "grid" or "staple"
    direction: str  # "vertical" or "horizontal"
    actual_pitch_nm: float
    width_nm: float
    thickness_nm: float
    z_bottom_nm: float
    z_top_nm: float
    resistance_per_square: float


def _compute_z_coordinates(config: Config) -> dict[str, tuple[float, float]]:
    """Walk the BEOL stack bottom-to-top and compute z coordinates."""
    layers_bottom_up = list(reversed(config.beol_stack.layers))
    z_coords: dict[str, tuple[float, float]] = {}
    z = 0.0

    for layer in layers_bottom_up:
        if layer.type == "substrate":
            z_coords[layer.name] = (-layer.thickness, 0.0)
            z = 0.0
            continue

        z_bottom = z
        z_top = z + layer.thickness
        z_coords[layer.name] = (z_bottom, z_top)

        # Vias occupy dielectric height and do not advance the running metal elevation.
        if layer.type != "via":
            z = z_top

    return z_coords


def _resolve_layers(config: Config) -> list[ResolvedLayer]:
    """Resolve configured/implicit metal layers in BEOL top-to-bottom order."""
    z_coords = _compute_z_coordinates(config)
    resolved: list[ResolvedLayer] = []

    metal_layers_in_beol = [l for l in config.beol_stack.layers if l.type == "metal"]

    # Direction alternates bottom-up, with the bottom metal always horizontal.
    directions: dict[str, str] = {}
    current_direction = "horizontal"
    for beol_layer in reversed(metal_layers_in_beol):
        directions[beol_layer.name] = current_direction
        current_direction = "vertical" if current_direction == "horizontal" else "horizontal"

    lowest_layer = config.lowest_metal_layer_name

    for beol_layer in metal_layers_in_beol:
        layer_name = beol_layer.name
        usage = config.grid.layer_usage.get(layer_name)

        if usage is None and layer_name != lowest_layer:
            continue

        if usage is None:
            # Implicit lowest-layer behavior.
            layer_type = "grid"
            width_nm = beol_layer.min_width or 0.0
            actual_pitch_nm = config.distance_to_nm(config.standard_cell_placement.row_height)
        else:
            layer_type = "grid" if usage.type == "g" else "staple"
            # Configured layer_usage width/pitch are multipliers:
            # width = factor * WMIN, pitch = factor * (WMIN + SMIN).
            base_width_nm = beol_layer.min_width or 0.0
            base_pitch_nm = beol_layer.pitch or 0.0
            width_nm = usage.width * base_width_nm
            actual_pitch_nm = usage.pitch * base_pitch_nm

        z_bottom, z_top = z_coords[layer_name]

        resolved.append(
            ResolvedLayer(
                name=layer_name,
                layer_type=layer_type,
                direction=directions[layer_name],
                actual_pitch_nm=actual_pitch_nm,
                width_nm=width_nm,
                thickness_nm=beol_layer.thickness,
                z_bottom_nm=z_bottom,
                z_top_nm=z_top,
                resistance_per_square=(beol_layer.resistance_per_square or 0.0),
            )
        )

    return resolved


def _node_name(layer: str, x: float, y: float) -> str:
    return f"{layer}_X_{int(round(x))}_Y_{int(round(y))}"


def _staple_key(layer: str, x: float, y: float, net: str) -> tuple[str, int, int, str]:
    return (layer, int(round(x)), int(round(y)), net)


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


def _create_named_node(
    nodes: dict[str, Node],
    name: str,
    layer: str,
    x: float,
    y: float,
    z: float,
    net: str,
) -> Node:
    if name not in nodes:
        nodes[name] = Node(name=name, x=x, y=y, z=z, layer=layer, net=net)
    return nodes[name]


def _get_resolved(resolved: list[ResolvedLayer], name: str) -> ResolvedLayer:
    for rl in resolved:
        if rl.name == name:
            return rl
    raise ValueError(f"Resolved layer '{name}' not found")


def _generate_stripes(
    resolved: list[ResolvedLayer],
    grid_size_x: float,
    grid_size_y: float,
    cycle: list[str],
) -> list[Stripe]:
    """Generate stripes for all grid-type layers."""
    stripes: list[Stripe] = []

    for rl in resolved:
        if rl.layer_type != "grid" or rl.actual_pitch_nm <= 0:
            continue

        if rl.direction == "vertical":
            n_stripes = int(grid_size_x / rl.actual_pitch_nm) + 1
            for i in range(n_stripes):
                pos = i * rl.actual_pitch_nm
                if pos > grid_size_x:
                    break
                net = cycle[i % len(cycle)]
                stripes.append(
                    Stripe(
                        layer=rl.name,
                        direction="vertical",
                        position=pos,
                        start=0.0,
                        end=grid_size_y,
                        width=rl.width_nm,
                        net=net,
                    )
                )
        else:
            n_stripes = int(grid_size_y / rl.actual_pitch_nm) + 1
            for i in range(n_stripes):
                pos = i * rl.actual_pitch_nm
                if pos > grid_size_y:
                    break
                net = cycle[i % len(cycle)]
                stripes.append(
                    Stripe(
                        layer=rl.name,
                        direction="horizontal",
                        position=pos,
                        start=0.0,
                        end=grid_size_x,
                        width=rl.width_nm,
                        net=net,
                    )
                )

    return stripes


def _generate_staples(
    resolved: list[ResolvedLayer],
    config: Config,
    stripes: list[Stripe],
) -> list[Staple]:
    """Generate staple pads at intersections of nearest grid layers above and below."""
    staples: list[Staple] = []
    by_layer: dict[str, list[Stripe]] = {}
    for s in stripes:
        by_layer.setdefault(s.layer, []).append(s)

    via_lookup = {(v.from_layer, v.to_layer): v for v in config.itf_stack.vias}
    via_lookup.update({(b, a): v for (a, b), v in list(via_lookup.items())})

    for idx, rl in enumerate(resolved):
        if rl.layer_type != "staple":
            continue

        above = next((x for x in reversed(resolved[:idx]) if x.layer_type == "grid"), None)
        below = next((x for x in resolved[idx + 1 :] if x.layer_type == "grid"), None)
        if above is None or below is None:
            continue

        above_stripes = by_layer.get(above.name, [])
        below_stripes = by_layer.get(below.name, [])
        if not above_stripes or not below_stripes:
            continue

        # Spec rule: staple side = 2x via side length of the via above this staple layer.
        staple_size = 0.0
        if idx > 0:
            upper_adjacent = resolved[idx - 1]
            itf_via = via_lookup.get((upper_adjacent.name, rl.name))
            if itf_via and itf_via.area_nm2 > 0.0:
                staple_size = 2.0 * math.sqrt(itf_via.area_nm2)

        if staple_size <= 0.0:
            # Fallback to preserve behavior if via metadata is missing.
            above_beol = config.get_beol_layer(above.name)
            staple_size = 2.0 * (above_beol.min_width or above.width_nm)

        seen: set[tuple[str, int, int, str]] = set()
        for sa in above_stripes:
            for sb in below_stripes:
                if sa.net != sb.net:
                    continue

                if sa.direction == "vertical" and sb.direction == "horizontal":
                    x, y = sa.position, sb.position
                elif sa.direction == "horizontal" and sb.direction == "vertical":
                    x, y = sb.position, sa.position
                else:
                    continue

                key = _staple_key(rl.name, x, y, sa.net)
                if key in seen:
                    continue
                seen.add(key)
                staples.append(Staple(layer=rl.name, x=x, y=y, size=staple_size, net=sa.net))

    return staples


def _create_staple_segments(
    config: Config,
    resolved: list[ResolvedLayer],
    staples: list[Staple],
    nodes: dict[str, Node],
) -> tuple[list[Segment], dict[tuple[str, int, int, str], tuple[Node, Node]]]:
    """Create one-square staple resistors and return staple-side nodes for via attachment."""
    segments: list[Segment] = []
    staple_nodes: dict[tuple[str, int, int, str], tuple[Node, Node]] = {}

    for staple in staples:
        rl = _get_resolved(resolved, staple.layer)
        beol = config.get_beol_layer(staple.layer)

        key = _staple_key(staple.layer, staple.x, staple.y, staple.net)
        up_name = f"{staple.layer}_X_{int(round(staple.x))}_Y_{int(round(staple.y))}_UP"
        dn_name = f"{staple.layer}_X_{int(round(staple.x))}_Y_{int(round(staple.y))}_DN"
        up_node = _create_named_node(nodes, up_name, staple.layer, staple.x, staple.y, rl.z_top_nm, staple.net)
        dn_node = _create_named_node(nodes, dn_name, staple.layer, staple.x, staple.y, rl.z_bottom_nm, staple.net)
        staple_nodes[key] = (up_node, dn_node)

        r = (beol.resistance_per_square or 0.0) * 1.0
        segments.append(
            Segment(
                node_a=up_node,
                node_b=dn_node,
                layer=staple.layer,
                width=staple.size,
                thickness=rl.thickness_nm,
                length=staple.size,
                net=staple.net,
                resistance=r,
                cap_plate=0.0,
                cap_fringe=0.0,
            )
        )

    return segments, staple_nodes


def _add_via_nodes(
    vias: list[ViaConnection],
    via_beol,
    itf_via,
    node_top: Node,
    node_bot: Node,
    net: str,
) -> None:
    vias.append(
        ViaConnection(
            node_top=node_top,
            node_bot=node_bot,
            via_layer=via_beol.name,
            width=via_beol.min_width,
            height=via_beol.thickness,
            net=net,
            resistance=itf_via.resistance_per_via,
        )
    )


def _place_vias(
    config: Config,
    resolved: list[ResolvedLayer],
    stripes: list[Stripe],
    staples: list[Staple],
    staple_nodes: dict[tuple[str, int, int, str], tuple[Node, Node]],
    nodes: dict[str, Node],
) -> list[ViaConnection]:
    """Place vias between adjacent resolved layers where shapes intersect."""
    vias: list[ViaConnection] = []

    stripes_by_layer: dict[str, list[Stripe]] = {}
    for s in stripes:
        stripes_by_layer.setdefault(s.layer, []).append(s)

    staples_by_layer: dict[str, list[Staple]] = {}
    for s in staples:
        staples_by_layer.setdefault(s.layer, []).append(s)

    via_lookup = {(v.from_layer, v.to_layer): v for v in config.itf_stack.vias}
    via_lookup.update({(b, a): v for (a, b), v in list(via_lookup.items())})

    for i in range(len(resolved) - 1):
        top_rl = resolved[i]
        bot_rl = resolved[i + 1]
        itf_via = via_lookup.get((top_rl.name, bot_rl.name))
        if not itf_via:
            continue

        via_beol = config.get_beol_layer(itf_via.name)

        # Grid-to-grid via placement.
        if top_rl.layer_type == "grid" and bot_rl.layer_type == "grid":
            top_stripes = stripes_by_layer.get(top_rl.name, [])
            bot_stripes = stripes_by_layer.get(bot_rl.name, [])
            for ts in top_stripes:
                for bs in bot_stripes:
                    if ts.net != bs.net or ts.direction == bs.direction:
                        continue
                    if ts.direction == "vertical":
                        x, y = ts.position, bs.position
                    else:
                        x, y = bs.position, ts.position

                    node_top = _get_or_create_node(nodes, top_rl.name, x, y, top_rl.z_bottom_nm, ts.net)
                    node_bot = _get_or_create_node(nodes, bot_rl.name, x, y, bot_rl.z_bottom_nm, ts.net)
                    _add_via_nodes(vias, via_beol, itf_via, node_top, node_bot, ts.net)
            continue

        # Grid-to-staple (top grid, bottom staple).
        if top_rl.layer_type == "grid" and bot_rl.layer_type == "staple":
            top_stripes = stripes_by_layer.get(top_rl.name, [])
            bot_staples = staples_by_layer.get(bot_rl.name, [])
            for ts in top_stripes:
                for st in bot_staples:
                    if ts.net != st.net:
                        continue
                    if ts.direction == "vertical":
                        intersects = abs(st.x - ts.position) <= (st.size + ts.width) / 2.0
                    else:
                        intersects = abs(st.y - ts.position) <= (st.size + ts.width) / 2.0
                    if not intersects:
                        continue

                    key = _staple_key(st.layer, st.x, st.y, st.net)
                    up_node = staple_nodes[key][0]
                    node_top = _get_or_create_node(nodes, top_rl.name, st.x, st.y, top_rl.z_bottom_nm, ts.net)
                    _add_via_nodes(vias, via_beol, itf_via, node_top, up_node, ts.net)
            continue

        # Staple-to-grid (top staple, bottom grid).
        if top_rl.layer_type == "staple" and bot_rl.layer_type == "grid":
            top_staples = staples_by_layer.get(top_rl.name, [])
            bot_stripes = stripes_by_layer.get(bot_rl.name, [])
            for st in top_staples:
                for bs in bot_stripes:
                    if st.net != bs.net:
                        continue
                    if bs.direction == "vertical":
                        intersects = abs(st.x - bs.position) <= (st.size + bs.width) / 2.0
                    else:
                        intersects = abs(st.y - bs.position) <= (st.size + bs.width) / 2.0
                    if not intersects:
                        continue

                    key = _staple_key(st.layer, st.x, st.y, st.net)
                    dn_node = staple_nodes[key][1]
                    node_bot = _get_or_create_node(nodes, bot_rl.name, st.x, st.y, bot_rl.z_bottom_nm, st.net)
                    _add_via_nodes(vias, via_beol, itf_via, dn_node, node_bot, st.net)
            continue

        # Staple-to-staple (adjacent staple layers).
        if top_rl.layer_type == "staple" and bot_rl.layer_type == "staple":
            top_staples = staples_by_layer.get(top_rl.name, [])
            bot_staples = staples_by_layer.get(bot_rl.name, [])
            bot_by_xy_net = {
                (int(round(st.x)), int(round(st.y)), st.net): st
                for st in bot_staples
            }
            for st_top in top_staples:
                key_xy_net = (int(round(st_top.x)), int(round(st_top.y)), st_top.net)
                st_bot = bot_by_xy_net.get(key_xy_net)
                if st_bot is None:
                    continue

                top_key = _staple_key(st_top.layer, st_top.x, st_top.y, st_top.net)
                bot_key = _staple_key(st_bot.layer, st_bot.x, st_bot.y, st_bot.net)
                top_nodes = staple_nodes.get(top_key)
                bot_nodes = staple_nodes.get(bot_key)
                if top_nodes is None or bot_nodes is None:
                    continue

                dn_node = top_nodes[1]
                up_node = bot_nodes[0]
                _add_via_nodes(vias, via_beol, itf_via, dn_node, up_node, st_top.net)

    return vias


def _segment_stripes(
    resolved: list[ResolvedLayer],
    config: Config,
    stripes: list[Stripe],
    nodes: dict[str, Node],
) -> list[Segment]:
    """Break stripes into segments at via/cell/PLOC connection points."""
    segments: list[Segment] = []

    stripe_points: dict[tuple[str, float, str], list[float]] = {}
    for stripe in stripes:
        key = (stripe.layer, stripe.position, stripe.net)
        stripe_points.setdefault(key, [])

    for node in nodes.values():
        for stripe in stripes:
            if stripe.layer != node.layer or stripe.net != node.net:
                continue
            if stripe.direction == "vertical" and abs(stripe.position - node.x) < 1e-6:
                stripe_points[(stripe.layer, stripe.position, stripe.net)].append(node.y)
            elif stripe.direction == "horizontal" and abs(stripe.position - node.y) < 1e-6:
                stripe_points[(stripe.layer, stripe.position, stripe.net)].append(node.x)

    oxide_mat = config.get_material("oxide")

    for stripe in stripes:
        rl = _get_resolved(resolved, stripe.layer)
        key = (stripe.layer, stripe.position, stripe.net)
        points = sorted(set([stripe.start, *stripe_points.get(key, []), stripe.end]))

        for j in range(len(points) - 1):
            c_start, c_end = points[j], points[j + 1]
            length = c_end - c_start
            if length <= 1e-9:
                continue

            if stripe.direction == "vertical":
                x = stripe.position
                node_a = _get_or_create_node(nodes, stripe.layer, x, c_start, rl.z_bottom_nm, stripe.net)
                node_b = _get_or_create_node(nodes, stripe.layer, x, c_end, rl.z_bottom_nm, stripe.net)
            else:
                y = stripe.position
                node_a = _get_or_create_node(nodes, stripe.layer, c_start, y, rl.z_bottom_nm, stripe.net)
                node_b = _get_or_create_node(nodes, stripe.layer, c_end, y, rl.z_bottom_nm, stripe.net)

            if rl.width_nm > 0 and rl.resistance_per_square > 0:
                r = rl.resistance_per_square * (length / rl.width_nm)
            else:
                r = 0.0
            cp, cf, _ = total_segment_capacitance(
                rl.width_nm,
                length,
                rl.thickness_nm,
                rl.z_bottom_nm,
                oxide_mat.get("relative_permittivity", 2.5),
            )
            segments.append(
                Segment(
                    node_a=node_a,
                    node_b=node_b,
                    layer=stripe.layer,
                    width=rl.width_nm,
                    thickness=rl.thickness_nm,
                    length=length,
                    net=stripe.net,
                    resistance=r,
                    cap_plate=cp,
                    cap_fringe=cf,
                )
            )
    return segments


def _chain_standard_cells(cells: list[CellPlacement], config: Config, rng: random.Random) -> None:
    """Modify pin connections to chain cells together."""
    if not cells:
        return

    rng.shuffle(cells)

    cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
    input_pins = [p.name for p in cell_cfg.pins if p.direction == "input" and p.type == "signal"]
    output_pins = [p.name for p in cell_cfg.pins if p.direction == "output" and p.type == "signal"]
    if not input_pins or not output_pins:
        return

    input_pin = input_pins[0]
    output_pin = output_pins[0]

    max_chain_len = config.spice_netlist.cell_chains.max_instance_count_per_chain

    chain_num = 0
    for i in range(0, len(cells), max_chain_len):
        chain = cells[i : i + max_chain_len]
        chain[0].pin_connections[input_pin] = f"CHAIN_{chain_num}_IN"

        for j in range(len(chain) - 1):
            cell_a = chain[j]
            cell_b = chain[j + 1]
            connection_net = f"chain_{chain_num}_link_{j}"
            cell_a.pin_connections[output_pin] = connection_net
            cell_b.pin_connections[input_pin] = connection_net

        chain[-1].pin_connections[output_pin] = f"CHAIN_{chain_num}_OUT"
        chain_num += 1


def _place_standard_cells(
    config: Config,
    resolved: list[ResolvedLayer],
    stripes: list[Stripe],
    nodes: dict[str, Node],
    grid_size_x: float,
    rng: random.Random,
) -> list[CellPlacement]:
    """Place standard cells on the lowest available grid layer."""
    cells: list[CellPlacement] = []

    available_grid_layers = [rl for rl in resolved if rl.layer_type == "grid"]
    if not available_grid_layers:
        return cells
    placement_layer = available_grid_layers[-1].name
    if not any(rl.name == placement_layer for rl in resolved):
        return cells

    layer_rl = _get_resolved(resolved, placement_layer)
    layer_stripes = sorted([s for s in stripes if s.layer == placement_layer], key=lambda s: s.position)
    if len(layer_stripes) < 2:
        return cells

    power_net_name = config.pg_nets.power.name
    ground_net_name = config.pg_nets.ground.name

    cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
    cell_width_nm = config.distance_to_nm(cell_cfg.width)
    half_w = cell_width_nm / 2.0

    site_w_nm = config.distance_to_nm(config.standard_cell_placement.site_width)
    min_cc_x_nm = config.distance_to_nm(config.standard_cell_placement.min_space["x"])
    step_nm = max(min_cc_x_nm, cell_width_nm)

    stagger_cfg = config.standard_cell_placement.stagger_row_start
    stagger_range_nm = config.distance_to_nm(float(stagger_cfg.get("range", 0.0)))
    stagger_random = bool(stagger_cfg.get("random", False))

    nc_counter = 0
    logical_row = 0

    for row_idx in range(len(layer_stripes) - 1):
        stripe_a = layer_stripes[row_idx]
        stripe_b = layer_stripes[row_idx + 1]
        if {stripe_a.net, stripe_b.net} != {power_net_name, ground_net_name}:
            continue

        flipped = logical_row % 2 == 1
        if flipped:
            vdd_stripe = stripe_a if stripe_a.net == power_net_name else stripe_b
            vss_stripe = stripe_b if vdd_stripe is stripe_a else stripe_a
        else:
            vdd_stripe = stripe_b if stripe_b.net == power_net_name else stripe_a
            vss_stripe = stripe_a if vdd_stripe is stripe_b else stripe_b

        cell_y_nm = (vdd_stripe.position + vss_stripe.position) / 2.0

        if stagger_range_nm > 0:
            if stagger_random:
                row_offset_nm = rng.uniform(0.0, stagger_range_nm)
            else:
                row_offset_nm = stagger_range_nm if logical_row % 2 else 0.0
        else:
            row_offset_nm = 0.0

        left_edge_nm = math.floor(row_offset_nm / site_w_nm) * site_w_nm
        col_idx = 0

        while left_edge_nm + cell_width_nm <= grid_size_x + 1e-6:
            cell_x_nm = left_edge_nm + half_w
            instance_name = f"XCELL_R{logical_row}_C{col_idx}"

            pin_connections: dict[str, str] = {}
            for pin in cell_cfg.pins:
                if pin.type in {"power", "ground"}:
                    net_name = power_net_name if pin.type == "power" else ground_net_name
                    stripe = vdd_stripe if pin.type == "power" else vss_stripe
                    node = _get_or_create_node(
                        nodes,
                        placement_layer,
                        cell_x_nm,
                        stripe.position,
                        layer_rl.z_bottom_nm,
                        net_name,
                    )
                    pin_connections[pin.name] = node.name
                else:
                    pin_connections[pin.name] = f"NC_{pin.name}_{nc_counter}"
                    nc_counter += 1

            cells.append(
                CellPlacement(
                    instance_name=instance_name,
                    cell_name=cell_cfg.name,
                    x=cell_x_nm,
                    y=cell_y_nm,
                    row=logical_row,
                    flipped=flipped,
                    pin_connections=pin_connections,
                )
            )

            next_left = left_edge_nm + step_nm
            left_edge_nm = math.ceil(next_left / site_w_nm) * site_w_nm
            col_idx += 1

        logical_row += 1

    _chain_standard_cells(cells, config, rng)
    return cells


def _generate_plocs(
    config: Config,
    resolved: list[ResolvedLayer],
    stripes: list[Stripe],
    nodes: dict[str, Node],
    grid_size_x: float,
    grid_size_y: float,
) -> list[PlocPoint]:
    """Generate top-layer PLOC points and snap to valid PG stripes."""
    top_grid = next((rl for rl in resolved if rl.layer_type == "grid"), None)
    if top_grid is None:
        return []

    top_stripes = [s for s in stripes if s.layer == top_grid.name]
    if not top_stripes:
        return []

    power_net = config.pg_nets.power.name
    ground_net = config.pg_nets.ground.name

    pitch_nm = config.distance_to_nm(config.ploc.pitch)
    x0 = config.distance_to_nm(config.ploc.offset_from_origin["x"])
    y0 = config.distance_to_nm(config.ploc.offset_from_origin["y"])
    if pitch_nm <= 0:
        return []

    plocs: list[PlocPoint] = []
    seen_nodes: set[str] = set()

    row = 0
    y = y0
    while y <= grid_size_y + 1e-6:
        x = x0 + (pitch_nm / 2.0 if config.ploc.staggered and row % 2 else 0.0)
        col = 0
        while x <= grid_size_x + 1e-6:
            desired_net = power_net if col % 2 == 0 else ground_net
            candidates = [s for s in top_stripes if s.net == desired_net]
            if candidates:
                if top_grid.direction == "vertical":
                    snap = min(candidates, key=lambda s: abs(s.position - x))
                    sx, sy = snap.position, y
                else:
                    snap = min(candidates, key=lambda s: abs(s.position - y))
                    sx, sy = x, snap.position

                node = _get_or_create_node(nodes, top_grid.name, sx, sy, top_grid.z_bottom_nm, desired_net)
                if node.name not in seen_nodes:
                    seen_nodes.add(node.name)
                    plocs.append(PlocPoint(node_name=node.name, x=sx, y=sy, layer=top_grid.name, net=desired_net))

            x += pitch_nm
            col += 1

        y += pitch_nm
        row += 1

    return plocs


def _place_dcap_cells(
    config: Config,
    resolved: list[ResolvedLayer],
    stripes: list[Stripe],
    nodes: dict[str, Node],
    grid_size_x: float,
    grid_size_y: float,
    chain_cells: list[CellPlacement],
    rng: random.Random,
) -> list[CellPlacement]:
    """Place decoupling capacitance cells in available sites."""
    dcap_cfg = config.standard_cell_placement.dcap_cells
    if dcap_cfg is None or not dcap_cfg.enabled:
        return []

    dcap_cell_cfg = config.get_cell_by_name(dcap_cfg.cell)
    dcap_w_nm = config.distance_to_nm(dcap_cell_cfg.width)
    dcap_h_nm = config.distance_to_nm(dcap_cell_cfg.height)

    # Compute max dcap count from density
    grid_area = grid_size_x * grid_size_y
    dcap_area = dcap_w_nm * dcap_h_nm
    if dcap_area <= 0:
        return []
    max_count = int(math.floor(dcap_cfg.max_density_pct / 100.0 * grid_area / dcap_area))
    if max_count < 1:
        return []

    # Find the placement layer (lowest grid layer)
    available_grid_layers = [rl for rl in resolved if rl.layer_type == "grid"]
    if not available_grid_layers:
        return []
    placement_layer = available_grid_layers[-1].name
    layer_rl = _get_resolved(resolved, placement_layer)
    layer_stripes = sorted(
        [s for s in stripes if s.layer == placement_layer], key=lambda s: s.position
    )
    if len(layer_stripes) < 2:
        return []

    power_net_name = config.pg_nets.power.name
    ground_net_name = config.pg_nets.ground.name
    site_w_nm = config.distance_to_nm(config.standard_cell_placement.site_width)

    chain_cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
    chain_w_nm = config.distance_to_nm(chain_cell_cfg.width)

    # Build occupied intervals per row from chain cells
    occupied_by_row: dict[int, list[tuple[float, float]]] = {}
    for cell in chain_cells:
        half_w = chain_w_nm / 2.0
        left = cell.x - half_w
        right = cell.x + half_w
        occupied_by_row.setdefault(cell.row, []).append((left, right))
    for intervals in occupied_by_row.values():
        intervals.sort()

    def _overlaps(left: float, right: float, intervals: list[tuple[float, float]]) -> bool:
        for il, ir in intervals:
            if left < ir and right > il:
                return True
        return False

    # Collect all legal dcap sites
    legal_sites: list[tuple[int, float, float, bool, Stripe, Stripe]] = []
    logical_row = 0
    for row_idx in range(len(layer_stripes) - 1):
        stripe_a = layer_stripes[row_idx]
        stripe_b = layer_stripes[row_idx + 1]
        if {stripe_a.net, stripe_b.net} != {power_net_name, ground_net_name}:
            continue

        flipped = logical_row % 2 == 1
        if flipped:
            vdd_stripe = stripe_a if stripe_a.net == power_net_name else stripe_b
            vss_stripe = stripe_b if vdd_stripe is stripe_a else stripe_a
        else:
            vdd_stripe = stripe_b if stripe_b.net == power_net_name else stripe_a
            vss_stripe = stripe_a if vdd_stripe is stripe_b else stripe_b

        row_intervals = occupied_by_row.get(logical_row, [])

        left_edge = 0.0
        while left_edge + dcap_w_nm <= grid_size_x + 1e-6:
            right_edge = left_edge + dcap_w_nm
            if not _overlaps(left_edge, right_edge, row_intervals):
                legal_sites.append(
                    (logical_row, left_edge, right_edge, flipped, vdd_stripe, vss_stripe)
                )
            left_edge = math.ceil((left_edge + site_w_nm) / site_w_nm) * site_w_nm

        logical_row += 1

    if not legal_sites:
        return []

    # Greedy random selection: shuffle, then pick non-overlapping sites.
    # Track occupied intervals per row to prevent dcap-to-dcap overlap.
    rng.shuffle(legal_sites)
    placed_by_row: dict[int, list[tuple[float, float]]] = {}
    chosen: list[tuple[int, float, float, bool, Stripe, Stripe]] = []
    for site in legal_sites:
        if len(chosen) >= max_count:
            break
        row, left, right = site[0], site[1], site[2]
        row_placed = placed_by_row.get(row, [])
        if _overlaps(left, right, row_placed):
            continue
        chosen.append(site)
        placed_by_row.setdefault(row, []).append((left, right))

    dcap_cells: list[CellPlacement] = []
    for idx, (row, left, right, flipped, vdd_stripe, vss_stripe) in enumerate(chosen):
        cell_x_nm = (left + right) / 2.0
        cell_y_nm = (vdd_stripe.position + vss_stripe.position) / 2.0
        instance_name = f"XDCAP_{idx}"

        pin_connections: dict[str, str] = {}
        for pin in dcap_cell_cfg.pins:
            if pin.type == "power":
                node = _get_or_create_node(
                    nodes, placement_layer, cell_x_nm, vdd_stripe.position,
                    layer_rl.z_bottom_nm, power_net_name,
                )
                pin_connections[pin.name] = node.name
            elif pin.type == "ground":
                node = _get_or_create_node(
                    nodes, placement_layer, cell_x_nm, vss_stripe.position,
                    layer_rl.z_bottom_nm, ground_net_name,
                )
                pin_connections[pin.name] = node.name

        dcap_cells.append(
            CellPlacement(
                instance_name=instance_name,
                cell_name=dcap_cell_cfg.name,
                x=cell_x_nm,
                y=cell_y_nm,
                row=row,
                flipped=flipped,
                pin_connections=pin_connections,
            )
        )

    return dcap_cells


def build_grid(config: Config) -> Grid:
    """Build the complete power grid from configuration."""
    rng = config.make_rng()
    resolved = _resolve_layers(config)
    nodes: dict[str, Node] = {}

    grid_size_y = config.grid.size["rows"] * config.distance_to_nm(config.standard_cell_placement.row_height)
    grid_size_x = config.grid.size["sites"] * config.distance_to_nm(config.standard_cell_placement.site_width)

    nets = [config.pg_nets.power.name, config.pg_nets.ground.name]

    stripes = _generate_stripes(resolved, grid_size_x, grid_size_y, nets)
    staples = _generate_staples(resolved, config, stripes)
    staple_segments, staple_nodes = _create_staple_segments(config, resolved, staples, nodes)
    vias = _place_vias(config, resolved, stripes, staples, staple_nodes, nodes)
    cells = _place_standard_cells(config, resolved, stripes, nodes, grid_size_x, rng)
    dcap_cells = _place_dcap_cells(
        config, resolved, stripes, nodes, grid_size_x, grid_size_y, cells, rng,
    )
    plocs = _generate_plocs(config, resolved, stripes, nodes, grid_size_x, grid_size_y)
    segments = _segment_stripes(resolved, config, stripes, nodes)
    segments.extend(staple_segments)

    power_net = config.pg_nets.power.name
    ground_net = config.pg_nets.ground.name
    power_plocs = sum(1 for p in plocs if p.net == power_net)
    ground_plocs = sum(1 for p in plocs if p.net == ground_net)
    if power_plocs < 1 or ground_plocs < 1:
        raise ValueError(
            "PLOC generation failed: expected at least one power and one ground PLOC "
            f"(power={power_plocs}, ground={ground_plocs}). "
            "Adjust ploc.pitch and/or ploc.offset_from_origin to place valid points in-grid."
        )

    return Grid(
        stripes=stripes,
        staples=staples,
        segments=segments,
        vias=vias,
        cells=cells,
        dcap_cells=dcap_cells,
        plocs=plocs,
        nodes=nodes,
    )
