"""Plotly 2D visualization of the power grid."""

from __future__ import annotations

import math
from pathlib import Path

import plotly.graph_objects as go
import plotly.offline as pyo

from pg_grid_netlist_gen.config import Config
from pg_grid_netlist_gen.geometry import CellPlacement, Grid

# Layer color map
LAYER_COLORS: dict[str, str] = {
    "M0": "rgba(40, 120, 220, 0.70)",
    "M1": "rgba(40, 120, 220, 0.70)",
    "M2": "rgba(40, 120, 220, 0.70)",
    "M3": "rgba(40, 120, 220, 0.70)",
    "M4": "rgba(40, 170, 90, 0.70)",
    "M5": "rgba(40, 170, 90, 0.70)",
    "M6": "rgba(40, 170, 90, 0.70)",
    "M7": "rgba(220, 150, 50, 0.70)",
    "M8": "rgba(220, 150, 50, 0.70)",
    "M9": "rgba(220, 150, 50, 0.70)",
    "M10": "rgba(245, 175, 80, 0.70)",
    "M11": "rgba(245, 175, 80, 0.70)",
    "M12": "rgba(210, 70, 70, 0.70)",
    "M13": "rgba(210, 70, 70, 0.70)",
    "RDL": "rgba(190, 50, 50, 0.70)",
}
VIA_COLOR = "rgba(150, 150, 150, 0.85)"
CELL_COLOR = "rgba(220, 200, 60, 0.55)"
DCAP_COLOR = "rgba(100, 200, 180, 0.55)"
FLIGHT_COLOR = "rgba(120, 30, 30, 0.85)"
OXIDE_COLOR = "rgba(236, 236, 210, 0.75)"
PLOC_DEFAULT_COLOR = "rgba(30, 30, 30, 0.80)"


def _get_color(layer_name: str) -> str:
    if layer_name.startswith("V"):
        return VIA_COLOR
    return LAYER_COLORS.get(layer_name, "rgba(130, 130, 130, 0.7)")


def _in_region_segment(seg, region):
    x_min, y_min, x_max, y_max = region
    sx_min, sx_max = min(seg.node_a.x, seg.node_b.x), max(seg.node_a.x, seg.node_b.x)
    sy_min, sy_max = min(seg.node_a.y, seg.node_b.y), max(seg.node_a.y, seg.node_b.y)
    return sx_max >= x_min and sx_min <= x_max and sy_max >= y_min and sy_min <= y_max


def _in_region_point(x, y, region):
    return region[0] <= x <= region[2] and region[1] <= y <= region[3]


def _stripe_intersections_for_layers(
    grid: Grid, layer_a: str, layer_b: str
) -> list[tuple[float, float, str]]:
    """Build via display points from intersections of two stripe layers by net."""
    stripes_a = [s for s in grid.stripes if s.layer == layer_a]
    stripes_b = [s for s in grid.stripes if s.layer == layer_b]
    if not stripes_a or not stripes_b:
        return []

    points: list[tuple[float, float, str]] = []
    seen: set[tuple[int, int, str]] = set()
    for sa in stripes_a:
        for sb in stripes_b:
            if sa.net != sb.net:
                continue
            if sa.direction == "horizontal" and sb.direction == "vertical":
                x, y = sb.position, sa.position
            elif sa.direction == "vertical" and sb.direction == "horizontal":
                x, y = sa.position, sb.position
            else:
                continue

            key = (int(round(x)), int(round(y)), sa.net)
            if key in seen:
                continue
            seen.add(key)
            points.append((x, y, sa.net))

    return points


def _via_points_from_itf_connection(
    grid: Grid, config: Config, via_layer: str
) -> list[tuple[float, float, str]]:
    """Derive via display points from ITF FROM/TO metals for this via layer."""
    itf_via = next((v for v in config.itf_stack.vias if v.name == via_layer), None)
    if itf_via is None:
        return []

    return _stripe_intersections_for_layers(grid, itf_via.from_layer, itf_via.to_layer)


def _is_grid_layer(config: Config, layer_name: str) -> bool:
    """Check if a layer is a grid-type layer (not staple)."""
    usage = config.grid.layer_usage.get(layer_name)
    if usage is not None:
        return usage.type == "grid"
    # The implicit lowest layer is always grid-type.
    return layer_name == config.lowest_metal_layer_name


def _multi_via_rects(
    x_center: float,
    y_center: float,
    top_stripe_w: float,
    bot_stripe_w: float,
    via_side: float,
    min_space_factor: float,
) -> list[tuple[float, float, float, float]]:
    """Return list of (x0, y0, x1, y1) for each via in a 2D array at a crossing."""
    rect_w = top_stripe_w
    rect_h = bot_stripe_w
    via_pitch = via_side * min_space_factor
    if via_pitch <= 0:
        half = via_side / 2.0
        return [(x_center - half, y_center - half, x_center + half, y_center + half)]

    n_x = max(1, int(rect_w / via_pitch))
    n_y = max(1, int(rect_h / via_pitch))

    # Center the array within the overlap rectangle.
    array_w = (n_x - 1) * via_pitch if n_x > 1 else 0.0
    array_h = (n_y - 1) * via_pitch if n_y > 1 else 0.0
    start_x = x_center - array_w / 2.0
    start_y = y_center - array_h / 2.0

    rects: list[tuple[float, float, float, float]] = []
    half = via_side / 2.0
    for ix in range(n_x):
        for iy in range(n_y):
            cx = start_x + ix * via_pitch
            cy = start_y + iy * via_pitch
            rects.append((cx - half, cy - half, cx + half, cy + half))
    return rects


def _stripe_intersections_with_widths(
    grid: Grid, config: Config, layer_a: str, layer_b: str
) -> list[tuple[float, float, str, float, float]]:
    """Like _stripe_intersections_for_layers but also returns stripe widths.

    Returns list of (x, y, net, stripe_a_width, stripe_b_width).
    """
    stripes_a = [s for s in grid.stripes if s.layer == layer_a]
    stripes_b = [s for s in grid.stripes if s.layer == layer_b]
    if not stripes_a or not stripes_b:
        return []

    points: list[tuple[float, float, str, float, float]] = []
    seen: set[tuple[int, int, str]] = set()
    for sa in stripes_a:
        for sb in stripes_b:
            if sa.net != sb.net:
                continue
            if sa.direction == "horizontal" and sb.direction == "vertical":
                x, y = sb.position, sa.position
            elif sa.direction == "vertical" and sb.direction == "horizontal":
                x, y = sa.position, sb.position
            else:
                continue

            key = (int(round(x)), int(round(y)), sa.net)
            if key in seen:
                continue
            seen.add(key)
            points.append((x, y, sa.net, sa.width, sb.width))

    return points


def _pin_side_for_cell(cell: CellPlacement, pin_location: str) -> str:
    """Return effective pin side after row flipping."""
    if not cell.flipped:
        return pin_location
    if pin_location == "top":
        return "bottom"
    if pin_location == "bottom":
        return "top"
    return pin_location


def _pin_xy(
    cell: CellPlacement, pin_location: str, cell_w_nm: float, cell_h_nm: float
) -> tuple[float, float]:
    side = _pin_side_for_cell(cell, pin_location)
    if side == "left":
        return cell.x - (cell_w_nm / 2.0), cell.y
    if side == "right":
        return cell.x + (cell_w_nm / 2.0), cell.y
    if side == "top":
        return cell.x, cell.y + (cell_h_nm / 2.0)
    return cell.x, cell.y - (cell_h_nm / 2.0)


def _draw_order_from_itf(config: Config) -> list[str]:
    """Return legend/draw order as top metal, via to next metal, next metal, ..."""
    metals_top_to_bottom = [c.name for c in config.itf_stack.conductors]
    if not metals_top_to_bottom:
        return []

    via_between: dict[frozenset[str], str] = {}
    conductor_names = set(metals_top_to_bottom)
    for via in config.itf_stack.vias:
        if via.from_layer not in conductor_names or via.to_layer not in conductor_names:
            continue
        via_between[frozenset((via.from_layer, via.to_layer))] = via.name

    order: list[str] = [metals_top_to_bottom[0]]
    for i in range(len(metals_top_to_bottom) - 1):
        upper = metals_top_to_bottom[i]
        lower = metals_top_to_bottom[i + 1]
        via_name = via_between.get(frozenset((upper, lower)))
        if via_name:
            order.append(via_name)
        order.append(lower)
    return order


def _ploc_color(config: Config, net_name: str) -> str:
    if net_name == config.pg_nets.power.name:
        return "rgba(190, 50, 50, 0.85)"
    if net_name == config.pg_nets.ground.name:
        return "rgba(40, 120, 220, 0.85)"
    return PLOC_DEFAULT_COLOR


def _build_cross_section(config: Config) -> go.Figure:
    """Build a standalone cross-section figure of the ITF stack.

    Requirements from AGENTS.md:
    - Render FEOL as a layer; bottom of FEOL = y=0
    - Show substrate as a layer 3x the FEOL thickness, below FEOL
    - Add an in-plot "Substrate" label on the substrate layer
    - VIA shapes labeled with just VIA<N>
    - No legend panel — hover-data is sufficient
    - Do not render dielectric layers (skip <via_layer>_diel)
    - Layers with the same thickness share the same color
    """
    fig = go.Figure()

    conductors_bottom_up = list(reversed(config.itf_stack.conductors))
    # Top-to-bottom order for legend ordering.
    conductors_top_to_bottom = list(config.itf_stack.conductors)

    feol_thickness_nm = config.distance_to_nm(config.feol_thickness)

    # Collect all layer thicknesses to build a thickness-to-color map.
    all_thicknesses: list[float] = []
    for metal in conductors_bottom_up:
        all_thicknesses.append(metal.thickness_nm)
    all_thicknesses.append(feol_thickness_nm)  # FEOL and substrate share this thickness

    unique_thicknesses = sorted(set(round(t, 6) for t in all_thicknesses if t > 0))
    # Assign a distinct color per unique thickness.
    _palette = [
        "rgba(40, 120, 220, 0.70)",
        "rgba(40, 170, 90, 0.70)",
        "rgba(220, 150, 50, 0.70)",
        "rgba(210, 70, 70, 0.70)",
        "rgba(130, 80, 190, 0.70)",
        "rgba(245, 175, 80, 0.70)",
        "rgba(80, 190, 190, 0.70)",
        "rgba(190, 50, 50, 0.70)",
        "rgba(100, 200, 100, 0.70)",
        "rgba(200, 100, 200, 0.70)",
        "rgba(160, 160, 60, 0.70)",
        "rgba(60, 160, 160, 0.70)",
    ]
    thickness_color: dict[float, str] = {}
    for i, t in enumerate(unique_thicknesses):
        thickness_color[t] = _palette[i % len(_palette)]

    def _color_for_thickness(t: float) -> str:
        return thickness_color.get(round(t, 6), "rgba(130, 130, 130, 0.70)")

    def _add_rect(
        x0: float,
        x1: float,
        z_bottom: float,
        z_top: float,
        color: str,
        legend_name: str,
    ) -> None:
        x = [x0, x1, x1, x0, x0, None]
        y = [
            z_bottom / 1000.0,
            z_bottom / 1000.0,
            z_top / 1000.0,
            z_top / 1000.0,
            z_bottom / 1000.0,
            None,
        ]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                fill="toself",
                fillcolor=color,
                line=dict(color=color, width=0.5),
                name=legend_name,
                legendgroup=legend_name,
                showlegend=False,
                hoverinfo="text",
                hovertext=legend_name,
                hovertemplate="%{hovertext}<extra></extra>",
            ),
        )

    # --- First pass: compute all z-coordinates bottom-up ---
    # y=0 = bottom of FEOL.  Substrate sits below at (-3*feol_thickness, 0).
    # FEOL sits at (0, feol_thickness). Metals start at feol_thickness.
    z = feol_thickness_nm  # bottom of lowest metal
    dielectrics = {d.name: d for d in config.itf_stack.dielectrics}
    metal_z: dict[str, tuple[float, float]] = {}
    for metal in conductors_bottom_up:
        m_bottom = z
        m_top = m_bottom + metal.thickness_nm
        metal_z[metal.name] = (m_bottom, m_top)
        z = m_top

        diel_name = f"{metal.name}_diel"
        dielectric = dielectrics.get(diel_name)
        if dielectric is not None and dielectric.thickness_nm > 0:
            z += dielectric.thickness_nm

    # Compute via z-coordinates.
    via_entries: list[tuple[object, float, float, float, int]] = []
    for idx, via in enumerate(config.itf_stack.vias):
        from_layer = via.from_layer
        to_layer = via.to_layer
        if from_layer not in metal_z or to_layer not in metal_z:
            continue

        from_z = metal_z[from_layer]
        to_z = metal_z[to_layer]
        if from_z[0] <= to_z[0]:
            lower_name, upper_name = from_layer, to_layer
        else:
            lower_name, upper_name = to_layer, from_layer

        lower_top = metal_z[lower_name][1]
        upper_bottom = metal_z[upper_name][0]
        via_bottom = lower_top
        via_top = upper_bottom

        via_side = max(via.area_nm2, 0.0) ** 0.5
        if via_top <= via_bottom:
            via_top = via_bottom + max(via_side, 1.0)

        via_entries.append((via, via_bottom, via_top, via_side, idx))

    # Normalize via widths for rendering.
    ref_side = 0.0
    if via_entries:
        top_entry = max(via_entries, key=lambda t: t[1])
        ref_side = top_entry[3]
    if ref_side <= 0.0:
        ref_side = max((entry[3] for entry in via_entries), default=1.0)
    if ref_side <= 0.0:
        ref_side = 1.0

    # --- Second pass: add traces in top-to-bottom order for correct legend ---
    # Build a via lookup: keyed by (from_layer, to_layer) or (to_layer, from_layer).
    via_by_upper_lower: dict[
        tuple[str, str], tuple[object, float, float, float, int]
    ] = {}
    for entry in via_entries:
        via_obj, vb, vt, vs, vi = entry
        from_layer = via_obj.from_layer
        to_layer = via_obj.to_layer
        from_z_val = metal_z.get(from_layer, (0, 0))
        to_z_val = metal_z.get(to_layer, (0, 0))
        if from_z_val[0] >= to_z_val[0]:
            upper, lower = from_layer, to_layer
        else:
            upper, lower = to_layer, from_layer
        via_by_upper_lower[(upper, lower)] = entry

    for i, metal in enumerate(conductors_top_to_bottom):
        # Add metal layer.
        m_bottom, m_top = metal_z[metal.name]
        _add_rect(
            0.0,
            1.0,
            m_bottom,
            m_top,
            _color_for_thickness(metal.thickness_nm),
            metal.name,
        )

        # If there's a via below this metal (connecting to next metal down), add it.
        if i < len(conductors_top_to_bottom) - 1:
            next_metal = conductors_top_to_bottom[i + 1]
            entry = via_by_upper_lower.get((metal.name, next_metal.name))
            if entry is not None:
                via_obj, via_bottom, via_top, via_side, via_idx = entry
                width_norm = via_side / ref_side
                half_w = max(0.001, min(0.45, 0.06 * width_norm))
                via_x0, via_x1 = 0.5 - half_w, 0.5 + half_w
                via_label = f"VIA{via_idx}"
                _add_rect(via_x0, via_x1, via_bottom, via_top, VIA_COLOR, via_label)

    # FEOL layer: from y=0 to y=feol_thickness.
    feol_color = "rgba(180, 160, 120, 0.70)"
    _add_rect(0.0, 1.0, 0.0, feol_thickness_nm, feol_color, "FEOL")

    # Substrate layer: 3x FEOL thickness, below FEOL.
    substrate_color = "rgba(140, 140, 160, 0.70)"
    substrate_bottom = -3 * feol_thickness_nm
    _add_rect(0.0, 1.0, substrate_bottom, 0.0, substrate_color, "Substrate")

    # In-plot label for substrate.
    substrate_mid_y = (substrate_bottom / 1000.0 + 0.0) / 2.0
    fig.add_annotation(
        x=0.5,
        y=substrate_mid_y,
        text="Substrate",
        showarrow=False,
        font=dict(size=14, color="white"),
        xanchor="center",
        yanchor="middle",
    )

    fig.update_layout(
        title="ITF Layer Cross-Section",
        xaxis=dict(showticklabels=False, title="", range=[0, 1]),
        yaxis=dict(title="Height (um)"),
        width=1050,
        height=800,
    )

    return fig


def render_grid(
    grid: Grid,
    config: Config,
    output_path: str | Path | None,
    viz_region: tuple[float, float, float, float] | None = None,
    open_browser: bool = False,
    save_image_layer: str | None = None,
    output_dir: Path | None = None,
) -> None:
    """Render combined 2D grid and layer-stack cross-section to HTML.

    Uses two independent Plotly figures combined into a single HTML file
    via plotly.offline.plot() with include_plotlyjs=True for the first
    div and include_plotlyjs=False for the second.
    """
    if not (output_path or save_image_layer):
        return

    output_dir = output_dir or (
        Path(output_path).parent if output_path else Path("output")
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    region_nm = tuple(c * 1000 for c in viz_region) if viz_region else None

    fig = go.Figure()

    generated_metal_layers = {seg.layer for seg in grid.segments}
    bottom_to_top = [c.name for c in reversed(config.itf_stack.conductors)]
    bottom_four_visible = [
        name for name in bottom_to_top if name in generated_metal_layers
    ][:4]
    bottom_four_set = set(bottom_four_visible)
    visible_via_layers = {
        v.name
        for v in config.itf_stack.vias
        if v.from_layer in bottom_four_set or v.to_layer in bottom_four_set
    }
    configured_visible_objects: set[str] | None = None
    if config.visualizer and config.visualizer.initial_visible_objects is not None:
        configured_visible_objects = {
            str(name).strip()
            for name in config.visualizer.initial_visible_objects
            if str(name).strip()
        }

    def _legend_visibility(name: str, default_visible: bool) -> bool | str:
        if configured_visible_objects is None:
            return True if default_visible else "legendonly"
        return True if name in configured_visible_objects else "legendonly"

    # Draw metal and via traces in ITF top-down connectivity order.
    layer_draw_order = _draw_order_from_itf(config)
    beol_by_name = {layer.name: layer for layer in config.beol_stack.layers}
    for layer_name in layer_draw_order:
        beol_layer = beol_by_name.get(layer_name)
        if beol_layer is None:
            continue
        color = _get_color(layer_name)

        if beol_layer.type == "metal":
            by_net: dict[str, tuple[list[float | None], list[float | None]]] = {}
            for seg in grid.segments:
                if seg.layer != layer_name:
                    continue
                if region_nm and not _in_region_segment(seg, region_nm):
                    continue

                half_w = seg.width / 2.0
                dx = abs(seg.node_a.x - seg.node_b.x)
                dy = abs(seg.node_a.y - seg.node_b.y)
                if dx < 1e-6 and dy < 1e-6:
                    # Staple through-segment: same XY node, different Z in netlist model.
                    # Render as a square in 2D.
                    x0, x1 = seg.node_a.x - half_w, seg.node_a.x + half_w
                    y0, y1 = seg.node_a.y - half_w, seg.node_a.y + half_w
                elif dx < 1e-6:
                    x0, x1 = seg.node_a.x - half_w, seg.node_a.x + half_w
                    y0, y1 = (
                        min(seg.node_a.y, seg.node_b.y),
                        max(seg.node_a.y, seg.node_b.y),
                    )
                elif dy < 1e-6:
                    x0, x1 = (
                        min(seg.node_a.x, seg.node_b.x),
                        max(seg.node_a.x, seg.node_b.x),
                    )
                    y0, y1 = seg.node_a.y - half_w, seg.node_a.y + half_w
                else:
                    # Non-orthogonal segment (e.g., staple through-segment). Draw as small square.
                    half = max(seg.width / 2.0, 1.0)
                    x0, x1 = seg.node_a.x - half, seg.node_a.x + half
                    y0, y1 = seg.node_a.y - half, seg.node_a.y + half

                xs, ys = by_net.setdefault(seg.net, ([], []))
                xs.extend(
                    [
                        x0 / 1000.0,
                        x1 / 1000.0,
                        x1 / 1000.0,
                        x0 / 1000.0,
                        x0 / 1000.0,
                        None,
                    ]
                )
                ys.extend(
                    [
                        y0 / 1000.0,
                        y0 / 1000.0,
                        y1 / 1000.0,
                        y1 / 1000.0,
                        y0 / 1000.0,
                        None,
                    ]
                )

            if by_net:
                visible = _legend_visibility(
                    layer_name, layer_name in bottom_four_visible
                )
                for idx, (net_name, (shapes_x, shapes_y)) in enumerate(
                    sorted(by_net.items())
                ):
                    fig.add_trace(
                        go.Scatter(
                            x=shapes_x,
                            y=shapes_y,
                            mode="lines",
                            fill="toself",
                            fillcolor=color,
                            line=dict(color=color, width=0.5),
                            name=layer_name,
                            legendgroup=layer_name,
                            showlegend=(idx == 0),
                            hoveron="fills",
                            hoverinfo="text",
                            hovertext=f"{layer_name}:{net_name}",
                            hovertemplate="%{hovertext}<extra></extra>",
                            visible=visible,
                        ),
                    )

        elif beol_layer.type == "via":
            by_net: dict[str, tuple[list[float | None], list[float | None]]] = {}
            via_w = beol_layer.min_width or 0.0

            # Determine if this via connects two grid layers for multi-via rendering.
            itf_via_obj = next(
                (v for v in config.itf_stack.vias if v.name == layer_name), None
            )
            is_grid_to_grid = (
                itf_via_obj is not None
                and _is_grid_layer(config, itf_via_obj.from_layer)
                and _is_grid_layer(config, itf_via_obj.to_layer)
            )

            if is_grid_to_grid and itf_via_obj is not None:
                # Multi-via array rendering for grid-to-grid crossings.
                via_side = (
                    math.sqrt(itf_via_obj.area_nm2)
                    if itf_via_obj.area_nm2 > 0
                    else via_w
                )
                cross_pts = _stripe_intersections_with_widths(
                    grid, config, itf_via_obj.from_layer, itf_via_obj.to_layer
                )
                for x_nm, y_nm, net_name, w_a, w_b in cross_pts:
                    if region_nm and not _in_region_point(x_nm, y_nm, region_nm):
                        continue
                    rects = _multi_via_rects(
                        x_nm, y_nm, w_a, w_b, via_side, config.grid.via_min_space_factor
                    )
                    xs, ys = by_net.setdefault(net_name, ([], []))
                    for vx0, vy0, vx1, vy1 in rects:
                        xs.extend(
                            [
                                vx0 / 1000.0,
                                vx1 / 1000.0,
                                vx1 / 1000.0,
                                vx0 / 1000.0,
                                vx0 / 1000.0,
                                None,
                            ]
                        )
                        ys.extend(
                            [
                                vy0 / 1000.0,
                                vy0 / 1000.0,
                                vy1 / 1000.0,
                                vy1 / 1000.0,
                                vy0 / 1000.0,
                                None,
                            ]
                        )

            elif itf_via_obj is not None:
                # Single-via rendering for non-grid-to-grid connections.
                derived_points = _via_points_from_itf_connection(
                    grid, config, layer_name
                )
                if derived_points:
                    for x_nm, y_nm, net_name in derived_points:
                        if region_nm and not _in_region_point(x_nm, y_nm, region_nm):
                            continue
                        half_w = via_w / 2.0
                        x0, y0 = x_nm - half_w, y_nm - half_w
                        x1, y1 = x_nm + half_w, y_nm + half_w
                        xs, ys = by_net.setdefault(net_name, ([], []))
                        xs.extend(
                            [
                                x0 / 1000.0,
                                x1 / 1000.0,
                                x1 / 1000.0,
                                x0 / 1000.0,
                                x0 / 1000.0,
                                None,
                            ]
                        )
                        ys.extend(
                            [
                                y0 / 1000.0,
                                y0 / 1000.0,
                                y1 / 1000.0,
                                y1 / 1000.0,
                                y0 / 1000.0,
                                None,
                            ]
                        )

            if not by_net:
                for via in grid.vias:
                    if via.via_layer != layer_name:
                        continue
                    if region_nm and not _in_region_point(
                        via.node_top.x, via.node_top.y, region_nm
                    ):
                        continue

                    half_w = via.width / 2.0
                    x0, y0 = via.node_top.x - half_w, via.node_top.y - half_w
                    x1, y1 = via.node_top.x + half_w, via.node_top.y + half_w
                    xs, ys = by_net.setdefault(via.net, ([], []))
                    xs.extend(
                        [
                            x0 / 1000.0,
                            x1 / 1000.0,
                            x1 / 1000.0,
                            x0 / 1000.0,
                            x0 / 1000.0,
                            None,
                        ]
                    )
                    ys.extend(
                        [
                            y0 / 1000.0,
                            y0 / 1000.0,
                            y1 / 1000.0,
                            y1 / 1000.0,
                            y0 / 1000.0,
                            None,
                        ]
                    )

            if by_net:
                visible = _legend_visibility(
                    layer_name, layer_name in visible_via_layers
                )
                for idx, (net_name, (shapes_x, shapes_y)) in enumerate(
                    sorted(by_net.items())
                ):
                    fig.add_trace(
                        go.Scatter(
                            x=shapes_x,
                            y=shapes_y,
                            mode="lines",
                            fill="toself",
                            fillcolor=color,
                            line=dict(color=color, width=0.5),
                            name=layer_name,
                            legendgroup=layer_name,
                            showlegend=(idx == 0),
                            hoveron="fills",
                            hoverinfo="text",
                            hovertext=f"{layer_name}:{net_name}",
                            hovertemplate="%{hovertext}<extra></extra>",
                            visible=visible,
                        ),
                    )

    # Cells (chain cells).
    if grid.cells:
        cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
        cell_w_nm = config.distance_to_nm(cell_cfg.width)
        cell_h_nm = config.distance_to_nm(cell_cfg.height)

        cell_shapes_x: list[float | None] = []
        cell_shapes_y: list[float | None] = []
        for cell in grid.cells:
            if region_nm and not _in_region_point(cell.x, cell.y, region_nm):
                continue
            x0, y0 = cell.x - cell_w_nm / 2.0, cell.y - cell_h_nm / 2.0
            x1, y1 = cell.x + cell_w_nm / 2.0, cell.y + cell_h_nm / 2.0
            cell_shapes_x.extend(
                [x0 / 1000.0, x1 / 1000.0, x1 / 1000.0, x0 / 1000.0, x0 / 1000.0, None]
            )
            cell_shapes_y.extend(
                [y0 / 1000.0, y0 / 1000.0, y1 / 1000.0, y1 / 1000.0, y0 / 1000.0, None]
            )

        if cell_shapes_x:
            fig.add_trace(
                go.Scatter(
                    x=cell_shapes_x,
                    y=cell_shapes_y,
                    mode="lines",
                    fill="toself",
                    fillcolor=CELL_COLOR,
                    line=dict(color="rgba(0,0,0,0.5)", width=0.6),
                    name="Cells",
                    hoverinfo="name",
                    visible=_legend_visibility("Cells", True),
                ),
            )

    # Dcap Cells.
    if grid.dcap_cells:
        dcap_cfg = config.get_cell_by_name(
            config.standard_cell_placement.dcap_cells.cell
        )
        dcap_w_nm = config.distance_to_nm(dcap_cfg.width)
        dcap_h_nm = config.distance_to_nm(dcap_cfg.height)

        dcap_shapes_x: list[float | None] = []
        dcap_shapes_y: list[float | None] = []
        for cell in grid.dcap_cells:
            if region_nm and not _in_region_point(cell.x, cell.y, region_nm):
                continue
            x0, y0 = cell.x - dcap_w_nm / 2.0, cell.y - dcap_h_nm / 2.0
            x1, y1 = cell.x + dcap_w_nm / 2.0, cell.y + dcap_h_nm / 2.0
            dcap_shapes_x.extend(
                [x0 / 1000.0, x1 / 1000.0, x1 / 1000.0, x0 / 1000.0, x0 / 1000.0, None]
            )
            dcap_shapes_y.extend(
                [y0 / 1000.0, y0 / 1000.0, y1 / 1000.0, y1 / 1000.0, y0 / 1000.0, None]
            )

        if dcap_shapes_x:
            fig.add_trace(
                go.Scatter(
                    x=dcap_shapes_x,
                    y=dcap_shapes_y,
                    mode="lines",
                    fill="toself",
                    fillcolor=DCAP_COLOR,
                    line=dict(color="rgba(0,0,0,0.5)", width=0.6),
                    name="Dcap Cells",
                    legendgroup="dcap_cells",
                    hoverinfo="name",
                    visible=_legend_visibility("Dcap Cells", True),
                ),
            )

    # Flight lines.
    if grid.cells:
        cell_cfg = config.get_cell_by_name(config.spice_netlist.cell_chains.chain_cell)
        cell_w_nm = config.distance_to_nm(cell_cfg.width)
        cell_h_nm = config.distance_to_nm(cell_cfg.height)
        in_pin = next(
            (p for p in cell_cfg.pins if p.type == "signal" and p.direction == "input"),
            None,
        )
        out_pin = next(
            (
                p
                for p in cell_cfg.pins
                if p.type == "signal" and p.direction == "output"
            ),
            None,
        )
        if in_pin and out_pin:
            in_net_to_cell: dict[str, CellPlacement] = {}
            for c in grid.cells:
                n = c.pin_connections.get(in_pin.name)
                if n:
                    in_net_to_cell[n] = c

            fx: list[float | None] = []
            fy: list[float | None] = []
            for src in grid.cells:
                out_net = src.pin_connections.get(out_pin.name)
                if not out_net:
                    continue
                dst = in_net_to_cell.get(out_net)
                if not dst:
                    continue

                x0, y0 = _pin_xy(src, out_pin.location, cell_w_nm, cell_h_nm)
                x1, y1 = _pin_xy(dst, in_pin.location, cell_w_nm, cell_h_nm)
                if region_nm and not (
                    _in_region_point(x0, y0, region_nm)
                    or _in_region_point(x1, y1, region_nm)
                ):
                    continue

                fx.extend([x0 / 1000.0, x1 / 1000.0, None])
                fy.extend([y0 / 1000.0, y1 / 1000.0, None])

            if fx:
                fig.add_trace(
                    go.Scatter(
                        x=fx,
                        y=fy,
                        mode="lines",
                        line=dict(color=FLIGHT_COLOR, width=1.0, dash="dot"),
                        name="Flight Lines",
                        legendgroup="flight_lines",
                        hoverinfo="name",
                        visible=_legend_visibility("Flight Lines", False),
                    ),
                )

    # PLOC points. Add last so they appear at the bottom of legend stack.
    if grid.plocs:
        diameter_nm = config.distance_to_nm(config.ploc.visualizer_render_diameter)
        radius_um = max(diameter_nm / 2000.0, 1e-6)
        by_net: dict[str, list[tuple[float, float]]] = {}
        for p in grid.plocs:
            if region_nm and not _in_region_point(p.x, p.y, region_nm):
                continue
            by_net.setdefault(p.net, []).append((p.x / 1000.0, p.y / 1000.0))

        n_circle_pts = 24
        for net_name, centers in sorted(by_net.items()):
            cx: list[float | None] = []
            cy: list[float | None] = []
            for x0, y0 in centers:
                for i in range(n_circle_pts + 1):
                    theta = (2.0 * math.pi * i) / n_circle_pts
                    cx.append(x0 + radius_um * math.cos(theta))
                    cy.append(y0 + radius_um * math.sin(theta))
                cx.append(None)
                cy.append(None)

            color = _ploc_color(config, net_name)
            fig.add_trace(
                go.Scatter(
                    x=cx,
                    y=cy,
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=dict(color=color, width=0.8),
                    name=f"PLOC:{net_name}",
                    legendgroup=f"PLOC:{net_name}",
                    showlegend=True,
                    hoveron="fills",
                    hoverinfo="text",
                    hovertext=f"PLOC:{net_name}",
                    hovertemplate="%{hovertext}<extra></extra>",
                    visible=(
                        _legend_visibility(f"PLOC:{net_name}", True)
                        if configured_visible_objects is None
                        or "PLOC" not in configured_visible_objects
                        else True
                    ),
                ),
            )

    fig.update_layout(
        title=f"Power Grid View - {config.beol_stack.technology}",
        width=1050,
        height=1050,
        legend=dict(orientation="v"),
        xaxis=dict(title="X (um)"),
        yaxis=dict(title="Y (um)", scaleanchor="x", scaleratio=1),
    )

    # Build the cross-section as an independent figure.
    fig2 = _build_cross_section(config)

    if output_path:
        # Two independent plots in a single HTML file per AGENTS.md.
        div1 = pyo.plot(fig, include_plotlyjs=True, output_type="div")
        div2 = pyo.plot(fig2, include_plotlyjs=False, output_type="div")
        html_content = f"<html><body>{div1}{div2}</body></html>"
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(html_content)

        if open_browser:
            import webbrowser

            webbrowser.open(f"file://{output_path.resolve()}")

    if save_image_layer:
        # Reserved for future static-image extraction implementation.
        pass
