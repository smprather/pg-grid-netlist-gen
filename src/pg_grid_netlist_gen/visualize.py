"""Plotly 2D visualization of the power grid."""

from __future__ import annotations

import math
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots

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


def _add_cross_section(fig: go.Figure, config: Config) -> None:
    """Add a stack cross-section in subplot row 2."""

    def _add_rect(
        x0: float, x1: float, z_bottom: float, z_top: float, color: str
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
                showlegend=False,
                hoverinfo="skip",
            ),
            row=2,
            col=1,
        )

    def _add_label(text: str, z_bottom: float, z_top: float) -> None:
        fig.add_annotation(
            x=0.5,
            y=((z_bottom + z_top) / 2.0) / 1000.0,
            text=text,
            showarrow=False,
            xref="x2",
            yref="y2",
            font=dict(size=10),
        )

    dielectrics = {d.name: d for d in config.itf_stack.dielectrics}
    conductors_bottom_up = list(reversed(config.itf_stack.conductors))

    # Draw full-width metals and dielectrics.
    z = 0.0
    metal_z: dict[str, tuple[float, float]] = {}
    for metal in conductors_bottom_up:
        m_bottom = z
        m_top = m_bottom + metal.thickness_nm
        metal_z[metal.name] = (m_bottom, m_top)
        _add_rect(0.0, 1.0, m_bottom, m_top, _get_color(metal.name))
        _add_label(metal.name, m_bottom, m_top)
        z = m_top

        diel_name = f"{metal.name}_diel"
        dielectric = dielectrics.get(diel_name)
        if dielectric is not None and dielectric.thickness_nm > 0:
            d_bottom = z
            d_top = d_bottom + dielectric.thickness_nm
            _add_rect(0.0, 1.0, d_bottom, d_top, OXIDE_COLOR)
            _add_label(diel_name, d_bottom, d_top)
            z = d_top

    # Build sample via "plugs" centered between connected metal layers.
    via_entries: list[tuple[object, float, float, float]] = []
    for via in config.itf_stack.vias:
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
            # Ensure a visible plug even if layer spacing data collapses.
            via_top = via_bottom + max(via_side, 1.0)

        via_entries.append((via, via_bottom, via_top, via_side))

    # Normalize via widths to the top-most via layer's side length.
    ref_side = 0.0
    if via_entries:
        top_entry = max(
            via_entries, key=lambda t: t[1]
        )  # highest via by stack position
        ref_side = top_entry[3]
    if ref_side <= 0.0:
        ref_side = max((entry[3] for entry in via_entries), default=1.0)
    if ref_side <= 0.0:
        ref_side = 1.0

    for via, via_bottom, via_top, via_side in via_entries:
        # Top-most via renders at width 0.12; others scale relative to it.
        width_norm = via_side / ref_side
        half_w = max(0.001, min(0.45, 0.06 * width_norm))
        via_x0, via_x1 = 0.5 - half_w, 0.5 + half_w
        _add_rect(via_x0, via_x1, via_bottom, via_top, _get_color(via.name))
        _add_label(via.name, via_bottom, via_top)


def render_grid(
    grid: Grid,
    config: Config,
    output_path: str | Path | None,
    viz_region: tuple[float, float, float, float] | None = None,
    open_browser: bool = False,
    save_image_layer: str | None = None,
    output_dir: Path | None = None,
) -> None:
    """Render combined 2D grid and layer-stack cross-section to HTML."""
    if not (output_path or save_image_layer):
        return

    output_dir = output_dir or (
        Path(output_path).parent if output_path else Path("output")
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    region_nm = tuple(c * 1000 for c in viz_region) if viz_region else None

    fig = make_subplots(
        rows=2,
        cols=1,
        vertical_spacing=0.08,
        row_heights=[0.72, 0.84],
        subplot_titles=("2D Grid and Cells", "ITF Layer Cross-Section"),
    )

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
                visible = True if layer_name in bottom_four_visible else "legendonly"
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
                        row=1,
                        col=1,
                    )

        elif beol_layer.type == "via":
            by_net: dict[str, tuple[list[float | None], list[float | None]]] = {}
            via_w = beol_layer.min_width or 0.0
            derived_points = _via_points_from_itf_connection(grid, config, layer_name)
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
            else:
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
                visible = True if layer_name in visible_via_layers else "legendonly"
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
                        row=1,
                        col=1,
                    )

    # Cells.
    if grid.cells:
        cell_cfg = config.standard_cells[0]
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
                    visible=True,
                ),
                row=1,
                col=1,
            )

    # Flight lines.
    if grid.cells:
        cell_cfg = config.standard_cells[0]
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
                        visible="legendonly",
                    ),
                    row=1,
                    col=1,
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
                    visible=True,
                ),
                row=1,
                col=1,
            )

    _add_cross_section(fig, config)

    fig.update_layout(
        title=f"Power Grid View - {config.beol_stack.technology}",
        width=1050,
        height=1872,
        legend=dict(orientation="v"),
    )

    fig.update_xaxes(title_text="X (um)", row=1, col=1)
    fig.update_yaxes(title_text="Y (um)", row=1, col=1, scaleanchor="x", scaleratio=1)

    fig.update_xaxes(showticklabels=False, title_text="", range=[0, 1], row=2, col=1)
    fig.update_yaxes(title_text="Height (um)", row=2, col=1)

    if output_path:
        # fig.write_html(str(output_path), include_plotlyjs="cdn")
        fig.write_html(str(output_path))
        if open_browser:
            import webbrowser

            webbrowser.open(f"file://{Path(output_path).resolve()}")

    if save_image_layer:
        # Reserved for future static-image extraction implementation.
        pass
