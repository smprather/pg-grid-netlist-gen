"""Plotly 2D visualization of the power grid."""

from __future__ import annotations

from pathlib import Path

import plotly.graph_objects as go

from pg_grid_netlist_gen.config import Config
from pg_grid_netlist_gen.geometry import Grid

# Color scheme per layer tier (using rgba for fill color with some transparency)
LAYER_COLORS: dict[str, str] = {
    "M1": "rgba(30, 100, 220, 0.7)",
    "MINT1": "rgba(40, 160, 80, 0.7)",
    "MINT2": "rgba(50, 170, 90, 0.7)",
    "MINT3": "rgba(60, 180, 100, 0.7)",
    "MINT4": "rgba(70, 190, 110, 0.7)",
    "MINT5": "rgba(80, 200, 120, 0.7)",
    "MSMG1": "rgba(220, 140, 40, 0.7)",
    "MSMG2": "rgba(230, 150, 50, 0.7)",
    "MSMG3": "rgba(240, 160, 60, 0.7)",
    "MSMG4": "rgba(250, 170, 70, 0.7)",
    "MSMG5": "rgba(255, 180, 80, 0.7)",
    "MG1": "rgba(200, 50, 50, 0.7)",
    "MG2": "rgba(220, 70, 70, 0.7)",
    # Explicit colors for via layers
    "VG1": "rgba(180, 40, 40, 0.7)",
    "VINT1": "rgba(180, 40, 40, 0.7)",
    "VINT2": "rgba(180, 40, 40, 0.7)",
    "VINT3": "rgba(180, 40, 40, 0.7)",
    "VINT4": "rgba(180, 40, 40, 0.7)",
    "VINT5": "rgba(180, 40, 40, 0.7)",
    "VSMG1": "rgba(180, 40, 40, 0.7)",
    "VSMG2": "rgba(180, 40, 40, 0.7)",
    "VSMG3": "rgba(180, 40, 40, 0.7)",
    "VSMG4": "rgba(180, 40, 40, 0.7)",
    "VSMG5": "rgba(180, 40, 40, 0.7)",
    "V1": "rgba(180, 40, 40, 0.7)",
    "V0": "rgba(180, 40, 40, 0.7)",
}

VIA_COLOR = "rgba(160, 160, 160, 0.8)" # This is a fallback if a via layer is not in LAYER_COLORS
CELL_COLOR = "rgba(220, 200, 50, 0.5)"


def render_grid(
    grid: Grid,
    config: Config,
    output_path: str | Path,
    viz_region: tuple[float, float, float, float] | None = None,
    open_browser: bool = False,
) -> None:
    """Render the grid as a 2D interactive HTML visualization."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    region_nm: tuple[float, float, float, float] | None = None
    if viz_region:
        region_nm = (
            viz_region[0] * 1000, viz_region[1] * 1000,
            viz_region[2] * 1000, viz_region[3] * 1000,
        )

    fig = go.Figure()

    # List to store traces as they are generated, maintaining BEOL stack order
    ordered_traces = []
    
    # Process layers based on BEOL stack order (highest to lowest)
    for beol_layer in config.beol_stack.layers:
        layer_name = beol_layer.name
        
        if beol_layer.type == "metal":
            # Check if this metal layer is actually part of the grid (has segments or staples)
            has_segments = any(seg.layer == layer_name for seg in grid.segments)
            has_staples = any(staple.layer == layer_name for staple in grid.staples)

            if not (has_segments or has_staples):
                continue # Skip if no components for this metal layer in the grid

            current_layer_xs = []
            current_layer_ys = []
            color = LAYER_COLORS.get(layer_name, "rgba(128, 128, 128, 0.7)")

            # Segments
            for seg in grid.segments:
                if seg.layer == layer_name:
                    if region_nm and not _in_region_segment(seg, region_nm):
                        continue
                    half_w = seg.width / 2.0
                    if seg.node_a.x == seg.node_b.x: # Vertical
                        x0, x1 = (seg.node_a.x - half_w), (seg.node_a.x + half_w)
                        y0, y1 = min(seg.node_a.y, seg.node_b.y), max(seg.node_a.y, seg.node_b.y)
                    else: # Horizontal
                        x0, x1 = min(seg.node_a.x, seg.node_b.x), max(seg.node_a.x, seg.node_b.x)
                        y0, y1 = (seg.node_a.y - half_w), (seg.node_a.y + half_w)
                    
                    current_layer_xs.extend([x0/1000, x1/1000, x1/1000, x0/1000, x0/1000, None])
                    current_layer_ys.extend([y0/1000, y0/1000, y1/1000, y1/1000, y0/1000, None])

            # Staples
            for staple in grid.staples:
                if staple.layer == layer_name:
                    if region_nm and not _in_region_point(staple.x, staple.y, region_nm):
                        continue
                    half_s = staple.size / 2.0
                    x0, y0 = staple.x - half_s, staple.y - half_s
                    x1, y1 = staple.x + half_s, staple.y + half_s

                    current_layer_xs.extend([x0/1000, x1/1000, x1/1000, x0/1000, x0/1000, None])
                    current_layer_ys.extend([y0/1000, y0/1000, y1/1000, y1/1000, y0/1000, None])
            
            if current_layer_xs:
                ordered_traces.append(
                    go.Scatter(
                        x=current_layer_xs, y=current_layer_ys,
                        mode='lines',
                        fill='toself',
                        fillcolor=color,
                        line=dict(color=color, width=0.5),
                        name=layer_name,
                        hoverinfo='name',
                        visible=False # Hidden by default
                    )
                )

        elif beol_layer.type == "via":
            # Check if this via layer has any vias placed
            has_vias = any(via.via_layer == layer_name for via in grid.vias)
            if not has_vias:
                continue # Skip if no vias for this layer

            current_via_xs = []
            current_via_ys = []
            color = LAYER_COLORS.get(layer_name, VIA_COLOR)

            for via in grid.vias:
                if via.via_layer == layer_name:
                    if region_nm and not _in_region_point(via.node_top.x, via.node_top.y, region_nm):
                        continue
                    half_w = via.width / 2.0
                    x0, y0 = via.node_top.x - half_w, via.node_top.y - half_w
                    x1, y1 = via.node_top.x + half_w, via.node_top.y + half_w
                    current_via_xs.extend([x0/1000, x1/1000, x1/1000, x0/1000, x0/1000, None])
                    current_via_ys.extend([y0/1000, y0/1000, y1/1000, y1/1000, y0/1000, None])
            
            if current_via_xs:
                ordered_traces.append(
                    go.Scatter(
                        x=current_via_xs, y=current_via_ys,
                        mode='lines',
                        fill='toself',
                        fillcolor=color,
                        line=dict(color=color, width=0.5),
                        name=layer_name,
                        hoverinfo='name',
                        visible=False # Hidden by default
                    )
                )

    # Process Cells (always added last if present)
    if grid.cells:
        cell_cfg = config.standard_cells[0]
        cell_w, cell_h = cell_cfg.size["x"], cell_cfg.size["y"]
        
        cell_xs = []
        cell_ys = []
        for cell in grid.cells:
            if region_nm and not _in_region_point(cell.x, cell.y, region_nm):
                continue
            x0, y0 = cell.x - cell_w / 2, cell.y - cell_h / 2
            x1, y1 = cell.x + cell_w / 2, cell.y + cell_h / 2
            cell_xs.extend([x0/1000, x1/1000, x1/1000, x0/1000, x0/1000, None])
            cell_ys.extend([y0/1000, y0/1000, y1/1000, y1/1000, y0/1000, None])
        
        if cell_xs:
            ordered_traces.append(
                go.Scatter(
                    x=cell_xs, y=cell_ys,
                    mode='lines',
                    fill='toself',
                    fillcolor=CELL_COLOR,
                    line=dict(color="rgba(0,0,0,0.5)", width=0.5),
                    name="Cells",
                    hoverinfo='name',
                    visible=False # Hidden by default
                )
            )

    # Add all traces to the figure in their determined order
    for trace in ordered_traces:
        fig.add_trace(trace)

    # Create dropdown buttons
    buttons = [
        dict(
            label="None",
            method="restyle",
            args=["visible", [False] * len(fig.data)],
        ),
        dict(
            label="All",
            method="restyle",
            args=["visible", [True] * len(fig.data)],
        ),
    ]

    # Add buttons for individual layers based on the order in ordered_traces
    for i, trace in enumerate(fig.data):
        visibility = [False] * len(fig.data)
        visibility[i] = True # Set only this trace to visible
        
        buttons.append(
            dict(
                label=trace.name,
                method="restyle",
                args=["visible", visibility],
            )
        )

    updatemenus = [
        dict(
            active=0, # "None" is active by default
            buttons=buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.01,
            xanchor="left",
            y=1.1,
            yanchor="top",
        )
    ]

    fig.update_layout(
        title=f"2D View - {config.beol_stack.technology} ({config.beol_stack.node})",
        xaxis_title="X (μm)",
        yaxis_title="Y (μm)",
        updatemenus=updatemenus,
        width=800,
        height=800,
        xaxis=dict(scaleanchor="y", scaleratio=1),
        yaxis=dict(autorange="reversed"),
    )

    fig.write_html(str(output_path), include_plotlyjs="cdn")

    if open_browser:
        import webbrowser
        webbrowser.open(f"file://{output_path.resolve()}")

def _in_region_segment(seg, region):
    x_min, y_min, x_max, y_max = region
    sx_min = min(seg.node_a.x, seg.node_b.x)
    sx_max = max(seg.node_a.x, seg.node_b.x)
    sy_min = min(seg.node_a.y, seg.node_b.y)
    sy_max = max(seg.node_a.y, seg.node_b.y)
    return sx_max >= x_min and sx_min <= x_max and sy_max >= y_min and sy_min <= y_max

def _in_region_point(x, y, region):
    return region[0] <= x <= region[2] and region[1] <= y <= region[3]
