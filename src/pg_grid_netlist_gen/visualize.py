"""Plotly 3D visualization of the power grid."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import plotly.graph_objects as go

from pg_grid_netlist_gen.config import Config
from pg_grid_netlist_gen.geometry import Grid, Segment, Staple, ViaConnection

# Color scheme per layer tier
LAYER_COLORS: dict[str, str] = {
    "M1": "rgb(30, 100, 220)",
    "MINT1": "rgb(40, 160, 80)",
    "MINT2": "rgb(50, 170, 90)",
    "MINT3": "rgb(60, 180, 100)",
    "MINT4": "rgb(70, 190, 110)",
    "MINT5": "rgb(80, 200, 120)",
    "MSMG1": "rgb(220, 140, 40)",
    "MSMG2": "rgb(230, 150, 50)",
    "MSMG3": "rgb(240, 160, 60)",
    "MSMG4": "rgb(250, 170, 70)",
    "MSMG5": "rgb(255, 180, 80)",
    "MG1": "rgb(200, 50, 50)",
    "MG2": "rgb(220, 70, 70)",
    "VG1": "rgb(180, 40, 40)",
}

VIA_COLOR = "rgb(160, 160, 160)"
CELL_COLOR = "rgb(220, 200, 50)"


def _box_mesh(
    x0: float, y0: float, z0: float,
    x1: float, y1: float, z1: float,
) -> tuple[list[float], list[float], list[float], list[int], list[int], list[int]]:
    """Generate vertices and triangle indices for a 3D box."""
    # 8 vertices
    vx = [x0, x1, x1, x0, x0, x1, x1, x0]
    vy = [y0, y0, y1, y1, y0, y0, y1, y1]
    vz = [z0, z0, z0, z0, z1, z1, z1, z1]

    # 12 triangles (2 per face)
    ti = [0, 0, 4, 4, 0, 0, 1, 1, 0, 0, 3, 3]
    tj = [1, 2, 5, 6, 1, 4, 2, 5, 3, 4, 2, 6]
    tk = [2, 3, 6, 7, 4, 5, 5, 6, 4, 7, 6, 7]

    return vx, vy, vz, ti, tj, tk


def _merge_meshes(
    boxes: list[tuple[list[float], list[float], list[float], list[int], list[int], list[int]]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Merge multiple box meshes into a single mesh."""
    all_x, all_y, all_z = [], [], []
    all_i, all_j, all_k = [], [], []
    offset = 0

    for vx, vy, vz, ti, tj, tk in boxes:
        all_x.extend(vx)
        all_y.extend(vy)
        all_z.extend(vz)
        all_i.extend(idx + offset for idx in ti)
        all_j.extend(idx + offset for idx in tj)
        all_k.extend(idx + offset for idx in tk)
        offset += len(vx)

    return (
        np.array(all_x), np.array(all_y), np.array(all_z),
        np.array(all_i), np.array(all_j), np.array(all_k),
    )


def render_grid(
    grid: Grid,
    config: Config,
    output_path: str | Path,
    z_exaggeration: float = 50.0,
    viz_region: tuple[float, float, float, float] | None = None,
    open_browser: bool = False,
) -> None:
    """Render the grid as a 3D HTML visualization."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert viz_region from microns to nm if specified
    region_nm: tuple[float, float, float, float] | None = None
    if viz_region:
        region_nm = (
            viz_region[0] * 1000, viz_region[1] * 1000,
            viz_region[2] * 1000, viz_region[3] * 1000,
        )

    fig = go.Figure()

    # Group segments by layer
    layer_segments: dict[str, list[Segment]] = {}
    for seg in grid.segments:
        if region_nm and not _in_region_segment(seg, region_nm):
            continue
        layer_segments.setdefault(seg.layer, []).append(seg)

    # Render segment boxes per layer
    for layer_name, segs in layer_segments.items():
        boxes = []
        for seg in segs:
            box = _segment_to_box(seg, z_exaggeration)
            if box:
                boxes.append(box)

        if not boxes:
            continue

        x, y, z, i, j, k = _merge_meshes(boxes)
        color = LAYER_COLORS.get(layer_name, "rgb(128, 128, 128)")

        fig.add_trace(go.Mesh3d(
            x=x / 1000, y=y / 1000, z=z,  # convert to μm for display
            i=i, j=j, k=k,
            color=color,
            opacity=0.6,
            name=f"{layer_name} (segments)",
            hoverinfo="name",
            showlegend=True,
        ))

    # Render staples per layer
    layer_staples: dict[str, list[Staple]] = {}
    for staple in grid.staples:
        if region_nm and not _in_region_point(staple.x, staple.y, region_nm):
            continue
        layer_staples.setdefault(staple.layer, []).append(staple)

    for layer_name, stps in layer_staples.items():
        beol_layer = config.get_beol_layer(layer_name)
        z_coords = _compute_z_for_layer(config, layer_name)
        if z_coords is None:
            continue
        z_bot, z_top = z_coords

        boxes = []
        for staple in stps:
            half = staple.size / 2.0
            boxes.append(_box_mesh(
                staple.x - half, staple.y - half, z_bot * z_exaggeration,
                staple.x + half, staple.y + half, z_top * z_exaggeration,
            ))

        if not boxes:
            continue

        x, y, z, i, j, k = _merge_meshes(boxes)
        color = LAYER_COLORS.get(layer_name, "rgb(128, 128, 128)")

        fig.add_trace(go.Mesh3d(
            x=x / 1000, y=y / 1000, z=z,
            i=i, j=j, k=k,
            color=color,
            opacity=0.5,
            name=f"{layer_name} (staples)",
            hoverinfo="name",
            showlegend=True,
        ))

    # Render vias
    via_boxes: dict[str, list] = {}
    for via in grid.vias:
        if region_nm and not _in_region_point(via.node_top.x, via.node_top.y, region_nm):
            continue
        via_boxes.setdefault(via.via_layer, []).append(via)

    for via_layer, via_list in via_boxes.items():
        z_coords = _compute_z_for_layer(config, via_layer)
        if z_coords is None:
            continue
        z_bot, z_top = z_coords

        boxes = []
        for via in via_list:
            half_w = via.width / 2.0
            boxes.append(_box_mesh(
                via.node_top.x - half_w, via.node_top.y - half_w, z_bot * z_exaggeration,
                via.node_top.x + half_w, via.node_top.y + half_w, z_top * z_exaggeration,
            ))

        if not boxes:
            continue

        x, y, z, i, j, k = _merge_meshes(boxes)

        fig.add_trace(go.Mesh3d(
            x=x / 1000, y=y / 1000, z=z,
            i=i, j=j, k=k,
            color=VIA_COLOR,
            opacity=0.4,
            name=f"{via_layer} (vias)",
            hoverinfo="name",
            showlegend=True,
        ))

    # Render standard cells as blocks
    if grid.cells:
        cell_cfg = config.standard_cells[0]
        cell_w = cell_cfg.size["x"]
        cell_h = cell_cfg.size["y"]
        # Cells sit below M1 — use a small z range below the BEOL
        cell_z_bot = -100 * z_exaggeration
        cell_z_top = 0

        boxes = []
        for cell in grid.cells:
            if region_nm and not _in_region_point(cell.x, cell.y, region_nm):
                continue
            boxes.append(_box_mesh(
                cell.x - cell_w / 2, cell.y - cell_h / 2, cell_z_bot,
                cell.x + cell_w / 2, cell.y + cell_h / 2, cell_z_top,
            ))

        if boxes:
            x, y, z, i, j, k = _merge_meshes(boxes)
            fig.add_trace(go.Mesh3d(
                x=x / 1000, y=y / 1000, z=z,
                i=i, j=j, k=k,
                color=CELL_COLOR,
                opacity=0.5,
                name="Standard Cells",
                hoverinfo="name",
                showlegend=True,
            ))

    fig.update_layout(
        title=f"Power Grid - {config.beol_stack.technology} ({config.beol_stack.node})",
        scene=dict(
            xaxis_title="X (μm)",
            yaxis_title="Y (μm)",
            zaxis_title="Z (nm, exaggerated)",
            aspectmode="data",
        ),
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )

    fig.write_html(str(output_path), include_plotlyjs="cdn")

    if open_browser:
        import webbrowser
        webbrowser.open(f"file://{output_path.resolve()}")


def _segment_to_box(
    seg: Segment,
    z_exaggeration: float,
) -> tuple[list[float], list[float], list[float], list[int], list[int], list[int]] | None:
    """Convert a segment to a 3D box."""
    half_w = seg.width / 2.0
    z_bot = seg.node_a.z * z_exaggeration
    z_top = (seg.node_a.z + seg.thickness) * z_exaggeration

    if seg.node_a.x == seg.node_b.x:
        # Vertical segment (runs along Y)
        x0 = seg.node_a.x - half_w
        x1 = seg.node_a.x + half_w
        y0 = min(seg.node_a.y, seg.node_b.y)
        y1 = max(seg.node_a.y, seg.node_b.y)
    else:
        # Horizontal segment (runs along X)
        x0 = min(seg.node_a.x, seg.node_b.x)
        x1 = max(seg.node_a.x, seg.node_b.x)
        y0 = seg.node_a.y - half_w
        y1 = seg.node_a.y + half_w

    if x0 == x1 or y0 == y1:
        return None

    return _box_mesh(x0, y0, z_bot, x1, y1, z_top)


def _compute_z_for_layer(
    config: Config,
    layer_name: str,
) -> tuple[float, float] | None:
    """Compute z_bottom, z_top for a layer from the BEOL stack."""
    layers_bottom_up = list(reversed(config.beol_stack.layers))
    z = 0.0
    for layer in layers_bottom_up:
        if layer.type == "substrate":
            z = 0.0
            continue
        z_bottom = z
        z_top = z + layer.thickness
        if layer.name == layer_name:
            return z_bottom, z_top
        z = z_top
    return None


def _in_region_segment(
    seg: Segment,
    region: tuple[float, float, float, float],
) -> bool:
    """Check if a segment overlaps with the region."""
    x_min, y_min, x_max, y_max = region
    sx_min = min(seg.node_a.x, seg.node_b.x)
    sx_max = max(seg.node_a.x, seg.node_b.x)
    sy_min = min(seg.node_a.y, seg.node_b.y)
    sy_max = max(seg.node_a.y, seg.node_b.y)
    return sx_max >= x_min and sx_min <= x_max and sy_max >= y_min and sy_min <= y_max


def _in_region_point(x: float, y: float, region: tuple[float, float, float, float]) -> bool:
    """Check if a point is within the region."""
    return region[0] <= x <= region[2] and region[1] <= y <= region[3]
