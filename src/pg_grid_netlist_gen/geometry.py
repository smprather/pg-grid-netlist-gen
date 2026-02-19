"""Geometry data model for the power grid."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


@dataclass
class Node:
    """A point in the grid where elements connect."""

    name: str
    x: float  # nm
    y: float  # nm
    z: float  # nm (from BEOL stack, bottom of metal)
    layer: str
    net: str


@dataclass
class Segment:
    """A wire segment between two nodes on the same layer."""

    node_a: Node
    node_b: Node
    layer: str
    width: float  # nm
    thickness: float  # nm
    length: float  # nm
    net: str
    resistance: float = 0.0  # Ω
    cap_plate: float = 0.0  # F
    cap_fringe: float = 0.0  # F


@dataclass
class ViaConnection:
    """A via connecting two nodes on adjacent layers."""

    node_top: Node
    node_bot: Node
    via_layer: str
    width: float  # nm
    height: float  # nm (via thickness)
    net: str
    resistance: float = 0.0  # Ω


@dataclass
class Stripe:
    """A full metal stripe on a grid layer before segmentation."""

    layer: str
    direction: Literal["vertical", "horizontal"]
    position: float  # nm — X for vertical, Y for horizontal
    start: float  # nm
    end: float  # nm
    width: float  # nm
    net: str


@dataclass
class Staple:
    """A square metal pad on a staple layer."""

    layer: str
    x: float  # nm (center)
    y: float  # nm (center)
    size: float  # nm (side length)
    net: str


@dataclass
class CellPlacement:
    """A placed standard cell instance."""

    instance_name: str
    cell_name: str
    x: float  # nm
    y: float  # nm
    row: int
    flipped: bool
    pin_connections: dict[str, str] = field(default_factory=dict)


@dataclass
class PlocPoint:
    """A snapped PLOC supply connection point."""

    node_name: str
    x: float  # nm
    y: float  # nm
    layer: str
    net: str


@dataclass
class Grid:
    """Complete generated grid, the central data structure."""

    stripes: list[Stripe] = field(default_factory=list)
    staples: list[Staple] = field(default_factory=list)
    segments: list[Segment] = field(default_factory=list)
    vias: list[ViaConnection] = field(default_factory=list)
    cells: list[CellPlacement] = field(default_factory=list)
    dcap_cells: list[CellPlacement] = field(default_factory=list)
    plocs: list[PlocPoint] = field(default_factory=list)
    nodes: dict[str, Node] = field(default_factory=dict)
