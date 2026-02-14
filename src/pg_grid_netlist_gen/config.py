"""Pydantic models for YAML configuration parsing."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, model_validator


class Material(BaseModel):
    name: str
    conductivity: float  # MS/m
    relative_permittivity: float


class StandardCellConfig(BaseModel):
    cell: str
    size: dict[str, float]  # x, y in nm
    instance_pins: list[str]
    distance_apart_um: float = 5000.0


class BeolLayer(BaseModel):
    name: str
    type: Literal["metal", "via", "oxide", "substrate"]
    thickness: float  # nm
    min_width: float | None = None
    pitch: float | None = None


class GridLayerConfig(BaseModel):
    type: Literal["grid", "staple"]
    direction: Literal["vertical", "horizontal"]
    pitch: float  # factor (multiplied by BEOL pitch)
    width: float | None = None  # nm, defaults to BEOL min_width
    material: str | None = None


class GridSize(BaseModel):
    x: float  # nm
    y: float  # nm


class GridConfig(BaseModel):
    nets: list[str]
    cycle: list[str] | None = None
    grid_size: GridSize
    layers: dict[str, GridLayerConfig]

    @model_validator(mode="before")
    @classmethod
    def merge_common_and_remove(cls, data: dict[str, Any]) -> dict[str, Any]:
        layers = data.get("layers", {})
        common = layers.pop("COMMON", {})
        for name, layer_cfg in layers.items():
            for key, val in common.items():
                if key not in layer_cfg:
                    layer_cfg[key] = val
        return data

    @model_validator(mode="after")
    def set_cycle_default(self) -> GridConfig:
        if self.cycle is None:
            self.cycle = list(self.nets)
        return self


class BeolStack(BaseModel):
    technology: str
    node: str
    units: float  # 1e-9 for nm
    layers: list[BeolLayer]


class Config(BaseModel):
    materials: list[Material]
    standard_cells: list[StandardCellConfig]
    grid: GridConfig
    beol_stack: BeolStack

    def get_material(self, name: str) -> Material:
        for m in self.materials:
            if m.name == name:
                return m
        raise ValueError(f"Material '{name}' not found")

    def get_beol_layer(self, name: str) -> BeolLayer:
        for layer in self.beol_stack.layers:
            if layer.name == name:
                return layer
        raise ValueError(f"BEOL layer '{name}' not found")

    def grid_layer_order(self) -> list[str]:
        """Return grid layer names in top-to-bottom order (matching BEOL stack order)."""
        beol_names = [l.name for l in self.beol_stack.layers]
        grid_names = list(self.grid.layers.keys())
        return [n for n in beol_names if n in grid_names]


def load_config(path: str | Path) -> Config:
    """Load and validate a YAML configuration file."""
    with open(path) as f:
        raw = yaml.safe_load(f)
    # The top-level 'units' key is informational, not part of the model
    raw.pop("units", None)
    return Config.model_validate(raw)
