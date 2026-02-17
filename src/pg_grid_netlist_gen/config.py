"""Pydantic models for YAML configuration parsing."""

from __future__ import annotations

from pathlib import Path
import random
from typing import Any, Literal

import yaml
from pydantic import BaseModel, PrivateAttr, model_validator

from pg_grid_netlist_gen.itf_parser import ItfStack, parse_itf


class BeolLayer(BaseModel):
    name: str
    type: Literal["metal", "via", "oxide", "substrate"]
    thickness: float  # nm
    min_width: float | None = None
    pitch: float | None = None
    resistance_per_square: float | None = None


class BeolStack(BaseModel):
    technology: str
    layers: list[BeolLayer]


class ItfConfig(BaseModel):
    file: Path
    units: str


class UnitsConfig(BaseModel):
    distance: str
    capacitance: str
    resistance: str
    time: str


class PinConfig(BaseModel):
    name: str
    type: Literal["signal", "power", "ground"]
    direction: Literal["input", "output", "inout"] | None = None
    location: Literal["left", "right", "top", "bottom"]


class StandardCellConfig(BaseModel):
    name: str
    height: float
    width: float
    pins: list[PinConfig]
    spice_port_order: str
    unateness: Literal["positive", "negative"] = "positive"
    spice_netlist_file: str | None = None


class PlocConfig(BaseModel):
    pitch: float
    staggered: bool
    offset_from_origin: dict[str, float]
    visualizer_render_diameter: float = 1.0


class VisualizerConfig(BaseModel):
    initial_visible_objects: list[str] | None = None


class SpiceNetlistConfig(BaseModel):
    standard_cell_output_load: dict[str, Any]
    transient_simulation: dict[str, Any]
    chain_input_stimulus: dict[str, Any]
    instance_chains: dict[str, Any]
    ir_drop_measurement: dict[str, Any]


class StandardCellPlacementConfig(BaseModel):
    row_height: float
    site_width: float
    min_space: dict[str, float]
    stagger_row_start: dict[str, Any]


class PowerNetConfig(BaseModel):
    name: str
    voltage: float


class GroundNetConfig(BaseModel):
    name: str


class PgNetsConfig(BaseModel):
    power: PowerNetConfig
    ground: GroundNetConfig


class LayerUsageConfig(BaseModel):
    type: Literal["g", "s"]
    width: float
    pitch: float


class GridConfig(BaseModel):
    size: dict[str, int]
    layer_usage: dict[str, LayerUsageConfig]


class Config(BaseModel):
    random_seed_value: float | int | Literal["random"]
    itf: ItfConfig
    units: UnitsConfig
    standard_cells: list[StandardCellConfig]
    ploc: PlocConfig
    visualizer: VisualizerConfig | None = None
    spice_netlist: SpiceNetlistConfig
    standard_cell_placement: StandardCellPlacementConfig
    pg_nets: PgNetsConfig
    grid: GridConfig
    beol_thickness: float

    _beol_stack: BeolStack | None = PrivateAttr(default=None)
    _itf_stack: ItfStack | None = PrivateAttr(default=None)
    _config_path: Path | None = PrivateAttr(default=None)
    _effective_seed: int | float | None = PrivateAttr(default=None)

    @property
    def itf_stack(self) -> ItfStack:
        if self._itf_stack is None:
            _ = self.beol_stack
        return self._itf_stack

    @property
    def beol_stack(self) -> BeolStack:
        if self._beol_stack is None:
            self._itf_stack = parse_itf(self.itf.file, self.itf.units)
            self._beol_stack = self._create_beol_stack_from_itf(self._itf_stack)
        return self._beol_stack

    @property
    def effective_seed(self) -> int | float | None:
        return self._effective_seed

    def make_rng(self) -> random.Random:
        """Return a run RNG and capture the effective seed for reporting."""
        if self._effective_seed is None:
            if self.random_seed_value == "random":
                self._effective_seed = random.SystemRandom().randrange(0, 2**63)
            else:
                self._effective_seed = float(self.random_seed_value)
        return random.Random(self._effective_seed)

    @property
    def lowest_metal_layer_name(self) -> str:
        if not self.itf_stack.conductors:
            raise ValueError("ITF file contains no conductors")
        # ITF conductors are defined top-to-bottom in this flow.
        return self.itf_stack.conductors[-1].name

    def _create_beol_stack_from_itf(self, itf_data: ItfStack) -> BeolStack:
        layers: list[BeolLayer] = []

        dielectrics = {d.name: d for d in itf_data.dielectrics}
        vias_from = {v.from_layer: v for v in itf_data.vias}

        # ITF conductors are listed from highest to lowest.
        for cond in itf_data.conductors:
            diel_name = f"{cond.name}_diel"
            if diel_name in dielectrics:
                d = dielectrics[diel_name]
                layers.append(BeolLayer(name=d.name, type="oxide", thickness=d.thickness_nm))

            layers.append(
                BeolLayer(
                    name=cond.name,
                    type="metal",
                    thickness=cond.thickness_nm,
                    min_width=cond.min_width_nm,
                    pitch=cond.min_width_nm + cond.min_spacing_nm,
                    resistance_per_square=cond.resistance_per_square,
                )
            )

            if cond.name in vias_from:
                v = vias_from[cond.name]
                dielectric_thickness = 0.0
                if len(layers) > 1 and layers[-2].type == "oxide":
                    dielectric_thickness = layers[-2].thickness

                layers.append(
                    BeolLayer(
                        name=v.name,
                        type="via",
                        thickness=dielectric_thickness,
                        min_width=(v.area_nm2) ** 0.5,
                    )
                )

        layers.append(BeolLayer(name="substrate", type="substrate", thickness=1.0))
        return BeolStack(technology=itf_data.technology_name, layers=layers)

    def get_material(self, name: str) -> dict[str, float]:
        # Placeholder material model for capacitance calculations.
        if name == "oxide":
            return {"relative_permittivity": 2.5}
        return {}

    def get_beol_layer(self, name: str) -> BeolLayer:
        for layer in self.beol_stack.layers:
            if layer.name == name:
                return layer
        raise ValueError(f"BEOL layer '{name}' not found")

    def grid_layer_order(self) -> list[str]:
        """Return grid layer names in top-to-bottom order."""
        beol_names = [l.name for l in self.beol_stack.layers if l.type == "metal"]
        grid_layer_names = self.grid.layer_usage.keys()
        return [name for name in beol_names if name in grid_layer_names]

    def distance_to_nm(self, value: float) -> float:
        """Convert configured distance units to nanometers."""
        unit = self.units.distance.strip().lower()
        if unit == "nm":
            return value
        if unit == "um":
            return value * 1000.0
        if unit == "mm":
            return value * 1_000_000.0
        raise ValueError(f"Unsupported units.distance '{self.units.distance}'. Supported: nm, um, mm")

    @model_validator(mode="after")
    def validate_semantics(self) -> Config:
        # Trigger ITF load once so semantic checks can reference stack content.
        _ = self.itf_stack

        valid_metal_names = {c.name for c in self.itf_stack.conductors}
        configured_layers = set(self.grid.layer_usage.keys())
        invalid_layers = sorted(configured_layers - valid_metal_names)
        if invalid_layers:
            raise ValueError(
                f"grid.layer_usage contains layer(s) not present in ITF conductors: {', '.join(invalid_layers)}"
            )

        lowest_layer = self.lowest_metal_layer_name
        if lowest_layer in configured_layers:
            raise ValueError(
                "grid.layer_usage must not include the implicit lowest layer "
                f"'{lowest_layer}'. Remove it from YAML."
            )

        for layer_name, usage in self.grid.layer_usage.items():
            if usage.width <= 0 or usage.pitch <= 0:
                raise ValueError(
                    f"grid.layer_usage.{layer_name} must have width > 0 and pitch > 0"
                )
            conductor = next((c for c in self.itf_stack.conductors if c.name == layer_name), None)
            if conductor is None:
                continue
            if conductor.min_width_nm <= 0:
                raise ValueError(
                    f"ITF conductor {layer_name} has non-positive WMIN; "
                    "cannot apply layer_usage width multiplier"
                )
            if (conductor.min_width_nm + conductor.min_spacing_nm) <= 0:
                raise ValueError(
                    f"ITF conductor {layer_name} has non-positive (WMIN+SMIN); "
                    "cannot apply layer_usage pitch multiplier"
                )

        if self.grid.size.get("rows", 0) < 1 or self.grid.size.get("sites", 0) < 1:
            raise ValueError("grid.size.rows and grid.size.sites must both be >= 1")

        chains_cfg = self.spice_netlist.instance_chains
        if "max_instance_count_per_chain" not in chains_cfg and "max_count" in chains_cfg:
            # Backward-compatible alias support.
            chains_cfg["max_instance_count_per_chain"] = chains_cfg["max_count"]

        max_per_chain = int(chains_cfg.get("max_instance_count_per_chain", 0))
        if max_per_chain < 1:
            raise ValueError("spice_netlist.instance_chains.max_instance_count_per_chain must be >= 1")

        win = self.spice_netlist.ir_drop_measurement.get("averaging_window", {})
        start = float(win.get("start", -1))
        end = float(win.get("end", -1))
        if not (0.0 <= start < end <= 100.0):
            raise ValueError(
                "spice_netlist.ir_drop_measurement.averaging_window must satisfy 0 <= start < end <= 100"
            )

        if self.standard_cell_placement.row_height <= 0 or self.standard_cell_placement.site_width <= 0:
            raise ValueError("standard_cell_placement.row_height and site_width must be > 0")

        if self.standard_cell_placement.min_space.get("x", -1) < 0 or self.standard_cell_placement.min_space.get("y", -1) < 0:
            raise ValueError("standard_cell_placement.min_space.x and min_space.y must be >= 0")

        if self.itf.units.strip().lower() not in {"nm", "um", "mm"}:
            raise ValueError("itf.units must be one of: nm, um, mm")
        if self.units.capacitance.strip().lower() not in {"f", "pf", "ff", "nf", "uf"}:
            raise ValueError("units.capacitance must be one of: f, pf, ff, nf, uf")

        if self.ploc.pitch <= 0:
            raise ValueError("ploc.pitch must be > 0")
        if self.ploc.visualizer_render_diameter <= 0:
            raise ValueError("ploc.visualizer_render_diameter must be > 0")
        if self.ploc.offset_from_origin.get("x", 0) < 0 or self.ploc.offset_from_origin.get("y", 0) < 0:
            raise ValueError("ploc.offset_from_origin.x and y must be >= 0")

        if self.beol_thickness <= 0:
            raise ValueError("beol_thickness must be > 0")

        pi_segments = int(self.spice_netlist.standard_cell_output_load.get("number_pi_segments", 0))
        if pi_segments < 1:
            raise ValueError("spice_netlist.standard_cell_output_load.number_pi_segments must be >= 1")

        for cell in self.standard_cells:
            for pin in cell.pins:
                if pin.type == "signal" and pin.direction is None:
                    raise ValueError(
                        f"standard_cells[{cell.name}].pins[{pin.name}] must set direction for signal pins"
                    )

            pin_names = [p.name for p in cell.pins]
            port_names = cell.spice_port_order.split()
            if len(port_names) != len(pin_names) or set(port_names) != set(pin_names):
                raise ValueError(
                    f"standard_cells[{cell.name}].spice_port_order must list each declared pin exactly once"
                )

            if not any(p.type == "power" for p in cell.pins):
                raise ValueError(
                    f"standard_cells[{cell.name}] must include at least one pin with type=power"
                )
            if not any(p.type == "ground" for p in cell.pins):
                raise ValueError(
                    f"standard_cells[{cell.name}] must include at least one pin with type=ground"
                )

        if self.visualizer and self.visualizer.initial_visible_objects is not None:
            bad_entries = [v for v in self.visualizer.initial_visible_objects if not str(v).strip()]
            if bad_entries:
                raise ValueError("visualizer.initial_visible_objects entries must be non-empty strings")

        # Connectivity is checked for explicitly configured routing layers.
        # The implicit lowest ITF metal is not required in grid.layer_usage.
        active_layers_top_to_bottom = [
            c.name
            for c in self.itf_stack.conductors
            if c.name in configured_layers
        ]
        via_pairs = {(v.from_layer, v.to_layer) for v in self.itf_stack.vias}
        via_pairs |= {(b, a) for (a, b) in via_pairs}
        for i in range(len(active_layers_top_to_bottom) - 1):
            upper = active_layers_top_to_bottom[i]
            lower = active_layers_top_to_bottom[i + 1]
            if (upper, lower) not in via_pairs:
                raise ValueError(
                    f"Configured metal stack is not vertically connectable: missing ITF via between {upper} and {lower}"
                )

        return self


def load_config(path: str | Path) -> Config:
    """Load and validate a YAML configuration file."""
    path = Path(path)
    with open(path) as f:
        raw = yaml.safe_load(f)

    # Resolve ITF path relative to the config location before model validation.
    itf_file = Path(raw["itf"]["file"])
    if not itf_file.is_absolute():
        raw["itf"]["file"] = str((path.parent / itf_file).resolve())

    config = Config.model_validate(raw)
    config._config_path = path
    return config
