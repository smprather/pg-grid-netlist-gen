"""Parser for Interconnect Technology Format (ITF) files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re


@dataclass
class ItfConductor:
    name: str
    thickness_nm: float
    min_width_nm: float
    min_spacing_nm: float
    resistance_per_square: float


@dataclass
class ItfDielectric:
    name: str
    thickness_nm: float
    relative_permittivity: float


@dataclass
class ItfVia:
    name: str
    from_layer: str
    to_layer: str
    resistance_per_via: float
    area_nm2: float

@dataclass
class ItfStack:
    technology_name: str
    conductors: list[ItfConductor] = field(default_factory=list)
    dielectrics: list[ItfDielectric] = field(default_factory=list)
    vias: list[ItfVia] = field(default_factory=list)


def _parse_itf_line(line: str) -> dict[str, str]:
    """Parses a braced key-value section of an ITF line."""
    props = {}
    # Use regex to find all KEY=VALUE pairs
    for match in re.finditer(r"(\w+)\s*=\s*([\w.-]+)", line):
        props[match.group(1)] = match.group(2)
    return props


def parse_itf(path: str | Path, distance_unit: str = "um") -> ItfStack:
    """Parses an ITF file and returns a structured representation."""
    path = Path(path)
    stack = ItfStack(technology_name="Unknown")
    unit_key = distance_unit.strip().lower()
    unit_to_nm = {"nm": 1.0, "um": 1000.0, "mm": 1_000_000.0}
    if unit_key not in unit_to_nm:
        raise ValueError(f"Unsupported ITF unit '{distance_unit}'. Supported: nm, um, mm")
    scale = unit_to_nm[unit_key]
    
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("$"):
                continue

            parts = line.split()
            if not parts:
                continue
            
            keyword = parts[0]
            
            if keyword == "TECHNOLOGY":
                if len(parts) >= 3 and parts[1] == "=":
                    stack.technology_name = parts[2]
            elif keyword == "CONDUCTOR":
                name = parts[1]
                props = _parse_itf_line(line)
                stack.conductors.append(ItfConductor(
                    name=name,
                    thickness_nm=float(props.get("THICKNESS", 0.0)) * scale,
                    min_width_nm=float(props.get("WMIN", 0.0)) * scale,
                    min_spacing_nm=float(props.get("SMIN", 0.0)) * scale,
                    resistance_per_square=float(props.get("RPSQ", 0.0)),
                ))
            elif keyword == "DIELECTRIC":
                name = parts[1]
                props = _parse_itf_line(line)
                stack.dielectrics.append(ItfDielectric(
                    name=name,
                    thickness_nm=float(props.get("THICKNESS", 0.0)) * scale,
                    relative_permittivity=float(props.get("ER", 1.0)),
                ))
            elif keyword == "VIA":
                name = parts[1]
                props = _parse_itf_line(line)
                stack.vias.append(ItfVia(
                    name=name,
                    from_layer=props.get("FROM", ""),
                    to_layer=props.get("TO", ""),
                    resistance_per_via=float(props.get("RPV", 0.0)),
                    area_nm2=float(props.get("AREA", 0.0)) * (scale * scale),
                ))

    return stack
