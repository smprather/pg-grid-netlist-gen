"""Physics calculations for resistance and capacitance."""

from __future__ import annotations

import math

# Vacuum permittivity (F/m)
EPSILON_0 = 8.854187817e-12


def plate_capacitance(
    width_nm: float,
    length_nm: float,
    distance_to_substrate_nm: float,
    relative_permittivity: float,
) -> float:
    """Calculate plate-to-substrate capacitance.

    C = ε₀ × εᵣ × W × L / d
    """
    width_m = width_nm * 1e-9
    length_m = length_nm * 1e-9
    d_m = distance_to_substrate_nm * 1e-9
    if d_m <= 0:
        return 0.0
    return EPSILON_0 * relative_permittivity * width_m * length_m / d_m


def fringe_capacitance(
    length_nm: float,
    thickness_nm: float,
    distance_to_substrate_nm: float,
    relative_permittivity: float,
) -> float:
    """Calculate fringing capacitance using Sakurai-Tamaru model.

    C_fringe = ε₀ × εᵣ × L × (2/π) × ln(1 + 2T/d + 2√(T/d × (1 + T/d)))

    Both edges included.
    """
    length_m = length_nm * 1e-9
    thickness_m = thickness_nm * 1e-9
    d_m = distance_to_substrate_nm * 1e-9
    if d_m <= 0:
        return 0.0
    t_over_d = thickness_m / d_m
    arg = 1.0 + 2.0 * t_over_d + 2.0 * math.sqrt(t_over_d * (1.0 + t_over_d))
    c_per_edge = (
        EPSILON_0 * relative_permittivity * length_m * (2.0 / math.pi) * math.log(arg)
    )
    return 2.0 * c_per_edge  # both edges


def total_segment_capacitance(
    width_nm: float,
    length_nm: float,
    thickness_nm: float,
    distance_to_substrate_nm: float,
    relative_permittivity: float,
) -> tuple[float, float, float]:
    """Calculate total capacitance for a segment.

    Returns (cap_plate, cap_fringe, cap_total).
    """
    cp = plate_capacitance(
        width_nm, length_nm, distance_to_substrate_nm, relative_permittivity
    )
    cf = fringe_capacitance(
        length_nm, thickness_nm, distance_to_substrate_nm, relative_permittivity
    )
    return cp, cf, cp + cf
