# AGENTS.md

## Project Overview

This will be a script that generates a SPICE netlist of resistors, capacitors, and standard cell instances.
A yaml file (pg_grid_netlist_gen.yaml) will provide a description of the metal layer stackup, including vias, materials, and standard
cells.
The yaml file will also describe the grid configuration, including:

- Whether the layer is a "grid" layer, or a "staple" layer.
    - A "grid" layer contains stripes of metal that span the space of the grid
      either top-to-bottom for vertical layers, or side-to-side for horizontal layers.
    - A "staple" layer should just be a square shape of metal with sides 2x as thick as the min-width of the layer above it.
- The directions of the metal stripes on the layer, vertical or horizontal
- Standard cells should be connected to the bottom-most via layer
- The power and ground nets in the grid must match the standard cell power and ground pin names (this is a limitation that might be
  removed in the future).
- Any pin that is not power or ground can be connected to no-connect nets. This is temporary. The next phase in the poject will
  specify how to connect the inputs and outputs of the cells.
- Standard cells should be distributed according to the spacing definition between the rows of the bottom-most metal layer. Use
  your knowledge of how standard cells are typically placed, including a vertical-flip for every other row.

The output is:

1. A 3D rendering of the grid
    - Web browser display is fine. Or other. The agent is free to suggest alternatives.
2. A SPICE netlist model of the grid
    - Calculate the resistance between either vias, or standard cell tap-points (on the bottom layer of metal)
    - Calculate the capacitance of each resistance segment. Caclulate the plate-to-plate capacitance
      to the substrate layer. Also calculate fringing capacitance to the substrate layer.
      Do not calculate capacitance to other conductors. When calculating fringe capacitance, assume
      the segment of wire stands alone in the space.

The beol_stack/layers are listed in order from highest to lowest. Layers of type "via" are used to connect the "metal" layer above
and below the via. The via shape should be a square with x and y lengths of min_width. A via should be placed anywhere the layer
above and layer below intersect.

## Tech Stack

- Python 3.14
- Astral uv to manage the project
- Rich-Click for CLI
