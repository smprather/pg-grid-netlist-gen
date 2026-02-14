"""Rich-Click CLI for pg_grid_netlist_gen."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

click.rich_click.USE_RICH_MARKUP = True


@click.group()
def app() -> None:
    """Power/Ground Grid Netlist Generator."""


@app.command()
@click.argument("config_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=Path("output"),
    help="Output directory for generated files.",
)
@click.option("--netlist/--no-netlist", default=True, help="Generate SPICE netlist.")
@click.option("--viz/--no-viz", default=True, help="Generate 3D HTML visualization.")
@click.option(
    "--viz-region",
    type=str,
    default=None,
    help="Render subregion in microns: X1,Y1,X2,Y2",
)
@click.option(
    "--z-exaggeration",
    type=float,
    default=50.0,
    help="Z axis scale factor for visualization.",
)
@click.option("--open-browser", is_flag=True, help="Auto-open HTML after generation.")
def generate(
    config_file: Path,
    output_dir: Path,
    netlist: bool,
    viz: bool,
    viz_region: str | None,
    z_exaggeration: float,
    open_browser: bool,
) -> None:
    """Generate power grid netlist and visualization from CONFIG_FILE."""
    from pg_grid_netlist_gen.config import load_config
    from pg_grid_netlist_gen.grid_builder import build_grid

    click.echo(f"Loading config: {config_file}")
    config = load_config(config_file)

    click.echo("Building grid...")
    grid = build_grid(config)

    # Print stats
    click.echo(f"  Stripes:  {len(grid.stripes)}")
    click.echo(f"  Staples:  {len(grid.staples)}")
    click.echo(f"  Segments: {len(grid.segments)}")
    click.echo(f"  Vias:     {len(grid.vias)}")
    click.echo(f"  Cells:    {len(grid.cells)}")
    click.echo(f"  Nodes:    {len(grid.nodes)}")

    if netlist:
        from pg_grid_netlist_gen.netlist import write_netlist

        netlist_path = output_dir / "power_grid.sp"
        click.echo(f"Writing netlist: {netlist_path}")
        write_netlist(grid, config, netlist_path)

    if viz:
        from pg_grid_netlist_gen.visualize import render_grid

        viz_path = output_dir / "power_grid.html"
        region = None
        if viz_region:
            parts = [float(x) for x in viz_region.split(",")]
            if len(parts) != 4:
                raise click.BadParameter("viz-region must be X1,Y1,X2,Y2")
            region = (parts[0], parts[1], parts[2], parts[3])

        click.echo(f"Rendering visualization: {viz_path}")
        render_grid(
            grid, config, viz_path,
            z_exaggeration=z_exaggeration,
            viz_region=region,
            open_browser=open_browser,
        )

    click.echo("Done!")
