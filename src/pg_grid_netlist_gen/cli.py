"""Rich-Click CLI for pg_grid_netlist_gen."""

from __future__ import annotations

from pathlib import Path

import rich_click as click

click.rich_click.USE_RICH_MARKUP = True


@click.command()
@click.argument("config_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default=Path("output"),
    show_default=True,
    help="Output directory for generated files.",
)
@click.option("--netlist/--no-netlist", default=True, show_default=True, help="Generate SPICE netlist.")
@click.option("--viz/--no-viz", default=True, show_default=True, help="Generate 2D HTML visualization.")
@click.option(
    "--viz-region",
    type=str,
    default=None,
    show_default="None",
    help="Render subregion in microns: X1,Y1,X2,Y2",
)
@click.option("--open-browser", is_flag=True, default=False, show_default="False", help="Auto-open HTML after generation.")
@click.option("--save-image", type=str, default=None, show_default="None", help="Save a specific layer as a static PNG image.")
@click.option("--report/--no-report", default=True, show_default=True, help="Generate and print an ASCII summary report.")
def generate(
    config_file: Path,
    output_dir: Path,
    netlist: bool,
    viz: bool,
    viz_region: str | None,
    open_browser: bool,
    save_image: str | None,
    report: bool,
) -> None:
    """Generate power grid netlist and visualization from CONFIG_FILE."""
    from pg_grid_netlist_gen.config import load_config
    from pg_grid_netlist_gen.grid_builder import build_grid

    click.echo(f"Loading config: {config_file}")
    config = load_config(config_file)

    click.echo("Building grid...")
    grid = build_grid(config)
    output_dir.mkdir(parents=True, exist_ok=True)

    if report:
        from pg_grid_netlist_gen.reporter import generate_report

        summary_text = generate_report(grid, config)
        click.echo(summary_text)
        summary_path = output_dir / "pg_grid_summary.txt"
        summary_path.write_text(summary_text + "\n")
        click.echo(f"Writing summary report: {summary_path}")

    if netlist:
        from pg_grid_netlist_gen.netlist import write_netlist

        netlist_path = output_dir / "pg_grid_netlist.sp"
        click.echo(f"Writing netlist: {netlist_path}")
        write_netlist(grid, config, netlist_path)

    if viz or save_image:
        from pg_grid_netlist_gen.visualize import render_grid

        viz_path = output_dir / "pg_grid_visualization.html" if viz else None
        region = None
        if viz_region:
            parts = [float(x) for x in viz_region.split(",")]
            if len(parts) != 4:
                raise click.BadParameter("viz-region must be X1,Y1,X2,Y2")
            region = (parts[0], parts[1], parts[2], parts[3])

        if viz:
            click.echo(f"Rendering visualization: {viz_path}")
        
        render_grid(
            grid, config, 
            output_path=viz_path,
            viz_region=region,
            open_browser=open_browser,
            save_image_layer=save_image,
            output_dir=output_dir,
        )

    click.echo("Done!")


# Keep the public CLI symbol name unchanged for __main__/entry points.
app = generate
