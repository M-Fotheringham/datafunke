import pandas as pd
from astropathdb import AstroDB
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
import colorcet as cc
import holoviews as hv
from holoviews.element.tiles import EsriImagery
from holoviews.operation.datashader import rasterize

# , datashade
# import datashader as ds

hv.extension("bokeh")


def plot_cells(
    pheno=None,
    sample=None,
    database="wsi02",
    colour="red",
    x=None,
    y=None,
    mpl=False,
    geo=False,
    shader=False,
    save_label=False,
):

    if x is None and y is None:
        pheno = "CD8"
        sample = 125
        print(f"No coords provided, querying {database}: {sample} {pheno}.")

        db = AstroDB(database=database)

        sql = f"""
        select ct.px, ct.py
        from dbo.celltag ct, dbo.phenotype p
        where ct.ptype = p.ptype
        and p.phenotype = '{pheno}'
        and ct.sampleid = {sample}
        """
        cell_coords = db.query(sql)

        x = cell_coords["px"]
        y = cell_coords["py"]

    # Optional annotation overlay - to do

    # Create color_dict for phenotypes - to do

    if not any([mpl, geo, shader]):
        print(
            """Select a plotting method as True from: mpl, geo, or shader.
            Defaulting to mpl."""
        )
        mpl = True

    if mpl:
        plt.rcParams["font.family"] = ["Arial"]
        plt.rcParams["font.weight"] = "bold"

        plt.scatter(x, y, c=colour, s=1, edgecolor="black", linewidth=0.1)

        plt.gca().set_aspect("equal")
        plt.gca().set_axis_off()

        if save_label:
            plt.savefig(
                f"../Data/{save_label}_sampleid_{sample}_{pheno}.png",
                dpi=1200,
                transparent=True,
            )

    elif geo:
        plt.rcParams["font.family"] = ["Arial"]
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["figure.dpi"] = 1200
        # Alternatively, geopandas df for speed
        cells_plot = gpd.GeoSeries([Point(i, j) for i, j in zip(x, y)])
        cells_plot.plot(
            ax=plt.gca(),
            color=colour,
            edgecolor="black",
            markersize=1,
            linewidth=0.1,
            aspect="equal",
        )
        if save_label:
            plt.savefig(
                f"../Data/{save_label}_sampleid_{sample}_{pheno}.png",
                dpi=1200,
                transparent=True,
            )

    elif shader:
        coord_df = pd.DataFrame({"x": x, "y": y})

        map_tiles = EsriImagery().opts(
            alpha=0.5, width=900, height=480, bgcolor="black"
        )
        points = hv.Points(coord_df, ["x", "y"])
        # cells_plot = datashade(points, x_sampling=50, y_sampling=50,
        # cmap=cc.fire, width=900, height=480)
        ropts = dict(
            tools=["hover"],
            colorbar=True,
            colorbar_position="bottom",
            cmap=cc.fire,
            cnorm="linear",
        )
        cells_plot = rasterize(points, x_sampling=200, y_sampling=200).opts(
            **ropts
        )

        # cvs = ds.Canvas(plot_width=900, plot_height=480)
        # agg = cvs.points(coord_df, "x", "y")
        # return ds.tf.set_background(ds.tf.shade(agg, cmap=cc.fire), "black")

        plotted = map_tiles * cells_plot

        if save_label:
            hv.save(
                plotted, f"../Data/{save_label}_sampleid_{sample}_{pheno}.html"
            )

        return plotted
