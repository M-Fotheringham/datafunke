"""Plots data generated from histogrammer or any x, y format."""

import matplotlib.pyplot as plt


def hist_plotter(
    df=None,
    x=None,
    y=None,
    sample=None,
    pheno=None,
    limits=None,
    colour=None,
    r=1000,
    prob=True,
    save=True,
    panel=None,
    expr=False,
    dpi=900,
):
    """
    df: provide dataframe when plotting multiple phenotypes
    x, y: array-like or names of x and y cols when providing a df
    sample: only for labelling and saving
    pheno: for labelling and subsetting plot data
    limits: tuple. Defaults to full x domain
    colour: optional dict for colour-coding phenotypes
    r: magic variable for plotting text (will edit aesthetics later)
    prob: True. Formats y-axis as probability
    save: True. Sometimes you don't want to save/want to tweak the plot
    panel: a temp variable for determining which default color dict to use
    expr: combines lineage and expr phenotype labels to plot separately
    """

    # Preprocess inputs ====================================================
    if isinstance(pheno, str):
        pheno = [pheno]

    # Default colour dict based off wsi02 when multiple phenos provided
    # To do: add color col to df during query to automatically colour-code
    if (colour is None) and (df is not None):
        if panel is None:
            colour = {
                "CD8": "yellow",
                "FoxP3": "red",
                "FoxP3CD8": "lightsteelblue",
                "CD163": "magenta",
                "Tumor": "orange",
                "Other": "blue",
                "CD79a": "lime",
                "ERG": "yellow",
                "CD3": "red",
            }
        elif panel == "tbet":
            colour = {
                "Gzmb": "lime",
                "Lag3": "yellow",
                "Eomes": "red",
                "PD1": "orange",
                "Tbet": "lightblue",
                "CD8": "purple",
            }
    elif colour is None:
        colour = "lightsteelblue"
    # =======================================================================

    # Graphpad-ify ==========================================================
    # Print then clear fig to apply aesthetic changes
    # Change to not affect rcparams
    plt.rcParams["figure.dpi"] = dpi
    plt.rcParams["font.family"] = ["Arial"]
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    # plt.gca().spines["top"].set_visible(False)
    # plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_linewidth(1.75)
    plt.gca().spines["right"].set_linewidth(1.75)
    plt.gca().spines["left"].set_linewidth(1.75)
    plt.gca().spines["bottom"].set_linewidth(1.75)
    plt.gca().tick_params(width=1.75)
    # =======================================================================

    # Plot the data =========================================================
    if df is not None:
        df = df.copy()
        if expr:  # combine lineage and expr to filter/plot as needed
            df["phenotype"] = (
                df["phenotype"] + "_" + df["exprphenotype"].astype(str)
            )
        for p in pheno:
            d = df[df["phenotype"] == p]
            plt.plot(d[x], d[y], color=colour[p], label=p)
    else:
        plt.plot(x, y, color=colour)
    # =======================================================================

    # Additional plot formatting ============================================
    # Get height of plot from y ticks
    # first and last ticks are beyond plot range
    max_y_tick = plt.yticks()[0][-1]
    # For plotting tumor boundary text
    pen_y_tick = plt.yticks()[0][-2]

    plt.xticks(color="dimgrey", fontweight="normal")
    plt.yticks(color="dimgrey", fontweight="normal")

    # Set limits
    if limits is not None:
        if not isinstance(limits, tuple):
            print(
                f"""Limits ({limits}) not provided as tuple.
                Defaulting to full x domain."""
            )
        else:
            plt.xlim(limits[0], limits[1])
    else:
        plt.ylim(0, max_y_tick)

    # Tumour boundary line
    plt.axvline(x=0, color="green", linewidth=1.5)

    # Modify to make dynamic
    # Tumor/Stroma arrows =================================================
    # plt.gca().annotate(
    #     "",
    #     xy=(-1.8 * r, max_y_tick * 1.025),
    #     xytext=(0, max_y_tick * 1.025),
    #     arrowprops=dict(arrowstyle="-|>", color="black"),
    #     annotation_clip=False,
    # )
    # plt.gca().annotate(
    #     "",
    #     xy=(1.8 * r, max_y_tick * 1.025),
    #     xytext=(0, max_y_tick * 1.025),
    #     arrowprops=dict(arrowstyle="-|>", color="black"),
    #     annotation_clip=False,
    # )

    # # Tumor/stroma labels
    # plt.text(
    #     -r,
    #     max_y_tick * 1.05,
    #     "Tumor",
    #     fontsize=12,
    #     fontweight="normal",
    #     ha="center",
    #     va="center",
    # )
    # plt.text(
    #     r,
    #     max_y_tick * 1.05,
    #     "Stroma",
    #     fontsize=12,
    #     fontweight="normal",
    #     ha="center",
    #     va="center",
    # )
    # plt.text(-120, pen_y_tick, "Tumor", rotation=90, ha="center", va="center")
    # plt.text(140, pen_y_tick, "Boundary", rotation=90, ha="center", va="center")
    # ========================================================================

    # rvt = clin[clin["sampleid"] == int(sample)]["rvt"].to_list()[0]
    # \n(RVT: {rvt})
    # \nPhenotype: {', '.join(map(str, pheno))}"
    plt.title(
        f"Sample: {sample}",
        fontweight="bold",
        fontsize=16,
        pad=50,
    )
    plt.xlabel("Distance from Tumor Boundary (µm)", labelpad=15)
    if prob:
        plt.ylabel("Probability", labelpad=15)
    else:
        plt.ylabel(f"Density (cells/mm\u00b2)")

    # labels, colours = zip(*colour.items())
    plt.legend(
        bbox_to_anchor=(1.0, 0.75), loc="center left", title="Cell Phenotype"
    )
    # =======================================================================

    if save:
        plt.savefig(
            f"../Data/{sample}_{pheno[0]}_dist.png",
            dpi=dpi,
            bbox_inches="tight",
        )
