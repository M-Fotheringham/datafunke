import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


def scatterPlotter(
    df,
    x,
    y,
    x_label,
    y_label,
    point_labels=False,
    ref_line=False,
    label="_",
    log2=False,
    zero_val=0.0000005,
):

    plt.cla()
    plt.clf()

    if df is not None:

        # Remove NaNs, print which samples are excluded
        dropped = df.iloc[x[np.isnan(x)].index, :]["sampleid"].values
        if len(dropped) > 0:
            print(
                f"The following {len(dropped)} samples are dropped because they are missing x values: {dropped}."
            )
            y = y[~np.isnan(x)]
            x = x[~np.isnan(x)]
    #         drop_list = [*np.where(np.isnan(x))[0]]
    #         x = x[~np.isnan(x)]
    #         y_drop = y.values
    #         for idx in drop_list:
    #             y_drop = np.delete(y_drop, idx)

    if log2:
        # Side-step log2(0) error by parsing as -inf and replacing with small value
        with np.errstate(divide="ignore"):
            x = np.log2(x)
            y = np.log2(y)
        x[np.isneginf(x)] = zero_val
        y[np.isneginf(y)] = zero_val
        log_label = "True"
    else:
        log_label = "False"

    if not isinstance(x_label, (tuple, list)) or not isinstance(y_label, (tuple, list)):
        print(
            "Please provide x, y labels as tuples: (descriptive label, label suffix), ie ('CD8', ' Density (cells/mm\u00b2))."
        )

    # Scatter plot
    plt.figure(dpi=600)
    plt.scatter(x, y, c="crimson", s=4)
    plt.xlabel(x_label[0] + x_label[1], fontweight="bold", fontsize=14)
    #  + " Density (structures/mm\u00b2)"
    plt.ylabel(y_label[0] + y_label[1], fontweight="bold", fontsize=14)

    if ref_line:
        # .reshape(-1, 1)
        slope, intercept, pearson_r, pearson_p, std_err = stats.linregress(x, y)
        # Data is not normally distributed, use spearmanr
        r, p = stats.spearmanr(x, y)
        if p < 0.001:
            p = "<0.001"
        else:
            p = round(p, 3)
        r = round(r, 3)
        r2 = round(r**2, 3)
        # model = LinearRegression(fit_intercept=False)
        # model.fit(x, y)
        # slope = model.coef_[0][0]
        # Spearman R\u00B2
        plt.axline(
            (0, intercept),
            slope=slope,
            c="gainsboro",
            linestyle="--",
            linewidth=1,
            alpha=0.5,
        )
        plt.text(
            max(x) - (max(x) - min(y)) * 0.4,
            max(y) - (max(y) - min(y)) * 0.2,
            f"Spearman coef: {r}\np: {p}",
            fontweight="bold",
            fontsize=14,
            ha="left",
            va="center",
        )
        # plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), c="gainsboro", linestyle="--", linewidth=1, alpha=0.5)

    if point_labels:
        # Get sampleids from dataframe using indices from individual arrays
        point_labs = df.iloc[x.index, :]["sampleid"].reset_index(drop=True).values
        # Add sampleids to points
        for idx, (xi, yi) in enumerate(zip(x, y)):
            # pending -- offset overlapping points
            plt.text(
                xi,
                yi,
                f"{point_labs[idx]}",
                fontweight="bold",
                fontsize=6,
                ha="left",
                va="top",
                alpha=0.5,
            )

    # Graphpad-ify
    plt.rcParams["font.family"] = ["Arial"]
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.gca().tick_params(width=2)
    # plt.gca().set_box_aspect(1)

    plt.tight_layout()

    plt.savefig(
        f"../Data/{label}_{x_label[0]}_{y_label[0]}_scatter_log2={log_label}_08.2024.png".replace(
            " ", "_"
        ),
        dpi=600,
        transparent=True,
    )
