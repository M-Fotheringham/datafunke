"""Consistently create forest plot for any number of covariates."""

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def cox_plot(
    cohort=None,
    variable_cols=None,
    outcome_col=None,
    event_col=None,
    x_lim=5,
    y_min=-0.5,
    multivariable=True,
    cohort_label="None",
):
    if not isinstance(variable_cols, list) and variable_cols is not None:
        variable_cols = [variable_cols]

    if multivariable:
        multi = "multivariable_"
        # Drop irrelevant columns and any Nan vals
        data = cohort[
            ~cohort[[*variable_cols, outcome_col, event_col]].isna().any(axis=1)
        ][[*variable_cols, outcome_col, event_col]]
        cph = CoxPHFitter()
        cph.fit(df=data, duration_col=outcome_col, event_col=event_col)

        res = cph.summary

        cph.plot(hazard_ratios=True)

    else:
        multi = ""
        cap_height = 1 / 48
        # Perform CoxPH per variable
        res_list = []
        y_labels = []
        for idx, var in enumerate(variable_cols):
            # Drop irrelevant columns and any Nan vals
            data = cohort[~cohort[[var, outcome_col, event_col]].isna().any(axis=1)][
                [var, outcome_col, event_col]
            ]
            cph = CoxPHFitter()
            cph.fit(df=data, duration_col=outcome_col, event_col=event_col)

            res = cph.summary

            res_list.append(res)

            y_labels.append(var)

            plt.scatter(
                x=res["exp(coef)"].iloc[0],
                y=idx,
                marker="s",
                color="white",
                edgecolor="black",
                zorder=100,
            )

            plt.axvline(
                x=1, linestyle="--", color="black", alpha=0.90, linewidth=0.5, zorder=0
            )

            plt.axhline(
                y=idx,
                xmin=res["exp(coef) lower 95%"].iloc[0] / x_lim,
                xmax=res["exp(coef) upper 95%"].iloc[0] / x_lim,
                color="black",
                zorder=1,
            )

            # Plot end caps for CI bars
            plt.axvline(
                x=res["exp(coef) lower 95%"].iloc[0],
                ymin=(idx - y_min - cap_height) / len(variable_cols),
                ymax=(idx - y_min + cap_height) / len(variable_cols),
                color="black",
                zorder=1,
            )
            plt.axvline(
                x=res["exp(coef) upper 95%"].iloc[0],
                ymin=(idx - y_min - cap_height) / len(variable_cols),
                ymax=(idx - y_min + cap_height) / len(variable_cols),
                color="black",
                zorder=1,
            )

        res = pd.concat(res_list)

        plt.gca().set_yticks([*range(len(variable_cols))], y_labels)
        plt.xlabel("HR (95% CI)")

    res_plot = res[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    res_plot.columns = ["HR", "Lower 95% CI", "Upper 95% CI", "p-value"]
    res_plot = res_plot.round(3)
    res_plot.insert(loc=0, column="Covariate", value=res_plot.index)
    # Reverse order to match plotted HRs
    res_plot = res_plot.iloc[::-1].reset_index(drop=True)

    # minHR = res["exp(coef) lower 95%"].min()
    # maxHR = res["exp(coef) upper 95%"].max()

    # print(f"""Min HR C: {minHR}, Max HR C:{maxHR}.""")

    # Graphpadify
    plt.rcParams["figure.dpi"] = 900
    plt.rcParams["font.family"] = ["Arial"]
    # plt.rcParams["axes.labelweight"] = "bold"
    # plt.rcParams["font.weight"] = "bold"
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.xlim((0, x_lim))
    plt.ylim((y_min, y_min + len(variable_cols)))
    plt.title("CoxPH")

    # Plot legend
    tab = plt.table(
        cellText=res_plot.values,
        colLabels=res_plot.columns,
        cellLoc="left",
        loc="right",
        edges="horizontal",
    )
    tab.auto_set_column_width(col=list(range(len(res_plot.columns))))
    tab.scale(1, 1.25)
    # Bold column titles a la Matplotlib documentation
    for (row, col), cell in tab.get_celld().items():
        if (row == 0) or (col == -1):
            cell.set_text_props(fontproperties=FontProperties(weight="bold"))

    plt.savefig(
        f"""../Data/{cohort_label}_{outcome_col}_{multi}08.2024.png""",
        dpi=900,
        transparent=True,
        bbox_inches="tight",
    )

    return res_plot
