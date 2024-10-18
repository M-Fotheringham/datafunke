"""Compute tdist histogram with user-defined bin widths."""

import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from .get_area import get_area
from .get_cell_counts import get_cell_counts


def tdistogram(
    sampleid=None,
    database=None,
    phenos=None,
    tdist_filter=None,
    rdist_filter=None,
    all_reg=False,
    excl_ln=False,
    t_hist_step=50,
    prob=False,
    samplewise=True,
    smoothed=True,
    altdb=None,
    prop=False,
    t_hist_type=None,
):
    """Builds histogram of cell counts by tdist."""

    print("Counting cells...")
    cells = get_cell_counts(
        sampleid,
        database,
        phenos,
        tdist_filter,
        rdist_filter,
        all_reg,
        excl_ln,
        t_hist_step,
        altdb=altdb,
        t_hist_type=t_hist_type,
    )

    # Cells are assigned a dist of 32700+ when annotation absent, exclude
    # cells = cells[cells["dist_bin_um"] != 16350]

    print("Counting area...")
    areas = get_area(
        sampleid,
        database,
        phenos,
        tdist_filter,
        rdist_filter,
        all_reg,
        excl_ln,
        t_hist_step,
        altdb=altdb,
        t_hist_type=t_hist_type,
    )

    data = pd.merge(cells, areas, on=["sampleid", "tdist_bin"], how="left")

    # Drop rows without area (extra rows can be added during get_cell_counts)
    data = data[~data["area_mm"].isna()]

    if not samplewise:
        data = (
            data.groupby(["phenotype", "exprphenotype", "tdist_bin"])[
                ["c", "area_mm"]
            ]
            .sum()
            .reset_index()
        )

    data["density_mm"] = data["c"] / data["area_mm"]

    # Exclude bins with 0.4mm^2 or less area
    data = data[data["area_mm"] > 0.4**2].reset_index(drop=True)

    if prop:
        if samplewise:
            # Calculate proportion of each phenotype per sample
            data["grouped_c"] = data.groupby(["sampleid", "tdist_bin"])[
                "c"
            ].transform(lambda x: x[data["exprphenotype"] != "Total"].sum())
            data["prop"] = data["c"] / data["grouped_c"]
        else:
            # Calculate proportion of each phenotype per sample
            data["grouped_c"] = data.groupby("tdist_bin")["c"].transform(
                lambda x: x[data["exprphenotype"] != "Total"].sum()
            )
            data["prop"] = data["c"] / data["grouped_c"]

    if prob:

        if samplewise:
            # Calculate overall density per expr, per pheno, per sample
            grouped_den = data.groupby(
                ["sampleid", "phenotype", "exprphenotype"]
            )[["c", "area_mm"]].sum()

            grouped_den["grouped_den"] = (
                grouped_den["c"] / grouped_den["area_mm"]
            )

            grouped_den = grouped_den["grouped_den"].reset_index()

            data = pd.merge(
                data,
                grouped_den,
                on=["sampleid", "phenotype", "exprphenotype"],
                how="left",
            )

            data["prob"] = data["density_mm"] / data["grouped_den"]

            # Normalize before smoothing
            data["prob"] = data["prob"] / data["prob"].sum()

            if smoothed:
                grouped = data.groupby(
                    ["sampleid", "phenotype", "exprphenotype"]
                )

                smoothed = []
                for name, group in grouped:
                    tdist = group["tdist_bin"].values
                    density = group["prob"].values

                    loess_result = lowess(density, tdist, frac=0.15)

                    xl = np.linspace(min(tdist), max(tdist), 1000)

                    smoothed_y = np.interp(
                        xl, loess_result[:, 0], loess_result[:, 1]
                    )

                    smoothed_df = pd.DataFrame(
                        {"tdist_smoothed": xl, "smoothed_y": smoothed_y}
                    )

                    smoothed_df["tdist_bin"] = (
                        np.floor(smoothed_df["tdist_smoothed"] / t_hist_step)
                        * t_hist_step
                    ).astype(int)

                    smoothed_df = pd.merge(
                        group, smoothed_df, on="tdist_bin", how="right"
                    )

                    smoothed.append(smoothed_df)

                data = pd.concat(smoothed, ignore_index=True)

        else:
            # Calculate overall density per expr, per pheno, NOT per sample
            grouped_den = data.groupby(["phenotype", "exprphenotype"])[
                ["c", "area_mm"]
            ].sum()

            grouped_den["grouped_den"] = (
                grouped_den["c"] / grouped_den["area_mm"]
            )

            grouped_den = grouped_den["grouped_den"].reset_index()

            data = pd.merge(
                data, grouped_den, on=["phenotype", "exprphenotype"], how="left"
            )

            data["prob"] = data["density_mm"] / data["grouped_den"]

            # Normalize before smoothing
            data["prob"] = data["prob"] / data["prob"].sum()

            if smoothed:
                grouped = data.groupby(["phenotype", "exprphenotype"])

                smoothed = []
                for name, group in grouped:
                    tdist = group["tdist_bin"].values
                    density = group["prob"].values

                    loess_result = lowess(density, tdist, frac=0.15)

                    xl = np.linspace(min(tdist), max(tdist), 1000)

                    smoothed_y = np.interp(
                        xl, loess_result[:, 0], loess_result[:, 1]
                    )

                    smoothed_df = pd.DataFrame(
                        {"tdist_smoothed": xl, "smoothed_y": smoothed_y}
                    )

                    smoothed_df["tdist_bin"] = (
                        np.floor(smoothed_df["tdist_smoothed"] / t_hist_step)
                        * t_hist_step
                    ).astype(int)

                    smoothed_df = pd.merge(
                        group, smoothed_df, on="tdist_bin", how="right"
                    )

                    smoothed.append(smoothed_df)

                data = pd.concat(smoothed, ignore_index=True)

    # Reformat, modularize ====================================================
    else:
        if smoothed:
            if samplewise:
                grouped = data.groupby(
                    ["sampleid", "phenotype", "exprphenotype"]
                )

                smoothed = []
                for name, group in grouped:
                    tdist = group["tdist_bin"].values

                    if prop:
                        y_vals = group["prop"].values
                    else:
                        y_vals = group["density_mm"].values

                    loess_result = lowess(y_vals, tdist, frac=0.15)

                    xl = np.linspace(min(tdist), max(tdist), 1000)

                    smoothed_y = np.interp(
                        xl, loess_result[:, 0], loess_result[:, 1]
                    )

                    smoothed_df = pd.DataFrame(
                        {"tdist_smoothed": xl, "smoothed_y": smoothed_y}
                    )

                    smoothed_df["tdist_bin"] = (
                        np.floor(smoothed_df["tdist_smoothed"] / t_hist_step)
                        * t_hist_step
                    ).astype(int)

                    smoothed_df = pd.merge(
                        group, smoothed_df, on="tdist_bin", how="right"
                    )

                    smoothed.append(smoothed_df)

                data = pd.concat(smoothed, ignore_index=True)

            else:
                grouped = data.groupby(["phenotype", "exprphenotype"])

                smoothed = []
                for name, group in grouped:
                    tdist = group["tdist_bin"].values

                    if prop:
                        y_vals = group["prop"].values
                    else:
                        y_vals = group["density_mm"].values

                    loess_result = lowess(y_vals, tdist, frac=0.15)

                    xl = np.linspace(min(tdist), max(tdist), 1000)

                    smoothed_y = np.interp(
                        xl, loess_result[:, 0], loess_result[:, 1]
                    )

                    smoothed_df = pd.DataFrame(
                        {"tdist_smoothed": xl, "smoothed_y": smoothed_y}
                    )

                    smoothed_df["tdist_bin"] = (
                        np.floor(smoothed_df["tdist_smoothed"] / t_hist_step)
                        * t_hist_step
                    ).astype(int)

                    smoothed_df = pd.merge(
                        group, smoothed_df, on="tdist_bin", how="right"
                    )

                    smoothed.append(smoothed_df)

                data = pd.concat(smoothed, ignore_index=True)
    # ========================================================================

    return data
