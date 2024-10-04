import pandas as pd
import numpy as np
from scipy.stats import percentileofscore
import matplotlib.pyplot as plt


def threshold_grabber(thresh_vals, outcome_vals, true_labels=None):
    # Label outcome as High or Low according to median by default
    if true_labels is None:
        outcome_threshold = outcome_vals.median()
    true_labels = outcome_vals.apply(lambda x: "High" if x >= outcome_threshold else "Low")

    # Combine values into df for easy sorting
    df = pd.DataFrame({"thresh_vals": thresh_vals, "outcome_vals": outcome_vals, "true_labels": true_labels}).sort_values("thresh_vals").reset_index(drop=True)

    # Get threshold before, after, and in between each value by list compin'
    thresholds = sorted(list(set(
        [x - 1 if idx == 0 else (x + df["thresh_vals"][idx - 1]) / 2
        for idx, x in 
        enumerate([*df["thresh_vals"], max(df["thresh_vals"]) + 1])]
        )))

    # Calculate False Positive Rate and True Positive Rate at each threhold
    fpr_list = []
    tpr_list = []
    for t in thresholds:
        fpr = len(df[(df["thresh_vals"] >= t) & (df["true_labels"] == "Low")]) / len(df[df["true_labels"] == "Low"])
        tpr = len(df[(df["thresh_vals"] >= t) & (df["true_labels"] == "High")]) / len(df[df["true_labels"] == "High"])
        fpr_list.append(fpr)
        tpr_list.append(tpr)

    # Calculate auc from tpr and fpr
    auc = abs(round(np.trapz(y=tpr_list, x=fpr_list), 3))

    # Calculate Youden Index
    youden_idx = [tpr - fpr for tpr, fpr in zip(tpr_list, fpr_list)]

    # Calculate quantiles of dataset for each threshold
    quantiles = [percentileofscore(df["thresh_vals"], thresh, kind="rank") for thresh in thresholds]

    thresh_df = pd.DataFrame({"fpr": fpr_list, "tpr": tpr_list, "auc": auc, "Youden_idx": youden_idx, "threshold": thresholds, "quantile": quantiles})

    return thresh_df


def roc_plotter(thresh_vals=None, outcome_vals=None, name=None, true_labels=None, colour="blue", youden=False):

    thresh_df = threshold_grabber(thresh_vals, outcome_vals, true_labels)

    fpr_list = thresh_df["fpr"]
    tpr_list = thresh_df["tpr"]
    auc = thresh_df["auc"][0]

    plt.clf()
    plt.cla()

    plt.rcParams["figure.dpi"] = 600

    # Reference line
    plt.plot([0, 1], [0, 1], color="gainsboro", linestyle="--", alpha=0.5, zorder=1)

    plt.scatter(fpr_list, tpr_list, marker="^", color=colour).set_clip_on(False)

    plt.plot(fpr_list, tpr_list, linestyle="-.", color=colour, alpha=0.5)

    if youden:
        youden_idx = thresh_df["Youden_idx"].max()
        youden_fpr = thresh_df[thresh_df["Youden_idx"] == youden_idx]["fpr"].reset_index(drop=True)
        youden_tpr = thresh_df[thresh_df["Youden_idx"] == youden_idx]["tpr"].reset_index(drop=True)
        youden_thresh = thresh_df[thresh_df["Youden_idx"] == youden_idx]["threshold"].reset_index(drop=True)
        plt.text(x=0.60, y=0.15, s=f"Youden Index: {round(youden_idx, 3)}", va="center")
        plt.scatter(0.55, 0.15, marker="^", color="black").set_clip_on(False)

        for idx, thresh in enumerate(youden_thresh):
            # plt.axvline(x=youden_fpr[idx], ymin=0, ymax=youden_tpr[idx], color="black", linestyle="-")
            # plt.plot([youden_fpr[idx], (youden_fpr[idx]+youden_tpr[idx])/2], [youden_tpr[idx],(youden_fpr[idx]+youden_tpr[idx])/2], c="black", zorder=1)
            # plt.text(x=0.7, y=0.1, s=f"Youden Index: {round(youden_idx, 3)}")
            plt.scatter(youden_fpr[idx], youden_tpr[idx], marker="^", color="black", zorder=100).set_clip_on(False)

    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")

    plt.xlim([0, 1])
    plt.ylim([0, 1])

    plt.text(x=0.60, y=0.2, s=f"AUC = {auc}")

    # Graphpad-ify
    plt.rcParams["font.family"] = ["Arial"]
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_linewidth(2)
    plt.gca().spines["bottom"].set_linewidth(2)
    plt.gca().tick_params(width=2)
    plt.gca().set_aspect("equal")

    plt.savefig(f"../Data/ROC_{name}_07.2024.png", dpi=600, transparent=True)
