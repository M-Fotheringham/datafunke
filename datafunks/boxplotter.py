import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def boxplotter(df, ylabel="Density (cells/mm\u00b2)"):

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

    plt.boxplot(df)

    plt.ylabel(ylabel)
