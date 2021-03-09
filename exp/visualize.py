import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd
import os
import glob
import sys


def visualize_trials(file_name):
    df = pd.read_csv(file_name)

    sns.set()
    sns.set_context() # "paper"

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1)

    xlim_min = df["#Iters"][0]
    xlim_max = df["#Iters"][df["#Iters"].shape[0] - 1]
    ylim_min = 0.00 if "residual" in file_name else None
    ylim_max = None # magic number { 0.50, 0.70, 0.90 }

    ax.set_title(os.path.basename(file_name))

    num_trials = df.shape[1] - 5
    for i in range(num_trials):
        ax.plot(df["#Iters"], df["Trial #" + str(i + 1)])

    ax.set_xlim([xlim_min, xlim_max])
    ax.set_ylim([ylim_min, ylim_max])

    ax.set_xticks(np.arange(xlim_min, xlim_max + 1, 1))

    fig.tight_layout()

    if False:
        plt.show()
    elif False:
        plt.savefig(file_name + "-trials.png")
        plt.savefig(file_name + "-trials.pdf")
    else:
        plt.savefig(file_name + "-trials.pdf")

    plt.close()


def visualize_stats(file_name):
    df = pd.read_csv(file_name)

    sns.set()
    sns.set_context() # "paper"

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1)

    xlim_min = df["#Iters"][0]
    xlim_max = df["#Iters"][df["#Iters"].shape[0] - 1]
    ylim_min = 0.00 if "residual" in file_name else None
    ylim_max = None # magic number { 0.50, 0.70, 0.90 }

    ax.set_title(os.path.basename(file_name))

    ax.fill_between(df["#Iters"], df["Upper"], df["Lower"], alpha=0.2)
    ax.plot(df["#Iters"], df["Mean"])

    ax.set_xlim([xlim_min, xlim_max])
    ax.set_ylim([ylim_min, ylim_max])

    ax.set_xticks(np.arange(xlim_min, xlim_max + 1, 1))

    fig.tight_layout()

    if False:
        plt.show()
    elif False:
        plt.savefig(file_name + ".png")
        plt.savefig(file_name + ".pdf")
    else:
        plt.savefig(file_name + ".pdf")

    plt.close()


if __name__ == "__main__":
    assert len(sys.argv) >= 2

    data_dir = sys.argv[1]

    csvs = glob.glob(data_dir + "/**/*.csv")

    for csv in csvs:
        visualize_stats(csv)
        visualize_trials(csv)
