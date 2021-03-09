import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd
import os
import glob
import sys


legend_names = [
    "SLS",
    "SPS (Random)",
    "SPS (Ours)",
]

markers = [
    "s",
    "D",
    "o",
]

marker_scales = [
    4,
    4,
    5,
]

def visualize_value_together(problem_dir):
    csv_names = [
        "sls-value.csv",
        "sps-fixed-center-free-cross-vecs-value.csv",
        "sps-first-ei-then-plane-integral-value.csv",
    ]

    sns.set()
    sns.set_context()

    plt.rcParams['font.sans-serif'] = ["Linux Biolinum", "Linux Biolinum O"]

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(1, 1, 1)

    xlim_min = 1
    xlim_max = 20

    ylim_min = 0.0
    ylim_max = None

    optimal_value = 0.0 if "rosenbrock" in problem_dir else 1.0

    file_paths = list(map(lambda name: problem_dir + "/" + name, csv_names))

    for index, file_path in enumerate(file_paths):
        print(file_path)

        df = pd.read_csv(file_path)

        ax.fill_between(df["#Iters"], optimal_value - df["Upper"], optimal_value - df["Lower"], alpha=0.2, cmap="inferno")
        ax.plot(df["#Iters"], optimal_value - df["Mean"], marker=markers[index], ms=marker_scales[index])

    ax.set_xlim([xlim_min, xlim_max])
    ax.set_ylim([ylim_min, ylim_max])

    ax.set_xticks(np.arange(xlim_min, xlim_max + 1, 1))

    ax.legend(legend_names)

    if False:
        ax.set_title("Optimality Gap")
    ax.set_xlabel("#Iterations")
    ax.set_ylabel("Optimality gap")

    fig.tight_layout()

    plt.savefig(problem_dir + "/value.pdf")

def visualize_residuals_together(problem_dir):
    csv_names = [
        "sls-residual.csv",
        "sps-first-ei-then-plane-integral-residual.csv",
        "sps-fixed-center-free-cross-vecs-residual.csv",
    ]

    sns.set()
    sns.set_context("paper")

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    xlim_min = 1
    xlim_max = 20

    ax.set_title("Residual")

    file_paths = list(map(lambda name: problem_dir + "/" + name, csv_names))

    for index, file_path in enumerate(file_paths):
        print(file_path)

        df = pd.read_csv(file_path)

        ax.fill_between(df["#Iters"], df["Upper"], df["Lower"], alpha=0.2)
        ax.plot(df["#Iters"], df["Mean"])

    ax.set_xlim([xlim_min, xlim_max])
    # ax.set_ylim([ylim_min, ylim_max])

    ax.set_xticks(np.arange(xlim_min, xlim_max + 1, 1))

    fig.tight_layout()

    plt.savefig(problem_dir + "/residuals.pdf")

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
    else:
        plt.savefig(file_name + "-trials.png")
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
    else:
        plt.savefig(file_name + ".png")
        plt.savefig(file_name + ".pdf")

    plt.close()


if __name__ == "__main__":
    assert len(sys.argv) >= 2

    data_dir = sys.argv[1]

    problem_dirs = glob.glob(data_dir + "/*-*d")

    for problem_dir in problem_dirs:
        visualize_residuals_together(problem_dir)
        visualize_value_together(problem_dir)
