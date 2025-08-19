# python3 statistics.py </path/to/data/dir>

import numpy as np
import pandas as pd
import sys
import scipy.stats as st


file_names = [
    "sps-first-ei-then-plane-integral-value.csv",
    "sps-fixed-center-free-cross-vecs-value.csv",
    "sls-value.csv",
]

problem_dir_names = [
    "rosenbrock-5d",
    "rosenbrock-10d",
    "rosenbrock-15d",
    "rosenbrock-20d",
    "isometric-5d",
    "isometric-10d",
    "isometric-15d",
    "isometric-20d",
]

num_trials = 50
num_iters = 25

if __name__ == "__main__":
    assert len(sys.argv) >= 2

    data_dir = sys.argv[1]

    problem_dirs = list(map(lambda name: data_dir + "/" + name, problem_dir_names))

    for problem_dir in problem_dirs:
        print(problem_dir)

        file_paths = list(map(lambda name: problem_dir + "/" + name, file_names))

        raw_data_list = []
        for file_path in file_paths:
            data_frame = pd.read_csv(file_path)
            raw_data = data_frame.to_numpy()[:, 1 : num_trials + 1]  # hard coding
            raw_data_list.append(raw_data)

        print("-- ours vs. random -- ")

        for i in range(num_iters):
            stat, p_value = st.mannwhitneyu(
                raw_data_list[0][i, :], raw_data_list[1][i, :], alternative="two-sided"
            )
            print(
                "iter "
                + "{:2d}".format(i + 1)
                + ": "
                + ("*" if p_value < 0.05 else " ")
                + " "
                + ("*" if p_value < 0.05 / 3.0 else " ")
                + " (p = {:.6f}, U = {:.1f}, f = {:.6f})".format(
                    p_value, stat, stat / float(num_trials * num_trials)
                )
            )

        print("-- ours vs. line -- ")

        for i in range(num_iters):
            stat, p_value = st.mannwhitneyu(
                raw_data_list[0][i, :], raw_data_list[2][i, :], alternative="two-sided"
            )
            print(
                "iter "
                + "{:2d}".format(i + 1)
                + ": "
                + ("*" if p_value < 0.05 else " ")
                + " "
                + ("*" if p_value < 0.05 / 3.0 else " ")
                + " (p = {:.6f}, U = {:.1f}, f = {:.6f})".format(
                    p_value, stat, stat / float(num_trials * num_trials)
                )
            )

        print("-- random vs. line -- ")

        for i in range(num_iters):
            stat, p_value = st.mannwhitneyu(
                raw_data_list[1][i, :], raw_data_list[2][i, :], alternative="two-sided"
            )
            print(
                "iter "
                + "{:2d}".format(i + 1)
                + ": "
                + ("*" if p_value < 0.05 else " ")
                + " "
                + ("*" if p_value < 0.05 / 3.0 else " ")
                + " (p = {:.6f}, U = {:.1f}, f = {:.6f})".format(
                    p_value, stat, stat / float(num_trials * num_trials)
                )
            )
