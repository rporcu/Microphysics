#!/usr/bin/env python

import math
from pathlib import Path

import matplotlib.pyplot as plt

from rtl2.post import (
    append_avg_dp_value,
    read_avg_values,
    read_historic_pressure,
    read_two_pressures,
    read_velocity,
)


def plot(refdata: Path) -> Path:
    """Generate plot image file with matplotlib
    Returns: Path of the newly created image plot file"""

    running_avg = refdata
    x_avg_vals, y_avg_vals = read_avg_values(
        running_avg
    )  # includes latest value that was just appended

    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath} \usepackage{bm}")

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.0, 3.0), dpi=150)
    fig.tight_layout()

    ax1.set_xlabel(r"$t(s)$")
    ax1.set_ylabel(r"${dp}^\ast$")

    ax1.set_xlim(0.0, 1.2)
    ax1.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    ax1.set_ylim(bottom=0)

    ax2.plot([-10, 0], [1.0028038, 1.0028038], color="k", linewidth=0.5)
    ax2.plot([-10, 0], [0.99936198, 0.99936198], color="k", linewidth=0.5)
    ax2.plot(x_avg_vals[-10:], y_avg_vals[-10:], color="r", linewidth=0.5, marker=".")
    ax2.set_title(r"Recent $\bar{dp}^\ast$ over time")
    # ax2.set_xlabel("Previous runs")

    ax3.plot([-len(y_avg_vals), 0], [1.0028038, 1.0028038], color="k", linewidth=0.5)
    ax3.plot([-len(y_avg_vals), 0], [0.99936198, 0.99936198], color="k", linewidth=0.5)
    ax3.plot(x_avg_vals, y_avg_vals, color="r", linewidth=0.5, marker=".")
    ax3.set_title(r"All ${dp}^\ast$ over time")

    output_fname = "output.png"
    fig.savefig(output_fname)

    return Path(output_fname)
