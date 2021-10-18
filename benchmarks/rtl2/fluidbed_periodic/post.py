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

PG_1_FNAME = "uio_p_g_1.dat"
PG_2_FNAME = "uio_p_g_2.dat"
VEL_P_FNAME = "uio_vel_p_0.dat"

DENSITY = 1000.0
GRAVITY = 9.81
PARTICLE_VOLUME = math.pi * 0.0001 ** 3
CROSS_AREA = 6.0 * 0.002 ** 2
WEIGHT_OVER_AREA_PER_PARTICLE = DENSITY * PARTICLE_VOLUME * GRAVITY / CROSS_AREA


def plot(refdata: Path) -> Path:
    """Generate plot image file with matplotlib
    Returns: Path of the newly created image plot file"""

    running_avg = refdata
    refdata_dat = refdata.parent / "mfix.R2016.1.wdf.011619.dat"

    hist_vals = read_historic_pressure(refdata_dat)
    pg_data = read_two_pressures(Path(PG_1_FNAME), Path(PG_2_FNAME))
    vel_data = read_velocity(Path(VEL_P_FNAME))

    append_avg_dp_value(pg_data, vel_data, running_avg, WEIGHT_OVER_AREA_PER_PARTICLE)
    x_avg_vals, y_avg_vals = read_avg_values(running_avg)  # includes latest value that was just appended

    pressure_values = [(1.5 * datum['pg1'] - 0.5 * datum['pg2']) / 13.091 for datum in pg_data]
    pressure_times = [datum['t'] for datum in pg_data]
    hist_times = [datum['t'] for datum in hist_vals]
    hist_values = [datum['pg'] for datum in hist_vals]

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{bm}')

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.0, 3.0), dpi=150)
    fig.tight_layout()
    # fig.suptitle('Horizontally stacked subplots')

    ax1.plot(hist_times, hist_values, color='k', linewidth=0.5)
    ax1.plot(pressure_times, pressure_values, color='r', linewidth=0.5)
    # ax1.set_title("Pressure compared to historicvalues")
    ax1.set_xlabel(r'$t(s)$')
    ax1.set_ylabel(r'${dp}^\ast$')

    ax1.set_xlim(0.0, 1.2)
    ax1.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    ax1.set_ylim(bottom=0)

    ax2.plot([-10, 0], [0.99761206, 0.99761206], color='k', linewidth=0.5)
    ax2.plot(x_avg_vals[-10:], y_avg_vals[-10:], color='r', linewidth=0.5, marker='.')
    ax2.set_title(r'Recent $\bar{dp}^\ast$ over time')
    # ax2.set_xlabel("Previous runs")

    ax3.plot([-len(y_avg_vals), 0], [0.99761206, 0.99761206], color='k', linewidth=0.5)
    ax3.plot(x_avg_vals, y_avg_vals, color='r', linewidth=0.5, marker='.')
    ax3.set_title(r'All $\bar{dp}^\ast$ over time')
    # ax3.set_xlabel("Previous runs")

    output_fname = "output.png"
    fig.savefig(output_fname)
    return Path(output_fname)
