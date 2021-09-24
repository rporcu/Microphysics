from pathlib import Path
from typing import List, Mapping, Optional

import matplotlib.pyplot as plt

from rtl2.post import (
    # append_avg_meanvp_value,
    get_mean_vp,
    read_avg_values,
    read_historic_velocity,
    read_two_pressures,
    read_velocity,
)

VEL_P_FNAME = Path("uio_vel_p_0.dat")


def plot(refdata: Path) -> Path:
    """Generate plot image file with matplotlib
    Returns: Path of the newly created image plot file"""

    running_avg = refdata
    refdata_dat = refdata.parent / "mfix.R2016.1.wdf.011619.dat"

    refdata = read_historic_velocity(refdata_dat, 0.01)
    meanvp = get_mean_vp(VEL_P_FNAME)
    mean_vp = "{:16.8}".format(meanvp)

    # append new mean_vp value to running_val file
    assert running_avg.is_file()
    with open(running_avg, "a") as run_avg:
        run_avg.write(mean_vp)
        run_avg.write("\n")
    avg_vals = read_avg_values(running_avg)  # includes latest value that was just appended

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{bm}')

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.0, 3.0), dpi=150)
    fig.tight_layout()
    # fig.suptitle('Horizontally stacked subplots')

    vel_data = read_velocity(Path(VEL_P_FNAME))
    curr_times = [datum['t'] for datum in vel_data]
    curr_values = [datum['up'] for datum in vel_data]
    hist_times = [datum['t'] for datum in refdata]
    hist_values = [datum['velocity'] for datum in refdata]

    ax1.plot(hist_times, hist_values, color='k', linewidth=0.5)
    ax1.plot(curr_times, curr_values, color='r', linewidth=0.5)
    ax1.set_xlabel(r'$t(s)$')
    ax1.set_ylabel("mean $V_p$ (m/s)")

    ax1.set_xlim(0.0, 3.0)
    # ax1.set_xticks([0.,0.2,0.4,0.6,0.8,1.0,1.2])

    ax2.plot(list(range(-10, 0)), avg_vals[-10:], color='r', linewidth=0.5, marker='.')
    ax2.set_title(r'Recent $\bar{V_p}$ over time')
    # ax2.set_xlabel("Previous runs")

    ax3.plot(list(range(-len(avg_vals), 0)), avg_vals, color='r', linewidth=0.5, marker='.')
    ax3.set_title(r'All $\bar{V_p}$ over time')
    # ax3.set_xlabel("Previous runs")

    output_fname = "output.png"
    fig.savefig(output_fname)

    return Path(output_fname)
