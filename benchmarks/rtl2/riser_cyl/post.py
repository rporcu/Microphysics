from pathlib import Path

import matplotlib.pyplot as plt

from rtl2 import suite
from rtl2.test_util import git_commit
from datetime import datetime

from rtl2.post import (
    get_mean_vp,
    read_avg_values,
    read_historic_velocity,
    read_velocity,
)

VEL_P_FNAME = Path("uio_vel_p_0.dat")


def plot(refdata: Path) -> Path:
    """Generate plot image file with matplotlib
    Returns: Path of the newly created image plot file"""

    running_avg = refdata
    stl_16 = refdata.parent / "mfix.R2016.1.stl16.wdf.011419.dat"
    stl_24 = refdata.parent / "mfix.R2016.1.stl24.wdf.011519.dat"

    hist_16 = read_historic_velocity(stl_16, 0.01)
    hist_24 = read_historic_velocity(stl_24, 0.01)
    meanvp = get_mean_vp(VEL_P_FNAME)
    mean_vp = "{:16.8}".format(meanvp)

    # append new mean_vp value to running_val file
    assert running_avg.is_file()
    if not suite.get_suite().post_only:
        time = datetime.now().astimezone().isoformat()
        branch, sha = git_commit(suite.get_suite().sourceDir)
        with open(running_avg, "a") as run_avg:
            run_avg.write(f"{mean_vp}\t{time}\t{sha}\t{branch}\n")

    x_avg_vals, y_avg_vals = read_avg_values(running_avg)  # includes latest value that was just appended

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{bm}')

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.0, 3.0), dpi=150)
    fig.tight_layout()
    # fig.suptitle('Horizontally stacked subplots')

    vel_data = read_velocity(Path(VEL_P_FNAME))
    curr_times = [datum['t'] for datum in vel_data]
    curr_values = [datum['up'] for datum in vel_data]
    hist_1_times = [datum['t'] for datum in hist_16]
    hist_1_values = [datum['velocity'] for datum in hist_16]
    hist_2_times = [datum['t'] for datum in hist_24]
    hist_2_values = [datum['velocity'] for datum in hist_24]

    ax1.plot(hist_1_times, hist_1_values, color='k', linewidth=0.5)
    ax1.plot(hist_2_times, hist_2_values, color='k', linewidth=0.5)
    ax1.plot(curr_times, curr_values, color='r', linewidth=0.5)
    ax1.set_xlabel(r'$t(s)$')
    ax1.set_ylabel("mean $V_p$ (m/s)")

    ax1.set_xlim(0.0, 3.0)
    # ax1.set_xticks([0.,0.2,0.4,0.6,0.8,1.0,1.2])

    ax2.plot(x_avg_vals[-10:], y_avg_vals[-10:], color='r', linewidth=0.5, marker='.')
    ax2.set_title(r'Recent $\bar{V_p}$ over time')
    # ax2.set_xlabel("Previous runs")

    ax3.plot(x_avg_vals, y_avg_vals, color='r', linewidth=0.5, marker='.')
    ax3.set_title(r'All $\bar{V_p}$ over time')
    # ax3.set_xlabel("Previous runs")

    output_fname = "output.png"
    fig.savefig(output_fname)

    return Path(output_fname)
