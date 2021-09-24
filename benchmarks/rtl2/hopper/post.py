#!/usr/bin/env python

from pathlib import Path

import matplotlib.pyplot as plt
import numpy

from rtl2.post import read_avg_values, read_historic_np, read_velocity

NUM_PARTICLES = 12359
VEL_P_FNAME = "uio_vel_p_0.dat"


def plot(refdata: Path) -> Path:
    running_avg = refdata
    hist_fname = refdata.parent / "mfix.R2016.1.gran.wdf.102219.dat"
    avg_dat = refdata.parent / "runningave.dat"

    refdata = read_historic_np(Path(hist_fname))
    tzero = t_zero()
    tz = "{:16.8}".format(tzero)

    assert avg_dat.is_file()
    with open(avg_dat, "a") as run_avg:
        run_avg.write(tz)
        run_avg.write("\n")
    avg_vals = read_avg_values(avg_dat)  # includes latest value that was just appended

    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{bm}')

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.0, 3.0), dpi=150)
    fig.tight_layout()
    # fig.suptitle('Horizontally stacked subplots')

    np_data = read_velocity(Path(VEL_P_FNAME))
    curr_times = [datum['t'] for datum in np_data]
    curr_values = [datum['np'] for datum in np_data]
    hist_times = [datum['t'] for datum in refdata]
    hist_values = [NUM_PARTICLES * datum['normalized_num_of_particles'] for datum in refdata]

    ax1.plot(hist_times, hist_values, color='k', linewidth=0.5)
    ax1.plot(curr_times, curr_values, color='r', linewidth=0.5)
    # ax1.set_title("Inventory (Np) compared to historic values")
    ax1.set_xlabel(r'$t(s)$')
    ax1.set_ylabel(r'$N_p$ (\#)')

    avg_x = list(range(-len(avg_vals), 0))
    avg_arr = numpy.asarray(avg_vals)
    avg_masked = numpy.ma.masked_where(avg_arr > 0.2, avg_arr)

    ax2.plot([-10, 0], [0.117, 0.117], color='k', linewidth=0.5)
    ax2.plot(list(range(-10, 0)), avg_masked[-10:], color='r', linewidth=0.5, marker='.')
    ax2.set_title(r'Recent $T_0$ over time')
    # ax2.set_xlabel("Previous runs")

    ax3.plot([-len(avg_vals), 0], [0.117, 0.117], color='k', linewidth=0.5)
    ax3.plot(avg_x, avg_masked, color='r', linewidth=0.5, marker='.')
    ax3.set_title(r'All $T_0$ over time')
    # ax3.set_xlabel("Previous runs")

    output_fname = "output.png"
    fig.savefig(output_fname)
    return Path(output_fname)


def t_zero() -> float:
    data = []
    with open(VEL_P_FNAME) as f:
        next(f)  # skip header
        for line in f:
            t, np, up, vp, wp, ke = line.split()
            data.append(
                {
                    "t": float(t),
                    "np": float(np),
                    "up": float(up),
                    "vp": float(vp),
                    "wp": float(wp),
                    "ke": float(ke),
                }
            )

    for ii in range(1, len(data)):
        if data[ii]["np"] == 0:
            return data[ii]["t"]
    return 0.21
