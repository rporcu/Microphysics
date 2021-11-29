#!/usr/bin/env python

def process_velp(fname, vel_index):

    from math import sqrt
    from scipy import stats

    To = 5 # truncate time from startup
    Tb = 5 # bin averaging window size

    try:
        with open(fname) as f:
            lines = f.readlines()
    except:
        lines = []

    X = []
    S = []
    Te = To + Tb # bin end time
    bigT = 0.
    lX = 0.
    lS = 0.
    t0=0.
    for line in lines[1:]:
        t_str, np_str, u_str, v_str, w_str, ke_str = line.split()
        time = float(t_str)
        np = int(np_str)

        if vel_index == 0:
            x = float(u_str)
        elif vel_index == 1:
            x = float(v_str)
        elif vel_index == 2:
            x = float(w_str)
        else:
            raise NotImplementedError

        if time <= To:
            t0 = time
            continue

        dt = (time-t0)*0.5

        if time > Te:
            lX /= bigT
            lS /= bigT
            lS = sqrt(lS-(lX*lX))
            X.append(lX)
            S.append(lS)
            Te += Tb
            bigT = 0.
            lX=0.
            lS=0.

        if np > 0:
            bigT += dt
            lX += x*dt
            lS += x*x*dt

        t0 = time

    if bigT > 0.:
        lX /= bigT
        lS /= bigT

        if (lS-(lX*lX)) >=0:
            lS = sqrt(lS-(lX*lX))

            X.append(lX)
            S.append(lS)

    if(len(X) > 1):

        Fs = stats.t.ppf(1-0.05, len(X)-1)

        return {
            'mu_hatX':    stats.tmean(X),
            'sigma_hatX': stats.tstd(X)*Fs/len(X),
            'mu_hatS':    stats.tmean(S),
            'sigma_hatS': stats.tstd(S)*Fs/len(S),
        }

    else:

        return {
            'mu_hatX':    stats.tmean(X),
            'sigma_hatX': None,
            'mu_hatS':    stats.tmean(S),
            'sigma_hatS': None,
        }
