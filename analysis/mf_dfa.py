# Multifractal Detrended Fluctuation Analysis
# calculates DFA results and plots them with the help of display package
#
# author: lucaswzyk
import math
import numpy as np
import matplotlib.pyplot as plt

import analysis.analysis_settings as settings
import analysis.reader as rdr
import display.disp_mf_dfa as disp
import display.log_fit as lf

# constant settings for analysis, easy to change here
# for value description check analysis_settings.py
STD_SETTINGS = settings.AnalysisSettings(
    min_beat=0,
    max_beat=24462,
    min_beat_resolution=24463,
    min_scale=10,
    max_scale=5000,
    min_scale_resolution=70,
    exp_scales=True,
    q_range=5,
    min_q_resolution=20,
    deg_pol=2,
    file='../material/healthy_adult_25k.txt',
    int_scales=True
)

# will be set later on (could probably be done more nicely)
sett = None
data = None


def analyze_file(conf=STD_SETTINGS):
    global data, sett
    sett = conf
    data = rdr.RRData(conf)

    # calculate F^2(nu, s) as two-dimensional list
    f2s = f2()
    # calculate F_q
    f = fq(f2s)
    disp.plt_fq_over_s(f, sett.qs)
    disp.plt_h_over_q([h(f[q]) for q in sett.qs], sett.qs)


# returns F^2(nu, s) as two-dimensional list
def f2():
    f2_vs = {}

    for scale in sett.scales:
        n_s = n_scale(scale)
        if n_s == 0:
            print('Scale larger than time series: ', scale)
            break

        f2_vs[scale] = []

        for nu in range(0, data.num_rr - scale, scale):
            f2_vs[scale].append(var_segment(nu, nu + scale))
        for nu in range(data.num_rr - 1, scale - 1, -scale):
            f2_vs[scale].append(var_segment(nu - scale, nu))

    return f2_vs


# returns the variance in a given segment, F^2 for specific nu and s
def var_segment(min_beat, max_beat):
    scale = max_beat - min_beat
    rrs_cumul_loc = data.rrs_cumul[min_beat:max_beat]
    p = np.poly1d(np.polyfit(range(0, scale), rrs_cumul_loc, sett.deg_pol))
    detrended = sum([math.pow(rrs_cumul_loc[i] - p(i), 2) for i in range(0, scale)])

    return detrended / scale


# returns the number of segments that fit into the sample using segment size scale
def n_scale(scale):
    return int(data.num_rr / scale)


# returns F_q as a map with scales as keys
def fq(f2s):
    f = {}

    for q in sett.qs:
        f[q] = {}
        for scale in sett.scales:
            f2s_scale = f2s[scale]
            f[q][scale] = math.pow(
                            sum(
                                map(lambda x: 0 if x == 0 else math.pow(x, q / 2), f2s_scale)
                            ) / (2 * n_scale(scale)),
                          1 / q)
    return f


# returns local hurst exponent h(q) from F_q
def h(f_of_q):
    fit, p = lf.best_loglog_fit(sett.scales, [f_of_q[scale] for scale in sett.scales])
    return p[0]


analyze_file()
