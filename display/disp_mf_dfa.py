# display methods for DFA
#
# author: lucaswzyk

import matplotlib.pyplot as plt
import numpy as np
import display.log_fit as lf


# plots F_q(s) over s
def plt_fq_over_s(f, qs):
    plt.figure(1)
    legend = []

    for q in qs:
        scales = list(f[q].keys())
        fq_s = list(f[q].values())
        plt.loglog(scales, fq_s, 'o-', markersize=2)

        log_fit = lf.best_loglog_fit(scales, fq_s)[0]
        plt.loglog(scales, log_fit)
        legend.append('q = ' + str(q))
    plt.legend(legend)
    plt.show()


# plots original RR values, best fit and detrended values
def plt_detrending(rrs, detr, ps, ts):
    fig, (normal, pol, detrended) = plt.subplots(3)
    normal.plot(ts, rrs)
    pol.plot(ts, ps)
    detrended.plot(ts, detr)
    plt.show()


# plots h(q) over q
def plt_h_over_q(hqs, qs):
    plt.figure(2)

    new_qs = []
    tqs = []
    for i in range(len(qs)):
        if abs(qs[i]) <= 5:
            new_qs.append(qs[i])
            tqs.append(hqs[i])

    [m, k] = np.polyfit(new_qs, tqs, 1)
    lqs = np.array([m*q + k for q in new_qs])
    tqs = np.array(tqs)
    var = np.var(np.abs(tqs - lqs))
    print('var =', var)
    plt.plot(new_qs, lqs, "o-", markersize=2)
    plt.plot(new_qs, tqs, "o-", markersize=2)

    # plt.plot(qs, hqs, "o-")
    plt.show()
