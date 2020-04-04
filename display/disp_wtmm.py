# display methods for WTMM
#
# author: lucaswzyk

import numpy as np
import matplotlib.pyplot as plt


# plots tau(q) over q
def plot_tau_over_q(qs, taus):
    plt.plot(qs, taus)
    plt.show()


# plots loglog fit from which tau is derived directly into section of Z_q(s) over s
def plot_tau_fit(scales, zs, fit):
    plt.loglog(scales, zs, scales, fit)
    plt.show()


# plots wavelet transformation result with height as color
def plot_wav_transform(wt_root, settings):
    # coordinate sets for scatterplot
    xs = list(settings.beats) * len(settings.scales)
    ys = [[i]*len(settings.beats) for i in settings.scales]
    ys = list(np.ravel(ys))
    all_vals = list(wt_root.values())
    all_vals = list(np.ravel(all_vals))
    m = max(all_vals)

    # colors generated from wavelet transformation result as [r, g, b] with g limited to 1
    colors = [[0, val / m, 0] for val in all_vals]

    ax = plt.gca()
    ax.scatter(xs, ys, c=colors, s=.5, alpha=1)
    ax.set_yscale('log')
    # ax.colorbar()
    plt.show()

