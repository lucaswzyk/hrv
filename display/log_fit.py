# functions for working with semilog / loglog graphs
#
# author: lucaswzyk

import numpy as np
import math


# returns best polynomial fit of degree deg for log xs and log ys
def best_loglog_fit(xs, ys, deg=1):
    [m, k] = np.polyfit(np.log10(xs), np.log10(ys), deg)
    log_ys = [math.pow(10, m * np.log10(x) + k) for x in xs]

    return log_ys, [m, k]
