# container for settings used in analysis algorithms
#
# author: lucaswzyk

import math
import numpy as np


# Direct input values
# min_beat, max_beat      first/last beat to consider
# min_beat_resolution     min resolution for beat axis, actual value may be greater due to rounding or smaller
#                         due to insufficient sample size (max_beat - min_beat)
# min_scale, max_scale    smallest/largest applied wavelet scale
# min_scale_resolution    min resolution for scale axis (actual value may be greater due to rounding)
# exp_scales              boolean whether to use exponentially growing scales instead of linearly growing scales
# q_range                 boundary for relevant values of q (DFA), lower bound is -q_range
# min_q_resolution        min resolution of q values (DFA)
# wavelet_bound           boundary for relevant wavelet support, lower bound is -wavelet_bound
# deg_pol                 degree of polynomial fit (DFA)
# file                    .txt file with RR interval data
# int_scales              boolean whether to round scales to integers (relevant for WTMM)

# Generated values
# beat_step               distance between considered beats
# beats                   range determined by beat parameters
# scale_step              distance between considered scales (can be evaluated exponentially)
# scales                  list of considered scales
# q_step                  distance between considered qs (DFA)
# qs                      list of considered qs (will not contain 0)
class AnalysisSettings:
    def __init__(self,
                 min_beat, max_beat, min_beat_resolution,
                 min_scale, max_scale, min_scale_resolution, exp_scales,
                 min_q_resolution, q_range=20,
                 wavelet_bound=6,
                 deg_pol=2,
                 file='../material/healthy_adult_25k.txt',
                 int_scales=False):
        self.min_beat = min_beat
        self.max_beat = max_beat
        self.min_beat_resolution = min(min_beat_resolution, max_beat - min_beat)
        self.beat_step = max(math.floor((max_beat - min_beat) / min_beat_resolution), 1)
        self.beats = range(0, (max_beat - min_beat), self.beat_step)

        if min_scale == 0:
            print('WARNING: min_scale = 0')
        self.min_scale = min_scale
        self.max_scale = max_scale
        self.min_scale_resolution = min(min_scale_resolution, max_scale - min_scale)
        self.scale_step = float(max_scale - min_scale) / self.min_scale_resolution
        if exp_scales:
            self.exp_step = (math.log(max_scale, min_scale) - 1) / (min_scale_resolution - 1)
            if int_scales:
                self.scales = [int(math.pow(min_scale, 1 + e * self.exp_step)) for e in range(min_scale_resolution)]
            else:
                self.scales = [math.pow(min_scale, 1 + e * self.exp_step) for e in range(min_scale_resolution)]
        else:
            self.scales = np.linspace(min_scale, max_scale, num=self.min_scale_resolution)

        self.min_q = -q_range
        self.max_q = q_range
        self.min_q_resolution = min_q_resolution
        self.q_step = 2 * q_range / self.min_q_resolution
        qs = list(np.linspace(self.min_q, self.max_q, num=self.min_q_resolution+1))
        self.qs = qs
        # uncomment the following for quadratic q values normed by q_range
        # self.qs = [q*abs(q) / q_range for q in qs]
        if 0.0 in self.qs:
            self.qs.remove(0.0)

        self.wavelet_bound = wavelet_bound
        self.deg_pol = deg_pol
        self.file = file
