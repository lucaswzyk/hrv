# Wavelet Transform Modulus Maxima Method
# adjust result figure size so you get a clear picture
#
# author: lucaswzyk

import math

import analysis.analysis_settings as settings
import analysis.reader as rdr
import display.log_fit as lf
import display.disp_wtmm as disp


# constant settings for analysis, easy to change here
# for value description check analysis_settings.py
SETTINGS = settings.AnalysisSettings(
    min_beat=0,
    max_beat=24000,
    min_beat_resolution=1000,
    min_scale=2000,
    max_scale=70000,
    min_scale_resolution=100,
    exp_scales=True,
    q_range=5,
    min_q_resolution=40,
    file='../material/healthy_adult_25k.txt'
)

# will be set later on (could probably be done more nicely)
data = None
wt = None
wt_root = None


# returns fourth derivative of Gaussian function at value x
# most commonly used wavelet function, so we will use this as well
def third_gaussian(x):
    return - math.exp(-math.pow(x, 2) / 2) * x * (math.pow(x, 2) - 3) / math.sqrt(2*math.pi)


# returns wavelet data and min and max relevant beat
def gen_wavelet_vals(t0, scale):
    wavelet = []
    beat_min = 0
    beat_max = data.num_rr - 1

    # find first relevant beat (if) as beat_min, meaning first x coordinate with meaningful wavelet value
    # calculate relevant wavelet values (else)
    # break when wavelet becomes irrelevantly small again (elif) and saves last relevant beat as beat_max
    for i in range(0, data.num_rr):
        x = (data.rrs_cumul[i] - t0) / scale

        if x < -SETTINGS.wavelet_bound:
            beat_min = beat_min + 1
        elif x > SETTINGS.wavelet_bound:
            beat_max = i-1
            break
        else:
            wavelet.append(third_gaussian(x))

    return wavelet, beat_min, beat_max


# multiplies wavelet with signal and integrates relevant part
# uses trapezoid formula for integration
def wav_integration(wavelet, min_relev_beat, max_relev_beat, scale):
    integral = 0

    for j in range(min_relev_beat, max_relev_beat):
        # time between current data points (one RR interval)
        dx = data.rrs[j + 1]
        # product of series and wavelet at left beat
        left_val = data.rrs[j] * wavelet[j - min_relev_beat]
        # product of series and wavelet at right beat
        right_val = data.rrs[j + 1] * wavelet[(j + 1) - min_relev_beat]
        # trapezoid formula for integration
        integral += dx * (left_val + right_val)

    return abs(integral)/scale, math.pow(abs(integral)/scale, 1/2)


# returns fluctuation function Z_q(s) for specific q and scale
def mod_max(scale, q):
    wt_line = wt[scale]
    z_of_q = 0

    for i in range(1, len(wt_line) - 1):
        val = abs(wt_line[i])
        if abs(wt_line[i-1]) < val <= abs(wt_line[i+1]):
            z_of_q += math.pow(val, q)

    return z_of_q
    # the following yields simpler results (interpretation?)
    # return math.pow(z_of_q, 1/q)


# returns scaling exponent tau for specific q
def tau(q):
    scales = SETTINGS.scales
    zs = [mod_max(scale, q) for scale in scales]
    fit, p = lf.best_loglog_fit(scales, zs)
    disp.plot_tau_fit(scales, zs, fit)

    return p[0]


# setup and main loop
def analyze_file(conf=SETTINGS):
    global data, wt, wt_root
    data = rdr.RRData(conf)

    wt = {}
    wt_root = {}
    for scale in conf.scales:
        wt[scale] = []
        wt_root[scale] = []
        for beat0 in conf.beats:
            # get wavelet values at time of beat0 with scale scale
            wavelet, min_relev_beat, max_relev_beat = gen_wavelet_vals(data.rrs_cumul[beat0], scale)
            # calculate next result value and append to result wt
            integral = wav_integration(wavelet, min_relev_beat, max_relev_beat, scale)
            wt[scale].append(integral[0])
            wt_root[scale].append(integral[1])

        # print progress in %
        print("{:.0%}".format(math.pow((scale - conf.min_scale) / (conf.max_scale - conf.min_scale), 2)))

    disp.plot_tau_over_q(conf.qs, [tau(q) for q in conf.qs])
    disp.plot_wav_transform(wt_root, conf)


analyze_file()
