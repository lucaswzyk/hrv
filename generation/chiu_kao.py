# Chiu Kao Model (incl. Warner Cox Model)
#
# author: lucaswzyk

import numpy as np
from scipy.integrate import odeint
from matplotlib import rc, pyplot as plt

# Font formatting for output picture
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Constants
k1 = 6.356
k2 = 1.7#2.05
k3 = 2.25
k4 = 4
k5 = 29.13
k6 = 10# addition to hr0
k7 = 5.75
k8 = 5#.69
k9 = 7#3.87
k10 = 1# addition to rr0
k11 = 2.4
k12 = 1.6

a0 = 1
c1 = 1
v1 = 1
v2 = 10
nm = 1
ns = 1
nv = 1

starting_val = [.12, 1, 1, .88, 1, 0]#[.36, 4.1, 4.1, .64, .89, .16]
t_end = 280
t_res = 10*t_end
s_to_ms = 1000
t_step = s_to_ms / t_res

T = 1
m0 = 1.38


# comment one of the following function so it only returns 0
# to get a result as in 2.8
#
# f_s(t) an in 2.8
def fs(t):
    if t < 5:
        return 0
    elif t < 40:
        return 5
    elif t < 80:
        return 0
    elif t < 120:
        return 2
    elif t < 160:
        return 0
    elif t < 185:
        return 20
    elif t < 225:
        return 0
    elif t < 255:
        return 10
    return 0


# f_v(t) as in 2.8
def fv(t):
    if t < 7:
        return 0
    elif t < 30:
        return 5
    elif t < 65:
        return 0
    elif t < 85:
        return 10
    return 0


# right side of warner cox model
def warner_model(y, t):
    ab, a1, a2, b, n, c2 = y
    d_ab = k4 * a2 * b - k5 * ab
    d_a1 = (k1 * ns * fs(t) + k2 * (a0 - a1) + k3 * (a2 - a1)) / v1
    d_a2 = (k3 * (a1 - a2) - d_ab) / v2
    d_b = -d_ab
    d_n = k7 * (nm - n) - k8 * n * fv(t)
    d_c2 = (nv * k8 * n * c1 * fv(t) - k9 * c2) / v2
    return [d_ab, d_a1, d_a2, d_b, d_n, d_c2]


# Integral Pulse Frequency Modulation
# of linear combination of NE and ACh to get a discrete RR series
def ipfm(sol):
    xs = []
    hrs = []
    integral = t = 0

    for sample in sol:
        integral += (m0 + k11 * sample[0] - k12 * sample[5]) / t_res
        t += t_step
        if integral > T:
            xs.append(t)
            hrs.append(60 * s_to_ms / t)
            integral = t = 0

    return xs, hrs


# method for generating and displaying RR series based on f_s and f_v
def solve_and_display(u0, t):
    ts = [(i/t_res) for i in range(0, t_res*t)]
    sol = odeint(warner_model, u0, ts)

    ab, a1, a2, b, n, c2 = list(zip(*sol))
    fig, ([ne_model, ach_model, rr_plot], [ab_model, vesicle_model, hr_plot]) = plt.subplots(2, 3)
    # fig, (ach_model, vesicle_model, hr_plot) = plt.subplots(3)

    ne_model.plot(ts, a1, 'b', ts, a2, 'g')
    ne_model.legend(['$ A_1 $', '$ A_2 $'])
    ne_model.set(title='Norepinephrin', ylabel='Stoffkonzentration\n relativ zu $ A_0 $')
    ab_model.plot(ts, ab, 'r', linewidth=3)#, ts, b, 'y')
    # ab_model.legend(['$ AB $', '$ B $'])
    ab_model.set(title='Rezeptor (NE)', ylabel='Stoffkonzentration\n relativ zu $ B + AB $')
    ach_model.plot(ts, c2, 'b')
    ach_model.legend(['$ C_2 $'])
    ach_model.set(title='Acethylcholin', ylabel='Stoffkonzentration\n relativ zu $ C_1 $')
    vesicle_model.plot(ts, n, 'g')
    vesicle_model.legend(['$ N $'])
    vesicle_model.set(title='Vesikel', ylabel='Anzahl geladener Vesikel\n relativ zu $ N_{m} $')

    xs, hrs = ipfm(sol)
    print(hrs[3], min(hrs))
    rr_plot.plot(np.cumsum(xs)/s_to_ms, xs, linestyle='-', marker='o', markersize=2)
    rr_plot.set(title='RR-Intervalle', ylabel='Zeit in ms')
    # hr_plot.plot(ts, ab, 'r', linewidth=3)
    hr_plot.plot(np.cumsum(xs)/s_to_ms, hrs, linestyle='-', marker='o', markersize=2)
    # hr_plot.legend(['$ C_2 $', '$ HR $'])
    hr_plot.legend(['$ HR $'])
    hr_plot.set(title='Herzfrequenz (je Schlag)', ylabel='Herzfrequenz in Hz')

    fig.suptitle('Abszisse jeweils $ t $ in s', y='.07', size='10')
    plt.subplots_adjust(bottom=.15, wspace=.5, hspace=.4)
    plt.show()


solve_and_display(starting_val, t_end)
