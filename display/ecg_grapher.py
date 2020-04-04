# simple example of BioSPPy function
# plots width_t seconds of the given ecg file (.txt) starting at offset_t with
# the given resolution
#
# author: lucaswzyk

import matplotlib.pyplot as plt
import numpy as np
from biosppy import ecg, storage

file = 'ecg.txt'
width_t = 2
resolution = 1000
offset_t = 0

ts = [(i/250) for i in range(0+offset_t, width_t+offset_t)]
signal, mdata = storage.load_txt(file)
# ecg = [(np.sin(ts[i]*2*np.pi*39)*10**19*.25-ecg[i]) for i in range(0, width_t)]

out = ecg.ecg(signal[0:width_t*resolution], resolution, show=True)
# print(out['filtered'])
# plt.plot(range(0, len(out['filtered'])), out['filtered'])
# plt.show()

# plt.plot(ts, ecg)
# plt.show()
