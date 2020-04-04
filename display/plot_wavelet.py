# script for displaying the third derivative of the Gaussian function
#
# author: lucaswzyk

import math
import matplotlib.pyplot as plt

xs = [x / 10 for x in range(-60, 60)]
ys = [- math.exp(-math.pow(x, 2) / 2) * x * (math.pow(x, 2) - 3) / math.sqrt(2*math.pi) for x in xs]
plt.plot(xs, ys)
plt.show()
