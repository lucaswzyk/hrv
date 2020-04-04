# writes values of sin wave similar to RR values into a file
#
# author: lucaswzyk

import math
import matplotlib.pyplot as plt

f = open('../material/sin.txt', 'w')
data = [600 + 100 * math.sin(i * math.pi / 1200) for i in range(25000)]

for i in range(25000):
    f.write(str(600 + 100 * math.sin(i * math.pi / 1200)) + "\n")

plt.plot(list(range(25000)), data)
plt.show()
