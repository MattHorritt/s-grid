# Quick check of tide level calculations - pasted in from ltis_slr as want to do
# independently of grid

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

slr = 0.0
tidePeriod = 12.42  # Period between high tides
t1 = 2.91
t200 = 3.41
t_list = np.arange(0,72 * 3600, 30)
h_list = []

for t in t_list:
    theta = 2 * np.pi * t / (3600. * tidePeriod)

    if theta <= 2 * np.pi or theta > 3 * np.pi:
        wl = np.sin(theta) * t1
    else:
        wl = np.sin(theta) * t200

    wl += slr

    h_list.append(wl)

matplotlib.use("Qt5Agg")
plt.plot(t_list/3600, h_list)
plt.show()

