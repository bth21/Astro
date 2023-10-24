
import numpy as np
import matplotlib.pyplot as plt

brightness = np.loadtxt('brightness.csv', delimiter = ',', skiprows = 1)
eccent = np.loadtxt('eccent.csv', delimiter = ',', skiprows = 1)

plt.scatter(eccent, brightness)