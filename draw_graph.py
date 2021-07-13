import numpy as np
import matplotlib.pyplot as plt

cold_energies = np.load('cold_energies.npy')
hot_energies = np.load('hot_energies.npy')
plt.plot(cold_energies)
plt.plot(hot_energies)
plt.show()

