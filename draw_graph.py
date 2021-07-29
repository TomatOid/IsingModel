import numpy as np
import matplotlib.pyplot as plt

cold_energies = np.load('cold_energies.npy')
hot_energies = np.load('hot_energies.npy')
xs = range(0, 12 * 64 * 64, 64 * 64)
plt.plot(xs, cold_energies, marker = 'o', label = "cold start")
plt.plot(xs, hot_energies, marker = 'o', label = "hot start")
plt.xlabel("Iterations")
plt.ylabel("H / N");
plt.legend()

plt.savefig('hot_cold.pdf')

plt.clf()
with np.load('correlation_b03.npy') as data:
    correlations = data['correlations'].real
    error = abs(data['error'])
    print(error)
    plt.errorbar(range(len(correlations)), correlations, yerr = error, fmt = 'o')
    plt.show();
