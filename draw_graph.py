import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager

plt.rcParams.update({
    "text.usetex": True,
    "mathtext.fontset": 'stixsans',
    'text.latex.preamble': [r'\usepackage{cmbright}', r'\usepackage{amsmath}']
    })

cold_energies = np.load('cold_energies.npy')
hot_energies = np.load('hot_energies.npy')
plt.figure(figsize = (5, 3))
xs = range(0, 12 * 64 * 64, 64 * 64)
plt.plot(xs, cold_energies, marker = 'o', label = "cold start")
plt.plot(xs, hot_energies, marker = 'o', label = "hot start")
plt.xlabel("Iterations")
plt.ylabel("$\hat{H} / N_{sites}$");
plt.legend()

plt.tight_layout()
plt.savefig('hot_cold.pdf')

plt.clf()
with np.load('correlation_b03.npy') as data:
    plt.figure(figsize = (5, 3))
    correlations = data['correlations'].real[0:10]
    error = abs(data['error'])[0:10]
    print(error)
    plt.errorbar(range(len(correlations)), np.log(correlations), yerr = np.log((correlations + error) / correlations), fmt = '--.', capsize = 3)
    plt.xlabel("$t$")
    plt.ylabel("log $C_{0}(t)$")
    plt.tight_layout()
    plt.savefig('correlation_graph.pdf')
