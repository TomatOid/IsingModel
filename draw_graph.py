import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
import scipy.optimize as opt

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
plt.ylabel("$\hat{H} / N_{sites}$")
plt.legend()

plt.tight_layout()
plt.savefig('hot_cold.pdf')

plt.clf()
with np.load('corr_all_0.6.npy') as data:
    #plt.plot(data['correlations'].real)
    #plt.show()
    plt.figure(figsize = (5, 3))
    corr_raw = data['correlations'].real[1:6, 0:12]
    err_raw  = data['error'][1:6, 0:12]
    for correlations, error, n in zip(corr_raw, err_raw, range(len(corr_raw))):
        print(correlations)
        print(error)
        x_range  = np.arange(len(correlations))
        log_corr = np.log(correlations)
        errorbar = np.log((correlations + error) / correlations)
        nan_filt = np.logical_not(np.logical_or(np.isnan(log_corr), np.isnan(errorbar)))
        bars = plt.errorbar(x_range[nan_filt], log_corr[nan_filt], yerr = errorbar[nan_filt], fmt = '.', capsize = 3)
        #fit, fit_cov = opt.curve_fit(lambda x, m: m * x, x_range[nan_filt][1:], log_corr[nan_filt][1:], sigma = errorbar[nan_filt][1:])
        fit, fit_cov = opt.curve_fit(lambda x, E: np.exp(E * x), x_range[1:], correlations[1:], sigma = error[1:])
        plt.plot(range(len(correlations)), np.linspace(0, fit[0] * (len(correlations) - 1), len(correlations)),
                 color = bars.lines[0].get_color(), label = f'$E_{{{n}}} = {-fit[0]:.4f}$')
        print(fit)
    plt.legend()
    plt.xlabel("$t$")
    plt.ylabel("$\\log\\frac{C_{n}(t)}{\\mathcal{A}_n}$")
    plt.tight_layout()
    plt.savefig('correlation_graph_all.pdf')
