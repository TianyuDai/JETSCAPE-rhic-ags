import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d

fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(4, 1, hspace=0, wspace=0)
data_th = np.loadtxt('../../JETSCAPE-output/pp2760/pp2760_chargedHadron.txt')
data_ex = np.loadtxt('pp2760_chargedHadron_data.txt')

th_x = data_th.T[0]
th_y = data_th.T[1]
th_err = data_th.T[2]
th_interp = interp1d(th_x, th_y)

ex_x = (data_ex.T[0] + data_ex.T[1]) / 2
ex_y = data_ex.T[2]
ex_xerr = (data_ex.T[1]-data_ex.T[0])/2
ex_yerr = np.sqrt(data_ex.T[3]**2+data_ex.T[4]**2)
th_interp_y = [th_interp(pt) for pt in ex_x]

for i in range(2): 
    if (i == 0): 
        ax = fig.add_subplot(gs[:-1])
        ax.fill_between(th_x, th_y-th_err, th_y+th_err, alpha=0.2, edgecolor='none', label='JETSCAPE')
        ax.plot(th_x, th_y, linewidth=0.5)
        ax.errorbar(ex_x, ex_y, xerr=ex_xerr, yerr=ex_yerr, fmt='.', markersize=5, elinewidth=0.1, label='CMS')
        ax.set_yscale('log')
        ax.set_ylim(1e-13, 1e-5)
        ax.set_ylabel('$\\frac{d^3N}{d\eta d^2p_T} \left(\\frac{1}{GeV^2\cdot c^2}\\right)$')
        ax.tick_params(axis='x', which='both', bottom=False)
        ax.legend()
    else: 
        ax = fig.add_subplot(gs[-1])
        ax.errorbar(ex_x, ex_y/th_interp_y, xerr=ex_xerr, yerr=np.sqrt(ex_yerr**2/th_y[:-1]**2+th_err[:-1]**2/ex_y**2), fmt='.', label='JETSCAPE ratio', markersize=5, elinewidth=0.1)
        ratio_baseline = [1 for x in ex_x]
        ax.plot(ex_x, ratio_baseline, '--', color='black', markersize=1)
        ax.set_yscale('linear')
        ax.set_ylim(0., 2.)
        ax.set_ylabel('data/theory')
        ax.legend(loc='lower right')
    ax.set_xlim(10, 100)
    ax.tick_params(direction="in", which='both')

fig.suptitle('pp, 2760GeV, charged hadron, $|\eta|<1$, JETSCAPE vs. CMS12')
plt.savefig('../../JETSCAPE-output/pp2760/pp2760_chargedHadron_data_compare.pdf')
