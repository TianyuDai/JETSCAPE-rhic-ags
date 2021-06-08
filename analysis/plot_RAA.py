import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d

data_pp = np.loadtxt("../../JETSCAPE-output/AuAu200/pp200_chargedHadron.txt")
data_AA = np.loadtxt("../../JETSCAPE-output/AuAu200/AuAu200_chargedHadron.txt")

xpp = np.array(data_pp.T[0])
ypp = np.array(data_pp.T[1])
zpp = np.array(data_pp.T[2])

xAA = np.array(data_AA.T[0])
yAA = np.array(data_AA.T[1])
zAA = np.array(data_AA.T[2])

bins = np.linspace(15, 55, 5)
result = yAA/ypp

error = np.sqrt(result**2*((zAA/yAA)**2+(zpp/ypp)**2))

plt.figure(figsize=(10, 8))
plt.fill_between(bins, result-error, result+error, facecolor='cornflowerblue', alpha=0.2, edgecolor="none")
plt.plot(bins, result, color='cornflowerblue')
plt.text(16.5, 1.08, "$\sqrt{s_{NN}}=200GeV$, charged hadron, \n$10-20 \%$, $|y| \leq 1$", size=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
plt.tick_params(axis="x", labelsize=20)
plt.tick_params(axis="y", labelsize=20)
plt.text(16.5, 0.6, "single hydrodynamic event\n$\\alpha_s = 0.3$", ha='left', size=12)
plt.xlabel("$p_T$(GeV)", fontsize=20)
plt.ylabel("$R_{AA}^{h^{\pm}}$", fontsize=20)
plt.xlim(15, 55)
plt.savefig('../../JETSCAPE-output/AuAu200/AuAu200_chargedHadron_RAA_oneBin.pdf')

