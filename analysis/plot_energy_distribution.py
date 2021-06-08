import numpy as np
import matplotlib.pyplot as plt

energy_list = []
i_energy = 3
n_bins = 100
E0 = 120

with open("../../JETSCAPE-output/my_final_state_partons.txt", "r") as data_file: 
    lines = data_file.readlines()
    for line in lines: 
        if '#' in line:
            continue
        arr = line.split()
        energy_list.append(float(arr[i_energy]))

plt.figure()
plt.hist(energy_list, n_bins, density=True, histtype='step', stacked=True, fill=False)
plt.title("Energy distribution of final partons, $E_0 = 120$ GeV")
plt.xlabel('E (GeV)')
plt.ylabel('dN/dE')
plt.yscale('log')
plt.xlim(1, 120)
plt.savefig("../../JETSCAPE-output/parton_energy_distribution.pdf")

