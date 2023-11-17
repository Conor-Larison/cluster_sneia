from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

r500_data = ascii.read('inner_cluster_data.csv')
twor500_data = ascii.read('outer_cluster_data.csv')

r500_data = r500_data[np.where(r500_data['Host SFR'] == 'Q')]
twor500_data = twor500_data[np.where(twor500_data['Host SFR'] == 'Q')]

quiescent_data = ascii.read('quiescent_field.csv')

r500_x1s = r500_data['x1']
twor500_x1s = twor500_data['x1']

quiescent_x1s = quiescent_data['x1']

fig, axs = plt.subplots(1,1,figsize=(6,4))

w = 0.5
axs.hist(r500_x1s,color='k',zorder=5,label=r'quiescent inner cluster',histtype='step',density=False,bins=np.arange(-3, 3, w),linewidth=3)
axs.hist(twor500_x1s,color='limegreen',zorder=1,label=r'quiescent outer cluster',linestyle='dashdot',histtype='step',density=False,bins=np.arange(-3,3, w),linewidth=3)
axs.set_xlabel(r'$x_1$',fontsize=16)
axs.set_ylabel('number',fontsize=14)
axs.set_xlim(-3,3)

factor = 5
(quiescent_counts, quiescent_bins) = np.histogram(quiescent_x1s, bins=np.arange(-3,3, w))

axs.hist(quiescent_bins[:-1], quiescent_bins, weights=quiescent_counts/factor,color='orangered',zorder=1,label='field quiescent / 5',histtype='step',linestyle='dashed',density=False,linewidth=3)

axs.legend(loc=1,prop={'size': 12})
axs.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

plt.tight_layout()
plt.show()
plt.close()