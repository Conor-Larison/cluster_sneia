from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),
     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]


r500_data = ascii.read('inner_cluster_data.csv')
twor500_data = ascii.read('outer_cluster_data.csv')

quiescent_data = ascii.read('quiescent_field.csv')
sf_data = ascii.read('sf_field.csv')

r500_cluster_redshift = r500_data['Cluster z']
twor500_cluster_redshift = twor500_data['Cluster z']

quiescent_redshifts = quiescent_data['Host z']
sf_redshifts = sf_data['Host z']

fig, axs = plt.subplots(1,1,figsize=(6,4))

w = 0.01
axs.hist(r500_cluster_redshift,color='k',zorder=5,label=r'inner cluster',histtype='step',density=False,bins=np.arange(0,0.1, w),linewidth=3)
axs.hist(twor500_cluster_redshift,color='limegreen',zorder=1,label=r'outer cluster',histtype='step',linestyle='dashdot',density=False,bins=np.arange(0,0.1, w),linewidth=3)
axs.set_xlabel('redshift',fontsize=12)
axs.set_ylabel('number',fontsize=12)

(quiescent_counts, quiescent_bins) = np.histogram(quiescent_redshifts, bins=np.arange(0,0.1, w))
(sf_counts, sf_bins) = np.histogram(sf_redshifts, bins=np.arange(0,0.1, w))

factor = 5
axs.hist(quiescent_bins[:-1], quiescent_bins, weights=quiescent_counts/factor,color='orangered',zorder=1,label='field quiescent / 5',histtype='step',linestyle='dashed',density=False,linewidth=3)
axs.hist(sf_bins[:-1], sf_bins, weights=sf_counts/factor,color='slateblue',zorder=1,label='field star-forming / 5',histtype='step',linestyle=linestyle_tuple[3][-1],density=False,linewidth=3)
axs.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

plt.tight_layout()
plt.legend(prop={'size':11},loc=2)
plt.ylim(0,20)
plt.show()
plt.close()