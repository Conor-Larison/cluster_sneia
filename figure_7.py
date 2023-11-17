import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from scipy.stats import binned_statistic

quiescent_data = ascii.read('quiescent_field.csv')
r500_data = ascii.read('inner_cluster_data.csv')

r500_data = r500_data[np.where(r500_data['Host SFR'] == 'Q')]

r500_x1s = r500_data['x1']
quiescent_x1s = quiescent_data['x1']

r500_x1errs = r500_data['x1_err']
quiescent_x1errs = quiescent_data['x1_err']

r500_redshifts = r500_data['Host/SN z']
quiescent_redshifts = quiescent_data['Host z']


bin_stats = binned_statistic(quiescent_redshifts,quiescent_x1s,statistic='median',bins=7)
bin_counts = binned_statistic(quiescent_redshifts,quiescent_x1s,statistic='count',bins=7)
bin_stds = binned_statistic(quiescent_redshifts,quiescent_x1s,statistic='std',bins=7)

bin_errs = bin_stds[0][0:7]/np.sqrt(bin_counts[0][0:7])

r500_bin_stats = binned_statistic(r500_redshifts,r500_x1s,statistic='median',bins=5)
r500_bin_counts = binned_statistic(r500_redshifts,r500_x1s,statistic='count',bins=5)
r500_bin_stds = binned_statistic(r500_redshifts,r500_x1s,statistic='std',bins=5)

r500_bin_errs = r500_bin_stds[0][0:5]/np.sqrt(r500_bin_counts[0][0:5])

bin_width = bin_stats[1][1] - bin_stats[1][0]
r500_bin_width = r500_bin_stats[1][1] - r500_bin_stats[1][0]

plt.errorbar(quiescent_redshifts,quiescent_x1s,yerr=quiescent_x1errs,alpha=0.25,fmt='o',color='orangered',label='field quiescent')
plt.errorbar(r500_redshifts,r500_x1s,yerr=r500_x1errs,color='k',label='inner cluster quiescent',fmt='d')
plt.errorbar(bin_stats[1][0:7] + bin_width/2,bin_stats[0][0:7],yerr=bin_errs,color='orangered',linestyle='--',label='field quiescent trend',linewidth=3)
plt.errorbar(r500_bin_stats[1][0:5] + r500_bin_width/2,r500_bin_stats[0][0:5],yerr=r500_bin_errs,color='k',label=r'inner cluster quiescent trend',linewidth=3)


plt.xlabel('redshift',fontsize=14)
plt.ylabel(r'$x_1$',fontsize=16)
plt.xlim(0,0.1)
plt.ylim(-3.5,3.5)
plt.legend(prop={'size':9})
plt.show()
plt.close()