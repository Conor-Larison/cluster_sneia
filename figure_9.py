import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0 = 70 * u.km/u.s/u.Mpc,Tcmb0 = 2.725*u.K,Om0=0.3)


r500_data = ascii.read('inner_cluster_data.csv')
twor500_data = ascii.read('outer_cluster_data.csv')

r500_data = r500_data[np.where(r500_data['Cluster z'] > 0.01)]
r500_data = r500_data[np.where(r500_data['Outlier'] == 'n')]

twor500_data = twor500_data[np.where(twor500_data['Cluster z'] > 0.01)]
twor500_data = twor500_data[np.where(twor500_data['Outlier'] == 'n')]

r500_moduli = r500_data['mu_1']
twor500_moduli = twor500_data['mu_1']

r500_muerrs = r500_data['mu_err_1']
twor500_muerrs = twor500_data['mu_err_1']

r500_residuals = r500_data['HR_1']
twor500_residuals = twor500_data['HR_1']

r500_hd_redshifts = r500_data['Cluster z (HD)']
twor500_hd_redshifts = twor500_data['Cluster z (HD)']

z = np.arange(0.001,0.5,0.00001)
lumdistance = [cosmo.luminosity_distance(i) for i in z]
mu_pred = [5*np.log10(i / u.Mpc) + 25 for i in lumdistance]


fig,(ax1,ax2) = plt.subplots(2,1,figsize=(7.5,6),sharex=True,gridspec_kw={'height_ratios': [2.5, 1]})
fig.subplots_adjust(hspace=0)
ax1.errorbar(r500_hd_redshifts,r500_moduli,yerr=r500_muerrs,fmt='o',color='k',label=r'inner cluster $\mu_{obs}$')
ax1.errorbar(twor500_hd_redshifts,twor500_moduli,yerr=twor500_muerrs,fmt='*',color='limegreen',label=r'outer cluster $\mu_{obs}$')
ax1.plot(z,mu_pred,label=r'$\mu_{cosmo}$')
ax1.legend(prop={'size':12})
ax1.set_ylabel(r'$\mu$ (mag)',fontsize=14)
ax1.set_ylim(32.5,39)
ax1.set_xlim(0.008,0.1)
ax1.set_yticks(np.arange(33,40,1))
ax1.tick_params(axis='y', labelsize=12)
ax1.text(.70, .16, 'inner cluster RMS residual: 0.163 mag', ha='center', va='top', transform=ax1.transAxes,color = 'k',fontsize=12)
ax1.text(.70, .10, 'outer cluster RMS residual: 0.163 mag', ha='center', va='top', transform=ax1.transAxes,color = 'k',fontsize=12)


ax2.errorbar(r500_hd_redshifts,r500_residuals,yerr=r500_muerrs,fmt='o',color='k')
ax2.errorbar(twor500_hd_redshifts,twor500_residuals,yerr=twor500_muerrs,fmt='*',color='limegreen')
ax2.plot(z,[0]*len(z))
ax2.set_ylabel(r'HR (mag)',fontsize=14)
ax2.set_xlabel('redshift',fontsize=14)
ax2.set_xticks(np.arange(0.01,0.11,0.01))
ax2.set_yticks(np.arange(-0.4,0.5,0.4))
ax2.tick_params(axis='x', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)


plt.show()
plt.close()