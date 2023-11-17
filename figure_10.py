import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

r500_data = ascii.read('inner_cluster_data.csv')

quiescent_data = ascii.read('quiescent_field.csv')
sf_data = ascii.read('sf_field.csv')

r500_data = r500_data[np.where(r500_data['Outlier'] == 'n')]
r500_data = r500_data[np.where(r500_data['Ks/W1 log(M_stellar)'] != 'N/A')]
r500_data = r500_data[np.where(r500_data['Host/SN z'] > 0.01)]


quiescent_data = quiescent_data[np.where(quiescent_data['log(M_stellar)_err'] != 0)]
sf_data = sf_data[np.where(sf_data['log(M_stellar)_err'] != 0)]

quiescent_data = quiescent_data[np.where(quiescent_data['Outlier'] == 'n')]
sf_data = sf_data[np.where(sf_data['Outlier'] == 'n')]


r500_residuals = r500_data['HR_2']
r500_residual_errs = r500_data['mu_err_2']
r500_masses = r500_data['Ks/W1 log(M_stellar)']
r500_masserrs = r500_data['Ks/W1 log(M_stellar)_err']


quiescent_residuals = quiescent_data['HR']
quiescent_residual_errs = quiescent_data['mu_err']
sf_residuals = sf_data['HR']
sf_residual_errs = sf_data['mu_err']


quiescent_masses = quiescent_data['log(M_stellar)']
quiescent_mass_errs = quiescent_data['log(M_stellar)_err']
sf_masses = sf_data['log(M_stellar)']
sf_mass_errs = sf_data['log(M_stellar)_err']


masses = np.concatenate([quiescent_masses,sf_masses])
mass_errs = np.concatenate([quiescent_mass_errs,sf_mass_errs])
residuals = np.concatenate([quiescent_residuals,sf_residuals])
muerrs = np.concatenate([quiescent_residual_errs,sf_residual_errs])


avg_masses = 10
lower_bin = np.average(residuals[np.where(masses<avg_masses)],weights=1/muerrs[np.where(masses<avg_masses)]**2)
upper_bin = np.average(residuals[np.where(masses>avg_masses)],weights=1/muerrs[np.where(masses>avg_masses)]**2)

lower_bin_err = np.sqrt(1/np.sum(1/muerrs[np.where(masses<avg_masses)]**2))
upper_bin_err = np.sqrt(1/np.sum(1/muerrs[np.where(masses>avg_masses)]**2))


lower_avg = np.average(masses[np.where(masses<avg_masses)],weights=1/mass_errs[np.where(masses<avg_masses)]**2)
upper_avg = np.average(masses[np.where(masses>avg_masses)],weights=1/mass_errs[np.where(masses>avg_masses)]**2)

lower_std = np.sqrt(1/np.sum(1/mass_errs[np.where(masses<avg_masses)]**2))
upper_std = np.sqrt(1/np.sum(1/mass_errs[np.where(masses>avg_masses)]**2))

r500_masses = r500_masses.astype('float64')
r500_masserrs = r500_masserrs.astype('float64')


r500_lower_bin = np.average(r500_residuals[np.where(r500_masses<avg_masses)],weights=1/r500_residual_errs[np.where(r500_masses<avg_masses)]**2)
r500_upper_bin = np.average(r500_residuals[np.where(r500_masses>avg_masses)],weights=1/r500_residual_errs[np.where(r500_masses>avg_masses)]**2)

r500_lower_bin_err = np.sqrt(1/np.sum(1/r500_residual_errs[np.where(r500_masses<avg_masses)[0]]**2))
r500_upper_bin_err = np.sqrt(1/np.sum(1/r500_residual_errs[np.where(r500_masses>avg_masses)[0]]**2))


r500_lower_med = np.average(r500_masses[np.where(r500_masses<avg_masses)],weights=1/r500_masserrs[np.where(r500_masses<avg_masses)]**2)
r500_upper_med = np.average(r500_masses[np.where(r500_masses>avg_masses)],weights=1/r500_masserrs[np.where(r500_masses>avg_masses)]**2)

r500_lower_std = np.sqrt(1/np.sum(1/r500_masserrs[np.where(r500_masses<avg_masses)]**2))
r500_upper_std = np.sqrt(1/np.sum(1/r500_masserrs[np.where(r500_masses>avg_masses)]**2))


q_cluster_color = (0.6, 0.0, 0.0, 0.9)
sf_cluster_color = (0.0, 0.0, 0.6, 0.9)

plt.errorbar([lower_avg,upper_avg],[lower_bin,upper_bin],yerr=[lower_bin_err,upper_bin_err],xerr=[lower_std,upper_std],color='k',fmt='o',zorder=3,markersize=10)
plt.errorbar(quiescent_masses,quiescent_residuals,yerr=quiescent_residual_errs,xerr=quiescent_mass_errs,elinewidth=0.25,label='field quiescent',color='orangered',fmt='o',alpha=0.25,zorder=1)
plt.errorbar(sf_masses,sf_residuals,yerr=sf_residual_errs,xerr=sf_mass_errs,elinewidth=0.25,label='field star-forming',color='slateblue',fmt='d',alpha=0.25,zorder=1)
plt.errorbar([r500_lower_med,r500_upper_med],[r500_lower_bin,r500_upper_bin],yerr=[r500_lower_bin_err,r500_upper_bin_err],xerr=[r500_lower_std,r500_upper_std],color='gold',ecolor='k',mec='k',fmt='*',zorder=4,markersize=12,capsize=5)
plt.errorbar(r500_masses[(r500_data['Host SFR'] == 'Q')], r500_residuals[(r500_data['Host SFR'] == 'Q')], yerr = r500_residual_errs[(r500_data['Host SFR'] == 'Q')],xerr=r500_masserrs[(r500_data['Host SFR'] == 'Q')], color=q_cluster_color,fmt='*',alpha=0.6,markersize=8,zorder=4,label = 'inner cluster quiescent')
plt.errorbar(r500_masses[(r500_data['Host SFR'] == 'SF')], r500_residuals[(r500_data['Host SFR'] == 'SF')], yerr = r500_residual_errs[(r500_data['Host SFR'] == 'SF')],xerr=r500_masserrs[(r500_data['Host SFR'] == 'SF')],color=sf_cluster_color,fmt='*',alpha=0.6,markersize=8,zorder=4,label= 'inner cluster star-forming')


plt.axvline(avg_masses,linewidth='2',color='seagreen',zorder=5)

x = np.linspace(8,avg_masses,100)
nx = np.linspace(12,avg_masses,100)

plt.plot(x,[lower_bin]*len(x),color='k',linestyle='dashdot',zorder=3,linewidth=2)
plt.plot(nx,[upper_bin]*len(x),color='k',linestyle='dashdot',zorder=3,linewidth=2)
    
plt.ylabel(r'Hubble residual (mag)',fontsize=12)
plt.xlabel(r'log($M_{\ast}/M_{\odot}$)',fontsize=12)
plt.xlim(8,12)
plt.ylim(-0.8,0.6)
plt.legend(loc=3,)
plt.show()
plt.close()