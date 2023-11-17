import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import astropy.units as u
import matplotlib.colors as c
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0 = 70 * u.km/u.s/u.Mpc,Tcmb0 = 2.725*u.K,Om0=0.3)

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

r500_x1s = r500_data['x1']
r500_cs = r500_data['c']
r500_mBs = r500_data['mB']
r500_hd_redshifts = r500_data['Cluster z (HD)']

r500_x1_errs = r500_data['x1_err']
r500_c_errs = r500_data['c_err']
r500_mB_errs = r500_data['mB_err']

quiescent_x1s = quiescent_data['x1']
quiescent_cs = quiescent_data['c']
quiescent_mBs = quiescent_data['mB']
quiescent_hd_redshifts = quiescent_data['Host z (HD)']

quiescent_x1_errs = quiescent_data['x1_err']
quiescent_c_errs = quiescent_data['c_err']
quiescent_mB_errs = quiescent_data['mB_err']

sf_x1s = sf_data['x1']
sf_cs = sf_data['c']
sf_mBs = sf_data['mB']
sf_hd_redshifts = sf_data['Host z (HD)']

sf_x1_errs = sf_data['x1_err']
sf_c_errs = sf_data['c_err']
sf_mB_errs = sf_data['mB_err']

alpha = 0.127
beta = 2.395


q_cluster_color = (0.6, 0.0, 0.0, 0.9)
sf_cluster_color = (0.0, 0.0, 0.6, 0.9)

r500_mus = 5*np.log10(cosmo.luminosity_distance(r500_hd_redshifts)/ u.Mpc) + 25
quiescent_mus = 5*np.log10(cosmo.luminosity_distance(quiescent_hd_redshifts)/ u.Mpc) + 25
sf_mus = 5*np.log10(cosmo.luminosity_distance(sf_hd_redshifts)/ u.Mpc) + 25


quiescent_MBs = -(quiescent_mus - quiescent_mBs + beta*quiescent_cs)
quiescent_MB_errs = np.sqrt(quiescent_mB_errs**2 + beta**2*quiescent_c_errs**2 + alpha**2 * quiescent_x1_errs**2)

sf_MBs = -(sf_mus - sf_mBs + beta*sf_cs)
sf_MB_errs = np.sqrt(sf_mB_errs**2 + beta**2*sf_c_errs**2 + alpha**2 * sf_x1_errs**2)

r500_MBs = -(r500_mus - r500_mBs + beta*r500_cs)
r500_MB_errs = np.sqrt(r500_mB_errs**2 + beta**2*r500_c_errs**2 + alpha**2 * r500_x1_errs**2)

alpha_lower,beta_lower = -0.20251304994516123,-19.38988959131312
alpha_lower_err,beta_lower_err = 0.019792953791067422,0.03731504400093384

alpha_upper,beta_upper = -0.08737900127943338,-19.284595641164078

lower_x_dummy = np.arange(-3,-0.7,0.001)
lower_MB_dummy = alpha_lower*lower_x_dummy + beta_lower

upper_x_dummy = np.arange(-1.3,3,0.001)
upper_MB_dummy = alpha_upper*upper_x_dummy + beta_upper

plt.errorbar(quiescent_x1s, quiescent_MBs, xerr = quiescent_x1_errs,yerr=np.sqrt(0.14**2 + quiescent_MB_errs**2),color='orangered',label = 'field quiescent',alpha=0.25,elinewidth=0.25,fmt = 'o',zorder=1)
plt.errorbar(sf_x1s,sf_MBs,xerr = sf_x1_errs, yerr=np.sqrt(0.14**2 + sf_MB_errs**2),color='slateblue',label='field star-forming',alpha=0.25,elinewidth=0.25,fmt='d',zorder=1)
plt.errorbar(r500_x1s[(r500_data['Host SFR'] == 'Q')],r500_MBs[(r500_data['Host SFR'] == 'Q')],xerr = r500_x1_errs[(r500_data['Host SFR']== 'Q')],yerr=np.sqrt(0.14**2 + r500_MB_errs[(r500_data['Host SFR'] == 'Q')]**2),color=q_cluster_color,label='inner cluster quiescent',fmt='*',alpha=0.6,zorder=4,markersize=8)
plt.errorbar(r500_x1s[(r500_data['Host SFR'] == 'SF')],r500_MBs[(r500_data['Host SFR'] == 'SF')],xerr = r500_x1_errs[(r500_data['Host SFR'] == 'SF')],yerr=np.sqrt(0.14**2 + r500_MB_errs[(r500_data['Host SFR'] == 'SF')]**2),color=sf_cluster_color,label='inner cluster star-forming',fmt='*',alpha=0.6,zorder=4,markersize=8)


plt.plot(lower_x_dummy, lower_MB_dummy,color='k',zorder=3)
plt.plot(upper_x_dummy,upper_MB_dummy,color='k',linestyle='dashed',zorder=3)
plt.xlim(-3.7,3.7)
plt.ylim(-18.5,-20)
plt.gca().invert_xaxis()
plt.ylabel(r'$M_{corr} = m_B - \beta\,c - \mu_{\rm cosmo}$ (mag)')
plt.xlabel(r'$x_1$')
plt.legend(loc = 1)
plt.show()
plt.close()