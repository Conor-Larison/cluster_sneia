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
r500_masses = r500_data['Ks/W1 log(M_stellar)'].astype('float64')
r500_sfr_cats = r500_data['Host SFR']

r500_x1_errs = r500_data['x1_err']
r500_c_errs = r500_data['c_err']
r500_mB_errs = r500_data['mB_err']

quiescent_x1s = quiescent_data['x1']
quiescent_cs = quiescent_data['c']
quiescent_mBs = quiescent_data['mB']
quiescent_hd_redshifts = quiescent_data['Host z (HD)']
quiescent_masses = quiescent_data['log(M_stellar)']

quiescent_x1_errs = quiescent_data['x1_err']
quiescent_c_errs = quiescent_data['c_err']
quiescent_mB_errs = quiescent_data['mB_err']

sf_x1s = sf_data['x1']
sf_cs = sf_data['c']
sf_mBs = sf_data['mB']
sf_hd_redshifts = sf_data['Host z (HD)']
sf_masses = sf_data['log(M_stellar)']

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


fig,axs = plt.subplots(3,1,figsize=(5,11),sharex=True,sharey=True)
fig.subplots_adjust(hspace=0,wspace=0)
plt.tight_layout

lower_alpha_upper,lower_intercept_upper = -0.21014713990253098,-19.383121692434912

lower_x_upper = np.arange(-3,-0.7,0.1)
lower_y_upper = lower_x_upper * lower_alpha_upper + lower_intercept_upper

lower_alpha_lower,lower_intercept_lower = -0.09060476120815121,-19.2618179676212

lower_x_lower = np.arange(3,-1.3,-0.1)
lower_y_lower = lower_x_lower * lower_alpha_lower + lower_intercept_lower


axs[0].errorbar(sf_x1s[(sf_masses < 10.5)], sf_MBs[(sf_masses < 10.5)], yerr = np.sqrt(0.14**2 + sf_MB_errs[(sf_masses < 10.5)]**2),xerr=sf_x1_errs[(sf_masses < 10.5)],color='slateblue',fmt='d',alpha=0.25,elinewidth=0.25,zorder = 1)
axs[0].errorbar(quiescent_x1s[(quiescent_masses < 10.5)], quiescent_MBs[(quiescent_masses < 10.5)], yerr = np.sqrt(0.14**2 + quiescent_MB_errs[(quiescent_masses < 10.5)]**2),xerr=quiescent_x1_errs[(quiescent_masses < 10.5)],color='orangered',fmt='o',alpha=0.25,elinewidth=0.25,zorder = 1)
axs[0].errorbar(r500_x1s[(r500_masses < 10.5) & (r500_sfr_cats == 'SF')], r500_MBs[(r500_masses < 10.5) & (r500_sfr_cats == 'SF')], yerr = np.sqrt(0.14**2 + r500_MB_errs[(r500_masses < 10.5) & (r500_sfr_cats == 'SF')]**2),xerr=r500_x1_errs[(r500_masses < 10.5) & (r500_sfr_cats == 'SF')],color=sf_cluster_color,fmt='*',alpha=0.6,markersize=8, zorder=4)
axs[0].errorbar(r500_x1s[(r500_masses < 10.5) & (r500_sfr_cats == 'Q')], r500_MBs[(r500_masses < 10.5) & (r500_sfr_cats == 'Q')], yerr = np.sqrt(0.14**2 + r500_MB_errs[(r500_masses < 10.5) & (r500_sfr_cats == 'Q')]**2),xerr=r500_x1_errs[(r500_masses < 10.5) & (r500_sfr_cats == 'Q')],color=q_cluster_color,fmt='*',alpha=0.6,markersize=8, zorder = 4)
axs[0].plot(lower_x_upper, lower_y_upper,color='k',zorder=3)
axs[0].plot(lower_x_lower, lower_y_lower,color='k',linestyle='dashed',zorder=3)
axs[0].set_xlim(3.7,-3.7)
axs[0].set_ylim(-18.5,-19.9)
axs[0].set_title(r'log($M_{\ast}/M_{\odot}$)$ < 10.5$',fontsize=14,y=0)
axs[0].tick_params(axis='x', labelsize=14)
axs[0].tick_params(axis='y', labelsize=14)


mid_alpha_upper,mid_intercept_upper = -0.1888605072946966,-19.37625873765778


mid_x_upper = np.arange(-3,-0.7,0.1)
mid_y_upper = mid_x_upper * mid_alpha_upper + mid_intercept_upper


mid_alpha_lower,mid_intercept_lower = -0.12208388947703709,-19.31013458076518


mid_x_lower = np.arange(3,-1.3,-0.1)
mid_y_lower = mid_x_lower * mid_alpha_lower + mid_intercept_lower

axs[1].errorbar(sf_x1s[(sf_masses < 11) & (sf_masses >= 10.5)], sf_MBs[(sf_masses < 11) & (sf_masses >= 10.5)], yerr = np.sqrt(0.14**2 + sf_MB_errs[(sf_masses < 11) & (sf_masses >= 10.5)]**2),xerr=sf_x1_errs[(sf_masses < 11) & (sf_masses >= 10.5)],color='slateblue',fmt='d',alpha=0.25,elinewidth=0.25, zorder = 1)
axs[1].errorbar(quiescent_x1s[(quiescent_masses < 11) & (quiescent_masses >= 10.5)], quiescent_MBs[(quiescent_masses < 11) & (quiescent_masses >= 10.5)], yerr = np.sqrt(0.14**2 + quiescent_MB_errs[(quiescent_masses < 11) & (quiescent_masses >= 10.5)]**2),xerr=quiescent_x1_errs[(quiescent_masses < 11) & (quiescent_masses >= 10.5)],color='orangered',fmt='o',alpha=0.25,elinewidth=0.25, zorder = 1)
axs[1].errorbar(r500_x1s[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'SF')], r500_MBs[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'SF')], yerr = np.sqrt(0.14**2 + r500_MB_errs[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'SF')]**2),xerr=r500_x1_errs[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'SF')],color=sf_cluster_color,fmt='*',alpha=0.6,markersize=8,zorder = 4)
axs[1].errorbar(r500_x1s[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'Q')], r500_MBs[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'Q')], yerr = np.sqrt(0.14**2 + r500_MB_errs[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'Q')]**2),xerr=r500_x1_errs[(r500_masses < 11) & (r500_masses >= 10.5) & (r500_sfr_cats == 'Q')],color=q_cluster_color,fmt='*',alpha=0.6,markersize=8, zorder=4)
axs[1].plot(mid_x_upper, mid_y_upper,color='k')
axs[1].plot(mid_x_lower, mid_y_lower,color='k',linestyle='dashed')
axs[1].set_xlim(3.7,-3.7)
axs[1].set_ylim(-18.5,-19.9)
axs[1].set_title(r'$10.5 < $log($M_{\ast}/M_{\odot}$)$ < 11.0$',fontsize=14,y=0)
axs[1].tick_params(axis='x', labelsize=14)
axs[1].tick_params(axis='y', labelsize=14)

upper_alpha_upper,upper_intercept_upper = -0.21925764948029253,-19.436024115953945

upper_x_upper = np.arange(-3,-0.7,0.1)
upper_y_upper = upper_x_upper * upper_alpha_upper + upper_intercept_upper

upper_alpha_lower,upper_intercept_lower = -0.06971999005850965,-19.357773226266417

upper_x_lower = np.arange(3,-1.3,-0.1)
upper_y_lower = upper_x_lower * upper_alpha_lower + upper_intercept_lower

axs[2].errorbar(sf_x1s[(sf_masses >= 11)], sf_MBs[(sf_masses >= 11)], yerr = np.sqrt(0.14**2 + sf_MB_errs[(sf_masses >= 11)]**2),xerr=sf_x1_errs[(sf_masses >= 11)],color='slateblue',fmt='d',alpha=0.25,elinewidth=0.25,zorder=1)
axs[2].errorbar(quiescent_x1s[(quiescent_masses >= 11)], quiescent_MBs[(quiescent_masses >= 11)], yerr = np.sqrt(0.14**2 + quiescent_MB_errs[(quiescent_masses >= 11)]**2),xerr=quiescent_x1_errs[(quiescent_masses >= 11)],color='orangered',fmt='o',alpha=0.25,elinewidth=0.25,zorder=1)
axs[2].errorbar(r500_x1s[(r500_masses >= 11) & (r500_sfr_cats == 'SF')], r500_MBs[(r500_masses >= 11) & (r500_sfr_cats == 'SF')], yerr = np.sqrt(0.14**2 + r500_MB_errs[(r500_masses >= 11) & (r500_sfr_cats == 'SF')]**2),xerr=r500_x1_errs[(r500_masses >= 11) & (r500_sfr_cats == 'SF')],color=sf_cluster_color,fmt='*',alpha=0.6,markersize=8,zorder=4)
axs[2].errorbar(r500_x1s[(r500_masses >= 11) & (r500_sfr_cats == 'Q')], r500_MBs[(r500_masses >= 11) & (r500_sfr_cats == 'Q')], yerr = np.sqrt(0.14**2 + r500_MB_errs[(r500_masses >= 11) & (r500_sfr_cats == 'Q')]**2),xerr=r500_x1_errs[(r500_masses >= 11) & (r500_sfr_cats == 'Q')],color=q_cluster_color,fmt='*',alpha=0.6,markersize=8,zorder=4)
axs[2].plot(upper_x_upper, upper_y_upper,color='k')
axs[2].plot(upper_x_lower, upper_y_lower,color='k',linestyle='dashed')
axs[2].set_xlim(3.7,-3.7)
axs[2].set_ylim(-18.5,-19.9)
axs[2].set_title(r'log$(M_{\ast}/M_{\odot}$)$ > 11.0$',fontsize=14,y=0.0)
axs[2].tick_params(axis='x', labelsize=14)
axs[2].tick_params(axis='y', labelsize=14)
axs[2].set_xticks(np.arange(-3,4))
axs[2].set_xlabel(r'$x_1$',fontsize=14)


plt.show()
plt.close()