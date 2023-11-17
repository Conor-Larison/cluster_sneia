import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

r500_data = ascii.read('inner_cluster_data.csv')
twor500_data = ascii.read('outer_cluster_data.csv')
quiescent_data = ascii.read('quiescent_field.csv')
sf_data = ascii.read('sf_field.csv')

r500_x1s = r500_data['x1']
twor500_x1s = twor500_data['x1']
quiescent_x1s = quiescent_data['x1']
sf_x1s = sf_data['x1']

combined_cluster_x1s = np.concatenate([r500_x1s,twor500_x1s])

mu1, sig1, mu2, sig2= [-2.013,0.428,0.286,0.590]

data = np.arange(-3.2,3.2,0.1)
w = 0.5

r500_frac1 = 0.76

r500_model1 = (len(r500_x1s)*w)*r500_frac1/np.sqrt(2.0*np.pi*(sig1**2)) \
            * np.exp(-(data-mu1)**2/2.0/(sig1**2))

r500_model2 = (len(r500_x1s)*w)*(1.0 - r500_frac1)/np.sqrt(2.0*np.pi*(sig2**2)) \
            * np.exp(-(data-mu2)**2/2.0/(sig2**2))

r500_model = (r500_model1 + r500_model2)

cluster_frac1 = 0.59

cluster_model1 = (len(combined_cluster_x1s)*w)*cluster_frac1/np.sqrt(2.0*np.pi*(sig1**2)) \
            * np.exp(-(data-mu1)**2/2.0/(sig1**2))

cluster_model2 = (len(combined_cluster_x1s)*w)*(1.0 - cluster_frac1)/np.sqrt(2.0*np.pi*(sig2**2)) \
            * np.exp(-(data-mu2)**2/2.0/(sig2**2))

cluster_model = (cluster_model1 + cluster_model2)

quiescent_frac1 = 0.44

quiescent_model1 = (len(quiescent_x1s)*w)*quiescent_frac1/np.sqrt(2.0*np.pi*(sig1**2)) \
            * np.exp(-(data-mu1)**2/2.0/(sig1**2))

quiescent_model2 = (len(quiescent_x1s)*w)*(1.0 - quiescent_frac1)/np.sqrt(2.0*np.pi*(sig2**2)) \
            * np.exp(-(data-mu2)**2/2.0/(sig2**2))

quiescent_model = (quiescent_model1 + quiescent_model2)

sf_frac1 = 0.08


sf_model1 = (len(sf_x1s)*w)*sf_frac1/np.sqrt(2.0*np.pi*(sig1**2)) \
            * np.exp(-(data-mu1)**2/2.0/(sig1**2))

sf_model2 = (len(sf_x1s)*w)*(1.0 - sf_frac1)/np.sqrt(2.0*np.pi*(sig2**2)) \
            * np.exp(-(data-mu2)**2/2.0/(sig2**2))

sf_model = (sf_model1 + sf_model2)


fig,axs = plt.subplots(2,2,figsize = (8,6),sharex=True)
fig.subplots_adjust(hspace=0,wspace=0)
axs[0][0].get_shared_y_axes().join(axs[0][0], axs[0][1])
axs[1][0].get_shared_y_axes().join(axs[1][0], axs[1][1])

axs[0][0].hist(combined_cluster_x1s,density=False,bins=np.arange(-3,3,w),color='k',histtype='step',linewidth=3,label='full cluster')
axs[0][0].plot(data,cluster_model1,color='lightcoral',zorder=0,linewidth=3)
axs[0][0].plot(data,cluster_model2,color='dodgerblue',zorder=0,linewidth=3)
axs[0][0].fill_between(data,cluster_model1,color='lightcoral',alpha=0.5)
axs[0][0].fill_between(data,cluster_model2,color='dodgerblue',alpha=0.5)
axs[0][0].text(.31, .7, r'$\mathbf{f_1 = 0.59^{+0.05}_{-0.05}}$', ha='left', va='top', transform=axs[0][0].transAxes,color = 'coral',fontsize=12,weight='bold')

axs[0][1].hist(r500_x1s,bins=np.arange(-3,3,w),color='slategrey',histtype='step',linewidth=3,label = 'inner cluster')
axs[0][1].plot(data,r500_model1,color='lightcoral',zorder=0,linewidth=3)
axs[0][1].plot(data,r500_model2,color='dodgerblue',zorder=0,linewidth=3)
axs[0][1].fill_between(data,r500_model1,color='lightcoral',alpha=0.5)
axs[0][1].fill_between(data,r500_model2,color='dodgerblue',alpha=0.5)
axs[0][1].text(.02, .74, r'$\mathbf{f_1 = 0.76^{+0.06}_{-0.06}}$', ha='left', va='top', transform=axs[0][1].transAxes,color = 'coral',fontsize=12,weight='bold')

axs[1][0].hist(quiescent_x1s,bins=np.arange(-3,3,w),color='orangered',histtype='step',linewidth=3,label = 'field quiescent')
axs[1][0].plot(data,quiescent_model1,color='lightcoral',zorder=0,linewidth=3)
axs[1][0].plot(data,quiescent_model2,color='dodgerblue',zorder=0,linewidth=3)
axs[1][0].fill_between(data,quiescent_model1,color='lightcoral',alpha=0.5)
axs[1][0].fill_between(data,quiescent_model2,color='dodgerblue',alpha=0.5)
axs[1][0].text(.02, .62, r'$\mathbf{f_1 = 0.44^{+0.03}_{-0.03}}$', ha='left', va='top', transform=axs[1][0].transAxes,color = 'coral',fontsize=12,weight='bold')

axs[1][1].hist(sf_x1s,bins=np.arange(-3,3,w),color='slateblue',histtype='step',linewidth=3,label = 'field star-forming')
axs[1][1].plot(data,sf_model1,color='lightcoral',zorder=0,linewidth=3)
axs[1][1].plot(data,sf_model2,color='dodgerblue',zorder=0,linewidth=3)
axs[1][1].fill_between(data,sf_model1,color='lightcoral',alpha=0.5)
axs[1][1].fill_between(data,sf_model2,color='dodgerblue',alpha=0.5)
axs[1][1].text(.02, .38, r'$\mathbf{f_1 = 0.08^{+0.02}_{-0.01}}$', ha='left', va='top', transform=axs[1][1].transAxes,color = 'coral',fontsize=12,weight='bold')


axs[0][1].yaxis.set_tick_params(labelleft=False)
axs[1][1].yaxis.set_tick_params(labelleft=False)

axs[0][1].set_yticks([])
axs[1][1].set_yticks([])
axs[0][0].set_ylabel('number',fontsize=16)
axs[1][0].set_ylabel('number',fontsize=16)
axs[1][0].set_xlabel(r'$x_1$',fontsize=16)
axs[1][1].set_xlabel(r'$x_1$',fontsize=16)
axs[0][0].set_yticks(np.arange(0,35,5))

axs[0][0].set_ylim(0,31)
axs[0][0].set_xlim(-3.2,3.2)
axs[1][1].set_ylim(0,150)

axs[0][0].legend(fontsize=11)
axs[0][1].legend(fontsize=11,loc=1)
axs[1][0].legend(fontsize=11)
axs[1][1].legend(fontsize=11,loc=1)

axs[0][0].tick_params(axis='y', labelsize=12)
axs[1][0].tick_params(axis='y', labelsize=12)
axs[1][0].tick_params(axis='x', labelsize=12)
axs[1][1].tick_params(axis='x', labelsize=12)

plt.show()
plt.close()