from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


r500_data = ascii.read('inner_cluster_data.csv')
twor500_data = ascii.read('outer_cluster_data.csv')

quiescent_data = ascii.read('quiescent_field.csv')
sf_data = ascii.read('sf_field.csv')

r500_x1s = r500_data['x1']
r500_cs = r500_data['c']
twor500_x1s = twor500_data['x1']
twor500_cs = twor500_data['c']

quiescent_x1s = quiescent_data['x1']
quiescent_cs = quiescent_data['c']
sf_x1s = sf_data['x1']
sf_cs = sf_data['c']


fig, axs = plt.subplots(2,2,figsize=(7.5,6))

w = 0.5
axs[0][0].hist(r500_x1s,color='k',zorder=5,label=r'inner cluster',histtype='step',density=False,bins=np.arange(-3,3, w),linewidth=2)
axs[0][0].hist(twor500_x1s,color='limegreen',zorder=1,label=r'outer cluster',linestyle='dashdot',histtype='step',density=False,bins=np.arange(-3,3, w),linewidth=2)
axs[0][0].set_xlabel(r'$x_1$',fontsize=10)
axs[0][0].set_ylabel('number',fontsize=10)
axs[0][0].set_xlim(-3,3)
axs[0][0].set_ylim(0,20)
axs[0][0].yaxis.set_major_locator(plt.MaxNLocator(integer=True))
axs[0][0].legend(prop={'size': 9})

w1 = 0.05
axs[0][1].hist(r500_cs,color='black',zorder=5,label=r'inner cluster',histtype='step',density=False,bins=np.arange(-0.3,0.3, w1),linewidth=2)
axs[0][1].hist(twor500_cs,color='limegreen',zorder=1,label=r'outer cluster',linestyle='dashdot',histtype='step',density=False,bins=np.arange(-0.3,0.3, w1),linewidth=2)
axs[0][1].set_xlabel(r'$c$',fontsize=10)
axs[0][1].set_ylabel('number',fontsize=10)
axs[0][1].set_xlim(-0.3,0.3)
axs[0][1].set_ylim(0,20)
axs[0][1].yaxis.set_major_locator(plt.MaxNLocator(integer=True))
axs[0][1].legend(loc=1,prop={'size': 9})

axs[1][0].hist(quiescent_x1s,color='orangered',zorder=5,label='field quiescent',histtype='step',density=False,bins=np.arange(-3,3, w),linewidth=2)
axs[1][0].hist(sf_x1s,color='slateblue',zorder=1,label='field star-forming',linestyle='dashdot',histtype='step',density=False,bins=np.arange(-3,3, w),linewidth=2)
axs[1][0].set_xlabel(r'$x_1$',fontsize=10)
axs[1][0].set_ylabel('number',fontsize=10)
axs[1][0].set_xlim(-3,3)
axs[1][0].set_ylim(0,120)
axs[1][0].yaxis.set_major_locator(plt.MaxNLocator(integer=True))
axs[1][0].legend(prop={'size': 9})

axs[1][1].hist(quiescent_cs,color='orangered',zorder=5,label='field quiescent',histtype='step',density=False,bins=np.arange(-0.3,0.3, w1),linewidth=2)
axs[1][1].hist(sf_cs,color='slateblue',zorder=1,label='field star-forming',linestyle='dashdot',histtype='step',density=False,bins=np.arange(-0.3,0.3, w1),linewidth=2)
axs[1][1].set_xlabel(r'$c$',fontsize=10)
axs[1][1].set_ylabel('number',fontsize=10)
axs[1][1].set_xlim(-0.3,0.3)
axs[1][1].set_ylim(0,120)
axs[1][1].yaxis.set_major_locator(plt.MaxNLocator(integer=True))
axs[1][1].legend(prop={'size':9})

plt.tight_layout()
plt.show()
plt.close()