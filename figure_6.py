import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

r500_data = ascii.read('inner_cluster_data.csv')
twor500_data = ascii.read('outer_cluster_data.csv')

r500_data = r500_data[np.where(r500_data['Host SFR'] == 'Q')]
twor500_data = twor500_data[np.where(twor500_data['Host SFR'] == 'Q')]


r500_x1s = r500_data['x1']
twor500_x1s = twor500_data['x1']

r500_x1errs = r500_data['x1_err']
twor500_x1errs = twor500_data['x1_err']

r500_cluster_sep = r500_data['Projected r/r500']
twor500_cluster_sep = twor500_data['Projected r/r500']

combined_x1s = np.concatenate([r500_x1s,twor500_x1s])


frac1,mu1,sig1,mu2,sig2 = [0.69,-2.03,0.44,0.12,0.54]

data = np.arange(-4,4,0.01)
model1 = frac1/np.sqrt(2.0*np.pi*(sig1**2)) \
         * np.exp(-(data-mu1)**2/2.0/(sig1**2))

model2 = (1.0 - frac1)/np.sqrt(2.0*np.pi*(sig2**2)) \
         * np.exp(-(data-mu2)**2/2.0/(sig2**2))


pop1_color,pop2_color = ['lightcoral','dodgerblue']

pop1alpha = 0.5
pop2alpha = 0.5

fig,axs = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios': [1, 2.5]},figsize=(6,6))
fig.subplots_adjust(hspace=0)

w = 0.5
axs[0].hist(combined_x1s,color='grey',zorder=5,histtype='step',density=False,bins=np.arange(-3, 3, w),linewidth=2)
axs[0].plot(data,len(combined_x1s)*w*model1,color=pop1_color,zorder=0)
axs[0].plot(data,len(combined_x1s)*w*model2,color=pop2_color,zorder=0)
axs[0].axvline(-2.03,linestyle = 'dashed',color=pop1_color,zorder=1)
axs[0].fill_between(data,len(combined_x1s)*w*model1,color=pop1_color,alpha=pop1alpha)
axs[0].axvline(0.12,linestyle = 'dashed',color=pop2_color,alpha=pop2alpha,zorder=1)
axs[0].fill_between(data,len(combined_x1s)*w*model2,color=pop2_color,alpha=pop2alpha)

axs[1].errorbar(r500_x1s,r500_cluster_sep,xerr=r500_x1errs,color='k',fmt='d')
axs[1].errorbar(twor500_x1s,twor500_cluster_sep,xerr=twor500_x1errs,color='limegreen',fmt='o')
axs[1].set_xlabel(r'$x_1$',fontsize=16)
x = np.arange(-0.1,2.1,0.1)
y1 = np.array([-2.03] * len(x))
y2 = np.array([0.3] * len(x))
axs[1].axvline(-2.03,linestyle = 'dashed',color=pop1_color)
axs[1].axvline(0.12,linestyle = 'dashed',color=pop2_color,alpha=pop2alpha)
axs[1].set_ylabel(r'projected  $\mathrm{r}/\mathrm{r}_{500}$',fontsize=14)

axs[0].set_ylabel('number',fontsize=14)
axs[1].set_xlim(-3.5,2)
axs[1].set_ylim(-0.1,2.15)

plt.show()
plt.close()