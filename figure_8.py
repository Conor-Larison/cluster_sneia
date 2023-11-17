import matplotlib.pyplot as plt
import numpy as np
import pickle
import corner


r500_pkl_file = open('inner_cluster_mcmc_samples.pkl', 'rb')

r500_sample = np.array(pickle.load(r500_pkl_file))

twor500_pkl_file = open('outer_cluster_mcmc_samples.pkl', 'rb')

twor500_sample = np.array(pickle.load(twor500_pkl_file))

r500_pkl_file.close()
twor500_pkl_file.close()

tlabels = [r"$\alpha$", 
           r"$\beta$",
           r"$M_B$",
           r"$\sigma_{\mathrm{int}}$" ]

fig = plt.figure(figsize=(10,10))

corner.corner(twor500_sample, labels=tlabels, 
                    show_titles=True, title_fmt=".3f",verbose=True,
                    title_kwargs={"fontsize": 14}, label_kwargs={"fontsize": 20},color='limegreen',fig=fig)

corner.corner(r500_sample, labels=tlabels, 
                    show_titles=True, title_fmt=".3f",verbose=True,
                    title_kwargs={"fontsize": 14}, label_kwargs={"fontsize": 20},color='k',fig=fig)


axes = np.array(fig.axes)

r500_alpha = np.percentile(r500_sample[:, 0], [16, 50, 84])
r500_beta = np.percentile(r500_sample[:, 1], [16, 50, 84])
r500_MB = np.percentile(r500_sample[:, 2], [16, 50, 84])
r500_sigma = np.percentile(r500_sample[:, 3], [16, 50, 84])

twor500_alpha = np.percentile(twor500_sample[:, 0], [16, 50, 84])
twor500_beta = np.percentile(twor500_sample[:, 1], [16, 50, 84])
twor500_MB = np.percentile(twor500_sample[:, 2], [16, 50, 84])
twor500_sigma = np.percentile(twor500_sample[:, 3], [16, 50, 84])


axes[0].text(0.5, 1.225,r'$\mathbf{\alpha}$ = $\mathbf{{' + "%.3f" % r500_alpha[1] + '}_{' + "%.3f" % (r500_alpha[0] - r500_alpha[1]) + '}^{+' + "%.3f" % (r500_alpha[2] - r500_alpha[1])  + '}}$', ha='center', va='top', transform=axes[0].transAxes,color = 'k',fontsize=12,weight='bold')
axes[0].text(0.5, 1.11,r'$\mathbf{\alpha}$ = $\mathbf{{' + "%.3f" % twor500_alpha[1] + '}_{' + "%.3f" % (twor500_alpha[0] - twor500_alpha[1]) + '}^{+' + "%.3f" % (twor500_alpha[2] - twor500_alpha[1])  + '}}$', ha='center', va='top', transform=axes[0].transAxes,color = 'limegreen',fontsize=12,weight='bold')
axes[0].set_title('')
axes[5].text(0.5, 1.225,r'$\mathbf{\beta}$ = $\mathbf{{' + "%.3f" % r500_beta[1] + '}_{' + "%.3f" % (r500_beta[0] - r500_beta[1]) + '}^{+' + "%.3f" % (r500_beta[2] - r500_beta[1])  + '}}$', ha='center', va='top', transform=axes[5].transAxes,color = 'k',fontsize=12,weight='bold')
axes[5].text(0.5, 1.11,r'$\mathbf{\beta}$ = $\mathbf{{' + "%.3f" % twor500_beta[1] + '}_{' + "%.3f" % (twor500_beta[0] - twor500_beta[1]) + '}^{+' + "%.3f" % (twor500_beta[2] - twor500_beta[1])  + '}}$', ha='center', va='top', transform=axes[5].transAxes,color = 'limegreen',fontsize=12,weight='bold')
axes[5].set_title('')
axes[10].text(0.5, 1.225,r'$\mathbf{M_B}$ = $\mathbf{{' + "%.3f" % r500_MB[1] + '}_{' + "%.3f" % (r500_MB[0] - r500_MB[1]) + '}^{+' + "%.3f" % (r500_MB[2] - r500_MB[1])  + '}}$',ha='center', va='top', transform=axes[10].transAxes,color = 'k',fontsize=12,weight='bold')
axes[10].text(0.5, 1.11,r'$\mathbf{M_B}$ = $\mathbf{{' + "%.3f" % twor500_MB[1] + '}_{' + "%.3f" % (twor500_MB[0] - twor500_MB[1]) + '}^{+' + "%.3f" % (twor500_MB[2] - twor500_MB[1])  + '}}$', ha='center', va='top', transform=axes[10].transAxes,color = 'limegreen',fontsize=12,weight='bold')
axes[10].set_title('')
axes[15].text(0.5, 1.225,r'$\mathbf{\sigma_{\mathrm{int}}}$ = $\mathbf{{' + "%.3f" % r500_sigma[1] + '}_{' + "%.3f" % (r500_sigma[0] - r500_sigma[1]) + '}^{+' + "%.3f" % (r500_sigma[2] - r500_sigma[1])  + '}}$', ha='center', va='top', transform=axes[15].transAxes,color = 'k',fontsize=12,weight='bold')
axes[15].text(0.5, 1.11,r'$\mathbf{\sigma_{\mathrm{int}}}$ = $\mathbf{{' + "%.3f" % twor500_sigma[1] + '}_{' + "%.3f" % (twor500_sigma[0] - twor500_sigma[1]) + '}^{+' + "%.3f" % (twor500_sigma[2] - twor500_sigma[1])  + '}}$', ha='center', va='top', transform=axes[15].transAxes,color = 'limegreen',fontsize=12,weight='bold')
axes[15].set_title('')


plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

plt.show()
plt.close()