# Environmental SNe Ia
Data from "Environmental Dependence of Type Ia Supernovae in Low-Redshift Galaxy Clusters", [NASA ADS](https://ui.adsabs.harvard.edu/abs/2023arXiv230601088L/abstract)

inner_cluster_data.csv and outer_cluster_data.csv include the SALT3 mB, x1, and c parameter values, distance moduli and Hubble residuals (with _1 referring to Figure 9 and _2 referring to Figure 10), outlier designation from MCMC procedure, host cluster, host cluster redshift (with Hubble diagram version converted to frame of CMB), host cluster r500, projected separation from cluster center, NED Host galaxy name, photometrically-derived estimate for host mass, host or SN redshift used in analysis, and the Host SFR category (Q: quiescent, SF: star-forming, GV: green valley) for our cluster SNe Ia.

sf_field.csv and quiescent_field.csv contain SALT parameter values, distance moduli and Hubble residuals (from Figure 10), host galaxy sSFR and mass measurements, and host redshifts (all spectroscopic, also with Hubble diagram converted values) for SNe Ia in our field samples.

full_cluster.csv contains the data from the table in the appendix of the paper.

The inner_cluster_/outer_cluster_mcmc_samples.csv files contain the samples needed to reproduce the corner plot for Figure 10.

The Python scripts recreate the figures from the paper given the above data. The details for which columns and constraints needed to reproduce the figures are included in these files.
