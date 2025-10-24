# netcurve


This repository is used for replicating the results and figures in the paper "(Asymptotically Normal Estimation of Local Latent Network Curvature)[https://arxiv.org/abs/2211.11673]".  This is used in conjunction with the package (lolaR)[https://github.com/SteveJWR/lolaR] which contain the main methods for the paper. 


Within the R directory the following files are included: 

00_functions.R - Base functions for simulations and plotting. 
01a_plot_heatmaps.R - Includes plots for for biases as a function of triangle shape
01b_full_model_estimation.R - Simulation results for latent spaces of different curvature.
02a_snap_network_application.R  - Results for the geometry of citation networks application. 
02b_LANL_sequence.R - Results for the cybersecurity application from LANL clusters
03a_power_simulations_adjacent_sphere.R - Simulation results for testing against a non-constant curvature latent manifold with two different sized adjacent spheres. 
03b_power_simulations_multi_view.R - Simulation results for multi-view network problem
03c_changepoints_simulations_time_series.R - Simulation results for time series changepoints of curvature. 
03d_plot_sims.R - Plotting additional simulation results 
03e_graph_summary_statistics.R - Aggregating summary statistics of simulation results for the appendix. 
04b_misc_plots.R - Additional plots for the appendix. 
04c_variance_theory_plot.R - Theoretical variance as a function of triangle shape heatmaps. 
04d_midpoint_bias_theory.R - Theoretical bias as a function of midpoint misalignment in triangles.
clique_finder.R - A set of functions used for finding cliques within the networks for use in constructing distance matrix estimates. 

