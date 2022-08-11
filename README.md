# covR
An evaluation of how COVID-19 responses impact antibiotic resistance, built upon a deterministic, compartmental transmission model describing co-circulation of SARS-CoV-2 and competing strains of commensal bacteria in a healthcare facility population.

This R code supports Smith et al., (2022). Collateral impacts of pandemic COVID-19 drive the nosocomial spread of antibiotic resistance.

All code was developed, tested and run using Rv3.6.0

# about
Ordinary differential equations are written, and functions are provided to run numerical simulations and calculate epidemiological indicators. Three main sets of analyses are conducted, representing three main sets of numerical simulation. Figures are generated.

Analysis 1 describes model simulations over a baseline parameter set (with parameters for SARS-CoV-2 transmission and healthcare facility characteristics largely derived from the hospital outbreak described by Shirreff et al., Emerging Infect Dis 2022). 

Analysis 2 describes Monte Carlo simulations accounting for broad parameter uncertainty across bacterial parameters ("generic MRB") and healthcare facility parameters ("generic hospitals"). 

Analysis 3 describes Monte Carlo simulations applied to case studies of specific bacterial species, hospital wards and COVID-19 response scenarios. 

# Rproject
* (1) covR.Rproj 
  * Files are associated with an R project (for stable working directory)

# model
* (2) ODEs.R:
  *  File containing ODEs

# functions
* (3) functions.R
  * contains all functions for simulation and analysis

# baseline parameter values
* (4) pars_states.R
  * contains a range of parameters and objects used in simulations and analysis, including baseline model parameters and initial state variables for numerical ODE integration, and vectors used for data management and plotting 

# main analysis pipeline
* (5) pipeline_analyses.R 
  * file for execution of all main analyses, including parameter sampling, ODE integration, calculation of epidemiological indicators, bootstrap resampling of  indicators, calculation of partial rank correlation coefficients

# figure rendering pipelines
* (6.1) pipeline_figures_analysis1.R
* (6.2) pipeline_figures_analysis1.R
* (6.3) pipeline_figures_analysis1.R
  * files for rendering figures for respective analyses (described above); it is necessary for users to create a sub-directory "plots" in which to save figures rendered
  
# data
* (7.1) outputs_analysis1/
  * parameter values and simulation outputs for analysis 1
* (7.2) outputs_analysis2/
  * parameter values and simulation outputs for analysis 2, as well as summarized final indicators calculated from simulation outputs
* (7.3) outputs_analysis3/
  * parameter values and simulations outputs for analyis 3 (raw data are only provided for t_policy_in = 21 days, as space exceeded), as well as summarized final indicators

# contact
David Smith \
Institut Pasteur / Inserm / UVSQ \
david.smith@pasteur.fr \
davidrobertmundysmith@gmail.com
