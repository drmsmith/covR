library(deSolve)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(scico)
library(epiR)
library(RColorBrewer)

source("ODEs.R")
source("functions.R")
source("pars_states.R")


################
### PIPELINE ###
################
### run code for various analyses

### ANALYSIS1: for baseline parameter values, vary values of tau over range and evaluate impacts on dynamics
### ANALYSIS2: Monte Carlo simulations across broad parameter space, evaluate indicators
### ANALYSIS3: Case studies for E. coli and MRSA in a selection of healthcare facilities

### NB: always save data to folders 1+ levels above where they are later stored/called from (pour assurer que ca ne s'ecrase pas)

##################
### ANALYSIS 1 ###
##################

### (1.1) loop through pandemic impacts with identical values, solve ODEs and save output to folder /data/

f_dynamics_loop_tau(model_loop = ODEs_covR, pars_loop = pars_m, states_loop = states_m_num, time_out = 180, list_impacts = list_impacts, save_data = 'YES')

### (1.2) loop through sets of tau with random values

# first, generate and save vectors of random tau
df_tau_random = data.frame(handrub = runif(1000,0,0.9999),
                           masks = runif(1000,0,0.9999),
                           prophylaxis = runif(1000,0,0.9999),
                           distancing = runif(1000,0,0.9999),
                           disorg = runif(1000,0,0.9999),
                           comm_denom = runif(1000,0,0.9999),
                           a_surge = runif(1000,0,0.9999),
                           chi = runif(1000,0,0.9999),
                           covid_stay = runif(1000,0,0.9999),
                           adm_reduc = runif(1000,0,0.9999))

# save(df_tau_random, file = "data/df_tau_random_vals.Rdata")


# second, run a loop sampling the desired tau and running ODEs
f_dynamics_loop_random_tau(model_loop = ODEs_covR, pars_loop = pars_m, states_loop = states_m_num, time_out = 180, n_parsets = 100, save_data = 'YES')



##################
### ANALYSIS 2 ###
##################

### (2.1)  generate parameter sets for Monte Carlo simulation
pars_analysis2 = f_generate_parsets_uncertain_facility_bacteria(500, pars_m)

#save(pars_analysis2, file = "data/parsets_analysis2.Rdata")

### (2.2)  run MC simulations and save individually
f_runODEs_calculateIndicators(model_loop = ODEs_covR,
                              value_tau_loop = 0.5,
                              time_out_loop = 180,
                              n_parset_start = 289,
                              n_parset_end = 289,
                              filepath_parsets = "data/outputs_analysis2/parsets_analysis2.Rdata", 
                              save_data = 'YES')


### (2.3) combine data into a single dataframe
# load data
filepath_indicators = paste0("data/outputs_analysis2/simu_indicators_raw/")
files_indicators = list.files(filepath_indicators); 
# concatenate in list
list_indicators = list()
for(n_row in 1:length(files_indicators)){
  list_indicators[[n_row]] = loadRData(paste0(filepath_indicators, files_indicators[n_row]))
}
# combine
df_indicators = do.call(rbind.data.frame, list_indicators)
# save(df_indicators, file = paste0("data/indicators_summarized_raw.Rdata"))

### (2.4) calculate "deltas": difference between indicators for each tau and baseline simulations without any pandemic impacts
f_indicator_deltas(filepath_indicators_summarized = "data/outputs_analysis2/indicators_summarized_raw.Rdata", save_data = "YES")

### (2.5) bootstrap "deltas" to assess whether sufficient # of simulations
f_bootstrap_deltas(n_bootstraps = 50, save_data = 'YES')

### (2.6) PRCC for key outcomes
pars_analysis2 = loadRData("data/outputs_analysis2/parsets_analysis2.Rdata")
df_indicators = loadRData("data/outputs_analysis2/indicators_summarized_raw.Rdata")
f_prcc(m_pars = pars_analysis2, 
       df_outputs = df_indicators, 
       vec_pars_varied = vec_pars_varied, 
       vec_tau = as.character(df_pandI_withNone$par), 
       vec_outcomes = c("V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe", "V_inc_cumul_all",
                        "Cr_pd", "R_rate_pd",
                        "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog", "Cr_inc_cumul_all",
                        "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe", "Tr_inc_cumul_all"), 
       save.data = "YES")


##################
### ANALYSIS 3 ###
##################

### (3.1)  generate parameter sets for Monte Carlo simulation
pars_analysis3_rehab = f_generate_parsets_casestudies_setting(500, pars_m, "rehab")
pars_analysis3_geriatric = f_generate_parsets_casestudies_setting(500, pars_m, "geriatric")
pars_analysis3_paediatric = f_generate_parsets_casestudies_setting(500, pars_m, "paediatric")
pars_analysis3_organized = f_generate_parsets_casestudies_scenario(500, pars_m, "organized")
pars_analysis3_intermediate = f_generate_parsets_casestudies_scenario(500, pars_m, "intermediate")
pars_analysis3_overwhelmed = f_generate_parsets_casestudies_scenario(500, pars_m, "overwhelmed")
pars_analysis3_s_aureus = f_generate_parsets_casestudies_species(500, pars_m, "s_aureus")
pars_analysis3_e_coli = f_generate_parsets_casestudies_species(500, pars_m, "e_coli")

pars_analysis3_rehab_organized_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_rehab, pars_analysis3_organized, pars_analysis3_s_aureus)
pars_analysis3_rehab_intermediate_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_rehab, pars_analysis3_intermediate, pars_analysis3_s_aureus)
pars_analysis3_rehab_overwhelmed_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_rehab, pars_analysis3_overwhelmed, pars_analysis3_s_aureus)

pars_analysis3_rehab_organized_e_coli = f_parsets_setting_scenario_species(pars_analysis3_rehab, pars_analysis3_organized, pars_analysis3_e_coli)
pars_analysis3_rehab_intermediate_e_coli = f_parsets_setting_scenario_species(pars_analysis3_rehab, pars_analysis3_intermediate, pars_analysis3_e_coli)
pars_analysis3_rehab_overwhelmed_e_coli = f_parsets_setting_scenario_species(pars_analysis3_rehab, pars_analysis3_overwhelmed, pars_analysis3_e_coli)

pars_analysis3_geriatric_organized_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_geriatric, pars_analysis3_organized, pars_analysis3_s_aureus)
pars_analysis3_geriatric_intermediate_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_geriatric, pars_analysis3_intermediate, pars_analysis3_s_aureus)
pars_analysis3_geriatric_overwhelmed_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_geriatric, pars_analysis3_overwhelmed, pars_analysis3_s_aureus)

pars_analysis3_geriatric_organized_e_coli = f_parsets_setting_scenario_species(pars_analysis3_geriatric, pars_analysis3_organized, pars_analysis3_e_coli)
pars_analysis3_geriatric_intermediate_e_coli = f_parsets_setting_scenario_species(pars_analysis3_geriatric, pars_analysis3_intermediate, pars_analysis3_e_coli)
pars_analysis3_geriatric_overwhelmed_e_coli = f_parsets_setting_scenario_species(pars_analysis3_geriatric, pars_analysis3_overwhelmed, pars_analysis3_e_coli)

pars_analysis3_paediatric_organized_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_paediatric, pars_analysis3_organized, pars_analysis3_s_aureus)
pars_analysis3_paediatric_intermediate_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_paediatric, pars_analysis3_intermediate, pars_analysis3_s_aureus)
pars_analysis3_paediatric_overwhelmed_s_aureus = f_parsets_setting_scenario_species(pars_analysis3_paediatric, pars_analysis3_overwhelmed, pars_analysis3_s_aureus)

pars_analysis3_paediatric_organized_e_coli = f_parsets_setting_scenario_species(pars_analysis3_paediatric, pars_analysis3_organized, pars_analysis3_e_coli)
pars_analysis3_paediatric_intermediate_e_coli = f_parsets_setting_scenario_species(pars_analysis3_paediatric, pars_analysis3_intermediate, pars_analysis3_e_coli)
pars_analysis3_paediatric_overwhelmed_e_coli = f_parsets_setting_scenario_species(pars_analysis3_paediatric, pars_analysis3_overwhelmed, pars_analysis3_e_coli)

save(pars_analysis3_rehab_organized_s_aureus, file = "data/parsets_analysis3_rehab_organized_s_aureus.Rdata")
save(pars_analysis3_rehab_intermediate_s_aureus, file = "data/parsets_analysis3_rehab_intermediate_s_aureus.Rdata")
save(pars_analysis3_rehab_overwhelmed_s_aureus, file = "data/parsets_analysis3_rehab_overwhelmed_s_aureus.Rdata")
save(pars_analysis3_rehab_organized_e_coli, file = "data/parsets_analysis3_rehab_organized_e_coli.Rdata")
save(pars_analysis3_rehab_intermediate_e_coli, file = "data/parsets_analysis3_rehab_intermediate_e_coli.Rdata")
save(pars_analysis3_rehab_overwhelmed_e_coli, file = "data/parsets_analysis3_rehab_overwhelmed_e_coli.Rdata")
save(pars_analysis3_geriatric_organized_s_aureus, file = "data/parsets_analysis3_geriatric_organized_s_aureus.Rdata")
save(pars_analysis3_geriatric_intermediate_s_aureus, file = "data/parsets_analysis3_geriatric_intermediate_s_aureus.Rdata")
save(pars_analysis3_geriatric_overwhelmed_s_aureus, file = "data/parsets_analysis3_geriatric_overwhelmed_s_aureus.Rdata")
save(pars_analysis3_geriatric_organized_e_coli, file = "data/parsets_analysis3_geriatric_organized_e_coli.Rdata")
save(pars_analysis3_geriatric_intermediate_e_coli, file = "data/parsets_analysis3_geriatric_intermediate_e_coli.Rdata")
save(pars_analysis3_geriatric_overwhelmed_e_coli, file = "data/parsets_analysis3_geriatric_overwhelmed_e_coli.Rdata")
save(pars_analysis3_paediatric_organized_s_aureus, file = "data/parsets_analysis3_paediatric_organized_s_aureus.Rdata")
save(pars_analysis3_paediatric_intermediate_s_aureus, file = "data/parsets_analysis3_paediatric_intermediate_s_aureus.Rdata")
save(pars_analysis3_paediatric_overwhelmed_s_aureus, file = "data/parsets_analysis3_paediatric_overwhelmed_s_aureus.Rdata")
save(pars_analysis3_paediatric_organized_e_coli, file = "data/parsets_analysis3_paediatric_organized_e_coli.Rdata")
save(pars_analysis3_paediatric_intermediate_e_coli, file = "data/parsets_analysis3_paediatric_intermediate_e_coli.Rdata")
save(pars_analysis3_paediatric_overwhelmed_e_coli, file = "data/parsets_analysis3_paediatric_overwhelmed_e_coli.Rdata")

### (3.2) run simulations
vec_beta_S = c(0.76, 1.28, 2.4)
vec_t_policy = c(7, 21, 35)
vec_scenarios = c("organized", "intermediate", "overwhelmed")

# (3.2.1) rehabilitation hospital, S. aureus
for(beta_S_i in vec_beta_S){
  for(t_policy_i in vec_t_policy){
    for(scenario_i in vec_scenarios){
      
      print(paste0("for S. aureus in rehab, on scenario ", scenario_i, " with beta_S = ", beta_S_i, " and t_policy = ", t_policy_i))
      
      f_case_studies_dynamics(model_loop = ODEs_covR, 
                              time_out_loop = 180,
                              setting = "rehab", 
                              scenario = scenario_i,
                              species = "s_aureus", 
                              in_beta_S = beta_S_i,
                              in_t_policy = t_policy_i,
                              n_parset_start = 1, 
                              n_parset_end = 100, 
                              save_data = 'YES')
    }
  }
}


# (3.2.2) rehabilitation hospital, E. coli
for(beta_S_i in vec_beta_S){
  for(t_policy_i in vec_t_policy){
    for(scenario_i in vec_scenarios){
      
      print(paste0("for E. coli in rehab, on scenario ", scenario_i, " with beta_S = ", beta_S_i, " and t_policy = ", t_policy_i))
      
      f_case_studies_dynamics(model_loop = ODEs_covR, 
                              time_out_loop = 180,
                              setting = "rehab", 
                              scenario = scenario_i,
                              species = "e_coli", 
                              in_beta_S = beta_S_i,
                              in_t_policy = t_policy_i,
                              n_parset_start = 1, 
                              n_parset_end = 100, 
                              save_data = 'YES')
    }
  }
}

# (3.2.3) geriatric ward, S. aureus
for(beta_S_i in vec_beta_S){
  for(t_policy_i in vec_t_policy){
    for(scenario_i in vec_scenarios){
      
      print(paste0("for S. aureus in geriatrc, on scenario ", scenario_i, " with beta_S = ", beta_S_i, " and t_policy = ", t_policy_i))
      
      f_case_studies_dynamics(model_loop = ODEs_covR, 
                              time_out_loop = 180,
                              setting = "geriatric", 
                              scenario = scenario_i,
                              species = "s_aureus", 
                              in_beta_S = beta_S_i,
                              in_t_policy = t_policy_i,
                              n_parset_start = 1, 
                              n_parset_end = 100, 
                              save_data = 'YES')
    }
  }
}

# (3.2.4) geriatric ward, E. coli
for(beta_S_i in vec_beta_S){
  for(t_policy_i in vec_t_policy){
    for(scenario_i in vec_scenarios){
      
      print(paste0("for E. coli in geriatric, on scenario ", scenario_i, " with beta_S = ", beta_S_i, " and t_policy = ", t_policy_i))
      
      f_case_studies_dynamics(model_loop = ODEs_covR, 
                              time_out_loop = 180,
                              setting = "geriatric", 
                              scenario = scenario_i,
                              species = "e_coli", 
                              in_beta_S = beta_S_i,
                              in_t_policy = t_policy_i,
                              n_parset_start = 1, 
                              n_parset_end = 100, 
                              save_data = 'YES')
    }
  }
}

# (3.2.5) paediatric ward, S. aureus
for(beta_S_i in vec_beta_S){
  for(t_policy_i in vec_t_policy){
    for(scenario_i in vec_scenarios){
      
      print(paste0("for S. aureus in paediatric, on scenario ", scenario_i, " with beta_S = ", beta_S_i, " and t_policy = ", t_policy_i))
      
      f_case_studies_dynamics(model_loop = ODEs_covR, 
                              time_out_loop = 180,
                              setting = "paediatric", 
                              scenario = scenario_i,
                              species = "s_aureus", 
                              in_beta_S = beta_S_i,
                              in_t_policy = t_policy_i,
                              n_parset_start = 1, 
                              n_parset_end = 100, 
                              save_data = 'YES')
    }
  }
}

# (3.2.6) paediatric ward, E. coli
for(beta_S_i in vec_beta_S){
  for(t_policy_i in vec_t_policy){
    for(scenario_i in vec_scenarios){
      
      print(paste0("for E. coli in paediatric, on scenario ", scenario_i, " with beta_S = ", beta_S_i, " and t_policy = ", t_policy_i))
      
      f_case_studies_dynamics(model_loop = ODEs_covR, 
                              time_out_loop = 180,
                              setting = "paediatric", 
                              scenario = scenario_i,
                              species = "e_coli", 
                              in_beta_S = beta_S_i,
                              in_t_policy = t_policy_i,
                              n_parset_start = 1, 
                              n_parset_end = 100, 
                              save_data = 'YES')
    }
  }
}
