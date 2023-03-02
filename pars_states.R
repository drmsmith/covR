#############################
### PARAMETERS AND STATES ###
#############################

######################
### Organizational ###
######################

### Prevalence categories from "plottable" data
vec_comparts_PA = c("S_pa", "E_pa", "I_pa", "R_pa", "U_pa", "Cs_pa", "Cr_pa", "Total_pa")
vec_comparts_PE = c("S_pe", "E_pe", "I_pe", "R_pe", "U_pe", "Cs_pe", "Cr_pe", "Total_pe")

### Names of prevalence columns from ODE output
cols_prevalence = c("S_U_pa", "S_Cs_pa", "S_Cr_pa", 
                    "E_U_pa", "E_Cs_pa", "E_Cr_pa", 
                    "I_U_pa", "I_Cs_pa", "I_Cr_pa", 
                    "R_U_pa", "R_Cs_pa", "R_Cr_pa",
                    "S_U_pe", "S_Cs_pe", "S_Cr_pe",
                    "E_U_pe", "E_Cs_pe", "E_Cr_pe", 
                    "I_U_pe", "I_Cs_pe", "I_Cr_pe", 
                    "R_U_pe", "R_Cs_pe", "R_Cr_pe",
                    "SL")

### Names of incidence columns from ODE output
cols_incidence = c('incS_pa_pa', 'incS_pe_pa', 'incS_pa_pe', 'incS_pe_pe',
                   'incCs_pa_pa', 'incCs_pe_pa', 'incCs_pa_pe', 'incCs_pe_pe',
                   'incCr_pa_pa', 'incCr_pe_pa', 'incCr_pa_pe', 'incCr_pe_pe', 'incCr_endog')

### Names of metadata columns from ODE output
cols_meta = c("abx", "kappa_patients", "kappa_staff", "hh", "staffing_ratio", "adm", "adm_R", "foi_S", "foi_Cr", "patientdays", "staffdays")

##################
### PARAMETERS ###
##################

#########################
### m_: Generic means ###
#########################

### for the default generic bacteria

# demography
m_mu = 1/80 # Rene Muret
m_f_Cs = 0.15 
m_f_Cr = 0.05
m_f_U = 1 - m_f_Cs - m_f_Cr 

# population size parameters
m_Nbeds = 350 # Rene Muret
m_Nhcw = 529 # Rene Muret

# contact rates
m_kappa_pa_pa = 5
m_kappa_pa_pe = 15

m_kappa_pe_pe = 10
m_kappa_pe_pa = m_kappa_pa_pe*(m_Nbeds/m_Nhcw)

# nosocomial parameters
m_rho = 3 #
m_hyg = 0.4 #

# antibiotics
m_a = 0.1 # previously 0.2 (before switch to Rene Muret)
m_theta = 0.25 #
m_r_S = 0.1 #
m_r_R = 0.9 #

# ARB epi pars (other than transmission)
m_gamma_Cs = 0.03 #
m_cost = 0.1 #
m_gamma_Cr = m_gamma_Cs*(1+m_cost)
m_alpha = 0.01 #

# SARS-CoV-2 epi pars (other than transmission)
m_upsilon = 1/7
m_eta = 1/5

# transmission rates
m_beta_S = 1.28 # Shirreff
m_beta_Cs = 0.2 # 
m_beta_Cr = 0.2 #

m_pi_S = m_beta_S/(m_kappa_pa_pa+m_kappa_pa_pe+m_kappa_pe_pe+m_kappa_pe_pa)

m_pi_Cs = m_beta_Cs/(m_kappa_pa_pa+m_kappa_pe_pa)
m_pi_Cr = m_beta_Cr/(m_kappa_pa_pa+m_kappa_pe_pa)


# sick leave
m_upsilon_SL = 1/7 # (assumed)
m_prop_symptomatic = 1-0.475

# dynamic covid impacts
m_disorg = 0
m_comm_denom = 0 
m_a_surge = 0 
m_chi = 0 # default none
m_covid_stay = 0
m_adm_reduc = 0

# policy covid impacts
m_handrub = 0 # m_ipc_C = 0
m_masks = 0 # m_ipc_S = 0
m_prophylaxis = 0 # m_a_proph = 0
m_distancing = 0 # m_lockd = 0

# when are policy parameters introduced? and how long is the burn-in period?
m_t_policy = 21
m_t_burnin = 7

### put together in a vector
pars_m <- c(
  
  # demography
  mu = m_mu, f_U = m_f_U, f_Cs = m_f_Cs, f_Cr = m_f_Cr, Nbeds = m_Nbeds, Nhcw = m_Nhcw,
  
  # contact rates
  kappa_pa_pa = m_kappa_pa_pa, kappa_pa_pe = m_kappa_pa_pe, 
  kappa_pe_pe = m_kappa_pe_pe, kappa_pe_pa = m_kappa_pe_pa,
  
  # nosocomial parameters
  rho = m_rho, hyg = m_hyg,
  
  # antibiotics
  a = m_a, theta = m_theta, r_S = m_r_S, r_R = m_r_R,
  
  # ARB epi pars (other than transmission)
  gamma_Cs = m_gamma_Cs, cost = m_cost, gamma_Cr = m_gamma_Cr, alpha = m_alpha,
  
  # SARS-CoV-2 epi pars (other than transmission)
  upsilon = m_upsilon, eta = m_eta,
  
  # overall transmission rates
  beta_S = m_beta_S, beta_Cs = m_beta_Cs, beta_Cr = m_beta_Cr,
  
  # Transmission rates per contact = overall transmission rate / sum of all contacts
  pi_S = m_pi_S, pi_Cs = m_pi_Cs, pi_Cr = m_pi_Cr,
  
  # sick leave
  chi = m_chi, upsilon_SL = m_upsilon_SL, prop_symptomatic = m_prop_symptomatic,
  
  # dynamic covid impacts
  disorg = m_disorg, comm_denom = m_comm_denom, a_surge = m_a_surge, covid_stay = m_covid_stay, adm_reduc = m_adm_reduc,
  
  # policy covid impacts
  handrub = m_handrub, masks = m_masks, prophylaxis = m_prophylaxis, distancing = m_distancing,
  t_policy = m_t_policy,
  t_burnin = m_t_burnin
)


##############
### STATES ###
##############

states_m_num = c(S_U_pa = m_Nbeds/3, S_Cs_pa = m_Nbeds/3, S_Cr_pa = m_Nbeds/3, E_U_pa = 0, E_Cs_pa = 0, E_Cr_pa = 0, I_U_pa = 0, I_Cs_pa = 0, I_Cr_pa = 0, R_U_pa = 0, R_Cs_pa = 0, R_Cr_pa = 0,
                 S_U_pe = m_Nhcw/3, S_Cs_pe = m_Nhcw/3, S_Cr_pe = m_Nhcw/3, E_U_pe = 0, E_Cs_pe = 0, E_Cr_pe = 0, I_U_pe = 0, I_Cs_pe = 0, I_Cr_pe = 0, R_U_pe = 0, R_Cs_pe = 0, R_Cr_pe = 0,
                 SL = 0,
                 incS_pa_pa = 0, incS_pe_pa = 0, incS_pa_pe = 0, incS_pe_pe = 0,
                 incCs_pa_pa = 0, incCs_pe_pa = 0, incCs_pa_pe = 0, incCs_pe_pe = 0,
                 incCr_pa_pa = 0, incCr_pe_pa = 0, incCr_pa_pe = 0, incCr_pe_pe = 0, incCr_endog = 0,
                 abx = 0, kappa_patients = 0, kappa_staff = 0, hh = 0, staffing_ratio = 0, adm = 0, adm_R = 0, foi_S = 0, foi_Cr = 0, patientdays = 0, staffdays = 0)


########################
### PANDEMIC IMPACTS ###
########################

### Vectors of combos and groups of pandemic impacts

vec_impacts_abx = c("prophylaxis", "a_surge")
vec_impacts_contact = c("distancing", "disorg")
vec_impacts_ipc = c("handrub", "masks")
vec_impacts_disease = c("chi", "covid_stay")
vec_impacts_admission = c("comm_denom", "adm_reduc")
vec_impacts_policy = c("prophylaxis", "distancing", "handrub", "masks")
vec_impacts_caseload = c("a_surge", "chi", "disorg", "comm_denom", "adm_reduc", "covid_stay")
vec_impacts_all = c("prophylaxis", "distancing", "handrub", "masks", "a_surge", "chi", "disorg", "comm_denom", "adm_reduc", "covid_stay")

list_impacts = list("prophylaxis", "distancing", "handrub", "masks", "a_surge", "chi", "disorg", "comm_denom", "adm_reduc", "covid_stay",
                   vec_impacts_abx, vec_impacts_contact, vec_impacts_ipc, vec_impacts_disease, vec_impacts_admission, 
                   vec_impacts_policy, vec_impacts_caseload, vec_impacts_all)

vec_impactnames = c("prophylaxis", "distancing", "handrub", "masks", "a_surge", "sickleave", "disorg", "comm_denom", "adm_reduc", "covid_stay",
                    "combo_antibiotics", "combo_contact", "combo_ipc", "combo_disease", "combo_admission", "combo_policy", "combo_caseload", "combo_all")

names(list_impacts) <- vec_impactnames

# vec_pars_and_combos = c(vec_pandI_parameters, vec_combonames)


### Data frame of characteristics of different pandemic impacts

# policy responses
df_pandI_prophylaxis = data.frame(impact = 'COVID-19 prescribing', par = 'prophylaxis', type = 'policy', category = "abx")
df_pandI_handrub = data.frame(impact = 'hand hygiene', par = 'handrub', type = 'policy', category = "ipc")
df_pandI_masks = data.frame(impact = 'universal masking', par = 'masks', type = 'policy', category = "ipc")
df_pandI_distancing = data.frame(impact = 'patient lockdown', par = 'distancing', type = 'policy', category = "contact")

# dynamic responses
df_pandI_disorg = data.frame(impact = 'care disorganization', par = 'disorg', type = 'burden', category = "contact")
df_pandI_a_surge = data.frame(impact = 'abandoned stewardship', par = 'a_surge', type = 'burden', category = "abx")
df_pandI_comm_denom = data.frame(impact = 'sicker casemix', par = 'comm_denom', type = 'burden', category = "admission")
df_pandI_sickleave = data.frame(impact = 'staff sick leave', par = 'sickleave', type = 'burden', category = "disease")
df_pandI_covid_stay = data.frame(impact = 'COVID-19 stays', par = 'covid_stay', type = 'burden', category = "disease")
df_pandI_adm_reduc = data.frame(impact = 'reduced admission', par = 'adm_reduc', type = 'burden', category = "admission")

# combos
df_pandI_combo_admission = data.frame(impact = 'admission', par = 'combo_admission', type = 'burden', category = "admission")
df_pandI_combo_all = data.frame(impact = 'all responses', par = 'combo_all', type = 'mixed', category = "mixed")
df_pandI_combo_abx = data.frame(impact = 'antibiotics', par = 'combo_antibiotics', type = 'mixed', category = "abx")
df_pandI_combo_caseload = data.frame(impact = 'caseload', par = 'combo_caseload', type = 'burden', category = "mixed")
df_pandI_combo_contact = data.frame(impact = 'contact', par = 'combo_contact', type = 'mixed', category = "contact")
df_pandI_combo_disease = data.frame(impact = 'disease', par = 'combo_disease', type = 'burden', category = "disease")
df_pandI_combo_ipc = data.frame(impact = 'IPC', par = 'combo_ipc', type = 'policy', category = "ipc")
df_pandI_combo_policy = data.frame(impact = 'policy change', par = 'combo_policy', type = 'policy', category = "mixed")


### all together
df_pandI = rbind(df_pandI_prophylaxis,df_pandI_handrub,df_pandI_masks,df_pandI_distancing,
                 df_pandI_sickleave,df_pandI_disorg,df_pandI_a_surge,df_pandI_comm_denom,df_pandI_covid_stay,df_pandI_adm_reduc,
                 df_pandI_combo_admission, df_pandI_combo_all, df_pandI_combo_abx, df_pandI_combo_caseload,
                 df_pandI_combo_contact,df_pandI_combo_disease, df_pandI_combo_ipc, df_pandI_combo_policy
                 )

vec_pandI_parameters = levels(df_pandI$par)

### all together with "none"
df_pandI_withNone = rbind(df_pandI,
                          data.frame(impact = "none", par = "none", type = "none", category = "none"))

#########################
### PRCC: pars varied ###
#########################

vec_pars_varied = c("mu", "kappa_pa_pa", "kappa_pa_pe", "kappa_pe_pe", "rho", "hyg", "a", "Nbeds", "Nhcw", 
                    "gamma_Cs", "cost", "r_S", "r_R", "theta", "alpha", "f_Cs", "f_Cr", "pi_S", "pi_Cr")

################
### PLOTTING ###
################

### set global themes, colours, etc.
th <- theme_classic()

col_pat_staff = c("#FEA200", "#02AFAE", "#B20093")
col_incidence = brewer.pal(n = 8, name = "Dark2")
col_policy_surge = c("#762a83", "#5aae61")
col_strains = c("#FF040A", "#0264B8")
col_staffing = pal_nejm()(1)
col_hh = pal_nejm()(2)[2]
col_abx = "#8C564B"
col_admission = pal_nejm()(5)[5]
col_foi = pal_nejm()(8)[8]
col_pathogen = c("#d95f02", "#e7298a")
col_species = c("#FFD200", "#FE492C")
col_settings = c("#FFD400", "#01A8B5", "#FF492A")
col_settings_scenario = c("#FFF594", "#FFD400", "#C8A803", #yellow
                          "#01E2F4", "#01A8B5", "#007A84", #cyan
                          "#FF877F", "#FF492A", "#CE3922") #red

cols_nejm_paired = pal_nejm()(5)[c(1,1,2,2,3,3,4,4,5,5)]
cols_10_quali =  c("#f94144", "#f3722c", "#f8961e", "#f9844a", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1")
cols_10_pairs = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
cols_policy_surge = c("#E4AA00", "#C00000", "#767171")
cols_sensitivity = rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))
cols_sensitivity = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

alpha_vline = 0.3
alpha_range_tau = c(0.2,1)

labels_5combos = c(expression(paste("antibiotics (", tau[cp]," and ", tau[as],")")), 
                  expression(paste("contact (", tau[pl]," and ", tau[cd],")")), 
                  expression(paste("IPC (", tau[um]," and ", tau[hh],")")), 
                  expression(paste("disease (", tau[cs]," and ", tau[ss],")")),
                  expression(paste("admission (", tau[ra]," and ", tau[sc],")")))

labels_all10 = c(expression(paste("abandoned stewardship (", tau[as],")")), 
                 expression(paste("COVID-19 prescribing (", tau[cp],")")), 
                 expression(paste("care disorganization (", tau[cd],")")), 
                 expression(paste("patient lockdown (", tau[pl],")")),
                 expression(paste("universal masking (", tau[um], ")")), 
                 expression(paste("hand hygiene (", tau[hh], ")")), 
                 expression(paste("COVID-19 stays (", tau[cs],")")),
                 expression(paste("staff sick leave (", tau[ss],")")),
                 expression(paste("reduced admissions (", tau[ra],")")),
                 expression(paste("sicker casemix (", tau[sc],")")))

vec_impact_clean_order = c("none",
                           "abandoned stewardship", "COVID-19 prescribing",
                           "care disorganization", "patient lockdown",
                           "universal masking", "hand hygiene",
                           "COVID-19 stays", "staff sick leave",
                           "reduced admission", "sicker casemix",
                           "antibiotics", "contact","IPC","disease","admission",
                           "policy change", "caseload",
                           "all responses")

### case studies
size_text = 3

vec_wards_clean = c("short-stay ward\n(geriatric)", "general ward\n(paediatric)", "rehabilitation ward\n(geriatric)")
vec_species_clean = c("MRSA", "ESBL-E. coli")
vec_wards_scenarios_clean = c("short-stay ward\n(geriatric) (organized)", "short-stay ward\n(geriatric) (intermediate)", "short-stay ward\n(geriatric) (overwhelmed)",
                              "general ward\n(paediatric) (organized)", "general ward\n(paediatric) (intermediate)", "general ward\n(paediatric) (overwhelmed)",
                              "rehabilitation ward\n(geriatric) (organized)", "rehabilitation ward\n(geriatric) (intermediate)", "rehabilitation ward\n(geriatric) (overwhelmed)")
