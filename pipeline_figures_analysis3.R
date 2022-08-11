library(deSolve)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(cowplot)
library(RColorBrewer)
library(ggallin)

source("ODEs.R")
source("functions.R")
source("pars_states.R")

#################
### LOAD DATA ###
#################
# create dataframes grouping all species and settings, for the 6 groups of outcomes simulated
# (1.) indicators no covid
# (2.) indicators covid
# (3.) indicators deltas
# (4.) metadata no covid
# (5.) metadata covid
# (6.) metadata deltas

list_casestudy_indicators_nocovid = list()
list_casestudy_indicators_covid = list()
list_casestudy_indicators_deltas = list()
list_casestudy_metadata_nocovid = list()
list_casestudy_metadata_covid = list()
list_casestudy_metadata_deltas = list()

qounter = 0

### NB: for space and size reasons, raw data only uploaded for baseline t_policy_i = 21
### users are invited to run simulations for 7, 14, 28 and 35 using parameters and code provided
#vec_t_policy_i = c(7, 14, 21, 28, 35)
vec_t_policy_i = 21
vec_beta_S_i = c(0.76, 1.28, 2.4, 3.6, 4.8)

for(species_i in c("s_aureus", "e_coli")){
  for(setting_i in c("rehab", "geriatric", "paediatric")){
    for(scenario_i in c("overwhelmed", "intermediate", "organized")){
      for(t_policy_i in vec_t_policy_i){
        for(beta_S_i in vec_beta_S_i){
          
          print(paste(species_i, setting_i, scenario_i, t_policy_i, beta_S_i))
          
          qounter = qounter + 1
          
          vec_setting_scenario_species_i = paste0(setting_i, "_", scenario_i, "_", species_i)
          
          filename_indicators_nocovid_i = paste0(vec_setting_scenario_species_i, "_indicators_nocovid_t_policy_", t_policy_i, "_beta_S_", beta_S_i, ".Rdata")
          filename_indicators_covid_i = paste0(vec_setting_scenario_species_i, "_indicators_covid_t_policy_", t_policy_i, "_beta_S_", beta_S_i,  ".Rdata")
          filename_indicators_deltas_i = paste0(vec_setting_scenario_species_i, "_indicators_deltas_t_policy_", t_policy_i, "_beta_S_", beta_S_i,  ".Rdata")
          filename_metadata_nocovid_i = paste0(vec_setting_scenario_species_i, "_metadata_nocovid_t_policy_", t_policy_i, "_beta_S_", beta_S_i,  ".Rdata")
          filename_metadata_covid_i = paste0(vec_setting_scenario_species_i, "_metadata_covid_t_policy_", t_policy_i, "_beta_S_", beta_S_i,  ".Rdata")
          filename_metadata_deltas_i = paste0(vec_setting_scenario_species_i, "_metadata_deltas_t_policy_", t_policy_i, "_beta_S_", beta_S_i,  ".Rdata")
          
          data_indicators_nocovid_i = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species_i, "/", filename_indicators_nocovid_i))%>%
            mutate(setting = setting_i,
                   species = species_i,
                   scenario = scenario_i)
          data_indicators_covid_i = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species_i, "/", filename_indicators_covid_i))%>%
            mutate(setting = setting_i,
                   species = species_i,
                   scenario = scenario_i)
          data_indicators_deltas_i = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species_i, "/", filename_indicators_deltas_i))%>%
            mutate(setting = setting_i,
                   species = species_i,
                   scenario = scenario_i)
          data_metadata_nocovid_i = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species_i, "/", filename_metadata_nocovid_i))%>%
            mutate(setting = setting_i,
                   species = species_i,
                   scenario = scenario_i)
          data_metadata_covid_i = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species_i, "/", filename_metadata_covid_i))%>%
            mutate(setting = setting_i,
                   species = species_i,
                   scenario = scenario_i)
          data_metadata_deltas_i = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species_i, "/", filename_metadata_deltas_i))%>%
            mutate(setting = setting_i,
                   species = species_i,
                   scenario = scenario_i)
          
          list_casestudy_indicators_nocovid[[qounter]] = data_indicators_nocovid_i
          list_casestudy_indicators_covid[[qounter]] = data_indicators_covid_i
          list_casestudy_indicators_deltas[[qounter]] = data_indicators_deltas_i
          list_casestudy_metadata_nocovid[[qounter]] = data_metadata_nocovid_i
          list_casestudy_metadata_covid[[qounter]] = data_metadata_covid_i
          list_casestudy_metadata_deltas[[qounter]] = data_metadata_deltas_i
        }
      }
    }
  }
}

# rbind lists into df and update factor labels
df_casestudy_indicators_nocovid = do.call(rbind, list_casestudy_indicators_nocovid)%>%
  mutate(setting = factor(setting,
                          levels = c("geriatric", "paediatric", "rehab"),
                          labels = vec_wards_clean),
         species = factor(species,
                          levels = c("s_aureus", "e_coli"),
                          labels = vec_species_clean))
df_casestudy_indicators_covid = do.call(rbind, list_casestudy_indicators_covid)%>%
  mutate(setting = factor(setting,
                          levels = c("geriatric", "paediatric", "rehab"),
                          labels = vec_wards_clean),
         species = factor(species,
                          levels = c("s_aureus", "e_coli"),
                          labels = vec_species_clean))
df_casestudy_indicators_deltas = do.call(rbind, list_casestudy_indicators_deltas)%>%
  mutate(setting = factor(setting,
                          levels = c("geriatric", "paediatric", "rehab"),
                          labels = vec_wards_clean),
         species = factor(species,
                          levels = c("s_aureus", "e_coli"),
                          labels = vec_species_clean))
df_casestudy_metadata_nocovid = do.call(rbind, list_casestudy_metadata_nocovid)%>%
  mutate(setting = factor(setting,
                          levels = c("geriatric", "paediatric", "rehab"),
                          labels = vec_wards_clean),
         species = factor(species,
                          levels = c("s_aureus", "e_coli"),
                          labels = vec_species_clean))
df_casestudy_metadata_covid = do.call(rbind, list_casestudy_metadata_covid)%>%
  mutate(setting = factor(setting,
                          levels = c("geriatric", "paediatric", "rehab"),
                          labels = vec_wards_clean),
         species = factor(species,
                          levels = c("s_aureus", "e_coli"),
                          labels = vec_species_clean))
df_casestudy_metadata_deltas = do.call(rbind, list_casestudy_metadata_deltas)%>%
  mutate(setting = factor(setting,
                          levels = c("geriatric", "paediatric", "rehab"),
                          labels = vec_wards_clean),
         species = factor(species,
                          levels = c("s_aureus", "e_coli"),
                          labels = vec_species_clean))



#########################
### NUMERICAL RESULTS ###
#########################

### RAW INDICATORS
df_casestudy_indicators_covid%>%
  group_by(species, setting, scenario, t_policy, beta_S)%>%
  summarise(lower = quantile(Cr_inc_cumul_all, 0.25),
            median = median(Cr_inc_cumul_all),
            upper = quantile(Cr_inc_cumul_all, 0.75))


### DELTAS: INDICATORS
# CHANGE IN INCIDENCE
df_casestudy_indicators_deltas%>%
  group_by(species, setting, scenario, t_policy, beta_S)%>%
  summarise(mean = (mean(Cr_inc_cumul_all)-1)*100,
            lower = (quantile(Cr_inc_cumul_all, 0.025)-1)*100,
            lower_iqr = (quantile(Cr_inc_cumul_all, 0.25)-1)*100,
            median = (median(Cr_inc_cumul_all)-1)*100,
            upper_iqr = (quantile(Cr_inc_cumul_all, 0.75)-1)*100,
            upper = (quantile(Cr_inc_cumul_all, 0.975)-1)*100)

# CHANGE IN PD COLONIZED
df_casestudy_indicators_deltas%>%
  group_by(species, setting, scenario, t_policy, beta_S)%>%
  summarise(mean = (mean(Cr_pd)-1)*100,
            lower = (quantile(Cr_pd, 0.025)-1)*100,
            lower_iqr = (quantile(Cr_pd, 0.25)-1)*100,
            median = (median(Cr_pd)-1)*100,
            upper_iqr = (quantile(Cr_pd, 0.75)-1)*100,
            upper = (quantile(Cr_pd, 0.975)-1)*100)

# CHANGE IN R RATE
df_casestudy_indicators_deltas%>%
  group_by(species, setting, scenario, t_policy, beta_S)%>%
  summarise(mean = (mean(R_rate_pd)-1)*100,
            lower = (quantile(R_rate_pd, 0.025)-1)*100,
            lower_iqr = (quantile(R_rate_pd, 0.25)-1)*100,
            median = (median(R_rate_pd)-1)*100,
            upper_iqr = (quantile(R_rate_pd, 0.75)-1)*100,
            upper = (quantile(R_rate_pd, 0.975)-1)*100)




##############################################################
### RAW A: ARB outbreak characteristics (without COVID-19) ###
##############################################################

# Average colonization prevalence
p_cs_nocovid_Cprev = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Cr_prev_max)%>%
  ggplot(aes(y = Cr_prev_max*100,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('patient colonization prevalence (%)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# Average HCW carriage
p_cs_nocovid_Tprev = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Tr_prev_max)%>%
  ggplot(aes(y = Tr_prev_max*100,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('HCW carriage prevalence (%)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Tr_prev_max*100, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# cumulative colonization incidence
p_cs_nocovid_Binc = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Cr_incRate_all)%>%
  ggplot(aes(y = Cr_incRate_all,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('colonization incidence rate')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_inc_cumul_all, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
        )

# colonization incidencerate (by route)
p_cs_nocovid_Binc_pa_pa = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Cr_incRate_pa_pa)%>%
  ggplot(aes(y = Cr_incRate_pa_pa,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('incidence rate (patient-to-patient)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_inc_cumul_pa_pa, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p_cs_nocovid_Binc_pe_pa = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Cr_incRate_pe_pa)%>%
  ggplot(aes(y = Cr_incRate_pe_pa,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('incidence rate (HCW-to-patient)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_inc_cumul_pe_pa, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p_cs_nocovid_Binc_endog = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Cr_incRate_endog)%>%
  ggplot(aes(y = Cr_incRate_endog,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('incidence rate (endogenous)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_inc_cumul_endog, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# carriage incidence rate
p_cs_nocovid_Binc_hcw = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Tr_incRate_all)%>%
  ggplot(aes(y = Tr_incRate_all,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('HCW carriage incidence rate')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Tr_inc_cumul_all, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# cumultaive patient days colonized
p_cs_nocovid_Bpd = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, Cr_pd)%>%
  ggplot(aes(y = Cr_pd,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('cumulative patient-days colonized')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_pd, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# cumulative resistance rate
p_cs_nocovid_Rrate = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  dplyr::select(n_sim, setting, species, R_rate_pd)%>%
  ggplot(aes(y = R_rate_pd*100,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('average resistance rate (%)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$R_rate_pd*100, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

p_cs_nocovid_summary_patients = ggarrange(p_cs_nocovid_Cprev+ggtitle(expression(paste(bold("a.  "), "patient prevalence"))),
                                 p_cs_nocovid_Binc+ggtitle(expression(paste(bold("b.  "), "patient incidence"))), 
                                 p_cs_nocovid_Bpd+ggtitle(expression(paste(bold("c.  "), "patient-days colonized"))), 
                                 p_cs_nocovid_Rrate+ggtitle(expression(paste(bold("d.  "), "resistance rate"))),
                                 ncol = 4, nrow = 1,
                                 common.legend = T, legend = 'bottom')

p_cs_nocovid_summary = ggarrange(p_cs_nocovid_Cprev+ggtitle(expression(paste(bold("a.  "), "patient prevalence"))),
          p_cs_nocovid_Binc+ggtitle(expression(paste(bold("b.  "), "patient incidence"))), 
          p_cs_nocovid_Bpd+ggtitle(expression(paste(bold("c.  "), "patient-days colonized"))), 
          p_cs_nocovid_Rrate+ggtitle(expression(paste(bold("d.  "), "resistance rate"))), 
          p_cs_nocovid_Tprev+ggtitle(expression(paste(bold("e.  "), "HCW carriage"))),
          p_cs_nocovid_Binc_hcw+ggtitle(expression(paste(bold("f.  "), "HCW incidence"))),
          ncol = 3, nrow = 2,
          common.legend = T, legend = 'bottom')

# colonization incidence --> relative contribution by route
p_cs_nocovid_Binc_routes_relative = df_casestudy_indicators_nocovid%>%
  filter(t_policy == 21, beta_S == 1.28, scenario == "organized")%>%
  mutate(`patient to patient` = Cr_inc_cumul_pa_pa/(Cr_inc_cumul_pa_pa+ Cr_inc_cumul_pe_pa + Cr_inc_cumul_endog),
         `HCW to patient` = Cr_inc_cumul_pe_pa/(Cr_inc_cumul_pa_pa+ Cr_inc_cumul_pe_pa + Cr_inc_cumul_endog),
         `endogenous` = Cr_inc_cumul_endog/(Cr_inc_cumul_pa_pa+ Cr_inc_cumul_pe_pa + Cr_inc_cumul_endog))%>%
  dplyr::select(n_sim, setting, species, `patient to patient`, `HCW to patient`, `endogenous`)%>%
  pivot_longer(-c(n_sim, setting, species), names_to = "names", values_to = "values")%>%
  ggplot(aes(y = values,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('relative incidence (share of patient colonization acquisition events)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_inc_cumul_endog, c(0.025, 0.975)))+
  facet_grid(cols = vars(species), rows = vars(names), scales = 'free_x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'bottom'
  )+guides(fill=guide_legend(nrow=1,byrow=TRUE))

ggsave(p_cs_nocovid_summary_patients, filename = "plots/cs_nocovid_arb_summary_patients.pdf", width = 30, height = 10, units = "cm")
ggsave(p_cs_nocovid_summary_patients, filename = "plots/cs_nocovid_arb_summary_patients.png", width = 30, height = 10, units = "cm")

ggsave(p_cs_nocovid_summary, filename = "plots/cs_nocovid_arb_summary.pdf", width = 30, height = 20, units = "cm")
ggsave(p_cs_nocovid_summary, filename = "plots/cs_nocovid_arb_summary.png", width = 30, height = 20, units = "cm")

ggsave(p_cs_nocovid_Binc_routes_relative, filename = "plots/cs_nocovid_arb_routes_relative.pdf", width = 15, height = 15, units = "cm")
ggsave(p_cs_nocovid_Binc_routes_relative, filename = "plots/cs_nocovid_arb_routes_relative.png", width = 15, height = 15, units = "cm")


#####################################################################
### RAW B: Healthcare facility characteristics (without COVID-19) ### "metadata"
#####################################################################


# Antibiotics
p_cs_nocovid_abx = df_casestudy_metadata_nocovid%>%
  filter(species == "MRSA", t_policy == 21, beta_S == 1.28, scenario == "intermediate")%>%
  dplyr::select(n_sim, setting, species, abx_daily)%>%
  ggplot(aes(y = abx_daily,fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('average daily # of patients\nexposed to antibiotics')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# kappa
p_cs_nocovid_contacts = df_casestudy_metadata_nocovid%>%
  filter(species == "MRSA", t_policy == 21, beta_S == 1.28, scenario == "intermediate")%>%
  dplyr::select(n_sim, setting, species, kappa_patients_daily, kappa_staff_daily)%>%
  pivot_longer(-c(n_sim, setting, species), names_to = "indicator", values_to = "value")%>%
  mutate(category = factor(indicator,
                           levels = c("kappa_patients_daily", "kappa_staff_daily"),
                           labels = c("patients", "HCWs")))%>%
  ggplot(aes(y = value, fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('average daily # of contacts')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )+
  facet_grid(cols = vars(category))


# hand hygiene
p_cs_nocovid_hh = df_casestudy_metadata_nocovid%>%
  filter(species == "MRSA", t_policy == 21, beta_S == 1.28, scenario == "intermediate")%>%
  dplyr::select(n_sim, setting, species, hh_daily)%>%
  ggplot(aes(y = 1/hh_daily, fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('average delay between compliant\nhand-washing events (hours)')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# staffing ratio
p_cs_nocovid_staffing_ratio = df_casestudy_metadata_nocovid%>%
  filter(species == "MRSA", t_policy == 21, beta_S == 1.28, scenario == "intermediate")%>%
  dplyr::select(n_sim, setting, species, staffing_ratio_daily)%>%
  ggplot(aes(y = staffing_ratio_daily, fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab("average HCW:patient ratio")+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# admissions
p_cs_nocovid_adm = df_casestudy_metadata_nocovid%>%
  filter(species == "MRSA", t_policy == 21, beta_S == 1.28, scenario == "intermediate")%>%
  dplyr::select(n_sim, setting, species, adm_daily)%>%
  ggplot(aes(y = adm_daily, fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab("average daily patient admissions")+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )

# number of individuals
p_cs_nocovid_totals = df_casestudy_metadata_nocovid%>%
  filter(species == "MRSA", t_policy == 21, beta_S == 1.28, scenario == "intermediate")%>%
  dplyr::select(n_sim, setting, species, patientdays_daily, staffdays_daily)%>%
  pivot_longer(-c(n_sim, setting, species), names_to = "indicator", values_to = "value")%>%
  mutate(category = factor(indicator,
                           levels = c("patientdays_daily", "staffdays_daily"),
                           labels = c("patients", "HCWs")))%>%
  ggplot(aes(y = value, fill = setting, x = species))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, width = 1)+
  theme_bw()+
  scale_fill_manual(name = "healthcare facility", values = col_settings)+
  ylab('average daily # of individuals')+xlab('')+
  #scale_y_continuous(limits = quantile(df_casestudy_indicators_nocovid$Cr_prev_max*100, c(0.025, 0.975)))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )+
  facet_grid(cols = vars(category))

p_cs_nocovid_metadata_summary = ggarrange(p_cs_nocovid_totals+
            ggtitle(expression(paste(bold("a.  "), "totals"))),
          p_cs_nocovid_staffing_ratio+
            ggtitle(expression(paste(bold("b.  "), "staffing"))),
          p_cs_nocovid_hh+
            ggtitle(expression(paste(bold("c.  "), "hand hygiene"))),
          p_cs_nocovid_contacts+
            ggtitle(expression(paste(bold("d.  "), "contact"))),
          p_cs_nocovid_abx+
            ggtitle(expression(paste(bold("e.  "), "antibiotics"))),
          p_cs_nocovid_adm+
            ggtitle(expression(paste(bold("f.  "), "admissions"))),
          nrow = 2, ncol = 3, common.legend = T, legend = 'bottom')


ggsave(p_cs_nocovid_metadata_summary, filename = "plots/cs_nocovid_metadata_summary.pdf", width = 20, height = 15, units = "cm")
ggsave(p_cs_nocovid_metadata_summary, filename = "plots/cs_nocovid_metadata_summary.png", width = 20, height = 15, units = "cm")


#################################
### RAW C: COVID-19 incidence ###
#################################


p_cs_covid_incidence = df_casestudy_indicators_covid%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean),
         `patient to patient` = V_inc_cumul_pa_pa,
         `patient to HCW` = V_inc_cumul_pa_pe,
         `HCW to patient` = V_inc_cumul_pe_pa,
         `HCW to HCW` = V_inc_cumul_pe_pe,
         `total` = V_inc_cumul_all)%>%
  filter(t_policy == 21, beta_S == 1.28, species == "MRSA")%>%
  dplyr::select(n_sim, setting, scenario, species, `patient to patient`, `patient to HCW`, `HCW to patient`, `HCW to HCW`, `total`, setting_scenario)%>%
  pivot_longer(-c(n_sim, setting, scenario, species,setting_scenario), names_to = "route", values_to = "V_inc")%>%
  ggplot(aes(x = V_inc, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('cumulative SARS-CoV-2 infection incidence')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(cols = vars(setting), rows = vars(route), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = 'none'
  )

ggsave(p_cs_covid_incidence, filename = "plots/cs_covid_incidence.pdf", width = 20, height = 18, units = "cm")
ggsave(p_cs_covid_incidence, filename = "plots/cs_covid_incidence.png", width = 20, height = 18, units = "cm")




##############################
### DELTAS: AMR INDICATORS ### "Figure 5"
##############################


### MRSA
### patient colonization incidence rate (all routes combined)
p_cs_deltas_s_aureus_Cinc = df_casestudy_indicators_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean))%>%
  filter(t_policy == 21, beta_S == 1.28, species == "MRSA")%>%
  dplyr::select(n_sim, setting, scenario, species, Cr_inc_cumul_all, setting_scenario)%>%
  mutate(delta = (Cr_inc_cumul_all-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change in cumulative colonization incidence')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )

### patient-days colonized
p_cs_deltas_s_aureus_pd = df_casestudy_indicators_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean))%>%
  filter(t_policy == 21, beta_S == 1.28, species == "MRSA")%>%
  dplyr::select(n_sim, setting, scenario, species, Cr_pd, setting_scenario)%>%
  mutate(delta = (Cr_pd-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change in cumulative patient-days colonized')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )

### resistance rate
p_cs_deltas_s_aureus_Rrate = df_casestudy_indicators_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean))%>%
  filter(t_policy == 21, beta_S == 1.28, species == "MRSA")%>%
  dplyr::select(n_sim, setting, scenario, species, R_rate_pd, setting_scenario)%>%
  mutate(delta = (R_rate_pd-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change in average resistance rate')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )

### E coli
### patient colonization incidence rate (all routes combined)
p_cs_deltas_e_coli_Cinc = df_casestudy_indicators_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean))%>%
  filter(t_policy == 21, beta_S == 1.28, species == "ESBL-E. coli")%>%
  dplyr::select(n_sim, setting, scenario, species, Cr_inc_cumul_all, setting_scenario)%>%
  mutate(delta = (Cr_inc_cumul_all-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change in cumulative colonization incidence')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )

### patient-days colonized
p_cs_deltas_e_coli_pd = df_casestudy_indicators_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean))%>%
  filter(t_policy == 21, beta_S == 1.28, species == "ESBL-E. coli")%>%
  dplyr::select(n_sim, setting, scenario, species, Cr_pd, setting_scenario)%>%
  mutate(delta = (Cr_pd-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change in cumulative patient-days colonized')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )

### resistance rate
p_cs_deltas_e_coli_Rrate = df_casestudy_indicators_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean))%>%
  filter(t_policy == 21, beta_S == 1.28, species == "ESBL-E. coli")%>%
  dplyr::select(n_sim, setting, scenario, species, R_rate_pd, setting_scenario)%>%
  mutate(delta = (R_rate_pd-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change in average resistance rate')+ylab('COVID-19 response scenario\n')+
  #scale_y_continuous(limits = quantile((df_casestudy_indicators_deltas$Cr_inc_cumul_all-1)*100, c(0.025, 0.975)))+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )



### plot limits
inc_max = df_casestudy_indicators_deltas%>%filter(t_policy == 21, beta_S == 1.28)%>%summarise(max((Cr_inc_cumul_all-1)*100))%>%as.numeric()
inc_min = df_casestudy_indicators_deltas%>%filter(t_policy == 21, beta_S == 1.28)%>%summarise(min((Cr_inc_cumul_all-1)*100))%>%as.numeric()
pd_max = df_casestudy_indicators_deltas%>%filter(t_policy == 21, beta_S == 1.28)%>%summarise(max((Cr_pd-1)*100))%>%as.numeric()
pd_min = df_casestudy_indicators_deltas%>%filter(t_policy == 21, beta_S == 1.28)%>%summarise(min((Cr_pd-1)*100))%>%as.numeric()
rrate_max = df_casestudy_indicators_deltas%>%filter(t_policy == 21, beta_S == 1.28)%>%summarise(max((R_rate_pd-1)*100))%>%as.numeric()
rrate_min = df_casestudy_indicators_deltas%>%filter(t_policy == 21, beta_S == 1.28)%>%summarise(min((R_rate_pd-1)*100))%>%as.numeric()



### INDIVIDUAL INDICATORS

### INCIDENCE
p_cs_deltas_incidence = plot_grid(p_cs_deltas_s_aureus_Cinc+theme(legend.position = "none")+
                                         scale_x_continuous(limits = c(inc_min,inc_max))+
                                         ggtitle(expression(paste(bold(a.), "  MRSA"))), 
                                       p_cs_deltas_e_coli_Cinc+theme(legend.position = "none")+
                                         scale_x_continuous(limits = c(inc_min,inc_max))+
                                         ggtitle(expression(paste(bold(b.), "  ESBL-E. coli"))), 
                                       ncol = 1, nrow = 2, align = 'v', axis = 'lr')

### PATIENT DAYS
p_cs_deltas_pd = plot_grid(p_cs_deltas_s_aureus_pd+theme(legend.position = "none")+
                                     scale_x_continuous(limits = c(pd_min,pd_max))+
                                     ggtitle(expression(paste(bold(a.), "  MRSA"))),
                                   p_cs_deltas_e_coli_pd+theme(legend.position = "none")+
                                     scale_x_continuous(limits = c(pd_min,pd_max))+
                                     ggtitle(expression(paste(bold(b.), "  ESBL-E. coli"))),
                                   ncol = 1, nrow = 2, align = 'v', axis = 'lr')

### RESISTANCE RATE
p_cs_deltas_resistance = plot_grid(p_cs_deltas_s_aureus_Rrate+theme(legend.position = "none")+
                                     scale_x_continuous(limits = c(rrate_min,rrate_max))+
                                     ggtitle(expression(paste(bold(a.), "  MRSA"))),
                                   p_cs_deltas_e_coli_Rrate+theme(legend.position = "none")+
                                     scale_x_continuous(limits = c(rrate_min,rrate_max))+
                                     ggtitle(expression(paste(bold(b.), "  ESBL-E. coli"))),
                                   ncol = 1, nrow = 2, align = 'v', axis = 'lr')

ggsave(p_cs_deltas_incidence, filename = "plots/cs_deltas_incidence.pdf", width = 15, height = 22, units = "cm")
ggsave(p_cs_deltas_incidence, filename = "plots/cs_deltas_incidence.png", width = 15, height = 22, units = "cm")

ggsave(p_cs_deltas_pd, filename = "plots/cs_deltas_pd.pdf", width = 15, height = 22, units = "cm")
ggsave(p_cs_deltas_pd, filename = "plots/cs_deltas_pd.png", width = 15, height = 22, units = "cm")

ggsave(p_cs_deltas_resistance, filename = "plots/cs_deltas_resistance.pdf", width = 15, height = 22, units = "cm")
ggsave(p_cs_deltas_resistance, filename = "plots/cs_deltas_resistance.png", width = 15, height = 22, units = "cm")


### COMBINED INDICATORS
p_cs_deltas_incidence_resistance = plot_grid(p_cs_deltas_s_aureus_Cinc+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(inc_min,inc_max))+
                                               ggtitle(label = "cumulative incidence",
                                                       subtitle = "a.  MRSA")+
                                               theme(axis.text.y = element_text(size = 11),
                                                     plot.title = element_text(hjust = 0.5, face = "bold"),
                                                     plot.subtitle = element_text(face = "bold")), 
                                             p_cs_deltas_s_aureus_Rrate+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(rrate_min,rrate_max))+
                                               ggtitle(label = "resistance rate",
                                                       subtitle = "b.  MRSA")+
                                               theme(axis.text.y = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     plot.title = element_text(hjust = 0.5, face = "bold"),
                                                     plot.subtitle = element_text(face = "bold")),
                                             p_cs_deltas_e_coli_Cinc+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(inc_min,inc_max))+
                                               ggtitle(label = " ",
                                                       subtitle = "c.  ESBL-E. coli")+
                                               theme(axis.text.y = element_text(size = 11),
                                                     plot.subtitle = element_text(face = "bold")), 
                                             p_cs_deltas_e_coli_Rrate+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(rrate_min,rrate_max))+
                                               ggtitle(label = " ",
                                                       subtitle = "d.  ESBL-E. coli")+
                                               theme(axis.text.y = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     plot.subtitle = element_text(face = "bold")),
                                             ncol = 2, nrow = 2, align = 'h', axis = 'tb', rel_widths = c(1.35,1,1.35,1))

ggsave(p_cs_deltas_incidence_resistance, filename = "plots/cs_deltas_incidence_resistance.pdf", width = 23, height = 23, units = "cm")
ggsave(p_cs_deltas_incidence_resistance, filename = "plots/cs_deltas_incidence_resistance.png", width = 23, height = 23, units = "cm")


p_cs_deltas_pd_resistance = plot_grid(p_cs_deltas_s_aureus_pd+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(inc_min,inc_max))+
                                               ggtitle(expression(paste(bold(a.), "  MRSA (patient-days colonized)"))), 
                                             p_cs_deltas_s_aureus_Rrate+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(rrate_min,rrate_max))+
                                               ggtitle(expression(paste(bold(b.), "  MRSA (average resistance rate)")))+
                                               theme(axis.text.y = element_blank(),
                                                     axis.title.y = element_blank()),
                                             p_cs_deltas_e_coli_pd+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(inc_min,inc_max))+
                                               ggtitle(expression(paste(bold(c.), "  ESBL-E. coli (patient-days colonized)"))), 
                                             p_cs_deltas_e_coli_Rrate+theme(legend.position = "none")+
                                               scale_x_continuous(limits = c(rrate_min,rrate_max))+
                                               ggtitle(expression(paste(bold(d.), "  ESBL-E. coli (average resistance rate)")))+
                                               theme(axis.text.y = element_blank(),
                                                     axis.title.y = element_blank()),
                                             ncol = 2, nrow = 2, align = 'h', axis = 'tb', rel_widths = c(1.2,1,1.2,1))

ggsave(p_cs_deltas_pd_resistance, filename = "plots/cs_deltas_pd_resistance.pdf", width = 23, height = 23, units = "cm")
ggsave(p_cs_deltas_pd_resistance, filename = "plots/cs_deltas_pd_resistance.png", width = 23, height = 23, units = "cm")


########################
### DELTAS: METADATA ### 
########################


### Antibiotic use
### patient colonization incidence rate (all routes combined)
p_cs_deltas_metadata = df_casestudy_metadata_deltas%>%
  mutate(setting = factor(setting, levels = vec_wards_clean),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")),
         setting_scenario = paste0(setting, " (", scenario, ")"),
         setting_scenario = factor(setting_scenario,
                                   levels = vec_wards_scenarios_clean),
         `antibiotic exposure` = abx_daily,
         `patient contacts` = kappa_patients_daily,
         `HCW contacts` = kappa_staff_daily,
         `hand hygiene` = hh_daily,
         `HCW:patient ratio` = staffing_ratio_daily,
         `patient admissions` = adm_daily,
         `patient-days` = patientdays_daily,
         `HCW-days` = staffdays_daily)%>%
  filter(t_policy == 21, beta_S == 1.28, species == "MRSA")%>%
  dplyr::select(n_sim, setting, scenario, species, `antibiotic exposure`, `patient contacts`, `HCW contacts`, 
                `hand hygiene`, `HCW:patient ratio`, `patient admissions`, `patient-days`, `HCW-days`, setting_scenario)%>%
  pivot_longer(-c(n_sim, setting, scenario, species, setting_scenario), names_to = "metadata", values_to = "value")%>%
  mutate(delta = (value-1)*100)%>%
  ggplot(aes(x = delta, y = scenario, fill = setting_scenario))+
  geom_vline(xintercept = 0)+
  #geom_boxplot(notch = F, outlier.alpha = 0.4, position = 'dodge')+
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), scale = "width")+
  theme_bw()+
  scale_fill_manual(values = col_settings_scenario)+
  xlab('% change due to COVID-19 responses')+ylab('COVID-19 response scenario\n')+
  facet_grid(rows = vars(setting), scales = 'free_y')+
  theme(axis.text.y = element_text(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = 'none'
  )+
  facet_grid(rows = vars(metadata), cols = vars(setting))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-10,-1,0,1,10,100))

ggsave(p_cs_deltas_metadata, filename = "plots/cs_deltas_metadata.pdf", width = 23, height = 25, units = "cm")
ggsave(p_cs_deltas_metadata, filename = "plots/cs_deltas_metadata.png", width = 23, height = 25, units = "cm")









#####################################
### IMPACT OF BETA_S AND T_POLICY ###
#####################################

### change in incidence due to COVID-19

df_sensitivity = df_casestudy_indicators_deltas%>%
  mutate(delta_incidence = (Cr_inc_cumul_all-1)*100,
         delta_pd = (Cr_pd-1)*100,
         delta_resistance = (R_rate_pd-1)*100)%>%
  group_by(setting, scenario, species, t_policy, beta_S)%>%
  summarise(delta_incidence_mean = mean(delta_incidence),
            delta_pd_mean = mean(delta_pd),
            delta_resistance_mean = mean(delta_resistance))%>%
  mutate(t_policy = factor(t_policy, levels = c(7,14,21,28,35)),
         beta_S = factor(beta_S, levels = c(0.76, 1.28, 2.4, 3.6, 4.8)),
         scenario = factor(scenario, levels = c("organized", "intermediate", "overwhelmed")))


### geriatric rehab ward
p_sa_rehab_mrsa_organized_incidence = df_sensitivity%>%
  filter(species == "MRSA", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_rehab_mrsa_organized_resistance = df_sensitivity%>%
  filter(species == "MRSA", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_rehab_mrsa_overwhelmed_incidence = df_sensitivity%>%
  filter(species == "MRSA", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_rehab_mrsa_overwhelmed_resistance = df_sensitivity%>%
  filter(species == "MRSA", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')


p_sa_rehab_esbl_organized_incidence = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_rehab_esbl_organized_resistance = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_rehab_esbl_overwhelmed_incidence = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_rehab_esbl_overwhelmed_resistance = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'rehabilitation ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

### general paediatric ward
p_sa_paediatric_mrsa_organized_incidence = df_sensitivity%>%
  filter(species == "MRSA", setting == 'general ward\n(paediatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_paediatric_mrsa_organized_resistance = df_sensitivity%>%
  filter(species == "MRSA", setting == 'general ward\n(paediatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_paediatric_mrsa_overwhelmed_incidence = df_sensitivity%>%
  filter(species == "MRSA", setting == 'general ward\n(paediatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_paediatric_mrsa_overwhelmed_resistance = df_sensitivity%>%
  filter(species == "MRSA", setting == 'general ward\n(paediatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')


p_sa_paediatric_esbl_organized_incidence = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'general ward\n(paediatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_paediatric_esbl_organized_resistance = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'general ward\n(paediatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_paediatric_esbl_overwhelmed_incidence = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'general ward\n(paediatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_paediatric_esbl_overwhelmed_resistance = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'general ward\n(paediatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

### short-stay geriatric ward
p_sa_geriatric_mrsa_organized_incidence = df_sensitivity%>%
  filter(species == "MRSA", setting == 'short-stay ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_geriatric_mrsa_organized_resistance = df_sensitivity%>%
  filter(species == "MRSA", setting == 'short-stay ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_geriatric_mrsa_overwhelmed_incidence = df_sensitivity%>%
  filter(species == "MRSA", setting == 'short-stay ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_geriatric_mrsa_overwhelmed_resistance = df_sensitivity%>%
  filter(species == "MRSA", setting == 'short-stay ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')


p_sa_geriatric_esbl_organized_incidence = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'short-stay ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_geriatric_esbl_organized_resistance = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'short-stay ward\n(geriatric)', scenario %in% c("organized"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_geriatric_esbl_overwhelmed_incidence = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'short-stay ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_incidence_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_incidence_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in\nincidence due to\nCOVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')

p_sa_geriatric_esbl_overwhelmed_resistance = df_sensitivity%>%
  filter(species == "ESBL-E. coli", setting == 'short-stay ward\n(geriatric)', scenario %in% c("overwhelmed"))%>%
  mutate(label = paste0(round(delta_resistance_mean,1), "%"))%>%
  ggplot(aes(x = beta_S, y = t_policy, fill = delta_resistance_mean))+
  geom_tile(colour = 'white')+
  theme_bw()+
  scale_fill_gradientn("% change in average resistance\nrate due to COVID-19 responses", colours = cols_sensitivity)+
  xlab(expression(paste('SARS-CoV-2 transmission rate (', beta[V], ')')))+ylab(expression(paste('days to COVID-19 policy implementation (', t[policy], ')')))+
  geom_text(aes(x = beta_S, y = t_policy, label = label), colour = 'black', size = size_text, fontface = 'bold')


p_sa_rehab = plot_grid(p_sa_rehab_mrsa_organized_incidence + ggtitle("rehab, mrsa, organized, incidence")+theme(legend.position = 'none'), 
                       p_sa_rehab_mrsa_organized_resistance  + ggtitle("rehab, mrsa, organized, resistance")+theme(legend.position = 'none'),
                       p_sa_rehab_mrsa_overwhelmed_incidence + ggtitle("rehab, mrsa, overwhelmed, incidence")+theme(legend.position = 'none'), 
                       p_sa_rehab_mrsa_overwhelmed_resistance + ggtitle("rehab, mrsa, overwhelmed, resistance")+theme(legend.position = 'none'),
                       p_sa_rehab_esbl_organized_incidence + ggtitle("rehab, esbl, organized, incidence")+theme(legend.position = 'none'), 
                       p_sa_rehab_esbl_organized_resistance + ggtitle("rehab, esbl, organized, resistance")+theme(legend.position = 'none'), 
                       p_sa_rehab_esbl_overwhelmed_incidence + ggtitle("rehab, esbl, overwhelmed, incidence")+theme(legend.position = 'none'), 
                       p_sa_rehab_esbl_overwhelmed_resistance + ggtitle("rehab, esbl, overwhelmed, resistance")+theme(legend.position = 'none'),
                       nrow = 4, ncol = 2,
                       align = 'hv', axis = 'tblr')

p_sa_paediatric = plot_grid(p_sa_paediatric_mrsa_organized_incidence + ggtitle("paediatric, mrsa, organized, incidence")+theme(legend.position = 'none'), 
                       p_sa_paediatric_mrsa_organized_resistance  + ggtitle("paediatric, mrsa, organized, resistance")+theme(legend.position = 'none'),
                       p_sa_paediatric_mrsa_overwhelmed_incidence + ggtitle("paediatric, mrsa, overwhelmed, incidence")+theme(legend.position = 'none'), 
                       p_sa_paediatric_mrsa_overwhelmed_resistance + ggtitle("paediatric, mrsa, overwhelmed, resistance")+theme(legend.position = 'none'),
                       p_sa_paediatric_esbl_organized_incidence + ggtitle("paediatric, esbl, organized, incidence")+theme(legend.position = 'none'), 
                       p_sa_paediatric_esbl_organized_resistance + ggtitle("paediatric, esbl, organized, resistance")+theme(legend.position = 'none'), 
                       p_sa_paediatric_esbl_overwhelmed_incidence + ggtitle("paediatric, esbl, overwhelmed, incidence")+theme(legend.position = 'none'), 
                       p_sa_paediatric_esbl_overwhelmed_resistance + ggtitle("paediatric, esbl, overwhelmed, resistance")+theme(legend.position = 'none'),
                       nrow = 4, ncol = 2,
                       align = 'hv', axis = 'tblr')

p_sa_geriatric = plot_grid(p_sa_geriatric_mrsa_organized_incidence + ggtitle("geriatric, mrsa, organized, incidence")+theme(legend.position = 'none'), 
                            p_sa_geriatric_mrsa_organized_resistance  + ggtitle("geriatric, mrsa, organized, resistance")+theme(legend.position = 'none'),
                            p_sa_geriatric_mrsa_overwhelmed_incidence + ggtitle("geriatric, mrsa, overwhelmed, incidence")+theme(legend.position = 'none'), 
                            p_sa_geriatric_mrsa_overwhelmed_resistance + ggtitle("geriatric, mrsa, overwhelmed, resistance")+theme(legend.position = 'none'),
                            p_sa_geriatric_esbl_organized_incidence + ggtitle("geriatric, esbl, organized, incidence")+theme(legend.position = 'none'), 
                            p_sa_geriatric_esbl_organized_resistance + ggtitle("geriatric, esbl, organized, resistance")+theme(legend.position = 'none'), 
                            p_sa_geriatric_esbl_overwhelmed_incidence + ggtitle("geriatric, esbl, overwhelmed, incidence")+theme(legend.position = 'none'), 
                            p_sa_geriatric_esbl_overwhelmed_resistance + ggtitle("geriatric, esbl, overwhelmed, resistance")+theme(legend.position = 'none'),
                            nrow = 4, ncol = 2,
                            align = 'hv', axis = 'tblr')


ggsave(p_sa_rehab, filename = "plots/cs_sa_rehab.pdf", width = 23, height = 30, units = "cm")
ggsave(p_sa_paediatric, filename = "plots/cs_sa_paediatric.pdf", width = 23, height = 30, units = "cm")
ggsave(p_sa_geriatric, filename = "plots/cs_sa_geriatric.pdf", width = 23, height = 30, units = "cm")

### 
p_sa_mrsa_rehab = plot_grid(p_sa_rehab_mrsa_organized_resistance+
                              theme(legend.position = 'none')+
                              ggtitle(expression(paste(bold(a.), "  organized COVID-19 response"))), 
                            p_sa_rehab_mrsa_overwhelmed_resistance+
                              theme(legend.position = 'none')+
                              ggtitle(expression(paste(bold(b.), "  overwhelmed COVID-19 response"))))


ggsave(p_sa_mrsa_rehab, filename = "plots/cs_sa_mrsa_rehab.pdf", width = 20, height = 10, units = "cm")
ggsave(p_sa_mrsa_rehab, filename = "plots/cs_sa_mrsa_rehab.png", width = 20, height = 10, units = "cm")

### The rest
### MRSA
p_sa_mrsa = plot_grid(p_sa_rehab_mrsa_organized_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(a.), "  organized rehabilitation ward (geriatric)"))), 
                      p_sa_geriatric_mrsa_organized_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(b.), "  organized short-stay ward (geriatric)"))),
                      p_sa_paediatric_mrsa_organized_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(c.), "  organized general ward (paediatric)"))),
                      p_sa_rehab_mrsa_overwhelmed_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(d.), "  overwhelmed rehabilitation ward (geriatric)"))),
                      p_sa_geriatric_mrsa_overwhelmed_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(e.), "  overwhelmed short-stay ward (geriatric)"))),
                      p_sa_paediatric_mrsa_overwhelmed_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(f.), "  overwhelmed general ward (paediatric)"))),
                      ncol = 3,
                      align = 'hv',
                      axis = 'tlbr')

ggsave(p_sa_mrsa, filename = "plots/cs_sa_mrsa.pdf", width = 33, height = 24, units = "cm")
ggsave(p_sa_mrsa, filename = "plots/cs_sa_mrsa.png", width = 33, height = 24, units = "cm")

p_sa_esbl = plot_grid(p_sa_rehab_esbl_organized_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(a.), "  organized rehabilitation ward (geriatric)"))), 
                      p_sa_geriatric_esbl_organized_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(b.), "  organized short-stay ward (geriatric)"))),
                      p_sa_paediatric_esbl_organized_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(c.), "  organized general ward (paediatric)"))),
                      p_sa_rehab_esbl_overwhelmed_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(d.), "  overwhelmed rehabilitation ward (geriatric)"))),
                      p_sa_geriatric_esbl_overwhelmed_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(e.), "  overwhelmed short-stay ward (geriatric)"))),
                      p_sa_paediatric_esbl_overwhelmed_resistance+
                        theme(legend.position = 'none')+
                        ggtitle(expression(paste(bold(f.), "  overwhelmed general ward (paediatric)"))),
                      ncol = 3,
                      align = 'hv',
                      axis = 'tlbr')

ggsave(p_sa_esbl, filename = "plots/cs_sa_esbl.pdf", width = 33, height = 24, units = "cm")
ggsave(p_sa_esbl, filename = "plots/cs_sa_esbl.png", width = 33, height = 24, units = "cm")



