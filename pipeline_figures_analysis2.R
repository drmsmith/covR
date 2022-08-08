library(deSolve)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(cowplot)
library(RColorBrewer)
library(ggallin)
library(ggExtra)
library(wesanderson)

source("ODEs.R")
source("functions.R")
source("pars_states.R")

###########################
### DATA FOR ANALYSIS 2 ###
###########################
data_indicators_raw = loadRData("data/outputs_analysis2/indicators_summarized_raw.Rdata")%>%
  mutate(par = tau)%>%
  left_join(df_pandI_withNone)

data_deltas_raw = loadRData("data/outputs_analysis2/indicators_summarized_deltas.Rdata")%>%
  mutate(par = tau)%>%
  left_join(df_pandI_withNone)%>%
  mutate(impact = factor(impact,
                         levels = vec_impact_clean_order))

data_deltas_means = data_deltas_raw%>%
  group_by(tau)%>%
  summarise(across(where(is.numeric), mean))%>%
  mutate(par = tau)%>%
  left_join(df_pandI_withNone)

#########################
### NUMERICAL RESULTS ###
#########################
vec_quantiles_iqr = c(0.25,0.5,0.75)
vec_quantiles_95 = c(0.025,0.5,0.975)

f_outs = function(data_in){return(c(mean(data_in),
                                    quantile(data_in,c(0.025,0.25,0.5,0.75,0.975))))}
### RAW INDICATORS

data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(Cs_inc_cumul_all, 
                Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog, Cr_inc_cumul_all, 
                Ts_inc_cumul_all, 
                Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe, Tr_inc_cumul_all,
                Cs_prev_max, Cr_prev_max, 
                Cs_prev_min, Cr_prev_min, 
                Ts_prev_max, Tr_prev_max, 
                Ts_prev_min, Tr_prev_min,
                R_rate_pd,
                V_prev_max_all_I,
                V_lag_peak_inc_all,
                V_inc_cumul_all,
                meta_abx_total, meta_kappa_patients_total, meta_kappa_staff_total, meta_hh_total, 
                meta_adm_total, meta_adm_R_total,meta_patientdays_total, meta_staffdays_total, meta_staffing_ratio_total,
                meta_foi_S_total, meta_foi_Cr_total)%>%
  mutate(Cr_prop_pa_pa = Cr_inc_cumul_pa_pa/Cr_inc_cumul_all,
         Cr_prop_pe_pa = Cr_inc_cumul_pe_pa/Cr_inc_cumul_all,
         Cr_prop_endog = Cr_inc_cumul_endog/Cr_inc_cumul_all,
         Tr_prop_pe_pe = Tr_inc_cumul_pe_pe/Tr_inc_cumul_all,
         Tr_prop_pa_pe = Tr_inc_cumul_pa_pe/Tr_inc_cumul_all,)%>%
  summarise(Cs_inc_cumul_all = f_outs(Cs_inc_cumul_all)/180,
            Cr_inc_cumul_pa_pa = f_outs(Cr_inc_cumul_pa_pa)/180,
            Cr_inc_cumul_pe_pa = f_outs(Cr_inc_cumul_pe_pa)/180,
            Cr_inc_cumul_endog = f_outs(Cr_inc_cumul_endog)/180,
            Cr_prop_pa_pa_mean = mean(Cr_prop_pa_pa),
            Cr_prop_pe_pa_mean = mean(Cr_prop_pe_pa),
            Cr_prop_endog_mean = mean(Cr_prop_endog),
            Cr_prop_pa_pa = f_outs(Cr_prop_pa_pa),
            Cr_prop_pe_pa = f_outs(Cr_prop_pe_pa),
            Cr_prop_endog = f_outs(Cr_prop_endog),
            Ts_inc_cumul_all = f_outs(Ts_inc_cumul_all)/180,
            Tr_inc_cumul_all = f_outs(Tr_inc_cumul_all)/180,
            Tr_prop_pe_pe_mean = mean(Tr_prop_pe_pe),
            Tr_prop_pa_pe_mean = mean(Tr_prop_pa_pe),
            Tr_prop_pe_pe = f_outs(Tr_prop_pe_pe),
            Tr_prop_pa_pe = f_outs(Tr_prop_pa_pe),
            Cs_prev_max = f_outs(Cs_prev_max),
            Cs_prev_min = f_outs(Cs_prev_min),
            Cr_prev_max = f_outs(Cr_prev_max),
            Cr_prev_min = f_outs(Cr_prev_min),
            Ts_prev_max = f_outs(Ts_prev_max),
            Ts_prev_min = f_outs(Ts_prev_min),
            Tr_prev_max = f_outs(Tr_prev_max),
            Tr_prev_min = f_outs(Tr_prev_min),
            R_rate_pd = f_outs(R_rate_pd),
            V_prev_max_all_I = f_outs(V_prev_max_all_I),
            V_lag_peak_inc_all = f_outs(V_lag_peak_inc_all),
            V_inc_cumul_all = f_outs(V_inc_cumul_all),
            meta_abx_total = f_outs(meta_abx_total)/180,
            meta_kappa_patients_total = f_outs(meta_kappa_patients_total)/180,
            meta_kappa_staff_total = f_outs(meta_kappa_staff_total)/180, 
            meta_hh_total = f_outs(meta_hh_total)/180, 
            meta_adm_total = f_outs(meta_adm_total)/180,
            meta_adm_R_total = f_outs(meta_adm_R_total)/180,
            meta_patientdays_total = f_outs(meta_patientdays_total)/180,
            meta_staffdays_total = f_outs(meta_staffdays_total)/180,
            meta_staffing_ratio_total = f_outs(meta_staffing_ratio_total)/180,
            meta_foi_S_total = f_outs(meta_foi_S_total)/180,
            meta_foi_Cr_total = f_outs(meta_foi_Cr_total)/180)

### DELTAS
# across all tau
df_deltas_csv = data_deltas_raw%>%
  dplyr::select(tau,
                Cs_inc_cumul_all, 
                Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog, Cr_inc_cumul_all, 
                Ts_inc_cumul_all, 
                Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe, Tr_inc_cumul_all,
                Cs_prev_max, Cr_prev_max, 
                Cs_prev_min, Cr_prev_min, 
                Ts_prev_max, Tr_prev_max, 
                Ts_prev_min, Tr_prev_min,
                Cr_pd, R_rate_pd,
                V_prev_max_all_I,
                V_lag_peak_inc_all,
                V_inc_cumul_all,
                meta_abx_total, meta_kappa_patients_total, meta_kappa_staff_total, meta_hh_total, 
                meta_adm_total, meta_adm_R_total,meta_patientdays_total, meta_staffdays_total, meta_staffing_ratio_total,
                meta_foi_S_total, meta_foi_Cr_total)%>%
  mutate(Cs_inc_cumul_all = (Cs_inc_cumul_all-1)*100, 
         Cr_inc_cumul_pa_pa = (Cr_inc_cumul_pa_pa-1)*100,
         Cr_inc_cumul_pe_pa = (Cr_inc_cumul_pe_pa-1)*100,
         Cr_inc_cumul_endog = (Cr_inc_cumul_endog-1)*100,
         Cr_inc_cumul_all = (Cr_inc_cumul_all-1)*100, 
         Ts_inc_cumul_all = (Ts_inc_cumul_all-1)*100, 
         Tr_inc_cumul_pe_pe = (Tr_inc_cumul_pe_pe-1)*100,
         Tr_inc_cumul_pa_pe = (Tr_inc_cumul_pa_pe-1)*100,
         Tr_inc_cumul_all = (Tr_inc_cumul_all-1)*100,
         Cs_prev_max = (Cs_prev_max-1)*100,
         Cr_prev_max = (Cr_prev_max-1)*100,
         Cs_prev_min = (Cs_prev_min-1)*100,
         Cr_prev_min = (Cr_prev_min-1)*100, 
         Ts_prev_max = (Ts_prev_max-1)*100, 
         Tr_prev_max = (Tr_prev_max-1)*100, 
         Ts_prev_min = (Ts_prev_min-1)*100, 
         Tr_prev_min = (Tr_prev_min-1)*100,
         Cr_pd = (Cr_pd-1)*100,
         R_rate_pd = (R_rate_pd-1)*100,
         V_prev_max_all_I = (V_prev_max_all_I-1)*100,
         V_lag_peak_inc_all = (V_lag_peak_inc_all-1)*100,
         V_inc_cumul_all = (V_inc_cumul_all-1)*100,
         meta_abx_total = (meta_abx_total-1)*100,
         meta_kappa_patients_total = (meta_kappa_patients_total-1)*100, 
         meta_kappa_staff_total = (meta_kappa_staff_total-1)*100,
         meta_hh_total = (meta_hh_total-1)*100, 
         meta_adm_total = (meta_adm_total-1)*100,
         meta_adm_R_total = (meta_adm_R_total-1)*100,
         meta_patientdays_total = (meta_patientdays_total-1)*100,
         meta_staffdays_total = (meta_staffdays_total-1)*100,
         meta_staffing_ratio_total = (meta_staffing_ratio_total-1)*100,
         meta_foi_S_total = (meta_foi_S_total-1)*100,
         meta_foi_Cr_total = (meta_foi_Cr_total-1)*100)%>%
  group_by(tau)%>%
  summarise(Cs_inc_cumul_all = f_outs(Cs_inc_cumul_all),
            Cr_inc_cumul_pa_pa = f_outs(Cr_inc_cumul_pa_pa),
            Cr_inc_cumul_pe_pa = f_outs(Cr_inc_cumul_pe_pa),
            Cr_inc_cumul_endog = f_outs(Cr_inc_cumul_endog),
            Cr_inc_cumul_all = f_outs(Cr_inc_cumul_all),
            Ts_inc_cumul_all = f_outs(Ts_inc_cumul_all),
            Tr_inc_cumul_all = f_outs(Tr_inc_cumul_all),
            Cs_prev_max = f_outs(Cs_prev_max),
            Cs_prev_min = f_outs(Cs_prev_min),
            Cr_prev_max = f_outs(Cr_prev_max),
            Cr_prev_min = f_outs(Cr_prev_min),
            Ts_prev_max = f_outs(Ts_prev_max),
            Ts_prev_min = f_outs(Ts_prev_min),
            Tr_prev_max = f_outs(Tr_prev_max),
            Tr_prev_min = f_outs(Tr_prev_min),
            Cr_pd = f_outs(Cr_pd),
            R_rate_pd = f_outs(R_rate_pd),
            V_prev_max_all_I = f_outs(V_prev_max_all_I),
            V_lag_peak_inc_all = f_outs(V_lag_peak_inc_all),
            V_inc_cumul_all = f_outs(V_inc_cumul_all),
            meta_abx_total = f_outs(meta_abx_total),
            meta_kappa_patients_total = f_outs(meta_kappa_patients_total),
            meta_kappa_staff_total = f_outs(meta_kappa_staff_total), 
            meta_hh_total = f_outs(meta_hh_total), 
            meta_adm_total = f_outs(meta_adm_total),
            meta_adm_R_total = f_outs(meta_adm_R_total),
            meta_patientdays_total = f_outs(meta_patientdays_total),
            meta_staffdays_total = f_outs(meta_staffdays_total),
            meta_staffing_ratio_total = f_outs(meta_staffing_ratio_total),
            meta_foi_S_total = f_outs(meta_foi_S_total),
            meta_foi_Cr_total = f_outs(meta_foi_Cr_total))

tau_csv = "combo_caseload"
indicator_csv = "Cr_inc_cumul_all"
df_deltas_csv%>%
  filter(tau == tau_csv)%>%select(indicator_csv)

# across all tau
data_deltas_raw%>%
  filter(tau == "combo_all")%>%
  dplyr::select(Cs_inc_cumul_all, 
                Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog, Cr_inc_cumul_all, 
                Ts_inc_cumul_all, 
                Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe, Tr_inc_cumul_all,
                Cs_prev_max, Cr_prev_max, 
                Cs_prev_min, Cr_prev_min, 
                Ts_prev_max, Tr_prev_max, 
                Ts_prev_min, Tr_prev_min,
                R_rate_pd,
                V_prev_max_all_I,
                V_lag_peak_inc_all,
                V_inc_cumul_all,
                meta_abx_total, meta_kappa_patients_total, meta_kappa_staff_total, meta_hh_total, 
                meta_adm_total, meta_adm_R_total,meta_patientdays_total, meta_staffdays_total, meta_staffing_ratio_total,
                meta_foi_S_total, meta_foi_Cr_total)%>%
  mutate(Cs_inc_cumul_all = (Cs_inc_cumul_all-1)*100, 
         Cr_inc_cumul_pa_pa = (Cr_inc_cumul_pa_pa-1)*100,
         Cr_inc_cumul_pe_pa = (Cr_inc_cumul_pe_pa-1)*100,
         Cr_inc_cumul_endog = (Cr_inc_cumul_endog-1)*100,
         Cr_inc_cumul_all = (Cr_inc_cumul_all-1)*100, 
         Ts_inc_cumul_all = (Ts_inc_cumul_all-1)*100, 
         Tr_inc_cumul_pe_pe = (Tr_inc_cumul_pe_pe-1)*100,
         Tr_inc_cumul_pa_pe = (Tr_inc_cumul_pa_pe-1)*100,
         Tr_inc_cumul_all = (Tr_inc_cumul_all-1)*100,
         Cs_prev_max = (Cs_prev_max-1)*100,
         Cr_prev_max = (Cr_prev_max-1)*100,
         Cs_prev_min = (Cs_prev_min-1)*100,
         Cr_prev_min = (Cr_prev_min-1)*100, 
         Ts_prev_max = (Ts_prev_max-1)*100, 
         Tr_prev_max = (Tr_prev_max-1)*100, 
         Ts_prev_min = (Ts_prev_min-1)*100, 
         Tr_prev_min = (Tr_prev_min-1)*100,
         R_rate_pd = (R_rate_pd-1)*100,
         V_prev_max_all_I = (V_prev_max_all_I-1)*100,
         V_lag_peak_inc_all = (V_lag_peak_inc_all-1)*100,
         V_inc_cumul_all = (V_inc_cumul_all-1)*100,
         meta_abx_total = (meta_abx_total-1)*100,
         meta_kappa_patients_total = (meta_kappa_patients_total-1)*100, 
         meta_kappa_staff_total = (meta_kappa_staff_total-1)*100,
         meta_hh_total = (meta_hh_total-1)*100, 
         meta_adm_total = (meta_adm_total-1)*100,
         meta_adm_R_total = (meta_adm_R_total-1)*100,
         meta_patientdays_total = (meta_patientdays_total-1)*100,
         meta_staffdays_total = (meta_staffdays_total-1)*100,
         meta_staffing_ratio_total = (meta_staffing_ratio_total-1)*100,
         meta_foi_S_total = (meta_foi_S_total-1)*100,
         meta_foi_Cr_total = (meta_foi_Cr_total-1)*100)%>%
  summarise(Cs_inc_cumul_all = quantile(Cs_inc_cumul_all, vec_quantiles_iqr),
            Cr_inc_cumul_pa_pa = quantile(Cr_inc_cumul_pa_pa, vec_quantiles_iqr),
            Cr_inc_cumul_pe_pa = quantile(Cr_inc_cumul_pe_pa, vec_quantiles_iqr),
            Cr_inc_cumul_endog = quantile(Cr_inc_cumul_endog, vec_quantiles_iqr),
            Ts_inc_cumul_all = quantile(Ts_inc_cumul_all, vec_quantiles_iqr),
            Tr_inc_cumul_all = quantile(Tr_inc_cumul_all, vec_quantiles_iqr),
            Cs_prev_max = quantile(Cs_prev_max, vec_quantiles_iqr),
            Cs_prev_min = quantile(Cs_prev_min, vec_quantiles_iqr),
            Cr_prev_max = quantile(Cr_prev_max, vec_quantiles_iqr),
            Cr_prev_min = quantile(Cr_prev_min, vec_quantiles_iqr),
            Ts_prev_max = quantile(Ts_prev_max, vec_quantiles_iqr),
            Ts_prev_min = quantile(Ts_prev_min, vec_quantiles_iqr),
            Tr_prev_max = quantile(Tr_prev_max, vec_quantiles_iqr),
            Tr_prev_min = quantile(Tr_prev_min, vec_quantiles_iqr),
            R_rate_pd = quantile(R_rate_pd, vec_quantiles_iqr),
            V_prev_max_all_I = quantile(V_prev_max_all_I, vec_quantiles_iqr),
            V_lag_peak_inc_all = quantile(V_lag_peak_inc_all, vec_quantiles_iqr),
            V_inc_cumul_all = quantile(V_inc_cumul_all, vec_quantiles_iqr),
            meta_abx_total = quantile(meta_abx_total, vec_quantiles_iqr),
            meta_kappa_patients_total = quantile(meta_kappa_patients_total, vec_quantiles_iqr),
            meta_kappa_staff_total = quantile(meta_kappa_staff_total, vec_quantiles_iqr), 
            meta_hh_total = quantile(meta_hh_total, vec_quantiles_iqr), 
            meta_adm_total = quantile(meta_adm_total, vec_quantiles_iqr),
            meta_adm_R_total = quantile(meta_adm_R_total, vec_quantiles_iqr),
            meta_patientdays_total = quantile(meta_patientdays_total, vec_quantiles_iqr),
            meta_staffdays_total = quantile(meta_staffdays_total, vec_quantiles_iqr),
            meta_staffing_ratio_total = quantile(meta_staffing_ratio_total, vec_quantiles_iqr),
            meta_foi_S_total = quantile(meta_foi_S_total, vec_quantiles_iqr),
            meta_foi_Cr_total = quantile(meta_foi_Cr_total, vec_quantiles_iqr))





data_deltas_raw%>%
  dplyr::select(tau, Cr_inc_cumul_all, Cr_pd, R_rate_pd, V_inc_cumul_all)%>%
  group_by(tau)%>%
  summarise(Cr_inc_cumul_all = (1-c(mean(Cr_inc_cumul_all), quantile(Cr_inc_cumul_all, vec_quantiles_95)))*100,
            Cr_pd = (1-c(mean(Cr_pd), quantile(Cr_pd, vec_quantiles_95)))*100,
            R_rate_pd = (1-c(mean(R_rate_pd), quantile(R_rate_pd, vec_quantiles_95)))*100,
            V_inc_cumul_all = (1-c(mean(V_inc_cumul_all), quantile(V_inc_cumul_all, vec_quantiles_95)))*100)



##############################
### FIGURES FOR ANALYSIS 2 ###
##############################

##################
### PARAMETERS ###
##################



### NB: will compare parameters for analyses 2 and 3, so load parameters from both
df_parsets_general_species = loadRData("data/outputs_analysis2/parsets_analysis2.Rdata")%>%
  as.data.frame()%>%
  mutate(bacteria = "generic MRB")
df_parsets_e_coli = loadRData("data/outputs_analysis3/geriatric_overwhelmed_e_coli/parsets_analysis3_geriatric_overwhelmed_e_coli.Rdata")%>%
  as.data.frame()%>%
  mutate(bacteria = "E. coli")
df_parsets_s_aureus = loadRData("data/outputs_analysis3/geriatric_overwhelmed_s_aureus/parsets_analysis3_geriatric_overwhelmed_s_aureus.Rdata")%>%
  as.data.frame()%>%
  mutate(bacteria = "S. aureus")

df_parsets_general_setting = loadRData("data/outputs_analysis2/parsets_analysis2.Rdata")%>%
  as.data.frame()%>%
  mutate(setting = "generic hospital")
df_parsets_geriatric = loadRData("data/outputs_analysis3/geriatric_overwhelmed_e_coli/parsets_analysis3_geriatric_overwhelmed_e_coli.Rdata")%>%
  as.data.frame()%>%
  mutate(setting = vec_wards_clean[1])
df_parsets_paediatric = loadRData("data/outputs_analysis3/paediatric_overwhelmed_e_coli/parsets_analysis3_paediatric_overwhelmed_e_coli.Rdata")%>%
  as.data.frame()%>%
  mutate(setting = vec_wards_clean[2])
df_parsets_rehab = loadRData("data/outputs_analysis3/rehab_overwhelmed_e_coli/parsets_analysis3_rehab_overwhelmed_e_coli.Rdata")%>%
  as.data.frame()%>%
  mutate(setting = vec_wards_clean[3])

df_parsets_bacteria = rbind(df_parsets_general_species, df_parsets_e_coli, df_parsets_s_aureus)
df_parsets_setting = rbind(df_parsets_general_setting, df_parsets_geriatric, df_parsets_paediatric, df_parsets_rehab)


p_pars_species = df_parsets_bacteria%>%
  dplyr::select(bacteria, f_Cr, pi_Cr, alpha, gamma_Cr, r_R)%>%
  pivot_longer(-bacteria, names_to = "parameter", values_to = "value")%>%
  mutate(bacteria = factor(bacteria, levels = c("generic MRB", "E. coli", "S. aureus")),
         parameter = factor(parameter, 
                            levels = c("f_Cr", "pi_Cr", "gamma_Cr", "alpha", "r_R"),
                            labels = c(expression(f[C^R]), expression(pi[C^R]), expression(gamma[C^R]), expression(paste(alpha,"'")), expression(r[R]))))%>%
  ggplot(aes(x = value, y = bacteria, colour = bacteria))+
  geom_jitter(alpha = 0.6)+
  geom_boxplot(alpha = 1, colour = 'black', fill = NA, outlier.shape = NA)+
  facet_wrap(facets = vars(parameter), scales = "free_x", nrow = 1, labeller = label_parsed)+
  theme_bw()+
  theme(legend.position = 'none',
        strip.text = element_text(size = 12))+
  ylab('')+xlab('parameter value')

p_pars_settings = df_parsets_setting%>%
  dplyr::select(setting, mu, a, kappa_pa_pe, kappa_pa_pa, kappa_pe_pe)%>%
  mutate(mu = 1/mu)%>%
  pivot_longer(-setting, names_to = "parameter", values_to = "value")%>%
  mutate(setting = factor(setting, levels = c("generic hospital", vec_wards_clean)),
         parameter = factor(parameter, 
                            levels = c("mu", "a", "kappa_pa_pe", "kappa_pa_pa", "kappa_pe_pe"),
                            labels = c(expression(mu^-1), expression(A[base]), expression(kappa^"pat->hcw"), expression(kappa^"pat->pat"), expression(kappa^"hcw->hcw"))))%>%
  ggplot(aes(x = value, y = setting, colour = setting))+
  geom_jitter(alpha = 0.6)+
  geom_boxplot(alpha = 1, colour = 'black', fill = NA, outlier.shape = NA)+
  facet_wrap(facets = vars(parameter), scales = "free_x", nrow = 1, labeller = label_parsed)+
  theme_bw()+
  theme(legend.position = 'none',
        strip.text = element_text(size = 12))+
  ylab('')+xlab('parameter value')
  

p_pars_species_settings = plot_grid(p_pars_species, p_pars_settings, labels = c('a.', 'b.'), nrow = 2, align = 'v', axis = 'lr')

ggsave(p_pars_species_settings, file = paste0("plots", "/pars_analysis2_vs_analysis3.png"), width = 25, height = 20, unit = 'cm')
ggsave(p_pars_species_settings, file = paste0("plots", "/pars_analysis2_vs_analysis3.pdf"), width = 25, height = 20, unit = 'cm')


##################################
### RAW INDICATORS: SARS-CoV-2 ###
##################################
p_V_prev_max = data_indicators_raw%>%dplyr::select(n_sim, tau, category, V_prev_max_pa, V_prev_max_pe, V_prev_max_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_prev_max_pa', "patient", ifelse(name == 'V_prev_max_pe', "HCW", "all")))%>%
  ggplot(aes(x = value*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('maximum SARS-CoV-2 infection prevalence (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Max daily incidence (patient/staff)
p_V_inc_daily_max = data_indicators_raw%>%dplyr::select(n_sim, tau, category, V_inc_daily_max_pa_pa, V_inc_daily_max_pa_pe, V_inc_daily_max_pe_pa, V_inc_daily_max_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_inc_daily_max_pa_pa', "patient -> patient", 
                        ifelse(name == 'V_inc_daily_max_pa_pe', "patient -> HCW",
                               ifelse(name == 'V_inc_daily_max_pe_pa', "HCW -> patient", "HCW -> HCW"))))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('maximum daily SARS-CoV-2 infection incidence\n(axis limited to x<=100 for ease of visual inspection)')+
  coord_cartesian(xlim = c(0, 100))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Cumulative incidence (patient/staff)
p_V_inc_cumul = data_indicators_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_pa_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pa, V_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_inc_cumul_pa_pa', "patient -> patient", 
                        ifelse(name == 'V_inc_cumul_pa_pe', "patient -> HCW",
                               ifelse(name == 'V_inc_cumul_pe_pa', "HCW -> patient", "HCW -> HCW"))))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('cumulative SARS-CoV-2 infection incidence\n(axis limited to x<=4,000 for ease of visual inspection)')+
  coord_cartesian(xlim = c(0, 4000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Lag to max prevalence
p_V_lag_peak_inc = data_indicators_raw%>%dplyr::select(n_sim, tau, category, V_lag_peak_inc_pa_pa, V_lag_peak_inc_pa_pe, V_lag_peak_inc_pe_pa, V_lag_peak_inc_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_lag_peak_inc_pa_pa', "patient -> patient", 
                        ifelse(name == 'V_lag_peak_inc_pa_pe', "patient -> HCW",
                               ifelse(name == 'V_lag_peak_inc_pe_pa', "HCW -> patient", "HCW -> HCW"))))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('lag to max SARS-CoV-2 infection incidence (days)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################
### RAW INDICATORS: BACTERIA ###
################################

### Max prevalence (stratified by strain)
p_B_prev_max = data_indicators_raw%>%dplyr::select(n_sim, tau, category, Cs_prev_max, Cr_prev_max, Ts_prev_max, Tr_prev_max)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_prev_max', "S strain (patients)",
                        ifelse(name == 'Cr_prev_max', "R strain (patients)",
                               ifelse(name == 'Ts_prev_max', "S strain (HCW)", "R strain (HCW)"))))%>%
  ggplot(aes(x = value*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('maximum bacterial colonization/carriage prevalence (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Max daily incidence (resistant strain, patient/staff)
p_B_inc_daily_max = data_indicators_raw%>%dplyr::select(n_sim, tau, category, Cr_inc_daily_max_pa_pa, Cr_inc_daily_max_pe_pa, Cr_inc_daily_max_endog, Tr_inc_daily_max_pe_pe, Tr_inc_daily_max_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cr_inc_daily_max_pa_pa', "patient -> patient (colonization)", 
                        ifelse(name == 'Cr_inc_daily_max_pe_pa', "HCW -> patient (colonization)",
                               ifelse(name == 'Cr_inc_daily_max_endog', "endogenous (colonization)", 
                                      ifelse(name == "Tr_inc_daily_max_pe_pe", "patient -> HCW (carriage)", "HCW -> HCW (carriage)")))))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('maximum daily bacterial colonization/carriage incidence\n(axis limited to x<=250 for ease of visual inspection)')+
  coord_cartesian(xlim = c(0, 250))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Cumulative incidence (stratified by strain, colonization only)
p_B_inc_cumul = data_indicators_raw%>%dplyr::select(n_sim, tau, category, Cs_inc_cumul_pa_pa, Cs_inc_cumul_pe_pa, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_inc_cumul_pa_pa', "patient -> patient (S strain)", 
                        ifelse(name == 'Cs_inc_cumul_pe_pa', "HCW -> patient (S strain)",
                               ifelse(name == 'Cr_inc_cumul_pa_pa', "patient -> patient (R strain)", 
                                      ifelse(name == "Cr_inc_cumul_pe_pa", "HCW -> patient (R strain)", "endogenous (R strain)")))))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('cumulative bacterial colonization incidence\n(axis limited to x<=1000 for ease of visual inspection)')+
  coord_cartesian(xlim = c(0, 1000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Cumulative incidence resistant strain, relative contribution of each route
p_B_inc_cumul_prop_routes = data_indicators_raw%>%
  mutate(Cr_inc_cumul_prop_pa_pa = Cr_inc_cumul_pa_pa/Cr_inc_cumul_all,
         Cr_inc_cumul_prop_pe_pa = Cr_inc_cumul_pe_pa/Cr_inc_cumul_all,
         Cr_inc_cumul_prop_endog = Cr_inc_cumul_endog/Cr_inc_cumul_all)%>%
  dplyr::select(n_sim, tau, category, Cr_inc_cumul_prop_pa_pa, Cr_inc_cumul_prop_pe_pa, Cr_inc_cumul_prop_endog)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cr_inc_cumul_prop_pa_pa', "patient -> patient", 
                        ifelse(name == 'Cr_inc_cumul_prop_pe_pa', "HCW -> patient", "endogenous")))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('proportion of bacterial colonization incidence')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Cumulative patient-days colonized (stratified by strain, colonization only)
p_B_pd = data_indicators_raw%>%dplyr::select(n_sim, tau, category, Cs_pd, Cr_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_pd', "S strain", "R strain"))%>%
  ggplot(aes(x = value, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('cumulative patient patient-days colonized')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Resistance rate
p_B_R_rate = data_indicators_raw%>%dplyr::select(n_sim, tau, category, R_rate_max, R_rate_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'R_rate_max', "resistance rate (max)", "resistance rate (patient-days colonized)"))%>%
  ggplot(aes(x = value*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('resistance rate (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######################
### BOTH PATHOGENS ###
######################
## side-by-side comparison of cumulative incidence
p_V_B_inc_cumul = data_indicators_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_all, Cr_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_inc_cumul_all', "SARS-CoV-2 infection", "resistant bacteria colonization"))%>%
  ggplot(aes(x = value/1000, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('cumulative incidence (x1,000)')+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


ggsave(p_V_prev_max, file = paste0("plots", "/V_prev_max.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_V_inc_daily_max, file = paste0("plots", "/V_inc_daily_max.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_V_inc_cumul, file = paste0("plots", "/V_inc_cumul.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_V_lag_peak_inc, file = paste0("plots", "/V_lag_peak_inc.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_B_prev_max, file = paste0("plots", "/B_prev_max.png"), width = 25, height = 14, unit = 'cm')
ggsave(p_B_inc_daily_max, file = paste0("plots", "/B_inc_daily_max.png"), width = 30, height = 14, unit = 'cm')
ggsave(p_B_inc_cumul, file = paste0("plots", "/B_inc_cumul.png"), width = 30, height = 14, unit = 'cm')
ggsave(p_B_inc_cumul_prop_routes, file = paste0("plots", "/B_inc_cumul_prop_routes.png"), width = 30, height = 14, unit = 'cm')
ggsave(p_B_pd, file = paste0("plots", "/B_pd.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_B_R_rate, file = paste0("plots", "/B_R_rate.png"), width = 25, height = 14, unit = 'cm')
ggsave(p_V_B_inc_cumul, file = paste0("plots", "/V_B_inc_cumul.png"), width = 25, height = 14, unit = 'cm')


####################################
### DELTA INDICATORS: SARS-CoV-2 ###
####################################

##################
### SARS-CoV-2 ###
##################

### Max prevalence (stratified by patient/staff)
p_V_prev_max = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_prev_max_pa, V_prev_max_pe, V_prev_max_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_prev_max_pa', "patient", ifelse(name == 'V_prev_max_pe', "HCW", "all")))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in maximum SARS-CoV-2 infection prevalence (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)

### Max daily incidence (patient/staff)
p_V_inc_daily_max = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_daily_max_pa_pa, V_inc_daily_max_pa_pe, V_inc_daily_max_pe_pa, V_inc_daily_max_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_inc_daily_max_pa_pa', "patient -> patient", 
                        ifelse(name == 'V_inc_daily_max_pa_pe', "patient -> HCW",
                               ifelse(name == 'V_inc_daily_max_pe_pa', "HCW -> patient", "HCW -> HCW"))))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in maximum daily SARS-CoV-2 infection incidence (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


### Cumulative incidence (patient/staff)
p_V_inc_cumul = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_pa_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pa, V_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_inc_cumul_pa_pa', "patient -> patient", 
                        ifelse(name == 'V_inc_cumul_pa_pe', "patient -> HCW",
                               ifelse(name == 'V_inc_cumul_pe_pa', "HCW -> patient", "HCW -> HCW"))))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative SARS-CoV-2 infection incidence (%)\n(axis limited to x>=-200 for ease of visual inspection)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)+
  coord_cartesian(xlim = c(-100, 200))

### Lag to max prevalence
p_V_lag_peak_inc = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_lag_peak_inc_pa_pa, V_lag_peak_inc_pa_pe, V_lag_peak_inc_pe_pa, V_lag_peak_inc_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_lag_peak_inc_pa_pa', "patient -> patient", 
                        ifelse(name == 'V_lag_peak_inc_pa_pe', "patient -> HCW",
                               ifelse(name == 'V_lag_peak_inc_pe_pa', "HCW -> patient", "HCW -> HCW"))))%>%
  ggplot(aes(x = (1-value)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in delay to max SARS-CoV-2 infection incidence (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


################
### BACTERIA ###
################

### Max prevalence (stratified by strain)
p_B_prev_max = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cs_prev_max, Cr_prev_max, Ts_prev_max, Tr_prev_max)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_prev_max', "S strain (patients)",
                        ifelse(name == 'Cr_prev_max', "R strain (patients)",
                               ifelse(name == 'Ts_prev_max', "S strain (HCW)", "R strain (HCW)"))))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in maximum bacterial colonization/carriage prevalence (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)+
  coord_cartesian(xlim = c(-20, 100))

### Max daily incidence (resistant strain, patient/staff)
### NB: transient carriage labels were off on first batch of simulations, will need to update accordingly (Tr_inc_daily_max_pa_pa to Tr_inc_daily_max_pe_pe AND Tr_inc_daily_max_pa_pe to Tr_inc_daily_max_pe_pa)
p_B_inc_daily_max = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_inc_daily_max_pa_pa, Cr_inc_daily_max_pe_pa, Cr_inc_daily_max_endog, Tr_inc_daily_max_pe_pe, Tr_inc_daily_max_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cr_inc_daily_max_pa_pa', "patient -> patient (colonization)", 
                        ifelse(name == 'Cr_inc_daily_max_pe_pa', "HCW -> patient (colonization)",
                               ifelse(name == 'Cr_inc_daily_max_endog', "endogenous (colonization)", 
                                      ifelse(name == "Tr_inc_daily_max_pa_pe", "patient -> HCW (carriage)", "HCW -> HCW (carriage)")))))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in maximum daily colonization/carriage incidence of resistant bacteria (%)\n(axis limited to x<=150 for ease of visual inspection)')+
  coord_cartesian(xlim = c(-50, 150))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


### Cumulative incidence (stratified by strain, colonization only)
p_B_inc_cumul = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cs_inc_cumul_pa_pa, Cs_inc_cumul_pe_pa, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_inc_cumul_pa_pa', "patient -> patient (S strain)", 
                        ifelse(name == 'Cs_inc_cumul_pe_pa', "HCW -> patient (S strain)",
                               ifelse(name == 'Cr_inc_cumul_pa_pa', "patient -> patient (R strain)", 
                                      ifelse(name == "Cr_inc_cumul_pe_pa", "HCW -> patient (R strain)", "endogenous (R strain)")))))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative bacterial colonization incidence (%)\n(axis limited to -100<=x<=100 for ease of visual inspection)')+
  coord_cartesian(xlim = c(-100, 100))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


### Cumulative patient-days colonized (stratified by strain, colonization only)
p_B_pd = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cs_pd, Cr_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_pd', "S strain", "R strain"))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative patient patient-days colonized (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


### Resistance rate
p_B_R_rate = data_deltas_raw%>%dplyr::select(n_sim, tau, category, R_rate_max, R_rate_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'R_rate_max', "resistance rate (max)", "resistance rate (patient-days colonized)"))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in resistance rate (%)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)+
  coord_cartesian(xlim = c(-20, 100))


### Incidence rate 
data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cs_incRate_pa_pa, Cs_incRate_pe_pa, Cr_incRate_pa_pa, Cr_incRate_pe_pa, Cr_incRate_endog)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'Cs_incRate_pa_pa', "patient -> patient (S strain)", 
                        ifelse(name == 'Cs_incRate_pe_pa', "HCW -> patient (S strain)",
                               ifelse(name == 'Cr_incRate_pa_pa', "patient -> patient (R strain)", 
                                      ifelse(name == "Cr_incRate_pe_pa", "HCW -> patient (R strain)", "endogenous (R strain)")))))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in bacterial colonization incidence rate (%)\n(axis limited to -100<=x<=100 for ease of visual inspection)')+
  coord_cartesian(xlim = c(-100, 100))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)


######################
### BOTH PATHOGENS ###
######################
## side-by-side comparison of cumulative incidence

p_V_B_inc_cumul = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_all, Cr_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_inc_cumul_all', "SARS-CoV-2", "Bacteria"))%>%
  ggplot(aes(x = (value-1)*100, y = tau, fill = category))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  facet_grid(rows = vars(category), cols = vars(group), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative incidence (%)\n(axis limited to x>=-200 for ease of visual inspection)')+
  geom_vline(xintercept = 0, linetype = 2, alpha = alpha_vline)+
  coord_cartesian(xlim = c(-100, 200))


ggsave(p_V_prev_max, file = paste0("plots", "/delta_V_prev_max.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_V_inc_daily_max, file = paste0("plots", "/delta_V_inc_daily_max.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_V_inc_cumul, file = paste0("plots", "/delta_V_inc_cumul.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_V_lag_peak_inc, file = paste0("plots", "/delta_V_lag_peak_inc.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_B_prev_max, file = paste0("plots", "/delta_B_prev_max.png"), width = 25, height = 14, unit = 'cm')
ggsave(p_B_inc_daily_max, file = paste0("plots", "/delta_B_inc_daily_max.png"), width = 30, height = 14, unit = 'cm')
ggsave(p_B_inc_cumul, file = paste0("plots", "/delta_B_inc_cumul.png"), width = 30, height = 14, unit = 'cm')
ggsave(p_B_pd, file = paste0("plots", "/delta_B_pd.png"), width = 20, height = 14, unit = 'cm')
ggsave(p_B_R_rate, file = paste0("plots", "/delta_B_R_rate.png"), width = 25, height = 14, unit = 'cm')
ggsave(p_V_B_inc_cumul, file = paste0("plots", "/delta_V_B_inc_cumul.png"), width = 25, height = 14, unit = 'cm')





###################
### CLEAN PLOTS ###
###################
### In absence of impacts, key indicators: prevalence (both), incidence (both), resistance rate (bacteria), delays (virus)

################
### RAW DATA ###
################

# infection prevalence
p_raw_prevalence_V = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, V_prev_max_pa_I, V_prev_max_pe_I, V_prev_max_all_I)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(group = ifelse(name == 'V_prev_max_pa_I', "patients", ifelse(name == 'V_prev_max_pe_I', "HCWs", "all")))%>%
  mutate(group = factor(group, levels = c("patients", "HCWs", "all")))%>%
  ggplot(aes(y = value*100, x = group, fill = group))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975))+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = col_pat_staff)+
  xlab('')+ylab('maximum prevalence, SARS-CoV-2 infection (%)')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  ylim(0,100)

# infection incidence
p_raw_incidence_V = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, V_inc_cumul_pa_pa, V_inc_cumul_pe_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(route = ifelse(name == 'V_inc_cumul_pa_pa', "patient to patient", 
                        ifelse(name == 'V_inc_cumul_pe_pa', "HCW to patient", 
                               ifelse(name == "V_inc_cumul_pa_pe", "patient to HCW", "HCW to HCW"))))%>%
  mutate(route = factor(route, levels = c("patient to patient", "patient to HCW", "HCW to patient", "HCW to HCW")))%>%
  ggplot(aes(y = value, x = route, fill = route))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975))+
  theme_bw()+
  scale_fill_manual(values = col_incidence)+
  xlab('')+ylab('cumulative incidence, SARS-CoV-2 infection')+
  scale_y_continuous(breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000), labels = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000), trans = "log10")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")

# lags to peak incidence
p_raw_lag_incidence_V = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, V_lag_peak_inc_all)%>%
  ggplot(aes(y = V_lag_peak_inc_all, x = 1))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = "antiquewhite2")+
  theme_bw()+
  xlab('')+ylab('lag to peak SARS-CoV-2 infection incidence (days)')+
  theme(axis.text.x = element_blank(), legend.position = "none",
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# colonization and carriage prevalence
p_raw_prevalence_C = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, Cs_prev_max, Cr_prev_max)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(strain = ifelse(name %in% c("Cs_prev_max"), "S strain", "R strain"))%>%
  mutate(strain = factor(strain, levels = c("S strain", "R strain")))%>%
  ggplot(aes(y = value*100, x = strain, fill = strain))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975))+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = rev(col_strains))+
  xlab('')+ylab('mean prevalence, patient colonization (%)')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")+
  ylim(0,100)

p_raw_prevalence_T = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, Ts_prev_max, Tr_prev_max)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(strain = ifelse(name %in% c("Ts_prev_max"), "S strain", "R strain"))%>%
  mutate(strain = factor(strain, levels = c("S strain", "R strain")))%>%
  ggplot(aes(y = value*100, x = strain, fill = strain))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975))+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = rev(col_strains))+
  xlab('')+ylab('mean prevalence, transient HCW carriage (%)')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")+
  ylim(0,100)

# colonization resistance rate
p_raw_R_rate = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, R_rate_pd)%>%
  ggplot(aes(y = R_rate_pd*100, x = 1))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = 'grey')+
  theme_bw()+
  xlab('')+ylab('mean resistance rate among colonized patients (%)')+
  theme(axis.text.x = element_blank(), legend.position = "none",
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# colonization incidence
p_raw_incidence_C = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(route = ifelse(name == 'Cr_inc_cumul_pa_pa', "patient\nto\npatient", 
                        ifelse(name == 'Cr_inc_cumul_pe_pa', "HCW\nto\npatient", "endogenous")))%>%
  mutate(route = factor(route, levels = c("patient\nto\npatient", "HCW\nto\npatient", "endogenous")))%>%
  ggplot(aes(y = value/180, x = route, fill = route))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975))+
  theme_bw()+
  scale_fill_manual(values = col_incidence[c(1,3,6)])+
  xlab('')+ylab('mean daily incidence, colonization (R strain)')+
  scale_y_continuous(limits = c(0.00001, 2000), breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000), labels = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000), trans = "log10")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")

# carriage incidence
p_raw_incidence_T = data_indicators_raw%>%
  filter(tau == "none")%>%
  dplyr::select(n_sim, tau, category, Tr_inc_cumul_pa_pe, Tr_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(route = ifelse(name == 'Tr_inc_cumul_pa_pe', "patient\nto\nHCW", "HCW\nto\nHCW"))%>%
  mutate(route = factor(route, levels = c("patient\nto\nHCW", "HCW\nto\nHCW")))%>%
  ggplot(aes(y = value/180, x = route, fill = route))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975))+
  theme_bw()+
  scale_fill_manual(values = col_incidence[c(2,4)])+
  xlab('')+ylab('mean daily incidence, transient carriage (R strain)')+
  scale_y_continuous(limits = c(0.00001, 2000), breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000), labels = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000), trans = "log10")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")



top_row = plot_grid(p_raw_prevalence_C+ggtitle(expression(paste(bold("a.")))), 
                       p_raw_prevalence_T+ggtitle(expression(paste(bold("b.")))), 
                       p_raw_R_rate+ggtitle(expression(paste(bold("c.")))), 
                       p_raw_incidence_C+ggtitle(expression(paste(bold("d.")))), 
                       p_raw_incidence_T+ggtitle(expression(paste(bold("e.")))), 
                       ncol = 5, rel_widths = c(1,1,0.7,1.3,1.1), 
                       align = 'h', axis = 'tb')

bottom_row = plot_grid(p_raw_prevalence_V+ggtitle(expression(paste(bold("f.")))),
                       p_raw_incidence_V+ggtitle(expression(paste(bold("g.")))), 
                       p_raw_lag_incidence_V+ggtitle(expression(paste(bold("h.")))), 
                       ncol = 3, rel_widths = c(1.2,2,0.6), 
                       align = 'h', axis = 'tb')

p_summary_no_tau = plot_grid(top_row, bottom_row, nrow = 2)

ggsave(p_summary_no_tau, filename = "plots/summary_no_tau.pdf", width = 30, height = 20, units = "cm")
ggsave(p_summary_no_tau, filename = "plots/summary_no_tau.png", width = 30, height = 20, units = "cm")


################
### METADATA ###
################

data_indicators_raw_metadata = data_indicators_raw%>%
  filter(tau == "none")%>%
  mutate(`antibiotic exposure` = meta_abx_total/180,
         `patient contacts` = meta_kappa_patients_total/180,
         `HCW contacts` = meta_kappa_staff_total/180,
         `hand hygiene` = meta_hh_total/180,
         `HCW:patient ratio` = meta_staffing_ratio_total/180,
         `patient admissions` = meta_adm_R_total/180,
         `patient-days` = meta_patientdays_total/180,
         `HCW-days` = meta_staffdays_total/180)%>%
  dplyr::select(n_sim, tau, category, `antibiotic exposure`, `patient contacts`, `HCW contacts`, `hand hygiene`,`HCW:patient ratio`,
                `patient admissions`, `patient-days`, `HCW-days`)%>%
  pivot_longer(-c(n_sim, tau, category))

p_raw_abx = data_indicators_raw_metadata%>%
  filter(name == 'antibiotic exposure')%>%
  ggplot(aes(y = value, x = name))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = col_abx)+
  theme_bw()+
  #scale_fill_manual(name = "pandemic impact", values = col_pat_staff)+
  xlab('')+ylab('average daily # of patients\nexposed to antibiotics')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  facet_grid(cols = vars(name), scales = 'free')

p_raw_kappa_pat = data_indicators_raw_metadata%>%
  filter(name == 'patient contacts')%>%
  ggplot(aes(y = value, x = name))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = col_pat_staff[1])+
  theme_bw()+
  xlab('')+ylab('average daily # of patient contacts')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  facet_grid(cols = vars(name), scales = 'free')

p_raw_kappa_hcw = data_indicators_raw_metadata%>%
  filter(name == 'HCW contacts')%>%
  ggplot(aes(y = value, x = name))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = col_pat_staff[2])+
  theme_bw()+
  xlab('')+ylab('average daily # of HCW contacts')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  facet_grid(cols = vars(name), scales = 'free')

p_raw_hh = data_indicators_raw_metadata%>%
  filter(name == 'hand hygiene')%>%
  ggplot(aes(y = 1/value, x = name))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = col_hh)+
  theme_bw()+
  xlab('')+ylab('average delay between compliant\nhand-washing events (hours)')+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  facet_grid(cols = vars(name), scales = 'free')

p_raw_staffing = data_indicators_raw_metadata%>%
  filter(name == 'HCW:patient ratio')%>%
  ggplot(aes(y = value, x = name))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = col_staffing)+
  theme_bw()+
  xlab('')+ylab(expression(paste('average HCW:patient ratio')))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  facet_grid(cols = vars(name), scales = 'free')

p_raw_adm = data_indicators_raw_metadata%>%
  filter(name == 'patient admissions')%>%
  ggplot(aes(y = value, x = name))+
  geom_violin(width = 1, draw_quantiles = c(0.025, 0.5, 0.975), fill = col_admission)+
  theme_bw()+
  xlab('')+ylab("average daily # of patients admitted\nalready colonized with MRB")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = 'none')+
  facet_grid(cols = vars(name), scales = 'free')

p_summary_no_tau_metadata = plot_grid(p_raw_abx+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                      p_raw_kappa_pat+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                      p_raw_kappa_hcw+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                      p_raw_hh+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                      p_raw_staffing+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                                      p_raw_adm+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
          ncol = 3, nrow = 2, align = "hv", axis = "tblr")

ggsave(p_summary_no_tau_metadata, filename = "plots/summary_no_tau_metadata.pdf", width = 16, height = 16, units = "cm")
ggsave(p_summary_no_tau_metadata, filename = "plots/summary_no_tau_metadata.png", width = 16, height = 16, units = "cm")

#####################################################
### CHANGE IN INDICATORS FOR EACH PANDEMIC IMPACT ###
#####################################################

### PREVALENCE ###

# SARS-CoV-2
p_delta_V_prev = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_prev_max_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(par %in% vec_pandI_parameters[c(1:10,12)])%>%
  mutate(category = factor(category, 
         levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
         labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  ggplot(aes(y = impact, x = value, fill = category))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(category), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in maximum infection prevalence, SARS-CoV-2  (%)')+
  coord_cartesian(xlim = c(-100, 100))
  
  #scale_x_continuous(limits = quantile((data_deltas_raw$V_prev_max_all-1)*100, c(0.025, 0.975)))

# Bacteria (max prev)
p_delta_B_prev = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_prev_max)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(par %in% vec_pandI_parameters[c(1:10,12)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  ggplot(aes(y = impact, x = value, fill = category))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = F, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(category), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in maximum colonization prevalence, resistant bacteria (%)')+
  coord_cartesian(xlim = c(-100, 100))
  #scale_x_continuous(limits = quantile((data_deltas_raw$Cr_prev_max-1)*100, c(0.025, 0.975)))

# Bacteria (patient-days colonized)
p_delta_B_pd = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(par %in% vec_pandI_parameters[c(1:10,12)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  ggplot(aes(y = impact, x = value, fill = category,))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(category), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative patient-days colonized, resistant bacteria (%)')+
  coord_cartesian(xlim = c(-100, 100))
  #scale_x_continuous(limits = quantile((data_deltas_raw$Cr_pd-1)*100, c(0.025, 0.975)))

### INCIDENCE

# SARS-CoV-2
p_delta_V_inc = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(par %in% vec_pandI_parameters[c(1:10,12)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  ggplot(aes(y = impact, x = value, fill = category))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(category), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative incidence, SARS-CoV-2 (%)')+
  coord_cartesian(xlim = c(-100, 100))
  #scale_x_continuous(limits = quantile((data_deltas_raw$V_inc_cumul_all-1)*100, c(0.025, 0.975)))

# Bacteria
p_delta_B_inc = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(par %in% vec_pandI_parameters[c(1:10,12)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  ggplot(aes(y = impact, x = value, fill = category))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(category), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in cumulative incidence, resistant bacteria (%)')+
  coord_cartesian(xlim = c(-100, 100))

### RESISTANCE RATE

p_delta_R_rate = data_deltas_raw%>%dplyr::select(n_sim, tau, category, R_rate_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(par %in% vec_pandI_parameters[c(1:10,12)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  ggplot(aes(y = impact, x = value, fill = category))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(category), scales = 'free_y', space = "free_y")+
  scale_fill_nejm(name = "pandemic impact")+
  scale_colour_nejm(name = "pandemic impact")+
  ylab('')+xlab('change in average resistance rate (%)')+
  coord_cartesian(xlim = c(-100, 100))
  #scale_x_continuous(limits = quantile((data_deltas_raw$R_rate_pd-1)*100, c(0.025, 0.975)))

### combine

# V incidence vs. B incidence
p_delta_V_inc_B_inc = ggarrange(p_delta_V_inc+ggtitle(expression(paste(bold("a.")))),
                                p_delta_B_inc+ggtitle(expression(paste(bold("b.")))),
                                ncol = 2,
                                legend = "none")
ggsave(p_delta_V_inc_B_inc, file = paste0("plots", "/deltaBinc_deltaVinc.pdf"), width = 29, height = 15, unit = 'cm')
ggsave(p_delta_V_inc_B_inc, file = paste0("plots", "/deltaBinc_deltaVinc.png"), width = 29, height = 15, unit = 'cm')

# V incidence vs. B incidence vs. resistance rate
p_delta_V_inc_B_inc_R_rate = plot_grid(p_delta_V_inc+
                                         ggtitle(expression(paste(bold("a."))))+
                                         theme(legend.position = 'none'),
                                       p_delta_B_inc+ggtitle(expression(paste(bold("b."))))+
                                         theme(legend.position = 'none'),
                                       p_delta_R_rate+ggtitle(expression(paste(bold("c."))))+
                                         theme(legend.position = 'none'),
                                       ncol = 3)
ggsave(p_delta_V_inc_B_inc_R_rate, file = paste0("plots", "/deltaBinc_deltaVinc_Rrate.pdf"), width = 40, height = 15, unit = 'cm')
ggsave(p_delta_V_inc_B_inc_R_rate, file = paste0("plots", "/deltaBinc_deltaVinc_Rrate.png"), width = 40, height = 15, unit = 'cm')


# B patient-days colonized vs B R rate
p_delta_B_pd_B_Rrate = ggarrange(p_delta_B_pd+ggtitle(expression(paste(bold("a.")))),
                                 p_delta_R_rate+ggtitle(expression(paste(bold("b.")))),
                                ncol = 2,
                                legend = "none")

ggsave(p_delta_B_pd_B_Rrate, file = paste0("plots", "/deltaBpd_deltaRrate.pdf"), width = 30, height = 15, unit = 'cm')
ggsave(p_delta_B_pd_B_Rrate, file = paste0("plots", "/deltaBpd_deltaRrate.png"), width = 30, height = 15, unit = 'cm')

# V incidence, B incidence, resistance rate
p_delta_V_inc_routes = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_pa_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pa, V_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[c(1:10,14,18)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')),
         name = factor(name, 
                       levels = c("V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe"),
                       labels = c("patient to patient", "patient to HCW", "HCW to patient", "HCW to HCW")),
         impact = factor(impact, 
                         levels = df_pandI$impact[c(7,1,6,4,3,2,9,5,10,8,18,14)]))%>%
  ggplot(aes(y = name, x = value, fill = category))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_wrap(facets = vars(impact), scales = 'free_y', ncol = 2)+
  scale_fill_nejm(name = "pandemic impact")+
  ylab('')+xlab('% change in cumulative SARS-CoV-2 infection incidence')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-100, 300))
  #scale_x_continuous(limits = quantile((data_deltas_raw$V_inc_cumul_pa_pa-1)*100, c(0.025, 0.975)))

p_delta_B_inc_routes = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog, Tr_inc_cumul_pa_pe, Tr_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[c(1:10,14,18)])%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')),
         name = factor(name, 
                       levels = c("Cr_inc_cumul_endog", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Tr_inc_cumul_pa_pe", "Tr_inc_cumul_pe_pe"),
                       labels = c("endogenous","patient to patient", "HCW to patient", "patient to HCW", "HCW to HCW")),
         impact = factor(impact, 
                         levels = df_pandI$impact[c(7,1,6,4,3,2,9,5,10,8,18,14)]))%>%
  ggplot(aes(y = name, x = value, fill = category))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_wrap(facets = vars(impact), scales = 'free_y', ncol = 2)+
  scale_fill_nejm(name = "pandemic impact")+
  ylab('')+xlab('% change in cumulative MRB incidence')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-100, 300))
  #scale_x_continuous(limits = quantile((data_deltas_raw$Cr_inc_cumul_endog-1)*100, c(0.001, 0.975)))

ggsave(p_delta_V_inc_routes, file = paste0("plots", "/deltaV_inc_routes.pdf"), width = 20, height = 25, unit = 'cm')
ggsave(p_delta_V_inc_routes, file = paste0("plots", "/deltaV_inc_routes.png"), width = 20, height = 25, unit = 'cm')

ggsave(p_delta_B_inc_routes, file = paste0("plots", "/deltaB_inc_routes.pdf"), width = 20, height = 28, unit = 'cm')
ggsave(p_delta_B_inc_routes, file = paste0("plots", "/deltaB_inc_routes.png"), width = 20, height = 28, unit = 'cm')


### MIXED INDICATORS FOR ALL TAU, AND FOR THE FOUR TAU HAVING THE GREATEST AND LEAST IMPACT

### SARS-CoV-2 incidence 
# total, then all 4 routes
p_delta_mixed_Vinc = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_all, V_inc_cumul_pa_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pa, V_inc_cumul_pe_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(route = factor(name, 
                        levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe"),
                        labels = c("total", "patient to patient", "patient to HCW", "HCW to patient", "HCW to HCW")),
         state = "SARS-CoV-2 infection")%>%
  ggplot(aes(y = route, x = value, fill = route))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#ef3b2c", "#fc9272", "#fc9272", "#fc9272", "#fc9272"))+
  ylab('')+xlab('% change in cumulative infection incidence (SARS-CoV-2)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(-100,-80,-60,-40,-20,0,20, 40, 60))+
  ggtitle(expression(paste(bold("a."))))



# total, for two greatest and two least
df_sort_Vinc = data_deltas_raw%>%
  filter(impact %in% df_pandI$impact[1:10])%>%
  select(impact, V_inc_cumul_all)%>%
  group_by(impact)%>%
  summarise(median = median(V_inc_cumul_all),
            mean = mean(V_inc_cumul_all))%>%
  arrange(mean)

toptwo_Vinc = as.character(df_sort_Vinc$impact[1:2])
bottomtwo_Vinc = as.character(df_sort_Vinc$impact[9:10])
  

p_delta_mixed_Vinc_topresponses = data_deltas_raw%>%dplyr::select(n_sim, tau, category, V_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% c(toptwo_Vinc, bottomtwo_Vinc))%>%
  mutate(state = "SARS-CoV-2 infection",
         direction = factor(impact,
                            levels = c(toptwo_Vinc, bottomtwo_Vinc),
                            labels = c("responses causing greatest decrease", "responses causing greatest decrease",
                                       "responses causing greatest increase", "responses causing greatest increase")))%>%
  ggplot(aes(y = impact, x = value, fill = direction))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_wrap(facets = vars(direction), nrow = 2, scales = 'free_y')+
  scale_fill_manual(values = c("#878787", "#b2182b"))+
  ylab('')+xlab('change in cumulative incidence (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(-100,-75,-50,-25, 0, 25,50,75, 100))+
  coord_cartesian(xlim = c(-100, 100))


p_grid_Vinc = plot_grid(p_delta_mixed_Vinc, p_delta_mixed_Vinc_topresponses, ncol = 2, align = "h", axis = "bt")


### Cr incidence (total, then all 3 routes)
p_delta_mixed_Binc = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_inc_cumul_all, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(name = factor(name, 
                       levels = c("Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog"),
                       labels = c("total", "patient to patient", "HCW to patient", "endogenous")),
         state = "patient colonization")%>%
  ggplot(aes(y = name, x = value, fill = name))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#2171b5", "#9ecae1", "#9ecae1", "#9ecae1"))+
  ylab('')+xlab('% change in cumulative patient colonization incidence (MRB)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600))+
  coord_cartesian(xlim = c(-100, 300))+
  ggtitle(expression(paste(bold("b. "))))

# total, for two greatest and two least
df_sort_Binc = data_deltas_raw%>%
  filter(impact %in% df_pandI$impact[1:10])%>%
  select(impact, Cr_inc_cumul_all)%>%
  group_by(impact)%>%
  summarise(median = median(Cr_inc_cumul_all),
            mean = mean(Cr_inc_cumul_all))%>%
  arrange(mean)

toptwo_Binc = as.character(df_sort_Binc$impact[1:2])
bottomtwo_Binc = as.character(df_sort_Binc$impact[9:10])

p_delta_mixed_Binc_topresponses = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Cr_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% c(toptwo_Binc, bottomtwo_Binc))%>%
  mutate(state = "patient colonization",
         direction = factor(impact,
                            levels = c(toptwo_Binc, bottomtwo_Binc),
                            labels = c("responses causing greatest decrease", "responses causing greatest decrease",
                                       "responses causing greatest increase", "responses causing greatest increase")))%>%
  ggplot(aes(y = impact, x = value, fill = direction))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_wrap(facets = vars(direction), nrow = 2, scales = 'free_y')+
  scale_fill_manual(values = c("#878787", "#b2182b"))+
  ylab('')+xlab('change in cumulative incidence (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(-100,-75,-50,-25, 0, 25,50,75, 100))+
  coord_cartesian(xlim = c(-100, 100))


p_grid_Binc = plot_grid(p_delta_mixed_Binc, p_delta_mixed_Binc_topresponses, ncol = 2, align = "h", axis = "bt")

### Tr incidence (total, then both routes)
p_delta_mixed_Tinc = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Tr_inc_cumul_all, Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(name = factor(name, 
                       levels = c("Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe"),
                       labels = c("total", "HCW to HCW", "patient to HCW")),
         state = "HCW carriage")%>%
  ggplot(aes(y = name, x = value, fill = name))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#6a51a3", "#bcbddc", "#bcbddc", "#bcbddc"))+
  ylab('')+xlab('% change in cumulative HCW carriage incidence (MRB)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  ggtitle(expression(paste(bold("c.  "))))

# total, for two greatest and two least
df_sort_Tinc = data_deltas_raw%>%
  filter(impact %in% df_pandI$impact[1:10])%>%
  select(impact, Tr_inc_cumul_all)%>%
  group_by(impact)%>%
  summarise(median = median(Tr_inc_cumul_all),
            mean = mean(Tr_inc_cumul_all))%>%
  arrange(mean)

toptwo_Tinc = as.character(df_sort_Tinc$impact[1:2])
bottomtwo_Tinc = as.character(df_sort_Tinc$impact[9:10])

p_delta_mixed_Tinc_topresponses = data_deltas_raw%>%dplyr::select(n_sim, tau, category, Tr_inc_cumul_all)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% c(toptwo_Tinc, bottomtwo_Tinc))%>%
  mutate(state = "HCW carriage (resistant bacteria)",
         direction = factor(impact,
                            levels = c(toptwo_Tinc, bottomtwo_Tinc),
                            labels = c("responses causing greatest decrease", "responses causing greatest decrease",
                                       "responses causing greatest increase", "responses causing greatest increase")))%>%
  ggplot(aes(y = impact, x = value, fill = direction))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_wrap(facets = vars(direction), nrow = 2, scales = 'free_y')+
  scale_fill_manual(values = c("#878787", "#b2182b"))+
  ylab('')+xlab('change in cumulative incidence (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = c(-100,-75,-50,-25, 0, 25,50,75, 100))+
  coord_cartesian(xlim = c(-100, 100))



p_delta_mixed = plot_grid(p_delta_mixed_Vinc+scale_x_continuous(breaks = c(-100,0,100,200,300))+coord_cartesian(xlim = c(-100,250)), 
                          p_delta_mixed_Binc+scale_x_continuous(breaks = c(-100,0,100,200,300))+coord_cartesian(xlim = c(-100,250)),
                          p_delta_mixed_Tinc+scale_x_continuous(breaks = c(-100,0,100,200,300))+coord_cartesian(xlim = c(-100,250)),
                          nrow = 3, rel_heights = c(6,5,4), align = 'v', axis = 'lr')

ggsave(p_delta_mixed, file = "plots/delta_inc_mixed.pdf", width = 14, height = 16, unit = 'cm')
ggsave(p_delta_mixed, file = "plots/delta_inc_mixed.png", width = 14, height = 16, unit = 'cm')



### all in one plot
p_delta_mixed_indicators = data_deltas_raw%>%dplyr::select(n_sim, tau, category, 
                                V_inc_cumul_all, V_inc_cumul_pa_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pa, V_inc_cumul_pe_pe,
                                R_rate_pd, Cr_pd,
                                Cr_inc_cumul_all, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog,
                                Tr_inc_cumul_all, Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(name = factor(name,
                       levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe",
                                  "Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                  "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe",
                                  "R_rate_pd", "Cr_pd")),
         route = factor(name, 
                       levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe",
                                  "Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                  "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe",
                                  "R_rate_pd", "Cr_pd"),
                       labels = c("total", "patient to patient", "patient to HCW", "HCW to patient", "HCW to HCW",
                                  "total", "patient to patient", "HCW to patient", "endogenous",
                                  "total", "HCW to HCW", "patient to HCW",
                                  "average resistance rate", "cumulative patient-days colonized")),
         indicator = factor(name,
                            levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe",
                                       "Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                       "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe",
                                       "R_rate_pd", "Cr_pd"),
                            labels = c(rep("SARS-CoV-2\ninfection\nincidence", 5),
                                       rep("patient\ncolonization\nincidence", 4),
                                       rep("HCW\ncarriage\nincidence", 3),
                                       rep("global\nresistance\nindicator",2))))%>%
  ggplot(aes(y = route, x = value, fill = name))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#ef3b2c", "#fc9272", "#fc9272", "#fc9272", "#fc9272",
                                                         "#2171b5", "#9ecae1", "#9ecae1", "#9ecae1",
                                                         "#6a51a3", "#bcbddc", "#bcbddc",
                                                         "#d94801", "#fec44f"))+
  ylab('')+xlab('change due to COVID-19 response (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  facet_grid(facets = vars(indicator), rows = 4, scales = 'free_y', space = 'free_y')+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-100,150))

ggsave(p_delta_mixed_indicators, file = "plots/delta_indicator_mixed.pdf", width = 20, height = 15, unit = 'cm')
ggsave(p_delta_mixed_indicators, file = "plots/delta_indicator_mixed.png", width = 20, height = 15, unit = 'cm')



### incidence (left) vs. resistance rate on its own (right)
p_delta_mixed_incidence = data_deltas_raw%>%dplyr::select(n_sim, tau, category, 
                                                           V_inc_cumul_all, V_inc_cumul_pa_pa, V_inc_cumul_pa_pe, V_inc_cumul_pe_pa, V_inc_cumul_pe_pe,
                                                           Cr_inc_cumul_all, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog,
                                                           Tr_inc_cumul_all, Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(name = factor(name,
                       levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe",
                                  "Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                  "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe")),
         route = factor(name, 
                        levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe",
                                   "Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                   "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe"),
                        labels = c("total", "patient to patient", "patient to HCW", "HCW to patient", "HCW to HCW",
                                   "total", "patient to patient", "HCW to patient", "endogenous",
                                   "total", "HCW to HCW", "patient to HCW")),
         indicator = factor(name,
                            levels = c("V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe",
                                       "Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                       "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe"),
                            labels = c(rep("SARS-CoV-2\ninfection\nincidence", 5),
                                       rep("patient\ncolonization\nincidence", 4),
                                       rep("HCW\ncarriage\nincidence", 3))))%>%
  ggplot(aes(y = route, x = value, fill = name))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#31a354", "#c7e9c0", "#c7e9c0", "#c7e9c0", "#c7e9c0",
                                                         "#2171b5", "#9ecae1", "#9ecae1", "#9ecae1",
                                                         "#6a51a3", "#bcbddc", "#bcbddc"))+
  ylab('')+xlab('change in incidence due to COVID-19 response (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  facet_grid(facets = vars(indicator), rows = 3, scales = 'free_y', space = 'free_y')+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-100,150))


p_delta_mixed_resistance = data_deltas_raw%>%dplyr::select(n_sim, tau, category, R_rate_pd)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  ggplot(aes(y = value, x = NA))+
  geom_boxplot(notch = T, outlier.alpha = 0.4, fill = "#ef3b2c")+
  theme_bw()+
  ylab('change in average resistance rate due to COVID-19 response (%)')+xlab('')+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-10,-1,0,1, 10, 100))
  
p_delta_mixed_indicators_incidence_resistance = plot_grid(p_delta_mixed_incidence+ggtitle(expression(paste(bold(a.)))), 
          p_delta_mixed_resistance+ggtitle(expression(paste(bold(b.)))),
                          ncol = 2, rel_widths = c(5,1), align = 'h', axis = 'tb')

ggsave(p_delta_mixed_indicators_incidence_resistance, file = "plots/delta_indicator_mixed_inc_res.pdf", width = 25, height = 15, unit = 'cm')
ggsave(p_delta_mixed_indicators_incidence_resistance, file = "plots/delta_indicator_mixed_inc_res.png", width = 25, height = 15, unit = 'cm')


############################################
### INCIDENCE ABSOLUTE VS INCIDENCE RATE ###
############################################

p_delta_incidence_absolute = data_deltas_raw%>%dplyr::select(n_sim, tau, category, 
                                                             Cr_inc_cumul_all, Cr_inc_cumul_pa_pa, Cr_inc_cumul_pe_pa, Cr_inc_cumul_endog,
                                                             Tr_inc_cumul_all, Tr_inc_cumul_pe_pe, Tr_inc_cumul_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(name = factor(name,
                       levels = c("Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                  "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe")),
         route = factor(name, 
                        levels = c("Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                   "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe"),
                        labels = c("total", "patient to patient", "HCW to patient", "endogenous",
                                   "total", "HCW to HCW", "patient to HCW")),
         indicator = factor(name,
                            levels = c("Cr_inc_cumul_all", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_inc_cumul_endog",
                                       "Tr_inc_cumul_all", "Tr_inc_cumul_pe_pe", "Tr_inc_cumul_pa_pe"),
                            labels = c(rep("patient colonization incidence", 4),
                                       rep("HCW carriage incidence", 3))))%>%
  ggplot(aes(y = route, x = value, fill = name))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#2171b5", "#9ecae1", "#9ecae1", "#9ecae1",
                                                         "#6a51a3", "#bcbddc", "#bcbddc"))+
  ylab('')+xlab('% change due to COVID-19 response (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  facet_grid(facets = vars(indicator), rows = 3, scales = 'free_y', space = 'free_y')+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-100,150))


p_delta_incidence_rate = data_deltas_raw%>%dplyr::select(n_sim, tau, category, 
                                Cr_incRate_all, Cr_incRate_pa_pa, Cr_incRate_pe_pa, Cr_incRate_endog,
                                Tr_incRate_all, Tr_incRate_pe_pe, Tr_incRate_pa_pe)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_pandI)%>%
  ungroup()%>%
  filter(impact %in% df_pandI$impact[12])%>%
  mutate(name = factor(name,
                       levels = c("Cr_incRate_all", "Cr_incRate_pa_pa", "Cr_incRate_pe_pa", "Cr_incRate_endog",
                                  "Tr_incRate_all", "Tr_incRate_pe_pe", "Tr_incRate_pa_pe")),
         route = factor(name, 
                        levels = c("Cr_incRate_all", "Cr_incRate_pa_pa", "Cr_incRate_pe_pa", "Cr_incRate_endog",
                                   "Tr_incRate_all", "Tr_incRate_pe_pe", "Tr_incRate_pa_pe"),
                        labels = c("total", "patient to patient", "HCW to patient", "endogenous",
                                   "total", "HCW to HCW", "patient to HCW")),
         indicator = factor(name,
                            levels = c("Cr_incRate_all", "Cr_incRate_pa_pa", "Cr_incRate_pe_pa", "Cr_incRate_endog",
                                       "Tr_incRate_all", "Tr_incRate_pe_pe", "Tr_incRate_pa_pe"),
                            labels = c(rep("patient colonization incidence", 4),
                                       rep("HCW carriage incidence", 3))))%>%
  ggplot(aes(y = route, x = value, fill = name))+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  scale_fill_manual(name = "pandemic impact", values = c("#2171b5", "#9ecae1", "#9ecae1", "#9ecae1",
                                                         "#6a51a3", "#bcbddc", "#bcbddc"))+
  ylab('')+xlab('% change due to COVID-19 response (%)')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  facet_grid(facets = vars(indicator), rows = 3, scales = 'free_y', space = 'free_y')+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-100,150))

p_delta_incidence_absolute_vs_rate = plot_grid(p_delta_incidence_absolute+ggtitle("cumulative incidence"),
                                               p_delta_incidence_rate+ggtitle("incidence rate"),
                                               ncol = 2)



# save just cumulative incidence
ggsave(p_delta_incidence_absolute, file = "plots/delta_inc_cumul_bacteria_allresponses.pdf", width = 20, height = 10, unit = 'cm')
ggsave(p_delta_incidence_absolute, file = "plots/delta_inc_cumul_bacteria_allresponses.png", width = 20, height = 10, unit = 'cm')


ggsave(p_delta_incidence_absolute_vs_rate, file = "plots/delta_inc_cumul_vs_rate_bacteria_allresponses.pdf", width = 35, height = 15, unit = 'cm')
ggsave(p_delta_incidence_absolute_vs_rate, file = "plots/delta_inc_cumul_vs_rate_bacteria_allresponses.png", width = 35, height = 15, unit = 'cm')




##################################################
### METADATA + INCIDENCE + RESISTANCE BOXPLOTS ###
##################################################

vec_metadata_totals = c("meta_abx_total", "meta_kappa_patients_total", "meta_kappa_staff_total", "meta_hh_total", 
                        "meta_adm_total", "meta_adm_R_total","meta_patientdays_total", "meta_staffdays_total", "meta_staffing_ratio_total",
                        "Tr_incRate_pe_pe", "Tr_incRate_pa_pe",
                        "Cr_incRate_pa_pa", "Cr_incRate_pe_pa", "Cr_incRate_endog",
                        "Cr_pd", "Cr_inc_cumul_all", "Cr_prev_max", "R_rate_pd", "R_rate_max")

df_metadata_labels = data.frame(name = vec_metadata_totals,
                                name_clean = c("patient-days of antibiotic exposure", "patient contacts", "HCW contacts","HCW hand decontamination rate", 
                                               "admissions","admissions colonized with MRB",
                                               "patient-days", "HCW-days", "HCW:patient ratio",
                                               "HCW to HCW", "patient to HCW",
                                               "patient to patient", "HCW to patient", "endogenous",
                                               "cumulative patient-days colonized", "cumulative colonization incidence", "colonization prevalence (daily max)", "average resistance rate", "resistance rate (daily max)"),
                                group = c(rep("behaviour",4),
                                          rep("demography",5),
                                          rep("MRB incidence rate", 5),
                                          rep("MRB burden", 5)))

vec_metadata_colours = c("#984ea3", "#377eb8", "#ff7f00","#e41a1c")

p_metadata_incidence_resistance = data_deltas_raw%>%
  filter(par == "combo_all")%>%
  dplyr::select(n_sim, tau, category, df_metadata_labels$name)%>%
  pivot_longer(-c(n_sim, tau, category))%>%
  mutate(value = (value - 1)*100,
         par = tau)%>%
  left_join(df_metadata_labels)%>%
  ungroup()%>%
  mutate(group = factor(group,
                        levels = c('behaviour', 'demography', 'FOI', 'MRB incidence rate', 'MRB burden')),
         name_clean = factor(name_clean,
                             levels = c("patient-days of antibiotic exposure", "HCW hand decontamination rate", "HCW contacts", "patient contacts",
                                        "admissions colonized with MRB", "admissions", "HCW:patient ratio", "HCW-days", "patient-days",
                                        "endogenous", "patient to patient", "HCW to patient", 
                                        "HCW to HCW", "patient to HCW",
                                        "resistance rate (daily max)", "average resistance rate", "colonization prevalence (daily max)", "cumulative colonization incidence", "cumulative patient-days colonized")))%>%
  ggplot(aes(y = name_clean, x = value, fill = group))+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_boxplot(notch = T, outlier.alpha = 0.4)+
  theme_bw()+
  facet_grid(rows = vars(group), scales = 'free_y', space = "free_y")+
  ylab('')+xlab('% change due to COVID-19 responses')+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100, -10,-1,0,1, 10, 100, 1000))+
  scale_fill_manual(values = vec_metadata_colours)+
  theme(legend.position = "none")

ggsave(p_metadata_incidence_resistance, file = "plots/delta_metadata_incidence_resistance.pdf", width = 20, height = 15, unit = 'cm')
ggsave(p_metadata_incidence_resistance, file = "plots/delta_metadata_incidence_resistance.png", width = 20, height = 15, unit = 'cm')



#############################################
### SIMULTANEOUS INDICATORS SCATTER PLOTS ###
#############################################

#############################################
### Delta incidence (Cr) by incidence (V) ###
#############################################
### policy vs. surge impacts
p_deltaB_deltaV_inc_cumul_policy_surge = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(12,14,18)])%>%
  dplyr::select(n_sim, tau, category, V_inc_cumul_all, Cr_inc_cumul_all)%>%
  ggplot(aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_colour_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_shape_manual(name = "COVID-19 response", values = c(22,25,23), labels = c("policy", "caseload", "all"))+
  ylab('% change in cumulative SARS-CoV-2 infection incidence')+
  xlab('% change in cumulative MRB colonization incidence')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100, 1000))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100,1000))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

p_deltaB_deltaV_inc_cumul_policy_surge_marginal = ggMarginal(p_deltaB_deltaV_inc_cumul_policy_surge+
                                                               theme(legend.position = 'bottom')+
                                                               scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10), limits = c(-100,50)), 
                                                             type = "histogram", groupFill = T, alpha = 1)



### 5 combos
p_deltaB_deltaV_inc_cumul_combos = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)])%>%
  dplyr::select(n_sim, tau, category, V_inc_cumul_all, Cr_inc_cumul_all)%>%
  ggplot(aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_nejm(name = "COVID-19 response", labels = labels_5combos)+
  scale_colour_nejm(name = "COVID-19 response", labels = labels_5combos)+
  scale_shape_manual(name = "COVID-19 response", values = c(21,22,23,24,25), labels = labels_5combos)+
  ylab('% change in cumulative SARS-CoV-2 infection incidence')+
  xlab('% change in cumulative MRB colonization incidence')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100, 1000))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(hjust = 0))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, fill = tau, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

p_deltaB_deltaV_inc_cumul_combos_marginal = ggMarginal(p_deltaB_deltaV_inc_cumul_combos+
                                                         guides(fill=guide_legend(nrow=3,byrow=TRUE))+
                                                         theme(legend.position = 'bottom')+
                                                         scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10,100), limits = c(-100,100)), 
                                                       type = "histogram", groupFill = T, alpha = 1)

### all ten
p_deltaB_deltaV_inc_cumul_allTau = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(1:10)])%>%
  dplyr::select(n_sim, impact, category, V_inc_cumul_all, Cr_inc_cumul_all)%>%
  ggplot(aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, colour = impact, fill = impact, shape = impact))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(name = "COVID-19 response", values = cols_nejm_paired, labels = labels_all10)+
  scale_colour_manual(name = "COVID-19 response", values = cols_nejm_paired, labels = labels_all10)+
  scale_shape_manual(name = "COVID-19 response", values = c(21,24,21,24,21,24,21,24,21,24), labels = labels_all10)+
  ylab('% change in cumulative SARS-CoV-2 infection incidence')+
  xlab('% change in cumulative MRB colonization incidence')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100, 1000))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100, 1000))+
  #guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(hjust = 0))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(1:10)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, fill = impact, shape = impact),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(1:10)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (V_inc_cumul_all-1)*100, fill = impact, shape = impact),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)


#################################################
### Delta incidence by resistance rate for Cr ###
#################################################

# policy vs. surge impacts
p_deltaBinccumul_deltaBrrate_policy_surge = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(12,14,18)])%>%
  dplyr::select(n_sim, tau, category, Cr_inc_cumul_all, R_rate_pd)%>%
  ggplot(aes(y = (R_rate_pd-1)*100, x = (Cr_inc_cumul_all-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_colour_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_shape_manual(name = "COVID-19 response", values = c(22,25,23), labels = c("policy", "caseload", "all"))+
  xlab('% change in cumulative MRB colonization incidence')+
  ylab('% change in average MRB resistance rate')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100,1000))+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, fill = tau, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

p_deltaBinccumul_deltaBrrate_policy_surge_marginal = ggMarginal(p_deltaBinccumul_deltaBrrate_policy_surge+
                                                                  theme(legend.position = 'bottom')+
                                                                  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-10,-1,0,1, 10, 100), limits = c(-35,275)),
                                                       type = "histogram", groupFill = T, alpha = 1)

# 5 combos
p_deltaBinccumul_deltaBrrate_combos = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)])%>%
  dplyr::select(n_sim, tau, category, Cr_inc_cumul_all, R_rate_pd)%>%
  ggplot(aes(y = (R_rate_pd-1)*100, x = (Cr_inc_cumul_all-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_nejm(name = "COVID-19 response", labels = labels_5combos)+
  scale_colour_nejm(name = "COVID-19 response", labels = labels_5combos)+
  scale_shape_manual(name = "COVID-19 response", values = c(21,22,23,24,25), labels = labels_5combos)+
  xlab('% change in cumulative MRB colonization incidence')+
  ylab('% change in average MRB resistance rate')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100,1000))+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(hjust = 0))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

p_deltaBinccumul_deltaBrrate_combos_marginal = ggMarginal(p_deltaBinccumul_deltaBrrate_combos+
                                                         guides(fill=guide_legend(nrow=3,byrow=TRUE))+
                                                         theme(legend.position = 'bottom'),
                                                       type = "histogram", groupFill = T, alpha = 1)

### all ten
p_deltaBinccumul_deltaBrrate_allTau = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(1:10)])%>%
  dplyr::select(n_sim, impact, category, Cr_inc_cumul_all, R_rate_pd)%>%
  ggplot(aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, colour = impact, fill = impact, shape = impact))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(name = "COVID-19 response", values = cols_nejm_paired, labels = labels_all10)+
  scale_colour_manual(name = "COVID-19 response", values = cols_nejm_paired, labels = labels_all10)+
  scale_shape_manual(name = "COVID-19 response", values = c(21,24,21,24,21,24,21,24,21,24), labels = labels_all10)+
  xlab('% change in cumulative MRB colonization incidence')+
  ylab('% change in average MRB resistance rate')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-10,-1,0,1, 10, 100))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100, 1000))+
  #guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(hjust = 0))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(1:10)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, fill = impact, shape = impact),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(1:10)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (R_rate_pd-1)*100, fill = impact, shape = impact),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

##############################################################
### Delta patient-days colonized by resistance rate for Cr ###
##############################################################

# policy vs. surge impacts
p_deltaBpd_deltaBrrate_policy_surge = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(12,14,18)])%>%
  dplyr::select(n_sim, tau, category, Cr_pd, R_rate_pd)%>%
  ggplot(aes(y = (R_rate_pd-1)*100, x = (Cr_pd-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_colour_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_shape_manual(name = "COVID-19 response", values = c(22,25,23), labels = c("policy", "caseload", "all"))+
  xlab('% change in cumulative patient-days colonized with MRB')+
  ylab('% change in average MRB resistance rate')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(y = (R_rate_pd-1)*100, x = (Cr_pd-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(y = (R_rate_pd-1)*100, x = (Cr_pd-1)*100, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

p_deltaBpd_deltaBrrate_policy_surge_marginal = ggMarginal(p_deltaBpd_deltaBrrate_policy_surge+
                                                            theme(legend.position = 'bottom')+
                                                            scale_y_continuous(trans=pseudolog10_trans, breaks = c(-10,-1,0,1, 10, 100), limits = c(-35,300)),
                                                          type = "histogram", groupFill = T, alpha = 1)

# 5 combos
p_deltaBpd_deltaBrrate_combos = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)])%>%
  dplyr::select(n_sim, tau, category, Cr_pd, R_rate_pd)%>%
  ggplot(aes(y = (R_rate_pd-1)*100, x = (Cr_pd-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_nejm(name = "COVID-19 response", labels = labels_5combos)+
  scale_colour_nejm(name = "COVID-19 response", labels = labels_5combos)+
  scale_shape_manual(name = "COVID-19 response", values = c(21,22,23,24,25), labels = labels_5combos)+
  xlab('% change in cumulative patient-days colonized with MRB')+
  ylab('% change in average MRB resistance rate')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(hjust = 0))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)]),
             mapping = aes(y = (R_rate_pd-1)*100, x = (Cr_pd-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(11,13,15,16,17)]),
             mapping = aes(y = (R_rate_pd-1)*100, x = (Cr_pd-1)*100, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

p_deltaBpd_deltaBrrate_combos_marginal = ggMarginal(p_deltaBpd_deltaBrrate_combos+
                                                            guides(fill=guide_legend(nrow=3,byrow=TRUE))+
                                                            theme(legend.position = 'bottom')+
                                                      scale_y_continuous(trans=pseudolog10_trans, breaks = c(-10,-1,0,1, 10, 100), limits = c(-30,300)),
                                                          type = "histogram", groupFill = T, alpha = 1)


### all ten
p_deltaBpd_deltaBrrate_allTau = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(1:10)])%>%
  dplyr::select(n_sim, impact, category, Cr_pd, R_rate_pd)%>%
  ggplot(aes(x = (Cr_pd-1)*100, y = (R_rate_pd-1)*100, colour = impact, fill = impact, shape = impact))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(name = "COVID-19 response", values = cols_nejm_paired, labels = labels_all10)+
  scale_colour_manual(name = "COVID-19 response", values = cols_nejm_paired, labels = labels_all10)+
  scale_shape_manual(name = "COVID-19 response", values = c(21,24,21,24,21,24,21,24,21,24), labels = labels_all10)+
  xlab('% change in cumulative patient-days colonized with MRB')+
  ylab('% change in average MRB resistance rate')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100, 1000))+
  #guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(hjust = 0))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(1:10)]),
             mapping = aes(x = (Cr_pd-1)*100, y = (R_rate_pd-1)*100, fill = impact, shape = impact),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(1:10)]),
             mapping = aes(x = (Cr_pd-1)*100, y = (R_rate_pd-1)*100, fill = impact, shape = impact),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)

#######################################################
### Delta patient-days colonized by incidence for Cr ###
#######################################################
p_deltaBpd_deltaBinccumul_policy_surge = data_deltas_raw%>%
  filter(tau %in% vec_pandI_parameters[c(12,14,18)])%>%
  dplyr::select(n_sim, tau, category, Cr_pd, Cr_inc_cumul_all)%>%
  ggplot(aes(x = (Cr_inc_cumul_all-1)*100, y = (Cr_pd-1)*100, colour = tau, fill = tau, shape = tau))+
  geom_point(alpha = 0.3)+
  theme_bw()+
  scale_fill_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_colour_manual(values = cols_policy_surge, name = "COVID-19 response", labels = c("policy", "caseload", "all"))+
  scale_shape_manual(name = "COVID-19 response", values = c(22,25,23), labels = c("policy", "caseload", "all"))+
  ylab('% change in cumulative patient-days colonized with MRB')+
  xlab('% change in cumulative MRB colonization incidence')+
  geom_vline(xintercept = 0, linetype = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = 1, alpha = 1)+
  scale_y_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  scale_x_continuous(trans=pseudolog10_trans, breaks = c(-100,-10,-1,0,1, 10, 100))+
  guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(alpha = 1)))+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (Cr_pd-1)*100, fill = tau, shape = tau),
             colour = 'black', size = 4, alpha = 1, stroke = 1)+
  geom_point(data_deltas_means%>%filter(tau %in% vec_pandI_parameters[c(12,14,18)]),
             mapping = aes(x = (Cr_inc_cumul_all-1)*100, y = (Cr_pd-1)*100, colour = tau, fill = tau, shape = tau),
             colour = 'black', fill = 'black', size = 0.5, alpha = 1)


### combine plots
p_BbyV_5combos = ggarrange(p_deltaB_deltaV_inc_cumul_combos+
                             coord_cartesian(ylim = c(-100,275), xlim = c(-100,450))+
                             ggtitle(expression(paste(bold(a.)))), 
                           p_deltaBinccumul_deltaBrrate_combos+
                             coord_cartesian(ylim = c(-100,275), xlim = c(-100,450))+
                             ggtitle(expression(paste(bold(b.)))), 
                           common.legend = T, legend = 'bottom', ncol = 2)

p_BbyV_policy_surge = ggarrange(p_deltaB_deltaV_inc_cumul_policy_surge+
                                  coord_cartesian(ylim = c(-100,220), xlim = c(-100,275))+
                                  ggtitle(expression(paste(bold(a.)))), 
                                p_deltaBpd_deltaBrrate_policy_surge+
                                  coord_cartesian(ylim = c(-100,220), xlim = c(-100,275))+
                                  ggtitle(expression(paste(bold(b.)))), 
                                common.legend = T, legend = 'bottom', ncol = 2)


p_deltaBandV_allTau = ggarrange(p_deltaB_deltaV_inc_cumul_allTau+
                                  coord_cartesian(ylim = c(-100,275), xlim = c(-100,500))+
                                  ggtitle(expression(paste(bold(a.)))), 
                                p_deltaBpd_deltaBrrate_allTau+
                                  coord_cartesian(ylim = c(-100,275), xlim = c(-100,500))+
                                  ggtitle(expression(paste(bold(b.)))), 
                                common.legend = T, legend = 'bottom', ncol = 2)

p_deltaBpd_deltaBrrate_alltau_policy_surge = plot_grid(p_deltaBpd_deltaBrrate_allTau+
                                                         coord_cartesian(ylim = c(-50,275), xlim = c(-100,250))+
                                                         ggtitle(expression(paste(bold(a.)))), 
                                                       p_deltaBpd_deltaBrrate_policy_surge+
                                                         coord_cartesian(ylim = c(-50,275), xlim = c(-100,250))+
                                                         ggtitle(expression(paste(bold(b.)))), 
                                                       nrow = 2,
                                                       align = "v", axis = "lr")

p_deltaBinc_deltaVinc_alltau_policy_surge = plot_grid(p_deltaB_deltaV_inc_cumul_allTau+
                                                         coord_cartesian(ylim = c(-100,60), xlim = c(-100,450))+
                                                         ggtitle(expression(paste(bold(a.)))), 
                                                      p_deltaB_deltaV_inc_cumul_policy_surge+
                                                         coord_cartesian(ylim = c(-100,60), xlim = c(-100,450))+
                                                         ggtitle(expression(paste(bold(b.)))), 
                                                       nrow = 2,
                                                       align = "v", axis = "lr")

p_5combos_allindicators = ggarrange(p_deltaB_deltaV_inc_cumul_combos+
                                      coord_cartesian(ylim = c(-100,275), xlim = c(-100,450))+
                                      ggtitle(expression(paste(bold(a.)))),
                                    p_deltaBpd_deltaBrrate_combos+
                                      coord_cartesian(ylim = c(-100,275), xlim = c(-100,450))+
                                      ggtitle(expression(paste(bold(b.)))), 
                                    common.legend = T, legend = 'right', nrow = 2)





### save ###

### individual plots
ggsave(p_deltaB_deltaV_inc_cumul_policy_surge_marginal, file = "plots/scatter_deltaB_deltaV_inc_cumul_policy_surge_marginal.pdf", width = 15, height = 16, unit = 'cm')
ggsave(p_deltaB_deltaV_inc_cumul_policy_surge_marginal, file = "plots/scatter_deltaB_deltaV_inc_cumul_policy_surge_marginal.png", width = 15, height = 16, unit = 'cm')

ggsave(p_deltaB_deltaV_inc_cumul_combos_marginal, file = "plots/scatter_deltaB_deltaV_inc_cumul_5combos_marginal.pdf", width = 15, height = 18, unit = 'cm')
ggsave(p_deltaB_deltaV_inc_cumul_combos_marginal, file = "plots/scatter_deltaB_deltaV_inc_cumul_5combos_marginal.png", width = 15, height = 18, unit = 'cm')

ggsave(p_deltaBinccumul_deltaBrrate_policy_surge_marginal, file = "plots/scatter_deltaBinccumul_deltaBrrate_policy_surge_marginal.pdf", width = 15, height = 16, unit = 'cm')
ggsave(p_deltaBinccumul_deltaBrrate_policy_surge_marginal, file = "plots/scatter_deltaBinccumul_deltaBrrate_policy_surge_marginal.png", width = 15, height = 16, unit = 'cm')

ggsave(p_deltaBinccumul_deltaBrrate_combos_marginal, file = "plots/scatter_deltaBinccumul_deltaBrrate_combos_marginal.pdf", width = 15, height = 18, unit = 'cm')
ggsave(p_deltaBinccumul_deltaBrrate_combos_marginal, file = "plots/scatter_deltaBinccumul_deltaBrrate_combos_marginal.png", width = 15, height = 18, unit = 'cm')

ggsave(p_deltaBpd_deltaBrrate_policy_surge_marginal, file = "plots/scatter_deltaBpd_deltaBrrate_policy_surge_marginal.pdf", width = 15, height = 16, unit = 'cm')
ggsave(p_deltaBpd_deltaBrrate_policy_surge_marginal, file = "plots/scatter_deltaBpd_deltaBrrate_policy_surge_marginal.png", width = 15, height = 16, unit = 'cm')

ggsave(p_deltaBpd_deltaBrrate_combos_marginal, file = "plots/scatter_deltaBpd_deltaBrrate_combos_marginal.pdf", width = 15, height = 18, unit = 'cm')
ggsave(p_deltaBpd_deltaBrrate_combos_marginal, file = "plots/scatter_deltaBpd_deltaBrrate_combos_marginal.png", width = 15, height = 18, unit = 'cm')


### combined plots
ggsave(p_BbyV_5combos, file = paste0("plots", "/scatter_deltaB_deltaV_5combos.pdf"), width = 25, height = 15, unit = 'cm')
ggsave(p_BbyV_5combos, file = paste0("plots", "/scatter_deltaB_deltaV_5combos.png"), width = 25, height = 15, unit = 'cm')

ggsave(p_BbyV_policy_surge, file = paste0("plots", "/scatter_deltaB_deltaV_policy_surge.pdf"), width = 25, height = 15, unit = 'cm')
ggsave(p_BbyV_policy_surge, file = paste0("plots", "/scatter_deltaB_deltaV_policy_surge.png"), width = 25, height = 15, unit = 'cm')

ggsave(p_deltaBandV_allTau, file = paste0("plots", "/scatter_deltaB_deltaV_allTau.pdf"), width = 25, height = 15, unit = 'cm')
ggsave(p_deltaBandV_allTau, file = paste0("plots", "/scatter_deltaB_deltaV_allTau.png"), width = 25, height = 15, unit = 'cm')

ggsave(p_deltaBpd_deltaBrrate_alltau_policy_surge, file = paste0("plots", "/scatter_deltaB_alltau_policy_surge.pdf"), width = 18, height = 25, unit = 'cm')
ggsave(p_deltaBpd_deltaBrrate_alltau_policy_surge, file = paste0("plots", "/scatter_deltaB_alltau_policy_surge.png"), width = 18, height = 25, unit = 'cm')

ggsave(p_deltaBinc_deltaVinc_alltau_policy_surge, file = paste0("plots", "/scatter_deltaB_deltaV_alltau_policy_surge.pdf"), width = 18, height = 25, unit = 'cm')
ggsave(p_deltaBinc_deltaVinc_alltau_policy_surge, file = paste0("plots", "/scatter_deltaB_deltaV_alltau_policy_surge.png"), width = 18, height = 25, unit = 'cm')

ggsave(p_5combos_allindicators, file = paste0("plots", "/scatter_allindicators_5combos.pdf"), width = 18, height = 25, unit = 'cm')
ggsave(p_5combos_allindicators, file = paste0("plots", "/scatter_allindicators_5combos.png"), width = 18, height = 25, unit = 'cm')

################################
### BOOTSTRAP DELTA OUTCOMES ###
################################

df_bootstrap = loadRData("data/outputs_analysis2/indicators_deltas_bootstrapped.Rdata")%>%
  mutate(par = tau)%>%
  left_join(df_pandI)%>%
  mutate(category = factor(category, 
                           levels = c('abx', 'contact', 'ipc', 'disease', 'admission', 'mixed'), 
                           labels = c('antibiotics', 'contact', 'IPC', 'disease', 'admission', 'all')))%>%
  mutate(impact = factor(impact,
                         levels = df_pandI$impact[c(1,7,13,
                                                    4,6,15,
                                                    2,3,17,
                                                    9,5,16,
                                                    8,10,11,
                                                    18,14,12)]))

p_bootstrap_V_inc_cumul_all = df_bootstrap%>%
  ggplot(aes(x = factor(n_sims), y = (mean_V_inc_cumul_all-1)*100, fill = category))+
  geom_boxplot()+
  facet_wrap(facets = vars(impact), scales = 'free_y', ncol = 3)+
  theme_bw()+
  xlab('number of simulations')+
  ylab('mean change in cumulative SARS-CoV-2 infection incidence (%)')+
  scale_fill_nejm()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1))


p_bootstrap_Cr_inc_cumul_all = df_bootstrap%>%
  ggplot(aes(x = factor(n_sims), y = (mean_Cr_inc_cumul_all-1)*100, fill = category))+
  geom_boxplot()+
  facet_wrap(facets = vars(impact), scales = 'free_y', ncol = 3)+
  theme_bw()+
  xlab('number of simulations')+
  ylab('mean change in cumulative MRB colonization incidence (%)')+
  scale_fill_nejm(name = "COVID-19 response")+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(p_bootstrap_V_inc_cumul_all, file = paste0("plots", "/bootstrap_Vinc.pdf"), height = 25, width = 25, unit = 'cm')
ggsave(p_bootstrap_V_inc_cumul_all, file = paste0("plots", "/bootstrap_Vinc.png"), height = 25, width = 25, unit = 'cm')

ggsave(p_bootstrap_Cr_inc_cumul_all, file = paste0("plots", "/bootstrap_Binc.pdf"), height = 25, width = 25, unit = 'cm')
ggsave(p_bootstrap_Cr_inc_cumul_all, file = paste0("plots", "/bootstrap_Binc.png"), height = 25, width = 25, unit = 'cm')


######################
### PRCC HEAT MAPS ###
######################
df_prcc = loadRData("data/outputs_analysis2/data_prcc.Rdata")%>%
  mutate(label = as.character(round(gamma, 2)))%>%
  mutate(label = ifelse(p.value <0.000005, paste0(label, "***"),
                        ifelse(p.value <0.0005, paste0(label, "**"),
                               ifelse(p.value<0.05, paste0(label, "*"), label))))%>%
  mutate(pathogen = case_when(grepl("V_", outcome) ~ "SARS-CoV-2 outcomes",
                              TRUE ~ "MRB outcomes"))%>%
  mutate(symbol = par_name)

df_prcc$outcome = factor(df_prcc$outcome,
                         levels = c("Cr_inc_cumul_all", "Cr_inc_cumul_endog", "Cr_inc_cumul_pa_pa", "Cr_inc_cumul_pe_pa", "Cr_pd", "R_rate_pd",
                                    "Tr_inc_cumul_all","Tr_inc_cumul_pa_pe", "Tr_inc_cumul_pe_pe",
                                    "V_inc_cumul_all", "V_inc_cumul_pa_pa", "V_inc_cumul_pa_pe", "V_inc_cumul_pe_pa", "V_inc_cumul_pe_pe"),
                         labels = c("colonization incidence (all)", "colonization incidence (endogenous)", "colonization incidence (patient to patient)", "colonization incidence (HCW to patient)", "patient-days colonized", "average resistance rate",
                                    "carriage incidence (all)", "carriage incidence (patient to HCW)", "carriage incidence (HCW to HCW)", 
                                    "infection incidence (all)", "infection incidence (patient to patient)", "infection incidence (patient to HCW)", "infection incidence (HCW to patient)", "infection incidence (HCW to HCW)"))

df_prcc$symbol  = factor(df_prcc$symbol,
                         levels = c("a", "alpha", "cost", "f_Cr", "f_Cs", "gamma_Cs", "hyg", 
                                    "kappa_pa_pa", "kappa_pa_pe", "kappa_pe_pe", "mu",
                                    "Nbeds", "Nhcw", "pi_Cr", "pi_S", "r_R", "r_S", "rho", "theta"),
                         labels = c(expression(A[base]), expression(paste(alpha,"'")), expression(c), expression(f_C^R), expression(f[C^S]), expression(gamma), expression(H[base]),
                                    expression(kappa^"pat->pat"), expression(kappa^"pat->hcw"), expression(kappa^"hcw->hcw"), expression(mu),
                                    expression(N^beds), expression(N^hcw), expression(pi[B]), expression(pi[V]), expression(r[R]), expression(r[S]), expression(rho), expression(theta)
                                    )
                         )

symbol_levels = c("a", "alpha", "cost", "f_Cr", "f_Cs", "gamma_Cs", "hyg", 
                  "kappa_pa_pa", "kappa_pa_pe", "kappa_pe_pe", "mu",
                  "Nbeds", "Nhcw", "pi_Cr", "pi_S", "r_R", "r_S", "rho", "theta")
symbol_labels = c(expression(A[base]), expression(paste(alpha,"'")), expression(c), expression(f[C^R]), expression(f[C^S]), expression(gamma), expression(H[base]),
                  expression(kappa^"pat->pat"), expression(kappa^"pat->hcw"), expression(kappa^"hcw->hcw"), expression(mu),
                  expression(N^beds), expression(N^hcw), expression(pi[B]), expression(pi[V]), expression(r[R]), expression(r[S]), expression(rho), expression(theta))

size_text = 3

### compare all main outcomes for COVID-19 response "none"
p_prcc_V_notau = df_prcc%>%
  filter(tau == 'none', pathogen == "SARS-CoV-2 outcomes")%>%
  ggplot(aes(x = par_name, y = outcome, fill = gamma))+
  geom_tile(color = "white")+
  geom_text(aes(x = par_name, y = outcome, label = label), colour = 'black', size = size_text, fontface = 'bold')+
  scale_fill_gradient2(low = wes_palette("Zissou1")[1],  mid = 'white', high = wes_palette("Zissou1")[5],
                       limits = c(-1,1), midpoint = 0, name = 'PRCC')+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 10, vjust = 0.5))+
  ylab('')+
  xlab('')+
  scale_x_discrete(labels = symbol_labels)

p_prcc_B_notau = df_prcc%>%
  filter(tau == 'none', pathogen == "MRB outcomes")%>%
  ggplot(aes(x = par_name, y = outcome, fill = gamma))+
  geom_tile(color = "white")+
  geom_text(aes(x = par_name, y = outcome, label = label), colour = 'black', size = size_text, fontface = 'bold')+
  scale_fill_gradient2(low = wes_palette("Zissou1")[1],  mid = 'white', high = wes_palette("Zissou1")[5],
                       limits = c(-1,1), midpoint = 0, name = 'PRCC')+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 10, vjust = 0.5))+
  ylab('')+
  xlab('')+
  scale_x_discrete(labels = symbol_labels)


# compare all main outcomes for combined COVID-19 response "combo_all"
p_prcc_V_alltau = df_prcc%>%
  filter(tau == 'combo_all', pathogen == "SARS-CoV-2 outcomes")%>%
  ggplot(aes(x = par_name, y = outcome, fill = gamma))+
  geom_tile(color = "white")+
  geom_text(aes(x = par_name, y = outcome, label = label), colour = 'black', size = size_text, fontface = 'bold')+
  scale_fill_gradient2(low = wes_palette("Zissou1")[1],  mid = 'white', high = wes_palette("Zissou1")[5],
                       limits = c(-1,1), midpoint = 0, name = 'PRCC')+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 10, vjust = 0.5))+
  ylab('')+
  xlab('')+
  scale_x_discrete(labels = symbol_labels)

p_prcc_B_alltau = df_prcc%>%
  filter(tau == 'combo_all', pathogen == "MRB outcomes")%>%
  ggplot(aes(x = par_name, y = outcome, fill = gamma))+
  geom_tile(color = "white")+
  geom_text(aes(x = par_name, y = outcome, label = label), colour = 'black', size = size_text, fontface = 'bold')+
  scale_fill_gradient2(low = wes_palette("Zissou1")[1],  mid = 'white', high = wes_palette("Zissou1")[5],
                       limits = c(-1,1), midpoint = 0, name = 'PRCC')+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 10, vjust = 0.5))+
  ylab('')+
  xlab('')+
  scale_x_discrete(labels = symbol_labels)


### combine and save
p_prcc_notau = plot_grid(p_prcc_V_notau+theme(legend.position = 'none')+ggtitle(expression(paste(bold("a.  "), "SARS-CoV-2 outcomes"))), 
                         p_prcc_B_notau+theme(legend.position = 'bottom')+ggtitle(expression(paste(bold("b.  "), "MRB outcomes"))), 
                         nrow = 2, ncol = 1, align = "v", axis = "lr", rel_heights = c(1,2))

p_prcc_alltau = plot_grid(p_prcc_V_alltau+theme(legend.position = 'none')+ggtitle(expression(paste(bold("a.  "), "SARS-CoV-2 outcomes"))), 
                          p_prcc_B_alltau+theme(legend.position = 'bottom')+ggtitle(expression(paste(bold("b.  "), "MRB outcomes"))), 
                          nrow = 2, ncol = 1, align = "v", axis = "lr", rel_heights = c(1,2))

ggsave(p_prcc_notau, file = paste0("plots", "/prcc_notau.pdf"), height = 25, width = 35, unit = 'cm')
ggsave(p_prcc_notau, file = paste0("plots", "/prcc_notau.png"), height = 25, width = 35, unit = 'cm')

ggsave(p_prcc_alltau, file = paste0("plots", "/prcc_alltau.pdf"), height = 25, width = 35, unit = 'cm')
ggsave(p_prcc_alltau, file = paste0("plots", "/prcc_alltau.png"), height = 25, width = 35, unit = 'cm')

