library(deSolve)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(shades)
library(RColorBrewer)
library(cowplot)

source("ODEs.R")
source("functions.R")
source("pars_states.R")

##################################################
### SUPPLEMENTARY FIGURE FOR MODEL ASSUMPTIONS ###
##################################################

# supplementary figure: visualizing these functions
df_hh_i = data.frame()
for(ratio_i in c(1,1.5,2)){
  print(ratio_i)
  for(H_i in seq(0.1,0.6,length.out = 21)){
    for(cc_i in seq(1,30,length.out = 21)){
      
      val_i=f_eta(H_i, cc_i, ratio_i)
      
      df_i=data.frame(ratio = ratio_i, H = H_i, cc = cc_i, val = val_i)
      df_hh_i = rbind(df_hh_i, df_i)
    }
  }
}

df_hh_i$ratio = factor(df_hh_i$ratio, 
                       levels = c(1,1.5,2),
                       labels = c("HCW:patient ratio = 1", 
                                  "HCW:patient ratio = 1.5", 
                                  "HCW:patient ratio = 2"))

p_hh = ggplot(df_hh_i, aes(x = cc, y = H, fill = val/24))+
  geom_tile()+
  theme_bw()+
  scico::scale_fill_scico(name = expression(paste(omega, " (/hour)")), palette = "lajolla")+
  ggtitle(expression(paste(bold("a.  "), "rate of HCW decontamination")))+
  ylab(expression(paste(H, ' (hand hygiene compliance)')))+
  xlab(expression(paste(kappa^"hcw->pat", ' (daily number of HCW contacts with patients)')))+
  theme(legend.title = element_text(hjust = 0))+
  facet_grid(rows = vars(ratio))+
  theme(legend.position = 'bottom')



df_tau_i = data.frame()
for(tau_i in seq(0,1,length.out = 21)){
  print(tau_i)
  for(infect_i in seq(0,0.3,length.out = 21)){
    val_i=f_surge_pars(tau_i, infect_i, N=1)
    df_i=data.frame(tau = tau_i, freq = infect_i, val = val_i)
    df_tau_i = rbind(df_tau_i, df_i)
  }
}

p_tau_i = ggplot(df_tau_i, aes(x = freq, y = tau, fill = val))+
  geom_tile()+
  theme_bw()+
  scico::scale_fill_scico(name = expression(paste(Tau[x](t))), palette = "tokyo")+
  ggtitle(expression(paste(bold("b.  "), 'prevalence-scaled COVID-19 responses (', tau[as],', ', tau[cd],', ', tau[ra],', ', tau[sc],')')))+
  ylab(expression(paste(tau, " (normalized strength of COVID-19 response)")))+
  xlab(expression(paste(x[t], ' [patient SARS-CoV-2 infection prevalence (', I^pat/N^pat, ')]')))+
  theme(legend.title = element_text(hjust = 0))+
  theme(legend.position = 'bottom')

p_hh_tau = plot_grid(p_hh, p_tau_i, ncol = 2, align = 'h', axis = 'tb', rel_widths = c(0.75,1))

ggsave("plots/assumptions_hh_tau.png", plot = p_hh_tau, width = 25, height = 15, unit = 'cm', )
ggsave("plots/assumptions_hh_tau.pdf", plot = p_hh_tau, width = 25, height = 15, unit = 'cm', )

##############################
### FIGURES FOR ANALYSIS 1 ###
##############################

###################################
### Dynamics over ranges of tau ###
###################################

### (1.1) For each pandemic impact, and combination of pandemic impacts, 
### how does varying tau impact SARS-CoV-2 prevalence, bacterial colonization prevalence and the resistance rate

### Load data
df_tau_impacts = loadRData("data/outputs_analysis1/df_tau_impacts.Rdata")%>%
  f_dynamics_plottable_tau()%>%
  left_join(df_pandI)%>%
  mutate(par_val = round(par_val, 2))

df_tau_impacts_metadata = loadRData("data/outputs_analysis1/df_tau_impacts.Rdata")%>%
  f_dynamics_metadata()%>%
  pivot_longer(-c(time, par, par_val), values_to = 'value', names_to = 'state')%>%
  left_join(df_pandI)%>%
  mutate(par_val = round(par_val, 2))%>%
  dplyr::filter(time >0) # daily metadata calculated as difference from previous day, so no data for time=0
  

### Vectors to select particular states to plot
vec_V_pa_pe = c('I_pa', 'I_pe')
vec_B_Cs_Cr = c('Cr_pa', 'Cs_pa')
vec_R_rate = c('R_rate')

### Extract min and max values
min_Rrate = min(filter(df_tau_impacts, state == 'R_rate', state %in% vec_R_rate)$value)
max_Rrate = max(filter(df_tau_impacts, state == 'R_rate', state %in% vec_R_rate)$value)

###  Pandemic impacts (in order rendered)
vec_impacts = levels(factor(df_tau_impacts$par))

vec_values_tau = as.numeric(levels(factor(df_tau_impacts$par_val)))


### (1.1.1) SARS-CoV-2 prevalence
# SARS-CoV-2 infection prevalence among patients and staff
# loop through all pandemic impacts and save each plot into an empty list

list_plots_V_prev = list()
df_tau_impacts_Vprev = df_tau_impacts%>%filter(organism == 'SARS-CoV-2', state %in% vec_V_pa_pe)
max_Vs = df_tau_impacts_Vprev%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  max_V = as.numeric(max(max_Vs[,2]))
  #max_V = as.numeric(max_Vs[i,2])
  
  dat_i = df_tau_impacts_Vprev%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_V_i = ggplot(dat_i, aes(x = time, y = value, colour = state, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    th+
    xlab('time (days)')+
    ylab('SARS-CoV-2\ninfection\nprevalence')+
    scale_colour_manual('type of individual', values = col_pat_staff, labels = c('patients', 'HCWs'))+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, guide = "none")+
    ylim(-0.001, max_V)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_V_prev[[i]] = p_V_i
}

### (1.1.2) Bacteria prevalence
# Bacterial colonizaton prevalence among patients
# loop through all pandemic impacts and save each plot into an empty list

list_plots_B_prev = list()
df_tau_impacts_Bprev = df_tau_impacts%>%filter(organism == 'ARB', state %in% vec_B_Cs_Cr)
max_Bs = df_tau_impacts_Bprev%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  max_B = as.numeric(max(max_Bs[,2]))
  #max_B = as.numeric(max_Bs[i,2])
  
  dat_i = df_tau_impacts_Bprev%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_B_i = ggplot(dat_i, aes(x = time, y = value, colour = state, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    th+
    xlab('time (days)')+
    ylab('bacterial\ncolonization\nprevalence')+
    scale_colour_manual('strain', values = col_strains, labels = c('drug-resistant', 'drug-sensitive'))+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, guide = "none")+
    ylim(0, max_B)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_B_prev[[i]] = p_B_i
}


### (1.1.3) Bacterial resistance
# Resistance rate across both strains
# loop through all pandemic impacts and save each plot into an empty list

list_plots_R_rate = list()
df_tau_impacts_R_rate = df_tau_impacts%>%filter(organism == 'ARB', state %in% vec_R_rate)

min_R_rates = df_tau_impacts_R_rate%>%group_by(par)%>%summarise(min = min(value))
max_R_rates = df_tau_impacts_R_rate%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_R_rate = as.numeric(min(min_R_rates[,2]))
  max_R_rate = as.numeric(max(max_R_rates[,2]))
  # min_R_rate = as.numeric(min_R_rates[i,2])
  # max_R_rate = as.numeric(max_R_rates[i,2])
  
  dat_i = df_tau_impacts_R_rate%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_R_i = ggplot(dat_i, aes(x = time, y = value, colour = state, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity', colour = 'black')+
    th+
    xlab('time (days)')+
    ylab('bacterial\nresistance\nrate')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    ylim(min_R_rate, max_R_rate)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(alpha=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_R_rate[[i]] = p_R_i
}


### (1.1.4) Population size 
# Npat and Nhcw through time
# loop through all pandemic impacts and save each plot into an empty list

list_plots_N = list()
df_tau_impacts_N = df_tau_impacts%>%filter(state %in% c("Total_pa", "Total_pe"))

min_Ns = df_tau_impacts_N%>%group_by(par)%>%summarise(min = min(value))
max_Ns = df_tau_impacts_N%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_N = as.numeric(min(min_Ns[,2]))
  max_N = as.numeric(max(max_Ns[,2]))
  # min_N = as.numeric(min_Ns[i,2])
  # max_N = as.numeric(max_Ns[i,2])
  
  dat_i = df_tau_impacts_N%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_N_i = ggplot(dat_i, aes(x = time, y = value, colour = state, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    th+
    xlab('time (days)')+
    ylab('number of\nindividuals\nin hospital')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    scale_colour_manual('type of individual', values = col_pat_staff, labels = c('patients', 'HCWs'))+
    ylim(min_N, max_N)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_N[[i]] = p_N_i
}


### (1.1.5) Staffing ratio
# Nhcw:Npat through time
# loop through all pandemic impacts and save each plot into an empty list

list_plots_staffing = list()
df_tau_impacts_staffing = df_tau_impacts_metadata%>%filter(state == 'staffing_ratio_daily')

min_staffings = df_tau_impacts_staffing%>%group_by(par)%>%summarise(min = min(value))
max_staffings = df_tau_impacts_staffing%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_staffing = as.numeric(min(min_staffings[,2]))
  max_staffing = as.numeric(max(max_staffings[,2]))
  # min_staffing = as.numeric(min_staffings[i,2])
  # max_staffing = as.numeric(max_staffings[i,2])
  
  dat_i = df_tau_impacts_staffing%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_staffing_i = ggplot(dat_i, aes(x = time, y = value, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity', colour = col_staffing)+
    th+
    xlab('time (days)')+
    ylab('\n\nHCW-to-patient ratio')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    ylim(min_staffing, max_staffing)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_staffing[[i]] = p_staffing_i
}


### (1.1.6) Contact rates (kappa)
# Daily contact rates for patients and staff through time
# loop through all pandemic impacts and save each plot into an empty list

list_plots_kappa = list()
df_tau_impacts_kappa = df_tau_impacts_metadata%>%filter(state %in% c("kappa_patients_daily", "kappa_staff_daily"))

min_kappas = df_tau_impacts_kappa%>%group_by(par)%>%summarise(min = min(value))
max_kappas = df_tau_impacts_kappa%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_kappa = as.numeric(min(min_kappas[,2]))
  max_kappa = as.numeric(max(max_kappas[,2]))
  # min_kappa = as.numeric(min_kappas[i,2])
  # max_kappa = as.numeric(max_kappas[i,2])
  
  dat_i = df_tau_impacts_kappa%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_kappa_i = ggplot(dat_i, aes(x = time, y = value, colour = state, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    th+
    xlab('time (days)')+
    ylab('daily number of\nclose-proximity\ninteractions')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    scale_colour_manual('type of individual', values = col_pat_staff, labels = c('patients', 'HCWs'))+
    ylim(min_kappa, max_kappa)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_kappa[[i]] = p_kappa_i
}


### (1.1.7) Hand hygiene
# Average delay between successful handwashing
# loop through all pandemic impacts and save each plot into an empty list

list_plots_hh = list()
df_tau_impacts_hh = df_tau_impacts_metadata%>%filter(state %in% c("hh_daily"))

min_hhs = df_tau_impacts_hh%>%group_by(par)%>%summarise(min = min(value))
max_hhs = df_tau_impacts_hh%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_hh = 1/as.numeric(min(min_hhs[,2]))*24
  max_hh = 1/as.numeric(max(max_hhs[,2]))*24
  # min_hh = 1/as.numeric(min_hhs[i,2])*24
  # max_hh = 1/as.numeric(max_hhs[i,2])*24
  
  dat_i = df_tau_impacts_hh%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_hh_i = ggplot(dat_i, aes(x = time, y = (1/value)*24, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity', colour = col_hh)+
    th+
    xlab('time (days)')+
    ylab('delay between\ncompliant hand-washing\nevents (hours)')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    ylim(min_hh, max_hh)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_hh[[i]] = p_hh_i
}

### (1.1.8) Antibiotics
# Daily number of patients on antibiotics
# loop through all pandemic impacts and save each plot into an empty list

list_plots_abx = list()
df_tau_impacts_abx = df_tau_impacts_metadata%>%filter(state %in% c("abx_daily"))

min_abxs = df_tau_impacts_abx%>%group_by(par)%>%summarise(min = min(value))
max_abxs = df_tau_impacts_abx%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_abx = as.numeric(min(min_abxs[,2]))
  max_abx = as.numeric(max(max_abxs[,2]))
  # min_abx = as.numeric(min_abxs[i,2])
  # max_abx = as.numeric(max_abxs[i,2])
  
  dat_i = df_tau_impacts_abx%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_abx_i = ggplot(dat_i, aes(x = time, y = value, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity', colour = col_abx)+
    th+
    xlab('time (days)')+
    ylab('number of patients\nexposed to\nantibiotics')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    ylim(min_abx, max_abx)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_abx[[i]] = p_abx_i
}



### (1.1.9) Admission
# Daily number of patients admitted colonized with resistant bacteria
# loop through all pandemic impacts and save each plot into an empty list

list_plots_admission = list()
df_tau_impacts_admission = df_tau_impacts_metadata%>%filter(state %in% c("adm_R_daily"))

min_admissions = df_tau_impacts_admission%>%group_by(par)%>%summarise(min = min(value))
max_admissions = df_tau_impacts_admission%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_admission = as.numeric(min(min_admissions[,2]))
  max_admission = as.numeric(max(max_admissions[,2]))
  # min_admission = as.numeric(min_admissions[i,2])
  # max_admission = as.numeric(max_admissions[i,2])
  
  dat_i = df_tau_impacts_admission%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_admission_i = ggplot(dat_i, aes(x = time, y = value, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity', colour = col_admission)+
    th+
    xlab('time (days)')+
    ylab('number of patients\nadmitted colonized with\nresistant strain')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    ylim(min_admission, max_admission)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_admission[[i]] = p_admission_i
}

### (1.1.10) SARS-CoV-2 force of infection
# SARS-CoV-2 force of infection FROM patients AND staff (i.e. not incidence)
# loop through all pandemic impacts and save each plot into an empty list

list_plots_foi_V = list()
df_tau_impacts_foi_V = df_tau_impacts_metadata%>%filter(state %in% c("foi_V_daily"))

min_foi_Vs = df_tau_impacts_foi_V%>%group_by(par)%>%summarise(min = min(value))
max_foi_Vs = df_tau_impacts_foi_V%>%group_by(par)%>%summarise(max = max(value))

for(i in 1:length(vec_impacts)){
  
  par_i = vec_impacts[i]
  min_foi_V = as.numeric(min(min_foi_Vs[,2]))
  max_foi_V = as.numeric(max(max_foi_Vs[,2]))
  # min_foi_V = as.numeric(min_foi_Vs[i,2])
  # max_foi_V = as.numeric(max_foi_Vs[i,2])
  
  dat_i = df_tau_impacts_foi_V%>%filter(par == par_i)
  if(levels(factor(dat_i$type)) == 'burden'){alpha_vline = 0}else{alpha_vline = 0.3}
  
  p_foi_V_i = ggplot(dat_i, aes(x = time, y = value, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity', colour = col_foi)+
    th+
    xlab('time (days)')+
    ylab('\nSARS-CoV-2\nforce of infection')+
    scale_alpha_continuous(expression(tau), range = alpha_range_tau, breaks = vec_values_tau)+
    ylim(min_foi_V, max_foi_V)+
    geom_vline(xintercept = pars_m['t_policy'], alpha = alpha_vline, size = 2)+
    guides(colour=guide_legend(nrow=2,byrow=TRUE))
  
  list_plots_foi_V[[i]] = p_foi_V_i
}


#####################
### COMBINE PLOTS ###
#####################

p_impacts_individual = plot_grid(list_plots_V_prev[[which(vec_impacts == "a_surge")]]+
                                   ggtitle(expression(paste(bold("a.  "), "abandoned stewardship (", tau[as], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "a_surge")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "a_surge")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "prophylaxis")]]+
                                   ggtitle(expression(paste(bold("b.  "), "COVID-19 prophylaxis (", tau[cp], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "prophylaxis")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "prophylaxis")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "disorg")]]+
                                   ggtitle(expression(paste(bold("c.  "), "cohorting disorganization (", tau[cd], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "disorg")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "disorg")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "distancing")]]+
                                   ggtitle(expression(paste(bold("d.  "), "patient lockdown (", tau[pl], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "distancing")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "distancing")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "masks")]]+
                                   ggtitle(expression(paste(bold("e.  "), "universal masking (", tau[um], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "masks")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "masks")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "handrub")]]+
                                   ggtitle(expression(paste(bold("f.  "), "hand hygiene (", tau[hh], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "handrub")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "handrub")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "covid_stay")]]+
                                   ggtitle(expression(paste(bold("g.  "), "COVID-19 stays (", tau[cs], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "covid_stay")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "covid_stay")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "sickleave")]]+
                                   ggtitle(expression(paste(bold("h.  "), "staff sick leave (", tau[ss], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "sickleave")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "sickleave")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "adm_reduc")]]+
                                   ggtitle(expression(paste(bold("i.  "), "reduced admission (", tau[ra], ")")))+
                                   theme(legend.position = "none"),
                                 list_plots_B_prev[[which(vec_impacts == "adm_reduc")]]+
                                   theme(legend.position = "none"),
                                 list_plots_R_rate[[which(vec_impacts == "adm_reduc")]]+
                                   theme(legend.position = "none"),
                                 list_plots_V_prev[[which(vec_impacts == "comm_denom")]]+
                                   ggtitle(expression(paste(bold("j.  "), "sicker casemix (", tau[sc], ")")))+
                                   theme(legend.position = "bottom"),
                                 list_plots_B_prev[[which(vec_impacts == "comm_denom")]]+
                                   theme(legend.position = "bottom"),
                                 list_plots_R_rate[[which(vec_impacts == "comm_denom")]]+
                                   theme(legend.position = "bottom"),
                                 ncol = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5))


ggsave("plots/tau_impacts_individual.png", plot = p_impacts_individual, width = 30, height = 50, unit = 'cm', )
ggsave("plots/tau_impacts_individual.pdf", plot = p_impacts_individual, width = 30, height = 50, unit = 'cm', )

### Main categories

p_impacts_categories_V = ggarrange(list_plots_V_prev[[which(vec_impacts == "combo_ipc")]], 
                                   list_plots_V_prev[[which(vec_impacts == "combo_antibiotics")]], 
                                   list_plots_V_prev[[which(vec_impacts == "combo_contact")]], 
                                   list_plots_V_prev[[which(vec_impacts == "combo_disease")]], 
                                   list_plots_V_prev[[which(vec_impacts == "combo_admission")]], 
                                   nrow = 5, ncol = 1,
                                   common.legend = T,
                                   legend = 'bottom')


p_impacts_categories_B = ggarrange(list_plots_B_prev[[which(vec_impacts == "combo_ipc")]], 
                                   list_plots_B_prev[[which(vec_impacts == "combo_antibiotics")]], 
                                   list_plots_B_prev[[which(vec_impacts == "combo_contact")]],
                                   list_plots_B_prev[[which(vec_impacts == "combo_disease")]], 
                                   list_plots_B_prev[[which(vec_impacts == "combo_admission")]], 
                                   nrow = 5, ncol = 1,
                                   common.legend = T,
                                   legend = 'bottom')

p_impacts_categories_R_rate = ggarrange(list_plots_R_rate[[which(vec_impacts == "combo_ipc")]], 
                                        list_plots_R_rate[[which(vec_impacts == "combo_antibiotics")]], 
                                        list_plots_R_rate[[which(vec_impacts == "combo_contact")]],
                                        list_plots_R_rate[[which(vec_impacts == "combo_disease")]], 
                                        list_plots_R_rate[[which(vec_impacts == "combo_admission")]],
                                        nrow = 5, ncol = 1,
                                        common.legend = T,
                                        legend = 'bottom')

p_impacts_categories = ggarrange(p_impacts_categories_V,
                                 p_impacts_categories_B,
                                 p_impacts_categories_R_rate,
                                 ncol = 3)
ggsave("plots/tau_impacts_categories.png", plot = p_impacts_categories, width = 25, height = 30, unit = 'cm', )
ggsave("plots/tau_impacts_categories.pdf", plot = p_impacts_categories, width = 25, height = 30, unit = 'cm', )


### IPC

p_impacts_IPC_V = ggarrange(list_plots_V_prev[[which(vec_impacts == "masks")]]+ggtitle(expression(paste(bold("a.  "), "universal masking (", tau[um],")"))), 
                            list_plots_V_prev[[which(vec_impacts == "handrub")]]+ggtitle(expression(paste(bold("d.  "), "hand hygiene (", tau[hh],")"))),
                            list_plots_V_prev[[which(vec_impacts == "combo_ipc")]]+ggtitle(expression(paste(bold("g.  "), "both (", tau[um]," and ", tau[hh],")"))),
                            nrow = 3, ncol = 1,
                            common.legend = T,
                            legend = 'bottom',
                            labels = c('A', 'B', 'C'))

p_impacts_IPC_B = ggarrange(list_plots_B_prev[[which(vec_impacts == "masks")]]+ggtitle(expression(paste(bold("b.  "), "universal masking (", tau[um],")"))), 
                            list_plots_B_prev[[which(vec_impacts == "handrub")]]+ggtitle(expression(paste(bold("e.  "), "hand hygiene (", tau[hh],")"))),
                            list_plots_B_prev[[which(vec_impacts == "combo_ipc")]]+ggtitle(expression(paste(bold("h.  "), "both (", tau[um]," and ", tau[hh],")"))),
                            nrow = 3, ncol = 1,
                            common.legend = T,
                            legend = 'bottom')

p_impacts_IPC_R_rate = ggarrange(list_plots_R_rate[[which(vec_impacts == "masks")]]+ggtitle(expression(paste(bold("c.  "), "universal masking (", tau[um],")"))), 
                            list_plots_R_rate[[which(vec_impacts == "handrub")]]+ggtitle(expression(paste(bold("f.  "), "hand hygiene (", tau[hh],")"))),
                            list_plots_R_rate[[which(vec_impacts == "combo_ipc")]]+ggtitle(expression(paste(bold("i.  "), "both (", tau[um]," and ", tau[hh],")"))),
                            nrow = 3, ncol = 1,
                            common.legend = T,
                            legend = 'bottom')

p_impacts_IPC = ggarrange(p_impacts_IPC_V,
                          p_impacts_IPC_B,
                          p_impacts_IPC_R_rate,
                          ncol = 3)

ggsave("plots/tau_impacts_IPC.png", plot = p_impacts_IPC, width = 25, height = 20, unit = 'cm', )
ggsave("plots/tau_impacts_IPC.pdf", plot = p_impacts_IPC, width = 25, height = 20, unit = 'cm', )


### ABX

p_impacts_abx_V = ggarrange(list_plots_V_prev[[which(vec_impacts == "prophylaxis")]]+ggtitle(expression(paste(bold("a.  "), "COVID-19 prophylaxis (", tau[cp],")"))), 
                            list_plots_V_prev[[which(vec_impacts == "a_surge")]]+ggtitle(expression(paste(bold("d.  "), "abandoned stewardship (", tau[as],")"))),
                            list_plots_V_prev[[which(vec_impacts == "combo_antibiotics")]]+ggtitle(expression(paste(bold("g.  "), "both (", tau[cp]," and ", tau[as], ")"))),
                            nrow = 3, ncol = 1,
                            common.legend = T,
                            legend = 'bottom')

p_impacts_abx_B = ggarrange(list_plots_B_prev[[which(vec_impacts == "prophylaxis")]]+ggtitle(expression(paste(bold("b.  "), "COVID-19 prophylaxis (", tau[cp],")"))), 
                            list_plots_B_prev[[which(vec_impacts == "a_surge")]]+ggtitle(expression(paste(bold("e.  "), "abandoned stewardship (", tau[as],")"))),
                            list_plots_B_prev[[which(vec_impacts == "combo_antibiotics")]]+ggtitle(expression(paste(bold("h.  "), "both (", tau[cp]," and ", tau[as], ")"))),
                            nrow = 3, ncol = 1,
                            common.legend = T,
                            legend = 'bottom')

p_impacts_abx_R_rate = ggarrange(list_plots_R_rate[[which(vec_impacts == "prophylaxis")]]+ggtitle(expression(paste(bold("c.  "), "COVID-19 prophylaxis (", tau[cp],")"))), 
                                 list_plots_R_rate[[which(vec_impacts == "a_surge")]]+ggtitle(expression(paste(bold("f.  "), "abandoned stewardship (", tau[as],")"))),
                                 list_plots_R_rate[[which(vec_impacts == "combo_antibiotics")]]+ggtitle(expression(paste(bold("i.  "), "both (", tau[cp]," and ", tau[as], ")"))),
                                 nrow = 3, ncol = 1,
                                 common.legend = T,
                                 legend = 'bottom')

p_impacts_abx = ggarrange(p_impacts_abx_V,
                          p_impacts_abx_B,
                          p_impacts_abx_R_rate,
                          ncol = 3)

ggsave("plots/tau_impacts_abx.png", plot = p_impacts_abx, width = 27, height = 20, unit = 'cm', )
ggsave("plots/tau_impacts_abx.pdf", plot = p_impacts_abx, width = 27, height = 20, unit = 'cm', )


### Contact

p_impacts_contact_V = ggarrange(list_plots_V_prev[[which(vec_impacts == "disorg")]]+ggtitle(expression(paste(bold("a.  "), "cohort disorganization (", tau[cd],")"))), 
                                list_plots_V_prev[[which(vec_impacts == "distancing")]]+ggtitle(expression(paste(bold("d.  "), "patient lockdown (", tau[pl],")"))),
                                list_plots_V_prev[[which(vec_impacts == "combo_contact")]]+ggtitle(expression(paste(bold("g.  "), "both (", tau[cd]," and ", tau[pl],")"))),
                                nrow = 3, ncol = 1,
                                common.legend = T,
                                legend = 'bottom')

p_impacts_contact_B = ggarrange(list_plots_B_prev[[which(vec_impacts == "disorg")]]+ggtitle(expression(paste(bold("b.  "), "cohort disorganization (", tau[cd],")"))), 
                                list_plots_B_prev[[which(vec_impacts == "distancing")]]+ggtitle(expression(paste(bold("e.  "), "patient lockdown (", tau[pl],")"))),
                                list_plots_B_prev[[which(vec_impacts == "combo_contact")]]+ggtitle(expression(paste(bold("h.  "), "both (", tau[cd]," and ", tau[pl],")"))),
                                nrow = 3, ncol = 1,
                                common.legend = T,
                                legend = 'bottom')

p_impacts_contact_R_rate = ggarrange(list_plots_R_rate[[which(vec_impacts == "disorg")]]+ggtitle(expression(paste(bold("c.  "), "cohort disorganization (", tau[cd],")"))), 
                                     list_plots_R_rate[[which(vec_impacts == "distancing")]]+ggtitle(expression(paste(bold("f.  "), "patient lockdown (", tau[pl],")"))),
                                     list_plots_R_rate[[which(vec_impacts == "combo_contact")]]+ggtitle(expression(paste(bold("i.  "), "both (", tau[cd]," and ", tau[pl],")"))),
                                     nrow = 3, ncol = 1,
                                     common.legend = T,
                                     legend = 'bottom')

p_impacts_contact = ggarrange(p_impacts_contact_V,
                              p_impacts_contact_B,
                              p_impacts_contact_R_rate,
                              ncol = 3)

ggsave("plots/tau_impacts_contact.png", plot = p_impacts_contact, width = 25, height = 20, unit = 'cm', )
ggsave("plots/tau_impacts_contact.pdf", plot = p_impacts_contact, width = 25, height = 20, unit = 'cm', )


### Disease

p_impacts_disease_V = ggarrange(list_plots_V_prev[[which(vec_impacts == "covid_stay")]]+ggtitle(expression(paste(bold("a.  "), "COVID-19 stays (", tau[cs],")"))), 
                                list_plots_V_prev[[which(vec_impacts == "sickleave")]]+ggtitle(expression(paste(bold("d.  "), "staff sick leave (", tau[ss],")"))),
                                list_plots_V_prev[[which(vec_impacts == "combo_disease")]]+ggtitle(expression(paste(bold("g.  "), "both (", tau[cc]," and ", tau[ss],")"))),
                                nrow = 3, ncol = 1,
                                common.legend = T,
                                legend = 'bottom')

p_impacts_disease_B = ggarrange(list_plots_B_prev[[which(vec_impacts == "covid_stay")]]+ggtitle(expression(paste(bold("b.  "), "COVID-19 stays (", tau[cs],")"))), 
                                list_plots_B_prev[[which(vec_impacts == "sickleave")]]+ggtitle(expression(paste(bold("e.  "), "staff sick leave (", tau[ss],")"))),
                                list_plots_B_prev[[which(vec_impacts == "combo_disease")]]+ggtitle(expression(paste(bold("h.  "), "both (", tau[cc]," and ", tau[ss],")"))),
                                nrow = 3, ncol = 1,
                                common.legend = T,
                                legend = 'bottom')

p_impacts_disease_R_rate = ggarrange(list_plots_R_rate[[which(vec_impacts == "covid_stay")]]+ggtitle(expression(paste(bold("c.  "), "COVID-19 stays (", tau[cs],")"))), 
                                     list_plots_R_rate[[which(vec_impacts == "sickleave")]]+ggtitle(expression(paste(bold("f.  "), "staff sick leave (", tau[ss],")"))),
                                     list_plots_R_rate[[which(vec_impacts == "combo_disease")]]+ggtitle(expression(paste(bold("i.  "), "both (", tau[cc]," and ", tau[ss],")"))),
                                     nrow = 3, ncol = 1,
                                     common.legend = T,
                                     legend = 'bottom')

p_impacts_disease = ggarrange(p_impacts_disease_V,
                              p_impacts_disease_B,
                              p_impacts_disease_R_rate,
                              ncol = 3)

ggsave("plots/tau_impacts_disease.png", plot = p_impacts_disease, width = 25, height = 20, unit = 'cm', )
ggsave("plots/tau_impacts_disease.pdf", plot = p_impacts_disease, width = 25, height = 20, unit = 'cm', )


### Admission

p_impacts_admission_V = ggarrange(list_plots_V_prev[[which(vec_impacts == "adm_reduc")]]+ggtitle(expression(paste(bold("a.  "), "reduced admission (", tau[ra],")"))), 
                                  list_plots_V_prev[[which(vec_impacts == "comm_denom")]]+ggtitle(expression(paste(bold("d.  "), "sicker casemix (", tau[sc],")"))),
                                  list_plots_V_prev[[which(vec_impacts == "combo_admission")]]+ggtitle(expression(paste(bold("g.  "), "both (", tau[ra]," and ", tau[rc],")"))),
                                  nrow = 3, ncol = 1,
                                  common.legend = T,
                                  legend = 'bottom')

p_impacts_admission_B = ggarrange(list_plots_B_prev[[which(vec_impacts == "adm_reduc")]]+ggtitle(expression(paste(bold("b.  "), "reduced admission (", tau[ra],")"))), 
                                  list_plots_B_prev[[which(vec_impacts == "comm_denom")]]+ggtitle(expression(paste(bold("e.  "), "sicker casemix (", tau[sc],")"))),
                                  list_plots_B_prev[[which(vec_impacts == "combo_admission")]]+ggtitle(expression(paste(bold("h.  "), "both (", tau[ra]," and ", tau[rc],")"))),
                                  nrow = 3, ncol = 1,
                                  common.legend = T,
                                  legend = 'bottom')

p_impacts_admission_R_rate = ggarrange(list_plots_R_rate[[which(vec_impacts == "adm_reduc")]]+ggtitle(expression(paste(bold("c.  "), "reduced admission (", tau[ra],")"))), 
                                       list_plots_R_rate[[which(vec_impacts == "comm_denom")]]+ggtitle(expression(paste(bold("f.  "), "sicker casemix (", tau[sc],")"))),
                                       list_plots_R_rate[[which(vec_impacts == "combo_admission")]]+ggtitle(expression(paste(bold("i.  "), "both (", tau[ra]," and ", tau[rc],")"))),
                                       nrow = 3, ncol = 1,
                                       common.legend = T,
                                       legend = 'bottom')

p_impacts_admission = ggarrange(p_impacts_admission_V,
                                p_impacts_admission_B,
                                p_impacts_admission_R_rate,
                                ncol = 3)

ggsave("plots/tau_impacts_admission.png", plot = p_impacts_admission, width = 25, height = 20, unit = 'cm', )
ggsave("plots/tau_impacts_admission.pdf", plot = p_impacts_admission, width = 25, height = 20, unit = 'cm', )



### All impacts

p_impacts_all = ggarrange(list_plots_V_prev[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("a."))))+
                            guides(colour=guide_legend(nrow=2,byrow=TRUE), alpha = FALSE)+theme(legend.position = c(0.75, 0.92))+ylab('SARS-CoV-2 infection prevalence'),
                          list_plots_B_prev[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("b."))))+
                            guides(colour=guide_legend(nrow=2,byrow=TRUE, alpha = FALSE))+theme(legend.position = c(0.75, 0.92))+ylab('bacterial colonization prevalence'),
                          list_plots_R_rate[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("c."))))+
                            guides(alpha=guide_legend(nrow=2,byrow=TRUE))+theme(legend.position = c(0.7, 0.92))+ylab('bacterial resistance rate'),
                          ncol = 3, nrow = 1)
ggsave("plots/tau_impacts_all.png", plot = p_impacts_all, width = 30, height = 10, unit = 'cm', )
ggsave("plots/tau_impacts_all.pdf", plot = p_impacts_all, width = 30, height = 10, unit = 'cm', )


### All impacts, including metadata

p_impacts_all_with_metadata = ggarrange(list_plots_N[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("a.")))),
                                        list_plots_staffing[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("b.")))),
                                        list_plots_hh[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("c.")))),
                                        list_plots_kappa[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("d.")))),
                                        list_plots_abx[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("e.")))),
                                        list_plots_admission[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("f.")))),
                                        list_plots_V_prev[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("g.")))),
                                        list_plots_B_prev[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("h.")))),
                                        list_plots_R_rate[[which(vec_impacts == "combo_all")]]+ggtitle(expression(paste(bold("i.")))),
                                        ncol = 3, nrow = 3, 
                                        legend = 'none')
ggsave("plots/tau_impacts_all_with_metadata.png", plot = p_impacts_all_with_metadata, width = 25, height = 20, unit = 'cm', )
ggsave("plots/tau_impacts_all_with_metadata.pdf", plot = p_impacts_all_with_metadata, width = 25, height = 20, unit = 'cm', )


### Relevant metadata for each impact
# tweak y-axes and ggtitles plot-by-plot
p_impacts_each_on_metadata = ggarrange(list_plots_foi_V[[which(vec_impacts == "masks")]]+ggtitle(expression(paste(bold("a.  "), "universal masking (", tau[um],")")))+
                                         guides(alpha=guide_legend(nrow=5,byrow=TRUE, override.aes = list(colour = 'black')))+theme(legend.position = c(0.8, 0.7)),
                                         list_plots_hh[[which(vec_impacts == "handrub")]]+ggtitle(expression(paste(bold("b.  "), "hand hygiene (", tau[hh],")")))+ theme(legend.position="none"),
                                         list_plots_abx[[which(vec_impacts == "prophylaxis")]]+ggtitle(expression(paste(bold("c.  "), "COVID-19 prophylaxis (", tau[cp],")")))+ theme(legend.position="none"),
                                         list_plots_abx[[which(vec_impacts == "a_surge")]]+ggtitle(expression(paste(bold("d.  "), "abandoned stewardship (", tau[as],")")))+ theme(legend.position="none"),
                                         list_plots_kappa[[which(vec_impacts == "distancing")]]+ggtitle(expression(paste(bold("e.  "), "patient lockdown (", tau[pl],")")))+ theme(legend.position="none"),
                                         list_plots_kappa[[which(vec_impacts == "disorg")]]+ggtitle(expression(paste(bold("f.  "), "cohorting disorganization (", tau[cd],")")))+ theme(legend.position="none")+
                                         guides(colour=guide_legend(nrow=1,byrow=TRUE), alpha = FALSE)+theme(legend.position = c(0.75, 0.15)),
                                         list_plots_staffing[[which(vec_impacts == "covid_stay")]]+ggtitle(expression(paste(bold("g.  "), "COVID-19 stays (", tau[cs],")")))+ theme(legend.position="none"),
                                         list_plots_staffing[[which(vec_impacts == "sickleave")]]+ggtitle(expression(paste(bold("h.  "), "staff sick leave (", tau[ss],")")))+ theme(legend.position="none"),
                                         list_plots_admission[[which(vec_impacts == "adm_reduc")]]+ggtitle(expression(paste(bold("i.  "), "reduced admission (", tau[ra],")")))+ theme(legend.position="none"),
                                         list_plots_admission[[which(vec_impacts == "comm_denom")]]+ggtitle(expression(paste(bold("j.  "), "sicker casemix (", tau[sc],")")))+ theme(legend.position="none"),
                                         ncol = 2, nrow = 5, align = 'hv')

ggsave("plots/tau_impacts_each_on_metadata.png", plot = p_impacts_each_on_metadata, width = 25, height = 40, unit = 'cm', )
ggsave("plots/tau_impacts_each_on_metadata.pdf", plot = p_impacts_each_on_metadata, width = 25, height = 40, unit = 'cm', )



################################
### Dynamics over random tau ###
################################

df_tau_random = loadRData("data/outputs_analysis1/df_tau_random.Rdata")%>%
  mutate(I_pa = I_U_pa + I_Cs_pa + I_Cr_pa,
         I_pe = I_U_pe + I_Cs_pe + I_Cr_pe,
         Cs_pa = S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa,
         Cr_pa = S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa)%>%
  mutate(R_rate = Cr_pa/(Cs_pa + Cr_pa))%>%
  dplyr::select(n_sim, time, I_pa, I_pe, Cs_pa, Cr_pa, R_rate)%>%
  pivot_longer(-c(n_sim, time), names_to = "state", values_to = "value")

df_tau_random0 = df_tau_random%>%
  filter(time == 0)

df_tau_random1 = df_tau_random0%>%mutate(time = -1)
df_tau_random2 = df_tau_random0%>%mutate(time = -2)
df_tau_random3 = df_tau_random0%>%mutate(time = -3)
df_tau_random4 = df_tau_random0%>%mutate(time = -4)
df_tau_random5 = df_tau_random0%>%mutate(time = -5)
df_tau_random6 = df_tau_random0%>%mutate(time = -6)
df_tau_random7 = df_tau_random0%>%mutate(time = -7)
df_tau_random8 = df_tau_random0%>%mutate(time = -8)
df_tau_random9 = df_tau_random0%>%mutate(time = -9)
df_tau_random10 = df_tau_random0%>%mutate(time = -10)

df_tau_random = rbind(df_tau_random, df_tau_random1, df_tau_random2, df_tau_random3, df_tau_random4, df_tau_random5,
                      df_tau_random6, df_tau_random7, df_tau_random8, df_tau_random9, df_tau_random10)

df_tau_random_mean = df_tau_random%>%group_by(time, state)%>%summarize(mean = mean(value))

p_impacts_random_V = ggplot(df_tau_random%>%filter(state %in% c('I_pa', 'I_pe')), aes(x = time, y = abs(value), colour = factor(state), group = interaction(state, factor(n_sim))))+
  geom_line(stat = 'identity', alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('SARS-CoV-2 infection prevalence')+
  ylim(0, 120)+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  theme(legend.position = c(0.8,0.95))+
  scale_colour_manual('type of individual', values = col_pat_staff, labels = c('patients', 'HCWs'))+
  scale_fill_manual('type of individual', values = col_pat_staff, labels = c('patients', 'HCWs'))+
  geom_point(df_tau_random_mean%>%filter(state %in% c('I_pa', 'I_pe'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state, fill = factor(state)), 
             alpha = 1, colour = 'black', shape = 21, size = 2)+
  guides(colour=guide_legend(override.aes = list(alpha = 1)))

p_impacts_random_B = ggplot(df_tau_random%>%filter(state %in% c('Cs_pa', 'Cr_pa')), aes(x = time, y = abs(value), colour = factor(state), group = interaction(state, factor(n_sim))))+
  geom_line(stat = 'identity', alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('bacterial colonization prevalence')+
  ylim(0, 120)+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  theme(legend.position = c(0.8,0.95))+
  scale_colour_manual('strain', values = col_strains, labels = c('drug-resistant', 'drug-sensitive'))+
  scale_fill_manual('strain', values = col_strains, labels = c('drug-resistant', 'drug-sensitive'))+
  geom_point(df_tau_random_mean%>%filter(state %in% c('Cs_pa', 'Cr_pa'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state, fill = factor(state)), 
             alpha = 1, colour = 'black', shape = 21, size = 2)+
  guides(colour=guide_legend(override.aes = list(alpha = 1)))

p_impacts_random_R_rate = ggplot(df_tau_random%>%filter(state %in% c('R_rate')), aes(x = time, y = abs(value), group = interaction(state, factor(n_sim))))+
  geom_line(stat = 'identity', alpha = 0.2, colour = 'black')+
  th+
  xlab('time (days)')+
  ylab('bacterial resistance rate')+
  ylim(0.55, 1)+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  geom_point(df_tau_random_mean%>%filter(state %in% c('R_rate'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state),
             fill = "white", alpha = 1, colour = 'black', shape = 21, size = 2)

p_impacts_random = ggarrange(p_impacts_random_V+ggtitle(expression(paste(bold("a.  ")))), 
                             p_impacts_random_B+ggtitle(expression(paste(bold("b.  ")))), 
                             p_impacts_random_R_rate+ggtitle(expression(paste(bold("c.  ")))), 
                             ncol = 3)

ggsave("plots/tau_impacts_random.png", plot = p_impacts_random, width = 30, height = 10, unit = 'cm', )
ggsave("plots/tau_impacts_random.pdf", plot = p_impacts_random, width = 30, height = 10, unit = 'cm', )

################################
### METADATA FROM RANDOM TAU ###
################################

### Load and prep data
df_tau_random_metadata = loadRData("data/outputs_analysis1/df_tau_random.Rdata")%>%
  f_dynamics_metadata_n_sim(par = 'all', par_val = 'random')%>%
  pivot_longer(-c(n_sim, time, par, par_val), values_to = 'value', names_to = 'state')%>%
  filter(time>0)

df_tau_random_metadata0a = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "abx_daily", value = as.numeric(pars_m['a']*pars_m['Nbeds']))
df_tau_random_metadata0b = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "adm_daily", value = as.numeric(pars_m['mu']*pars_m['Nbeds']))
df_tau_random_metadata0c = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "adm_R_daily", value = as.numeric(pars_m['mu']*pars_m['Nbeds']*pars_m['f_Cr']))
df_tau_random_metadata0d = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "foi_Cr_daily", value = 0.10622)
df_tau_random_metadata0e = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "foi_V_daily", value = 0)
df_tau_random_metadata0f = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "hh_daily", value = as.numeric(f_eta(pars_m['hyg'], pars_m['kappa_pe_pa'], pars_m['Nhcw']/pars_m['Nbeds'])))
df_tau_random_metadata0g = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "kappa_patients_daily", value = as.numeric(pars_m['kappa_pa_pe']+pars_m['kappa_pa_pa']))
df_tau_random_metadata0h = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "kappa_staff_daily", value = as.numeric(pars_m['kappa_pe_pa']+pars_m['kappa_pe_pe']))
df_tau_random_metadata0i = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "patientdays_daily", value = as.numeric(pars_m['Nbeds']))
df_tau_random_metadata0j = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "staffdays_daily", value = as.numeric(pars_m['Nhcw']))
df_tau_random_metadata0k = data.frame(n_sim = 1:100, time = 0, par = "all", par_val = "random", state = "staffing_ratio_daily", value = as.numeric(pars_m['Nhcw']/pars_m['Nbeds']))

df_tau_random_metadata0 = rbind(df_tau_random_metadata0a, df_tau_random_metadata0b, df_tau_random_metadata0c, df_tau_random_metadata0d, df_tau_random_metadata0e, df_tau_random_metadata0f,
                                df_tau_random_metadata0g, df_tau_random_metadata0h, df_tau_random_metadata0i, df_tau_random_metadata0j, df_tau_random_metadata0k)

df_tau_random_metadata1 = df_tau_random_metadata0%>%mutate(time = -1)
df_tau_random_metadata2 = df_tau_random_metadata0%>%mutate(time = -2)
df_tau_random_metadata3 = df_tau_random_metadata0%>%mutate(time = -3)
df_tau_random_metadata4 = df_tau_random_metadata0%>%mutate(time = -4)
df_tau_random_metadata5 = df_tau_random_metadata0%>%mutate(time = -5)
df_tau_random_metadata6 = df_tau_random_metadata0%>%mutate(time = -6)
df_tau_random_metadata7 = df_tau_random_metadata0%>%mutate(time = -7)
df_tau_random_metadata8 = df_tau_random_metadata0%>%mutate(time = -8)
df_tau_random_metadata9 = df_tau_random_metadata0%>%mutate(time = -9)
df_tau_random_metadata10 = df_tau_random_metadata0%>%mutate(time = -10)

df_tau_random_metadata = rbind(df_tau_random_metadata, df_tau_random_metadata0, df_tau_random_metadata1, df_tau_random_metadata2, df_tau_random_metadata3, df_tau_random_metadata4, 
                               df_tau_random_metadata5,df_tau_random_metadata6, df_tau_random_metadata7, df_tau_random_metadata8, df_tau_random_metadata9, df_tau_random_metadata10)

df_tau_random_metadata_mean = df_tau_random_metadata%>%
  group_by(time, par, par_val, state)%>%
  summarise(mean = mean(value))


### kappa
p_impacts_random_kappa = df_tau_random_metadata%>%
  filter(state %in% c('kappa_patients_daily', 'kappa_staff_daily'))%>%
  ggplot(aes(x = time, y = value, group = interaction(n_sim,state), colour = state))+
  geom_line(stat = 'identity', alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('daily number of\nclose-proximity\ncontacts')+
  scale_colour_manual('type of individual', 
                      values = col_pat_staff, 
                      labels = c("patients", "HCWs"))+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.position = c(0.7,0.85))+
  ylim(15,32)+
  geom_point(df_tau_random_metadata_mean%>%filter(state %in% c('kappa_patients_daily', 'kappa_staff_daily'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state, fill = state),
             alpha = 1, colour = 'black', shape = 21, size = 2)+
  scale_fill_manual('type of individual', 
                      values = col_pat_staff, 
                      labels = c("patients", "HCWs"))

### hand hygiene
p_impacts_random_hh = df_tau_random_metadata%>%
  filter(state == 'hh_daily')%>%
  ggplot(aes(x = time, y = 1/value*24, group = interaction(n_sim,state)))+
  geom_line(stat = 'identity', alpha = 0.2, colour = col_hh)+
  th+
  xlab('time (days)')+
  ylab('delay between\ncompliant hand-washing\nevents (hours)')+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  geom_point(df_tau_random_metadata_mean%>%filter(state %in% c('hh_daily'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = 1/abs(mean)*24, group = state),
             fill = brightness(col_hh, 1.5),alpha = 1, colour = 'black', shape = 21, size = 2)

### foi
p_impacts_random_foi = df_tau_random_metadata%>%
  filter(state %in% c('foi_Cr_daily', 'foi_V_daily'))%>%
  ggplot(aes(x = time, y = value, group = interaction(n_sim,state), colour = state))+
  geom_line(stat = 'identity', alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('\n\nforce of infection')+
  scale_colour_manual('pathogen', 
                      values = col_pathogen, 
                      labels = c("resistant bacteria", "SARS-CoV-2"))+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(alpha = 1)))+
  theme(legend.position = c(0.7,0.85))+
  geom_point(df_tau_random_metadata_mean%>%filter(state %in% c('foi_Cr_daily', 'foi_V_daily'),
                                                  time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state, fill = state),
             alpha = 1, colour = 'black', shape = 21, size = 2)+
  scale_fill_manual('pathogen', 
                      values = brightness(col_pathogen,1.5), 
                      labels = c("resistant bacteria", "SARS-CoV-2"))

### staffing ratio
p_impacts_random_staffing = df_tau_random_metadata%>%
  filter(state == 'staffing_ratio_daily')%>%
  ggplot(aes(x = time, y = value, group = n_sim))+
  geom_line(stat = 'identity', colour = col_staffing, alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('\n\nHCW-to-patient ratio')+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  geom_point(df_tau_random_metadata_mean%>%filter(state %in% c('staffing_ratio_daily'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state),
             fill = brightness(col_staffing, 1.5),alpha = 1, colour = 'black', shape = 21, size = 2)

### Abx 
p_impacts_random_abx = df_tau_random_metadata%>%
  filter(state == 'abx_daily')%>%
  ggplot(aes(x = time, y = value, group = n_sim))+
  geom_line(stat = 'identity', colour = col_abx, alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('number of patients\nexposed to\nantibiotics')+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  ylim(0,185)+
  geom_point(df_tau_random_metadata_mean%>%filter(state %in% c('abx_daily'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state),
             fill = col_abx,alpha = 1, colour = 'black', shape = 21, size = 2)

### Adm_R
p_impacts_random_admR = df_tau_random_metadata%>%
  filter(state == 'adm_R_daily')%>%
  ggplot(aes(x = time, y = value, group = n_sim))+
  geom_line(stat = 'identity', colour = col_admission, alpha = 0.2)+
  th+
  xlab('time (days)')+
  ylab('number of patients\nadmitted colonized with\nresistant strain')+
  geom_vline(xintercept = pars_m['t_policy'], alpha = 0.3, size = 2)+
  geom_vline(xintercept = 0, alpha = 0.8, size = 0.5, linetype = 2)+
  geom_point(df_tau_random_metadata_mean%>%filter(state %in% c('adm_R_daily'), time %in% seq(-10,180, by = 5)), 
             mapping = aes(x = time, y = abs(mean), group = state),
             fill = brightness(col_admission, 1.5),alpha = 1, colour = 'black', shape = 21, size = 2)
  
### combine only metadata
p_impacts_random_with_metadata_only = ggarrange(
  p_impacts_random_staffing+ggtitle(expression(paste(bold("a.")))),
  p_impacts_random_hh+ggtitle(expression(paste(bold("b.")))),
  p_impacts_random_kappa+ggtitle(expression(paste(bold("c.")))),
  p_impacts_random_abx+ggtitle(expression(paste(bold("d.")))),
  p_impacts_random_admR+ggtitle(expression(paste(bold("e.")))), 
  p_impacts_random_foi+ggtitle(expression(paste(bold("f.")))),
  ncol = 3, nrow = 2, align = 'hv')

ggsave("plots/tau_impacts_random_with_metadata_only.png", plot = p_impacts_random_with_metadata_only, width = 28, height = 18, unit = 'cm', )
ggsave("plots/tau_impacts_random_with_metadata_only.pdf", plot = p_impacts_random_with_metadata_only, width = 28, height = 18, unit = 'cm', )


### combine metadata with prevalence and R rate
p_impacts_random_with_metadata = ggarrange(
  p_impacts_random_staffing+ggtitle(expression(paste(bold("a."))))+
    annotate("text", x = -8, y = 1.8, label = "SARS-CoV-2 introduced", 
             angle = 90, hjust = 1, size = 3)+
    annotate("text", x = 21-8, y = 1.8, label = "COVID-19 policies instated", 
             angle = 90, hjust = 1, size = 3),
  p_impacts_random_hh+ggtitle(expression(paste(bold("b.")))),
  p_impacts_random_kappa+ggtitle(expression(paste(bold("c.")))),
  p_impacts_random_abx+ggtitle(expression(paste(bold("d.")))),
  p_impacts_random_admR+ggtitle(expression(paste(bold("e.")))), 
  p_impacts_random_foi+ggtitle(expression(paste(bold("f.")))),
  p_impacts_random_V+ylab('SARS-CoV-2\ninfection\nprevalence')+ggtitle(expression(paste(bold("g.")))), 
  p_impacts_random_B+ylab('bacterial\ncolonization\nprevalence')+ggtitle(expression(paste(bold("h.")))), 
  p_impacts_random_R_rate+ylab('bacterial\nresistance\nrate')+ggtitle(expression(paste(bold("i. ")))), 
  ncol = 3, nrow = 3, align = 'hv')

ggsave("plots/tau_impacts_random_with_metadata.png", plot = p_impacts_random_with_metadata, width = 28, height = 25, unit = 'cm', )
ggsave("plots/tau_impacts_random_with_metadata.pdf", plot = p_impacts_random_with_metadata, width = 28, height = 25, unit = 'cm', )

