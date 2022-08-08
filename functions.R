#################
### FUNCTIONS ###
#################

# last updated August 2022

###################
### ASSUMPTIONS ###
###################

### (0) ### Accessory functions (built independently into ODEs)

# Hand hygiene function
f_eta = function(H, cc, ratio){
  eta = (H*cc*ratio)
  return(eta)
}

# generic tau function
f_surge_pars = function(tau, infect, N){
  par_adjusted = qbeta(infect/N, shape1 = 2, shape2 = 1/tau)
  return(par_adjusted)
}


####################
### HOUSEKEEPING ###
####################

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



###################
### ODE SOLVING ###
###################

### (1.1) FUNCTION TO SOLVE ODEs
# (a) find endemic equilibrium of MRB given input parameters
# (b) set initial conditions to endemic equilibrium, set incidence and metadata to zero
# (c) introduce SARS-CoV-2 against backdrop of endemic bacteria, and simulate dynamics

f_dynamics_add_covid = function(ODE_model, states_base, pars_base, pars_covid, time_out=100){
  
  # (a) find endemic eqbm of ARB in absence of covid (both PA and PE)
  
  out_base = ode(y = states_base, times = c(0,10,100,1000,10000,10001),#seq(0,1000,100), 
                 func = ODE_model, 
                 parms = pars_base, method = 'bdf_d')%>% 
    as.data.frame()
  
  # number of states
  ncols_prev = 24
  
  # if model not at equilibrium, return warning
  if(sum(out_base[nrow(out_base),2:(ncols_prev+1)] - out_base[nrow(out_base)-1,2:(ncols_prev+1)])>0.01){
    warning("not at equilibrium"); print(paste0("t1: ", out_base[nrow(out_base)-1,2:(ncols_prev+1)],
                                                " and t2: ", out_base[nrow(out_base),2:(ncols_prev+1)]))
  }
  
  # (b) update input states using equilibrium values (replace 1 susceptible, uncolonized patient with 1 exposed, uncolonized patient) and set incidence to zero
  
  states_eqbm_first = round(as.numeric(out_base[nrow(out_base),2:ncol(out_base)]),6)
  
  # print(sum(out_base[nrow(out_base),2:(ncols_prev+1)]))
  
  if(sum(out_base[nrow(out_base),2:(ncols_prev+1)])>2.1){
    newCase = 1}else{
      newCase = 0.001; print("warning: N as proportion not total")
    }
  
  states_eqbm = c(S_U_pa = states_eqbm_first[1]-newCase,
                  S_Cs_pa = states_eqbm_first[2], 
                  S_Cr_pa = states_eqbm_first[3], 
                  E_U_pa = states_eqbm_first[4]+newCase, 
                  E_Cs_pa = states_eqbm_first[5], 
                  E_Cr_pa = states_eqbm_first[6], 
                  I_U_pa = states_eqbm_first[7], 
                  I_Cs_pa = states_eqbm_first[8], 
                  I_Cr_pa = states_eqbm_first[9], 
                  R_U_pa = states_eqbm_first[10], 
                  R_Cs_pa = states_eqbm_first[11], 
                  R_Cr_pa = states_eqbm_first[12],
                  S_U_pe = states_eqbm_first[13]-newCase, 
                  S_Cs_pe = states_eqbm_first[14], 
                  S_Cr_pe = states_eqbm_first[15], 
                  E_U_pe = states_eqbm_first[16]+newCase, 
                  E_Cs_pe = states_eqbm_first[17], 
                  E_Cr_pe = states_eqbm_first[18], 
                  I_U_pe = states_eqbm_first[19], 
                  I_Cs_pe = states_eqbm_first[20], 
                  I_Cr_pe = states_eqbm_first[21], 
                  R_U_pe = states_eqbm_first[22], 
                  R_Cs_pe = states_eqbm_first[23], 
                  R_Cr_pe = states_eqbm_first[24],
                  SL = 0,
                  incS_pa_pa = 0, incS_pe_pa = 0, incS_pa_pe = 0, incS_pe_pe = 0,
                  incCs_pa_pa = 0, incCs_pe_pa = 0, incCs_pa_pe = 0, incCs_pe_pe = 0,
                  incCr_pa_pa = 0, incCr_pe_pa = 0, incCr_pa_pe = 0, incCr_pe_pe = 0, incCr_endog = 0,
                  abx = 0, kappa_patients = 0, kappa_staff = 0, hh = 0, staffing_ratio = 0, adm = 0, adm_R = 0, foi_S = 0, foi_Cr = 0, patientdays = 0, staffdays = 0)
  
  # (c) simulate dynamics with introduction of covid and associated parset
  
  out_covid = ode(y = states_eqbm, 
                  times = c(0:time_out), 
                  func = ODE_model, 
                  parms = pars_covid, method = 'bdf_d')%>%
    as.data.frame()
  
  return(out_covid)
}


### (1.2) Function to make ODE ouptput data in long form, plottable via ggplot

f_dynamics_plottable = function(data_in){
  
  data_prev = data_in%>%
    dplyr::select(time, cols_prevalence)%>%
    mutate(S_pa = S_U_pa + S_Cs_pa + S_Cr_pa,
           E_pa = E_U_pa + E_Cs_pa + E_Cr_pa,
           I_pa = I_U_pa + I_Cs_pa + I_Cr_pa,
           R_pa = R_U_pa + R_Cs_pa + R_Cr_pa,
           
           U_pa = S_U_pa + E_U_pa + I_U_pa + R_U_pa,
           Cs_pa = S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa,
           Cr_pa = S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa,
           
           S_pe = S_U_pe + S_Cs_pe + S_Cr_pe,
           E_pe = E_U_pe + E_Cs_pe + E_Cr_pe,
           I_pe = I_U_pe + I_Cs_pe + I_Cr_pe,
           R_pe = R_U_pe + R_Cs_pe + R_Cr_pe,
           
           U_pe = S_U_pe + E_U_pe + I_U_pe + R_U_pe,
           Cs_pe = S_Cs_pe + E_Cs_pe + I_Cs_pe + R_Cs_pe,
           Cr_pe = S_Cr_pe + E_Cr_pe + I_Cr_pe + R_Cr_pe,
           
           Total_pa = S_U_pa + S_Cs_pa + S_Cr_pa + E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa,
           Total_pe = S_U_pe + S_Cs_pe + S_Cr_pe + E_U_pe + E_Cs_pe + E_Cr_pe + I_U_pe + I_Cs_pe + I_Cr_pe + R_U_pe + R_Cs_pe + R_Cr_pe,
           
           R_rate = (S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa)/(S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa + S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa))%>%
    pivot_longer(-c(time), values_to = "value", names_to = "state")%>%
    mutate(measure = "prevalence", 
           organism = ifelse(state %in% c('S_pa', 'E_pa', 'I_pa', 'R_pa', 'S_pe', 'E_pe', 'I_pe', 'R_pe'), 'SARS-CoV-2', ifelse(state == 'Total', 'Both', 'ARB')),
           route = NA)
  
  data_inc = data_in%>%
    dplyr::select(time, cols_incidence)%>%
    pivot_longer(-c(time), values_to = "value", names_to = "state")%>%
    mutate(measure = "incidence",
           organism = ifelse(state %in% c('incS_pa_pa', 'incS_pe_pa', 'incS_pa_pe', 'incS_pe_pe'), 'SARS-CoV-2',
                             ifelse(state %in% c('incCs_pa_pa', 'incCs_pe_pa', 'incCs_pa_pe', 'incCs_pe_pe'), 'ARB sensitive', 'ARB resistant')),
           route = ifelse(state %in% c('incS_pa_pa', 'incCs_pa_pa', 'incCr_pa_pa'), 'PA -> PA',
                          ifelse(state %in% c('incS_pa_pe', 'incCs_pa_pe', 'incCr_pa_pe'), 'PA -> PE',
                                 ifelse(state %in% c('incS_pe_pa', 'incCs_pe_pa', 'incCr_pe_pa'), 'PE -> PA',
                                        ifelse(state %in% c('incS_pe_pe', 'incCs_pe_pe', 'incCr_pe_pe'), 'PE -> PE',
                                               ifelse(state == 'incCr_endog', 'endogenous', 'ERROR'))))))
  
  data_meta = data_in%>%
    dplyr::select(time, cols_meta)%>%
    pivot_longer(-c(time), values_to = "value", names_to = "state")%>%
    mutate(measure = "metadata",
           organism = NA,
           route = NA)
  
  data_out = rbind(data_prev, data_inc, data_meta)
  
  return(data_out)
}

### (1.3) Function to simultaneously plot multiple indicators from said output

f_dynamics_plot = function(out_plottable){
  p1a = ggplot(out_plottable%>%filter(measure == 'prevalence', state %in% c("S_pa", "E_pa", "I_pa", "R_pa", "U_pa", "Cs_pa", "Cr_pa", "Total_pa")), 
               aes(x = time, y = value, colour = state, linetype = organism))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('prevalence')+
    ggtitle('patient prevalence')
  
  p1b = ggplot(out_plottable%>%filter(measure == 'prevalence', state %in% c("S_pe", "E_pe", "I_pe", "R_pe", "U_pe", "Cs_pe", "Cr_pe", "Total_pe")), 
               aes(x = time, y = value, colour = state, linetype = organism))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('prevalence')+
    ggtitle('staff prevalence')
  
  p1c = ggplot(out_plottable%>%filter(measure == 'incidence', organism == 'SARS-CoV-2'), 
               aes(x = time, y = value, colour = route))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('incidence')+
    ggtitle('SARS-CoV-2 incidence')
  
  
  p1d = ggplot(out_plottable%>%filter(measure == 'incidence', organism != 'SARS-CoV-2'), 
               aes(x = time, y = value, colour = route, linetype = organism))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('incidence')+
    scale_y_continuous(trans = 'log10')+
    ggtitle('ARB incidence')
  
  p = ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2)
  return(p)
}


#################################
### SOLVE ODEs OVER PAR RANGE ###
#################################

### (2.1) Use above functions to solve ODEs while varying a parameter value over a range, make data immediately plottable

f_dynamics_par_range = function(ODE_model, states_base, pars_base, time_out,
                                parameter, parameter_range){
  
  out = data.frame()
  
  for(i in parameter_range){
    print(paste0(parameter," ", i, " of ", last(parameter_range)))
    
    pars_covid_i = pars_base
    pars_covid_i[parameter] = i
    
    out_temp = f_dynamics_add_covid(ODE_model = ODE_model, 
                                    states_base = states_base, 
                                    pars_base = pars_base, 
                                    pars_covid = pars_covid_i, 
                                    time_out=time_out)%>%
      f_dynamics_plottable()%>%
      mutate(par_varied = parameter, par_val = i)
    
    out = rbind(out, out_temp)
  }
  
  return(out)
}



### (2.2) Plot outputs over parameter space

f_dynamics_plot_range = function(out_plottable){
  p1a = ggplot(out_plottable%>%filter(measure == 'prevalence', state %in% c("S_pa", "E_pa", "I_pa", "R_pa", "U_pa", "Cs_pa", "Cr_pa", "Total_pa")), 
               aes(x = time, y = value, colour = state, linetype = organism, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('prevalence')+
    ggtitle('patient prevalence')
  
  p1b = ggplot(out_plottable%>%filter(measure == 'prevalence', state %in% c("S_pe", "E_pe", "I_pe", "R_pe", "U_pe", "Cs_pe", "Cr_pe", "Total_pe")), 
               aes(x = time, y = value, colour = state, linetype = organism, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('prevalence')+
    ggtitle('staff prevalence')
  
  p1c = ggplot(out_plottable%>%filter(measure == 'incidence', organism == 'SARS-CoV-2'), 
               aes(x = time, y = value, colour = route, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('incidence')+
    ggtitle('SARS-CoV-2 incidence')
  
  
  p1d = ggplot(out_plottable%>%filter(measure == 'incidence', organism != 'SARS-CoV-2'), 
               aes(x = time, y = value, colour = route, linetype = organism, alpha = par_val, group = interaction(state, par_val)))+
    geom_line(stat = 'identity')+
    theme_bw()+
    ylab('incidence')+
    scale_y_continuous(trans = 'log10')+
    ggtitle('ARB incidence')
  
  p = ggarrange(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2)
  return(p)
}


#########################
### METADATA DYNAMICS ###
#########################

### For a given output dataset, calculate daily "incidence" of metadata indicators (abx, kappa_patients, kappa_staff, hh, staffing_ratio, adm_R)
## input data: wide ODE output (e.g. from f_dynamics_add_covid)

f_dynamics_metadata = function(ode_output, par = 'none', par_val = 'none'){
  
  ### if no par and par_val provided, assume these are already in the data
  
  if(par == "none" | par_val == "none"){  
    df_dynamics_metadata = ode_output%>%
      # calculate daily values of metadata indicators as row-wise difference using function diff()
      mutate(abx_daily = c(0, diff(abx)),
             kappa_patients_daily = c(0, diff(kappa_patients)),
             kappa_staff_daily = c(0, diff(kappa_staff)),
             hh_daily = c(0, diff(hh)),
             staffing_ratio_daily = c(0, diff(staffing_ratio)),
             adm_daily = c(0, diff(adm)),
             adm_R_daily = c(0, diff(adm_R)),
             foi_V_daily = c(0, diff(foi_S)),
             foi_Cr_daily = c(0, diff(foi_Cr)),
             patientdays_daily = c(0, diff(patientdays)),
             staffdays_daily = c(0, diff(staffdays))
      )%>%
      dplyr::select(time, par, par_val, abx_daily, kappa_patients_daily, kappa_staff_daily,hh_daily, staffing_ratio_daily, adm_daily, adm_R_daily, foi_V_daily, foi_Cr_daily, patientdays_daily, staffdays_daily)
    
  }else{
    
    ### Otherwise mutate the par and par_val provided by the function
    df_dynamics_metadata = ode_output%>%
      # calculate daily values of metadata indicators as row-wise difference using function diff()
      mutate(par = par, par_val = par_val)%>%
      mutate(abx_daily = c(0, diff(abx)),
             kappa_patients_daily = c(0, diff(kappa_patients)),
             kappa_staff_daily = c(0, diff(kappa_staff)),
             hh_daily = c(0, diff(hh)),
             staffing_ratio_daily = c(0, diff(staffing_ratio)),
             adm_daily = c(0, diff(adm)),
             adm_R_daily = c(0, diff(adm_R)),
             foi_V_daily = c(0, diff(foi_S)),
             foi_Cr_daily = c(0, diff(foi_Cr)),
             patientdays_daily = c(0, diff(patientdays)),
             staffdays_daily = c(0, diff(staffdays))
      )%>%
      dplyr::select(time, par, par_val, abx_daily, kappa_patients_daily, kappa_staff_daily,hh_daily, staffing_ratio_daily, adm_daily, adm_R_daily, foi_V_daily, foi_Cr_daily, patientdays_daily, staffdays_daily)
  }
  
  return(df_dynamics_metadata)
}


f_dynamics_metadata_n_sim = function(ode_output, par = 'none', par_val = 'none'){
  
  ### same as above, but when there are a number of distinct simulations to consider
  
  ### if no par and par_val provided, assume these are already in the data
  
  if(par == "none" | par_val == "none"){  
    df_dynamics_metadata = ode_output%>%
      # calculate daily values of metadata indicators as row-wise difference using function diff()
      mutate(abx_daily = c(0, diff(abx)),
             kappa_patients_daily = c(0, diff(kappa_patients)),
             kappa_staff_daily = c(0, diff(kappa_staff)),
             hh_daily = c(0, diff(hh)),
             staffing_ratio_daily = c(0, diff(staffing_ratio)),
             adm_daily = c(0, diff(adm)),
             adm_R_daily = c(0, diff(adm_R)),
             foi_V_daily = c(0, diff(foi_S)),
             foi_Cr_daily = c(0, diff(foi_Cr)),
             patientdays_daily = c(0, diff(patientdays)),
             staffdays_daily = c(0, diff(staffdays))
      )%>%
      dplyr::select(n_sim, time, par, par_val, abx_daily, kappa_patients_daily, kappa_staff_daily,hh_daily, staffing_ratio_daily, adm_daily, adm_R_daily, foi_V_daily, foi_Cr_daily, patientdays_daily, staffdays_daily)
    
  }else{
    
    ### Otherwise mutate the par and par_val provided by the function
    df_dynamics_metadata = ode_output%>%
      # calculate daily values of metadata indicators as row-wise difference using function diff()
      mutate(par = par, par_val = par_val)%>%
      mutate(abx_daily = c(0, diff(abx)),
             kappa_patients_daily = c(0, diff(kappa_patients)),
             kappa_staff_daily = c(0, diff(kappa_staff)),
             hh_daily = c(0, diff(hh)),
             staffing_ratio_daily = c(0, diff(staffing_ratio)),
             adm_daily = c(0, diff(adm)),
             adm_R_daily = c(0, diff(adm_R)),
             foi_V_daily = c(0, diff(foi_S)),
             foi_Cr_daily = c(0, diff(foi_Cr)),
             patientdays_daily = c(0, diff(patientdays)),
             staffdays_daily = c(0, diff(staffdays))
      )%>%
      dplyr::select(n_sim, time, par, par_val, abx_daily, kappa_patients_daily, kappa_staff_daily,hh_daily, staffing_ratio_daily, adm_daily, adm_R_daily, foi_V_daily, foi_Cr_daily, patientdays_daily, staffdays_daily)
  }
  
  return(df_dynamics_metadata)
}
########################################################
### LOOP THROUGH PANDEMIC IMPACTS WITH RANDOM VALUES ###
########################################################

f_dynamics_loop_random_tau = function(model_loop, pars_loop, states_loop, time_out, n_parsets, save_data = 'NO'){
  
  # load random parameter values
  df_tau_random_vals = loadRData("data/outputs_analysis1/df_tau_random_vals.Rdata")
  
  # initialize empty dataframe
  df_tau_out = data.frame()
  
  for(parset in 1:n_parsets){
    
    print(paste0("For random tau, evaluating parameter set ", parset, " of ", n_parsets))
    
    pars_loop_i = pars_loop
    pars_tau_i = df_tau_random_vals[parset,]
    pars_loop_i["handrub"] = as.numeric(pars_tau_i["handrub"])
    pars_loop_i["masks"] = as.numeric(pars_tau_i["masks"])
    pars_loop_i["prophylaxis"] = as.numeric(pars_tau_i["prophylaxis"])
    pars_loop_i["distancing"] = as.numeric(pars_tau_i["distancing"])
    pars_loop_i["disorg"] = as.numeric(pars_tau_i["disorg"])
    pars_loop_i["comm_denom"] = as.numeric(pars_tau_i["comm_denom"])
    pars_loop_i["a_surge"] = as.numeric(pars_tau_i["a_surge"])
    pars_loop_i["chi"] = as.numeric(pars_tau_i["chi"])
    pars_loop_i["covid_stay"] = as.numeric(pars_tau_i["covid_stay"])
    pars_loop_i["adm_reduc"] = as.numeric(pars_tau_i["adm_reduc"])
    
    df_out_loop = f_dynamics_add_covid(ODE_model = model_loop, 
                                       states_base = states_loop, 
                                       pars_base = pars_loop, 
                                       pars_covid = pars_loop_i,
                                       time_out = time_out)%>%
      mutate(n_sim = parset)
    
    df_tau_out = rbind(df_tau_out, df_out_loop)
  }
  
  if(save_data == 'YES'){
    print(paste0('saving random tau output data'))
    # don't save direct to raw data folder : save to initial greater folder so as not to accidentally overwrite
    save(df_tau_out, file = paste0("data/df_tau_random.Rdata"))
  }else{
    return(df_tau_out)
  }
}


###########################################
### SOLVE FOR GAMUT OF PANDEMIC IMPACTS ###
###########################################
# loop through "list_impacts" -- list of all tau and combinations of tau -- and solve for a given equal value of tau

f_dynamics_loop_tau = function(model_loop, pars_loop, states_loop, time_out, list_impacts, save_data = 'NO'){
  
  # initialize empty data.frame and loop counter
  df_tau = data.frame()
  qounter = 0
  
  # loop through impacts and solve ODEs 
  for(impact in list_impacts){
    
    # extract name of pandemic impact(s) evaluated
    qounter = qounter+1
    impact_name = names(list_impacts)[qounter]
    
    # loop through impact strengths
    for(par_val in seq(0,0.999,length.out = 7)){
      
      print(paste("For the parameter ", impact_name, ", evaluating ", par_val))
      
      pars_covid = pars_loop
      
      # update value for each pandemic impact parameter considered 
      for(par in impact){
        pars_covid[par] = par_val
      }
      
      df_loop = f_dynamics_add_covid(ODE_model = model_loop, 
                                     states_base = states_loop, 
                                     pars_base = pars_loop, 
                                     pars_covid = pars_covid,
                                     time_out = time_out)%>%
        mutate(par = impact_name, par_val = par_val)
      
      df_tau = rbind(df_tau, df_loop)
      
    }
  }
  
  
  if(save_data == 'YES'){
    print(paste0('saving ', impact_name))
    # don't save direct to raw data folder : save to initial greater folder so as not to accidentally overwrite
    save(df_tau, file = paste0("data/df_tau_impacts.Rdata"))
  }else{
    return(df_tau)
  }
}

## Make these data plottable, accounting for varying tau
f_dynamics_plottable_tau = function(data_in){
  
  data_prev = data_in%>%
    dplyr::select(time, cols_prevalence, par, par_val)%>%
    mutate(S_pa = S_U_pa + S_Cs_pa + S_Cr_pa,
           E_pa = E_U_pa + E_Cs_pa + E_Cr_pa,
           I_pa = I_U_pa + I_Cs_pa + I_Cr_pa,
           R_pa = R_U_pa + R_Cs_pa + R_Cr_pa,
           
           U_pa = S_U_pa + E_U_pa + I_U_pa + R_U_pa,
           Cs_pa = S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa,
           Cr_pa = S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa,
           
           S_pe = S_U_pe + S_Cs_pe + S_Cr_pe,
           E_pe = E_U_pe + E_Cs_pe + E_Cr_pe,
           I_pe = I_U_pe + I_Cs_pe + I_Cr_pe,
           R_pe = R_U_pe + R_Cs_pe + R_Cr_pe,
           
           U_pe = S_U_pe + E_U_pe + I_U_pe + R_U_pe,
           Cs_pe = S_Cs_pe + E_Cs_pe + I_Cs_pe + R_Cs_pe,
           Cr_pe = S_Cr_pe + E_Cr_pe + I_Cr_pe + R_Cr_pe,
           
           Total_pa = S_U_pa + S_Cs_pa + S_Cr_pa + E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa,
           Total_pe = S_U_pe + S_Cs_pe + S_Cr_pe + E_U_pe + E_Cs_pe + E_Cr_pe + I_U_pe + I_Cs_pe + I_Cr_pe + R_U_pe + R_Cs_pe + R_Cr_pe,
           
           R_rate = (S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa)/(S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa + S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa))%>%
    pivot_longer(-c(time, par, par_val), values_to = "value", names_to = "state")%>%
    mutate(measure = "prevalence", 
           organism = ifelse(state %in% c('S_pa', 'E_pa', 'I_pa', 'R_pa', 'S_pe', 'E_pe', 'I_pe', 'R_pe'), 'SARS-CoV-2', ifelse(state == 'Total', 'Both', 'ARB')),
           route = NA)
  
  data_inc = data_in%>%
    dplyr::select(time, cols_incidence, par, par_val)%>%
    pivot_longer(-c(time, par, par_val), values_to = "value", names_to = "state")%>%
    mutate(measure = "incidence",
           organism = ifelse(state %in% c('incS_pa_pa', 'incS_pe_pa', 'incS_pa_pe', 'incS_pe_pe'), 'SARS-CoV-2',
                             ifelse(state %in% c('incCs_pa_pa', 'incCs_pe_pa', 'incCs_pa_pe', 'incCs_pe_pe'), 'ARB sensitive', 'ARB resistant')),
           route = ifelse(state %in% c('incS_pa_pa', 'incCs_pa_pa', 'incCr_pa_pa'), 'PA -> PA',
                          ifelse(state %in% c('incS_pa_pe', 'incCs_pa_pe', 'incCr_pa_pe'), 'PA -> PE',
                                 ifelse(state %in% c('incS_pe_pa', 'incCs_pe_pa', 'incCr_pe_pa'), 'PE -> PA',
                                        ifelse(state %in% c('incS_pe_pe', 'incCs_pe_pe', 'incCr_pe_pe'), 'PE -> PE',
                                               ifelse(state == 'incCr_endog', 'endogenous', 'ERROR'))))))
  
  data_meta = data_in%>%
    dplyr::select(time, cols_meta, par, par_val)%>%
    pivot_longer(-c(time, par, par_val), values_to = "value", names_to = "state")%>%
    mutate(measure = "metadata",
           organism = NA,
           route = NA)
  
  data_out = rbind(data_prev, data_inc, data_meta)
  
  return(data_out)
}



##################
### INDICATORS ###
##################

### Take ODE output and calculate key indicators, returning one row per simulation

# specify input data that are key for tracking simulation output (n_sim, tau, tau_val)
f_indicators = function(n_sim, tau, tau_val, ode_output){
  
  # calculate raw outputs for each t in simulation time by grouping compartments as needed
  ode_output_grouped = ode_output%>%
    mutate(S_pa = S_U_pa + S_Cs_pa + S_Cr_pa,
           E_pa = E_U_pa + E_Cs_pa + E_Cr_pa,
           I_pa = I_U_pa + I_Cs_pa + I_Cr_pa,
           R_pa = R_U_pa + R_Cs_pa + R_Cr_pa,
           
           U_pa = S_U_pa + E_U_pa + I_U_pa + R_U_pa,
           Cs_pa = S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa,
           Cr_pa = S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa,
           
           S_pe = S_U_pe + S_Cs_pe + S_Cr_pe,
           E_pe = E_U_pe + E_Cs_pe + E_Cr_pe,
           I_pe = I_U_pe + I_Cs_pe + I_Cr_pe,
           R_pe = R_U_pe + R_Cs_pe + R_Cr_pe,
           
           U_pe = S_U_pe + E_U_pe + I_U_pe + R_U_pe,
           Cs_pe = S_Cs_pe + E_Cs_pe + I_Cs_pe + R_Cs_pe,
           Cr_pe = S_Cr_pe + E_Cr_pe + I_Cr_pe + R_Cr_pe,
           
           N_pa = S_U_pa + S_Cs_pa + S_Cr_pa + E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa,
           N_pe = S_U_pe + S_Cs_pe + S_Cr_pe + E_U_pe + E_Cs_pe + E_Cr_pe + I_U_pe + I_Cs_pe + I_Cr_pe + R_U_pe + R_Cs_pe + R_Cr_pe,
           
           R_rate = (S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa)/(S_Cs_pa + E_Cs_pa + I_Cs_pa + R_Cs_pa + S_Cr_pa + E_Cr_pa + I_Cr_pa + R_Cr_pa),
           
           incS_all = incS_pa_pa + incS_pa_pe + incS_pe_pa + incS_pe_pe,
           incCs_all = incCs_pa_pa + incCs_pe_pa, ### ONLY COLONIZATION (not transient carriage)
           incCr_all = incCr_pa_pa + incCr_pe_pa + incCr_endog,
           incTs_all = incCs_pa_pe + incCs_pe_pe, ### ONLY TRANSIENT CARRIAGE (not colonization)
           incTr_all = incCr_pa_pe + incCr_pe_pe)%>%
    ### calculate prevalence (proportional)
    mutate(V_prev_pa = (I_pa + E_pa)/N_pa,
           V_prev_pa_I = I_pa/N_pa,
           V_prev_pe = (I_pe + E_pe)/N_pe,
           V_prev_pe_I = I_pe/N_pe,
           V_prev_all = (I_pa + E_pa + I_pe + E_pe)/(N_pa + N_pe),
           V_prev_all_I = (I_pa + I_pe)/(N_pa + N_pe),
           Cr_prev = Cr_pa/N_pa,
           Cs_prev = Cs_pa/N_pa,
           Tr_prev = Cr_pe/N_pe,
           Ts_prev = Cs_pe/N_pe
           )%>%
    ### calculate daily incidence as row-wise difference in cumulative incidence using function diff()
    mutate(V_inc_daily_pa_pa = c(0, diff(incS_pa_pa)),
           V_inc_daily_pa_pe = c(0, diff(incS_pa_pe)),
           V_inc_daily_pe_pa = c(0, diff(incS_pe_pa)),
           V_inc_daily_pe_pe = c(0, diff(incS_pe_pe)),
           V_inc_daily_all = c(0, diff(incS_all)),
           Cs_inc_daily_pa_pa = c(0, diff(incCs_pa_pa)),
           Cs_inc_daily_pe_pa = c(0, diff(incCs_pe_pa)), 
           Cs_inc_daily_all = c(0, diff(incCs_all)), 
           Cr_inc_daily_pa_pa = c(0, diff(incCr_pa_pa)),
           Cr_inc_daily_pe_pa = c(0, diff(incCr_pe_pa)), 
           Cr_inc_daily_endog = c(0, diff(incCr_endog)),
           Cr_inc_daily_all = c(0, diff(incCr_all)), 
           Ts_inc_daily_pa_pe = c(0, diff(incCs_pa_pe)), # number of patient contaminations of HCWs with Bs
           Ts_inc_daily_pe_pe = c(0, diff(incCs_pe_pe)), # number of HCW contaminations of other HCWs with Bs
           Ts_inc_daily_all = c(0, diff(incTs_all)),
           Tr_inc_daily_pa_pe = c(0, diff(incCr_pa_pe)), # number of patient contaminations of HCWs with Br
           Tr_inc_daily_pe_pe = c(0, diff(incCr_pe_pe)), # number of HCW contaminations of other HCWs with Br
           Tr_inc_daily_all = c(0, diff(incTr_all))
           )%>%
    ### calculate metadata indicators as row-wise difference in cumulative indicators using diff() 
    mutate(abx = c(0, diff(abx)),
           kappa_patients = c(0, diff(kappa_patients)),
           kappa_staff = c(0, diff(kappa_staff)),
           hh = c(0, diff(hh)),
           staffing_ratio = c(0, diff(staffing_ratio)),
           adm_daily = c(0, diff(adm)),
           adm_R = c(0, diff(adm_R)),
           foi_S = c(0, diff(foi_S)),
           foi_Cr = c(0, diff(foi_Cr)),
           patientdays = c(0, diff(patientdays)),
           staffdays = c(0, diff(staffdays))
           )
  

  # data.frame containing all key indicators
  df_indicators = data.frame(n_sim = n_sim,
                             tau = tau,
                             tau_val = tau_val,
                             ### SARS-CoV-2 indicators
                             # max prevalence
                             V_prev_max_pa = max(ode_output_grouped$V_prev_pa),
                             V_prev_max_pa_I = max(ode_output_grouped$V_prev_pa_I),
                             V_prev_max_pe = max(ode_output_grouped$V_prev_pe),
                             V_prev_max_pe_I = max(ode_output_grouped$V_prev_pe_I),
                             V_prev_max_all = max(ode_output_grouped$V_prev_all),
                             V_prev_max_all_I = max(ode_output_grouped$V_prev_all_I),
                             # cumulative incidence
                             V_inc_cumul_pa_pa = last(ode_output_grouped$incS_pa_pa),
                             V_inc_cumul_pa_pe = last(ode_output_grouped$incS_pa_pe),
                             V_inc_cumul_pe_pa = last(ode_output_grouped$incS_pe_pa),
                             V_inc_cumul_pe_pe = last(ode_output_grouped$incS_pe_pe),
                             V_inc_cumul_all = last(ode_output_grouped$incS_all),
                             # max daily incidence
                             V_inc_daily_max_pa_pa = max(ode_output_grouped$V_inc_daily_pa_pa),
                             V_inc_daily_max_pa_pe = max(ode_output_grouped$V_inc_daily_pa_pe),
                             V_inc_daily_max_pe_pa = max(ode_output_grouped$V_inc_daily_pe_pa),
                             V_inc_daily_max_pe_pe = max(ode_output_grouped$V_inc_daily_pe_pe),
                             V_inc_daily_max_all = max(ode_output_grouped$V_inc_daily_all),
                             # lag to peak prevalence
                             V_lag_peak_prev_pa = which.max(ode_output_grouped$V_prev_pa),
                             V_lag_peak_prev_pa_I = which.max(ode_output_grouped$V_prev_pa_I),
                             V_lag_peak_prev_pe = which.max(ode_output_grouped$V_prev_pe),
                             V_lag_peak_prev_pe_I = which.max(ode_output_grouped$V_prev_pe_I),
                             V_lag_peak_prev_all = which.max(ode_output_grouped$V_prev_all),
                             V_lag_peak_prev_all_I = which.max(ode_output_grouped$V_prev_all_I),
                             # lag to peak daily incidence
                             V_lag_peak_inc_pa_pa = which.max(ode_output_grouped$V_inc_daily_pa_pa),
                             V_lag_peak_inc_pa_pe = which.max(ode_output_grouped$V_inc_daily_pa_pe),
                             V_lag_peak_inc_pe_pa = which.max(ode_output_grouped$V_inc_daily_pe_pa),
                             V_lag_peak_inc_pe_pe = which.max(ode_output_grouped$V_inc_daily_pe_pe),
                             V_lag_peak_inc_all = which.max(ode_output_grouped$V_inc_daily_all),
                             
                             ### Bacterial colonization indicators
                             # max prevalence
                             Cs_prev_max = max(ode_output_grouped$Cs_prev),
                             Cr_prev_max = max(ode_output_grouped$Cr_prev),
                             Ts_prev_max = max(ode_output_grouped$Ts_prev),
                             Tr_prev_max = max(ode_output_grouped$Tr_prev),
                             # min prevalence
                             Cs_prev_min = min(ode_output_grouped$Cs_prev),
                             Cr_prev_min = min(ode_output_grouped$Cr_prev),
                             Ts_prev_min = min(ode_output_grouped$Ts_prev),
                             Tr_prev_min = min(ode_output_grouped$Tr_prev),
                             # person-days colonized
                             Cs_pd = sum(ode_output_grouped$Cs_pa),
                             Cr_pd = sum(ode_output_grouped$Cr_pa),
                             # cumulative incidence (colonization)
                             Cs_inc_cumul_pa_pa = last(ode_output_grouped$incCs_pa_pa),
                             Cs_inc_cumul_pe_pa = last(ode_output_grouped$incCs_pe_pa),
                             Cs_inc_cumul_all = last(ode_output_grouped$incCs_all),
                             Cr_inc_cumul_pa_pa = last(ode_output_grouped$incCr_pa_pa),
                             Cr_inc_cumul_pe_pa = last(ode_output_grouped$incCr_pe_pa),
                             Cr_inc_cumul_endog = last(ode_output_grouped$incCr_endog),
                             Cr_inc_cumul_all = last(ode_output_grouped$incCr_all),
                             # cumulative incidence (transient carriage)
                             Ts_inc_cumul_pa_pe = last(ode_output_grouped$incCs_pa_pe),
                             Ts_inc_cumul_pe_pe = last(ode_output_grouped$incCs_pe_pe),
                             Ts_inc_cumul_all = last(ode_output_grouped$incTs_all),
                             Tr_inc_cumul_pa_pe = last(ode_output_grouped$incCr_pa_pe),
                             Tr_inc_cumul_pe_pe = last(ode_output_grouped$incCr_pe_pe),
                             Tr_inc_cumul_all = last(ode_output_grouped$incTr_all),
                             # max daily incidence (colonization)
                             Cs_inc_daily_max_pa_pa = max(ode_output_grouped$Cs_inc_daily_pa_pa),
                             Cs_inc_daily_max_pe_pa = max(ode_output_grouped$Cs_inc_daily_pe_pa),
                             Cs_inc_daily_max_all = max(ode_output_grouped$Cs_inc_daily_all),
                             Cr_inc_daily_max_pa_pa = max(ode_output_grouped$Cr_inc_daily_pa_pa),
                             Cr_inc_daily_max_pe_pa = max(ode_output_grouped$Cr_inc_daily_pe_pa),
                             Cr_inc_daily_max_endog = max(ode_output_grouped$Cr_inc_daily_endog),
                             Cr_inc_daily_max_all = max(ode_output_grouped$Cr_inc_daily_all),
                             # max daily incidence (transient carriage)
                             Ts_inc_daily_max_pa_pe = max(ode_output_grouped$Ts_inc_daily_pa_pe),
                             Ts_inc_daily_max_pe_pe = max(ode_output_grouped$Ts_inc_daily_pe_pe),
                             Ts_inc_daily_max_all = max(ode_output_grouped$Ts_inc_daily_all),
                             Tr_inc_daily_max_pa_pe = max(ode_output_grouped$Tr_inc_daily_pa_pe),
                             Tr_inc_daily_max_pe_pe = max(ode_output_grouped$Tr_inc_daily_pe_pe),
                             Tr_inc_daily_max_all = max(ode_output_grouped$Tr_inc_daily_all),
                             # R_rate
                             R_rate_max = max(ode_output_grouped$R_rate),
                             R_rate_min = min(ode_output_grouped$R_rate),
                             # R_rate, overall across person-days
                             R_rate_pd = sum(ode_output_grouped$Cr_pa)/(sum(ode_output_grouped$Cr_pa)+sum(ode_output_grouped$Cs_pa)),
                             # incidence rate, over patient days (colonization)
                             Cs_incRate_pa_pa = last(ode_output_grouped$incCs_pa_pa)/last(ode_output_grouped$patientdays),
                             Cs_incRate_pe_pa = last(ode_output_grouped$incCs_pe_pa)/last(ode_output_grouped$patientdays),
                             Cs_incRate_all = last(ode_output_grouped$incCs_all)/last(ode_output_grouped$patientdays),
                             Cr_incRate_pa_pa = last(ode_output_grouped$incCr_pa_pa)/last(ode_output_grouped$patientdays),
                             Cr_incRate_pe_pa = last(ode_output_grouped$incCr_pe_pa)/last(ode_output_grouped$patientdays),
                             Cr_incRate_endog = last(ode_output_grouped$incCr_endog)/last(ode_output_grouped$patientdays),
                             Cr_incRate_all = last(ode_output_grouped$incCr_all)/last(ode_output_grouped$patientdays),
                             # incidence rate, over staff days (carriage)
                             Ts_incRate_pa_pe = last(ode_output_grouped$incCs_pa_pe)/last(ode_output_grouped$staffdays),
                             Ts_incRate_pe_pe = last(ode_output_grouped$incCs_pe_pe)/last(ode_output_grouped$staffdays),
                             Ts_incRate_all = last(ode_output_grouped$incTs_all)/last(ode_output_grouped$staffdays),
                             Tr_incRate_pa_pe = last(ode_output_grouped$incCr_pa_pe)/last(ode_output_grouped$staffdays),
                             Tr_incRate_pe_pe = last(ode_output_grouped$incCr_pe_pe)/last(ode_output_grouped$staffdays),
                             Tr_incRate_all = last(ode_output_grouped$incTr_all)/last(ode_output_grouped$staffdays),
                             # metadata indicators
                             time = last(ode_output_grouped$time),
                             
                             meta_abx_total = sum(ode_output_grouped$abx),
                             meta_abx_max = max(ode_output_grouped$abx),
                             meta_abx_min = min(ode_output_grouped$abx[2:nrow(ode_output_grouped)]),
                             
                             meta_kappa_patients_total = sum(ode_output_grouped$kappa_patients),
                             meta_kappa_patients_max = max(ode_output_grouped$kappa_patients),
                             meta_kappa_patients_min = min(ode_output_grouped$kappa_patients[2:nrow(ode_output_grouped)]),
                             
                             meta_kappa_staff_total = sum(ode_output_grouped$kappa_staff),
                             meta_kappa_staff_max = max(ode_output_grouped$kappa_staff),
                             meta_kappa_staff_min = min(ode_output_grouped$kappa_staff[2:nrow(ode_output_grouped)]),
                             
                             meta_hh_total = sum(ode_output_grouped$hh),
                             meta_hh_max = max(ode_output_grouped$hh),
                             meta_hh_min = min(ode_output_grouped$hh[2:nrow(ode_output_grouped)]),
                             
                             meta_staffing_ratio_total = sum(ode_output_grouped$staffing_ratio),
                             meta_staffing_ratio_max = max(ode_output_grouped$staffing_ratio),
                             meta_staffing_ratio_min = min(ode_output_grouped$staffing_ratio[2:nrow(ode_output_grouped)]),
                             
                             meta_adm_total = sum(ode_output_grouped$adm),
                             meta_adm_max = max(ode_output_grouped$adm),
                             meta_adm_min = min(ode_output_grouped$adm[2:nrow(ode_output_grouped)]),
                             
                             meta_adm_R_total = sum(ode_output_grouped$adm_R),
                             meta_adm_R_max = max(ode_output_grouped$adm_R),
                             meta_adm_R_min = min(ode_output_grouped$adm_R[2:nrow(ode_output_grouped)]),
                             
                             meta_foi_S_total = sum(ode_output_grouped$foi_S),
                             meta_foi_S_max = max(ode_output_grouped$foi_S),
                             meta_foi_S_min = min(ode_output_grouped$foi_S[2:nrow(ode_output_grouped)]),
                             
                             meta_foi_Cr_total = sum(ode_output_grouped$foi_Cr),
                             meta_foi_Cr_max = max(ode_output_grouped$foi_Cr),
                             meta_foi_Cr_min = min(ode_output_grouped$foi_Cr[2:nrow(ode_output_grouped)]),
                             
                             meta_patientdays_total = sum(ode_output_grouped$patientdays),
                             meta_patientdays_max = max(ode_output_grouped$patientdays),
                             meta_patientdays_min = min(ode_output_grouped$patientdays[2:nrow(ode_output_grouped)]),
                             
                             meta_staffdays_total = sum(ode_output_grouped$staffdays),
                             meta_staffdays_max = max(ode_output_grouped$staffdays),
                             meta_staffdays_min = min(ode_output_grouped$staffdays[2:nrow(ode_output_grouped)])
                             
  )
  
  # and return this huge df of indicators
  return(df_indicators)
}




####################################
### PARAMETER SAMPLING FUNCTIONS ###
####################################
### Functions to create vectors of states and parameters


##############
### STATES ###
##############

f_states = function(pars_r){
  
  Nbeds_r = pars_r['Nbeds']
  Nhcw_r = pars_r['Nhcw']
  
  states_r <- c(S_U_pa = as.numeric(Nbeds_r)/3, S_Cs_pa = as.numeric(Nbeds_r)/3, S_Cr_pa = as.numeric(Nbeds_r)/3, E_U_pa = 0, E_Cs_pa = 0, E_Cr_pa = 0, I_U_pa = 0, I_Cs_pa = 0, I_Cr_pa = 0, R_U_pa = 0, R_Cs_pa = 0, R_Cr_pa = 0,
                S_U_pe = as.numeric(Nhcw_r)/3, S_Cs_pe = as.numeric(Nhcw_r)/3, S_Cr_pe = as.numeric(Nhcw_r)/3, E_U_pe = 0, E_Cs_pe = 0, E_Cr_pe = 0, I_U_pe = 0, I_Cs_pe = 0, I_Cr_pe = 0, R_U_pe = 0, R_Cs_pe = 0, R_Cr_pe = 0,
                SL = 0,
                incS_pa_pa = 0, incS_pe_pa = 0, incS_pa_pe = 0, incS_pe_pe = 0,
                incCs_pa_pa = 0, incCs_pe_pa = 0, incCs_pa_pe = 0, incCs_pe_pe = 0,
                incCr_pa_pa = 0, incCr_pe_pa = 0, incCr_pa_pe = 0, incCr_pe_pe = 0, incCr_endog = 0,
                abx = 0, kappa_patients = 0, kappa_staff = 0, hh = 0, staffing_ratio = 0, adm = 0, adm_R = 0, foi_S = 0, foi_Cr = 0, patientdays = 0, staffdays = 0)
  
  ### return updated par set
  return(states_r)
}

#################
## PARAMETERS ###
#################


### ANALYSIS2 ###
### Uncertainty in both facility and bacterial characteristics parameters, built upon baseline par set pars_m

f_params_uncertain_facility_bacteria = function(pars_m){
  
  ### healthcare facility parameters
  r_mu = 1/rlnorm(1,log(10),log(3))
  r_Nbeds = runif(1,100,1000)
  r_Nhcw = r_Nbeds*runif(1,0.25,3)
  r_kappa_pa_pa = rlnorm(1,log(5),log(2))
  r_kappa_pa_pe = rlnorm(1,log(15),log(2))
  r_kappa_pe_pe = rlnorm(1,log(15),log(2))
  r_kappa_pe_pa = r_kappa_pa_pe*(r_Nbeds/r_Nhcw)
  r_rho = runif(1,1,5)
  r_hyg = runif(1,0.2,0.6)
  r_a = rbeta(1,1.5,5)

  
  ### bacteria parameters
  r_gamma_Cs = rlnorm(1, log(0.03), log(3))
  r_cost = rbeta(1, 0.5, 2)
  r_gamma_Cr = r_gamma_Cs*(1+r_cost)
  r_r_S = runif(1,0,0.25)
  r_r_R = runif(1,0.75,1)
  r_beta_Cs = rlnorm(1, log(0.2), log(1.5))
  r_beta_Cr = r_beta_Cs
  r_theta = rbeta(1, 0.8, 2)
  r_alpha = rgamma(1, 0.4, 10)
  r_f_Cs = runif(1,0.05,0.3)
  r_f_Cr = runif(1,0.001,0.15)
  r_f_U = 1 - r_f_Cs - r_f_Cr
  
  
  ### transmission rates per contact
  r_beta_S = as.numeric(pars_m['beta_S'])
  r_pi_S = r_beta_S/(r_kappa_pa_pa+r_kappa_pa_pe+r_kappa_pe_pe+r_kappa_pe_pa)
  r_pi_Cs = r_beta_Cs/(r_kappa_pa_pa+r_kappa_pe_pa)
  r_pi_Cr = r_beta_Cr/(r_kappa_pa_pa+r_kappa_pe_pa)
  
  
  ### all other parameters use baseline values given in pars_states, i.e. "m_par"
  
  pars_r = pars_m;
  
  ### update parameters
  pars_r['mu'] = r_mu;
  pars_r['kappa_pa_pa'] = r_kappa_pa_pa;
  pars_r['kappa_pa_pe'] = r_kappa_pa_pe;
  pars_r['kappa_pe_pe'] = r_kappa_pe_pe;
  pars_r['rho'] = r_rho;
  pars_r['hyg'] = r_hyg;
  pars_r['a'] = r_a;
  pars_r['Nbeds'] = r_Nbeds;
  pars_r['Nhcw'] = r_Nhcw;
  pars_r['gamma_Cs'] = r_gamma_Cs;
  pars_r['cost'] = r_cost;
  pars_r['gamma_Cr'] = r_gamma_Cr;
  pars_r['r_S'] = r_r_S;
  pars_r['r_R'] = r_r_R;
  pars_r['theta'] = r_theta;
  pars_r['alpha'] = r_alpha;
  pars_r['f_Cs'] = r_f_Cs;
  pars_r['f_Cr'] = r_f_Cr;
  pars_r['f_U'] = r_f_U;
  pars_r['pi_S'] = r_pi_S;
  pars_r['pi_Cs'] = r_pi_Cs;
  pars_r['pi_Cr'] = r_pi_Cr;
  
  ### return updated par set
  return(pars_r)
}



### ANALYSIS3 ###
### Healthcare setting parameters for rehabilitation hospital
f_params_rehab = function(pars_in){
  
  ### healthcare facility parameters
  r_mu = 1/rnorm(1,49,4.9)
  r_Nbeds = rnorm(1,30,3)
  r_Nhcw = rnorm(1,20,2)
  r_kappa_pa_pa = rnorm(1,7.1,0.71)
  r_kappa_pa_pe = rnorm(1,3.9,0.39)
  r_kappa_pe_pe = rnorm(1,4.9, 0.49)
  r_kappa_pe_pa = r_kappa_pa_pe*(r_Nbeds/r_Nhcw)
  r_rho = runif(1,1,3)
  r_hyg = rnorm(1,0.12,0.012)
  r_a = rnorm(1,0.049,0.0049)
  
  ### base
  pars_rehab = pars_in;
  
  ### update parameters
  pars_rehab['mu'] = r_mu;
  pars_rehab['kappa_pa_pa'] = r_kappa_pa_pa;
  pars_rehab['kappa_pa_pe'] = r_kappa_pa_pe;
  pars_rehab['kappa_pe_pe'] = r_kappa_pe_pe;
  pars_rehab['kappa_pe_pa'] = r_kappa_pe_pa;
  pars_rehab['rho'] = r_rho;
  pars_rehab['hyg'] = r_hyg;
  pars_rehab['a'] = r_a;
  pars_rehab['Nbeds'] = r_Nbeds;
  pars_rehab['Nhcw'] = r_Nhcw;
  
  ### return updated par set
  return(pars_rehab)
}

### Healthcare setting parameters for short-stay geriatric unit of tertiary hospital
f_params_geriatric = function(pars_in){
  
  ### healthcare facility parameters
  r_mu = 1/rnorm(1,2.5,0.25)
  r_Nbeds = rnorm(1,19,1.9)
  r_Nhcw = rnorm(1,38,3.8)
  r_kappa_pa_pa = rnorm(1,1,0.1)
  r_kappa_pa_pe = rnorm(1,30,3)
  r_kappa_pe_pe = rnorm(1,62, 6.2)
  r_kappa_pe_pa = r_kappa_pa_pe*(r_Nbeds/r_Nhcw)
  r_rho = runif(1,2,4)
  r_hyg = rnorm(1,0.46,0.046)
  r_a = rnorm(1,0.305,0.0305)
  
  ### base
  pars_geriatric = pars_in;
  
  ### update parameters
  pars_geriatric['mu'] = r_mu;
  pars_geriatric['kappa_pa_pa'] = r_kappa_pa_pa;
  pars_geriatric['kappa_pa_pe'] = r_kappa_pa_pe;
  pars_geriatric['kappa_pe_pe'] = r_kappa_pe_pe;
  pars_geriatric['kappa_pe_pa'] = r_kappa_pe_pa;
  pars_geriatric['rho'] = r_rho;
  pars_geriatric['hyg'] = r_hyg;
  pars_geriatric['a'] = r_a;
  pars_geriatric['Nbeds'] = r_Nbeds;
  pars_geriatric['Nhcw'] = r_Nhcw;
  
  ### return updated par set
  return(pars_geriatric)
}

### Healthcare setting parameters for general ward of tertiary paediatric hospital
f_params_paediatric = function(pars_in){
  
  ### healthcare facility parameters
  r_mu = 1/rnorm(1,7,0.7)
  r_Nbeds = rnorm(1,37,3.7)
  r_Nhcw = rnorm(1,36,3.6)
  r_kappa_pa_pa = rnorm(1,0.1,0.01)
  r_kappa_pa_pe = rnorm(1,1.3,0.13)
  r_kappa_pe_pe = rnorm(1,34,3.4)
  r_kappa_pe_pa = r_kappa_pa_pe*(r_Nbeds/r_Nhcw)
  r_rho = runif(1,3,5)
  r_hyg = rnorm(1,0.46,0.046)
  r_a = rnorm(1,0.395,0.0395)
  
  ### base
  pars_paediatric = pars_in;
  
  ### update parameters
  pars_paediatric['mu'] = r_mu;
  pars_paediatric['kappa_pa_pa'] = r_kappa_pa_pa;
  pars_paediatric['kappa_pa_pe'] = r_kappa_pa_pe;
  pars_paediatric['kappa_pe_pe'] = r_kappa_pe_pe;
  pars_paediatric['kappa_pe_pa'] = r_kappa_pe_pa;
  pars_paediatric['rho'] = r_rho;
  pars_paediatric['hyg'] = r_hyg;
  pars_paediatric['a'] = r_a;
  pars_paediatric['Nbeds'] = r_Nbeds;
  pars_paediatric['Nhcw'] = r_Nhcw;
  
  ### return updated par set
  return(pars_paediatric)
}


### Pandemic parameter for organized response
f_params_organized = function(pars_in){
  
  ### pandemic impacts
  r_prophylaxis = runif(1, 0.047, 0.152);
  r_a_surge = 0;
  r_masks = 1; while(r_masks>=1){r_masks = rnorm(1,0.99,0.003) * runif(1,0.5,1)}
  r_handrub = runif(1, 0.25, 0.3);
  r_disorg = runif(1, 0.4, 0.6);
  r_distancing = runif(1, 0.8, 1);
  r_covid_stay = runif(1, 0.8, 1);
  r_chi = runif(1, 0, 0.5);
  r_adm_reduc = runif(1, 0.4, 0.6);
  r_comm_denom = runif(1, 0.4, 0.6);
  
  ### base
  pars_organized = pars_in;
  
  ### update parameters
  pars_organized['prophylaxis'] = r_prophylaxis;
  pars_organized['a_surge'] = r_a_surge;
  pars_organized['masks'] = r_masks;
  pars_organized['handrub'] = r_handrub;
  pars_organized['disorg'] = r_disorg;
  pars_organized['distancing'] = r_distancing;
  pars_organized['covid_stay'] = r_covid_stay;
  pars_organized['chi'] = r_chi;
  pars_organized['adm_reduc'] = r_adm_reduc;
  pars_organized['comm_denom'] = r_comm_denom;
  
  ### return updated par set
  return(pars_organized)
}

### Pandemic parameter for intermediate response
f_params_intermediate = function(pars_in){
  
  ### pandemic impacts
  r_prophylaxis = runif(1, 0.68, 0.81);
  r_a_surge = 0;
  r_masks = 1; while(r_masks>=1){r_masks = rnorm(1,0.59,0.069) * runif(1,0.5,1)}
  r_handrub = runif(1, 0.15, 0.2);
  r_disorg = runif(1, 0.4, 0.6);
  r_distancing = runif(1, 0.4, 0.6);
  r_covid_stay = runif(1, 0.8, 1);
  r_chi = runif(1, 0, 0.5);
  r_adm_reduc = runif(1, 0.4, 0.6);
  r_comm_denom = runif(1, 0.4, 0.6);
  
  ### base
  pars_intermediate = pars_in;
  
  ### update parameters
  pars_intermediate['prophylaxis'] = r_prophylaxis;
  pars_intermediate['a_surge'] = r_a_surge;
  pars_intermediate['masks'] = r_masks;
  pars_intermediate['handrub'] = r_handrub;
  pars_intermediate['disorg'] = r_disorg;
  pars_intermediate['distancing'] = r_distancing;
  pars_intermediate['covid_stay'] = r_covid_stay;
  pars_intermediate['chi'] = r_chi;
  pars_intermediate['adm_reduc'] = r_adm_reduc;
  pars_intermediate['comm_denom'] = r_comm_denom;
  
  ### return updated par set
  return(pars_intermediate)
}

### Pandemic parameter for overwhelmed response
f_params_overwhelmed = function(pars_in){
  
  ### pandemic impacts
  r_prophylaxis = runif(1, 0.9, 1);
  r_a_surge = 0;
  r_masks = 0;
  r_handrub = runif(1, 0, 0.05);
  r_disorg = runif(1, 0.4, 0.6);
  r_distancing = runif(1, 0, 0.2);
  r_covid_stay = runif(1, 0.8, 1);
  r_chi = runif(1, 0, 0.5);
  r_adm_reduc = runif(1, 0.4, 0.6);
  r_comm_denom = runif(1, 0.4, 0.6);
  
  ### base
  pars_overwhelmed = pars_in;
  
  ### update parameters
  pars_overwhelmed['prophylaxis'] = r_prophylaxis;
  pars_overwhelmed['a_surge'] = r_a_surge;
  pars_overwhelmed['masks'] = r_masks;
  pars_overwhelmed['handrub'] = r_handrub;
  pars_overwhelmed['disorg'] = r_disorg;
  pars_overwhelmed['distancing'] = r_distancing;
  pars_overwhelmed['covid_stay'] = r_covid_stay;
  pars_overwhelmed['chi'] = r_chi;
  pars_overwhelmed['adm_reduc'] = r_adm_reduc;
  pars_overwhelmed['comm_denom'] = r_comm_denom;
  
  ### return updated par set
  return(pars_overwhelmed)
}

### Bacterial parameters for S. aureus
f_params_s_aureus = function(pars_in){
  
  ### bacteria parameters
  r_gamma_Cs = 1/rnorm(1, 287, 17.9)
  r_cost = rnorm(1, 0.2, 0.02)
  r_gamma_Cr = r_gamma_Cs*(1+r_cost)
  r_r_S = runif(1,0.172,0.489)
  r_r_R = runif(1,0.908,0.982)
  r_beta_Cs = rnorm(1, 0.057, 0.0057)
  r_beta_Cr = r_beta_Cs
  r_theta = 1/runif(1, 1, 10)
  r_alpha = -1; while(r_alpha<0){r_alpha=rnorm(1, 0.0016, 0.0008) * rcauchy(1, 2.97, 0.28)}
  prop_col = rnorm(1,0.0757, 0.00364)
  prop_r = rnorm(1, 0.16, 0.016)
  r_f_Cs = prop_col - (prop_col*prop_r)
  r_f_Cr = prop_col*prop_r
  r_f_U = 1 - r_f_Cs - r_f_Cr
  
  ### all other parameters use baseline values given in pars_states, i.e. "m_par"
  
  pars_s_aureus = pars_in;
  
  ### update parameters
  pars_s_aureus['gamma_Cs'] = r_gamma_Cs;
  pars_s_aureus['cost'] = r_cost;
  pars_s_aureus['gamma_Cr'] = r_gamma_Cr;
  pars_s_aureus['r_S'] = r_r_S;
  pars_s_aureus['r_R'] = r_r_R;
  pars_s_aureus['beta_Cs'] = r_beta_Cs
  pars_s_aureus['beta_Cr'] = r_beta_Cr
  pars_s_aureus['theta'] = r_theta;
  pars_s_aureus['alpha'] = r_alpha;
  pars_s_aureus['f_Cs'] = r_f_Cs;
  pars_s_aureus['f_Cr'] = r_f_Cr;
  pars_s_aureus['f_U'] = r_f_U;
  
  ### return updated par set
  return(pars_s_aureus)
}


### Bacterial parameters for E. coli
f_params_e_coli = function(pars_in){
  
  ### bacteria parameters
  r_gamma_Cs = rnorm(1, 0.00269, 0.000216)
  r_cost = 0
  r_gamma_Cr = r_gamma_Cs
  r_r_S = runif(1,0.096,0.365)
  r_r_R = runif(1,0.774,0.922)
  r_beta_Cs = -1; while(r_beta_Cs<0){r_beta_Cs=rnorm(1, 0.0078, 0.00334)}
  r_beta_Cr = r_beta_Cs
  r_theta = 1/runif(1, 1, 10)
  r_alpha = -1; while(r_alpha<0){r_alpha = rnorm(1, 0.0024, 0.000663) * rcauchy(1, 11.80, 0.80)}
  prop_col = rnorm(1,0.275, 0.0140)
  prop_r = -1; while(prop_r<0){prop_r = rnorm(1, 0.119, 0.0413)}
  r_f_Cs = prop_col - (prop_col*prop_r)
  r_f_Cr = prop_col*prop_r
  r_f_U = 1 - r_f_Cs - r_f_Cr
  
  ### all other parameters use baseline values given in pars_states, i.e. "m_par"
  
  pars_e_coli = pars_in;
  
  ### update parameters
  pars_e_coli['gamma_Cs'] = r_gamma_Cs;
  pars_e_coli['cost'] = r_cost;
  pars_e_coli['gamma_Cr'] = r_gamma_Cr;
  pars_e_coli['r_S'] = r_r_S;
  pars_e_coli['r_R'] = r_r_R;
  pars_e_coli['beta_Cs'] = r_beta_Cs
  pars_e_coli['beta_Cr'] = r_beta_Cr
  pars_e_coli['theta'] = r_theta;
  pars_e_coli['alpha'] = r_alpha;
  pars_e_coli['f_Cs'] = r_f_Cs;
  pars_e_coli['f_Cr'] = r_f_Cr;
  pars_e_coli['f_U'] = r_f_U;
  
  ### return updated par set
  return(pars_e_coli)
}





#########################################
### PARAMETER SET GENERATOR FUNCTIONS ###
#########################################

### ANALYSIS2 ###
## Generate n parameter sets
## not all uncertainty is always accounted for, so must provide the "base" par set "pars_m" as a scaffold

f_generate_parsets_uncertain_facility_bacteria = function(n_parsets, pars_m){
  
  mat_pars = matrix(data = NA, nrow = n_parsets, ncol = length(pars_m)+1, byrow = F)
  
  for(n in 1:n_parsets){
    if(n_parsets %% 20 ==0){print(paste0("generating parset ", n, " of ", n_parsets))}
    pars_n = c(parset = n, f_params_uncertain_facility_bacteria(pars_m))
    mat_pars[n,] = pars_n
  }
  colnames(mat_pars) = c('parset',names(pars_m))
  
  return(mat_pars)
}


### ANALYSIS3 ###
### CASE STUDIES
## Generate n parameter sets
## species parameters (transmission) depend on setting (contact rates) and should be matched to setting
## (e.g. any given sim x for E. coli and S. aureus in a particular setting should have same setting parameters)



### function 1: generate parameters for setting (on scaffold of pars_in)
f_generate_parsets_casestudies_setting = function(n_parsets, pars_in, setting){
  
  if(!setting %in% c("rehab", "geriatric", "paediatric")){warning("setting not available");return()}
  
  mat_pars_setting = matrix(data = NA, nrow = n_parsets, ncol = length(pars_in)+1, byrow = F)
  
  for(n in 1:n_parsets){
    if(n_parsets %% 20 ==0){print(paste0("generating parset ", n, " of ", n_parsets))}
    
    # SETTING: 
    if(setting == "rehab"){pars_setting = f_params_rehab(pars_in)}
    if(setting == "geriatric"){pars_setting = f_params_geriatric(pars_in)}
    if(setting == "paediatric"){pars_setting = f_params_paediatric(pars_in)}
    
    pars_out = c(parset = n, pars_setting)
    
    mat_pars_setting[n,] = pars_out
  }
  
  colnames(mat_pars_setting) = c('parset',names(pars_in))
  
  return(mat_pars_setting)
}

### function 2: generate parameters for scenario (on scaffold of pars_in)
f_generate_parsets_casestudies_scenario = function(n_parsets, pars_in, scenario){
  
  if(!scenario %in% c("organized", "intermediate", "overwhelmed")){warning("scenario not available");return()}
  
  mat_pars_scenario = matrix(data = NA, nrow = n_parsets, ncol = length(pars_in)+1, byrow = F)
  
  for(n in 1:n_parsets){
    if(n_parsets %% 20 ==0){print(paste0("generating parset ", n, " of ", n_parsets))}
    
    # SETTING: 
    if(scenario == "organized"){pars_scenario = f_params_organized(pars_in)}
    if(scenario == "intermediate"){pars_scenario = f_params_intermediate(pars_in)}
    if(scenario == "overwhelmed"){pars_scenario = f_params_overwhelmed(pars_in)}
    
    pars_out = c(parset = n, pars_scenario)
    
    mat_pars_scenario[n,] = pars_out
  }
  
  colnames(mat_pars_scenario) = c('parset',names(pars_in))
  
  return(mat_pars_scenario)
}

### function 3: generate parameters for species (on scaffold of pars_in)
# NB: because pi is calculated using both kappa and beta, this will have to be updated when species and setting parameters are combined
f_generate_parsets_casestudies_species = function(n_parsets, pars_in, species){
  
  if(!species %in% c("e_coli", "s_aureus")){warning("species not available");return()}
  
  mat_pars_species = matrix(data = NA, nrow = n_parsets, ncol = length(pars_in)+1, byrow = F)
  
  for(n in 1:n_parsets){
    if(n_parsets %% 20 ==0){print(paste0("generating parset ", n, " of ", n_parsets))}
    
    # SPECIES: update pars_in to reflect desired species
    if(species == "s_aureus"){pars_species = f_params_s_aureus(pars_in)}
    if(species == "e_coli"){pars_species = f_params_e_coli(pars_in)}
    
    pars_out = c(parset = n, pars_species)
    
    mat_pars_species[n,] = pars_out
  }
  
  colnames(mat_pars_species) = c('parset',names(pars_in))
  
  return(mat_pars_species)
}

### function 3: combine setting, scenario and species parameters
f_parsets_setting_scenario_species = function(pars_setting, pars_scenario, pars_species){
  
  colnames_remainder = c("parset", "upsilon", "eta",  "pi_S", "upsilon_SL", "prop_symptomatic")
  colnames_setting = c("mu", "Nbeds", "Nhcw", "kappa_pa_pa", "kappa_pa_pe", "kappa_pe_pe", "kappa_pe_pa", "rho", "hyg", "a", 
                       "t_policy", "t_burnin")
  colnames_scenario = c("chi", "disorg", "comm_denom", "a_surge", "covid_stay", "adm_reduc", "handrub", "masks", "prophylaxis", "distancing")
  colnames_species = c("beta_S", "f_U", "f_Cs", "f_Cr", "theta", "r_S", "r_R", "gamma_Cs", "cost", "gamma_Cr", "alpha", "beta_Cs", "beta_Cr",
                       "pi_Cs", "pi_Cr")
  
  if(nrow(pars_setting) != nrow(pars_scenario)){warning("nrows of pars_setting and pars_scenario don't match")}
  if(nrow(pars_setting) != nrow(pars_species)){warning("nrows of pars_setting and pars_species don't match")}
  
  if(ncol(pars_setting) != ncol(pars_scenario)){warning("ncols of pars_setting and pars_scenario don't match")}
  if(ncol(pars_setting) != ncol(pars_species)){warning("ncols of pars_setting and pars_species don't match")}
  
  
  mat_pars_setting_scenario_species = matrix(data = NA, nrow = nrow(pars_setting), ncol = ncol(pars_setting), byrow = F)
  
  for(n in 1:nrow(mat_pars_setting_scenario_species)){
    
    mat_pars_remainder_n = pars_species[n,colnames_remainder]
    mat_pars_setting_n = pars_setting[n,colnames_setting]
    mat_pars_scenario_n = pars_scenario[n,colnames_scenario]
    mat_pars_species_n = pars_species[n,colnames_species]
    
    ### combine parameters
    mat_pars_setting_scenario_species[n,] = c(mat_pars_remainder_n, mat_pars_setting_n, mat_pars_scenario_n, mat_pars_species_n)
    colnames(mat_pars_setting_scenario_species) = c(colnames_remainder, colnames_setting, colnames_scenario, colnames_species)
    
    ### update rates of transmission per contact
    # transmission rates per contact
    r_beta_S = as.numeric(mat_pars_species_n['beta_S'])
    r_beta_Cs = as.numeric(mat_pars_species_n['beta_Cs'])
    r_beta_Cr = as.numeric(mat_pars_species_n['beta_Cr'])
    
    # numbers of contacts
    r_kappa_pa_pa = as.numeric(mat_pars_setting_n['kappa_pa_pa'])
    r_kappa_pa_pe = as.numeric(mat_pars_setting_n['kappa_pa_pe'])
    r_kappa_pe_pe = as.numeric(mat_pars_setting_n['kappa_pe_pe'])
    r_kappa_pe_pa = as.numeric(mat_pars_setting_n['kappa_pe_pa'])
    
    ### transmission rates per contact
    
    r_pi_S = r_beta_S/(r_kappa_pa_pa+r_kappa_pa_pe+r_kappa_pe_pe+r_kappa_pe_pa)
    r_pi_Cs = r_beta_Cs/(r_kappa_pa_pa+r_kappa_pe_pa)
    r_pi_Cr = r_beta_Cr/(r_kappa_pa_pa+r_kappa_pe_pa)
    
    mat_pars_setting_scenario_species[n,'pi_S'] = r_pi_S
    mat_pars_setting_scenario_species[n,'pi_Cs'] = r_pi_Cs
    mat_pars_setting_scenario_species[n,'pi_Cr'] = r_pi_Cr
  }
  
  return(mat_pars_setting_scenario_species)
}



#####################################################
### Function to run ODEs and calculate indicators ###
#####################################################
# requires loading of pre-generated parameter sets

f_runODEs_calculateIndicators = function(model_loop,
                                         value_tau_loop,
                                         time_out_loop,
                                         n_parset_start,
                                         n_parset_end,
                                         filepath_parsets, 
                                         save_data = 'NO'){
  
  
  # load parameter sets
  m_pars_loop = loadRData(filepath_parsets)
  
  # create empty dataframe 
  df_loop_all = data.frame()
  
  # initiate loop
  for(n in n_parset_start:n_parset_end){
    
    print(paste0("on parset ",n, " of ", n_parset_end))
    
    # set pars and states for loop
    pars_loop = m_pars_loop[n,]
    pars_covid_loop = pars_loop
    states_loop = f_states(pars_loop)
    
    # run once without any tau
    ode_output_no_tau = f_dynamics_add_covid(ODE_model = model_loop, 
                                             states_base = states_loop, 
                                             pars_base = pars_loop, 
                                             pars_covid = pars_covid_loop,
                                             time_out = time_out_loop)
    
    df_indicators = f_indicators(n_sim = n, tau = "none", tau_val = value_tau_loop, ode_output_no_tau)
    
    qounter = 0
    for(impact in list_impacts){
      
      # set pars and states again for this loop (otherwise parameter updates accumulate)
      pars_loop = m_pars_loop[n,]
      pars_covid_loop = pars_loop
      states_loop = f_states(pars_loop)
      
      print(paste0("on impact ",impact))
      
      # update counter to implement new impact
      qounter = qounter+1
      impact_name = names(list_impacts)[qounter]
      
      # 
      for(par in impact){
        pars_covid_loop[par] = value_tau_loop
      }
      
      
      ### Run ODEs
      ode_output_tau = f_dynamics_add_covid(ODE_model = model_loop, 
                                            states_base = states_loop, 
                                            pars_base = pars_loop, 
                                            pars_covid = pars_covid_loop,
                                            time_out = time_out_loop)
      
      ### Calculate indicators and bind to dataframe 
      df_indicators_tau_temp = f_indicators(n_sim = n, tau = impact_name, tau_val = value_tau_loop, ode_output_tau)
      
      df_indicators = rbind(df_indicators, df_indicators_tau_temp)
    }
    
    #print(df_indicators)
    if(save_data == 'YES'){
      print(paste0("saving data for parset ", n, " of ", n_parset_end))
      save(df_indicators, file = paste0("data/simu_montecarlo_uncertainty_parset", n,".Rdata"))}else{
      return(df_indicators)
    }
  }
}

##############################################
### Function to calculate indicator deltas ###
##############################################
### Differences in indicators as resulting form pandemic impacts
### i.e. comparing simulations for each tau to baseline simulations with no tau

f_indicator_deltas = function(filepath_indicators_summarized, save_data = 'NO'){
  
  # initialize final df
  df_deltas = data.frame()
  
  # load df_indicators
  df_indicators_summarized_raw = loadRData(filepath_indicators_summarized)
  
  # for each simulation, separate simulations into baseline (data_sim_no_tau) vs pandemic impact (data_sim_tau)
  for(n_sim_i in 1:max(df_indicators_summarized_raw$n_sim)){
    
    print(paste0("on Monte Carlo simulation ",n_sim_i, " of ", max(df_indicators_summarized_raw$n_sim)))
    
    data_indicators_raw_i = df_indicators_summarized_raw%>%filter(n_sim == n_sim_i)
    
    cols_data_indicators_raw = colnames(data_indicators_raw_i)
    cols_indicators = cols_data_indicators_raw[!cols_data_indicators_raw %in% c("n_sim", "tau", "tau_val", "tau_group")]
    
    data_sim_no_tau = data_indicators_raw_i%>%filter(tau == 'none')
    data_sim_tau = data_indicators_raw_i%>%filter(tau != 'none')
    
    # isolate each pandemic impact (row_i) in data_tau
    for(row_i in 1:nrow(data_sim_tau)){
      
      data_sim_no_tau_i = data_sim_no_tau
      data_sim_tau_i = data_sim_tau[row_i,]
      
      # and for each indicator, calculate relative difference from baseline
      delta_wide = data_sim_tau_i[,cols_indicators]/data_sim_no_tau_i[,cols_indicators]
      
      df_deltas_i = data.frame(n_sim = n_sim_i,
                               tau = data_sim_tau$tau[row_i],
                               tau_val = data_sim_tau$tau_val[row_i],
                               delta_wide)
      
      df_deltas = rbind(df_deltas, df_deltas_i)
      
    }
  }
  
  if(save_data == 'YES'){
    save(df_deltas, file = paste0("data/indicators_summarized_deltas.Rdata"))
  }else{
    return(df_deltas)
  }
}

####################################
### Function to bootstrap deltas ###
####################################

f_bootstrap_deltas = function(n_bootstraps, save_data = 'NO'){
  
  # load deltas
  data_deltas_raw = loadRData("data/outputs_analysis2/indicators_summarized_deltas.Rdata")
  
  # initialize data.frame
  df_bootstrap = data.frame()
  
  # loop through an increasing # of simulations to resample from delta outcomes, 
  for(sims_i in c(5, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500)){
    
    sims_sampled_i = sample(as.numeric(levels(factor(data_deltas_raw$n_sim))), sims_i)
    
    for(tau_i in levels(factor(data_deltas_raw$tau))){
      print(tau_i)
      # select data
      dat_i = data_deltas_raw%>%filter(n_sim %in% sims_sampled_i, tau == tau_i)
      
      # empty vectors to store bootstrap means and variances
      mean_V_inc_cumul_all = c()
      var_V_inc_cumul_all = c()
      
      mean_Cr_inc_cumul_all = c()
      var_Cr_inc_cumul_all = c()
      
      # empty data.frame to bind results for each tau
      
      df_bootstrap_i = c()
      
      for(bootstrap_i in 1:n_bootstraps){
        
        # calculate boostrap mean and boostrap variance for SARS-CoV-2 incidence and bacterial colonization incidence
        sample_i_V_inc_cumul_all = sample(dat_i$V_inc_cumul_all, replace = T)
        sample_i_Cr_inc_cumul_all = sample(dat_i$Cr_inc_cumul_all, replace = T)
        
        mean_i_V_inc_cumul_all = mean(sample_i_V_inc_cumul_all)
        var_i_V_inc_cumul_all = var(sample_i_V_inc_cumul_all)
        
        mean_i_Cr_inc_cumul_all = mean(sample_i_Cr_inc_cumul_all)
        var_i_Cr_inc_cumul_all = var(sample_i_Cr_inc_cumul_all)
        
        # append vectors to build distribution
        mean_V_inc_cumul_all = c(mean_V_inc_cumul_all, mean_i_V_inc_cumul_all)
        var_V_inc_cumul_all = c(var_V_inc_cumul_all, var_i_V_inc_cumul_all)
        
        mean_Cr_inc_cumul_all = c(mean_Cr_inc_cumul_all, mean_i_Cr_inc_cumul_all)
        var_Cr_inc_cumul_all = c(var_Cr_inc_cumul_all, var_i_Cr_inc_cumul_all)
      }
      
      df_bootstrap_i = data.frame(n_sims = sims_i,
                                  tau = tau_i,
                                  mean_V_inc_cumul_all = mean_V_inc_cumul_all,
                                  var_V_inc_cumul_all = var_V_inc_cumul_all,
                                  mean_Cr_inc_cumul_all = mean_Cr_inc_cumul_all,
                                  var_Cr_inc_cumul_all = var_Cr_inc_cumul_all)
      
      df_bootstrap = rbind(df_bootstrap, df_bootstrap_i)
    }
  }
  
  if(save_data == 'YES'){
    save(df_bootstrap, file = paste0("data/indicators_deltas_bootstrapped.Rdata"))
  }else{
    return(df_bootstrap)
  }
}


#########################################################
### Function to calculate PRCC for analysis2 outcomes ###
#########################################################

### PRCC for key outcomes
f_prcc = function(m_pars, df_outputs, vec_pars_varied, vec_tau, vec_outcomes, save.data = "NO"){
  
  # initialize output df
  df_prcc_out = data.frame()
  
  # turn par matrix into par df
  df_pars = m_pars%>%as.data.frame()%>%
    mutate(n_sim = parset)%>%
    dplyr::select(n_sim, vec_pars_varied)
  
  # loop through selected tau and outcomes and conduct distinct PRCC for each
  for(tau_i in vec_tau){
    for(outcome_i in vec_outcomes){
      
      # combine input parameters and output data
      df_prcc_i = left_join(df_pars, 
                             df_outputs%>%
                              dplyr::filter(tau == tau_i)%>%
                              dplyr::select(n_sim, outcome_i), 
                             by = "n_sim")%>%
        dplyr::select(-n_sim)

      # calculate PRCC
      df_prcc_out_i = epi.prcc(df_prcc_i, sided.test = 2)%>%
        mutate(tau = tau_i,
               outcome = outcome_i)
      df_prcc_out_i['par_name'] = vec_pars_varied
      
      df_prcc_out = rbind(df_prcc_out, df_prcc_out_i)
    }
  }
  
  if(save.data == "YES"){
    save(df_prcc_out, file = paste0("data/data_prcc.Rdata"))
  }else{
    return(df_prcc_out)
  }
}


####################################
### Function to run case studies ###
####################################


f_case_studies_dynamics = function(model_loop, 
                                   time_out_loop,
                                   setting, 
                                   scenario,
                                   species,
                                   in_beta_S,
                                   in_t_policy,
                                   n_parset_start, 
                                   n_parset_end, 
                                   save_data = 'NO'){
  
  if(!setting %in% c("rehab", "geriatric", "paediatric")){warning("setting not available");return()}
  if(!scenario %in% c("organized", "intermediate", "overwhelmed")){warning("scenario not available");return()}
  if(!species %in% c("e_coli", "s_aureus")){warning("species not available");return()}
  
  vec_setting_scenario_species = paste0(setting, "_", scenario, "_", species)
  
  ### load parameter set
  m_pars_loop = loadRData(paste0("data/outputs_analysis3/", vec_setting_scenario_species, "/parsets_analysis3_", vec_setting_scenario_species, ".Rdata"))
  
  # create empty dataframe 
  df_loop_indicators_nocovid = data.frame()
  df_loop_indicators_covid = data.frame()
  df_loop_metadata_nocovid = data.frame()
  df_loop_metadata_covid = data.frame()
  df_loop_indicators_deltas = data.frame()
  df_loop_metadata_deltas = data.frame()
  
  for(n in n_parset_start:n_parset_end){
    
    print(paste0("on parset ",n, " of ", n_parset_end))
    
    ### set pars and states for loop
    pars_loop = m_pars_loop[n,]
    states_loop = f_states(pars_loop)
    
    ### loop part 1: 180 days without COVID-19
    # NB: still use 'add_covid' function but just use pre-pandemic parameters i.e. 
    
    # assign pars
    pars_loop_nocovid = pars_loop
    
    # no COVID-19 transmission
    pars_loop_nocovid["pi_S"] = 0
    
    # no COVID-19 response parameters
    pars_loop_nocovid['prophylaxis'] = 0;
    pars_loop_nocovid['a_surge'] = 0;
    pars_loop_nocovid['masks'] = 0;
    pars_loop_nocovid['handrub'] = 0;
    pars_loop_nocovid['disorg'] = 0;
    pars_loop_nocovid['distancing'] = 0;
    pars_loop_nocovid['covid_stay'] = 0;
    pars_loop_nocovid['chi'] = 0;
    pars_loop_nocovid['adm_reduc'] = 0;
    pars_loop_nocovid['comm_denom'] = 0;
    
    
    # run dynamics
    ode_output_nocovid = f_dynamics_add_covid(ODE_model = model_loop, 
                                              states_base = states_loop,
                                              pars_base = pars_loop_nocovid,
                                              pars_covid = pars_loop_nocovid,
                                              time_out = time_out_loop)
    
    indicators_nocovid = f_indicators(n_sim = n, 
                                      tau = vec_setting_scenario_species, 
                                      tau_val = 0, 
                                      ode_output_nocovid)
    
    metadata_nocovid0 = f_dynamics_metadata(ode_output_nocovid%>%mutate(n_sim = n),
                                           par = vec_setting_scenario_species,
                                           par_val = 0)%>%
      filter(time>0)%>%
      group_by(par)%>%
      summarise(across(where(is.numeric), sum))%>%
      mutate(n_sim = n)%>%
      relocate(n_sim, par, par_val)
    
    metadata_nocovid = as.data.frame(c(metadata_nocovid0[1:4], as.data.frame(do.call(cbind,lapply(metadata_nocovid0[5:ncol(metadata_nocovid0)], function(x) x/time_out_loop)))))
    
    
    ### loop part 2: 180 days with COVID-19
    
    # assign pars (including COVID-19 response parameters)
    pars_loop_covid = pars_loop
    
    # calculate updated pi_S based on given contact rates and particular input beta_S
    pi_S_i = as.numeric(in_beta_S/(pars_loop_covid['kappa_pa_pa']+pars_loop_covid['kappa_pa_pe']+pars_loop_covid['kappa_pe_pe']+pars_loop_covid['kappa_pe_pa']))
    
    # update SARS-CoV-2 transmission rate and t_policy
    pars_loop_covid['pi_S'] = pi_S_i
    pars_loop_covid['t_policy'] = in_t_policy
    
    # run dynamics
    ode_output_covid = f_dynamics_add_covid(ODE_model = model_loop, 
                                              states_base = states_loop,
                                              pars_base = pars_loop_nocovid,
                                              pars_covid = pars_loop_covid,
                                              time_out = time_out_loop)
    
    indicators_covid = f_indicators(n_sim = n, 
                                      tau = vec_setting_scenario_species, 
                                      tau_val = 0, 
                                      ode_output_covid)
    
    metadata_covid0 = f_dynamics_metadata(ode_output_covid%>%mutate(n_sim = n),
                                           par = vec_setting_scenario_species,
                                           par_val = 0)%>%
      filter(time>0)%>%
      group_by(par)%>%
      summarise(across(where(is.numeric), sum))%>%
      mutate(n_sim = n)%>%
      relocate(n_sim, par, par_val)
    
    metadata_covid = as.data.frame(c(metadata_covid0[1:4], as.data.frame(do.call(cbind,lapply(metadata_covid0[5:ncol(metadata_covid0)], function(x) x/time_out_loop)))))
    
    ### loop part 3: find differences between simulations with COVID-19 and without
    
    # indicators
    indicators_deltas = merge(indicators_covid[1:3], indicators_covid[4:length(indicators_covid)]/indicators_nocovid[4:length(indicators_nocovid)])
    
    # metadata
    metadata_deltas = merge(metadata_covid[1:3], metadata_covid[4:length(metadata_covid)]/metadata_nocovid[4:length(metadata_nocovid)])
    
    
    ### bind output data
    df_loop_indicators_nocovid = rbind(df_loop_indicators_nocovid, indicators_nocovid)
    df_loop_indicators_covid = rbind(df_loop_indicators_covid, indicators_covid)
    df_loop_indicators_deltas = rbind(df_loop_indicators_deltas, indicators_deltas)
    df_loop_metadata_nocovid = rbind(df_loop_metadata_nocovid, metadata_nocovid)
    df_loop_metadata_covid = rbind(df_loop_metadata_covid, metadata_covid)
    df_loop_metadata_deltas = rbind(df_loop_metadata_deltas, metadata_deltas)
  }
  
  # update outputs to include a column for variable pars (t_policy and beta_S)
  df_loop_indicators_nocovid = df_loop_indicators_nocovid%>%mutate(t_policy = in_t_policy,beta_S = in_beta_S)
  df_loop_indicators_covid = df_loop_indicators_covid%>%mutate(t_policy = in_t_policy,beta_S = in_beta_S)
  df_loop_indicators_deltas = df_loop_indicators_deltas%>%mutate(t_policy = in_t_policy,beta_S = in_beta_S)
  df_loop_metadata_nocovid = df_loop_metadata_nocovid%>%mutate(t_policy = in_t_policy,beta_S = in_beta_S)
  df_loop_metadata_covid = df_loop_metadata_covid%>%mutate(t_policy = in_t_policy,beta_S = in_beta_S)
  df_loop_metadata_deltas = df_loop_metadata_deltas%>%mutate(t_policy = in_t_policy,beta_S = in_beta_S)
  
  ### save outputs
  if(save_data == 'YES'){
    print(paste0('saving case study data for ', vec_setting_scenario_species))
    # don't save direct to raw data folder : save to initial greater folder so as not to accidentally overwrite
    save(df_loop_indicators_nocovid, file = paste0("data/outputs_analysis3/", vec_setting_scenario_species, "_indicators_nocovid_t_policy_", in_t_policy, "_beta_S_", in_beta_S, ".Rdata"))
    save(df_loop_indicators_covid, file = paste0("data/outputs_analysis3/", vec_setting_scenario_species, "_indicators_covid_t_policy_", in_t_policy, "_beta_S_", in_beta_S, ".Rdata"))
    save(df_loop_indicators_deltas, file = paste0("data/outputs_analysis3/", vec_setting_scenario_species, "_indicators_deltas_t_policy_", in_t_policy, "_beta_S_", in_beta_S, ".Rdata"))
    save(df_loop_metadata_nocovid, file = paste0("data/outputs_analysis3/", vec_setting_scenario_species, "_metadata_nocovid_t_policy_", in_t_policy, "_beta_S_", in_beta_S, ".Rdata"))
    save(df_loop_metadata_covid, file = paste0("data/outputs_analysis3/", vec_setting_scenario_species, "_metadata_covid_t_policy_", in_t_policy, "_beta_S_", in_beta_S, ".Rdata"))
    save(df_loop_metadata_deltas, file = paste0("data/outputs_analysis3/", vec_setting_scenario_species, "_metadata_deltas_t_policy_", in_t_policy, "_beta_S_", in_beta_S, ".Rdata"))
  }else{
    return(df_loop_indicators_deltas)
  }
}
