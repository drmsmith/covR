###################################
### Hospital coinfection models ###
###################################

### SARS-CoV-2 notation
# S = susceptible
# E = exposed
# I = infectious 
# R = recovered

### ARB notation
# U = uncolonized
# Cs = colonized sensitive 
# Cr = colonized resistant 

### NB: D/Ts/Tr is used for indicators and in manuscript for transient bacterial carriage among staff
###     U/Cs/Cr is used in ODEs for both patients and staff

ODEs_covR <- function(t, states, params){
  
  ### write derivatives
  with(as.list(c(states, params)),{
    
    ### DEFINE CUT-OFFS FOR POLICY CHANGE PARAMETERS
    # if there is no COVID (e.g. when finding equilibrium for MRB), none are in place
    # else if there is COVID:
    ### if t<t_policy, then policy parameters at zero
    ### else if t>= t_policy, then update policy parameters
    if(E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa < 0.001){
      ipc_C = 0;
      ipc_S = 0;
      a_proph = 0;
      lockd = 0
    }else{if(t<t_policy){
      ipc_C = 0;
      ipc_S = 0;
      a_proph = 0;
      lockd = 0
    }else{if((t - t_policy) < t_burnin){
      ipc_C = handrub*(t - t_policy)/t_burnin;
      ipc_S = masks*(t - t_policy)/t_burnin;
      a_proph = prophylaxis*(t - t_policy)/t_burnin;
      lockd = distancing*(t - t_policy)/t_burnin;
    }else{
      ipc_C = handrub;
      ipc_S = masks;
      a_proph = prophylaxis;
      lockd = distancing}
    }
    }
    
    
    # ### Introduce a new case weekly (1 new HCW infected weekly)
    # if(E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa > 0.001){
    #   if(t > 1 & t %% 7 == 0){
    #     if(S_U_pe > R_U_pe){
    #       S_U_pe = S_U_pe - 1
    #       E_U_pe = E_U_pe + 1
    #     }else{
    #       R_U_pe = R_U_pe - 1
    #       E_U_pe = E_U_pe + 1
    #     }
    #   }
    # }
    
    # ### DEFINE CUT-OFFS FOR SARS-CoV-2 outbreak
    # # if new infections too few (eg with all COVID-19 responses), solver stops before t_out (state variables too miniscule)
    # # so if the number of recovered > 1 while infection prevalence very small
    if(R_U_pa + R_Cs_pa + R_Cr_pa + R_U_pe + R_Cs_pe + R_Cr_pe > 0.8 &
       E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + E_U_pe + E_Cs_pe + E_Cr_pe + I_U_pe + I_Cs_pe + I_Cr_pe < 0.2){
      E_U_pa = 0;
      E_Cs_pa = 0;
      E_Cr_pa = 0;
      I_U_pa = 0;
      I_Cs_pa = 0;
      I_Cr_pa = 0;
      E_U_pe = 0;
      E_Cs_pe = 0;
      E_Cr_pe = 0;
      I_U_pe = 0;
      I_Cs_pe = 0;
      I_Cr_pe = 0;
    }
    
    ### DEFINE COMPOUND PARAMETERS
    
    ### sums and proportions of compartments
    # total pa/pe
    N_pa <- S_U_pa + S_Cs_pa + S_Cr_pa + E_U_pa + E_Cs_pa + E_Cr_pa + I_U_pa + I_Cs_pa + I_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa
    N_pe <- S_U_pe + S_Cs_pe + S_Cr_pe + E_U_pe + E_Cs_pe + E_Cr_pe + I_U_pe + I_Cs_pe + I_Cr_pe + R_U_pe + R_Cs_pe + R_Cr_pe
    
    # patient infectious (number, frequency)
    N_I_pa = I_U_pa + I_Cs_pa + I_Cr_pa
    freq_I_pa = N_I_pa*(1/N_pa)
    
    # staff infectious (number, frequency)
    N_I_pe = I_U_pe + I_Cs_pe + I_Cr_pe
    freq_I_pe = N_I_pe*(1/N_pe)
    
    # patient colonized Cs (number, frequency)
    N_Cs_pa = S_Cs_pa+E_Cs_pa+I_Cs_pa+R_Cs_pa
    freq_Cs_pa = N_Cs_pa*(1/N_pa)
    
    # patient colonized Cr (number, frequency)
    N_Cr_pa = S_Cr_pa+E_Cr_pa+I_Cr_pa+R_Cr_pa
    freq_Cr_pa = N_Cr_pa*(1/N_pa)
    
    # staff transiently carrying Cs (number, frequency)
    N_Cs_pe = S_Cs_pe+E_Cs_pe+I_Cs_pe+R_Cs_pe
    freq_Cs_pe = N_Cs_pe*(1/N_pe)
    
    # staff transiently carrying Cr (number, frequency)
    N_Cr_pe = S_Cr_pe+E_Cr_pe+I_Cr_pe+R_Cr_pe
    freq_Cr_pe = N_Cr_pe*(1/N_pe)
    
    ### contact rates, including dynamic pandemic impacts
    
    # dynamic disorganization
    dyn_disorg = qbeta(abs(freq_I_pa), shape1 = 2, shape2 = 1/disorg)
    #print(dyn_disorg)
    
    # dynamic contact rates (k instead of kappa); calculate k_pe_pa after adjusting k_pa_pe:
    k_pa_pa = kappa_pa_pa*(1-lockd)
    k_pa_pe = kappa_pa_pe*(1+dyn_disorg*2)
    k_pe_pa = k_pa_pe*(N_pa/N_pe)
    k_pe_pe = kappa_pe_pe
    
    ### hand hygiene
    hyg_handrub = (hyg + ipc_C*(1-hyg))
    omega = (hyg_handrub * k_pe_pa * N_pe / N_pa)
    
    ### antibiotic prescribing
    # dynamic antibiotic prescribing, for infectious vs non-infectious patients
    A_base = a
    A_surge = (1-A_base)*qbeta(abs(freq_I_pa), shape1 = 2, shape2 = 1/a_surge)
    A_proph = a_proph*prop_symptomatic
    
    a_nonI = A_base+A_surge
    a_I = A_base+A_surge+(1-A_base-A_surge)*A_proph
    
    # print(paste0("A_base is ", A_base, ", A_surge is ", A_surge, "and A_proph is ", A_proph, " so A_nonI is ", a_nonI, " and A_I is ", a_I))
    
    ### proportion on sick-leave
    chi_eff = chi*prop_symptomatic
    
    ### forces of infection
    # Cs, patient -> patient
    foi_Cs_pa_pa = k_pa_pa*(1/rho)*pi_Cs*(freq_Cs_pa)
    # Cs, patient -> staff
    foi_Cs_pa_pe = k_pe_pa*rho*pi_Cs*(freq_Cs_pa)
    # Cs, staff -> patient
    foi_Cs_pe_pa = k_pa_pe*rho*pi_Cs*(freq_Cs_pe)
    # Cs, staff -> staff
    foi_Cs_pe_pe = k_pe_pe*(1/rho)*pi_Cs*(freq_Cs_pe)
    
    # Cr, patient -> patient
    foi_Cr_pa_pa = k_pa_pa*(1/rho)*pi_Cr*(freq_Cr_pa)
    # Cr, patient -> staff
    foi_Cr_pa_pe = k_pe_pa*rho*pi_Cr*(freq_Cr_pa)
    # Cr, staff -> patient
    foi_Cr_pe_pa = k_pa_pe*rho*pi_Cr*(freq_Cr_pe)
    # Cr, staff -> staff
    foi_Cr_pe_pe = k_pe_pe*(1/rho)*pi_Cr*(freq_Cr_pe)
    
    # S, patient -> patient
    foi_S_pa_pa <- k_pa_pa*pi_S*(1-ipc_S)*(freq_I_pa)
    # S, patient -> staff
    foi_S_pa_pe <- k_pe_pa*pi_S*(1-ipc_S)*(freq_I_pa)
    # S, staff -> patient
    foi_S_pe_pa <- k_pa_pe*pi_S*(1-ipc_S)*(freq_I_pe)
    # S, staff -> staff
    foi_S_pe_pe <- k_pe_pe*pi_S*(1-ipc_S)*(freq_I_pe)
    
    ### entry fractions: 
    # impact of comm_denom
    
    dyn_admission = mu*Nbeds*(1-qbeta(abs(freq_I_pa), shape1 = 2, shape2 = 1/adm_reduc))
    dyn_comm_denom = qbeta(abs(freq_I_pa), shape1 = 2, shape2 = 1/comm_denom)
    
    #print(dyn_comm_denom)
    
    ### DERIVATIVES: PATIENTS
    
    d_S_U_pa <- dyn_admission*(f_U - (f_Cr)*dyn_comm_denom) - 
      S_U_pa*(a_nonI)*(alpha*r_R) - 
      S_U_pa*(foi_Cs_pa_pa + foi_Cr_pa_pa) -
      S_U_pa*(foi_Cs_pe_pa + foi_Cr_pe_pa) + 
      S_Cs_pa*(gamma_Cs + (a_nonI)*theta*(1-r_S)) + 
      S_Cr_pa*(gamma_Cr + (a_nonI)*theta*(1-r_R)) - 
      S_U_pa*(foi_S_pa_pa) - 
      S_U_pa*(foi_S_pe_pa) -
      S_U_pa*mu
    
    d_S_Cs_pa <-  dyn_admission*f_Cs +  #dyn_admission*(f_Cs + (f_U*dyn_comm_denom*(f_Cs/(f_Cs+f_Cr)))) + 
      S_U_pa*(foi_Cs_pa_pa) + 
      S_U_pa*(foi_Cs_pe_pa) - 
      S_Cs_pa*(gamma_Cs + (a_nonI)*theta*(1-r_S)) - 
      S_Cs_pa * (foi_S_pa_pa) -
      S_Cs_pa * (foi_S_pe_pa) - 
      S_Cs_pa*mu
    
    d_S_Cr_pa <- dyn_admission*(f_Cr + (f_Cr)*dyn_comm_denom) +   #dyn_admission*(f_Cr + (f_U*dyn_comm_denom*(f_Cr/(f_Cs+f_Cr)))) + 
      S_U_pa*(a_nonI)*(alpha*r_R) + 
      S_U_pa*(foi_Cr_pa_pa) +
      S_U_pa*(foi_Cr_pe_pa) - 
      S_Cr_pa*(gamma_Cr + (a_nonI)*theta*(1-r_R)) - 
      S_Cr_pa * (foi_S_pa_pa) - 
      S_Cr_pa * (foi_S_pe_pa) -
      S_Cr_pa*mu
    
    d_E_U_pa <- - E_U_pa*(a_nonI)*(alpha*r_R) - 
      E_U_pa*(foi_Cs_pa_pa + foi_Cr_pa_pa) -
      E_U_pa*(foi_Cs_pe_pa + foi_Cr_pe_pa) + 
      E_Cs_pa*(gamma_Cs + (a_nonI)*theta*(1-r_S)) + 
      E_Cr_pa*(gamma_Cr + (a_nonI)*theta*(1-r_R)) + 
      S_U_pa*(foi_S_pa_pa) +
      S_U_pa*(foi_S_pe_pa) - 
      E_U_pa*eta - 
      E_U_pa*mu
    
    d_E_Cs_pa <- E_U_pa*(foi_Cs_pa_pa) +
      E_U_pa*(foi_Cs_pe_pa) - 
      E_Cs_pa*(gamma_Cs + (a_nonI)*theta*(1-r_S)) + 
      S_Cs_pa * (foi_S_pa_pa) +
      S_Cs_pa * (foi_S_pe_pa) - 
      E_Cs_pa*eta - 
      E_Cs_pa*mu
    
    d_E_Cr_pa <- E_U_pa*(a_nonI)*(alpha*r_R) + 
      E_U_pa*(foi_Cr_pa_pa) +
      E_U_pa*(foi_Cr_pe_pa) - 
      E_Cr_pa*(gamma_Cr + (a_nonI)*theta*(1-r_R)) + 
      S_Cr_pa * (foi_S_pa_pa) +
      S_Cr_pa * (foi_S_pe_pa) - 
      E_Cr_pa*eta - 
      E_Cr_pa*mu
    
    d_I_U_pa <- E_U_pa*eta - 
      I_U_pa*(a_I)*(alpha*r_R) - 
      I_U_pa*(foi_Cs_pa_pa + foi_Cr_pa_pa) -
      I_U_pa*(foi_Cs_pe_pa + foi_Cr_pe_pa) + 
      I_Cs_pa*(gamma_Cs + (a_I)*theta*(1-r_S)) + 
      I_Cr_pa*(gamma_Cr + (a_I)*theta*(1-r_R)) - I_U_pa*(upsilon) - 
      I_U_pa*mu*(1-covid_stay*prop_symptomatic)
    
    d_I_Cs_pa <- E_Cs_pa*eta + 
      I_U_pa*(foi_Cs_pa_pa) +
      I_U_pa*(foi_Cs_pe_pa) - 
      I_Cs_pa*(gamma_Cs + (a_I)*theta*(1-r_S)) - 
      I_Cs_pa*(upsilon) - 
      I_Cs_pa*mu*(1-covid_stay*prop_symptomatic)
    
    d_I_Cr_pa <- E_Cr_pa*eta + 
      I_U_pa*(a_I)*(alpha*r_R) + 
      I_U_pa*(foi_Cr_pa_pa) +
      I_U_pa*(foi_Cr_pe_pa) - 
      I_Cr_pa*(gamma_Cr + (a_I)*theta*(1-r_R)) - 
      I_Cr_pa*(upsilon) - 
      I_Cr_pa*mu*(1-covid_stay*prop_symptomatic)
    
    d_R_U_pa <- - R_U_pa*(a_nonI)*(alpha*r_R) - 
      R_U_pa*(foi_Cs_pa_pa + foi_Cr_pa_pa) -
      R_U_pa*(foi_Cs_pe_pa + foi_Cr_pe_pa) + 
      R_Cs_pa*(gamma_Cs + (a_nonI)*theta*(1-r_S)) + 
      R_Cr_pa*(gamma_Cr + (a_nonI)*theta*(1-r_R)) + 
      I_U_pa*(upsilon) - 
      R_U_pa*mu
    
    d_R_Cs_pa <- R_U_pa*(foi_Cs_pa_pa) +
      R_U_pa*(foi_Cs_pe_pa) - 
      R_Cs_pa*(gamma_Cs + (a_nonI)*theta*(1-r_S)) + 
      I_Cs_pa*(upsilon) - 
      R_Cs_pa*mu
    
    d_R_Cr_pa <- R_U_pa*(a_nonI)*(alpha*r_R) + 
      R_U_pa*(foi_Cr_pa_pa) +
      R_U_pa*(foi_Cr_pe_pa) - 
      R_Cr_pa*(gamma_Cr + (a_nonI)*theta*(1-r_R)) + 
      I_Cr_pa*(upsilon) - 
      R_Cr_pa*mu
    
    ### DERIVATIVES: PERSONNEL (STAFF)
    d_S_U_pe <- - S_U_pe*(foi_Cs_pe_pe + foi_Cr_pe_pe) -
      S_U_pe*(foi_Cs_pa_pe + foi_Cr_pa_pe) + 
      S_Cs_pe*(omega) + 
      S_Cr_pe*(omega) - 
      S_U_pe*(foi_S_pe_pe) - 
      S_U_pe*(foi_S_pa_pe)
    
    d_S_Cs_pe <- S_U_pe*(foi_Cs_pe_pe) +
      S_U_pe*(foi_Cs_pa_pe) -
      S_Cs_pe*(omega) - 
      S_Cs_pe*(foi_S_pe_pe) - 
      S_Cs_pe*(foi_S_pa_pe)
    
    d_S_Cr_pe <- S_U_pe*(foi_Cr_pe_pe) +
      S_U_pe*(foi_Cr_pa_pe) - 
      S_Cr_pe*(omega) - 
      S_Cr_pe*(foi_S_pe_pe) - 
      S_Cr_pe*(foi_S_pa_pe)
    
    d_E_U_pe <- - E_U_pe*(foi_Cs_pe_pe + foi_Cr_pe_pe) -
      E_U_pe*(foi_Cs_pa_pe + foi_Cr_pa_pe) +
      E_Cs_pe*(omega) + 
      E_Cr_pe*(omega) +
      S_U_pe*(foi_S_pe_pe) + 
      S_U_pe*(foi_S_pa_pe) - 
      E_U_pe*eta
    
    d_E_Cs_pe <- E_U_pe*(foi_Cs_pe_pe) +
      E_U_pe*(foi_Cs_pa_pe) -
      E_Cs_pe*(omega) + 
      S_Cs_pe*(foi_S_pe_pe) +
      S_Cs_pe*(foi_S_pa_pe) - 
      E_Cs_pe*eta
    
    d_E_Cr_pe <- E_U_pe*(foi_Cr_pe_pe) +
      E_U_pe*(foi_Cr_pa_pe) - 
      E_Cr_pe*(omega) + 
      S_Cr_pe*(foi_S_pe_pe) +
      S_Cr_pe*(foi_S_pa_pe) -
      E_Cr_pe*eta
    
    d_I_U_pe <- - I_U_pe*(foi_Cs_pe_pe + foi_Cr_pe_pe) -
      I_U_pe*(foi_Cs_pa_pe + foi_Cr_pa_pe) +
      I_Cs_pe*(omega) + 
      I_Cr_pe*(omega) +
      E_U_pe*eta - 
      I_U_pe*upsilon - 
      I_U_pe*chi_eff
    
    d_I_Cs_pe <- I_U_pe*(foi_Cs_pe_pe) + 
      I_U_pe*(foi_Cs_pa_pe) -
      I_Cs_pe*(omega) +
      E_Cs_pe*eta - 
      I_Cs_pe*upsilon - 
      I_Cs_pe*chi_eff
    
    d_I_Cr_pe <- I_U_pe*(foi_Cr_pe_pe) +
      I_U_pe*(foi_Cr_pa_pe) -
      I_Cr_pe*(omega) +
      E_Cr_pe*eta - 
      I_Cr_pe*upsilon - 
      I_Cr_pe*chi_eff
    
    d_R_U_pe <- - R_U_pe*(foi_Cs_pe_pe + foi_Cr_pe_pe) -
      R_U_pe*(foi_Cs_pa_pe + foi_Cr_pa_pe) + 
      R_Cs_pe*(omega) + 
      R_Cr_pe*(omega) +
      I_U_pe*upsilon + 
      SL*(upsilon_SL)
    
    d_R_Cs_pe <- R_U_pe*(foi_Cs_pe_pe) +
      R_U_pe*(foi_Cs_pa_pe) -
      R_Cs_pe*(omega) + 
      I_Cs_pe*upsilon
    
    d_R_Cr_pe <- R_U_pe*(foi_Cr_pe_pe) +
      R_U_pe*(foi_Cr_pa_pe) -
      R_Cr_pe*(omega) +
      I_Cr_pe*upsilon
    
    d_SL = I_U_pe*chi_eff + 
      I_Cs_pe*chi_eff + 
      I_Cr_pe*chi_eff - 
      SL*(upsilon_SL)
    
    d_incS_pa_pa <- (S_U_pa+S_Cs_pa+S_Cr_pa)*(foi_S_pa_pa)
    d_incS_pe_pa <- (S_U_pa+S_Cs_pa+S_Cr_pa)*(foi_S_pe_pa)
    d_incS_pa_pe <- (S_U_pe+S_Cs_pe+S_Cr_pe)*(foi_S_pa_pe)
    d_incS_pe_pe <- (S_U_pe+S_Cs_pe+S_Cr_pe)*(foi_S_pe_pe)
    d_incCs_pa_pa <- (S_U_pa+E_U_pa+I_U_pa+R_U_pa)*(foi_Cs_pa_pa)
    d_incCs_pe_pa <- (S_U_pa+E_U_pa+I_U_pa+R_U_pa)*(foi_Cs_pe_pa)
    d_incCs_pa_pe <- (S_U_pe+E_U_pe+I_U_pe+R_U_pe)*(foi_Cs_pa_pe)
    d_incCs_pe_pe <- (S_U_pe+E_U_pe+I_U_pe+R_U_pe)*(foi_Cs_pe_pe)
    d_incCr_pa_pa <- (S_U_pa+E_U_pa+I_U_pa+R_U_pa)*(foi_Cr_pa_pa)
    d_incCr_pe_pa <- (S_U_pa+E_U_pa+I_U_pa+R_U_pa)*(foi_Cr_pe_pa)
    d_incCr_pa_pe <- (S_U_pe+E_U_pe+I_U_pe+R_U_pe)*(foi_Cr_pa_pe)
    d_incCr_pe_pe <- (S_U_pe+E_U_pe+I_U_pe+R_U_pe)*(foi_Cr_pe_pe)
    d_incCr_endog <-  (S_U_pa+E_U_pa+I_U_pa+R_U_pa)*(a_nonI)*(alpha*r_R)
    
    d_abx = a_nonI*(S_U_pa + S_Cs_pa + S_Cr_pa + E_U_pa + E_Cs_pa + E_Cr_pa + R_U_pa + R_Cs_pa + R_Cr_pa) + a_I*(I_U_pa + I_Cs_pa + I_Cr_pa)
    d_kappa_patients = (k_pa_pa + k_pa_pe)
    d_kappa_staff = (k_pe_pe + k_pe_pa)
    d_hh = omega
    d_staffing_ratio = N_pe/N_pa
    d_adm = dyn_admission
    d_adm_R =  dyn_admission*(f_Cr + f_Cr*dyn_comm_denom)
    d_foi_S = foi_S_pa_pa + foi_S_pa_pe + foi_S_pe_pa + foi_S_pe_pe
    d_foi_Cr = foi_Cr_pa_pa + foi_Cr_pa_pe + foi_Cr_pe_pa + foi_Cr_pe_pe
    d_patientdays = N_pa
    d_staffdays = N_pe
    
    
    # return derivatives
    der <- c(d_S_U_pa, d_S_Cs_pa, d_S_Cr_pa, d_E_U_pa, d_E_Cs_pa, d_E_Cr_pa, d_I_U_pa, d_I_Cs_pa, d_I_Cr_pa, d_R_U_pa, d_R_Cs_pa, d_R_Cr_pa,
             d_S_U_pe, d_S_Cs_pe, d_S_Cr_pe, d_E_U_pe, d_E_Cs_pe, d_E_Cr_pe, d_I_U_pe, d_I_Cs_pe, d_I_Cr_pe, d_R_U_pe, d_R_Cs_pe, d_R_Cr_pe,
             d_SL,
             d_incS_pa_pa, d_incS_pe_pa, d_incS_pa_pe, d_incS_pe_pe,
             d_incCs_pa_pa, d_incCs_pe_pa, d_incCs_pa_pe, d_incCs_pe_pe,
             d_incCr_pa_pa, d_incCr_pe_pa, d_incCr_pa_pe, d_incCr_pe_pe, d_incCr_endog,
             d_abx, d_kappa_patients, d_kappa_staff, d_hh, d_staffing_ratio, d_adm, d_adm_R, d_foi_S, d_foi_Cr, d_patientdays, d_staffdays)
    list(der)
  }) 
}
