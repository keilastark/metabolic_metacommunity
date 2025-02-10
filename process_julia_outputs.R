library(dplyr)
library(data.table)
library(cowplot)
library(scales)
reps <- 30
disp.results.df <- data.frame()
alpha.results.df <- data.frame()
alpha.disp.results.df <- data.frame()
asymmetry.results.df2 <- data.frame()
asymmetry.quick.results.df<-data.frame()
#dispersal outputs

for (warming_level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
for(i in 1:reps) {
  print(i)
  load(paste("./disp_outputfile_",warming_level, "_C_", i, ".RData", sep = ""))
  
  model_summary<- disp_model.df %>% 
    group_by(dispersal, topt, d_temp_func, warming, Time, Patch) %>% 
    dplyr::summarise(alpha_div = sum(N>0),
                     alpha_N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, topt,  d_temp_func, warming,Patch) %>% 
    dplyr::summarise(alpha_div = mean(alpha_div),
                     stability = mean(alpha_N, na.rm = TRUE)/sd(alpha_N,na.rm = TRUE),
                     alpha_N = mean(alpha_N))
                      
  
  model_summary_gamma<- disp_model.df %>% 
    group_by(dispersal, topt, d_temp_func,warming, Time, Species) %>% 
    dplyr::summarise(N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, topt, d_temp_func,warming, Time) %>% 
    dplyr::summarise(gamma_div = sum(N>0),
                     N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, topt,  d_temp_func,warming) %>% 
    dplyr::summarise(stability_gamma = mean(N)/sd(N),
                     gamma_div = mean(gamma_div))
  
  model_summary_t_gamma<- disp_model.df %>% 
    group_by(dispersal, topt, d_temp_func, warming,Patch, Species) %>% 
    dplyr::summarise(N = sum(N)) %>% 
    ungroup() %>% 
    group_by(dispersal, topt,  d_temp_func,warming, Patch) %>% 
    dplyr::summarise(gamma_div_time = sum(N>0)) %>% 
    ungroup() %>% 
    group_by(dispersal, topt,d_temp_func,warming) %>% 
    dplyr::summarise(gamma_div_time = mean(gamma_div_time))
  
  
  all_model_summary <- left_join(left_join(model_summary,model_summary_gamma),
                                           model_summary_t_gamma)
  
  all_model_summary <- all_model_summary %>% 
    mutate(beta_time = gamma_div_time/alpha_div, 
           beta_space = gamma_div/alpha_div)
  
  all_model_summary$rep <- i
  
  disp.results.df <- bind_rows(disp.results.df, all_model_summary)
  }
}
save(disp.results.df, file = "./final_dispersal_outputs.RData")

#alpha results

alpha.results.df <- data.frame()


for (level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for(i in 1:reps) {
    print(i)
    load(paste("./alpha_outputfile_", level, "_C_", i, ".RData", sep = ""))
    
    model_summary<- alpha_model.df %>% 
      group_by(alpha_temp_func, topt, alpha, warming,Time, Patch) %>% 
      dplyr::summarise(alpha_div = sum(N>0),
                       alpha_N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, alpha, warming,Patch) %>% 
      dplyr::summarise(alpha_div = mean(alpha_div),
                       stability = mean(alpha_N, na.rm = TRUE)/sd(alpha_N,na.rm = TRUE),
                       alpha_N = mean(alpha_N)) %>% 
      ungroup() %>% 
      mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
      group_by(alpha_temp_func, topt, warming,alpha) %>% 
      dplyr::summarise(alpha_div = mean(alpha_div),
                       stability = mean(stability, na.rm = TRUE),
                       alpha_N = mean(alpha_N))
    
    model_summary_gamma<- alpha_model.df %>% 
      group_by(alpha_temp_func, topt, alpha, warming,Time, Species) %>% 
      dplyr::summarise(N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, alpha,warming, Time) %>% 
      dplyr::summarise(gamma_div = sum(N>0),
                       N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, warming,alpha) %>% 
      dplyr::summarise(stability_gamma = mean(N)/sd(N),
                       gamma_div = mean(gamma_div))
    
    model_summary_t_gamma<- alpha_model.df %>% 
      group_by(alpha_temp_func, topt, alpha,warming, Patch, Species) %>% 
      dplyr::summarise(N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, alpha,warming, Patch) %>% 
      dplyr::summarise(gamma_div_time = sum(N>0)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, alpha,warming) %>% 
      dplyr::summarise(gamma_div_time = mean(gamma_div_time))
    
    all_model_summary <- left_join(left_join(model_summary,model_summary_gamma),
                                   model_summary_t_gamma)
    
    all_model_summary <- all_model_summary %>% 
      mutate(beta_time = gamma_div_time/alpha_div, 
             beta_space = gamma_div/alpha_div)
    
    all_model_summary$rep <- i
    
    alpha.results.df <- bind_rows(alpha.results.df, all_model_summary)
    
  }
}

save(alpha.results.df, file = "./final_alpha_outputs.RData")

###

for (level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for(i in 1:reps) {
    print(i)
    load(paste("./alpha_disp_outputfile_", level, "_C_", i, ".RData", sep = ""))
    
    model_summary<- alpha_disp_model.df %>% 
      group_by(alpha_temp_func, topt, d_temp_func, warming,Time, Patch) %>% 
      dplyr::summarise(alpha_div = sum(N>0),
                       alpha_N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func, warming,Patch) %>% 
      dplyr::summarise(alpha_div = mean(alpha_div),
                       stability = mean(alpha_N, na.rm = TRUE)/sd(alpha_N,na.rm = TRUE),
                       alpha_N = mean(alpha_N)) %>% 
      ungroup() %>% 
      mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
      group_by(alpha_temp_func, topt, warming,d_temp_func) %>% 
      dplyr::summarise(alpha_div = mean(alpha_div),
                       stability = mean(stability, na.rm = TRUE),
                       alpha_N = mean(alpha_N))
    
    model_summary_gamma<- alpha_disp_model.df %>% 
      group_by(alpha_temp_func, topt, d_temp_func, warming,Time, Species) %>% 
      dplyr::summarise(N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming, Time) %>% 
      dplyr::summarise(gamma_div = sum(N>0),
                       N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming) %>% 
      dplyr::summarise(stability_gamma = mean(N)/sd(N),
                       gamma_div = mean(gamma_div))
    
    model_summary_t_gamma<- alpha_disp_model.df %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming, Patch, Species) %>% 
      dplyr::summarise(N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming, Patch) %>% 
      dplyr::summarise(gamma_div_time = sum(N>0)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming) %>% 
      dplyr::summarise(gamma_div_time = mean(gamma_div_time))
    
    all_model_summary <- left_join(left_join(model_summary,model_summary_gamma),
                                   model_summary_t_gamma)
    
    all_model_summary <- all_model_summary %>% 
      mutate(beta_time = gamma_div_time/alpha_div, 
             beta_space = gamma_div/alpha_div)
    
    all_model_summary$rep <- i
    
    alpha.disp.results.df <- bind_rows(alpha.disp.results.df, all_model_summary)
    
  }
}

save(alpha.disp.results.df, file = "./final_alpha_disp_outputs.RData")




for (level in c("0.0", "2.5", "5.0", "7.5", "10.0")) {
  for(i in 1:15) {
    print(i)
    load(paste("./asymmetry_outputfile_", level, "_C_", i, ".RData", sep = ""))
    
    model_summary<- asymmetry_model.df2 %>% 
      group_by(alpha_temp_func, topt, d_temp_func, warming,Time, Patch) %>% 
      dplyr::summarise(alpha_div = sum(N>0),
                       alpha_N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func, warming,Patch) %>% 
      dplyr::summarise(alpha_div = mean(alpha_div),
                       stability = mean(alpha_N, na.rm = TRUE)/sd(alpha_N,na.rm = TRUE),
                       alpha_N = mean(alpha_N)) %>% 
      ungroup() %>% 
      mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
      group_by(alpha_temp_func, topt, warming,d_temp_func) %>% 
      dplyr::summarise(alpha_div = mean(alpha_div),
                       stability = mean(stability, na.rm = TRUE),
                       alpha_N = mean(alpha_N))
    
    model_summary_gamma<- asymmetry_model.df2 %>% 
      group_by(alpha_temp_func, topt, d_temp_func, warming,Time, Species) %>% 
      dplyr::summarise(N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming, Time) %>% 
      dplyr::summarise(gamma_div = sum(N>0),
                       N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming) %>% 
      dplyr::summarise(stability_gamma = mean(N)/sd(N),
                       gamma_div = mean(gamma_div))
    
    model_summary_t_gamma<- asymmetry_model.df2 %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming, Patch, Species) %>% 
      dplyr::summarise(N = sum(N)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming, Patch) %>% 
      dplyr::summarise(gamma_div_time = sum(N>0)) %>% 
      ungroup() %>% 
      group_by(alpha_temp_func, topt, d_temp_func,warming) %>% 
      dplyr::summarise(gamma_div_time = mean(gamma_div_time))
    
    all_model_summary <- left_join(left_join(model_summary,model_summary_gamma),
                                   model_summary_t_gamma)
    
    all_model_summary <- all_model_summary %>% 
      mutate(beta_time = gamma_div_time/alpha_div, 
             beta_space = gamma_div/alpha_div)
    
    all_model_summary$rep <- i
    
    asymmetry.results.df <- bind_rows(asymmetry.results.df, all_model_summary)
    
  }
}

save(asymmetry.results.df, file = "./final_asymmetry_outputs.RData")
