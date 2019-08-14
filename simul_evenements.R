rm(list=ls())
setwd("~/Bureau/VEM_for_estimated_networks/Codes/V2")


############################################################################################

#############################################################################################
source('Simul_data.R')
source('Main_function_last.R')
source('VEM_last.R')
#############################################################################################
# GammaArrete <- c(0.70, 0.30, 0.05, 0.1, 0.1, 0.5)
# PiMelange <- c(1/6, 1/3, 1/2)
# paramsim <- list(GammaArrete = GammaArrete, PiMelange  = PiMelange)
# n  =50; 
# val_p <- c(20,50,100,500)
# for (p in val_p){
#   for (sim in 1:100){ 
#     
#     datasim <- SIMULDATA(GammaArrete, PiMelange, n,p)
#     file_name = paste('res_simu/datasim/datasim_n',n,'_p',p,'_sim',sim,'.Rdata',sep='')
#     save(datasim,paramsim,file = file_name)
#     
#     }
# }

#------------------------------- REcupScore + VEM

n  = 50; 
val_p <- c(500,100,50,20)
K  = 3
N=1225
for (p in val_p) {
  for (sim in 1:100) { 
     
     print(paste("p=", p,", sim = ",sim))
     file_name_result = paste('res_simu/ressim/ressim_n',n,'_p',p,'_sim',sim,'.Rdata',sep='')
    
     if(!file.exists(file_name_result)){
      file_name_data = paste('res_simu/datasim/datasim_n',n,'_p',p,'_sim',sim,'.Rdata',sep='')
      load(file_name_data)
      score <- recup_scores(datasim$Y)
      output <- VEM(S = score,K, niter=1000, epsilon_tau=1e-4, epsilon_eta = 1e-4,verbose = FALSE)
      save(score,output,file = file_name_result)
       datasim = c(); paramsim = c(); output = c(); score = c()
      }
    }
  }







