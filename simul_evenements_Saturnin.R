rm(list=ls())
setwd("~/Bureau/VEM_for_estimated_networks/Codes/V2")


############################################################################################

#############################################################################################
source('Simul_data_last.R')
source('Main_function_last.R')
source('VEM_last.R')


library(saturnin)


n  = 50; 
val_p <- c(500,100,50,20)
K  = 3
N=1225
for (p in val_p) {
  for (sim in 1:100) { 
    
    print(paste("p=", p,", sim = ",sim))
    file_name_data = paste('res_simu/datasim/datasim_n',n,'_p',p,'_sim',sim,'.Rdata',sep='')
    load(file_name_data)
    
    
    weights = lweights_gaussian(datasim$Y)
    score = edge.prob(weights, log = TRUE)
    a <- 1
    while( min( min(score, (1-score)) ) < 1e-10 ){
  
      a <- a + 0.1
      
      weights <- lweights_gaussian(datasim$Y) / a
      score <- edge.prob(weights, log = TRUE)
      score <- mat_vect_low(score)
      
    }
    
    
    score_logit = log(score/(1-score))
    file_name_result = paste('res_simu/ressimSaturnin/ressimSaturnin_n',n,'_p',p,'_sim',sim,'.Rdata',sep='')
    save(score,file = file_name_result)
    output <- VEM(S = score_logit,K, niter=1000, epsilon_tau=1e-4, epsilon_eta = 1e-4,verbose = FALSE)
    save(score_logit, a, output, file=file_name_result)
    datasim = c(); paramsim = c(); output = c(); score = c()
    
  }
        }
    
    
