
############################################################
############ Fonctions de bases à utiliser dans R ####################
############################################################


#librairies nécessaires
library(sna)
library(MASS)
library(huge)
library(glasso)
library(blockmodels)
library(mclust)
library(ggplot2)
library(dplyr)


#############################################################
#Construit une matrice symétrique à partir d'un vecteur
#############################################################
vect_mat_low <- function(V, diag=F)
{
  N <- length(V)
  if(diag == F) {n <- (1+sqrt(1+8*N))/2} else{ n <- (-1+sqrt(1+8*N))/2}
  M <- matrix(0,n,n) 
  M[lower.tri(M, diag=diag)] <- V
  B <- t(M)
  diag(B) <- 0
  M <- B + M
  return(M)
}


###########################################################################
#Construit un vecteur à partir d'une matrice triangulaire inférieure
###########################################################################
mat_vect_low <- function(M, diag=F)
{
  V <- M[lower.tri(M, diag = diag)]
  return(V)
}


############################################################################
#Passage d'une matrice binaire à un vecteur des indices 
############################################################################
mat_bin_to_vect_ind <- function(n,pi_melange)
{
  Z_m <- rmultinom(n, size = 1, prob = pi_melange)
  Z_v <- which(Z_m == 1, 1)
  Z_v <- Z_v[,1] 
  return(Z_v, Z_m)
}


###########################################################################
#recuperer les indices de la matrice triangulaire inferieure 
###########################################################################

indices <- function(n, diag=F)
{
  N  <- (n * (n-1)/2)*(diag==F) +(n * (n+1)/2)*(diag==T)
  S <- vect_mat_low(c(1:N), diag=diag)
  S[upper.tri(S)]=0
  return(which(S!=0,arr.ind=T))
  
}

############################################################################
#passer de n à N et inversement:
############################################################################

#N <- n*(n-1)/2
#n <-(1+sqrt(1+8*N))/2

#############################################################################
#############################################################################
fun_test<-function(a_i_b_i){
  
  return(exp(seq(log(a_i_b_i["a"]),log(a_i_b_i["b"]),length.out = a_i_b_i["k"])))
}
#############################################################################
#############################################################################

borne_inf <- function(tau_hat,Pi_hat,A,eta0,eta1){
borne_inf <- sum(tau_hat %*% log(Pi_hat))  - sum(tau_hat * log(tau_hat)) +
  .5*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
    t(tau_hat[,k]) %*% as.matrix(vect_mat_low(A[,k,l])) %*% tau_hat[,l]})
  }))) - 
  .5*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
    t(tau_hat[,k]) %*% 
      as.matrix(vect_mat_low(eta0[, k, l]*log(eta0[, k, l]) +
                               eta1[, k, l]*log(eta1[, k, l]))) %*% 
      tau_hat[,l]})
  }))) 
return(borne_inf)
}
