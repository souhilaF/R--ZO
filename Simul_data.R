library(MASS)

#######################################################################################################
#-----------------------------------------Simulation donnees

SIMULDATA <- function(GammaArrete, PiMelange, n,P) {
  N = n*(n-1)/2
  #fonction de classes des observations
  Z_m <- rmultinom(n, size = 1, prob = PiMelange)
  Z_v <- which(Z_m == 1, 1)
  Z_v <- Z_v[,1] 
  #matrice des proba d'arrete gamma entre observations de chaque classe
  M <- matrix(0, ncol = length(PiMelange), nrow = length(PiMelange))
  colnames(M) = paste("Z", 1:length(PiMelange), sep = "")
  rownames(M) = paste("Z", 1:length(PiMelange), sep = "")
  
  M[lower.tri(M, diag=T)] <- GammaArrete
  B <- t(M)
  diag(B) <- 0
  M <- B + M
  
  #G la matrice d'adjacence
  G <- matrix(rbinom(n*n, size = 1, prob = M[Z_v,Z_v]), n, n)
  
  
  #On construit le graph d'adjacence non dirigé:
  #Arbitrairement on récupère la matrice superieure de la matrice d'adjacence R et on la symetrise.
  #ATTENTION ne pas faire tourner l'algo suivant si on veut un graph dirigé
  
  
  G[lower.tri(G, diag=T)] <- 0
  G <- t(G)+G
  
  
  
  #On construit la matrice de précision Omega = G + nu_var*D
  #et on trouve le nu_var le plus petit possible pour rendre la matrice Omega inversible
  #D est la matrice de degré
  
  
  
  D2 <- rowSums(G) 
  #on fait la somme sur les lignes [trouver un autre moyen avec les matrice d'adjacence dirigée(sum sur ligne different que sur colonne)]
  #on met ce vecteur sur une matrice diagonale pour récuperer la matrice de degré D
  D <- matrix(0, ncol=n, nrow=n)
  diag(D) <- D2 + (D2 == 0) * (min(D2[which(D2 != 0)])*0.01)
  
  
  
  
  
  
  #construction de la matrice Omega
  nu_var <- 0
  mu_var <- 0
  Omega <- (nu_var + 0.1) * D + (mu_var + 0.3) * G 
  
  while((min(eigen(Omega)$values)) <= 1e-4) {
    nu_var <- nu_var + 0.1
    mu_var <- mu_var + 0.1
    Omega <- (nu_var + 0.1) * D + (mu_var + 0.3) * G 
    
  }
  Omega <- Omega / mean(Omega)
  Sigma <- solve(Omega)
  
  
  #On génère des simulation de gaussiennes centrées de variance sigma sur 4 sites:
  
  Y <- mvrnorm(P, mu = rep(0,n), Sigma = Sigma)
  OutSimul <- list(Y = Y, G = G, Z = Z_v)
  return(OutSimul)
}




################################################################################################
#--------------------------bootstrap
bootst<-function(data){
  b<-sample(1:nrow(data),nrow(data),replace=T)
  result<-data[b,]
  return(result)
  
}

