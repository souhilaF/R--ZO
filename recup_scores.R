

################################################################################################
################################################################################################
fun_test<-function(a_i_b_i){
  
  return(exp(seq(log(a_i_b_i["a"]),log(a_i_b_i["b"]),length.out = a_i_b_i["k"])))
}

#######################################################################################################
#######################################################################################################



recup_scores <- function(Y){
  
  n <- length(Y[1,])
  N <- n*(n-1)/2
  
  Huge.edge = huge(Y, method='glasso')
  lambda.step = exp(mean(diff(log(Huge.edge$lambda))))
  lambda.min = min(Huge.edge$lambda)
  lambda.max = max(Huge.edge$lambda)
  lambda.nb = length(Huge.edge$lambda)
  
  #On veut qu'aucune aretes soient présentent au début donc que df=0 donc lambda max assez grand 
  while(min(Huge.edge$df)>0){
    lambda.max = lambda.max/(lambda.step)^5
    Huge.edge = huge(Y, method='glasso', lambda=lambda.max)
    
  }
  
  
  while(max(Huge.edge$df)<N){
    lambda.min = lambda.min*(lambda.step)^5
    Huge.edge = huge(Y, method='glasso', lambda=lambda.min)
    #plot(Huge.edge$lambda, Huge.edge$df, log='x'); abline(h=c(0, N))
    
  }
  
  
  
  
  lambda.seq = unique(sort(unlist(new.lambda.seq), decreasing = T))
  Huge.edge = huge(Y, method = 'glasso', lambda=lambda.seq)
  plot(Huge.edge$lambda, Huge.edge$df, log = 'x')
  
  longueur <- length(lambda.seq)
  
  ##construction matrice de scores
  S <- matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
     if(i != j){
       S[i,j] <- lambda.seq[min(which(lapply(Huge.edge$path, function(P) P[i,j])==1))]
     } 
    }
  }
  
  
  diag(S) = 0
  score <- log(S)
  diag(score) <- 0
  return(score)
}
#######################################################################################################

