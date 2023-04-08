l2_dist <- function(x1,x2){
  return(sqrt(sum((x1-x2)^2)))
}

summary_stat <- function(x){
  AA=sum(x=='AA')
  aa=sum(x=='aa')
  return(c(aa,AA))
}


sample_prior_beta<-function(N,a1,b1,a2,b2){
  return(list(rbeta(N,a1,b1),rbeta(N,a2,b2)))
}

ABC <- function(X,N,sample_prior,epsilon,dist,...){
  
  n=length(X)
  S=summary_stat(X)
  
  param=sample_prior(N,...)
  p=unlist(param[1])
  theta=unlist(param[2])
  keep=logical(N)
  
  p_aa=p*theta +(1-p)*theta^2
  p_AA=p*(1-theta)+(1-p)*(1-theta)^2
  p_aA=1-(p_aa+p_AA)
  
  for (i in 1:N){
    X_=sample(c('aa','AA','aA'), n, replace = TRUE, prob = c(p_aa[i],p_AA[i],p_aA[i]))
    S_=summary_stat(X_)
    if(dist(S,S_)<epsilon) keep[i]=TRUE
    if(dist(S,S_)<epsilon){
      AA=sum(X_=='AA')
      aa=sum(X_=='aa')
      aA=sum(X_=='aA')
    }
    
  }
  
  return(list(p[keep],theta[keep],sum(keep)))
}


sample_prior_model_unif<-function(N){
  return(sample(c(0,1),N,replace=TRUE))
}


ABC_model_selec <- function(X,N,sample_prior,epsilon,dist){
  
  n=length(X)
  S=summary_stat(X)
  
  model=sample_prior_model_unif(N)
  param=sample_prior(N,1,1,1,1)
  
  p=unlist(param[1])
  theta=unlist(param[2])
  keep=logical(N)
  
  ind_m0=which(model==0)  
  ind_m1=which(model==1)
  
  p_aa=double(N)
  p_AA=double(N)
  p_aA=double(N)
  
  p_aa[ind_m0]=(1-theta[ind_m0])^2
  p_AA[ind_m0]=theta[ind_m0]^2
  p_aA[ind_m0]=1-(p_aa[ind_m0]+p_AA[ind_m0])
  
  p_aa[ind_m1]=p[ind_m1]*theta[ind_m1] +(1-p[ind_m1])*theta[ind_m1]^2
  p_AA[ind_m1]=p[ind_m1]*(1-theta[ind_m1])+(1-p[ind_m1])*(1-theta[ind_m1])^2
  p_aA[ind_m1]=1-(p_aa[ind_m1]+p_AA[ind_m1])
  
  for (i in 1:N){
    X_=sample(c('aa','AA','aA'), n, replace = TRUE, prob = c(p_aa[i],p_AA[i],p_aA[i]))
    S_=summary_stat(X_)
    if(dist(S,S_)<epsilon) keep[i]=TRUE
  }
  
  return(list(model[keep],sum(keep)))
}
