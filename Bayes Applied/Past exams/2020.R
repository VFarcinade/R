rm(list=ls())
data=read.csv('wikipedia.csv',header=F)

#---------------------------------------------------------
#model selec

d1=data[data$V2=="TFA" | data$V2=="OUT",]

hist(d1[,1],freq=T)
par(mfrow = c(1,2))
hist(data[data$V2=="TFA" ,1],freq=TRUE)
hist(data[data$V2=="OUT" ,1],freq=TRUE)

N_TFA=sum(data$V2=="TFA")
N_OUT=sum(data$V2=="OUT")

S_TFA=sum(data$V1[data$V2=="TFA"])
S_OUT=sum(data$V1[data$V2=="OUT"])


N=N_TFA+N_OUT
S=S_TFA+S_OUT

a=1
b=100

a1=1
b1=100

a2=1
b2=50


like1 = function(lambda){
  return(lambda^N * exp(-lambda*S))
}

like2 = function(lambda1, lambda2){
  return(lambda1^N_TFA * exp(-lambda1*S_TFA) * lambda2^N_OUT * exp(-lambda2*S_OUT))
}

N_mcmc=1e5

lambda=rgamma(N_mcmc,a,b)
lambda1=rgamma(N_mcmc,a1,b1)
lambda2=rgamma(N_mcmc,a2,b2)

mcmc1 = 0 
mcmc2 = 0 


for (i in seq(1,N_mcmc)){
  mcmc1 = mcmc1 + like1(lambda[i])
  mcmc2 = mcmc2 + like2(lambda1[i],lambda2[i])
}

print(log(mcmc1/N_mcmc))

analyt1 = a*log(b) + lgamma(N+a) - lgamma(a) - (N+a)*log(S+b)
print(log(exp(analyt1)))

print(log(mcmc2/N_mcmc))

print(-log(mcmc1/N_mcmc)+log(mcmc2/N_mcmc))

print("Strong evidence pour model 2")

#---------------------------------------------------------
#gibbs

X=data[,1]
n=length(X)


gibbs = function(niter, a0, b0, a1, b1) {
  lbd0 = rep(1, niter)
  lbd1 = rep(1, niter)
  z = array(0, dim=c(niter,n))
  pb = txtProgressBar(min = 0, max = niter, style = 3)
  
  
  for (i in 2:niter) {
    S0=sum(which(z[i-1,]==0))
    N0=sum(z[i-1,]==0)
    S1=sum(which(z[i-1,]==1))
    N1=sum(z[i-1,]==1)
    
    lbd0[i] = rgamma(1, a0+N0, b0+S0)
    lbd1[i] = rgamma(1, a1+N1, b1+S1)
    
    p=(lbd0[i]+lbd1[i])*X/(1-(lbd0[i]+lbd1[i])*X)
    ind = which(p>0 & p<1)
    z[i,ind] = rbinom(length(ind), 1, p[ind])
    setTxtProgressBar(pb, i)
  }
  return(list(lbd0=lbd0, lbd1=lbd1, z=z))
}


niter=1e5
nburn=1e3

library(MCMCpack)
run = gibbs(niter, 1, 100, 2, 200)
effectiveSize(run$lbd0)
effectiveSize(run$lbd1)
plot(run$lbd0[nburn:niter], type="l")
plot(run$lbd1[nburn:niter], type="l")

#---------------------------------------------------------
#IS

like1 = function(X, lbd0, lbd1, z){
  p=(lbd0*exp(-lbd0*X))^(1-z)*(lbd1*exp(-lbd1*X))^z
  return (prod(p))
}

like2 = function(X, lbd){
  p=lbd*exp(-lbd*X)
  return (prod(p))
}

N_mc=1e4
marg_like1 = 0
marg_like2 = 0

lbd0 = rgamma(N_mc, 1, 100)
lbd1 = rgamma(N_mc, 2, 200)
lbd = rgamma(N_mc, 1, 150)

z = rbinom(N_mc*n, 1, 0.5)
z = matrix(z, nrow = N_mc, byrow = TRUE)

for (i in seq(1:N_mc)){
  marg_like1 = marg_like1 + like1(X, lbd0[i], lbd1[i], z[i,])
  marg_like2 = marg_like2 + like2(X, lbd[i])
}

marg_like1 = marg_like1/N_mc
marg_like2 = marg_like2/N_mc


print(log(marg_like1) -log(marg_like2))



#---------------------------------------------------------



