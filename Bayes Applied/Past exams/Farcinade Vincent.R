rm(list=ls())

data = read.csv("RoadFatalities.csv", sep=',', header=T)
Y = data$y
t = as.double(data$t)

#Question 2.b)

a=1
b=1
qgamma(c(.025, .975), a + sum(Y), 1/(b + length(Y)) )

a=2
b=2
qgamma(c(.025, .975), a + sum(Y), 1/(b + length(Y)) )


#----------------------------------------
#Question 5
sigma_a = 1e4
sigma_b = 1

logprior0 = function(param){
  a=param[1]
  b=param[2]
  return(dnorm(a, 0, sigma_a, log = T) + dnorm(b, 0, sigma_b, log = T))
}

loglike0 = function(param){
  a=param[1]
  b=param[2]
  return(sum(-(a+b*t) + Y*log(a+b*t)))
}


logpost0 = function(param){
  return( logprior0(param) + loglike0(param) )
}

param0 = glm(data$y~data$t, family = poisson(link=identity))$coefficients
param0=c(as.double(param0[1]),as.double(param0[2]))


#----------------------------------------
#Question 6 & 7

library(MCMCpack)
out = MCMCmetrop1R(logpost0, param0, mcmc=1e4, burnin=5e1,tune=0.5)
#plot(out)

out = MCMCmetrop1R(logpost0, param0, mcmc=1e4, burnin=5e1,tune=2)
#plot(out)


#----------------------------------------
#Question 8

#A DECOMMENTER

out0 = MCMCmetrop1R(logpost0, param0, mcmc=1e5, burnin=5e1,tune=2)
#plot(out0)
#effectiveSize(out0)
#acf(out0)

#----------------------------------------
#Question 9

alpha=out0[,1]
beta=out0[,2]

quantile(alpha, probs =c(0.025,0.975) )
quantile(beta, probs =c(0.025,0.975) )


#----------------------------------------
#Question 10
Tc=1992
Y1 = Y[t<=Tc]
Y2 = Y[t>Tc]
t1 = t[t<=Tc]
t2 = t[t>Tc]

logprior = function(param){
  a1=param[1]
  b1=param[2]
  a2=param[3]
  b2=param[4]
  
  return(dnorm(a1, 0, sigma_a, log = T) + dnorm(b1, 0, sigma_b, log = T) 
         +dnorm(a2, 0, sigma_a, log = T) + dnorm(b2, 0, sigma_b, log = T))
}

loglike = function(param){
  a1=param[1]
  b1=param[2]
  a2=param[3]
  b2=param[4]
  
  return(sum(-(a1+b1*t1) + Y1*log(a1+b1*t1)) + sum(-(a2+b2*t2) + Y2*log(a2+b2*t2)))
}

logpost = function(param){
  return( logprior(param) + loglike(param) )
}

param0=double(4)
param0[1:2] = glm(Y1~t1, family = poisson(link=identity))$coefficients
param0[3:4] = glm(Y2~t2, family = poisson(link=identity))$coefficients

#A DECOMMENTER

out = MCMCmetrop1R(logpost, param0, mcmc=1e5, burnin=5e1,tune=1.5)

#plot(out[,1:2])
#plot(out[,3:4])
#effectiveSize(out)
#acf(out)

#----------------------------------------
#Question 12

beta1 = out[,2]
beta2 = out[,4]
N_mc=length(beta1)

mc_estim=cumsum(beta2<beta1)/seq(1,N_mc)

plot(mc_estim[1e4:N_mc],type='l')

cat("MCMC(p)=",mc_estim[N_mc])

#----------------------------------------
#Question 15
alpha1 = out[,1]
alpha2 = out[,3]
beta1 = out[,2]
beta2 = out[,4]

print(mean(alpha1))
print(var(alpha1))

print(mean(alpha2))
print(var(alpha2))

print(mean(beta1))
print(var(beta2))

print(mean(beta2))
print(var(beta2))

#----------------------------------------
#Question 16

N_mc=1e4

margLike0 = 0
margLike1 = 0

Llike0 = function(param){
  a=param[1]
  b=param[2]
  return(sum(-(a+b*t)+Y*log(a+b*t) - lgamma(Y+1)))  
}


Llike1 = function(param){
  a1=param[1]
  b1=param[2]
  a2=param[3]
  b2=param[4]
  return( sum(-(a1+b1*t1) + Y1 *log(a1+b1*t1) - lgamma(Y1+1) ) 
          + sum(-(a2+b2*t2) + Y2 *log(a2+b2*t2) - lgamma(Y2+1) ) ) 
}

for (i in seq(1,N_mc)){
  param0=c(out0[i,1],out0[i,2])
  param1=c(out[i,1],out[i,2],out[i,3],out[i,4])  
  margLike0 = margLike0 + exp(Llike0(param0))
  margLike1 = margLike1 + exp(Llike1(param1))
}

margLike0=margLike0/N_mc
margLike1=margLike1/N_mc

cat("log(BF(2,1))=", log(margLike1)-log(margLike0))


