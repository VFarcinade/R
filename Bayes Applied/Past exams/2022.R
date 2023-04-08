rm(list=ls())
data = read.csv("examdata2022.csv")
name=unique(data$white)


loglike = function(z,gamma,delta){
  res=0
  
  for (i in name){
    for (j in name){
      if(i!=j) {
        
        zi=(z[name==i])
        zj=(z[name==j])
        
        if(data$result[data$white==i & data$black==j] == 1) res = res + zi
        if(data$result[data$white==i & data$black==j] == 0) res = res + zj + gamma
        if(data$result[data$white==i & data$black==j] == 0) res = res + delta + 0.5*(gamma+zi+zj)
        
        res = res - log(exp(zi)+exp(zj + gamma)+ exp(delta + 0.5*(gamma+zi+zj)))
      }
    }
  }
  return (res)
}

logprior = function(z){
  sum((dnorm(z,0,1,log=T)))
}

logpost = function(param, d = data){
  name=unique(data$white)
  z = param[1:8]
  gamma = param[9]
  delta = param[10]
  res = logprior(z) + loglike(z,gamma,delta)
  return (res)
}



library(MCMCpack)
out = MCMCmetrop1R(logpost, c(rep(0., 8), 0., 0.), mcmc=1e3, burnin=5e1,tune=0.75)

summary(out)
plot(out[,9:10])
effectiveSize(out)
acf(out[,7:10])

par(mfrow=c(2,2))
hist(out[, 1])
hist(out[, 2])
hist(out[, 9])
hist(out[, 10])

z = out[,1:8]
gamma = out[, 9]
delta = out[, 10]

#Q20
print(mean(gamma<0))

#Q21
for(i in (1:8)){
  for (j in (1:8)){
    if(i != j) {
      prob = mean(z[,i]>z[,j])
      if(prob > 0.95) cat(name[i], "VS", name[j], ": ", mean(z[,i]>z[,j]), "\n")
    }
  }
}

#Q22

argMax=apply(z,1,which.max)
print(mean(argMax==which(name=="Karjakin")))
