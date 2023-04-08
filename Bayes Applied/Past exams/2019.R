rm(list=ls())
data=read.csv('examdata2019.csv', header=T)

K=33500000

Y=data$deaths
T=length(Y)

p_obj=2242/8000000000

b = 1
a = b*p_obj/(1-p_obj)
b = 1

x=seq(0,1,0.01)
plot(x,dbeta(x,a,b),type='l')

B=rbeta(1000000,a,b)
print(mean(B))
print(var(B))

#____________________
#MargLike
like = function(Y,p){
  return(prod(dbinom(Y,K,p)))
}

Lbinom = function(N,n){
  return(lgamma(N+1)-lgamma(N-n+1)-lgamma(n+1))
}

loglike = function(Y,p){
  res = Lbinom(K,Y) + Y*log(p) + (K-Y)*log(1-p) 
  return(sum(res))
}

N_mc=1e4
MargLike=0
p=rbeta(N_mc,a,b)

for (i in seq(1,N_mc)) MargLike = MargLike + exp(loglike(Y,p[i]))
  
MargLike = MargLike/N_mc

print(log(MargLike))


#------------------------------------

D2015_summer=data$deaths[data$year==2015 & (data$month == 7 | data$month == 8 | data$month == 9)]
D2015_autumn=data$deaths[data$year==2015 & (data$month == 10 | data$month == 11 | data$month == 12)]

b1 = 1
a1 = 1
b2 = 1
a2 = 0.5

print(a1/(a1+b1))
print(a2/(a2+b2))

N_mc=1e4
MargLike0 = 0
MargLike1 = 0

p=rbeta(N_mc,a,b)
p1=rbeta(N_mc,a1,b1)
p2=rbeta(N_mc,a2,b2)

for (i in seq(1,N_mc)){
  MargLike0 = MargLike + exp(loglike(Y,p[i]))
  MargLike1 = MargLike + exp(loglike(D2015_summer,p1[i])+loglike(D2015_autumn,p2[i]))
} 

MargLike0 = MargLike0/N_mc
MargLike1 = MargLike1/N_mc

print(log(MargLike0))
print(log(MargLike1))

print(log(MargLike0) - log(MargLike1))

b = 10
a = b*p_obj/(1-p_obj)

T=length(Y)
S=sum(Y)
LogMargLike0 = sum ( Lbinom(K,Y) )  - lgamma(K*T+a+b) + lgamma(S+a) + lgamma(K*T-S+b)



print(LogMargLike0)
exp(LogMargLike0)
