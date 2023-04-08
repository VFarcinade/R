y = c(28, 8, -3, 7, -1, 1, 18, 12)
s = c(15, 10, 16, 11, 9, 11, 10, 18)
n = length(y)

# Q2
# constant prior on the mean
# plot the posterior for each theta_i
curve(dnorm(x, y[1], s[1]), 
      xlim=c(-50, 80), ylim = c(0, .1))
for(i in 2:8){
  curve(dnorm(x, y[i], s[i]), add=T, col=i)
}

# Q3
# With a constant prior on the mean, the posterior is normal N(a, b^2) with
b = sqrt(1 / sum(1 / s ^ 2))
a = b^2 * sum(y / s ^ 2)
curve(dnorm(x, a, b), add=T, lty=2)
# dashed line: posterior for mu in the model of question 3
# coloured, full lines: posterior for each theta_i in model of question 2

curve(dnorm(x, a, b), xlim=c(-5, 30), lty=2)

niter = 1e4
r = rep(NA, niter)
mu = rep(NA, niter)
for(i in 1:niter){
  mu[i] = rnorm(1, a, b)
  ysim = rnorm(8, mu[i], s)
  r[i] = ysim[1] > ysim[3]
}

hist(r)
mean(r)

# Q5
require(mvtnorm)
loglkd = function(theta, mu, tau, y, s) {
  if (tau < 0) {
    return(-Inf)
  }
  return(sum(dnorm(theta, mu, tau, log = T) + 
               sum(dnorm(y, theta, s, log = T))))
}

logprior = function(theta, mu, tau){
  return(0)
}

logpost = function(theta, mu, tau, y, s){
  return(loglkd(theta, mu, tau, y, s) + 
           logprior(theta, mu, tau))
}

logpost2 = function(par){
  # theta is par[1:8]
  # mu is par[9]
  # tau is par[10]
  return(loglkd(par[1:8], par[9], par[10], y, s) + 
           logprior(par[1:8], par[9], par[10]))
}

library(MCMCpack)
out = MCMCmetrop1R(logpost2, c(rep(0, 8),0, 10), mcmc=2e5, burnin=1e3,tune=1)
summary(out)
plot(out)
effectiveSize(out)
raftery.diag(out)
acf(out)
hist(out[, 1])

# Alternative package:
library(mcmc)
out = metrop(logpost2, c(rep(0, 8),0,10), 5e5, scale=.3)
out$accept
plot(ts(out$batch[, 10]))
hist(out$batch[, 9])
effectiveSize(out$batch)


# or a MH by hand:
MH = function(niter, y, s, omega = 1) {
  n = length(y)
  theta = matrix(NA, niter, n)
  mu = rep(NA, niter)
  tau = rep(NA, niter)
  
  #Initialization
  theta[1, ] = 0
  mu[1] = 0
  tau[1] = 10
  acc = 0
  
  for (t in 2:niter) {
    thetaprop = rmvnorm(1, theta[t - 1, ], diag(rep(omega ^ 2, n)))
    muprop = rnorm(1, mu[t-1], omega)
    tauprop = rnorm(1, tau[t-1], omega)
    alpha = logpost(thetaprop, muprop, tauprop, y, s) - 
      logpost(theta[t - 1, ], mu[t-1], tau[t-1], y, s)
    if (runif(1) < exp(alpha)) {
      theta[t,] = thetaprop
      mu[t] = muprop
      tau[t] = tauprop
      acc = acc + 1
    }
    else{
      theta[t,] = theta[t - 1, ]
      mu[t] = mu[t-1]
      tau[t] = tau[t-1]
    }
  }
  print(acc / niter) #acceptance ratio
  return(list(theta, mu, tau))
}

niter = 5e4
run1 = MH(niter, y, s, 3)
#ESS
niter / (2 * sum(acf(run1[[2]], lag = 300, plot=F)$acf))

par(mfrow=c(3, 1))
plot(run1[[3]], type="l")
acf(run1[[3]], lag=300)
hist(run1[[3]])

# Q6


theta.out = out$batch[, 1:8]

thetamax = apply(theta.out, 1, max)
hist(thetamax)

mean(theta.out[, 1] > theta.out[, 3])

#Q8: posterior predictive checks
nrep = 200
thetasim = theta.out[sample.int(niter, nrep), ]

yrep = matrix(NA, nrep, n)
for (i in 1:nrep) {
  yrep[i, ] = rnorm(n, thetasim[i, ], s)
}



hist(yrep[, 1], breaks=30)
mean(yrep[, 1] > 28) # posterior predictive p-value for the observation y1=28.

yrep.min = apply(yrep, 1, min)
yrep.max = apply(yrep, 1, max)
hist(yrep.min, breaks=30)
mean(yrep.max > 28)
mean(yrep.min < -3)
mean(rowMeans(yrep) < mean(y))


# posterior predictive probability of the event that Y_1 > Y_3
pred1 = rnorm(nrow(out$batch), out$batch[, 1], s[1])
pred3 = rnorm(nrow(out$batch), out$batch[, 3], s[3])
mean(pred1>pred3)

# Q9
# Use a Cauchy distribution
loglkd2 = function(theta, mu, tau, y, s) {
  if (tau < 0) {
    return(-Inf)
  }
  return(sum(dcauchy(theta, mu, tau, log = T) + sum(dnorm(y, theta, s, log = T))))
}

logpost2 = function(theta, mu, tau, y, s){
  return(loglkd2(theta, mu, tau, y, s) + 
           logprior(theta, mu, tau))
}

MH2 = function(niter, y, s, omega = 1) {
  n = length(y)
  theta = matrix(NA, niter, n)
  mu = rep(NA, niter)
  tau = rep(NA, niter)
  
  #Initialization
  theta[1, ] = 0
  mu[1] = 0
  tau[1] = 10
  acc = 0
  
  for (t in 2:niter) {
    thetaprop = rmvnorm(1, theta[t - 1, ], diag(rep(omega ^ 2, n)))
    muprop = rnorm(1, mu[t-1], omega)
    tauprop = rnorm(1, tau[t-1], omega)
    alpha = logpost2(thetaprop, muprop, tauprop, y, s) - 
      logpost2(theta[t -  1, ], mu[t-1], tau[t-1], y, s)
    if (runif(1) < exp(alpha)) {
      theta[t,] = thetaprop
      mu[t] = muprop
      tau[t] = tauprop
      acc = acc + 1
    }
    else{
      theta[t,] = theta[t - 1, ]
      mu[t] = mu[t-1]
      tau[t] = tau[t-1]
    }
  }
  print(acc / niter) # acceptance ratio
  return(list(theta, mu, tau))
}

run2 = MH2(5e4, y, s, 3)
par(mfrow = c(2, 1))
hist(run2[[2]], breaks=30)
hist(run1[[2]], breaks=30)
summary(run2[[2]])
summary(run1[[2]])

# Q10
AMH = function(niter, y, s) {
  n = length(y)
  param = matrix(NA, niter, n + 2) #The first n columns correspond to theta, then mu, then tau
  acc = 0
  
  #Initialization
  param[1, ] = c(rep(0, n + 1), 10)
  param[2, ] = rmvnorm(1, param[1, ], diag(rep(1, n + 2)))
  
  for (t in 3:niter) {
    if (t < 100) {
      paramprop = rmvnorm(1, param[t - 1, ], diag(1, n + 2, n + 2))
    }
    else{
      paramprop = rmvnorm(1, param[t - 1, ], 2.38 ^ 2 / (n + 2) * 
                            cov(param[max(1, t - 1e4):(t - 1), ]))
    }
    alpha = logpost(paramprop[1:n], paramprop[n + 1], paramprop[n + 2], y, s) -
      logpost(param[t - 1, 1:n], param[t - 1, n + 1], param[t - 1, n + 2], y, s)
    if (runif(1) < exp(alpha)) {
      param[t, ] = paramprop
      acc = acc + 1
    }
    else{
      param[t, ] = param[t - 1, ]
    }
  }
  print(acc / niter) #acceptance rate
  return(param)
}
niter = 5e4
run3 = AMH(niter, y, s)
niter / (2 * sum(acf(run3[1e4:niter, 9], lag = 100)$acf))

system.time(AMH(5e4, y, s))

# A faster version, which updates the covariance matrix less and less often
AMH2 = function(niter, y, s) {
  n = length(y)
  param = matrix(NA, niter, n + 2) #The first n columns correspond to theta, then mu, then tau
  lambda = rep(1, n + 2)
  mu = rep(0, n + 2)
  Sigma = diag(rep(1, n + 2))
  Omega = diag(sqrt(lambda)) %*% Sigma %*% diag(sqrt(lambda))
  acc = 0
  
  #Initialization
  param[1, ] = c(rep(0, n + 1), 10)
  Sigma = diag(rep(1, n + 2))
  
  for (t in 2:niter) {
    paramprop = rmvnorm(1, param[t - 1, ], Sigma)
    
    alpha = logpost(paramprop[1:n], paramprop[n + 1], paramprop[n + 2], y, s) -
      logpost(param[t - 1, 1:n], param[t - 1, n + 1], param[t - 1, n + 2], y, s)
    if (runif(1) < exp(alpha)) {
      param[t, ] = paramprop
      acc = acc + 1
    }
    else{
      param[t, ] = param[t - 1, ]
    }
    if (t > 100 &
        runif(1) < 1 / t) {
      Sigma = 2.38 ^ 2 / (n + 2) * cov(param[1:t, ])
    }
  }
  print(acc / niter) #acceptance rate
  return(param)
}
run4 = AMH2(niter, y, s)
niter / (2 * sum(acf(run4[1e4:niter, 9], lag = 300)$acf))

system.time(AMH2(5e4, y, s))

library(rstan)
rstan_options(auto_write = T)

schools_dat = list(J = 8, 
                   y = c(28,  8, -3,  7, -1,  1, 18, 12),
                   sigma = c(15, 10, 16, 11,  9, 11, 10, 18))


fit = stan(file = 'schools.stan', data = schools_dat)

print(fit)
plot(fit)
traceplot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))
