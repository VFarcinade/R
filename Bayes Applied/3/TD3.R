rm(list=ls())
d = read.csv("deathrate2.csv")[, -1]
y = d[ , 1]
X = as.matrix(d[ , -1])

n = length(y)
p = ncol(X)

library(corrplot)
corrplot(cor(d))


#Q1
reg.f = lm(deathrate ~ ., data=d)
summary(reg.f)
reg.f.aic = step(reg.f)
summary(reg.f.aic)


betahat = reg.f$coefficients #MLE for beta
residuals = reg.f$residuals
s2 = summary(reg.f)$sigma^2 # MLE for sigma^2

X = cbind(1, X) # add a column of 1s for the intercept

#Q2a
g = c(.1, 1, 10, 100, 1000)
#mean beta0:
betahat[1] * g/(g+1)
#mean for sigma^2
a = n/2
b = s2/2 + 1/(2*g+2)*t(betahat)%*%t(X)%*%X%*%betahat
b/(a-1)

#Q2b
g = n
q = 2
X0 = X[, -(8:9)]
BF = (g+1)^(+q/2) * 
  ((t(y)%*%y-g/(g+1)*t(y)%*%X0%*%solve(t(X0)%*%X0)%*%t(X0)%*%y)/
  (t(y)%*%y-g/(g+1)*t(y)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%y))^(n/2)
log10(BF) # Strong evidence in favour of H1

# Q3
marglkd = function(gamma, mat=X, g=nrow(mat)){
  q = sum(gamma)
  X1 = mat[, c(T, gamma)]
  if(q==0){return(-(q+1)/2*log(g+1) -n/2*log(t(y)%*%y))}
  m = -q/2*log(g+1) -
    n/2*log(t(y)%*%y - g/(g+1)* t(y)%*% X1 %*%solve(t(X1)%*%X1) %*%
              t(X1)%*%y)
return(m)
}

X1 = X[, 1:4]
marglkd(c(F, F, F), X1)
marglkd(c(F, F, T), X1)
marglkd(c(F, T, F), X1)
marglkd(c(F, T, T), X1)
marglkd(c(T, F, F), X1)
marglkd(c(T, F, T), X1)
marglkd(c(T, T, F), X1)
marglkd(c(T, T, T), X1)


#Q4
# install.packages("BMS")
library(BMS)
reg = bms(d, burn=1e4, iter=5e4)
coef(reg)
image(reg)
topmodels.bma(reg)[, 1:5]


# Q5
# Caution! For pred.density to work, the explained variable Y must be in the first column of the data.frame.
idx = 21:26
d2 = d[-idx, ]
d3 = d[idx, ]

reg2.f = step(lm(deathrate ~ ., data=d2))
pred.f = predict(reg2.f, d3, se.fit=T)


reg2.b = bms(d2, burn=1e4, iter=5e4, nmodel=2e3)
topmodels.bma(reg2.b)[, 1]

reg2.b1 = bms(d2, burn=1e4, iter=5e4, nmodel = 1)

pdens1 = pred.density(reg2.b1, d3)
pdens1
plot(pdens1, 2)

pdens.all = pred.density(reg2.b, d3)
pdens.all
plot(pdens.all, 2)

plot(d3$deathrate, d3$deathrate,col=1, xlim=c(800, 1100), ylim=c(800, 1100))
abline(0, 1)

# frequentist
points(d3$deathrate, pred.f$fit, col=2)
for(i in 1:6){
  lines(c(d3$deathrate[i], d3$deathrate[i]), pred.f$fit[i]+c(-2,2)*pred.f$se.fit[i], col=2)
}


# mixture of models
for(i in 1:6){
  lines(c(d3$deathrate[i], d3$deathrate[i])+6, quantile(pdens.all, c(.025, .975))[i, ], col=4)
}

