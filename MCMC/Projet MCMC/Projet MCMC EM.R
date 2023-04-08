#Exercice 1

log_likelihood <- function(x,p,theta){
  p_AA <- p*(1-theta) + (1-p)*(1-theta)**2
  p_aA <- 2*(1-p)*theta*(1-theta)
  p_aa <-p*theta+(1-p)*theta**2
  n_AA <- sum(x=="AA")
  n_aa <- sum(x=="aa")
  n_aA <- sum(x=="aA")
  return(n_AA*log(p_AA)+n_aa*log(p_aa)+n_aA*log(p_aA))
}

est_EM <- function(x,n_em,epsilon,p_0, theta_0) {
  
  # initialization
  n <- length(x)
  p <- rep(0,n_em+1)
  theta <- rep(0,n_em+1)
  like <- rep(0,n_em)
  p[1] <- p_0
  theta[1] <- theta_0
  like[1] <- log_likelihood(x,p[1],theta[1])
  
  n_AA <- sum(x=="AA")
  n_aa <- sum(x=="aa")
  n_aA <- sum(x=="aA")
  
  # Maximization
  for (t in 1:n_em){
    
    
    p_AA1 <- p[t]*(1-theta[t])/(p[t]*(1-theta[t])+(1-p[t])*(1-theta[t])**2)
    p_AA0 <- (1-p[t])*(1-theta[t])**2/(p[t]*(1-theta[t])+(1-p[t])*(1-theta[t])**2)
    p_aa1 <- p[t]*theta[t]/(p[t]*theta[t] + (1-p[t])*theta[t]**2)
    p_aa0 <- (1-p[t])*theta[t]**2/(p[t]*theta[t] + (1-p[t])*theta[t]**2)

    
    p[t+1] <- (n_AA*p_AA1 + n_aa*p_aa1)/n
    theta[t+1] <- ((1+p_aa0)*n_aa + n_aA)/(n + n_aa*p_aa0 + n_AA*p_AA0 + n_aA)
    
    like[t+1] <- log_likelihood(x,p[t+1],theta[t+1])
    
    if (abs(like[t+1]-like[t]) < epsilon) {
      n_stop <- t
      break}

  }
  return(list(p=p[1:n_stop], theta=theta[1:n_stop],like=like[1:n_stop]))
}


#vector of observations
xobs<-c(rep("AA",302),rep("aa",125),rep("aA",73))

#initialization of the parameters
p_0 <- 1/2
theta_0 <- 1/2
n_iter <- 10000
epsilon <- 1e-16

out_EM <- est_EM(xobs, n_iter, epsilon, p_0, theta_0)

par(mfrow=c(1,2))
plot(out_EM$p,main="Convergence of p_hat, estimator of p", 
     xlab="number of iterations",ylab="p_hat",panel.first = grid(),pch=20,
     col ='lightblue3',type='o',font.main=4,col.main="gray24",font.lab=3,
     col.lab="gray24")

plot(out_EM$theta,main="Convergence of theta_hat, estimator of theta", 
     xlab="number of iterations",ylab="theta_hat" ,panel.first = grid(),pch=20, 
     col ='lightblue3',type='o',font.main=4,col.main="gray24",font.lab=3,
     col.lab="gray24")



theta_star <- out_EM$theta[length(out_EM$theta)]
p_star <- out_EM$p[length(out_EM$p)]
theta_test <- seq(0,1,0.01)
p_test <- seq(0,1,0.01)

par(mfrow=c(1,2))
plot(theta_test, log_likelihood(xobs,p_star,theta_test),type='l',
     xlab='theta_test',ylab="l(xobs,p*,theta_test)",xlim=c(0,1),
     panel.first = grid(),font.lab=3,col.lab="gray24",col="gray47")
title(main=paste("Evolution of the log-vraisemblance","\n",sep=""),
      font.main=4,col.main="gray24")
title(main=paste("\n","depending on theta_test",sep=""),font.main=4,
      col.main="gray24",)
abline(v=theta_star,col='mediumpurple1',lwd=3, lty=2)
abline(h=log_likelihood(xobs,p_star,theta_star),col='mediumvioletred',
       lwd=3, lty=2)
legend(x="bottomright",legend=c("theta*","l(xobs,p*,theta*)"),lwd = 1,
       col=c("mediumpurple1","mediumvioletred"),cex = 0.8,bg="white", lty=2)

plot(p_test, log_likelihood(xobs,p_test,theta_star),type='l',
     xlab='p_test',ylab="l(xobs,p_test,theta*)",xlim=c(0,1),
     panel.first = grid(),font.lab=3,col.lab="gray24",col="gray47")
title(main=paste("Evolution of the log-vraisemblance","\n",sep=""),font.main=4,
      col.main="gray24")
title(main=paste("\n","depending on p_test",sep=""),font.main=4,
      col.main="gray24")
abline(v=p_star,col='mediumpurple1',lwd=3, lty=2)
abline(h=log_likelihood(xobs,p_star,theta_star),col='mediumvioletred',lwd=3, 
       lty=2)
legend(x="bottomleft",legend=c("p*","l(xobs,p*,theta*)"),
       col=c("mediumpurple1","mediumvioletred"),lwd=1,cex = 0.8,bg="white", 
       lty=2)


