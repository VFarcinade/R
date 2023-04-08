# For Beta prior 

prior_beta<- function (theta_new,p_new,theta_current,p_current,a_p,b_p,a_theta,b_theta){
  temp1 = (p_new/p_current)**(a_p-1)*((1-p_new)/(1-p_current))**(b_p-1)
  temp2 = (theta_new/theta_current)**(a_theta-1)*((1-theta_new)/(1-theta_current))**(b_theta-1)
  return((temp1*temp2)*(theta_new>=0)*(theta_new<=1)*(p_new>=0)*(p_new<=1)*(theta_current>=0)*(theta_current<=1)*(p_current>=0)*(p_current<=1))
}

# For uniform prior 
prior_unif<- function (theta_new,p_new,theta_current,p_current){
  return((theta_new>=0)*(theta_new<=1)*(p_new>=0)*(p_new<=1)*(theta_current>=0)*(theta_current<=1)*(p_current>=0)*(p_current<=1))
}


#Compute de conditionnal distribution z | x_obs, theta, p 
pcond_z_<-function(x_obs,p,theta){
  if(x_obs=='AA'){
    prob = p*(1-theta)/(p*(1-theta)+(1-p)*(1-theta)**2)
    return (runif(1)<prob)
  }
  if(x_obs=='aa'){
    prob = p*theta/(p*theta+(1-p)*theta**2)
    return (runif(1)<prob)
  }
  return(0)
}

pcond_z <- Vectorize(pcond_z_, c("x_obs"))


#Compute de conditionnal distribution p,theta |z, x_obs, 
pcond_p_theta <-function(x_obs, z_current,prior,a_p=1,b_p=1,a_theta=1,b_theta=1){
  N = sum(z_current)
  NaA = sum(x_obs=="aA")
  Naa_1 = sum(x_obs=="aa" & z_current == 1 )
  NAA_1 = sum(x_obs=="AA" & z_current == 1 )
  Naa_0 = sum(x_obs=="aa" & z_current == 0 )
  NAA_0 = sum(x_obs=="AA" & z_current == 0 )
 
  
  #with a uniform prior
  if(prior=="unif"){ 
    p = rbeta(1,N+1,500-N+1)
    theta = rbeta(1,Naa_1+2*Naa_0+NaA+1,NAA_1+2*NAA_0+NaA+1)
  }
  #with a beta prior
  else{
    p = rbeta(1,N+a_p,500-N+b_p)
    theta = rbeta(1,Naa_1+2*Naa_0+NaA+a_theta,NAA_1+2*NAA_0+NaA+b_theta)
  }
  return(c(p,theta))
}