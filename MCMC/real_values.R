n=302+125+73
paa=125/n
paA=73/n
pAA=302/n

theta_=paa+0.5*paA
p_=(paa-theta_^2)/(theta_*(1-theta_))

paa_=p_*theta_+(1-p_)*theta_^2
paA_=2*(1-p_)*theta_*(1-theta_)
pAA_=p_*(1-theta_) + (1-p_)*(1-theta_)^2

verif=sum(abs(c(paa,paA,pAA)-c(paa_,paA_,pAA)))

print("Méthode des moments:")
cat("theta_hat: ",theta_,"\n")
cat("p_hat: ",p_,"\n")
