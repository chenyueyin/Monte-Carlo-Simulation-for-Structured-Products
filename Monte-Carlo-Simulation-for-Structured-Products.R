###define parameters
S0<-2614.45
R0<-1512.155
S<-vector()
S[1]<-S0
R<-vector()
R[1]<-R0
sigma1<-0.17
sigma2<-0.19
div1<-0.0189
div2<-0.0133
r0<-0.02337
r1<-0.0268
rho<-0.9
dt<-1/252
#T<-1362/252
T<-1986/365
covm<-matrix(c(dt,rho*dt,rho*dt,dt),nrow = 2,ncol = 2)
library(MASS)
nsim<-100000
V<-rep(0,nsim)
for (k in 1:nsim){
  phi1<-rnorm(1)
  phi2<-rnorm(1)
  S[2]<-S[1]*exp((r1-div1-0.5*sigma1^2)*(1298/252)+sigma1*sqrt(1298/252)*phi1)
  R[2]<-R[1]*exp((r1-div2-0.5*sigma2^2)*(1298/252)+sigma2*rho*sqrt(1298/252)*phi1+sigma2*sqrt(1-rho^2)*sqrt(1298/252)*phi2)
  dW<-mvrnorm(64,mu=c(0,0),Sigma = covm)
  for (i in 3:66){
    S[i]<-S[i-1]+(r1-div1)*S[i-1]*dt+sigma1*S[i-1]*dW[(i-2),1]
    R[i]<-R[i-1]+(r1-div2)*R[i-1]*dt+sigma2*R[i-1]*dW[(i-2),2]
  }
  Smean<-mean(S[2:66])
  Rmean<-mean(R[2:66])
  if (min(Smean/S0,Rmean/R0)>=1.21){
    payoff=min((1000+1000*(min(Smean/S0,Rmean/R0)-1.21)*3.34+415),2116.4)
  } else if ((min(Smean/S0,Rmean/R0)<1.21)&(min(Smean/S0,Rmean/R0)>=1)){
    payoff=1000+1000*(min(Smean/S0,Rmean/R0)-1)*1.5+100
  } else if ((min(Smean/S0,Rmean/R0)<1)&(min(Smean/S0,Rmean/R0)>=0.95)){
    payoff=1000+1000*(min(Smean/S0,Rmean/R0)-0.95)*2
  } else if (min(Smean/S0,Rmean/R0)<0.95){
    payoff=1000*min(Smean/S0,Rmean/R0)+50
  }
  V[k]<-payoff*exp(-r1*T)
}
V0<-mean(V)
V0