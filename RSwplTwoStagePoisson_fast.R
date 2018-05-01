
#setting: notation
N1=30 ## number of sampling cluster in the first stage (population level) 
N2=40 ##number of elements in each sampling cluster (population level)
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)
population$PSU<-population$long
overlap=ceiling(N2*3/4)

model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}

population<-model_cluster(population, overlap)
T=length(unique(population$cluster))


##check
population<-model_cluster(population, overlap)
table(population$cluster==population$PSU)
table(table(population$cluster))
table(table(population$PSU))

#Model: parameter from random intercept model
truebeta1=1
truebeta2=3
truesigma2=2
truetau2_11=1
truetau_12=0.4
truetau2_22=1
PairCov<-matrix(c(truetau2_11, truetau_12, truetau_12, truetau2_22), nrow=2, byrow=T)
###check positive definite 
#install.packages("matrixcalc")
library("matrixcalc")
is.positive.definite(PairCov)


truevalue<-c(truebeta1,truebeta2, truesigma2, truetau2_11, truetau_12, truetau2_22)
names(truevalue)<-c("beta1", "beta2", "sigma2", "tau2_11", "tau_12", "tau2_22")

##Population data
#install.packages("MASS")
library("MASS")
#install.packages("rockchalk")
library(rockchalk)
re=mvrnorm(n=T, mu = c(0,0), Sigma = PairCov) #generate vector of random effect (a, b)
population$a<-re[,1][population$cluster]
population$b<-re[,2][population$cluster]


population$x<-rnorm(N1*N2)+rnorm(T)[population$cluster]
population$y<-with(population, truebeta1+a+truebeta2*x+b*x+rnorm(N1*N2,s=sqrt(truesigma2)))
population$r=with(population, x*(y-truebeta1-truebeta2*x))
population$ID_unit=with(population, 1:(N1*N2))

#uninformative two-stage sampling design (first-stage: Poisson, Second-stage:SRSWOR)
#install.packages("sampling")
library("sampling")

##uninformative Poisson first-stage 
pi1=runif(N1) ## sampling inclusion probability for sampling cluster
FirststagePoisson=UPpoisson(pi1)
n1=sum(FirststagePoisson) # number of PSU in the first-stage sample

FirststagePoissonSample=subset(population, population$PSU%in% which(FirststagePoisson==1))

#uninformative SRSWOR second-stage  
n2=ceiling(N2/10)##number of elements in each sampling cluster ( sample level)
SecondstageSRSWOR=unlist(lapply(rep(n2,n1), function(v) return(srswor(v, N2))))
TwostagePoissonSample<-FirststagePoissonSample[c(which(SecondstageSRSWOR==1)), ]


#informative two-stage sampling design (first-stage: Poisson, Second-stage:SRSWOR)

#informative Poisson first-stage sample
param1=c(0.5, 8)
pi1informative= function(r, sc,  param, N1){
   a=rep(NA, N1)
   b=rep(NA, N1)
   for (i in 1: N1){
      a[i]=mean(r[sc==i])
      b[i]=(param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))
   }
   b}
pi1is=pi1informative(population$r,population$PSU, param1 ,N1)
FirststagePoissonis=UPpoisson(pi1is)
n1is=sum(FirststagePoissonis) # number of PSU in the first-stage sample

FirststagePoissonSampleis=subset(population, population$PSU%in% which(FirststagePoissonis==1))


#informative second-stage sample(SRSWOR)
param2=c(0.05, 3.5)
n2informative= function(r, sc, param, N2){
   a=rep(NA, length=length(unique(population$sc)))
   b=rep(NA, length=length(unique(population$sc)))
   for (i in unique(sc)){
      a[i]=mean(r[sc==i])
      b[i]=2*ceiling((param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))*N2/2)
   }
   b
}

n2pop=n2informative(population$r,population$PSU, param2 ,N2)
n2is=n2pop*FirststagePoissonis
SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
TwostagePoissonSampleis=FirststagePoissonSampleis[c(which(SecondstageSRSWORis==1)), ]

#Estimation: full-likelihood
#install.packages("lme4")
library(lme4)

##Census estimator 
fit_NML=lmer(y~(1+x|cluster)+x,data=population)
###summary(fit_NML)

### get the fixed effect for alpha and beta 
fixef(fit_NML)

### get the sigma2 for the error term 
sigma(fit_NML)^2

### get the variance for the random effect \tau2_11, \tau_12, \tau2_22  
VarCorr(fit_NML)$cluster[c(1, 2,4)] # this is the variance not standard deviation

##uninformative two-stage sampling design (Poisson)
lmer(y~(1+x|cluster)+x,data=TwostagePoissonSample)

##informative two-stage sampling design (Poisson)
lmer(y~(1+x|cluster)+x,data=TwostagePoissonSampleis)

# Estimation: pairwise likelihood (without weight)
l2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+(x1^2)*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+(x2^2)*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   -log(det)/2-(1/2)*(1/det)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)
   
}	


dalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dr1<- -1
   dr2<- -1
   

   (-1/2)*(1/det)*(2*r1*dr1*pc22-2*dr1*r2*pc12-2*r1*dr2*pc12+2*r2*dr2*pc11 )
   
}

dbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau_12, tau2_22){
   
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dr1<- -x1
   dr2<- -x2
   
   (-1/2)*(1/det)*(2*r1*dr1*pc22-2*dr1*r2*pc12-2*r1*dr2*pc12+2*r2*dr2*pc11)
}	

dsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-1
   dpc22<-1
   dpc12<-0
   
   ddet<-dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   (-1/2)*(ddet/det)-1/2*(-ddet)/(det)^2*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
}

dtau2_11<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-1
   dpc22<-1
   dpc12<-ifelse(g1==g2, 1, 0)
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   
   (-1/2)*(ddet/det)-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
 
}	


dtau_12<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-2*x1
   dpc22<-2*x2
   dpc12<-ifelse(g1==g2,x1+x2, 0)
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
 
}

dtau2_22<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22){
   pc11<-tau2_11+2*x1*tau_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau_12+x2*tau_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-x1^2
   dpc22<-x2^2
   dpc12<-ifelse(g1==g2, x1*x2, 0 )
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
  
}




#optimization problem for PL (without weight)
fit_PL<-function(y,g,x, pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                   tau2_22=exp(theta[6]))
      sum(increment)/T
   }
   gr<-function(theta){
      incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                         tau2_22=exp(theta[6]))
      incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                        sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                        tau2_22=exp(theta[6]))
      incrementds=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                        alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                        tau2_11=exp(theta[4]), tau_12=theta[5],
                                        tau2_22=exp(theta[6]))
      incrementdt_11=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                      sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],
                                      tau2_22=exp(theta[6]))
      incrementdt_12=dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],
                                            tau2_22=exp(theta[6]))
      incrementdt_22=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],
                                            tau2_22=exp(theta[6]))
      c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt_11),sum(incrementdt_12),
        sum(incrementdt_22))/T
   }
   optim(pars,func1, gr, method="BFGS",control=list(fnscale=-1,parscale=c(1/n,1/n,1/n,1/n, 1/n, 1/n)))
}


###uninformative (without weight)
estimator_PL<-fit_PL(TwostagePoissonSample$y, TwostagePoissonSample$cluster, TwostagePoissonSample$x, pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                                                                                             log(truevalue[6])))
estimator_PL

###informative sampling (without weight)
estimatoris_PL<- fit_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                                                                                                      log(truevalue[6])))
estimatoris_PL


##Define the pairwise score function and checking the pairwise score at PML (without weight)
pairscore_PL<-function(y,g,x, theta){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                      sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                      tau2_22=exp(theta[6]))
   incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                     sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                     tau2_22=exp(theta[6]))
   incrementds=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                     alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                                     tau2_22=exp(theta[6]))
   incrementdt_11=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                   sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                                   tau2_22=exp(theta[6]))
   incrementdt_12=dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                         sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                                         tau2_22=exp(theta[6]))
   incrementdt_22=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                         sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                                         tau2_22=exp(theta[6]))
   c(sum(incrementda), sum(incrementdb), sum(incrementds),sum(incrementdt_11), sum(incrementdt_12),
     sum(incrementdt_22))/T
}



###uninformative
pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x, estimator_PL[[1]])


###informative sampling
pairscore_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, estimatoris_PL[[1]])




dyn.load("RSFourOrdPiTwostagePoisson.so")
# second-order inclusion probability
C2<-function(pos1, pos2,sc1, sc2,fss, n2infor,N2){
   .C("SecOrdPi",as.integer(pos1), as.integer(pos2),as.integer(sc1), as.integer(sc2), as.double(fss), as.double(n2infor),as.double(N2),length(pos1),rval=numeric(length(pos1)))$rval
}



# fourth-order inclusion probability
C4<-function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor,N2){
   .C("FourOrdPi",as.integer(pos1), as.integer(pos2),as.integer(pos3), as.integer(pos4),as.integer(sc1), as.integer(sc2),
      as.integer(sc3), as.integer(sc4), as.double(fss), as.double(n2infor),as.double(N2),length(pos1),rval=numeric(length(pos1)))$rval	
   
}


##Define the second-order inclusion probability
SecOrdPi<-function(pos1, pos2,sc1, sc2,fss, n2infor,N2){
   #pi<-SecOrdPiInternal(pos1, pos2,sc1, sc2,n1, N1, n2infor,N2)	
   Cpi<-C2(pos1, pos2,sc1, sc2,fss, n2infor,N2)
   #if ((pi-Cpi)/(pi+Cpi)>1e-10) stop(paste(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,":",pi,Cpi,sep=","))
   Cpi
}


##Define the  fourth-order inclusion probability
FouOrdPi<-function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor,N2){
   #pi<-FouOrdPiInternal(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,n2infor,N2)	
   Cpi<-C4(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor,N2)
   #if ((pi-Cpi)/(pi+Cpi)>1e-10) stop(paste(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,":",pi,Cpi,sep=","))
   Cpi
}

#Define the fourth-order Delta
FouOrdDel=function(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4,fss, n2infor, N2){
   FouOrdPi(pos1, pos2,pos3, pos4,sc1, sc2,sc3, sc4, fss, n2infor,N2)-
      SecOrdPi(pos1, pos2,sc1, sc2, fss, n2infor,N2)*
      SecOrdPi(pos3, pos4,sc3, sc4, fss, n2infor,N2)
}

# Estimation: weighted pairwise likeliood 
wl2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2, fss,  n2infor,N2){
   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor, N2)*
      (l2(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
}	


#wdalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2,fss,  n2infor,N2){
#   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor,N2)*
#      (dalpha(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
#}	

#wdbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2,fss,  n2infor,N2){
#   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor, N2)*
#      (dbeta(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
#}	


#wdsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2, fss,  n2infor,N2){
#   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor,N2)*
#      (dsigma2(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
#}	

#wdtau2_11<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2, fss,  n2infor,N2) {
#   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor,N2)*
#      (dtau2_11(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
#}

#wdtau_12<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2, fss,  n2infor,N2) {
#   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor,N2)*
#      (dtau_12(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
#}

#wdtau2_22<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau_12, tau2_22, pos1, pos2,sc1, sc2, fss,  n2infor,N2) {
#   1/SecOrdPi(pos1, pos2,sc1, sc2, fss,  n2infor,N2)*
#      (dtau2_22(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau_12, tau2_22))
#}

#optimization (WPL)
fit_WPL<-function(y,g,x, pos, sc,fss,  n2infor, N2,  pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      wincrement=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                     sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]),  pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
      sum(wincrement)/T
   }
   gr<-function(theta){
      wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j],fss, n2infor,N2)
      wincrementda=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],
                           tau2_22=exp(theta[6]))
      wincrementdb=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                          sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
      wincrementds=exp(theta[3])*wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                          alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
      wincrementdt_11=exp(theta[4])*wij*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                        sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
      wincrementdt_12=wij*dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
      wincrementdt_22=exp(theta[6])*wij*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
      c(sum(wincrementda), sum(wincrementdb), sum(wincrementds), sum(wincrementdt_11), sum(wincrementdt_12), 
        sum(wincrementdt_22))/T
   }
   optim(pars,func1,gr,  method="BFGS",
         control=list(fnscale=-1,parscale=c(1/n,1/n,1/n,1/n,1/n, 1/n)))
}

##uninformative sampling WPL (with weight)
estimator_WPL=fit_WPL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                      pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2, 
                      pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                             log(truevalue[6])))
estimator_WPL
##informative sampling WPL (with weight)
estimatoris_WPL=fit_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                        TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2,  
                        pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], log(truevalue[6])))
estimatoris_WPL
##Define the  weighted pairwise score function for WPL and check the value of weighted pairwise score function at WPML
pairscore_WPL<-function(y,g,x, theta, pos, sc,fss,  n2infor, N2){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j],fss, n2infor,N2)
   wincrementda=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                        sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
   wincrementdb=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                       sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   wincrementds=exp(theta[3])*wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                       alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   wincrementdt_11=exp(theta[4])*wij*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   wincrementdt_12=wij*dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   wincrementdt_22=exp(theta[6])*wij*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                           sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]))
   c(sum(wincrementda), sum(wincrementdb), sum(wincrementds), sum(wincrementdt_11), sum(wincrementdt_12), 
     sum(wincrementdt_22))/T
}

##uninformative sampling WPL (with weight)
estimator_WPL=fit_WPL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                      pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2, 
                      pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], log(truevalue[6])))
estimator_WPL
##informative sampling WPL (with weight)
estimatoris_WPL=fit_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                        TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2,  
                        pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], log(truevalue[6])))
estimatoris_WPL
#variance estimation for PL under stratified SRSWORS
#define the pairwise likelihood (without weight)
#install.packages("numDeriv")
library("numDeriv")

#uninformative 
#Calculate Hessian matrix H for PL (bread for uninformative sampling design)
#pl=function(theta,y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x){
#   n<-length(y)
#   ij=expand.grid(1:n,1:n)
#   ij<-ij[ij[,1]<ij[,2],]
#   ij<-ij[g[ij[,1]]==g[ij[,2]],]
#   i<-ij[,1]
#   j<-ij[,2]
#   increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
#                sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
#   sum(increment)/T
#}
#estH_PL=hessian(pl, estimator_PL[[1]])

estH_PL=-jacobian(function(theta){with(TwostagePoissonSample,
                                         pairscore_PL(y,cluster,x,theta))}, x=estimator_PL[[1]],method="simple")


#Calculate  variance matrix J  for PL (meat for uninformative sampling design)
fast_J_PL<-function(y,g,x,pos, sc,fss, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
         cat(".")
         incrementdaij=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                              tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdbij=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                             tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdsij=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                                             sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_11ij=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                 tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_12ij=dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                 tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_22ij=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                                 tau2_11=exp(theta[4]),tau_12=theta[5], tau2_22=exp(theta[6]) )
         psij=c(incrementdaij, incrementdbij, incrementdsij, incrementdt_11ij, incrementdt_12ij, incrementdt_22ij)
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
         jj<-rep(j,length(k))
         incrementdakl=dalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdbkl=dbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                             sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdskl=exp(theta[3])*dsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                             sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_11kl=exp(theta[4])*dtau2_11(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_12kl=dtau_12(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]))
         incrementdt_22kl=exp(theta[6])*dtau2_22(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]))
         pskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdt_11kl,incrementdt_12kl, incrementdt_22kl )
         sumpskl<-colSums(FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],fss,n2infor,N2)* pskl)
         psijkl<-tcrossprod(psij,sumpskl)
         sum=sum+psijkl
      }
   }
   rval<-sum/(T^2)
   ##attr(rval, "pairs")<-keep ##debug
   rval
}

estJ_PL=fast_J_PL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                  pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU,fss=pi1,
                  n2infor=FirststagePoisson*n2, N2, theta=estimator_PL[[1]] )

#sanwich estimator (uninformative sampling )
sanestimator_PL= solve(estH_PL)%*% estJ_PL%*% solve(t(estH_PL))

#Informative
#Calculate Hessian matrix H for PL (bread for informative sampling design)
#plis=function (theta, y=TwostagePoissonSampleis$y, g=TwostagePoissonSampleis$cluster, x=TwostagePoissonSampleis$x){
#   n<-length(y)
#   ij=expand.grid(1:n,1:n)
#   ij<-ij[ij[,1]<ij[,2],]
#   ij<-ij[g[ij[,1]]==g[ij[,2]],]
#   i<-ij[,1]
#   j<-ij[,2]
#   increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
#                sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]) )
#   sum(increment)/T
#}
#estHis_PL=hessian(plis, estimatoris_PL[[1]])


estHis_PL=-jacobian(function(theta){with(TwostagePoissonSampleis,
                                         pairscore_PL(y,cluster,x,theta))}, x=estimatoris_PL[[1]],method="simple")

#Calculate  variance matrix J  for PL (meat for informative sampling design)
estJis_PL=fast_J_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,pos=TwostagePoissonSampleis$ID_unit,  
                    sc=TwostagePoissonSampleis$PSU, fss=pi1is, n2infor=n2is, N2, theta=estimatoris_PL[[1]] )
#sanwich estimator (informative sampling )
sanestimatoris_PL = solve(estHis_PL)%*% estJis_PL%*% t(solve(estHis_PL))


#variance estimation for WPL under two-stage SRSWORS
##define H as in page 96 of my thesis as \hat{H}(\est)
#define weighted pairwise likelihood WPL 

##uninformative sampling
wpl=function (theta, y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster,x=TwostagePoissonSample$x,
              pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, 
              n2infor=FirststagePoisson*n2 , N2=length(unique(population$lat)) ){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
   sum(increment)/T
}
estH_WPL=hessian(wpl, estimator_WPL[[1]])


##define \hat{J}(\theta) as in page 97 of my thesis and  evaluate at the WPLE
fast_J_WPL<-function(y,g,x,  pos,  sc, fss, n2infor,N2, theta){
   n<-length(y)
   sum=0
   
   kl<-expand.grid(1:n,1:n)
   kl<-kl[kl[,1]<kl[,2],]
   kl<-kl[g[kl[,1]]==g[kl[,2,]],]
   k<-kl[,1]
   l<-kl[,2]
   
   for (i in 1:(n-1)){
      cat(i)
      js <- (i+1):n
      js <- js[g[js] %in% g[i]]
      for(j in js){
         cat(".")
         wij<-1/SecOrdPi(pos[i], pos[j],sc[i], sc[j],fss, n2infor,N2)
         incrementdaij=wij*dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                               tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdbij=wij*dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                              tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdsij=exp(theta[3])*wij*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdt_11ij=exp(theta[4])*wij*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                               tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdt_12ij=wij*dtau_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                               tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdt_22ij=exp(theta[6])*wij*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),
                                               tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         
         wpsij=c(incrementdaij, incrementdbij, incrementdsij, incrementdt_11ij, incrementdt_12ij, incrementdt_22ij)
         
         
         ## k,l vectorised: probably can't afford memory to do that for ijkl 
         ii <-rep(i, length(k))
         jj<-rep(j,length(k))
         wkl<-1/SecOrdPi(pos[k], pos[l],sc[k], sc[l],fss, n2infor,N2)
         incrementdakl=wkl*dalpha(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                               sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdbkl=wkl*dbeta(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                              sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdskl=exp(theta[3])*wkl*dsigma2(y[k],y[l],g[k],g[l],x[k],x[l],alpha=theta[1],beta=theta[2],
                                              sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdt_11kl=exp(theta[4])*wkl*dtau2_11(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdt_12kl=wkl*dtau_12(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         incrementdt_22kl=exp(theta[6])*wkl*dtau2_22(y[k],y[l],g[k],g[l],x[k],x[l], alpha=theta[1],beta=theta[2],
                                                  sigma2=exp(theta[3]), tau2_11=exp(theta[4]), tau_12=theta[5],  tau2_22=exp(theta[6]))
         wpskl=cbind(incrementdakl, incrementdbkl, incrementdskl, incrementdt_11kl, incrementdt_12kl, incrementdt_22kl)
         sumwpskl<-colSums( (1/FouOrdPi( pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],fss,  n2infor,N2))*FouOrdDel(pos[ii], pos[jj], pos[k], pos[l], sc[ii], sc[jj], sc[k], sc[l],fss,  n2infor,N2)* wpskl)
         wpsijkl<-tcrossprod(wpsij,sumwpskl)
         sum=sum+wpsijkl
      }
   }
   rval<-sum/(T^2)
   # attr(rval, "pairs")<-keep ##debug
   rval
}
estJ_WPL=fast_J_WPL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster,
                    x=TwostagePoissonSample$x, pos=TwostagePoissonSample$ID_unit,  sc=TwostagePoissonSample$PSU, fss=pi1, 
                    n2infor= FirststagePoisson*n2, N2, theta=estimator_WPL[[1]] )

estJ_WPL
# sanwich estimator H^{-1}J (H^{-1})^\T
##uninformaitve
sanestimator_WPL= solve(estH_WPL)%*% estJ_WPL%*% t(solve(estH_WPL))
sanestimator_WPL

##informative sampling
wplis=function (theta, y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,
                pos=TwostagePoissonSampleis$ID_unit, sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                n2infor=n2is , N2=length(unique(population$lat)) ){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                 sigma2=exp(theta[3]),tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]),  pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
   sum(increment)/T
}
estHis_WPL=hessian(wplis, estimatoris_WPL[[1]])


estJis_WPL=fast_J_WPL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,
                      x=TwostagePoissonSampleis$x, pos=TwostagePoissonSampleis$ID_unit,  sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                      n2infor= n2is, N2, theta=estimatoris_WPL[[1]] )
# sanwich estimator H^{-1}J (H^{-1})^\T
##informaitve
sanestimatoris_WPL= solve(estHis_WPL)%*% estJis_WPL%*% t(solve(estHis_WPL))
sanestimatoris_WPL

#simulation

LOTS=400
#Fit from NML,PL, WPL for uninformative sampling
Fit_NML<-matrix(0,nrow=LOTS,ncol=6)
Fit_PL<-matrix(0,nrow=LOTS,ncol=6)
Fit_WPL<-matrix(0,nrow=LOTS, ncol=6)

#Fit from NML, PL, WPL for informative sampling
Fitis_NML<-matrix(0,nrow=LOTS,ncol=6)
Fitis_PL<-matrix(0,nrow=LOTS,ncol=6)
Fitis_WPL<-matrix(0,nrow=LOTS, ncol=6)

#Hessian matrix for PL(without weight) for uninformative sampling
H_PL<-array(0, c(6,6, LOTS))

#Hessian matrix for PL(without weight) for informative sampling
His_PL<-array(0, c(6,6, LOTS))

#Hessian matrix for WPL for uninformative sampling
H_WPL<-array(0, c(6,6, LOTS))

#Hessian matrix for WPL for informative sampling
His_WPL<-array(0, c(6,6, LOTS))

#Variance matrix J for PL for uninformative sampling
J_PL<-array(0, c(6,6, LOTS))

#Variance matrix J for PL for informative sampling
Jis_PL<-array(0, c(6,6, LOTS))

#Variance matrix J for WPL for uninformative sampling
J_WPL<-array(0, c(6,6, LOTS))

#Variance matrix J for WPL for informative sampling
Jis_WPL<-array(0, c(6,6, LOTS))

#Sanwich variance estimator for PL for uninformative sampling
G_PL<-array(0, c(6,6, LOTS))

#Sanwich variance estimator for PL for informative sampling
Gis_PL<-array(0, c(6,6, LOTS))

#Sanwich variance estimator for WPL for uninformative sampling
G_WPL<-array(0, c(6,6, LOTS))

#Sanwich variance estimator for WPL for informative sampling
Gis_WPL<-array(0, c(6,6, LOTS))

#Pairwise score function for PL for informative sampling 
PS_PL<-matrix(0,nrow=LOTS,ncol=6)

#Pairwise score function for PL for  informative sampling
PSis_PL<-matrix(0,nrow=LOTS,ncol=6)

#Pairwise score function for WPL for informative sampling 
PS_WPL<-matrix(0,nrow=LOTS,ncol=6)

#Pairwise score function for WPL for  informative sampling
PSis_WPL<-matrix(0,nrow=LOTS,ncol=6)



##Estimation: NML, PL and WPL 
for(i in 1:LOTS){
   
   cat(i)
   ##Population data
   re=mvrnorm(n=T, mu = c(0,0), Sigma = PairCov) #generate vector of random effect (a, b)
   population$a<-re[,1][population$cluster]
   population$b<-re[,2][population$cluster]
   
   
   population$x<-rnorm(N1*N2)+rnorm(T)[population$cluster]
   population$y<-with(population, truebeta1+a+truebeta2*x+b*x+rnorm(N1*N2,s=sqrt(truesigma2)))
   population$r=with(population, x*(y-truebeta1-truebeta2*x))
   population$ID_unit=with(population, 1:(N1*N2))
   
   
   
   
   ##uninformative Poisson first-stage 
   pi1=runif(N1) ## sampling inclusion probability for sampling cluster
   FirststagePoisson=UPpoisson(pi1)
   n1=sum(FirststagePoisson) # number of PSU in the first-stage sample
   
   FirststagePoissonSample=subset(population, population$PSU%in% which(FirststagePoisson==1))
   
   #uninformative SRSWOR second-stage  
   n2=ceiling(N2/10)##number of elements in each sampling cluster ( sample level)
   SecondstageSRSWOR=unlist(lapply(rep(n2,n1), function(v) return(srswor(v, N2))))
   TwostagePoissonSample<-FirststagePoissonSample[c(which(SecondstageSRSWOR==1)), ]
   
   
   
   
   #informative two-stage sampling design (first-stage: Poisson, Second-stage:SRSWOR)
   
   #informative Poisson first-stage sample
   param1=c(0.5, 8)
   pi1informative= function(r, sc,  param, N1){
      a=rep(NA, N1)
      b=rep(NA, N1)
      for (i in 1: N1){
         a[i]=mean(r[sc==i])
         b[i]=(param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))
      }
      b}
   pi1is=pi1informative(population$r,population$PSU, param1 ,N1)
   FirststagePoissonis=UPpoisson(pi1is)
   n1is=sum(FirststagePoissonis) # number of PSU in the first-stage sample
   
   FirststagePoissonSampleis=subset(population, population$PSU%in% which(FirststagePoissonis==1))
   
   
   #informative second-stage sample(SRSWOR)
   param2=c(0.05, 3.5)
   n2informative= function(r, sc, param, N2){
      a=rep(NA, length=length(unique(population$sc)))
      b=rep(NA, length=length(unique(population$sc)))
      for (i in unique(sc)){
         a[i]=mean(r[sc==i])
         b[i]=2*ceiling((param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))*N2/2)
      }
      b
   }
   
   n2pop=n2informative(population$r,population$PSU, param2 ,N2)
   n2is=n2pop*FirststagePoissonis
   SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
   TwostagePoissonSampleis=FirststagePoissonSampleis[c(which(SecondstageSRSWORis==1)), ]
   
   
   #NML, PL and WPL (uninformative sampling)
   ra<-lmer(y~(1+x|cluster)+x,data=TwostagePoissonSample)
   rb<-fit_PL(TwostagePoissonSample$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                                                                                          log(truevalue[6])))
   rc<-fit_WPL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
               pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2, 
               pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], log(truevalue[6])))
   
   #NML, PL and WPL (informative sampling)
   rais<-lmer(y~(1|cluster)+x,data=TwostagePoissonSampleis)
   rbis<-fit_PL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster, TwostagePoissonSampleis$x, pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                 log(truevalue[6])))
   
   rcis<-fit_WPL(TwostagePoissonSampleis$y, TwostagePoissonSampleis$cluster,TwostagePoissonSampleis$x,
                 TwostagePoissonSampleis$ID_unit, TwostagePoissonSampleis$PSU, fss=pi1is,  n2infor=n2is , N2,  
                 pars=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                        log(truevalue[6])))
   
   #NML (uniformative sampling)
   Fit_NML[i,1:2]<-fixef(ra)-truevalue[1:2]
   Fit_NML[i,3]<-sigma(ra)^2-truevalue[3]
   Fit_NML[i,4:6]<-VarCorr(ra)$cluster[c(1, 2,4)] -truevalue[4:6]
   
   #PL (uninformative sampling)
   Fit_PL[i,1:2]<-rb$par[1:2]-truevalue[1:2]
   Fit_PL[i,3:4]<-exp(rb$par[3:4])-truevalue[3:4] #reparametrize by taking exponential 
   Fit_PL[i,5]<-rb$par[5]-truevalue[5]
   Fit_PL[i,6]<-exp(rb$par[6])-truevalue[6] #reparametrize by taking exponential 
   
   
   
   #WPL (uniformative sampling)
   Fit_WPL[i,1:2]<-rc$par[1:2]-truevalue[1:2]
   Fit_WPL[i,3:4]<-exp(rc$par[3:4])-truevalue[3:4] #reparametrize by taking exponential 
   Fit_WPL[i,5]<-rc$par[5]-truevalue[6]
   Fit_WPL[i,6]<-exp(rc$par[6])-truevalue[6]
   
   #NML (informative sampling)
   Fitis_NML[i,1:2]<-fixef(rais)-truevalue[1:2]
   Fitis_NML[i,3]<-sigma(rais)^2-truevalue[3]
   Fitis_NML[i,4:6]<-VarCorr(ra)$cluster[c(1, 2,4)] -truevalue[4:6]
   
   #PL (informative sampling)
   Fitis_PL[i,1:2]<-rbis$par[1:2]-truevalue[1:2]
   Fitis_PL[i,3:4]<-exp(rbis$par[3:4])-truevalue[3:4] #reparametrize by taking exponential 
   Fitis_PL[i,5]<-rbis$par[5]-truevalue[5]
   Fitis_PL[i,6]<-exp(rbis$par[6])-truevalue[6]
   
   #WPLE (informative sampling)
   Fitis_WPL[i,1:2]<-rcis$par[1:2]-truevalue[1:2]
   Fitis_WPL[i,3:4]<-exp(rcis$par[3:4])-truevalue[3:4] #reparametrize by taking exponential 
   Fitis_WPL[i,5]<-rcis$par[5]-truevalue[5]
   Fitis_WPL[i,6]<-exp(rcis$par[6])-truevalue[6]
   
   #Calculate Hessian matrix H for PL (bread for uninformative sampling design)
   pl=function(theta,y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x){
      n<-length(y)
      ij=expand.grid(1:n,1:n)
      ij<-ij[ij[,1]<ij[,2],]
      ij<-ij[g[ij[,1]]==g[ij[,2]],]
      i<-ij[,1]
      j<-ij[,2]
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]))
      sum(increment)/T
   }
   H_PL[,,i]=hessian(pl, rb[[1]])
   
   #Calculate  variance matrix J  for PL (meat for uniformative sampling design)
   J_PL[, , i]=fast_J_PL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                         pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU,fss=pi1,
                         n2infor=FirststagePoisson*n2, N2, theta=rb[[1]] )
   
   #sanwich estimator (uninformative sampling )
   G_PL[, ,i] = tryCatch(solve(H_PL[,,i])%*% J_PL[, , i]%*% t(solve(H_PL[,,i])),   error=function(e) matrix(NaN, 6,6))
   
   #Pairwise score function PL (uninformative sampling) 
   PS_PL[i, ]<- pairscore_PL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster, x=TwostagePoissonSample$x,
                             theta=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                          log(truevalue[6])))
   
   #Calculate Hessian matrix H for PL (bread for informative sampling design)
   plis=function (theta, y=TwostagePoissonSampleis$y, g=TwostagePoissonSampleis$cluster, x=TwostagePoissonSampleis$x){
      n<-length(y)
      ij=expand.grid(1:n,1:n)
      ij<-ij[ij[,1]<ij[,2],]
      ij<-ij[g[ij[,1]]==g[ij[,2]],]
      i<-ij[,1]
      j<-ij[,2]
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]) )
      sum(increment)/T
   }
   His_PL[,,i]=hessian(plis, rbis[[1]])
   
   #Calculate  variance matrix J  for PL (meat for informative sampling design)
   Jis_PL[, , i]=fast_J_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,pos=TwostagePoissonSampleis$ID_unit,  
                           sc=TwostagePoissonSampleis$PSU, fss=pi1is,    n2infor=n2is, N2, theta=rbis[[1]] )
   
   #sanwich estimator (informative sampling )
   Gis_PL[, ,i] = tryCatch(solve(His_PL[,,i])%*% Jis_PL[, , i]%*% t(solve(His_PL[,,i])),  error=function(e) matrix(NaN, 6,6))
   
   #Pairwise score function PL (informative sampling)
   PSis_PL[i, ]<- pairscore_PL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster, x=TwostagePoissonSampleis$x,
                               theta=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                       log(truevalue[6])))
   
   #Calculate Hessian matrix H for WPL (bread for uninformative sampling design)
   wpl=function (theta, y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster,x=TwostagePoissonSample$x,
                 pos=TwostagePoissonSample$ID_unit, sc=TwostagePoissonSample$PSU, fss=pi1, 
                 n2infor=FirststagePoisson*n2 , N2=length(unique(population$lat)) ){
      n<-length(y)
      ij=expand.grid(1:n,1:n)
      ij<-ij[ij[,1]<ij[,2],]
      ij<-ij[g[ij[,1]]==g[ij[,2]],]
      i<-ij[,1]
      j<-ij[,2]
      increment=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                    sigma2=exp(theta[3]),tau2_11=exp(theta[4]), tau_12=theta[5], tau2_22=exp(theta[6]), pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
      sum(increment)/T
   }
   H_WPL[,,i]=hessian(wpl, rc[[1]])
   
   #Calculate  variance matrix J  for WPL (meat for uniformative sampling design)
   J_WPL[, , i]=fast_J_WPL(y=TwostagePoissonSample$y,g=TwostagePoissonSample$cluster,
                           x=TwostagePoissonSample$x, pos=TwostagePoissonSample$ID_unit,  sc=TwostagePoissonSample$PSU, fss=pi1, 
                           n2infor= FirststagePoisson*n2, N2, theta=rc[[1]])
   
   #sanwich estimator (uninformative sampling )
   G_WPL[, ,i] =  tryCatch(solve(H_WPL[,,i])%*% J_WPL[, , i]%*% t(solve(H_WPL[,,i])),  error=function(e) matrix(NaN, 6,6))
   
   #Pairwise score function WPL (uninformative sampling)
   PS_WPL[i, ]<- pairscore_WPL(y=TwostagePoissonSample$y, g=TwostagePoissonSample$cluster,x= TwostagePoissonSample$x,
                               theta=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                       log(truevalue[6])),   TwostagePoissonSample$ID_unit, TwostagePoissonSample$PSU, fss=pi1, n2infor=FirststagePoisson*n2, N2)
   
   #Calculate Hessian matrix H  for WPL (bread for informative sampling design)
   ##informative sampling
   wplis=function (theta, y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,
                   pos=TwostagePoissonSampleis$ID_unit, sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                   n2infor=n2is , N2=length(unique(population$lat)) ){
      n<-length(y)
      ij=expand.grid(1:n,1:n)
      ij<-ij[ij[,1]<ij[,2],]
      ij<-ij[g[ij[,1]]==g[ij[,2]],]
      i<-ij[,1]
      j<-ij[,2]
      increment=wl2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                    sigma2=exp(theta[3]),tau2_11=exp(theta[4]),  tau_12=theta[5], tau2_22=exp(theta[6]),  pos[i], pos[j], sc[i], sc[j], fss,  n2infor,N2)
      sum(increment)/T
   }
   His_WPL[, , i]=hessian(wplis, rcis[[1]])
   
   #Calculate Variance matrix J  for WPL (meat for  informative sampling design)
   Jis_WPL[, , i]=fast_J_WPL(y=TwostagePoissonSampleis$y,g=TwostagePoissonSampleis$cluster,
                             x=TwostagePoissonSampleis$x, pos=TwostagePoissonSampleis$ID_unit,  sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                             n2infor= n2is, N2,  theta=rcis[[1]] )
   
   #sanwich estimator for WPL (informative sampling )
   Gis_WPL[,,i]= tryCatch(solve(His_WPL[, , i])%*% Jis_WPL[, , i]%*% t(solve(His_WPL[, , i])),  error=function(e) matrix(NaN, 6,6))
   
   #Pairwise score function WPL (informative sampling)
   PSis_WPL[i, ]<- pairscore_WPL(y=TwostagePoissonSampleis$y, g=TwostagePoissonSampleis$cluster,x=TwostagePoissonSampleis$x,
                                 theta=c(truevalue[1:2], log(truevalue[3:4]),  truevalue[5], 
                                         log(truevalue[6])),pos=TwostagePoissonSampleis$ID_unit, sc=TwostagePoissonSampleis$PSU, fss=pi1is, 
                                 n2infor=n2is , N2)
}	


# install.packages("RColorBrewer")
library(RColorBrewer)
color<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c" )

#boxplot for uninformative sampling (NML, PL and WPL)
color=c( rep(color, 3))
name=c("alpha_NML", "beta_NML", "sigma^2_NML", "tau^2_NML", "alpha_PL", "beta_PL", "sigma^2_PL", "tau^2_PL", "alpha_WPL", "beta_WPL", "sigma^2_WPL", "tau^2_WPL" )
boxplot(cbind(Fit_NML[,c(1:4)],Fit_PL[,c(1:4)], Fit_WPL[,c(1:4)]) ,   col=color)
abline(h=0)

#boxplot for informative sampling (NML,PL and WPL)
boxplot(cbind(Fitis_NML[,c(1:4)],Fitis_PL[,c(1:4)], Fitis_WPL[,c(1:4)]) ,   col=color)
abline(h=0)


#define the squre root of J
sqrt_diagJ_PL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagJ_PL[i,]= sqrt(diag(J_PL[,,i]))
}

sqrt_diagJis_PL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagJis_PL[i,]= sqrt(diag(Jis_PL[,,i]))
}

sqrt_diagJ_WPL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagJ_WPL[i,]= sqrt(diag(J_WPL[,,i]))
}

sqrt_diagJis_WPL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagJis_WPL[i,]= sqrt(diag(Jis_WPL[,,i]))
}

#define the squre root of G
sqrt_diagG_PL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagG_PL[i,]= sqrt(diag(G_PL[,,i]))
}

sqrt_diagGis_PL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagGis_PL[i,]= sqrt(diag(Gis_PL[,,i]))
}

sqrt_diagG_WPL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagG_WPL[i,]= sqrt(diag(G_WPL[,,i]))
}

sqrt_diagGis_WPL=matrix(0, nrow=LOTS, ncol=6)
for ( i in 1:LOTS){
    sqrt_diagGis_PL[i,]= sqrt(diag(Gis_WPL[,,i]))
}

#create a table for latex
#install.packages("xtable")
library(xtable)
construct_header <- function(df, grp_names, span, align = "c", draw_line = T) {
   if (length(align) == 1) align <- rep(align, length(grp_names))
   if (!all.equal(length(grp_names), length(span), length(align)))
      stop("grp_names and span have to have the same length!")
   
   if (ncol(df) < sum(span)) stop("Span has to be less or equal to the number of columns of df") 
   
   header <- mapply(function(s, a, grp) sprintf("\\multicolumn{%i}{%s}{%s}", s, a, grp),
                    span, align, grp_names)
   header <- paste(header, collapse = " & ")
   header <- paste0(header, " \\\\")
   
   if (draw_line) {
      # where do we span the lines:
      min_vals <- c(1, 1 + cumsum(span)[1:(length(span) - 1)])
      max_vals <- cumsum(span)
      line <- ifelse(grp_names == "", "", 
                     sprintf("\\cmidrule(lr){%i-%i}", min_vals, max_vals))
      line <- paste(line[line != ""], collapse = " ")
      
      header <- paste0(header, "  ", line, "\n  ")
   }
   
   addtorow <- list(pos = list(-1, -1, nrow(df)),
                    command = c("\\hline\n  ", header, "\\hline\n  "))
   return(addtorow)
}





#bias and sd for uninformative sampling (NML, PL, WPL)
df<- matrix(c(apply(Fit_NML, 2,  mean), apply(Fit_NML, 2, sd), apply(Fit_PL, 2, mean), apply(Fit_PL, 2, sd),
apply(sqrt_diagG_PL, 2, mean),  apply(Fit_WPL, 2, mean), apply(Fit_WPL, 2, sd),
apply(sqrt_diagG_WPL, 2, mean)),ncol=8)

df[, c(1,3, 5, 6, 8)]=df[, c(1, 3, 5, 6, 8 )]*100
df<-round(df, 2)
df<-cbind(c("beta_0", "beta_1", "sigma^2", "tau_11^2", "tau_12", "tau_22^2"), df)
colnames(df)<-c("",c("bias*100", "sd"), rep(c("bias*100", "sd","esd*100"), 2))
df           
df_header <- construct_header(
   # the data.frame or matrix that should be plotted
   df,
   # the labels of the groups that we want to insert
   grp_names = c("uninformtive", "NML", "PL", "WPL"),
   # the number of columns each group spans
   span = c(1, 2, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)
print(xtable(df), add.to.row = df_header, include.rownames = F,  hline.after = F, latex.environments=NULL,booktabs=TRUE)




#variance estimator for uninformative sampling (PL, WPL)          
vardf<- matrix(c( apply(PS_PL, 2, mean),apply(PS_PL, 2, sd),  apply(sqrt_diagJ_PL, 2, mean) ,
apply(PS_WPL, 2, mean),apply(PS_WPL, 2, sd),  apply(sqrt_diagJ_WPL, 2, mean)),ncol=6)

vardf<-round(vardf,2)
vardf<-cbind(c("beta_0", "beta_1", "sigma^2", "tau_11^2", "tau_12", "tau_22^2"), vardf)
colnames(vardf)<-c("parameter", rep(c("mean of PS", "sd of PS", "esd of PS"), 2))
vardf
vardf_header <- construct_header(
   # the data.frame or matrix that should be plotted
   vardf,
   # the labels of the groups that we want to insert
   grp_names = c("",  "PL", "WPL"),
   # the number of columns each group spans
   span = c(1, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)
print(xtable(vardf), add.to.row = vardf_header, include.rownames = F, hline.after = F)


#bias and sd for informative sampling (NML, PL, WPL)
dfis<- matrix(c(apply(Fitis_NML, 2,  mean), apply(Fitis_NML, 2, sd), apply(Fitis_PL, 2, mean), apply(Fitis_PL, 2, sd),
apply(sqrt_diagGis_PL, 2, mean),  apply(Fitis_WPL, 2, mean), apply(Fitis_WPL, 2, sd),
apply(sqrt_diagGis_WPL, 2, mean)),ncol=8)

df[, c(1,3, 5, 6, 8)]=df[, c(1, 3, 5, 6, 8 )]*100
dfis=round(dfis, 2)


dfis<-cbind(c("beta_0", "beta_1", "sigma^2", "tau_11^2", "tau_12", "tau_22^2"), dfis)
colnames(dfis)<-c("parameter", c("bias*100", "sd"), rep(c("bias*100", "sd", "esd*100"), 2))
dfis

dfis_header <- construct_header(
   # the data.frame or matrix that should be plotted
   dfis,
   # the labels of the groups that we want to insert
   grp_names = c("", "NML", "PL", "WPL"),
   # the number of columns each group spans
   span = c(1, 2, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)
print(xtable(dfis), add.to.row = dfis_header, include.rownames = F,  hline.after = F, latex.environments=NULL,booktabs=TRUE)


#variance estimator for informative sampling (PL, WPL)          
vardfis<-matrix(c( apply(PSis_PL, 2, mean),apply(PSis_PL, 2, sd),  apply(sqrt_diagJis_PL, 2, mean) ,
apply(PSis_WPL, 2, mean),apply(PSis_WPL, 2, sd),  apply(sqrt_diagJis_WPL, 2, mean)),ncol=6)
vardfis<-round(vardfis,2)
vardfis<-cbind(c("beta_0", "beta_1", "sigma^2", "tau_11^2", "tau_12", "tau_22^2"), vardfis)
colnames(vardfis)<-c("parameter", rep(c("mean of PS", "sd of PS", "esd of PS"), 2))
vardfis
vardfis_header <- construct_header(
   # the data.frame or matrix that should be plotted
   vardfis,
   # the labels of the groups that we want to insert
   grp_names = c("",  "PL", "WPL"),
   # the number of columns each group spans
   span = c(1, 3, 3),
   # the alignment of each group, can be a single character (lcr) or a vector
   align = "c"
)   
print(xtable(vardfis), add.to.row = vardfis_header, include.rownames = F, hline.after = F, latex.environments=NULL,booktabs=TRUE)


