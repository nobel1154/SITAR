require(splines)
require(MCMCpack)
require(multinomRob)

data<-read.csv("C:/Users/bhuiy/Downloads/FW%3a_update_with_the_SITAR_model_project/sitar.csv",header=TRUE)
head(data)
dim(data)

sigmaalpha=3
alphascale=1
sigmabeta=3
betascale=1
sigmagamma=3
gammascale=1

lambda=0.1
n=500
T=10
y<-matrix(0,nrow=n, ncol=T)
for (i in 1:500)
{
  for (j in 1:10)
  {
    y[i,j]<-data[(i-1)*10+j,14]
  }
  
}

time<-matrix(0,nrow=n, ncol=T)
for (i in 1:500)
{
  for (j in 1:10)
  {
    time[i,j]<-data[(i-1)*10+j,8]
  }
  
}

N=5000
k=2

p<-matrix(0,nrow=N,ncol=k)
pz<-matrix(0,nrow=n,ncol=k)
z<-matrix(0,nrow=N,ncol=n)  #nis number of person
lamda<-matrix(1,nrow=N,ncol=k)
alpha<-matrix(0,nrow=N,ncol=n)
beta<-matrix(0,nrow=N,ncol=n)
gamma<-matrix(0,nrow=N,ncol=n)
palpha<-matrix(0,nrow=N,ncol=k)
pbeta<-matrix(0,nrow=N,ncol=k)
pgamma<-matrix(0,nrow=N,ncol=k)
malpha<-matrix(1,nrow=N,ncol=k)
mbeta<-matrix(1,nrow=N,ncol=k)
mgamma<-matrix(1,nrow=N,ncol=k)
valpha<-matrix(1,nrow=N,ncol=k)
vbeta<-matrix(1,nrow=N,ncol=k)
vgamma<-matrix(1,nrow=N,ncol=k)


resi<-matrix(0,nrow=N,ncol=1)
meanalpha<-10
meanbeta<-0
meanalpha<-0
meangamma<-0
variancetheta<-100
resishape<-5
resiscale<-1


for (i in 2:N)
  
  
{
  
  
  s<-y[1,]-alpha[i-1,1]
  ti<-(time[1,]-beta[i-1,1])/exp(-gamma[i-1,1])
  
  
  
  h<-predict(lm(s ~ ns(ti, df = 5)))
  
  ss=sum((s-h)^2)
  
  for (l in 2:n)
  {
    
    s<-y[l,]-alpha[i-1,l]
    ti<-(time[l,]-beta[i-1,l])/exp(-gamma[i-1,l])
    
    
    h<-predict(lm(s ~ ns(ti, df = 5)))
    
    ss=ss+sum((s-h)^2)
    
  }
  
  
  resi[i]<-rinvgamma(1, shape=resishape+n*T/2, scale =resiscale+ss/2)
  
  
  for (j in 1:k)
  {
    lamda[i,j]=lamda[i-1,j]+sum(z[i-1,]==j)
  }
  
  p[i,]<-rdirichlet(1, lamda[i,])
  
  
  
  for (l in 1:n)
  {
    
    sum=dnorm(alpha[i-1,l],malpha[i-1,1],sqrt(valpha[i-1,1]))*dnorm(beta[i-1,l],mbeta[i-1,1],sqrt(vbeta[i-1,1]))*dnorm(gamma[i-1,l],mgamma[i-1,1],sqrt(vgamma[i-1,1]))*p[i,1]
    for (u in 2:k)
    {
      sum=sum+dnorm(alpha[i-1,l],malpha[i-1,u],sqrt(valpha[i-1,u]))*dnorm(beta[i-1,l],mbeta[i-1,u],sqrt(vbeta[i-1,u]))*dnorm(gamma[i-1,l],mgamma[i-1,u],sqrt(vgamma[i-1,u]))*p[i,u]
    }
    
    for (t in 1:k)
    {
      
      
      pz[l,t]<-dnorm(alpha[i-1,l],malpha[i-1,t],sqrt(valpha[i-1,t]))*dnorm(beta[i-1,l],mbeta[i-1,t],sqrt(vbeta[i-1,t]))*
        dnorm(gamma[i-1,l],mgamma[i-1,t],sqrt(vgamma[i-1,t]))*p[i,t]/sum
    }
    
    M<-rmultinom(n = 1,size=1,  c(pz[l,]))
    
    for (w in 1:k)
    {
      
      if (M[w]==1) { z[i,l]=w}
    }
    
    
    
  }
  
  
  
  
  
  
  
  
  for (l in 1:n)
  {
    s<-y[l,]-alpha[i-1,l]
    ti<-(time[l,]-beta[i-1,l])/exp(-gamma[i-1,l])
    
    
    h<-predict(lm(s ~ ns(ti, df = 5)))
    
    
    mean_alpha<-(valpha[i-1,z[i,l]]*sum(y[l,]-h)+malpha[i-1,z[i,l]]*resi[i])/(resi[i]+n*valpha[i-1,z[i,l]])
    var_alpha<-resi[i]*valpha[i-1,z[i,l]]/(resi[i]+n*valpha[i-1,z[i,l]])
    alpha[i,l]<-rnorm(1,mean_alpha,sqrt(var_alpha))
    
    
  }
  
  
  
  
  
  for (l in 1:n)
  {
    
    beta0<-beta[i-1,l]
    beta1<-rnorm(1,beta0,lambda*sqrt(2))
    s<-y[l,]-alpha[i,l]
    
    
    ti<-(time[l,]-beta[i-1,l])/exp(-gamma[i-1,l])
    
    
    h<-predict(lm(s ~ ns(ti, df = 5)))
    den.tor<-dnorm(y[l,1],alpha[i,l]+h[1],sqrt(resi[i]))
    
    for (ks in 2:T)
    {
      
      den.tor<-den.tor*dnorm(y[l,ks],alpha[i,l]+h[ks],sqrt(resi[i]))
      
    }
    
    den.tor<-den.tor*dnorm(beta0,mbeta[i-1,z[i,l]],sqrt(vbeta[i-1,z[i,l]]))
    ti<-(time[l,]-beta1)/exp(-gamma[i-1,l])
    
    
    h<-predict(lm(s ~ ns(ti, df = 5)))
    
    num.tor<-dnorm(y[l,1],alpha[i,l]+h[1],sqrt(resi[i]))
    
    for (ks in 2:T)
    {
      
      num.tor<-num.tor*dnorm(y[l,ks],alpha[i,l]+h[ks],sqrt(resi[i]))
      
    }
    
    num.tor<-num.tor*dnorm(beta1,mbeta[i-1,z[i,l]],sqrt(vbeta[i-1,z[i,l]]))
    
    u <- ( runif(1) < min(1,num.tor/den.tor))
    if (u ==1)
    {
      beta[i,l] = beta1
    } else {
      beta[i,l] = beta0
    }
    
    
  }
  
  
  
  
  for (l in 1:n)
  {
    
    gamma0<-gamma[i-1,l]
    gamma1<-rnorm(1,gamma0,lambda*sqrt(2))
    s<-y[l,]-alpha[i,l]
    
    
    ti<-(time[l,]-beta[i,l])/exp(-gamma[i-1,l])
    
    
    h<-predict(lm(s ~ ns(ti, df = 5)))
    den.tor<-dnorm(y[l,1],alpha[i,l]+h[1],sqrt(resi[i]))
    
    for (ks in 2:T)
    {
      
      den.tor<-den.tor*dnorm(y[l,ks],alpha[i,l]+h[ks],sqrt(resi[i]))
      
    }
    
    den.tor<-den.tor*dnorm(gamma0,mgamma[i-1,z[i-1,l]],sqrt(vgamma[i-1,z[i-1,l]]))
    ti<-(time[l,]-beta[i,l])/exp(-gamma1)
    
    
    h<-predict(lm(s ~ ns(ti, df = 5)))
    
    num.tor<-dnorm(y[l,1],alpha[i,l]+h[1],sqrt(resi[i]))
    
    for (ks in 2:T)
    {
      
      num.tor<-num.tor*dnorm(y[l,ks],alpha[i,l]+h[ks],sqrt(resi[i]))
      
    }
    
    num.tor<-num.tor*dnorm(gamma1,mgamma[i-1,z[i-1,l]],sqrt(vgamma[i-1,z[i-1,l]]))
    
    u <- ( runif(1) < min(1,num.tor/den.tor))
    if (u ==1)
    {
      gamma[i,l] = gamma1
    } else {
      gamma[i,l] = gamma0
    }
    
    
  }
  
  
  
  
  for (m in 1:k)
  {
    valpha[i,m]<-rinvgamma(1,sigmaalpha+sum(z[i,]==m)/2,scale=alphascale+sum((alpha[i,z[i,]==m]-malpha[z[i,]==m])^2)/2)
    
  }
  
  
  for (m in 1:k)
  {
    vbeta[i,m]<-rinvgamma(1,sigmabeta+sum(z[i,]==m)/2,scale=betascale+sum((beta[i,z[i,]==m]-mbeta[z[i,]==m])^2)/2)
    
  }
  
  
  
  
  for (m in 1:k)
  {
    vgamma[i,m]<-rinvgamma(1,sigmagamma+sum(z[i,]==m)/2,scale=gammascale+sum((gamma[i,z[i,]==m]-mgamma[z[i,]==m])^2)/2)
    
  }
  
  
  
  
  
  for (m in 1:k)
  {
    mean1=mean(alpha[i,z[i,]==m])*variancetheta/(variancetheta+valpha[i,m]/sum(z[i,]==m))+(meanalpha*valpha[i,m]/sum(z[i,]==m))/(variancetheta+valpha[i,m]/sum(z[i,]==m))
    variance1<-1/(1/variancetheta+ sum(z[i,]==m)/valpha[i,m])
    malpha[i,m]=rnorm(1,mean=mean1,sd=sqrt(variance1))
    
  }
  
  
  
  for (m in 1:k)
  {
    mean1=mean(beta[i,z[i,]==m])*variancetheta/(variancetheta+vbeta[i,m]/sum(z[i,]==m))+meanbeta*vbeta[i,m]/(variancetheta+vbeta[i,m]/sum(z[i,]==m))
    variance1<-1/(1/variancetheta+ sum(z[i,]==m)/vbeta[i,m])
    mbeta[i,m]=rnorm(1,mean=mean1,sd=sqrt(variance1))
    
  }
  
  
  
  
  for (m in 1:k)
  {
    mean1=mean(gamma[i,z[i,]==m])*variancetheta/(variancetheta+vgamma[i,m]/sum(z[i,]==m))+meangamma*vgamma[i,m]/(variancetheta+vgamma[i,m]/sum(z[i,]==m))
    variance1<-1/(1/variancetheta+ sum(z[i,]==m)/vgamma[i,m])
    mgamma[i,m]=rnorm(1,mean=mean1,sd=sqrt(variance1))
    
  }
  
  print (i)
  
}
