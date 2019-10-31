if (FALSE){
production=function(biomass,catch,r=.3,k=1000,p=1,cv=.1) 
   max(biomass*r/p*(1-(biomass/k)^p),0)*rlnorm(1,0,cv)
  

err=0.2
catch     =rep( 30,100)
biomass   =rep(1000,100)
catch     =rep(  80,100)
index     =rep(1000,100)
Error     =numeric(100)
integral  =numeric(100)
derivative=numeric(100)

for (i in 2:50){
  catch[i]  =80
  biomass[i]=biomass[i-1]-catch[i]+production(biomass[i-1])
  index[i]  =biomass[i]*rlnorm(1,0,err)
  }

setPt=mean(index[1:10])/2
plot(biomass[1:50],type="line",ylim=c(0,1500),xlab="Year")
points(index[1:50])
abline(h=setPt,col="red")

Error[1:50] = index[1:50] - setPt

k1=0.02
k2=0.05
k3=0.02
biomass[50]=biomass[50-1]-catch[50]+production(biomass[50-1])

for (i in 51:100){
  Error[i] = index[i-1] - setPt 
  integral[i] = sum(Error[49:(i-1)])
  derivative[i] = (Error[i-1]-Error[i-2])
  
  catch[i]  = catch[i-1]-k1*Error[i]+k2*integral[i]+k3*derivative[i]
  
  catch[i]=min(max(catch[i],0),100)
  
  biomass[i]=biomass[i-1]-catch[i]+production(biomass[i-1])
  index[i]  =biomass[i]*rlnorm(1,0,err)
  }


plot(biomass,type="line",ylim=c(0,1500),xlab="Year")
plot(catch,  type="line",xlab="Year")
points(index)
abline(h=setPt,col="red")
}

