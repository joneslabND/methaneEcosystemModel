# methane diffusivity
	#observations at 1 atm pressure
	#from 1973M1 in Data extract from Landolt-Bornstein - Group IV Physical Chemistry Volume 15A, Chapter 5.1.3: Diffusion of Gas/Vapor in Liquid; 2007, pp 1880-1882. J. WInkelmann
	
	#	Diffusion of methane (1); water (2)
tempK=c(273,278,283,288,293,298,303,308)	# K
D=c(1.05,1.1,1.22,1.33,1.47,1.64,1.74,1.96)	# *10^9 m2 s-1

tempK2=tempK*tempK
summary(lm(D~tempK+tempK2))
B0=2.053e1
B1=-1.579e-1
B2=3.167e-4

plot(tempK,D,xlab="temperature [K]",ylab="diffusivity [10^9 m2 s-1]")
pred=273:310
pred2=pred*pred
lines(pred,B0+B1*pred+B2*pred2,col='red')


#Diffusion+Production

Zmax=10
Zmix=2

temperature=8	#degreesC
K=273+temperature

prod=10

nsteps=200
c=matrix(0,2,nsteps+1)

D=(B0+B1*K+B2*K*K)/1e9

c[1,1]=100
c[2,1]=0

for(i in 1:nsteps){
	J=D*(c[1,i]-c[2,i])/(Zmax-Zmix)
	c[1,(i+1)]=c[1,i]-J+prod
	c[2,(i+1)]=c[2,i]+J
}

plot(1:(nsteps+1),c[1,],type='l',ylim=c(0,200))
lines(1:(nsteps+1),c[2,],col='red')



#Diffusion+Production w/ more than 2 depths
source("RKfix.R")

depths=Zmax-Zmix
dz=1	#m

prod=1

nsteps=1000
c=matrix(0,depths,nsteps+1)

D=(B0+B1*K+B2*K*K)/1e9*60*60*24	#m2 d-1		; bastviken 2002 used 60*60*24/1e8 to 60*60*24/1e5 m2 d-1

c[1,1]=1
c[2:depths,1]=0

dp<-function(x){
	dxdt=numeric(length(x))
	dxdt[1]=-D*(x[1]-x[2])/dz+prod
	for(i in 2:(length(dxdt)-1)){
		dxdt[i]=D*(x[(i-1)]-x[i])-D*(x[i]-x[(i+1)])
	}
	dxdt[length(dxdt)]=D*(x[(length(dxdt)-1)]-x[length(dxdt)])
	
	return(dxdt)
}


for(i in 1:nsteps){
	c[,(i+1)]=Ynext(dp,0,1,c[,i],5)
}

plot(1:(nsteps+1),c[1,],type='l',ylim=c(0,400))
lines(1:(nsteps+1),c[2,],col='red')



# from Rimmer et al. 2006
C0=rep(0,20)	
deltaZ=0.4	# m

DT=0.2		#m2 d-1
A=100	#mg m-2 d-1

deltaT=(0.5*deltaZ^2)/DT		# day

nsteps=400

C=matrix(NA,length(C0),(nsteps+1))
C[,1]=C0

for(i in 1:nsteps){
	C[1,(i+1)]=C[1,i]+(deltaT/(deltaZ^2))*((A*deltaZ)+DT*(C[2,i]-C[1,i]))
	for(j in 2:(length(C0)-1)){ 
		C[j,(i+1)]=C[j,i]+DT*deltaT/(deltaZ^2)*(C[j+1,i]+C[j-1,i]-2*C[j,i])
	}
	#static thermocline
	C[length(C0),(i+1)]=C[length(C0),i]+DT*(deltaT/(deltaZ^2))*(C[(length(C0)-1),i]-C[length(C0),i])
}

quartz()
plot(C[,1],seq(0,(length(C0)-1)*deltaZ,by=deltaZ),type='o',xlim=range(C))
lines(C[,ncol(C)],seq(0,(length(C0)-1)*deltaZ,by=deltaZ),col='red')
points(C[,ncol(C)],seq(0,(length(C0)-1)*deltaZ,by=deltaZ),col='red')
