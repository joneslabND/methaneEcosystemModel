##### Simplest CH4 ecosystem model
##### 11-5-14
##### WEW & SEJ

rm(list=ls())

library(deSolve)

# set up ordinary differential equations (ODEs) in function for deSolve to use for each time step

# State variables:
#		1. CH4 in profundal sediment
#		2. CH4 in littoral sediment
#		3. CH4 in epilimnion

# Daily time step

timeStepCH4<-function(times,y,params){
	with(as.list(params),{
		P=y[1]	# CH4 in profundal sediment
		L=y[2]	# CH4 in littoral sediment
		E=y[3]	# CH4 in epilimnion
	
		#dPdt=production-diffusion; assuming no ebullition
		#		production=f(GPP,settling,yield)
		#		diffusion=f(???)
		dP.dt=Ap*(GPP*zmix*sed*yield-Pdiff*(P/Vp-E/Ve))
		#dLdt=production-diffusion-ebullition
		#		production=f(GPP,settling,yield)
		#		diffusion=f(E,L,???)
		#		ebullition=f(L,???)
		dL.dt=Al*(GPP*zmix*sed*yield-Ldiff*(L/Vl-E/Ve))-ebull*L
		#dEdt=production+diffusion (littoral and profundal)-atm. diffusion-oxidation
		#		production=f(GPP)????
		#		diffusion=f(E,L,P,???)
		#		atm. diffusion=f(E,atm)
		#		oxidation=f(E)
		dE.dt=GPP*Eprod*Ve+Pdiff*(P/Vp-E/Ve)*Ap+Ldiff*(L/Vl-E/Ve)*Al-k*(E/Ve-atmCH4)*Ae-oxE*E
	
		list(c(dP.dt,dL.dt,dE.dt))
	})
}

# Parameters
Ap=1250		# Area of profundal sediments; [m2]
Al=250		# Area of littoral sediments; [m2]
Ae=Ap+Al		# Surface area of lake; [m2]
activeLayerDepth=0.2	# Active layer depth of sediments; [m]
Vp=Ap*activeLayerDepth		# Volume of "active" profundal sediments; [m3]
Vl=Al*activeLayerDepth		# Volume of "active" littoral sediments; [m3]
zmix=2	# Mixed layer depth of lake; [m]
Ve=Ae*zmix		# Volume of epilimnion; [m3]
sed=0.3		# Proportion of areal GPP that settles to sediments; []
yield=0.25	# Molar yield of CH4 from algal C; [mol CH4 (mol algal C)-1]
Pdiff=0.001	# Mass transfer coefficient/diffusivity???; [m d-1]
Ldiff=0.005	# Mass transfer coefficient/diffusivity???; [m d-1]
ebull=0.1		# Fraction of production released as ebullition [d-1] -> could be a more complex function that causes ebullition at a critical saturating concentration; could also be probabalistic
Eprod=0.0001	# Epilimnetic CH4 production per unit GPP; [mol CH4 (mol C)-1]
k=0.4	# Piston velocity; [m d-1]
oxE=0.9	# Fraction of epilimnion CH4 lost to oxidation; [d-1]
params=c(Vp=Vp,Vl=Vl,Ve=Ve,Ap=Ap,Al=Al,Ae=Ae,zmix=zmix,sed=sed,yield=yield,Pdiff=Pdiff,Ldiff=Ldiff,ebull=ebull,Eprod=Eprod,k,oxE)

init=c(P=0.2*Vp,L=0.2*Vl,E=0.2*Ve)		#starting all the same (sort of like after mixis); 200 uM
t.s=1:180

GPP=0.1		# [mol C m-3 d-1]Jake says this is average for morris --> Stuart will remember how to do this dynamically
atmCH4=0.1	# [mol CH4 m-3]Will thinks this is about right --> We can check this and use whatever the global average is

out=ode(y=init,times=t.s,func=timeStepCH4,parms=params)

dev.new()
par(mfrow=c(3,1))
plot(out[,1],out[,2],type='l',lwd=2,xlab="Time",ylab="Profundal CH4")
plot(out[,1],out[,3],type='l',lwd=2,xlab="Time",ylab="Littoral CH4",col='red')
plot(out[,1],out[,4],type='l',lwd=2,xlab="Time",ylab="Epilimnion CH4",col='green')