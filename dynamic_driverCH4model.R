##### Simplest CH4 ecosystem model
##### 11-5-14
##### WEW & SEJ

rm(list=ls())

library(deSolve)

#### set up ordinary differential equations (ODEs) in function for deSolve to use for each time step
	# State variables:
	#	1. CH4 in profundal sediment
	#	2. CH4 in littoral sediment
	#	3. CH4 in epilimnion

timeStepCH4<-function(time,y,params){
	with(as.list(params),{
		# state variables
		P=y[1]	# CH4 in profundal sediment
		Psom=y[2]	# sediment organic matter in profundal sediment
		L=y[3]	# CH4 in littoral sediment
		Lsom=y[4]	# sediment organic matter in littoral sediment
		E=y[5]	# CH4 in epilimnion
	
		# forcing data
		GPP=varied_GPP(time)
		atmCH4=varied_atmCH4(time)
		k=varied_k(time)
		
		#dP.dt=production-diffusion; assuming no ebullition
		#		production=f(Psom,yield)
		#		diffusion=f(???)
		dP.dt=maxPsomUptake*Psom/(kPprod*Vp+Psom)*Psom*yield-Pdiff*(P/Vp-E/Ve)*Ap
		#dPsom.dt=settling material-conversion to methane-permanent burial
		#		settling=f(GPP,settling)
		#		conversion=f(Psom)
		#		burial=f(Psom)
		dPsom.dt=Ap*GPP*zmix*sed-maxPsomUptake*Psom/(kPprod*Vp+Psom)*Psom-Psom*permBuryP
		#dL.dt=production-diffusion-ebullition
		#		production=f(Lsom,yield)
		#		diffusion=f(E,L,???)
		#		ebullition=f(L,???)
		dL.dt=maxLsomUptake*Lsom/(kLprod*Vl+Lsom)*Lsom*yield-Ldiff*(L/Vl-E/Ve)*Al-ebull*L
		#dLsom.dt=settling material - conversion to methane - permanent burial
		#		settling=f(GPP,settling)
		#		conversion=f(Lsom)
		#		burial=f(Lsom)
		dLsom.dt=Al*GPP*zmix*sed-maxLsomUptake*Lsom/(kLprod*Vl+Lsom)*Lsom-Lsom*permBuryL
		#dEdt=production+diffusion (littoral and profundal)-atm. diffusion-oxidation
		#		production=f(GPP)????
		#		diffusion=f(E,L,P,???)
		#		atm. diffusion=f(E,atm)
		#		oxidation=f(E)
		dE.dt=GPP*Eprod*Ve+Pdiff*(P/Vp-E/Ve)*Ap+Ldiff*(L/Vl-E/Ve)*Al-k*(E/Ve-atmCH4)*Ae-oxE*E
	
		list(c(dP.dt,dPsom.dt,dL.dt,dLsom.dt,dE.dt))
	})
}


#### Parameters
Ap=1250		# Area of profundal sediments; [m2]
Al=250		# Area of littoral sediments; [m2]
Ae=Ap+Al		# Surface area of lake; [m2]
activeLayerDepth=0.2	# Active layer depth of sediments; [m]
Vp=Ap*activeLayerDepth		# Volume of "active" profundal sediments; [m3]
Vl=Al*activeLayerDepth		# Volume of "active" littoral sediments; [m3]
zmix=2	# Mixed layer depth of lake; [m]
Ve=Ae*zmix		# Volume of epilimnion; [m3]
sed=0.3		# Proportion of areal GPP that settles to sediments; []
maxPsomUptake=0.004	# Maximum rate of Psom uptake; [d-1]
kPprod=10		# half-saturation constant for methane production in profundal; [mol algal C m-3]
maxLsomUptake=0.004	# Maximum rate of Psom uptake; [d-1]
kLprod=10		# half-saturation constant for methane production in profundal; [mol algal C m-3]
permBuryP=0.1	# permanent burial of Psom; [d-1]
permBuryL=0.1	# permanent burial of Lsom; [d-1]
yield=0.25		# Molar yield of CH4 from algal C; [mol CH4 (mol algal C)-1]
Pdiff=0.001	# Mass transfer coefficient/diffusivity???; [m d-1]
Ldiff=0.005	# Mass transfer coefficient/diffusivity???; [m d-1]
ebull=0.1		# Fraction of production released as ebullition [d-1] -> could be a more complex function that causes ebullition at a critical saturating concentration; could also be probabalistic
Eprod=0.0001	# Epilimnetic CH4 production per unit GPP; [mol CH4 (mol C)-1]
oxE=0.9	# Fraction of epilimnion CH4 lost to oxidation; [d-1]
params=c(Vp=Vp,Vl=Vl,Ve=Ve,Ap=Ap,Al=Al,Ae=Ae,zmix=zmix,sed=sed,maxPsomUptake=maxPsomUptake,kPprod=kPprod,maxLsomUptake=maxLsomUptake,kLprod=kLprod,permBuryP=permBuryP,permBuryL=permBuryL,yield=yield,Pdiff=Pdiff,Ldiff=Ldiff,ebull=ebull,Eprod=Eprod,oxE)


#### forcing data (GPP, atmCH4, k) and run parameters
nTimeSteps=180
t.s=1:nTimeSteps

# Could load observations, model output, or set to a constant
constant_GPP=0.1	# [mol C m-3 d-1]Jake says this is average for morris --> Stuart will remember how to do this dynamically
constant_atmCH4=0.1	# [mol CH4 m-3]Will thinks this is about right --> We can check this and use whatever the global average is
constant_k=0.4	# [m d-1]using a decent number for UNDERC lakes

obs_GPP=NULL	# can be set to two column matrix with simulation timestep in first column and observed/modeled GPP [mol C m-3 d-1] in second column
obs_atmCH4=NULL	# can be set to two column matrix with simulation timestep in first column and observed/modeled atmCH4 [mol CH4 m-3] in second column
obs_k=NULL	# can be set to two column matrix with simulation timestep in first column and observed/modeled k [m d-1] in second column

# generate approxfun functions for use in model integration
if(is.null(obs_GPP)){
	varied_GPP=approxfun(x=1:nTimeSteps,y=rep(constant_GPP,nTimeSteps),method="linear",rule=2)
}else{
	varied_GPP=approxfun(x=obsGPP[,1],y=obs_GPP[,2],method="linear",rule=2)
}
if(is.null(obs_atmCH4)){
	varied_atmCH4=approxfun(x=1:nTimeSteps,y=rep(constant_atmCH4,nTimeSteps),method="linear",rule=2)
}else{
	varied_atmCH4=approxfun(x=obs_atmCH4[,1],y=obs_atmCH4[,2],method="linear",rule=2)
}
if(is.null(obs_k)){
	varied_k=approxfun(x=1:nTimeSteps,y=rep(constant_k,nTimeSteps),method="linear",rule=2)
}else{
	varied_k=approxfun(x=obs_k[,1],y=obs_k[,2],method="linear",rule=2)
}


# initial values
init=c(P=0.2*Vp,Psom=0,L=0.2*Vl,Lsom=0,E=0.2*Ve)		# [mol CH4] starting all the same (sort of like after mixis); 200 uM
	
#### integrate model
out=ode(y=init,times=t.s,func=timeStepCH4,parms=params)

#### plot output
dev.new()
par(mfrow=c(5,1))
plot(out[,1],out[,2],type='l',lwd=2,xlab="Time",ylab="Profundal CH4 [mol CH4]",col='black')
plot(out[,1],out[,3],type='l',lwd=2,xlab="Time",ylab="Profundal OM [mol algal C]",col='black',lty=2)
plot(out[,1],out[,4],type='l',lwd=2,xlab="Time",ylab="Littoral CH4 [mol CH4]",col='red')
plot(out[,1],out[,5],type='l',lwd=2,xlab="Time",ylab="Littoral OM [mol algal C]",col='red',lty=2)
plot(out[,1],out[,6],type='l',lwd=2,xlab="Time",ylab="Epilimnion CH4",col='green')

print(round(out[nrow(out),],4))