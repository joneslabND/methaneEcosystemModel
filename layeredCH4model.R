##### Layered water column CH4 ecosystem model
##### 04-20-17
##### SEJ

rm(list=ls())

library(deSolve)

#### set up ordinary differential equations (ODEs) in function for deSolve to use for each time step
	# State variables:
	#	1 to (n-1): vector of CH4 in hypolimnion layers
	#	n: CH4 in epilimnion

timeStepCH4<-function(time,y,params,As,Vs){
	with(as.list(params),{
		# state variables
		H=y[1:N]	# CH4 in hypolimnion
		P=y[(N+1)]
		Psom=y[(N+2)]	# sediment organic matter in profundal sediment
		L=y[(N+3)]	# CH4 in littoral sediment
		Lsom=y[(N+4)]	# sediment organic matter in littoral sediment
		E=y[(N+5)]	# CH4 in epilimnion
	
		dH.dt=numeric(N)
	
		# forcing data
		GPP=varied_GPP(time)
		atmCH4=varied_atmCH4(time)
		k=varied_k(time)
		
		#dP.dt=production-diffusion; assuming no ebullition
		#		production=f(Psom,yield)
		#		diffusion=f(H,Pdiff)
		dP.dt=maxPsomUptake*Psom/(kPprod*(As[1]*activeLayerDepth)+Psom)*Psom*yield-(1/deltaZ)*Pdiff*(P/(As[1]*activeLayerDepth)-H[1]/Vs[1])*As[1]
		#dPsom.dt=settling material-conversion to methane-permanent burial
		#		settling=f(GPP,settling)
		#		conversion=f(Psom)
		#		burial=f(Psom)
		dPsom.dt=As[1]*GPP*zmix*sed-maxPsomUptake*Psom/(kPprod*(As[1]*activeLayerDepth)+Psom)*Psom-Psom*permBuryP
		#dH.dt=diffusion; assuming no ebullition
		#		diffusion=f(H,Dt)
		dH.dt[1]=(1/deltaZ)*(Pdiff*(P/(As[1]*activeLayerDepth)-H[1]/Vs[1])*As[1]+Hdiff*(H[2]/Vs[2]-H[1]/Vs[1])*As[2])
		for(i in 2:(N-1)){
			dH.dt[i]=(1/deltaZ)*Hdiff*((H[(i+1)]/Vs[(i+1)]-H[i]/Vs[i])*As[(i+1)]+(H[(i-1)]/Vs[(i-1)]-H[i]/Vs[i])*As[i])
		}
		dH.dt[N]=(1/deltaZ)*Hdiff*((E/Vs[(N+1)]-H[N]/Vs[N])*As[(N+1)]+(H[(N-1)]/Vs[N-1]-H[i]/Vs[i])*As[N])
		#dL.dt=production-diffusion-ebullition
		#		production=f(Lsom,yield)
		#		diffusion=f(E,L,???)
		#		ebullition=f(L,???)
		dL.dt=maxLsomUptake*Lsom/(kLprod*(As[(N+2)]*activeLayerDepth)+Lsom)*Lsom*yield-(1/deltaZ)*Ldiff*(L/(As[(N+2)]*activeLayerDepth)-E/Vs[(N+1)])*As[(N+2)]-ebull*L
		#dLsom.dt=settling material - conversion to methane - permanent burial
		#		settling=f(GPP,settling)
		#		conversion=f(Lsom)
		#		burial=f(Lsom)
		dLsom.dt=As[(N+2)]*GPP*zmix*sed-maxLsomUptake*Lsom/(kLprod*(As[(N+2)]*activeLayerDepth)+Lsom)*Lsom-Lsom*permBuryL
		#dEdt=production+diffusion (littoral and profundal)-atm. diffusion-oxidation
		#		production=f(GPP)????
		#		diffusion=f(E,L,P,???)
		#		atm. diffusion=f(E,atm)
		#		oxidation=f(E)
		dE.dt=GPP*Eprod*Vs[(N+1)]+(1/deltaZ)*Hdiff*(E/Vs[(N+1)]-H[N]/Vs[N])*As[(N+1)]+(1/deltaZ)*Ldiff*(L/(As[(N+2)]*activeLayerDepth)-E/Vs[(N+1)])*As[(N+2)]-k*(E/Vs[(N+1)]-atmCH4)*As[(N+3)]-oxE*E
		
		list(c(dH.dt,dP.dt,dPsom.dt,dL.dt,dLsom.dt,dE.dt))
	})
}


#### Parameters
zmax=6
zmix=2	# Mixed layer depth of lake; [m]
deltaZ=0.25					# hypolimnion slice thickness; [m]
N=(zmax-zmix)/deltaZ			# number of hypolimnion slices;
activeLayerDepth=0.2	# Active layer depth of sediments; [m]
Pdiff=0.05		# diffusivity; [m2 d-1]
Hdiff=0.15		# diffusivity; [m2 d-1]
sed=0.3		# Proportion of areal GPP that settles to sediments; []
maxPsomUptake=0.04	# Maximum rate of Psom uptake; [d-1]
kPprod=10		# half-saturation constant for methane production in profundal; [mol algal C m-3]
maxLsomUptake=0.04	# Maximum rate of Psom uptake; [d-1]
kLprod=10		# half-saturation constant for methane production in profundal; [mol algal C m-3]
permBuryP=0.1	# permanent burial of Psom; [d-1]
permBuryL=0.1	# permanent burial of Lsom; [d-1]
yield=0.25		# Molar yield of CH4 from algal C; [mol CH4 (mol algal C)-1]
Ldiff=0.05		# diffusivity; [m2 d-1]
ebull=0.1		# Fraction of production released as ebullition [d-1] -> could be a more complex function that causes ebullition at a critical saturating concentration; could also be probabalistic
Eprod=0.0001	# Epilimnetic CH4 production per unit GPP; [mol CH4 (mol C)-1]
oxE=0.99		# Fraction of epilimnion CH4 lost to oxidation; [d-1]

# for now an upside-down tophat shaped lake with radius of pelagic=900 m, and littoral zone 100 m wide
As=c(rep(pi*900^2,N+1),(pi*1000^2-pi*900^2),pi*1000^2)		# interface areas; [m2]
Vs=c(rep(pi*900^2*deltaZ,N),pi*1000^2*zmix)		# layer volumes; [m3]

params=c(activeLayerDepth=activeLayerDepth,deltaZ=deltaZ,N=16,Pdiff=Pdiff,Hdiff=Hdiff,zmix=zmix,sed=sed,maxPsomUptake=maxPsomUptake,kPprod=kPprod,maxLsomUptake=maxLsomUptake,kLprod=kLprod,permBuryP=permBuryP,permBuryL=permBuryL,yield=yield,Ldiff=Ldiff,ebull=ebull,Eprod=Eprod,oxE)



#### forcing data (GPP, atmCH4, k) and run parameters
nTimeSteps=360
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
init=rep(0,(N+5))		# [mol CH4] starting all the same (sort of like after mixis); 200 uM
	
#### integrate model
out=ode(y=init,times=t.s,func=timeStepCH4,parms=params,As=As,Vs=Vs)

#### plot output
profundalSedV=As[1]*activeLayerDepth
littoralSedV=As[(N+2)]*activeLayerDepth

# generate concentration data in uM CH4
outConc=out
outConc[,2:(N+1)]=out[,2:(N+1)]/Vs[1:N]*1000
outConc[,(N+2)]=out[,(N+2)]/profundalSedV*1000
outConc[,(N+3)]=out[,(N+3)]/profundalSedV*1000
outConc[,(N+4)]=out[,(N+4)]/littoralSedV*1000
outConc[,(N+5)]=out[,(N+5)]/littoralSedV*1000
outConc[,(N+6)]=out[,(N+6)]/Vs[(N+1)]*1000

dev.new()
par(mfrow=c(2,4))
plot(outConc[,1],outConc[,(N+2)],type='l',lwd=2,xlab="Time",ylab="Profundal sed. CH4 [uM CH4]",col='black')
plot(outConc[,1],outConc[,(N+3)],type='l',lwd=2,xlab="Time",ylab="Profundal OM [uM algal C]",col='black')
plot(outConc[,1],outConc[,2],type='l',lwd=2,xlab="Time",ylab="Hypo. CH4 [uM CH4]",col='black',ylim=range(outConc[,2:(N+1)]))
lines(outConc[,1],outConc[,(N+1)],lty=2,lwd=2)
plot(1,1,cex=0,axes=FALSE,xlab="",ylab="")
plot(outConc[,1],outConc[,(N+4)],type='l',lwd=2,xlab="Time",ylab="Littoral sed. CH4 [uM CH4]",col='black')
plot(outConc[,1],outConc[,(N+5)],type='l',lwd=2,xlab="Time",ylab="Littoral OM [uM algal C]",col='black',lty=2)
plot(outConc[,1],outConc[,(N+6)],type='l',lwd=2,xlab="Time",ylab="Epilimnion CH4 [uM CH4]",col='black')
plot(outConc[nTimeSteps,c(2:(N+1),(N+6))],c(seq((zmax-deltaZ),zmix,-deltaZ),zmix/2),type='o',ylim=c(zmax,0),ylab="Depth (m)",xlab="uM CH4")
