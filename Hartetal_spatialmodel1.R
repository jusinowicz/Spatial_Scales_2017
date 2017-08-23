################################################################################
#
# This is a spatially explicit simulation of N species competing in a 
# heterogeneous environment, numerical analysis of coexistence, and a simulated 
# sampling scheme to reproduce coexistence-area and species-area curves. Phase 1
# produces an equilibrium community, Phase 2 tests the coexistence of each 
# species in this community more formally by allowing it to invade against the 
# remaining residents, and Phase 3 repeatedly subsamples the equilibrium 
# community over larger extents and tests coexistence of the subsampled community.
#
# The model is a discrete space lottery model, where multiple 
# populations may occupy the same site in varying proportions. The model allows 
# dispersal of seeds based on a dispersal kernel, and competition from 
# neighboring sites according to a competition kernel. Sedentary adults do not 
# compete, but experience a  constant, density-independent rate of survival.   
#
# Currently set for short-range dispersal (and competition), and long range
# spatial correlation. 
# See section "Variables of model that can be tuned" to change these and other
# parameter values. 

#=========================================================================
## Load these libraries
#=========================================================================

library(fields)

#========================================================================
## Function definitions
#========================================================================
#R equivalent for matlab's Meshgrid
Meshgrid=function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 

#========================================================================
## Variables of model that can be tuned
#========================================================================
#Define basic simulation variables
nspp=30 #Number of starting species
ld=32 #Lattice (dimension) width/height -- This should be even
xx=matrix(seq((-ld/2),(ld/2))) #Grid coordinates
np=length(xx)
gens=100 #Number of generations for invasion 
burns=400 #Let the residents establish
genst=burns+gens

#Species characteristics
#Survival
sr=array(0.7,dim=c(np,np,nspp)) 

#interspecific competition
alphaij=1

#Mean reproduction: generated randomly such that each species has a similar 
#fitness, but with a chance of it being higher or lower. 
mn1=matrix(rnorm(nspp,mean=2, sd=0.01),nspp,1)

#Spatial variance and Correlation lengths. Effectively the autocorrelation. 
#As the correlation length becomes longer, more environmental variation is 
#found at relatively larger distances. This is symmetrical across species
vps=1
ziep=20 #Spatial correlation length
###

###Competition and Dispersal parameters for exponential (Laplacian) kernel: 
# values of the parameters produce shorter distances  
#Competition
#brr=1/(100*np) #Essentially global
brr=1 #short range

#Dispesal
#arr=1/(100*np) #essentially global
arr=1 #short range
###

#========================================================================
# Internal variabls
#========================================================================

# 2D 
#Make the full array of space-time coordinates (Here: for 2D space, but no 
#time dimension) Meshgrid will produce an R data frame
stc.temp=Meshgrid(xx,xx)
#Convert the coordinates into an array 
stc=array(c(matrix(stc.temp$x,np,np,byrow=T),matrix(stc.temp$y,np,np,byrow=T)),dim=c(np,np,2)) 
#Final object is two matrixes where entries in matrix 1 correspond to the x 
#coordinate of each site, and matrix 2 corresponds to y coordinates of each site. 


#Population matrixes. Each population gets to save only 2 time points 
#(the current, and the future). This is to save memory. Otherwise, we're 
#storing ld^2*nspp*genst entries. One matrix is for current population, 
#the other is for the next time step. 
nr1=array(data=0,dim=c(np,np,nspp))
nr1.nxt=array(data=0,dim=c(np,np,nspp))

#Random ICs 
#Perpetually initiate invasions across the lattice by letting each species 
#have a very low background invasion rate. Set this to be large 
#enough that invasion actually happens in a reasonable amount of time, but 
#does not influence coexistence.
beta=10^-6

#Reproduction. This is the environmentally-driven term (e.g.  germination in
#the annual plant model). There are many ways to do do this. Here we define
#an environmental template for species' intrinsic reproductive response.

Fr1= array(data=0, dim=c(np,np,nspp))
covFr1= array(data=0, dim=c(np,np,nspp))

#Calculate the covariance structure of Fr1 with the following. Note, 
#this assumes symmetrical x and y directions. 
for (i in 1:nspp){
covFr1[,,i]=vps*exp(-abs(stc[,,1])/ziep-abs(stc[,,2])/ziep)}

#Now make each species' spatial map of reproduction, Fr1, by choosing one 
#region in space for each species. This is region is defined by taking an 
#initial random seed and using the spatial covariance function to map out 
#the adjacent environmental response.  
ics=1 #1 for each species
ics.sites=matrix(0,nspp,3)
for (j in 1:nspp) { 
ics.sites[j,]=cbind(sample(matrix(1:np,np,np), ics, replace=F), sample(matrix(1:np,np,np), ics, replace=F),j)}
#Make these are unique so multiple species don't get the same starting site
ics.sites=ics.sites[!duplicated(ics.sites[,1:2]),]
fr1.tmp=array(data=0,dim=c(np,np,nspp))
fr1.tmp[ics.sites]=1
for (j in 1:nspp){
Fr1[,,j] = Re(fft(fft(fr1.tmp[,,j]*mn1[j])*fft(covFr1[,,1] ) ,inverse=T)/np^2)
}
#If you want to see the spatial map of species' reproduction
image.plot(apply(Fr1,c(1,2),sum))

#Resulting spatial means/medians. 
#b1 = matrix(0,np^2,nspp)
#for(r in 1:nspp) {b1[,r] = matrix(Fr1[,,r],np^2,1)} 
#cov_hist=matrix( cov(b1)[lower.tri(cov(b1))], sum(lower.tri(cov(b1))),1)
#mean(cov_hist)
#median(cov_hist)

### Competition and dispersal kernels. These are both exponential functions. 
# These are all symmetrical -- all intra/interspecific competition has the
# same kernels, all dispersal kernels identical

#Competition kernel 
kc= brr/2*exp(-brr*(abs(stc[,,1])+abs(stc[,,2])) )
kc=kc/(sum(kc)) #Normalize so that integral =1 
fkc=fft(kc)#/(np+1) #Fourier transform

#Dispersal kernel 
kd= arr/2*exp(-arr*(abs(stc[,,1])+abs(stc[,,2])))
#kd=matrix(1,dim(xx)[1], dim(xx)[2]) #Uniform dispersal
kd=kd/(sum(kd))
fkd=fft(kd)#/(np+1)

#========================================================================
#Phase1: Population dynamics to generate an equilibrium community 
#========================================================================

lam1 = array(0, dim=c(np,np,nspp) ) #Reproduction after competition
nr.disp1 = array(0, dim=c(np,np,nspp) ) #Dispersed seeds
plot.top = array(0, dim=c(np,np,nspp) ) #Plotting (most abundant species per site)
total.pop = matrix(0, genst, nspp) #Summary variable: total pop of species

### Burn first

for (n in 1:(burns)) {
nr1.nxt=array(data=0,dim=c(np,np,nspp))
if (n>=2) { beta = 0} #This is to cut off invasion at a certain point
#Determine the seed rain at each site across all species. 
#First, reproduction is reduced by competition. All heterospecifics have the 
#same effect, and thus are lumped as a second species. 
	for (s in 1:nspp) {

		
		fp1=fft(Fr1[,,s]*nr1[,,s]) #/(np+1) #fft of population
		fp2=fft(apply(Fr1*nr1,c(1,2),sum)-Fr1[,,s]*nr1[,,s]) #fft of all 
		#heterospecific populations.
		
		#Inverse fft of convolution of population and competition kernel -- 
		#gives competition species 1 experiences from itself
		Cr11=Re(fft((fp1*fkc),inverse=T)/(np+1))
		Cr11=cbind(Cr11[,ceiling(np/2):(np)],Cr11[,1:floor(np/2)]) #These have 
		#to get chopped up and rearranged due to R's default implementation of 
		#fft. 
		Cr11=rbind(Cr11[ceiling(np/2):(np),],Cr11[1:floor(np/2),]) 
		
		#From heterospecifics
		Cr12=Re(fft((fp2*fkc),inverse=T)/(np+1)) 
		Cr12=cbind(Cr12[,ceiling(np/2):(np)],Cr12[,1:floor(np/2)])
		Cr12=rbind(Cr12[ceiling(np/2):(np),],Cr12[1:floor(np/2),])
		
		#Total reproduction after competition
		lam1[,,s] = (nr1[,,s]*Fr1[,,s])/(Cr11+alphaij*Cr12) 
		#Clean up small negatives left over from fft
		lam1[,,s] = lam1[,,s]*(nr1[,,s]>1e-34) 
		#Nans are produced when populations are 0. 
		lam1[,,s][is.na(lam1[,,s])] = 0 
		
		#Seeds disperse	
		fd1 = fft(lam1[,,s])	
		nr.disp1.tmp=Re(fft((fd1*fkd),inverse=T)/(np+1)) #Convolution with dispersal kernel
		nr.disp1.tmp= cbind(nr.disp1.tmp[,ceiling(np/2):(np)],nr.disp1.tmp[,1:floor(np/2)])
		nr.disp1.tmp=rbind(nr.disp1.tmp[ceiling(np/2):(np),],nr.disp1.tmp[1:floor(np/2),])	
		
		nr.disp1[,,s]=nr.disp1.tmp+beta

		nr1.nxt[,,s] = nr.disp1[,,s]+nr1[,,s]*sr[1]
		total.pop[n,s] = sum(nr1.nxt[,,s])
		}

		#Wash, rinse, repeat
		nr1=nr1.nxt

		#Make a picture every 500 or so steps
		if (n %% 100 == 0){ print( paste( 'n = ', n))
		print(total.pop[n,]) }
		#nr1_pic=array(0, dim=c(np,np,nspp) )
		#for (s in 1:nspp) {
		#nr1_pic[,,s]=nr1[,,s]*s}
		#image.plot(apply(nr1_pic,c(1,2),sum)) }

}

#Use the final community as the equilibrium community to test for coexistence. 
#That is, this community is the one that will be used for the repeated invasion
#attempts below.

nr1.test=nr1

#==============================================================================
#Phase 2: Invasion for every species against community. 
#Species are taken in turn from the community generated above and allowed to 
#invade. Invasion is measured using the slope of the invasion growth rate (IGR).
#Roughly speaking, all species with positive IGR can coexist with  the final 
#community. 
#==============================================================================

#Variable to store slopes of invasion growth rates
igr.slope=matrix(0,nspp,1)

#Run this for all species
for (k in 1:nspp){

#Initialize the resident community
nr1=nr1.test
nr1[,,k]= 1e-6/np^2 #Reset invader to low density
total.pop[(burns),k] = sum(nr1[,,k])

#Run population dynamics for limited number of generations. This code is the 
#same as in the previous section of the code.  
for (n in (burns+1):(genst)) {
nr1.nxt=array(data=0,dim=c(np,np,nspp))

#Determine the seed rain at each site across all species. 
#First, reproduction reduced by competition. All heterospecifics have the same 
#effect, and thus are lumped as a second species. 
	for (s in 1:nspp) {

		
		fp1=fft(Fr1[,,s]*nr1[,,s]) #/(np+1) #fft of population
		fp2=fft(apply(Fr1*nr1,c(1,2),sum)-Fr1[,,s]*nr1[,,s]) #fft of all 
		#heterospecific populations.
		
		#Inverse fft of convolution of population and competition kernel -- 
		#gives competition species 1 experiences from itself
		Cr11=Re(fft((fp1*fkc),inverse=T)/(np+1))
		Cr11=cbind(Cr11[,ceiling(np/2):(np)],Cr11[,1:floor(np/2)])  #These have 
		#to get chopped up and rearranged due to R's default implementation of 
		#fft. 
		Cr11=rbind(Cr11[ceiling(np/2):(np),],Cr11[1:floor(np/2),])
		
		#1 from heterospecifics
		Cr12=Re(fft((fp2*fkc),inverse=T)/(np+1)) 
		Cr12=cbind(Cr12[,ceiling(np/2):(np)],Cr12[,1:floor(np/2)])
		Cr12=rbind(Cr12[ceiling(np/2):(np),],Cr12[1:floor(np/2),])
		
		#Total reproduction after competition
		lam1[,,s] = (nr1[,,s]*Fr1[,,s])/(Cr11+alphaij*Cr12) 
		#Clean up small negatives left over from fft
		lam1[,,s] = lam1[,,s]*(nr1[,,s]>1e-34) 
		#Nans are produced when populations are 0.
		lam1[,,s][is.na(lam1[,,s])] = 0  

		#Seeds disperse	
		fd1 = fft(lam1[,,s])	
		nr.disp1.tmp=Re(fft((fd1*fkd),inverse=T)/(np+1)) #Convolution with dispersal kernel
		nr.disp1.tmp= cbind(nr.disp1.tmp[,ceiling(np/2):(np)],nr.disp1.tmp[,1:floor(np/2)])
		nr.disp1.tmp=rbind(nr.disp1.tmp[ceiling(np/2):(np),],nr.disp1.tmp[1:floor(np/2),])	
		
		nr.disp1[,,s]=nr.disp1.tmp+beta

		nr1.nxt[,,s] = nr.disp1[,,s]+nr1[,,s]*sr[1]
		
		}
		total.pop[n,k] = sum(nr1.nxt[,,k])

		#Wash, rinse, repeat
		nr1=nr1.nxt

		#Update user and possibly make a picture every 100 or so steps
		if (n %% 100 == 0){ print( paste( 'n = ', n, 'k = ', k))
		print(total.pop[n,]) }
		#nr1_pic=array(0, dim=c(np,np,nspp) )
		#for (s in 1:nspp) {
		#nr1_pic[,,s]=nr1[,,s]*s}
		#image.plot(apply(nr1_pic,c(1,2),sum)) }

}

#Calculate the slopes of the IGRs 
t1=(burns+(genst-burns)/2):genst
a1=log(total.pop[t1,k])
igr.slope[k] =lm(a1~t1)$coefficients[2]

#k
}

plot(total.pop[(burns+1):genst,1], type="l")
cl=sample(rainbow(nspp)) #Choose a bunch of colors automatically
for (n in 2:nspp) {
	lines(total.pop[(burns+1):genst,n], col = cl[n])}

#Save useful data at the end of the run
save (file= "lott30spp32np1_lott_cont_shortlong4.var", "nr1","Fr1","covFr1","total.pop","igr.slope")

#===============================================================================
## Phase 3: Subsampling for species area and coexistence curves
#===============================================================================
#Load from file for this next run in order to keep sampling from same large grid
#load( "lott30spp32np1_lott_cont_shortlong4.var" )

#Loop over subsample sizes. These define a square "radius."" Samples are made by 
#choosing center coordinates, then counting out in each direction from there.

min.grid = 1 #smallest grid
max.grid = 15 #largest grid -- set according to spatial extent
grid.inc = 1 #increments
ngrids.max=2*floor(ld/(min.grid*2+1)) #Maximum number of grids to use 

#Species area curve
sa.curve=matrix(NaN,ngrids.max, (ceiling((max.grid-min.grid+1)/grid.inc))) 

#There are multiple ways to automatically test/summarize invasion success. These 
#three were compared and all give approximately equal results. 

#Simple yes or no based on sign of slope at a particular time point
ca.curve.yesno=matrix(NaN,ngrids.max, (ceiling((max.grid-min.grid+1)/grid.inc))) 
#Average slope through time
ca.curve.ave=matrix(NaN,ngrids.max, (ceiling((max.grid-min.grid+1)/grid.inc)))
#
ca.curve.inv=matrix(NaN,ngrids.max, (ceiling((max.grid-min.grid+1)/grid.inc)))

burns = 0 #Number to burn and equilibrate
ngens = 100 #Number of time steps for invasion attempts.
ngenst= burns+ngens #When to initialize first invasion
p=1 #Number of invasion attempts allowed at each subgrid size.

#Pare this down so that only species with successful invasion against the 
#community are used: 

nr1a= nr1[,,igr.slope>0]
Fr1a= Fr1[,,igr.slope>0]

nspp=dim(nr1a)[3] #Number of species that actually coexisted
gindex=1

#Outer loop over grid sample size. Note: much of the code below is a repeat of 
#code in Phase 1. 
for (g in seq(min.grid, max.grid, grid.inc) ){ 
	
	#Generate coordinates for subgrids to sample. 
	ngrids=2*floor(ld/(g*2+1)) #Use about 3/5 of total possible unique centers? 	
	#The sampling space is adjusted to avoid having to implement periodic bcs.
	g.centers=cbind(sample(matrix((1+g):(np-g),(np-2*g),(np-2*g)), ngrids, replace=F), sample(matrix((1+g):(np-g),(np-2*g),(np-2*g)), ngrids, replace=F))	
	###############################################################
	#Remake stc and the dispersal and competition kernels 
	###############################################################

	xx=matrix(seq(-g,g)) #Grid coordinates
	# 2D 
	#Make the full array of space-time coordinates (Here: for 2D space, no time)
	#Meshgrid will produce an R data frame
	stc.temp=Meshgrid(xx,xx)
	#Convert the coordinates into an array for easier access
	stc=array(c(matrix(stc.temp$x,(g*2+1),(g*2+1),byrow=T),matrix(stc.temp$y,(g*2+1),(g*2+1),byrow=T)),dim=c((g*2+1),(g*2+1),2)) 
	#Final object is two matrixes where entries in matrix 1 correspond to the x 
	#coordinate of each site, and matrix 2 corresponds to y coordinates of each 
	#site. 
	###################################################################

	### Competition and dispersal kernels. These are also exponential functions. 
	# Just make these all symmetrical for now -- all intra/interspecific competition with
	# same kernels, all dispersal kernels identical

	#Competition kernel 
	kc= brr/2*exp(-brr*(abs(stc[,,1])+abs(stc[,,2])) )
	kc=kc/(sum(kc)) #Normalize so that integral =1 
	fkc=fft(kc)#/(np+1) #Fourier transform

	#Dispersal kernel 
	kd= arr/2*exp(-arr*(abs(stc[,,1])+abs(stc[,,2])))
	#kd=matrix(1,dim(xx)[1], dim(xx)[2]) #Uniform dispersal
	kd=kd/(sum(kd))
	fkd=fft(kd)#/(np+1)

	#Second loop through each of the selected subgrids
	for (d in 1:ngrids) { 
		#Extract the sub grid
		nr1.samp = nr1a[(g.centers[d,1]-g):(g.centers[d,1]+g),(g.centers[d,2]-g):(g.centers[d,2]+g),]	
		#Count the number of species
		nr1.sa=nr1.samp
		#This loop is to identify the species that are present by tagging them
		for (s in 1:nspp) {
		nr1.sa[,,s]=(nr1.samp[,,s]>0.2)*s}
		spp.ids=unique ( matrix( nr1.sa,((g*2+1)^2*nspp),1)) #Species IDs
		spp.ids=spp.ids[-which(spp.ids == 0)]
		sa.curve [d,gindex] = length(spp.ids) #The sample for the SA curve

		#Now a loop to measure coexistence. Repeated invasion attempts for each 
		#of the species that are present in the subgrid.
			
		#First, loop over the invaders. 
		nspp.inv=length(spp.ids)

		#Keep track of invasion successes/failure.
		invasion.yesno=matrix(0,p,nspp.inv)
		#Average number of individuals following invasion 
		invasion.ave=matrix(0,p,nspp.inv)
		#Slope of IGR
		invasion.slope=matrix(0,p,nspp.inv) 
		#Keep track of total number of individuals of each species through time
		total.pop2 = matrix(0, ngenst, nspp) 

		#Third loop goes through each species in the subgrid and allows it to 
		#invade against the remaining community in that subgrid.
		for(k in 1:nspp.inv) {
		spp.inv= spp.ids[k]		
		
			#Repeated invasion attempts to build up stats, if desired
			for (c in 1:p) {		
		
			print( paste("gindex =", gindex, ", g = ", g, ", d =", d, ", k = ", k, " of ", nspp.inv, ", c =", c))

			#Take the Fr1 (reproduction as function of environment) subsample
			Fr1.samp = Fr1a[(g.centers[d,1]-g):(g.centers[d,1]+g),(g.centers[d,2]-g):(g.centers[d,2]+g),]
		
			#Initialize the matrix to the subsample of the population matrix, 
			#but remove invader		
			nr1.test=nr1.samp
			nr1.test[,,spp.inv]=matrix(0, (g*2+1),(g*2+1))
			
			
			################################################
			#Population dynamics 
			#
			################################################

			#Reproduction after competition
			lam1 = array(0, dim=c((g*2+1),(g*2+1),nspp) ) 
			#Dispersed seeds
			nr.disp1 = array(0, dim=c((g*2+1),(g*2+1),nspp) )
			#Dominant species per site
			plot.top = array(0, dim=c((g*2+1),(g*2+1),nspp) ) 
			
			
			#################################################
			#Invasion 
			#################################################

			#Across the lattice
			nr1.test[,,spp.inv]= 1e-6/((g*2+1)^2)
			total.pop2[1,k] = sum(nr1.test[,,k])

			########################################
			for (n in 2:(ngenst)) {
			nr1.nxt=array(data=0,dim=c((g*2+1),(g*2+1),nspp))
			
			
			#Determine the seed rain at each site across all species. 
			#First, reproduction reduced by competition. All heterospecifics have the same 
			#effect, and thus are lumped as a second species. 
				for (s in 1:nspp) {

					fp1=fft(Fr1.samp[,,s]*nr1.test[,,s]) #/(np+1) #fft of population
					fp2=fft(apply(Fr1.samp*nr1.test,c(1,2),sum)-Fr1.samp[,,s]*nr1.test[,,s]) #/(np+1) #fft of all heterospecific populations
					#fft of all 
					#heterospecific populations.
		
					#Inverse fft of convolution of population and competition kernel -- 
					#gives competition species 1 experiences from itself
					Cr11=Re(fft((fp1*fkc),inverse=T)/((g*2+1)+1))
					Cr11=cbind(Cr11[,ceiling((g*2+1)/2):((g*2+1))],Cr11[,1:floor((g*2+1)/2)])
					Cr11=rbind(Cr11[ceiling((g*2+1)/2):((g*2+1)),],Cr11[1:floor((g*2+1)/2),])  		
		
					#1 from heterospecifics
					Cr12=Re(fft((fp2*fkc),inverse=T)/((g*2+1)+1)) 
					Cr12=cbind(Cr12[,ceiling((g*2+1)/2):((g*2+1))],Cr12[,1:floor((g*2+1)/2)])
					Cr12=rbind(Cr12[ceiling((g*2+1)/2):((g*2+1)),],Cr12[1:floor((g*2+1)/2),])
		
					lam1[,,s] = (nr1.test[,,s]*Fr1.samp[,,s])/(Cr11+alphaij*Cr12) 
					lam1[,,s] = lam1[,,s]*(nr1.test[,,s]>1e-34) 
					lam1[,,s][is.na(lam1[,,s])] = 0 

					#Seeds disperse	
					fd1 = fft(lam1[,,s])	
					nr.disp1.tmp=Re(fft((fd1*fkd),inverse=T)/((g*2+1)+1)) #Convolution with dispersal kernel
					nr.disp1.tmp= cbind(nr.disp1.tmp[,ceiling((g*2+1)/2):((g*2+1))],nr.disp1.tmp[,1:floor((g*2+1)/2)])
					nr.disp1.tmp=rbind(nr.disp1.tmp[ceiling((g*2+1)/2):((g*2+1)),],nr.disp1.tmp[1:floor((g*2+1)/2),])	
					
					nr.disp1[,,s]=nr.disp1.tmp
					
					nr1.nxt[,,s] = nr.disp1[,,s]+nr1.test[,,s]*sr[1]
					
					}
					total.pop2[n,k] = sum(nr1.nxt[,,k])
					#Wash, rinse, repeat
					nr1.test=nr1.nxt
		
			
					#if (n %% 100 == 0){
					#nr1_pic=array(0, dim=c((g*2+1),(g*2+1),nspp) )
					#for (s in 1:nspp) {
					#nr1_pic[,,s]=nr1.test[,,s]*s}
					#image.plot(apply(nr1_pic,c(1,2),sum)) 
					

				}

			#Now use total.pop2 to decide if invasion was successful. We tested 
			#several ways to do this that produce similar results: 
			#1) Are there more than 0 individuals in population? 
			#1) Number of final indivduals > invasion number ?
			#2) log(Slope) after invasion > 0 ?
			thresh.inv = 1e-6 #How many individuals for (1)?
			invasion.yesno[c,k] = as.numeric( total.pop2[(ngenst-1),k] > thresh.inv )
			
			#Approximate slope			
			t1=(2+(ngenst-2)/2):ngenst
			a1=log(total.pop2[t1,k])
			invasion.slope[c,k] =lm(a1~t1)$coefficients[2]
			
			#More naive approximation
			invasion.ave[c,k] = total.pop2[(ngenst-1),k] - thresh.inv
			
			#c			
			}
		#k
		}
	ca.curve.yesno[d,gindex]= sum(colSums(invasion.yesno)>=1)
	ca.curve.ave[d,gindex]=sum(colMeans(invasion.ave)>0)
	ca.curve.inv[d,gindex]=sum(colMeans(invasion.ave)>0)
	#d

	}

	gindex = gindex+1

}

#Save useful data
save(file="lott30spp128np1_lott_cont_shortlong_subdata4.var", "ca.curve.yesno","ca.curve.inv", "ca.curve.ave","sa.curve")

#Generate a rough plot of SA and CA curves. 
#jpeg(file="sa_ca_30spp32np1.jpg",height=4.5, width=4.75, units="in", family='Helvetica', quality=100, res=300, pointsize=12)
par(mar=c(5,6,4,2))
plot((2*seq(1:max.grid)+1)^2,colMeans(sa.curve, na.rm=T),ylim=c(0,20), t="l", xlab="Area", ylab= "Species per Area, Coexistence per Area")
lines((2*seq(1:max.grid)+1)^2,colMeans(ca.curve.ave, na.rm=T),col="blue",xlab="Area", ylab= "Species per Area, Coexistence per Area")
#dev.off()

