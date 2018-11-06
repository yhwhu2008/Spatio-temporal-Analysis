## 1 generate simulated data with SNR 
################################################
# generate data for WinBUGS
source("lambda-fun-level.R")
library(sp)
library(maptools)
indiana <- readShapePoly("tl_2010_18_county10.shp"[1], proj4string=CRS("+proj=longlat +datum=NAD27"))
library(spdep)
innb<-poly2nb(indiana)
#get 0-1 adjacency matrix for Indiana
in.mat<-nb2mat(innb, style="B")
n.nodes <- length(innb)

#normalized adjacency matrix for IN
in.adjnorm<-nb2mat(innb)

## get IN population and area
inpop <- read.csv(file="INpop.csv", header = TRUE, sep = ",")#[1:92,]
popID = match(indiana@data$NAME10, inpop$county)
pop <- inpop[popID,2]
indiana@data$pop <- pop

# log transform
lpop <- log(pop)
#center them
Xpop <- lpop-mean(lpop)

###############################
#   Simulate delta process    #
###############################

n.iter <- 500
T <- 1000
#start day of outbreak
Tb <- 10
adj0 <- in.mat
ps <- 0.001
pc <- 0.1
pop1 <- lpop/10

# generating the Poisson counts, insert outbreak points
set.seed(110)
delta = read.csv("deltaNotAlign.csv",header = F)
# set.seed(110)
# delta = delta[sample(1:nrow(delta),replace = F),]
###############################
#Simulate from the full model #
###############################
# set initial values
n.iter <- 200
rho1 <- 0.5
rho2 <- 0.5
sigma.mu <- 0.1
sigma.lam <- 0.5
c <- 2
bmu1 <- 1
blam1 <- 1

# generating the Poisson counts
bslevel <- 5
snr <- 0.1  ################# small SNR
set.seed(106)

mpop <- matrix(rep(Xpop, T), nc=T)
#mland <- matrix(rep(Xland, T), nc=T)
beta.dow <- 0.2

# generate mu and lambda for fever
mu.o <- matrix(rnorm(n.nodes*T, mean=log(bslevel), sd=sigma.mu), nr=n.nodes, nc=T)+ bmu1*mpop #+ bmu2*mland
lambda.o <- lambda.fun(adj0=in.mat, n.nodes=n.nodes, n.iter=n.iter, T=T, rho1=rho1, rho2=rho2, sigma=sigma.lam, level=bslevel, c=c) + blam1*Xpop #+ blam2*Xland

dowtp <- rep(c(0,1,1,1,1,1,0), length.out=T)
plot.ts(dowtp, pch=2) #plot the DOW effect
dow <- matrix(rep(dowtp, n.nodes), nc=T, byrow=TRUE)
mu <- mu.o + beta.dow*dow
lambda <- lambda.o + beta.dow*dow

# y.mu <- apply(as.matrix(exp(mu)), c(1,2), rpois, n=1)
# y.lam <- apply(as.matrix(exp(lambda)), c(1,2), rpois, n=1)
# y <- y.mu + delta*y.lam*snr

y = apply(as.matrix(exp(mu+delta*snr*lambda)), c(1,2), rpois, n=1)
y.mu <- apply(as.matrix(exp(mu)), c(1,2), rpois, n=1)

## time series plot of y and y.mu for fever
plot.ts(colSums(y[,1:T]), ylab="", type="o", main=paste("beta.dow =",beta.dow))
lines(colSums(y.mu[,1:T]), lty=20, col="red")
legend("topleft", c("Y","Y.mu"), lty=c(19,20), col=c("black", "red"))
proc.time() - ptm



sim.name<-c(paste("Y",1:T,sep=""))
for (i in 2:T){
  indiana@data[,sim.name[i]] <- y[,i]/pop
  #nc@data[,sim.name[i]] <- y[,i]
}

completeSimulationData = indiana@data
pureSimulationData = 10^5*t(completeSimulationData[,(ncol(completeSimulationData)-T+2):ncol(completeSimulationData)])

save(completeSimulationData,file = "completeSimulationDataNotAlign.RData")
write.csv(pureSimulationData,file="pureSimulationData.csv")

#########################################################


## 2 fit DP, obtain clusters
T <- 1000
library("DPpackage")
load("completeSimulationData.RData")

indianaSimulationData = completeSimulationData[,c(5,(ncol(completeSimulationData)-T+2):ncol(completeSimulationData))]
simu = indianaSimulationData[order(indianaSimulationData$NAME10),]
completedata = matrix(0,nrow(simu)*(ncol(simu)-1),3)
colnames(completedata) = c("county","time","incidenceRate")
for (i in 1:(ncol(simu)-1)){
  completedata[((i-1)*nrow(simu)+1):(i*nrow(simu)),"county"] = as.character(simu$NAME10)
  completedata[((i-1)*nrow(simu)+1):(i*nrow(simu)),"time"] = i
  completedata[((i-1)*nrow(simu)+1):(i*nrow(simu)),"incidenceRate"] = round(as.numeric(simu[,i+1])*10^5)
}

library(prodlim)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

modeMonth=NULL


fit_save_state_byCountyM=array(0, dim=c(1,5000,100))

ptm=proc.time()

y=as.numeric(completedata[,"incidenceRate"])
# x1=cos(2*pi*completedata$days/365)
x1=as.numeric(completedata[,"time"])
z=completedata[,"county"]

numx=1

# Prior information

beta0<-rep(0,numx)
Sbeta0<-diag(1000,numx)
tinv<-diag(1,1)

prior<-list(a0=2,b0=100,nu0=var(y),tinv=tinv,mub=rep(0,1),Sb=diag(1000,1),beta0=beta0,Sbeta0=Sbeta0) # mode 24

# Initial state
state <- NULL
# MCMC parameters
nburn <- 5000
nsave <- 5000
nskip <- 0
ndisplay <- 1000
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

fit <- DPglmm(fixed=y~x1,random=~1|z, family=poisson(log),prior=prior,mcmc=mcmc, state=state,status=TRUE)
data=data.frame(fit$save.state)
fit_save_state_byCounty=matrix(unlist(data),nrow(data),ncol(data))
colnames(fit_save_state_byCounty)=colnames(data)

proc.time()-ptm

save(fit_save_state_byCounty, file = "fit_save_state_byCountySimulationData.RData")
load("fit_save_state_byCountySimulationData.RData")

# rescale color
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
library(maps)

# use match
Weekcolor=array(0,dim=c(1,5000,92))

w1=data.frame(fit_save_state_byCounty)
w1mode=w1[which(w1$thetasave.ncluster==Mode(w1$thetasave.ncluster)),((ncol(w1)-92)):(ncol(w1)-1)]
for (i in 1:nrow(w1mode)){
  initialLabel=as.numeric(w1mode[i,])
  sortLabel=sort(unique(as.numeric(unlist(initialLabel))))
  for (k in 1:length(sortLabel)){
    matchposition=match(initialLabel,sortLabel[k])
    cols=which(matchposition==1)
    Weekcolor[1,i,cols]=k
  }
}



totalColor=Weekcolor
temp=totalColor
rescaletotalColor=array(0,dim=dim(totalColor))

for (j in 1:5000){
  temp1=temp[1,j,]
  label=sort(unique(temp1))
  for (h in 1:length(label)){
    firstPoint=na.omit(temp1)[1]
    matchposition=match(temp1,firstPoint)
    cols=which(matchposition==1)
    rescaletotalColor[1,j,cols]=label[h]
    temp1[cols]=NA
  }
}

#########################
library(maps)
## Constructing a county neighborhood matrix for the state of indiana
##Step1: Create county ID
mn.county = map("county","indiana", fill=TRUE, plot=FALSE)
county.ID <- sapply(strsplit(mn.county$names, ","), function(x) x[2])

library(maptools)
mn.poly = map2SpatialPolygons(mn.county, IDs=county.ID)

##Convert polygon to nb object
# install.packages("spdep")
library(spdep)
mn.nb = poly2nb(mn.poly)
mn.adj.mat = nb2mat(mn.nb, style="B")
## The option style="B" produces the binary adjacency matrix
## Write the 0-1 adjacency matrix
W <- mn.adj.mat
W[(W>0)] = 1

colnames(W)=rownames(W)

# adjust W by the number of neighbours

adjW=W/rowSums(W)
sum(adjW)

clusterWeight=array(0,dim=c(1,5000,50))

temp=rescaletotalColor[1,,]
row_sub = apply(temp, 1, function(row) all(row !=0 ))
temp1=temp[row_sub,]

for (i in 1:nrow(temp1)){
  uniqueCluster=sort(unique(temp1[i,]))
  for (j in 1:length(uniqueCluster)){
    matchposition=match(temp1[i,],uniqueCluster[j])
    cols=which(matchposition==1)
    clusterWeight[1,i,j]=sum(adjW[,cols])
  }
}


vectorSelect=array(0,dim=c(1,5000,92))

i=1
for (j in 1:length(apply(clusterWeight[i,,], 1, function(row) all(row !=0 )))){
  for (h in 1:length(which(clusterWeight[i,j,]!=0))){
    matchposition=match(rescaletotalColor[i,j,],seq(1,length(which(clusterWeight[i,j,]!=0)))[h])
    cols=which(matchposition==1)
    vectorSelect[i,j,cols]=clusterWeight[i,j,h]
  }
}


sumWeightedCluster=matrix(0,5000,1)
maxVectorIndicator=NULL

clusterAfterSelection=matrix(0,1,92)

i=1
sumWeightedCluster[,i]=rowSums(vectorSelect[i,,])
maxVectorIndicator[i]=which.max(sumWeightedCluster[,i])
clusterAfterSelection[i,]=rescaletotalColor[i,which.max(sumWeightedCluster[,i]),]

maxMode=max(clusterAfterSelection)

##################################################
# find corresponding theta for the selected cluster result

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

temp = fit_save_state_byCounty[which(fit_save_state_byCounty[,"thetasave.ncluster"]==Mode(fit_save_state_byCounty[,"thetasave.ncluster"])),]
thetaSelected = temp[maxVectorIndicator,((ncol(temp)-92):(ncol(temp)-1))]
############################
write.csv(unique(thetaSelected),file = "thetaSelected.csv")
#################################



# mapAdjust = c(1:70,74,71,72,73,75:92) for simulated data no need they are the same

cluster = clusterAfterSelection

adjW1 = adjW
colnames(adjW1) = as.numeric(cluster)
spatialDiffusion = matrix(0,nrow(adjW1), length(unique(as.numeric(cluster))))

for(i in 1:ncol(spatialDiffusion)){
  if (length(which(colnames(adjW1)==i))>1){
    spatialDiffusion[,i] = rowMeans(adjW1[,which(colnames(adjW1)==i)])
  }else{
    spatialDiffusion[,i] = adjW1[,which(colnames(adjW1)==i)]
  }
}

rownames(spatialDiffusion) = as.numeric(cluster)

spatialDiffusion1 = matrix(0,length(unique(as.numeric(cluster))), length(unique(as.numeric(cluster))))
for(i in 1:nrow(spatialDiffusion1)){
  if (length(which(rownames(spatialDiffusion)==i))>1){
    spatialDiffusion1[i,] = colSums(spatialDiffusion[which(rownames(spatialDiffusion)==i),])
  }else{
    spatialDiffusion1[i,] = spatialDiffusion[which(rownames(spatialDiffusion)==i),]
  }
}


write.csv(spatialDiffusion1,file = "spatialDiffusion.csv")


simulatedObs1 = round(t(simu[,-1])*10^5)
colnames(simulatedObs1) = as.numeric(clusterAfterSelection)
simulatedObs = matrix(0,nrow(simulatedObs1), length(unique(as.numeric(clusterAfterSelection))))

for(i in 1:ncol(simulatedObs)){
  if (length(which(colnames(simulatedObs1)==i))>1){
    simulatedObs[,i] = rowMeans(simulatedObs1[,which(colnames(simulatedObs1)==i)])
  }else{
    simulatedObs[,i] = simulatedObs1[,which(colnames(simulatedObs1)==i)]
  }
}

write.csv(simulatedObs,file = "afterDP_Simulation.csv")


#####################################################

## use matlab to find MLE

# simulated data with MLE
set.seed(1234)
##############################################################################
##incidence rate by DP regions

data1 = read.csv("afterDP_Simulation.csv",header = F)
data = matrix(as.numeric(unlist(data1)),nrow(data1),ncol(data1))
##############################################
x1= 0.2; x2= 0; x3=5.5;x4=0.1;x5=8.6;x6=3.1;x7=0.00058;x8=1012;
diffusion1 = read.csv("spatialDiffusion.csv",header = F)
diffusion = x7*matrix(as.numeric(unlist(diffusion1)),nrow(diffusion1),ncol(diffusion1))

# specify the function in state equation
stateFunction = function(x){ # x: state
  y = diffusion%*%x + x*(1-x/x8)
  return(y)
}

mu = function(x, y){ # x: state, y: observation
  z = solve(sigmaInverse) %*% (solve(stateNoiseVariance) %*% stateFunction(x) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% y)
  return(z)
}

phi = function(x, y){ # x: state, y: observation
  z = 1/2 * t(y - observe2StateMatrix %*% stateFunction(x)) %*% solve(K) %*% (y - observe2StateMatrix %*% stateFunction(x))
  return(z)
}

##########################################################

observeDimension = ncol(data)
stateDimension = observeDimension
observe2StateMatrix = diag(x5,observeDimension)

diag(observe2StateMatrix[1:(observeDimension-1),2:observeDimension]) = x6
diag(observe2StateMatrix[2:observeDimension,1:(observeDimension-1)]) = x6 


stateNoiseMean = rep(0, stateDimension)
stateNoiseVariance = diag(1/nrow(data), stateDimension, stateDimension) # G
CQ = matrix(x2, stateDimension, stateDimension)
diag(CQ) = x1

stateNoiseVariance = t(CQ) %*% CQ


observeNoiseMean = rep(0, observeDimension)
CR = matrix(x4, observeDimension,observeDimension)
diag(CR) = x3
observeNoiseVariance = t(CR) %*% CR


sigmaInverse = solve(stateNoiseVariance) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% observe2StateMatrix
K = observe2StateMatrix %*% stateNoiseVariance %*% t(observe2StateMatrix) + observeNoiseVariance



library("Matrix")
L = chol(solve(sigmaInverse))

# specify number of state samples M (say draw M=1000 samples of x)

M=10

set.seed(1234)

# with resampling, better results
#######################################################################################################################
statePredict = array(0, dim=c(nrow(data)+1, M, ncol(data)))
library("MASS")
# statePredict[1,,] = mvrnorm(n = M, stateNoiseMean, stateNoiseVariance)
statePredict[1,,] = matrix(0, M, ncol(data))

statePredictResample = statePredict


weights = matrix(0, nrow(data)+1, M)
weights[1,] = 1/M


tempNegativePhi = matrix(0, nrow(data)+1, M)
tempNegativePhi[1,] = 1


# In case that -phi is too small
for (i in 2:(nrow(data)+1)){
  for (j in 1:M){
    statePredict[i,j,] = mu(statePredict[i-1,j,],data[i-1,]) + L %*% matrix(mvrnorm(n = 1, rep(0,stateDimension), diag(1,stateDimension, stateDimension)), stateDimension,1)
    tempNegativePhi[i,j] = -phi(statePredict[i-1,j,], data[i-1,])
  }
  updatePhi = tempNegativePhi[i,] + tempNegativePhi[i-1,]
  maxTemp = max(updatePhi)
  weights[i,] = exp(updatePhi-maxTemp)/sum(exp(updatePhi-maxTemp))
  resamples = sample(1:M, M, replace = TRUE, prob = weights[i,])
  statePredictResample[i,,] = statePredict[i,resamples,]
}

weightsRoundup = round(weights, digits = 3)

weightedStates = matrix(0,nrow(data),ncol(data))
for (i in 1:nrow(data)){
  for (j in 1:ncol(data)){
    weightedStates[i,j] = weighted.mean(statePredictResample[i+1,,j], weights[i+1,])
  }
}
weightedStatesSimulation = weightedStates
save(weightedStatesSimulation,file = "weightedStatesMsimulation.RData")



#############################################
#plot
indianaSimulatedPureData = read.csv("pureSimulationData.csv")
indianaSimulatedPureData = indianaSimulatedPureData[,-1]
indianaSimulatedPureData = matrix(as.numeric(unlist(indianaSimulatedPureData)),nrow(indianaSimulatedPureData),ncol(indianaSimulatedPureData))
data1 = read.csv("afterDP_Simulation.csv",header = F)
data = matrix(as.numeric(unlist(data1)),nrow(data1),ncol(data1))
colnames(data) = seq(1,ncol(data))

period = 12
colorAll = matrix(0,period,92)
for ( i in 1:period){
  interval = round(seq(1+(i-1)*nrow(data)/period, i*nrow(data)/period))
  for (j in 1:nrow(data)){
    cols = match(clusterAfterSelection,j)
    colorAll[i,!is.na(cols)] = mean(indianaSimulatedPureData[interval,!is.na(cols)])
  }
}

library("maps")
mapAdjust = c(1:14,16,17,15,18:43,45,46,44,47:70,74,71,72,73,75:92)
maxCOL=max(colorAll)
minCOL=min(colorAll)
colorkeypval<-list(labels=as.character(c(0.1, 0.2, 0.4, 0.6, 0.8 ,1)),
                   at=(minCOL:maxCOL)/maxCOL, height=1)

# pval <- colorRampPalette(c("white", "yellow","orange","red","darkred"), space = "rgb")
pval <- colorRampPalette(c("white", "grey","black"), space = "rgb")
pvalcols <- pval(maxCOL)

(palette(pvalcols))


plot(seq(1,maxCOL),seq(1,maxCOL),type="p",pch=16,col=seq(1,maxCOL))

par(mfrow=c(3,4),mar=c(0,0,0,0),mex=0.55)
for (i in 1:period){
  tempCol = (colorAll[i,][mapAdjust]-(minCOL-1)*maxCOL/(maxCOL-1))*(maxCOL-1)/(maxCOL-minCOL)
  map('county', 'indiana', fill = TRUE, col = tempCol,mar=c(5,0,0,0))
  #title(main=paste("Simulated Data, SNR=0.05,P",i),cex.main=1.2)
  title(main=paste("Period",i),cex.main=1.2,line=1)
}


##################################################################








# use adjacency matrix and theta values
# simulated data with MLE
set.seed(1234)
##############################################################################
##incidence rate by DP regions

data1 = read.csv("afterDP_Simulation.csv",header = F)
data = matrix(as.numeric(unlist(data1)),nrow(data1),ncol(data1))
##############################################

x1= 0.5; x2= 5.5; x3=0.1;x4=2.7;x5=0.0014;x6=3066;
diffusion1 = read.csv("spatialDiffusion.csv",header = F)
diffusion = x5*matrix(as.numeric(unlist(diffusion1)),nrow(diffusion1),ncol(diffusion1))
theta = read.csv("thetaSelected.csv", header = F)
theta = theta - min(theta) +1
# specify the function in state equation
stateFunction = function(x){ # x: state
  y = diffusion%*%x + x*(1-x/x6)
  return(y)
}

mu = function(x, y){ # x: state, y: observation
  z = solve(sigmaInverse) %*% (solve(stateNoiseVariance) %*% stateFunction(x) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% y)
  return(z)
}

phi = function(x, y){ # x: state, y: observation
  z = 1/2 * t(y - observe2StateMatrix %*% stateFunction(x)) %*% solve(K) %*% (y - observe2StateMatrix %*% stateFunction(x))
  return(z)
}

##########################################################

observeDimension = ncol(data)
stateDimension = observeDimension
observe2StateMatrix = x4*diag(as.numeric(unlist(theta)))

stateNoiseMean = rep(0, stateDimension)
stateNoiseVariance = diag(1/nrow(data), stateDimension, stateDimension) # G
CQ = diag(x1, stateDimension, stateDimension)

stateNoiseVariance = t(CQ) %*% CQ

observeNoiseMean = rep(0, observeDimension)
CR = matrix(x3, observeDimension,observeDimension)
diag(CR) = x2
observeNoiseVariance = t(CR) %*% CR

sigmaInverse = solve(stateNoiseVariance) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% observe2StateMatrix
K = observe2StateMatrix %*% stateNoiseVariance %*% t(observe2StateMatrix) + observeNoiseVariance


library("Matrix")
L = chol(solve(sigmaInverse))
# specify number of state samples M (say draw M=1000 samples of x)
M=10

set.seed(1234)

# with resampling, better results
#######################################################################################################################
statePredict = array(0, dim=c(nrow(data)+1, M, ncol(data)))
library("MASS")
# statePredict[1,,] = mvrnorm(n = M, stateNoiseMean, stateNoiseVariance)
statePredict[1,,] = mvrnorm(n = M, solve(observe2StateMatrix) %*% data[1,], stateNoiseVariance)

statePredictResample = statePredict

weights = matrix(0, nrow(data)+1, M)
weights[1,] = 1/M

tempNegativePhi = matrix(0, nrow(data)+1, M)
tempNegativePhi[1,] = 1


# In case that -phi is too ExSmall
for (i in 2:(nrow(data)+1)){
  for (j in 1:M){
    statePredict[i,j,] = mu(statePredict[i-1,j,],data[i-1,]) + L %*% matrix(mvrnorm(n = 1, rep(0,stateDimension), diag(1,stateDimension, stateDimension)), stateDimension,1)
    tempNegativePhi[i,j] = -phi(statePredict[i-1,j,], data[i-1,])
  }
  updatePhi = tempNegativePhi[i,] + tempNegativePhi[i-1,]
  maxTemp = max(updatePhi)
  weights[i,] = exp(updatePhi-maxTemp)/sum(exp(updatePhi-maxTemp))
  resamples = sample(1:M, M, replace = TRUE, prob = weights[i,])
  statePredictResample[i,,] = statePredict[i,resamples,]
}

weightsRoundup = round(weights, digits = 3)

weightedStates = matrix(0,nrow(data),ncol(data))
for (i in 1:nrow(data)){
  for (j in 1:ncol(data)){
    weightedStates[i,j] = weighted.mean(statePredictResample[i+1,,j], weights[i+1,])
  }
}
weightedStatesSimulation = weightedStates
save(weightedStatesSimulation,file = "weightedStatesMsimulationtheta.RData")








