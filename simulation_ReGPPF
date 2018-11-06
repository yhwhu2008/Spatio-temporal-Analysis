library(Rcpp)
T <- 1000

load("completeSimulationData0_01SNRNotAlign.RData")
indianaSimulationData = completeSimulationData0_01[,c(5,(ncol(completeSimulationData0_01)-T+2):ncol(completeSimulationData0_01))]

simu = indianaSimulationData[order(indianaSimulationData$NAME10),]

data = t(simu[,-1])
colnames(data) = simu[,1]
centeredData = matrix(0,nrow = nrow(data), ncol = ncol(data))
for (i in 1:ncol(data)){
  centeredData[,i] = (data[,i] - mean(data[,i]))/sd(data[,i])
}


shift = 0
centeredData = centeredData + shift

timeUnit = 7
numComp = 5 # number of selected components
lagNumClip = 0
N = ncol(centeredData)

numClip = floor(nrow(centeredData)/timeUnit) - lagNumClip
dataByUnit = array(0, dim = c(timeUnit, ncol(centeredData),numClip), dimnames = list(1:timeUnit, 1:ncol(centeredData),1:numClip))
for (i in 1:numClip){
  dataByUnit[,,i] = centeredData[((i-1)*timeUnit+1):(i*timeUnit),]
}

centerdataByUnit = array(0, dim = c(timeUnit, ncol(centeredData),numClip), dimnames = list(1:timeUnit, 1:ncol(centeredData),1:numClip))
for (i in 1:numClip){
  for (j in 1:N){
    centerdataByUnit[,j,i] = (dataByUnit[,j,i] - mean(dataByUnit[,j,i]))
  }
}

eigenvector = array(0, dim = c(N, numComp, numClip), dimnames = list(1:N, 1:numComp, 1:numClip))

for (i in 1:(numClip)){
  eigenvector[,,i] = prcomp(as.matrix(centerdataByUnit[,,i]))$rotation[,1:numComp]
}

compH = matrix(0,nrow = numComp, ncol = numClip)
for (i in 1:nrow(compH)){
  for (j in 1:ncol(compH)){
    compH[i,j] = eigenvector[,i,j] %*% colMeans(centerdataByUnit[,,j])
  }
}


initial_eigenvalue = (as.vector(prcomp(as.matrix(centerdataByUnit[,,1])))$sdev)^2
a_0 = initial_eigenvalue[1:numComp] + rnorm(5,0,1)


sum(initial_eigenvalue[1:numComp])/sum(initial_eigenvalue)
dd = matrix(0, nrow = numClip, ncol = ncol(centeredData))
for (i in 1:ncol(centeredData)){
  dd[,i] = colMeans(centerdataByUnit[,i,])
}
train_data = dd
D = nrow(train_data)  # days
N = ncol(train_data)  # counties
M = 10000 # 1086s
K = numComp

adjH = rep(0,N*K*D)
for (j in 0:(N-1)){
  for (k in 0:(K-1)){
    for (t in 0:(D-1)){
      adjH[K*D*j+k*D+t + 1] = eigenvector[j+1,k+1,t+1]
    }
  }
}

sourceCpp("gibbs_pca4.cpp")

ptm = proc.time()
set.seed(1234)
re = gibbs_pca4(M, K, train_data, adjH, a_0)  # 593s M = 10000  K = 4
sga = 1
proc.time() - ptm


aaa = re$a
sigma2_mean = colMeans(re$sigma2[5001:10000,])


xxi = re$xi
xi_final = array(0, dim = c(M,K,D))
for (m in 0:(M-1)){
  for (k in 0:(K-1)){
    xi_final[m+1,k+1,] = xxi[(K*D*(m)+D*k+1):(K*D*(m)+D*(k+1))]
  }
}

xi_final_1_mean = colMeans(xi_final[5001:10000,1,])
xi_final_2_mean = colMeans(xi_final[5001:10000,2,])
xi_final_3_mean = colMeans(xi_final[5001:10000,3,])
xi_final_4_mean = colMeans(xi_final[5001:10000,4,])
xi_final_5_mean = colMeans(xi_final[5001:10000,5,])



data_topred = centeredData[(dim(centerdataByUnit)[1]*dim(centerdataByUnit)[3]+1):nrow(centeredData),]
theta = rep(0,ncol(centeredData))
for (i in 1:length(theta)){
  theta[i] = sum(aaa[,numClip]*eigenvector[i,,numClip]*(1+colMeans(xi_final[5001:10000,,numClip])))
  # theta[i] = sum(aaa[,numClip]*eigenvector[i,,numClip]*(colMeans(xi_final[5001:10000,,numClip])))
}

data_pred = matrix(0, nrow(data_topred),ncol(data_topred))
for (i in 1:ncol(data_pred)){
  data_pred[,i] = rnorm(nrow(data_topred),theta[i], sqrt(sigma2_mean[i]))
}

J = nrow(data_pred)
scaler = 1
MSPE_ReGPCA = NULL
for (j in 1:J){
  MSPE_ReGPCA[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
}


## filter
##########################################################
dataObserv = t(aaa)
write.table(dataObserv, row.names = FALSE, col.names = FALSE, sep = ',', file = paste('ReGPCA_shift',timeUnit,'C',numComp,'.csv',sep = ''))
dataObserv = read.csv(file = paste('ReGPCA_shift',timeUnit,'C',numComp,'.csv',sep = ''),  header = FALSE)         
x1=2.5892;    x2=0.4329;   x3=-1.5486;   x4=-0.4575;     x5=1; x6=0; x7=1 # snr = 0.01

diffusion = x7*diag(1,ncol(dataObserv))
# specify the function in state equation
stateFunction = function(x){ # x: state
  # y = diffusion%*%x + x*(1-x/x8)
  y = x
  return(y)
}

mu = function(x, y){ # x: state, y: observation
  z = solve(sigmaInverse) %*% (solve(stateNoiseVariance) %*% stateFunction(x) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% as.numeric(y))
  return(z)
}

phi = function(x, y){ # x: state, y: observation
  z = 1/2 * t(as.numeric(y) - observe2StateMatrix %*% stateFunction(x)) %*% solve(K) %*% (as.numeric(y) - observe2StateMatrix %*% stateFunction(x))
  return(z)
}

##########################################################
data = dataObserv
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
# stateNoiseVariance = CQ

observeNoiseMean = rep(0, observeDimension)
CR = matrix(x4, observeDimension,observeDimension)
diag(CR) = x3
observeNoiseVariance = t(CR) %*% CR
# observeNoiseVariance = CR

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

# In case that -phi is too Small
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


aaa_f = t(weightedStates)

data_topred = centeredData[(dim(centerdataByUnit)[1]*dim(centerdataByUnit)[3]+1):nrow(centeredData),]
theta = rep(0,ncol(centeredData))
for (i in 1:length(theta)){
  for (k in 1:nrow(aaa_f)){
    theta[i] = theta[i] + aaa_f[k,ncol(aaa_f)]*eigenvector[i,k,ncol(aaa_f)]*(1+colMeans(xi_final[5001:10000,k,])[ncol(aaa_f)])
  }
}

data_pred = matrix(0, nrow(data_topred),ncol(data_topred))
for (i in 1:ncol(data_pred)){
  data_pred[,i] = rnorm(nrow(data_topred),theta[i], sqrt(sigma2_mean[i]))
  
}

J = nrow(data_pred)
scaler = 1
MSPE_ReGPCA_Filter = NULL
for (j in 1:J){
  MSPE_ReGPCA_Filter[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
}





# Full Bayesian, paramters and state same time
############################
ptm = proc.time()
set.seed(1234)
library(dlm)
data = read.csv("ReGPCA_shift7C5.csv", header = FALSE)
ncomp = 5
weightedStates = matrix(0,nrow(data)+1,ncomp)

VV = rep(0,ncomp)
WW = rep(0,ncomp)
for (icol in 1:ncomp){
  outGibbs <- dlmGibbsDIG(data[,icol], dlmModPoly(1, dV = var(data[,icol]), dW = var(data[,icol])),
                          a.y = var(data[,icol]), b.y = 1, a.theta = var(data[,icol]), b.theta = 1,
                          n.sample = 1000,
                          save.states = TRUE)
  burn <- 500
  
  dV <- outGibbs$dV[-(1:burn)]
  dW <- outGibbs$dW[-(1:burn)]
  
  postMean = mcmcMean(with(outGibbs, sqrt(cbind(dV, dW))))
  VV[icol] = as.vector(postMean)[1]
  WW[icol] = as.vector(postMean)[3]
  
  aa = matrix(outGibbs$theta,dim(outGibbs$theta)[1],dim(outGibbs$theta)[3]) 
  weightedStates[,icol] = rowMeans(aa)
}
proc.time() - ptm



## filter
##########################################################

dataObserv = read.csv("ReGPCA_shift7C5.csv", header = FALSE)
aaa = t(dataObserv)
aaa_f = t(weightedStates)

data_topred = centeredData[(dim(centerdataByUnit)[1]*dim(centerdataByUnit)[3]+1):nrow(centeredData),]
theta = rep(0,ncol(centeredData))
for (i in 1:length(theta)){
  for (k in 1:nrow(aaa_f)){
    theta[i] = theta[i] + aaa_f[k,numClip]*eigenvector[i,k,numClip]*(1+colMeans(xi_final[5001:10000,k,])[numClip])
  }
}

data_pred = matrix(0, nrow(data_topred),ncol(data_topred))
for (i in 1:ncol(data_pred)){
  data_pred[,i] = rnorm(nrow(data_topred),theta[i], sqrt(sigma2_mean[i]))
  
}

J = nrow(data_pred)
scaler = 1
MSPE_ReGPCA_Filter_B2 = NULL
for (j in 1:J){
  MSPE_ReGPCA_Filter_B2[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
}



plot(MSPE_ReGPCA,type = "b",xlab = "Days",ylab = "MSPE",ylim = c(min(MSPE_ReGPCA,MSPE_ReGPCA_Filter,MSPE_ReGPCA_Filter_B2), max(MSPE_ReGPCA,MSPE_ReGPCA_Filter,MSPE_ReGPCA_Filter_B2)+0.1))
lines(MSPE_ReGPCA_Filter, type = "b", col = "red")
lines(MSPE_ReGPCA_Filter_B2, type = "b", col = "green4")
legend("topright",c("ReGP","ReGPPF","ReGPPFB"),horiz = T,col = c("black","red","green"),lty = rep(1,3), pch=rep(1,3),lwd=rep(1,3),cex=1)

