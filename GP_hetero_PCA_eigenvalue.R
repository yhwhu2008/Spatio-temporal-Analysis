# # 20180402: gibbs_pca4 a_0 = eiganvalue of t_1
# T <- 1000
# 
# # load("R:\\THESIS\\IndianaData\\simulation\\completeSimulationDataExSmallSNRNotAlign.RData")
# # indianaSimulationData = completeSimulationDataExSmallSNR[,c(5,(ncol(completeSimulationDataExSmallSNR)-T+2):ncol(completeSimulationDataExSmallSNR))]
# 
# # load("R:\\THESIS\\IndianaData\\simulation\\completeSimulationDataSmallSNRNotAlign.RData")
# # indianaSimulationData = completeSimulationDataSmallSNR[,c(5,(ncol(completeSimulationDataSmallSNR)-T+2):ncol(completeSimulationDataSmallSNR))]
# 
# 
# # load("R:\\THESIS\\IndianaData\\simulation\\completeSimulationDataMediumSNRNotAlign.RData")
# # indianaSimulationData = completeSimulationDataMediumSNR[,c(5,(ncol(completeSimulationDataMediumSNR)-T+2):ncol(completeSimulationDataMediumSNR))]
# 
# 
# simu = indianaSimulationData[order(indianaSimulationData$NAME10),]
# 
# data = t(simu[,-1])
# colnames(data) = simu[,1]
# 
# # centeredData = data*10^5
# centeredData = matrix(0,nrow = nrow(data), ncol = ncol(data))
# for (i in 1:ncol(data)){
#   centeredData[,i] = (data[,i] - mean(data[,i]))/sd(data[,i])
# }
# 
# shift = 0
# centeredData = centeredData + shift
# 
# timeUnit = 7
# numComp = 5 # number of selected components
# lagNumClip = 0
# N = ncol(centeredData)
# 
# numClip = floor(nrow(centeredData)/timeUnit) - lagNumClip
# dataByUnit = array(0, dim = c(timeUnit, ncol(centeredData),numClip), dimnames = list(1:timeUnit, 1:ncol(centeredData),1:numClip))
# for (i in 1:numClip){
#   dataByUnit[,,i] = centeredData[((i-1)*timeUnit+1):(i*timeUnit),]
# }
# 
# eigenvector = array(0, dim = c(N, numComp, numClip), dimnames = list(1:N, 1:numComp, 1:numClip))
# 
# for (i in 1:(numClip)){
#   eigenvector[,,i] = prcomp(as.matrix(dataByUnit[,,i]))$rotation[,1:numComp]
# }
# 
# compH = matrix(0,nrow = numComp, ncol = numClip)
# for (i in 1:nrow(compH)){
#   for (j in 1:ncol(compH)){
#     compH[i,j] = eigenvector[,i,j] %*% colMeans(dataByUnit[,,j])
#   }
# }
# 
# 
# initial_eigenvalue = (as.vector(prcomp(as.matrix(dataByUnit[,,1])))$sdev)^2
# a_0 = initial_eigenvalue[1:numComp]
# 
# sum(initial_eigenvalue[1:numComp])/sum(initial_eigenvalue)
# 
# #install.packages("Rcpp")
# library(Rcpp)
# 
# # install.packages("RcppArmadillo")
# # library(RcppArmadillo)
# 
# sourceCpp("R:/FutureStudy/gibbs_pca4.cpp")
# 
# dd = matrix(0, nrow = numClip, ncol = ncol(centeredData))
# for (i in 1:ncol(centeredData)){
#   dd[,i] = colMeans(dataByUnit[,i,])
# }
# train_data = dd
# D = nrow(train_data)  # days
# N = ncol(train_data)  # counties
# M = 10000 # 1086s
# K = numComp
# 
# adjH = rep(0,N*K*D)
# for (j in 0:(N-1)){
#   for (k in 0:(K-1)){
#     for (t in 0:(D-1)){
#       adjH[K*D*j+k*D+t + 1] = eigenvector[j+1,k+1,t+1]
#     }
#   }
# }
# 
# ptm = proc.time()
# set.seed(1234)
# re = gibbs_pca4(M, K, train_data, adjH, a_0)  # 593s M = 10000  K = 4
# proc.time() - ptm
# 
# 
# aaa = re$a
# # aaa = array(0, dim = c(M,K,D))
# # for (m in 1:M){
# #   for (k in 1:K){
# #     aaa[m,k,] = aa[(K*D*(m-1)+D*(k-1)+1):(K*D*(m-1)+D*k)]
# #   }
# # }
# 
# 
# sigma2_mean = colMeans(re$sigma2[5001:10000,])
# # plot(sigma2_mean,type = "b",xlab = "days", ylab = "", main = "sigma^2 (ReGPCA)")
# # for (k in 1:K){
# #   plot(aaa[k,],type = "b",xlab = "days",ylab = "", main = paste("a_",k,"(ReGPCA)"))
# #   
# # }
# 
# xxi = re$xi
# xi_final = array(0, dim = c(M,K,D))
# for (m in 0:(M-1)){
#   for (k in 0:(K-1)){
#     xi_final[m+1,k+1,] = xxi[(K*D*(m)+D*k+1):(K*D*(m)+D*(k+1))]
#   }
# }
# 
# # xi_final_1_mean = colMeans(xi_final[5001:10000,1,])
# # plot(xi_final_1_mean,type = "b",xlab = "days",ylab = "", main = "xi_1 (ReGPCA)")
# # xi_final_2_mean = colMeans(xi_final[5001:10000,2,])
# # plot(xi_final_2_mean,type = "b",xlab = "days",ylab = "", main = "xi_2 (ReGPCA)")
# # xi_final_3_mean = colMeans(xi_final[5001:10000,3,])
# # plot(xi_final_3_mean,type = "b",xlab = "days",ylab = "", main = "xi_3 (ReGPCA)")
# # xi_final_4_mean = colMeans(xi_final[5001:10000,4,])
# # plot(xi_final_4_mean,type = "b",xlab = "days",ylab = "", main = "xi_4 (ReGPCA)")
# 
# 
# 
# # for( i in seq(10,80,10)){
# #   ts.plot(re$sigma2[5001:10000,i],main=paste("County",i),ylab=paste("sigma^2"))
# # }
# # 
# # 
# # ts.plot(xi_final[5001:10000,1,1],ylab = paste('xi',11))
# # ts.plot(xi_final[5001:10000,2,1],ylab = paste('xi',21))
# # ts.plot(xi_final[5001:10000,3,1],ylab = paste('xi',31))
# # ts.plot(xi_final[5001:10000,4,1],ylab = paste('xi',41))
# # 
# # ts.plot(aaa[5001:10000,1,1],ylab = paste('a',11))
# # ts.plot(aaa[5001:10000,2,1],ylab = paste('a',21))
# # ts.plot(aaa[5001:10000,3,1],ylab = paste('a',31))
# # ts.plot(aaa[5001:10000,4,1],ylab = paste('a',41))
# 
# 
# 
# data_topred = centeredData[(dim(dataByUnit)[1]*dim(dataByUnit)[3]+1):nrow(centeredData),]
# theta = rep(0,ncol(centeredData))
# for (i in 1:length(theta)){
#   theta[i] = sum(aaa[,numClip]*eigenvector[i,,numClip]*(1+colMeans(xi_final[5001:10000,,numClip])))
# }
# 
# data_pred = matrix(0, nrow(data_topred),ncol(data_topred))
# for (i in 1:ncol(data_pred)){
#   data_pred[,i] = rnorm(nrow(data_topred),theta[i], sqrt(sigma2_mean[i]))
# }
# 
# J = nrow(data_pred)
# scaler = 1
# MSPE_ReGPCA = NULL
# for (j in 1:J){
#   MSPE_ReGPCA[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
# }
# 
# # plot(MSPE_ReGPCA,type = "b",xlab = "days",ylab = "",main =paste("MSPE (ReGPCA, a_eigen), shift = ",shift))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## filter
# ##########################################################
# dataObserv = t(aaa)
# write.table(dataObserv, row.names = FALSE, col.names = FALSE, sep = ',', file = paste('R:\\FutureStudy\\matlab\\ReGPCA_shift',timeUnit,'C',numComp,'.csv',sep = ''))
#           
# # x1=5.7763;    x2=0.2826;   x3=-2.7445;   x4=-0.3047;  x5=1; x6=0; x7=1  # a_eigen numcomp = 5 Medium
# # x1=5.2496;    x2=0.2559;   x3=-2.6856;   x4=-0.3069;  x5=1; x6=0; x7=1  # a_eigen numcomp = 5 Small
# x1=5.2240;    x2=0.1751;   x3=-2.8729;   x4=-0.3026;  x5=1; x6=0; x7=1  # a_eigen numcomp = 5 ExSmall
# 
# 
# diffusion = x7*diag(1,ncol(dataObserv))
# # specify the function in state equation
# stateFunction = function(x){ # x: state
#   # y = diffusion%*%x + x*(1-x/x8)
#   y = diffusion%*%x
#   return(y)
# }
# 
# mu = function(x, y){ # x: state, y: observation
#   z = solve(sigmaInverse) %*% (solve(stateNoiseVariance) %*% stateFunction(x) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% y)
#   return(z)
# }
# 
# phi = function(x, y){ # x: state, y: observation
#   z = 1/2 * t(y - observe2StateMatrix %*% stateFunction(x)) %*% solve(K) %*% (y - observe2StateMatrix %*% stateFunction(x))
#   return(z)
# }
# 
# ##########################################################
# data = dataObserv
# observeDimension = ncol(data)
# stateDimension = observeDimension
# observe2StateMatrix = diag(x5,observeDimension)
# 
# diag(observe2StateMatrix[1:(observeDimension-1),2:observeDimension]) = x6
# diag(observe2StateMatrix[2:observeDimension,1:(observeDimension-1)]) = x6
# stateNoiseMean = rep(0, stateDimension)
# stateNoiseVariance = diag(1/nrow(data), stateDimension, stateDimension) # G
# CQ = matrix(x2, stateDimension, stateDimension)
# diag(CQ) = x1
# 
# stateNoiseVariance = t(CQ) %*% CQ
# # stateNoiseVariance = CQ
# 
# observeNoiseMean = rep(0, observeDimension)
# CR = matrix(x4, observeDimension,observeDimension)
# diag(CR) = x3
# observeNoiseVariance = t(CR) %*% CR
# # observeNoiseVariance = CR
# 
# sigmaInverse = solve(stateNoiseVariance) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% observe2StateMatrix
# K = observe2StateMatrix %*% stateNoiseVariance %*% t(observe2StateMatrix) + observeNoiseVariance
# 
# library("Matrix")
# L = chol(solve(sigmaInverse))
# 
# # specify number of state samples M (say draw M=1000 samples of x)
# 
# M=10
# 
# set.seed(1234)
# 
# # with resampling, better results
# #######################################################################################################################
# statePredict = array(0, dim=c(nrow(data)+1, M, ncol(data)))
# library("MASS")
# # statePredict[1,,] = mvrnorm(n = M, stateNoiseMean, stateNoiseVariance)
# statePredict[1,,] = matrix(0, M, ncol(data))
# statePredictResample = statePredict
# 
# weights = matrix(0, nrow(data)+1, M)
# weights[1,] = 1/M
# 
# tempNegativePhi = matrix(0, nrow(data)+1, M)
# tempNegativePhi[1,] = 1
# 
# # In case that -phi is too Small
# for (i in 2:(nrow(data)+1)){
#   for (j in 1:M){
#     statePredict[i,j,] = mu(statePredict[i-1,j,],data[i-1,]) + L %*% matrix(mvrnorm(n = 1, rep(0,stateDimension), diag(1,stateDimension, stateDimension)), stateDimension,1)
#     tempNegativePhi[i,j] = -phi(statePredict[i-1,j,], data[i-1,])
#   }
#   updatePhi = tempNegativePhi[i,] + tempNegativePhi[i-1,]
#   maxTemp = max(updatePhi)
#   weights[i,] = exp(updatePhi-maxTemp)/sum(exp(updatePhi-maxTemp))
#   resamples = sample(1:M, M, replace = TRUE, prob = weights[i,])
#   statePredictResample[i,,] = statePredict[i,resamples,]
# }
# 
# weightsRoundup = round(weights, digits = 3)
# 
# weightedStates = matrix(0,nrow(data),ncol(data))
# for (i in 1:nrow(data)){
#   for (j in 1:ncol(data)){
#     weightedStates[i,j] = weighted.mean(statePredictResample[i+1,,j], weights[i+1,])
#   }
# }
# 
# 
# aaa_f = t(weightedStates)
# 
# data_topred = centeredData[(dim(dataByUnit)[1]*dim(dataByUnit)[3]+1):nrow(centeredData),]
# theta = rep(0,ncol(centeredData))
# for (i in 1:length(theta)){
#   for (k in 1:nrow(aaa_f)){
#     theta[i] = theta[i] + aaa_f[k,numClip]*eigenvector[i,k,numClip]*(1+colMeans(xi_final[5001:10000,k,])[numClip])
#   }
# }
# 
# data_pred = matrix(0, nrow(data_topred),ncol(data_topred))
# for (i in 1:ncol(data_pred)){
#   data_pred[,i] = rnorm(nrow(data_topred),theta[i], sqrt(sigma2_mean[i]))
# }
# 
# J = nrow(data_pred)
# scaler = 1
# MSPE_ReGPCA_Filter = NULL
# for (j in 1:J){
#   MSPE_ReGPCA_Filter[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
# }
# 
# # plot(MSPE_ReGPCA_Filter,type = "b",xlab = "days",ylab = "",main = paste("MSPE (ReGPCA_Filter, a_eigen), shift = ",shift))
# 
# plot(MSPE_ReGPCA,type = "b",xlab = "days",ylab = "MSPE",main = "Simulation 3",ylim = c(min(MSPE_ReGPCA,MSPE_ReGPCA_Filter), max(MSPE_ReGPCA,MSPE_ReGPCA_Filter)))
# lines(MSPE_ReGPCA_Filter, type = "b", col = "red")
# legend("topright",c("ReGP","ReGPF"),horiz = F,col = c(1,2),lty = rep(1,2), pch=rep(1,2),lwd=rep(1,2),cex=1)
# 
# # plot.new()
# # par(xpd=TRUE)
# # legend("center",c("ReGP","ReGPF"),horiz = F,col = c(1,2),lty = rep(1,2), pch=rep(1,2),lwd=rep(2,2),cex=1.2)
# # par(xpd=FALSE)





## add randomness to a_0
T <- 1000

# load("R:\\THESIS\\IndianaData\\simulation\\completeSimulationDataExSmallSNRNotAlign.RData")
# indianaSimulationData = completeSimulationDataExSmallSNR[,c(5,(ncol(completeSimulationDataExSmallSNR)-T+2):ncol(completeSimulationDataExSmallSNR))]

# load("R:\\THESIS\\IndianaData\\simulation\\completeSimulationDataSmallSNRNotAlign.RData")
# indianaSimulationData = completeSimulationDataSmallSNR[,c(5,(ncol(completeSimulationDataSmallSNR)-T+2):ncol(completeSimulationDataSmallSNR))]


load("R:\\THESIS\\IndianaData\\simulation\\completeSimulationDataMediumSNRNotAlign.RData")
indianaSimulationData = completeSimulationDataMediumSNR[,c(5,(ncol(completeSimulationDataMediumSNR)-T+2):ncol(completeSimulationDataMediumSNR))]


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
# a_0 = initial_eigenvalue[1:numComp]
set.seed(1234)
a_0 = initial_eigenvalue[1:numComp] + rnorm(5,0,1)


sum(initial_eigenvalue[1:numComp])/sum(initial_eigenvalue)

#install.packages("Rcpp")
library(Rcpp)

# install.packages("RcppArmadillo")
# library(RcppArmadillo)


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

sourceCpp("R:/FutureStudy/gibbs_pca4.cpp")
# sourceCpp("R:/FutureStudy/gibbs_pca3.cpp")
ptm = proc.time()
set.seed(1234)
re = gibbs_pca4(M, K, train_data, adjH, a_0)  # 593s M = 10000  K = 4
sga = 1
# re = gibbs_pca3(M, K, train_data, adjH, sga)  # 593s M = 10000  K = 4
proc.time() - ptm


aaa = re$a
# aaa = array(0, dim = c(M,K,D))
# for (m in 1:M){
#   for (k in 1:K){
#     aaa[m,k,] = aa[(K*D*(m-1)+D*(k-1)+1):(K*D*(m-1)+D*k)]
#   }
# }


sigma2_mean = colMeans(re$sigma2[5001:10000,])
# plot(sigma2_mean,type = "b",xlab = "days", ylab = "", main = "sigma^2 (ReGPCA)")
# for (k in 1:K){
#   plot(aaa[k,],type = "b",xlab = "days",ylab = "", main = paste("a_",k,"(ReGPCA)"))
# 
# }

xxi = re$xi
xi_final = array(0, dim = c(M,K,D))
for (m in 0:(M-1)){
  for (k in 0:(K-1)){
    xi_final[m+1,k+1,] = xxi[(K*D*(m)+D*k+1):(K*D*(m)+D*(k+1))]
  }
}

# xi_final_1_mean = colMeans(xi_final[5001:10000,1,])
# plot(xi_final_1_mean,type = "b",xlab = "days",ylab = "", main = "xi_1 (ReGPCA)")
# xi_final_2_mean = colMeans(xi_final[5001:10000,2,])
# plot(xi_final_2_mean,type = "b",xlab = "days",ylab = "", main = "xi_2 (ReGPCA)")
# xi_final_3_mean = colMeans(xi_final[5001:10000,3,])
# plot(xi_final_3_mean,type = "b",xlab = "days",ylab = "", main = "xi_3 (ReGPCA)")
# xi_final_4_mean = colMeans(xi_final[5001:10000,4,])
# plot(xi_final_4_mean,type = "b",xlab = "days",ylab = "", main = "xi_4 (ReGPCA)")
# xi_final_5_mean = colMeans(xi_final[5001:10000,5,])
# plot(xi_final_5_mean,type = "b",xlab = "days",ylab = "", main = "xi_5 (ReGPCA)")


# countyName = colnames(pureCompleteInci)
# setwd("R:\\FutureStudy\\figure")
# for( i in seq(10,90,10)){
#   pdf(paste("tracePlotCounty_",i,".pdf",sep=''),width=6,height=3)
#   ts.plot(re$sigma2[5001:10000,i],main=paste(countyName[i]),ylab=paste("sigma^2"))
#   dev.off()
# }
# 
# 
# ts.plot(xi_final[5001:10000,1,1],ylab = paste('xi',11))
# ts.plot(xi_final[5001:10000,2,1],ylab = paste('xi',21))
# ts.plot(xi_final[5001:10000,3,1],ylab = paste('xi',31))
# ts.plot(xi_final[5001:10000,4,1],ylab = paste('xi',41))
# 
# ts.plot(aaa[5001:10000,1,1],ylab = paste('a',11))
# ts.plot(aaa[5001:10000,2,1],ylab = paste('a',21))
# ts.plot(aaa[5001:10000,3,1],ylab = paste('a',31))
# ts.plot(aaa[5001:10000,4,1],ylab = paste('a',41))



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

# plot(MSPE_ReGPCA,type = "b",xlab = "days",ylab = "",main =paste("MSPE (ReGPCA, a_eigen), shift = ",shift))
















## filter
##########################################################
dataObserv = t(aaa)
# write.table(dataObserv, row.names = FALSE, col.names = FALSE, sep = ',', file = paste('R:\\FutureStudy\\matlab\\ReGPCA_shift',timeUnit,'C',numComp,'.csv',sep = ''))

          
          
# x1=5.2170;    x2=0.1733;   x3=-2.8533;   x4=-0.3013;  x5=1; x6=0; x7=1  # numcomp = 5  ExSmall
# x1=5.2746;    x2=0.2743;   x3=-2.7089;   x4=-0.3253;  x5=1; x6=0; x7=1  # Small
x1=5.3017;    x2=0.2863;   x3=-2.6125;   x4=-0.2548;  x5=1; x6=0; x7=1  # eigen numcomp = 5  Medium

          

diffusion = x7*diag(1,ncol(dataObserv))
# specify the function in state equation
stateFunction = function(x){ # x: state
  # y = diffusion%*%x + x*(1-x/x8)
  y = diffusion%*%x
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
    theta[i] = theta[i] + aaa_f[k,numClip]*eigenvector[i,k,numClip]*(1+colMeans(xi_final[5001:10000,k,])[numClip])
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

# plot(MSPE_ReGPCA_Filter,type = "b",xlab = "days",ylab = "",main = paste("MSPE (ReGPCA_Filter, a_eigen), shift = ",shift))

# plot(MSPE_ReGPCA,type = "b",xlab = "days",ylab = "MSPE",main = "Simulation 3",ylim = c(min(MSPE_ReGPCA,MSPE_ReGPCA_Filter), max(MSPE_ReGPCA,MSPE_ReGPCA_Filter)))
# lines(MSPE_ReGPCA_Filter, type = "b", col = "red")
# legend("topright",c("ReGP","ReGPF"),horiz = F,col = c(1,2),lty = rep(1,2), pch=rep(1,2),lwd=rep(1,2),cex=1)

# plot.new()
# par(xpd=TRUE)
# legend("center",c("ReGP","ReGPF"),horiz = F,col = c(1,2),lty = rep(1,2), pch=rep(1,2),lwd=rep(2,2),cex=1.2)
# par(xpd=FALSE)






dataObserv = t(aaa)
library(dlm)
data = dataObserv
ncomp = 5
VV = rep(0,ncomp)
WW = rep(0,ncomp)
for (icol in 1:ncomp){
  outGibbs <- dlmGibbsDIG(data[,icol], dlmModPoly(1, dV = var(data[,icol]), dW = var(data[,icol])),
                          a.y = var(data[,icol]), b.y = 1, a.theta = var(data[,icol]), b.theta = 1,
                          n.sample = 1000,
                          
                          save.states = FALSE)
  burn <- 500
  
  dV <- outGibbs$dV[-(1:burn)]
  dW <- outGibbs$dW[-(1:burn)]
  
  postMean = mcmcMean(with(outGibbs, sqrt(cbind(dV, dW))))
  VV[icol] = as.vector(postMean)[1]
  WW[icol] = as.vector(postMean)[3]
}



###################################################################################################
###################################################################################################
ptm = proc.time()
set.seed(1234)
library(dlm)
data = read.csv("R:/FutureStudy/matlab/ReGPCA_shift7C5.csv", header = FALSE)
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
  
  aa = matrix(outGibbs$theta,156,1000) 
  weightedStates[,icol] = rowMeans(aa)
}
proc.time() - ptm
###################################################################################################
###################################################################################################
x5=1; x6=0; x7=1  




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
CQ = matrix(0, stateDimension, stateDimension)
diag(CQ) = WW

stateNoiseVariance = t(CQ) %*% CQ
# stateNoiseVariance = CQ

observeNoiseMean = rep(0, observeDimension)
CR = matrix(0, observeDimension,observeDimension)
diag(CR) = VV
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

for (i in 1:ncomp){
  ts.plot(aaa[i,],ylab = paste('a',i))
  lines(aaa_f[i,],col = "red")
  legend("topright",c("a","a_filter"),horiz = F,col = c(1,2),lty = rep(1,2),lwd=rep(1,2),cex=1)
}

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
MSPE_ReGPCA_Filter_B = NULL
for (j in 1:J){
  MSPE_ReGPCA_Filter_B[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
}

# plot(MSPE_ReGPCA_Filter_B,type = "b",xlab = "days",ylab = "",main = paste("MSPE (ReGPCA_Filter, a_eigen), shift = ",shift))

plot(MSPE_ReGPCA,type = "b",xlab = "days",ylab = "MSPE",main = "Simulation 3",ylim = c(min(MSPE_ReGPCA,MSPE_ReGPCA_Filter,MSPE_ReGPCA_Filter_B), max(MSPE_ReGPCA,MSPE_ReGPCA_Filter,MSPE_ReGPCA_Filter_B)))
lines(MSPE_ReGPCA_Filter, type = "b", col = "red")
lines(MSPE_ReGPCA_Filter_B, type = "b", col = "blue")

# legend("topright",c("ReGP","ReGPF","ReGPFB"),horiz = F,col = c("black","red","blue"),lty = rep(1,3), pch=rep(1,3),lwd=rep(1,3),cex=1)









aaa_f = t(weightedStates)
for (i in 1:ncomp){
  ts.plot(aaa[i,],ylab = paste('a',i))
  lines(aaa_f[i,],col = "red")
  legend("topright",c("a","a_filter"),horiz = F,col = c(1,2),lty = rep(1,2),lwd=rep(1,2),cex=1)
}

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
MSPE_ReGPCA_Filter_B = NULL
for (j in 1:J){
  MSPE_ReGPCA_Filter_B[j] = sum((data_pred[1:j,]/scaler - data_topred[1:j,]/scaler)^2)/(j*ncol(data_topred))
}

# plot(MSPE_ReGPCA_Filter_B,type = "b",xlab = "days",ylab = "",main = paste("MSPE (ReGPCA_Filter, a_eigen), shift = ",shift))

plot(MSPE_ReGPCA,type = "b",xlab = "days",ylab = "MSPE",ylim = c(min(MSPE_ReGPCA,MSPE_ReGPCA_Filter,MSPE_ReGPCA_Filter_B), max(MSPE_ReGPCA,MSPE_ReGPCA_Filter,MSPE_ReGPCA_Filter_B)))
lines(MSPE_ReGPCA_Filter, type = "b", col = "red")
lines(MSPE_ReGPCA_Filter_B, type = "b", col = "blue")

