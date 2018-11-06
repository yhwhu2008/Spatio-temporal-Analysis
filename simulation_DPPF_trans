# model 
#       y_it ~ Pois(exp(mu_i + a_1*cos(2*pi*t/C) + a_2*sin(2*pi*t/C)))
#       a_1 ~ U(0,1)
#       a_2 ~ U(0,1)
#       mu_i|G ~ G
#       G|alpha, G_0 ~ DP(alpha, G_0)
#       G_0 ~ N(0,1)
# posterior
#       y_t ~ \sum(j=1:J) w_j*Pois(exp(mu_j + a_1*cos(2*pi*t/C) + a_2*sin(2*pi*t/C)))
#       where w_j, mu_j, a_1, a_2 are all posterior estimated 

########################
########################
# monthN = readline("which month (2 digits): ")
########################
########################

monthN = '01'

set.seed(1234)

# input estimated parameters
a_1 = read.csv(paste('traceA1_k30_twoHour_auto_5e4_',monthN,'.csv',sep=''),header=F)
a_1 = a_1[nrow(a_1),]
a_2 = read.csv(paste('traceA2_k30_twoHour_auto_5e4_',monthN,'.csv',sep=''),header=F)
a_2 = a_2[nrow(a_2),]
weight = read.csv(paste('traceWEIGHT_k30_twoHour_auto_5e4_',monthN,'.csv',sep=''),header=F,sep = '')
weight = weight[nrow(weight),]
weight = weight[which(weight != 0)]
mu = read.csv(paste('traceMU_k30_twoHour_auto_5e4_',monthN,'.csv',sep=''),header=F, sep = '')
mu = mu[nrow(mu),]

yellowCabPro = read.csv(paste('yellowCabPartProcessData15M',monthN,'.csv',sep=''),header = T)
yellowCabPro = data.frame(yellowCabPro)
qtl = 0.9
realdata = yellowCabPro[,which(colSums(yellowCabPro) > quantile(colSums(yellowCabPro), probs = qtl))]

# genetate data for time to predicted

dd = c(31,28,31,30,31,30,31,31,30,31,30,31)

N = ncol(realdata)
days = dd[as.numeric(monthN)]
C = days*12 # Cycle
T = ceiling(C*0.8)
J = length(weight)# number of clusters

components = sample(1:J, prob = weight, size = N, replace = T)  
length(weight)
length(unique(components))

nsample = N
y_predict_1 = matrix(0,nrow = (C-T), ncol = N)

for (t in (T+1):C){
  y_predict_1[t-T,] = rpois(nsample, as.numeric(exp(mu[components]+ a_1*cos(2*pi*t/C) + a_2*sin(2*pi*t/C))))
}

y_predict_1 = data.frame(y_predict_1)
rownames(y_predict_1) = as.character(1:nrow(y_predict_1)) 
colnames(y_predict_1) = components
colnames(y_predict_1) <- paste(colnames(y_predict_1), 1:ncol(y_predict_1), sep = "_")
y_predict_1$ID <- rownames(y_predict_1)

library(reshape)
y_predict_1.m <- melt(y_predict_1, id = "ID")
y_predict_1.m <- cbind(y_predict_1.m, colsplit(y_predict_1.m$variable, split = "_", names = c("Measure", "N")))
y_predict_1.agg <- cast(y_predict_1.m, ID ~ Measure, fun = mean)

y_predict_1.agg$ID = as.numeric(y_predict_1.agg$ID)
y_predict_1.agg = y_predict_1.agg[order(y_predict_1.agg$ID),]
rownames(y_predict_1.agg) = y_predict_1.agg$ID

y_predict = y_predict_1[,1:(ncol(y_predict_1)-1)]


realdataToCompare_1 = realdata[(T+1):C,]

colnames(realdataToCompare_1) = components
colnames(realdataToCompare_1) <- paste(colnames(realdataToCompare_1), 1:ncol(realdataToCompare_1), sep = "_")
realdataToCompare_1$ID <- rownames(realdataToCompare_1)

realdataToCompare_1.m <- melt(realdataToCompare_1, id = "ID")
realdataToCompare_1.m <- cbind(realdataToCompare_1.m, colsplit(realdataToCompare_1.m$variable, split = "_", names = c("Measure", "N")))
realdataToCompare_1.agg <- cast(realdataToCompare_1.m, ID ~ Measure, fun = mean)

realdataToCompare = realdataToCompare_1[,1:(ncol(realdataToCompare_1)-1)]

y_predict_reorder = y_predict_1.agg[,-1]
y_predict_reorder <- y_predict_reorder[,order(colMeans(y_predict_reorder))]

realdataToCompare_reorder = realdataToCompare_1.agg[,-1]
realdataToCompare_reorder <- realdataToCompare_reorder[,order(colMeans(realdataToCompare_reorder))]

# use filter
set.seed(1234)

realdataTrain_1 = realdata[1:T,]

colnames(realdataTrain_1) = components
colnames(realdataTrain_1) <- paste(colnames(realdataTrain_1), 1:ncol(realdataTrain_1), sep = "_")
realdataTrain_1$ID <- rownames(realdataTrain_1)


realdataTrain_1.m <- melt(realdataTrain_1, id = "ID")
realdataTrain_1.m <- cbind(realdataTrain_1.m, colsplit(realdataTrain_1.m$variable, split = "_", names = c("Measure", "N")))
realdataTrain_1.agg <- cast(realdataTrain_1.m, ID ~ Measure, fun = mean)

realdataTrain = realdataTrain_1[,1:(ncol(realdataTrain_1)-1)]


y_predict_reorder = y_predict_1.agg[,-1]
y_predict_reorder <- y_predict_reorder[,order(colMeans(y_predict_reorder))]

realdataTrain_reorder = realdataTrain_1.agg[,-1]
realdataTrain_reorder <- realdataTrain_reorder[,order(colMeans(realdataTrain_reorder))]

# write.table(realdataTrain_reorder, file = paste('realdataTrain_reorder',monthN,'.csv',sep=''), row.names = F, col.names = F, sep=',')


# y_t = A*x_t + v_t
# x_t = D*x_{t-1} + w_t

## use matlab to find MLE
# simulated data with MLE
set.seed(1234)
##############################################################################
##incidence rate by DP regions

data1 = read.csv(paste('realdataTrain_reorder',monthN,'.csv',sep=''),header = F)
data = matrix(as.numeric(unlist(data1)),nrow(data1),ncol(data1))
observeDimension = ncol(data)
stateDimension = observeDimension
##############################################
x = matrix(0,nrow=12,ncol=6)
x[1,] = c(6.8216, 0.7260, -15.2793, -1.2095, 2.3566, -0.0266)     # M01
x[2,] = c(1.2373, 0.1582, -9.3074, -0.8456, 9.5240, -0.0214)      #M02
x[3,] = c(2.6044, 0.3816, -9.0397, 0.5596, 5.2881, -0.0193)      #M03
x[4,] = c(0.7740, 0.0968, -14.5960, -1.6566, 19.6617, -0.0191)    #M04
x[5,] = c(16.2064, 2.2956, -8.6438, -1.1795, 0.7367, -0.0155)   #M05
x[6,] = c(3.5991, 0.5038, -8.1626,  -1.0162, 2.9435, -0.0160)   #M06
x[7,] = c(73.4209, 12.5082, -5.0128,  2.5549, 0.1222, -0.0135)   #M07
x[8,] = c(9.8413, 1.4925, -7.5482, -0.9094, 1.0159, -0.0151)    #M08
x[9,] = c(0.4458, 0.0703, 4.7642, -2.5365, 20.6956, -0.0114)    #M09
x[10,] = c(11.9513, 1.6820, -8.1110,  -1.1565, 0.9186,  -0.0162)  #10
x[11,] = c(8.3549,  1.1731, -12.0603, -1.5235,  1.5141, -0.0155)  #M11
x[12,] = c(0.6088,  0.0865, -12.7479, -1.3253,  23.4970,  -0.0215)  #M12
D = diag(x[as.numeric(monthN),6],stateDimension)  # state to state
A = diag(x[as.numeric(monthN),5],stateDimension)
observe2StateMatrix = A

stateNoiseMean = rep(0, stateDimension)
stateNoiseVariance = diag(1/nrow(data), stateDimension, stateDimension) # G
CQ = matrix(x[as.numeric(monthN),2], stateDimension, stateDimension)
diag(CQ) = x[as.numeric(monthN),1]

stateNoiseVariance = t(CQ) %*% CQ

observeNoiseMean = rep(0, observeDimension)
CR = matrix(x[as.numeric(monthN),4], observeDimension,observeDimension)
diag(CR) = x[as.numeric(monthN),3]
observeNoiseVariance = t(CR) %*% CR

# specify the function in state equation
stateFunction = function(s){ # s: state
  y = D%*%s
  return(y)
}

mu = function(s, y){ # s: state, y: observation
  z = solve(sigmaInverse) %*% (solve(stateNoiseVariance) %*% stateFunction(s) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% y)
  return(z)
}

phi = function(s, y){ # s: state, y: observation
  z = 1/2 * t(y - observe2StateMatrix %*% stateFunction(s)) %*% solve(K) %*% (y - observe2StateMatrix %*% stateFunction(s))
  return(z)
}

##########################################################

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

set.seed(1234)
MM = 3
y_filterPredict = matrix(0,nrow = nrow(y_predict_reorder), ncol = ncol(y_predict_reorder))
x_filterPredict = y_filterPredict
x_filterPredict[1,] = weightedStates[nrow(weightedStates),]
# y_filterPredict[1,] = A %*% x_filterPredict[1,] + colMeans(mvrnorm(n=MM,observeNoiseMean,observeNoiseVariance))
y_filterPredict[1,] = A %*% x_filterPredict[1,]


i=2
while (i <= nrow(y_filterPredict)){
  # x_filterPredict[i,] = D %*% x_filterPredict[i-1,] + mvrnorm(n=1, stateNoiseMean, stateNoiseVariance)
  x_filterPredict[i,] = D %*% x_filterPredict[i-1,] + colMeans(mvrnorm(n=MM, stateNoiseMean, stateNoiseVariance))
  y_filterPredict[i,] = A %*% x_filterPredict[i,]
  if (all(x_filterPredict[i,] >= 0)){
    i = i+1
  }
}

y_filterPredict = data.frame(y_filterPredict)


# scaler = 100
scaler = 120  # counts/scaler = taxi call per minites
J = nrow(realdataToCompare_reorder)
MSPE1 = NULL
MSPE2 = NULL
for (j in 1:J){
  MSPE1[j] = sum(y_predict_reorder[1:j,]/scaler - realdataToCompare_reorder[1:j,]/scaler)^2/(j*ncol(realdataToCompare_reorder))
  MSPE2[j] = sum(y_filterPredict[1:j,]/scaler - realdataToCompare_reorder[1:j,]/scaler)^2/(j*ncol(realdataToCompare_reorder))
}


MSPE1 = NULL
MSPE2 = NULL
for (j in 1:J){
  MSPE1[j] = sum((y_predict_reorder[1:j,]/scaler - realdataToCompare_reorder[1:j,]/scaler)^2)/(j*ncol(realdataToCompare_reorder))
  MSPE2[j] = sum((y_filterPredict[1:j,]/scaler - realdataToCompare_reorder[1:j,]/scaler)^2)/(j*ncol(realdataToCompare_reorder))
}

plot(2*seq(1,length(MSPE1)),MSPE1,col=1,type='l', ylab = '', xlab = '',cex.axis=1.5)
lines(2*seq(1,length(MSPE2)),MSPE2,col=2,type='l')

########################################################
##########################################################

#  nonlinear model
set.seed(1234)
##############################################################################
##incidence rate by DP regions

data1 = read.csv(paste('realdataTrain_reorder',monthN,'.csv',sep=''),header = F)
data = matrix(as.numeric(unlist(data1)),nrow(data1),ncol(data1))
observeDimension = ncol(data)
stateDimension = observeDimension
##############################################
x = matrix(0,nrow=12,ncol=7)
x[1,] = c(12,           1.1,       -15,      -1.4,       1.3,       0.1,     1285.4)   #M01
x[2,] = c( 5.3,         0.9,      -7.8,       0.5,       2.3,     0.017,     1315.4)   #M02
x[3,] = c(49.7,           5,     -10.5,        -1,       0.3,      0.04,      10374)   #M03
x[4,] = c( 4.2179,    0.4779,  -14.3638,   -1.7879,    3.5936,    0.0496,  584.9565)   #M04
x[5,] = c( 2.1720 ,   0.4060,   -7.0468,    0.2756,    5.8368,    0.0221,  588.8904)   #M05
x[6,] = c(0.1754,    0.0231,   -8.1437,   -1.1010,   59.8272,    0.0480,    28.7836)   #M06
x[7,] = c(3.5005,    0.5783,   -8.8190,   -1.3290,    2.5559,    0.0307,   814.7901)   #M07
x[8,] = c(7.7,          1.1,      -7.6,     -0.97,       1.3,     0.016,   3049.8)    #M08
x[9,] = c(0.8003,    0.1213,   -8.5679,   -1.3639,   11.4651,    0.0343,  210.2451)    #M09
x[10,] = c(2.0164,   0.3925,   -6.3831,    0.4578,    5.8362,    0.0255,  441.3788)  #10
x[11,] = c(14,          1.9,     -11.8,      -1.3,       0.9,      0.02,  4670)  #M11
x[12,] = c(28.8,        3.8,     -12.7,      -1.5,       0.5,      0.04,    4240)  #M12
D = diag(x[as.numeric(monthN),6],stateDimension)  # state to state
A = diag(x[as.numeric(monthN),5],stateDimension)
observe2StateMatrix = A

stateNoiseMean = rep(0, stateDimension)
stateNoiseVariance = diag(1/nrow(data), stateDimension, stateDimension) # G
CQ = matrix(x[as.numeric(monthN),2], stateDimension, stateDimension)
diag(CQ) = x[as.numeric(monthN),1]

stateNoiseVariance = t(CQ) %*% CQ

observeNoiseMean = rep(0, observeDimension)
CR = matrix(x[as.numeric(monthN),4], observeDimension,observeDimension)
diag(CR) = x[as.numeric(monthN),3]
observeNoiseVariance = t(CR) %*% CR

# specify the function in state equation
stateFunction = function(s){ # s: state
  y = D%*%s + s*(1-s/x[as.numeric(monthN),7])
  return(y)
}

mu = function(s, y){ # s: state, y: observation
  z = solve(sigmaInverse) %*% (solve(stateNoiseVariance) %*% stateFunction(s) + t(observe2StateMatrix) %*% solve(observeNoiseVariance) %*% y)
  return(z)
}

phi = function(s, y){ # s: state, y: observation
  z = 1/2 * t(y - observe2StateMatrix %*% stateFunction(s)) %*% solve(K) %*% (y - observe2StateMatrix %*% stateFunction(s))
  return(z)
}

##########################################################

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


set.seed(1234)
MM = 3
y_filterPredict = matrix(0,nrow = nrow(y_predict_reorder), ncol = ncol(y_predict_reorder))
x_filterPredict = y_filterPredict
x_filterPredict[1,] = weightedStates[nrow(weightedStates),]
# y_filterPredict[1,] = A %*% x_filterPredict[1,] + colMeans(mvrnorm(n=MM,observeNoiseMean,observeNoiseVariance))
y_filterPredict[1,] = A %*% x_filterPredict[1,]


i=2
while (i <= nrow(y_filterPredict)){
  # x_filterPredict[i,] = D %*% x_filterPredict[i-1,] + mvrnorm(n=1, stateNoiseMean, stateNoiseVariance)
  x_filterPredict[i,] = D %*% x_filterPredict[i-1,] + colMeans(mvrnorm(n=MM, stateNoiseMean, stateNoiseVariance))
  y_filterPredict[i,] = A %*% x_filterPredict[i,]
  if (all(x_filterPredict[i,] >= 0)){
    i = i+1
  }
}

y_filterPredict = data.frame(y_filterPredict)


scaler = 120  # counts/scaler = taxi call per minites
J = nrow(realdataToCompare_reorder)
MSPE3 = NULL
for (j in 1:J){
  MSPE3[j] = sum(y_filterPredict[1:j,]/scaler - realdataToCompare_reorder[1:j,]/scaler)^2/(j*ncol(realdataToCompare_reorder))
}


save(MSPE1, file = paste('MSPE_SB',monthN,'.csv',sep=''))
save(MSPE2, file = paste('MSPE_linearDPPF',monthN,'.csv',sep=''))
save(MSPE3, file = paste('MSPE_nonlinearDPPF',monthN,'.csv',sep=''))
