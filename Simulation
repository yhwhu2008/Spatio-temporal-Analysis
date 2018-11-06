
# function for lambda
lambda.fun <- function(adj0=adj0, n.nodes=92, n.iter=500, T=22, rho1=0.5, rho2=0.5, sigma=0.5, level=5, c=1){

sims.xi <- matrix(0,nrow=n.nodes,ncol=T)
xi <- matrix(0,nrow=n.nodes,ncol=1)

#scale by row sums
rowtot <- rowSums(adj0)
adj <- diag(1/rowtot)%*%adj0

#simulate from time 2 to T
for (k in 2:T){
  #Gibbs sampler to update xi one at a time
  xi[,1]<-sims.xi[,k-1]
  for (j in 1:n.iter){
    for (i in 1:n.nodes){
      temp <- adj%*%xi
      xi[i,1] <- rnorm(n=1, mean=rho1*temp[i,]+rho2*sims.xi[i,k-1], sd=sigma/rowtot[i])
      }
  }
  sims.xi[,k]<-xi[,1]
}
return(sims.xi+log(level*c))
}













