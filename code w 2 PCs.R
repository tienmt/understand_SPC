###==============================================
###==============================================
### 1). TWO PCs, X1 IS UNCORRELATED W. X2 =======
###==============================================
sigma <- function(theta=0, lambda=c(1,1)) {
  cos.t <- cos(theta)
  sin.t <- sin(theta)
  a <- matrix(c(cos.t, sin.t, -sin.t, cos.t), ncol=2)
  t(a) %*% diag(lambda) %*% a
}
library(MASS)

p <- 1000    # dimensions
n1 <- 300   # First group population
n2 <- 300   # Second group population
n = n1+n2   # total sample
X <- rbind(mvrnorm(n1, c(-2,1), sigma(0, c(1/2,1)) ), 
           mvrnorm(n2, c(2,1), sigma(0, c(1,1/3)) ) ) 
eps <- 1  # Error SD should be small compared to the SDs for the blobs
X <- cbind(X, matrix(rnorm( n*(p-2),sd=eps), ncol=p-2))

fit <- prcomp(X)          # PCA
summary(fit)              # Brief summary showing two principal components
plot(fit$x[, 1:2], asp=1) # Display the first two components $
plot(X[, 1:2], asp=1,main = "data");    # Display the original data for comparison
points(X[1:n1,1:2], col="Blue") #...distinguish one of the blobs

#### data generation 1: sparsity, ideal for LASSO  #########
beta = rep(0,p)
beta[1] = 1
spa = 20      # nonzero
beta[2:spa] = sample(c(-1,1), spa-1,replace = T)
sd = 1         # standard deviation of the noise
Y = X%*%beta + rnorm(n,0,sd)

# try with lasso
require(glmnet)
lasso.fit <- cv.glmnet(X,Y,intercept=F)
head(coef(lasso.fit )[-1] )
beta[1:6]

###==============================================
####  simulate many times ###############
coefSPC.l = coefCPC.l = BiasSPC.l =
  BiasCPC.l = varSPC.l = varCPC.l = list()
Kmax = 30
iters = 100        # number of repeated run
pb <- txtProgressBar(min = 0, max = iters, style = 3)
for (ss in 1:iters) {
  X <- rbind(mvrnorm(n1, c(-2,1), sigma(0, c(1/2,1)) ), 
             mvrnorm(n2, c(2,1), sigma(0, c(1,1/3)) ) ) 
  X <- cbind(X, matrix(rnorm( n*(p-2),sd = eps), ncol= p-2))
  x.svd <- svd(X)
  xm.svd = svd(X[,-1])
  Z = Zm = cbind(X[,1], matrix(0, nr=n, nc=Kmax))
  Z[,2:(Kmax+1)] = x.svd$u[,1:Kmax]
  Zm[,2:(Kmax+1)] = xm.svd$u[,1:Kmax]
  gam.k =  diag(xm.svd$d) %*% t(xm.svd$v) %*% beta[-1]
  gam.bar = diag(x.svd$d) %*% t(x.svd$v) %*% beta
  coefSPC.l[[ss]] = coefCPC.l[[ss]] = BiasSPC.l[[ss]] =
    BiasCPC.l[[ss]] = varSPC.l[[ss]] = varCPC.l[[ss]] = rep(NA,Kmax)
  for (k in 1:Kmax){
    coefSPC.l[[ss]][k] = lm(Y~Z[,1:k]+0)$coefficients[1]
    coefCPC.l[[ss]][k] = lm(Y~Zm[,1:k]+0)$coefficients[1]
    BiasCPC.l[[ss]][k] = X[,1] %*% xm.svd$u[,k:min(p-1,n)] %*% gam.k[k:min(p-1,n)]/ 
      (X[,1] %*%X[,1] - sum(( X[,1] %*% xm.svd$u[,1:k])^2))
    varCPC.l[[ss]][k] = sd^2/ (X[,1] %*%X[,1] - sum(( X[,1] %*% xm.svd$u[,1:k])^2))
    BiasSPC.l[[ss]][k] = -beta[1] + X[,1] %*% x.svd$u[,k:min(p,n)] %*% gam.bar[k:min(p,n)]/ 
      (X[,1] %*%X[,1] - sum(( X[,1] %*% x.svd$u[,1:k])^2))
    varSPC.l[[ss]][k] = sd^2/ (X[,1] %*%X[,1] - sum(( X[,1] %*% x.svd$u[,1:k])^2))
  }
  setTxtProgressBar(pb, ss)
}







### calculating the average over multiple runs
coefSPC.l.mean = rowMeans(matrix(unlist(coefSPC.l), ncol = iters, byrow = F))
coefCPC.l.mean = rowMeans(matrix(unlist(coefCPC.l), ncol = iters, byrow = F))

plot(coefSPC.l.mean, type="l",lwd=3,col='blue',
     xlab = 'numbers of PCs k',ylab = expression(beta[1] ) )
grid()
lines(coefCPC.l.mean,col="red",lwd=3)
legend("topright", inset=c(0,0), cex = .5,
       c("SPC","CPC"),bg="grey96", lwd=c(3,3),col = c("blue","red") )


#require(matrixStats)
#coefSPC.l.sd = rowSds(matrix(unlist(coefSPC.l), ncol = iters, byrow = F))
#coefCPC.l.sd = rowSds(matrix(unlist(coefCPC.l), ncol = iters, byrow = F))





par(mfrow=c(1,4))         # Prepare to plot
plot(fit$x[, 1:2], asp=1) # Display the first two components $
plot(x[, 1:2], asp=1,main = "data");    # Display the original data for comparison
points(x[1:n1,1:2], col="Blue") #...distinguish one of the blobs
### SVD  ##########
x.svd <- svd(X); 
plot(x.svd$d, pch = '.',
     main = 'eigenvalues',ylab = 'X.svd$d', xlab = '')        # SPC method
points(x.svd$d[1:2], col="Blue",pch = 15)
xm.svd = svd(X[,-1])
plot(xm.svd$d, pch = '.',
     main = 'eigenvalues',ylab = 'X[,-1].svd$d', xlab = '')  # CPC method
points(xm.svd$d[1], col="Blue",pch = 15)

