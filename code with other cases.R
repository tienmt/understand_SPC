###==============================================
###==============================================
### 2) all variables are iid ===============####
###==============================================
p = 100
n = 1000

#### data generation  sparsity, #########
beta = rep(0,p)
beta[1] = 1
spa = p      # nonzero
beta[2:spa] = sample(c(-1,1), spa-1,replace = T)
sd = 1         # standard deviation of the noise
X = matrix(rnorm(n*p),nrow = n)
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
  X = matrix(rnorm(n*p),nrow = n)
  Y = X%*%beta + rnorm(n,0,sd)
  x.svd <- svd(X)        
  xm.svd = svd(X[,-1])   
  Z = Zm = cbind(X[,1], matrix(0, nr=n, nc=Kmax))
  Z[,2:(Kmax+1)] = x.svd$u[,1:Kmax]    # design for CPC
  Zm[,2:(Kmax+1)] = xm.svd$u[,1:Kmax]  # design for SPC
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







###==============================================
###==============================================
######   BINARY MATRIX ##############################
###==============================================
p = 100
n = 1000

#### data generation 1: sparsity, ideal for LASSO  #########
beta = rep(0,p)
beta[1] = 1
spa = 20      # nonzero
beta[2:spa] = sample(c(0,1), spa-1,replace = T)
sd = 1         # standard deviation of the noise

X = matrix(sample(c(-1,1), n*p,replace = T,prob = c(.5,.5)),nrow = n)
Y = X%*%beta + rnorm(n,0,sd)

###==============================================
####  simulate many times ###############
coefSPC.l = coefCPC.l = BiasSPC.l =
  BiasCPC.l = varSPC.l = varCPC.l = list()
Kmax = 30
iters = 100        # number of repeated run
pb <- txtProgressBar(min = 0, max = iters, style = 3)
for (ss in 1:iters) {
  X = matrix(sample(c(-1,1), n*p,replace =T,prob = c(.5,.5)),nrow = n)
  Y = X%*%beta + rnorm(n,0,sd)
  x.svd <- svd(X)        
  xm.svd = svd(X[,-1])   
  Z = Zm = cbind(X[,1], matrix(0, nr=n, nc=Kmax))
  Z[,2:(Kmax+1)] = x.svd$u[,1:Kmax]    # design for CPC
  Zm[,2:(Kmax+1)] = xm.svd$u[,1:Kmax]  # design for SPC
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


# try with lasso
require(glmnet)
lasso.fit <- cv.glmnet(X,Y,intercept=F)
head(coef(lasso.fit )[-1] )
beta[1:6]








###==============================================
###==============================================
### all variables are correlated =========== ####
####=============================================
p = 100
n = 1000

# generate covariance matrix  #
S = matrix(0, nrow = p, ncol = p)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    S[i, j] = 0.5^abs(i - j)
  }
}
S = S + t(S)
diag(S) = 1
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)

#### data generation 1: sparsity, ideal for LASSO  #########
beta = rep(0,p)
beta[1] = 1
spa = p     # nonzero
beta[2:spa] = sample(c(-1,1), spa-1,replace = T)

sd = 1         # standard deviation of the noise
X = matrix( rnorm(n*p), nr=n) %*% S.sqrt  
Y = X%*%beta + rnorm(n,0,sd)

#########################################3
####  simulate many times ###############
coefSPC.l = coefCPC.l = BiasSPC.l =
  BiasCPC.l = varSPC.l = varCPC.l = list()
iters = 100        # number of repeated run
Kmax = 30
pb <- txtProgressBar(min = 0, max = iters, style = 3)
for (ss in 1:iters) {
  X = matrix(rnorm(n*p),nrow = n)%*% S.sqrt
  Y = X%*%beta + rnorm(n,0,sd)
  x.svd <- svd(X)        
  xm.svd = svd(X[,-1])   
  Z = Zm = cbind(X[,1], matrix(0, nr=n, nc=Kmax))
  Z[,2:(Kmax+1)] = x.svd$u[,1:Kmax]    # design for CPC
  Zm[,2:(Kmax+1)] = xm.svd$u[,1:Kmax]  # design for SPC
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



# try with lasso
require(glmnet)
lasso.fit <- cv.glmnet(X,Y,intercept=F)
head(coef(lasso.fit )[-1] )
beta[1:6]
