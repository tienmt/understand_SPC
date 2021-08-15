require(BGLR)
require(irlba)
data(wheat)
Y <- wheat.Y[,1]
X <- wheat.X; rm(wheat.A,wheat.sets,wheat.X)
n = nrow(X)
p = ncol(X)


Kmax = 10    # maximal number of PCs used in procedure
x.svd <- irlba(X,Kmax)  # SPC method
coefCPC = coefSPC = c()
pb <- txtProgressBar(min = 0, max = p, style = 3)
for (jj in 1:p)  {
  xm.svd = irlba(X[,-jj],Kmax)    # CPC method
  coefCPC[jj] = coef(summary(lm(Y~ X[,jj] + xm.svd$u )))[2,1]
  coefSPC[jj] = coef(summary(lm(Y~ X[,jj] + x.svd$u )))[2,1]
  setTxtProgressBar(pb, jj)
}
par(mfrow=c(1,2))
hist(abs(coefCPC - coefSPC),breaks = 20,main ='',xlab = 'absolute error')
hist(abs(coefCPC - coefSPC)/abs(coefSPC),
     breaks = 30,main ='',xlab = 'relative error')
axis(1,at=seq(0,39,0.5),labels =seq(0,39,0.5) )

tam = abs(coefCPC - coefSPC)/abs(coefSPC)
tam.indx =  which(tam>0.5)



Kmax = 30    # maximal number of PCs used in procedure
x.svd <- irlba(X,Kmax)  # SPC method
coefCPC.602 = coefSPC.602 = c()
xm.svd = irlba(X[,-602],Kmax)    # CPC method
for (k in 1:Kmax)  {
  coefCPC.602[k] = coef(summary(lm(Y~ X[,602] + xm.svd$u[,1:k] )))[2,1]
  coefSPC.602[k] = coef(summary(lm(Y~ X[,602] + x.svd$u[,1:k] )))[2,1]
}
plot(coefCPC.602,type='b',col='blue')
lines(coefSPC.602,type='b',col='red')
grid()



sele.CPC = which(coefCPC < 1e-4)
pval.corrected.bonf <- p.adjust(coefCPC, method="bonferroni")
sele10 <- which(pval < 0.005)
pval <- apply(X, 2, function(x) coef(summary(lm(Y~x + x.svd$u[,1:30])))[2,4])
plot(-log10(pval.corrected.bonf)); abline(h=3,col='red')
pval.corrected.bonf <- p.adjust(pval, method="bonferroni")
sele10 <- which(pval < 0.005)

sele10 = which(apply(X, 2, function(x) coef(summary(lm(Y~x + x.svd$u[,1:10] )))[2,4]) < 1e-6)
sele20 = which(apply(X, 2, function(x) coef(summary(lm(Y~x + x.svd$u[,1:20] )))[2,4]) < 1e-5)
sele30 = which(apply(X, 2, function(x) coef(summary(lm(Y~x + x.svd$u[,1:30] )))[2,4]) < 1e-4)
sele50 = which(apply(X, 2, function(x) coef(summary(lm(Y~x + x.svd$u[,1:50] )))[2,4]) < 1e-3)
#
# try with lasso
lasso.fit <- glmnet::cv.glmnet(X,Y)
sum(coef(lasso.fit ) !=0 )
plot( abs(coef(lasso.fit)[-1] ) ); abline(h=0.2,col='red')
tam = coef(lasso.fit)[-1] ; names(tam) = colnames(X)
sele.lasso = which(abs(tam ) > 0.2)



  
plot(-log10(coefCPC.l[25:35]), type='l',col='red',ylim = c(2.6,4.3),xaxt = "n")
lines(-log10(coefSPC.l[25:35]), type ='l',col='blue',xaxt = "n")
abline(h=4)
grid()
axis(1,at=1:11,labels = 25:35)#,xaxt = "n"

#== with adegenet ====
library(adegenet)
pca1 <- dudi.pca(df = X, scale = FALSE, scannf = FALSE, nf = 50)
head(pca1$eig)
s.label(pca1$li, sub="PCA - PC 1 and 2")
add.scatter.eig(pca1$eig,4,1,2, ratio=.3, posi="topleft")
D <- dist(pca1$li[,1:2])^2
clust <- hclust(D, method="complete")
plot(clust, main="Clustering (complete linkage) based on the first 2 PCs", cex=.4)

pop <- factor(cutree(clust, k=2))
s.class(pca1$li, fac=pop, col=transp(funky(5)), cpoint=2,
         sub="PCA - axes 1 and 2")
s.class(pca1$li, xax=3, yax=4, fac=pop, col=transp(funky(5)),
         cpoint=2, sub="PCA - axes 3 and 4")
dapc1 <- dapc(X, Y)
loadingplot(dapc1$var.contr, threshold = 0.0008)



tam = apply(X, 2, function(x) coef(summary(lm(Y~x )))[2,4])
pvaltam.crcted.bonf <- p.adjust(tam, method="bonferroni")
res <- which(pvaltam.crcted.bonf < 0.05)
sigSNPs <- list()
sigSNPs[[1]] <- list(res)
names(sigSNPs)[[1]] <- "univariate"
names(sigSNPs$univariate)[[1]] <- "bonferroni"
pval.corrected.fdr <- p.adjust(tam, method="fdr")
res <- which(pval.corrected.fdr < 0.05)
sigSNPs$univariate[[2]] <- res
names(sigSNPs$univariate)[[2]] <- "fdr"
LASSO <- cv.glmnet(X, Y)
beta <- as.vector(t(coef(LASSO, s="lambda.min")))
res <- which(beta[-1] !=0)
coefs.LASSO <- beta[-1][res]
names(coefs.LASSO) <- colnames(X)[res]
sigSNPs[[2]] <- list(res)
names(sigSNPs)[[2]] <- "multivariate"
names(sigSNPs$multivariate)[[1]] <- "lasso"
yBin<-ifelse(Y>0,1,0)
xval1 <- xvalDapc(X, yBin, n.rep=20) 
dapc1 <- xval1[[7]]
result <- snpzip(X, dapc1, method="ward", xval.plot = FALSE,
                 plot = TRUE, loading.plot = TRUE)
res <- result$FS[[2]]
sigSNPs$multivariate[[2]] <- res
names(sigSNPs$multivariate)[[2]] <- "dapc"

snps.corrected <- apply(X, 2, function(e)
  residuals(lm(e~pca1$li[,1]+pca1$li[,2]+pca1$li[,3]+pca1$li[,4]+pca1$li[,5]))) 
pval2 <- numeric(0)
for(i in 1:ncol(snps.corrected)){
  foo <- suppressWarnings(glm(Y ~ snps.corrected[,i]) )
  ANOVA <- anova(foo, test="Chisq")
  pval2[i] <- ANOVA$"Pr(>Chi)"[2]
} #
pval.corrected.bonf <- p.adjust(pval2, method="bonferroni")
res <- which(pval.corrected.bonf < 0.05)

sigSNPs.PCA <- list()
sigSNPs.PCA[[1]] <- list(res)
names(sigSNPs.PCA)[[1]] <- "univariate"
names(sigSNPs.PCA$univariate)[[1]] <- "bonferroni"
pval.corrected.fdr <- p.adjust(pval2, method="fdr")
res <- which(pval.corrected.fdr < 0.05)
sigSNPs.PCA$univariate[[2]] <- res
names(sigSNPs.PCA$univariate)[[2]] <- "fdr"
LASSO <- cv.glmnet(snps.corrected, Y)
beta <- as.vector(t(coef(LASSO, s="lambda.min")))
res <- which(beta[-1] !=0)
coefs.LASSO <- beta[-1][res]
names(coefs.LASSO) <- colnames(snps.corrected)[res]
sigSNPs.PCA[[2]] <- list(res)
names(sigSNPs.PCA)[[2]] <- "multivariate"
names(sigSNPs.PCA$multivariate) <- "lasso"
