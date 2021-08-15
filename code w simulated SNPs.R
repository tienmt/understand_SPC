#####################################################
###### SIMULATED SNPs MATRIX ##############################
snps <- simGWAS$snps
phen <- factor(simGWAS$phen)

# try with lasso
require(glmnet)
LASSO <- cv.glmnet(snps,phen,family="binomial")
res <- which(as.vector(t(coef(LASSO, s="lambda.min")))[-1] !=0)
pval <- apply(snps, 2, function(e)
  fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
pval.corrected.bonf <- p.adjust(pval, method="bonferroni")
res <- which(pval.corrected.bonf < 0.05)

#########################################3
####  simulate many times ###############
coefSPC.l = coefCPC.l = list()
x.svd <- svd(snps);  # SPC method
pval <-apply(snps, 2, function(x) coef(summary(glm(phen~x + x.svd$u[,1:10], family = binomial )))[2,4]) #
pval.corrected.bonf <- p.adjust(pval, method="bonferroni")
res <- which(pval.corrected.bonf < 0.05)



Kmax = 100    # maximal number of PCs used in procedure
coefCPC = c()
pb <- txtProgressBar(min = 0, max = p, style = 3)
for (jj in 1:p)  {
  xm.svd = svd(snps[,-jj]);   # CPC method
  coefCPC[jj] = coef(summary(lm(Y~ snps[,jj] + xm.svd$u[,1:Kmax] )))[2,4]
  setTxtProgressBar(pb, jj)
}
names(coefCPC) = colnames(X)
sele.CPC = which(coefCPC < 1e-4)

