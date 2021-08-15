
mat.coefSPC = matrix(unlist(coefSPC.l), ncol = iters, byrow = F)
mat.coefCPC = matrix(unlist(coefCPC.l), ncol = iters, byrow = F)
mat.BiasSPC = matrix(unlist(BiasSPC.l), ncol = iters, byrow = F)
mat.BiasCPC = matrix(unlist(BiasCPC.l), ncol = iters, byrow = F)
mat.varSPC = matrix(unlist(varSPC.l), ncol = iters, byrow = F)
mat.varCPC = matrix(unlist(varCPC.l), ncol = iters, byrow = F)
require(matrixStats)
alpha = 0.05
my.coefSPC = data.frame('mean' = rowMeans(mat.coefSPC),'type' = 'PSC'
              ,'low' = rowMeans(mat.coefSPC) - qnorm(1-alpha/2)*sqrt(rowMeans(mat.varSPC)) #-rowMeans(mat.BiasSPC) 
              ,'high' = rowMeans(mat.coefSPC) + qnorm(1-alpha/2)*sqrt(rowMeans(mat.varSPC))
                 ,'PCs' = 1:Kmax)
my.coefCPC = data.frame('mean' = rowMeans(mat.coefCPC),'type' = 'CPC'
               ,'low' = rowMeans(mat.coefCPC) - qnorm(1-alpha/2)*sqrt(rowMeans(mat.varCPC))
               ,'high'= rowMeans(mat.coefCPC) + qnorm(1-alpha/2)*sqrt(rowMeans(mat.varCPC))
                        , 'PCs' = 1:Kmax)
A = rbind(my.coefSPC,my.coefCPC)
require(ggplot2)
mytitle = paste('Independent X, n=',n,', p=',p,', sparse=',spa)
ggplot(data=A, aes(x=PCs, y=mean, colour=type) )+ 
  geom_point() +   geom_line()+ 
  ggtitle( mytitle)+
  ylab(expression(beta[1] ))+ theme_bw()+
  geom_ribbon(aes(ymin=A$low, ymax=A$high), linetype=2, alpha=.1)

