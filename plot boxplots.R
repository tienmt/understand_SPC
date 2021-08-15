cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mat.coefSPC = matrix(unlist(coefSPC.l), ncol = iters, byrow = F)
mat.coefCPC = matrix(unlist(coefCPC.l), ncol = iters, byrow = F)
require(matrixStats)
#== summary data =====
my.coefSPC = data.frame('mean' = rowMeans(mat.coefSPC)
                        ,'low' = rowQuantiles(mat.coefSPC,probs = 0.025)
                        ,'high' = rowQuantiles(mat.coefSPC,probs = 0.975)
                        ,'pcs' = 1:Kmax ,'type' = 'PSC')
my.coefCPC = data.frame('mean' = rowMeans(mat.coefCPC)
                        ,'low' = rowQuantiles(mat.coefCPC,probs = 0.025)
                        ,'high'= rowQuantiles(mat.coefCPC,probs = 0.975)
                        , 'pcs' = 1:Kmax ,'type' = 'CPC')
A = rbind(my.coefSPC,my.coefCPC)
require(ggplot2)
#== bias n variances data =====
mat.BiasSPC = matrix(unlist(BiasSPC.l), ncol = iters, byrow = F)
mat.BiasCPC = matrix(unlist(BiasCPC.l), ncol = iters, byrow = F)
mat.varSPC = matrix(unlist(varSPC.l), ncol = iters, byrow = F)
mat.varCPC = matrix(unlist(varCPC.l), ncol = iters, byrow = F)
kmax = 20
mydat = data.frame(SPC = mat.coefSPC[kmax,],
                   CPC =  mat.coefCPC[kmax,])
mybiasdat = rbind(data.frame(bias = abs(mat.BiasSPC[kmax,]),'type' = 'PSC'),
                 data.frame(bias = abs(mat.BiasCPC[kmax,]),'type' = 'CPC') )
myvardat = rbind(data.frame(var = mat.varSPC[kmax,],'type' = 'PSC'),
                 data.frame(var =  mat.varCPC[kmax,],'type' = 'CPC') )
# main plots ================
bias<- ggplot(mybiasdat,aes(type, bias)) + geom_violin(trim=FALSE)+
   geom_boxplot(width=0.1) + xlab('')
vari <- ggplot(myvardat,aes(type, var)) + geom_violin(trim=FALSE)+
   geom_boxplot(width=0.1) + xlab('')+ylab('variance')
mytitle = paste('Structured X, n=',n,', p=',p,', sparse=',spa)
me.plot <- ggplot(data=A, aes(x=pcs, y=mean, colour=type) )+ 
  geom_point(lwd = 1) +   geom_line(lwd = 1)+ theme_bw()+
  ggtitle( mytitle)+ xlab('principal components')+
  scale_y_continuous(name=expression(beta[1]) )+ #,breaks=seq(0 , 1.2,.2)
  geom_ribbon(aes(ymin=A$low, ymax=A$high), linetype=6, alpha=.02, lwd = 1)+ 
  geom_hline(yintercept=1,linetype='dotdash')

library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
define_region <- function(row, col)  viewport(layout.pos.row = row, layout.pos.col = col)
# Arrange the plots
print(me.plot, vp = define_region(row = 1:2, col = 1:2)) 
print(bias, vp = define_region(row = 3, col = 1))
print(vari, vp = define_region(row = 3, col = 2))



