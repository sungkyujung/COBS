
rm(list = ls())
library(RColorBrewer) #Heatmap color
require(lattice)

##### Toy data, Check Figure 1 in Notes : Y=UV+E #####

n=100  # Sample size, signal
p=120 # Dim of Y,V
q=40   # Dim of X,B
r=4    # r is the intrinsic rank of Y

gv=3; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b

# U=XB+F, but just plot U for the toy data, no need to simulate B and F....
set.seed(1234)
U1=as.matrix(sort(rnorm(n,0,15),decreasing=F)); sd(U1)
for_u2=sort(rnorm(n,0,2),decreasing=F); U2=as.matrix(1.5*for_u2^2+for_u2); sd(U2)
for_u3=sort(runif(n,-1,1),decreasing=F); U3=as.matrix(100*(for_u3^4-for_u3^2)); sd(U3)
U4=as.matrix(rnorm(n,0,5)); sd(U4)
U=as.matrix(cbind(U1,U2,U3,U4))


#plot(U4,pch=20,xlab="",ylab="",cex.axis=0.5)
# quartz.save(file="U4.png", type = "png", device = dev.cur(),
#             dpi = 300, width=7, height=4)

# V
V1=rep(1,p); V1=V1/norm(V1,"2")     # Full-Joint
V2=c(rep(1,20),rep(-1,20),rep(1,20),rep(-1,20),rep(0,40)); V2=V2/norm(V2,"2")   # Partial-Joint in block 1 and 2
V3=c(rep(c(1,-1),20),rep(0,80)); V3=V3/norm(V3,"2")   # Individual in block 1
V4=c(rep(0,80),rep(c(1,-1),20)); V4=V4/norm(V4,"2")   # Individual in block 3
V=as.matrix(cbind(V1,V2,V3,V4)); apply(V,2,function(x) norm(x,"2"))

gp=ceiling(1:p/ size.groupv)

# E
set.seed(1234)
E=matrix(rnorm(n*p,0,1),nrow=n)


# Seperate V, Y into blocks;
Signal=U%*%t(V)
Signal1=Signal[,1:40]; Signal2=Signal[,41:80]; Signal3=Signal[,81:120]

toy_fullJ=U1%*%t(V1); toy_partialJ=U2%*%t(V2); toy_ind_1=U3%*%t(V3); toy_ind_2=U4%*%t(V4)

Toy_Y=Signal+E
Toy_Y1=as.matrix(Toy_Y[,1:40]); Toy_Y2=as.matrix(Toy_Y[,41:80]); Toy_Y3=as.matrix(Toy_Y[,81:120]) 

##### Heatmaps #####
par(mar=c(4,4,4,4),mgp=c(2,0.8,0),oma=c(0,0,2,0))
# quartz.save(file="./ToyDataPlot/Ind2.png", type = "png", device = dev.cur(),
#             dpi = 300, width=5, height=5)
# Error Heatmap
levelplot(t(E[n:1,]),  
          xlab="", ylab="", main="", 
          col.regions=colorRampPalette(brewer.pal(9, "RdGy"))(50), colorkey=F,
          scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
          panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = c(40,80),lty=2)
          })

# Y Heatmaps
levelplot(t(Toy_Y[n:1,]), 
        xlab="", ylab="", main="", 
        scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
        col.regions=colorRampPalette(brewer.pal(9, "RdBu"))(50), colorkey=F,
        panel = function(...){
                panel.levelplot(...)
                panel.abline(v = c(40,80),lty=2)
        })

# heatmap(Toy_Y1,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
#         xlab=NULL, ylab=NULL, margins=c(4,4), main="Y1", 
#         col= colorRampPalette(brewer.pal(9, "RdBu"))(50))
# 
# heatmap(Toy_Y2,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
#         xlab=NULL, ylab=NULL, margins=c(4,4), main="Y2", 
#         col= colorRampPalette(brewer.pal(9, "RdBu"))(50))
# 
# heatmap(Toy_Y3,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
#         xlab=NULL, ylab=NULL, margins=c(4,4), main="Y3", 
#         col= colorRampPalette(brewer.pal(9, "RdBu"))(50))

# Signal Heatmap

levelplot(t(Signal[n:1,]), 
          xlab="", ylab="", main="", 
          col.regions=colorRampPalette(brewer.pal(9, "RdBu"))(50), colorkey=F,
          scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
          panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = c(40,80),lty=2)
          })

# Signal Heatmap - Patterns
levelplot(t(toy_fullJ[n:1,]),
          xlab="", ylab="", main="", 
          col.regions=colorRampPalette(brewer.pal(9, "RdBu"))(50), colorkey=F,
          scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
          panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = c(40,80),lty=2)
          })

levelplot(t(toy_partialJ[n:1,]), 
          xlab="", ylab="", main="", 
          col.regions=colorRampPalette(brewer.pal(9, "RdBu"))(50), colorkey=F,
          scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
          panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = c(40,80),lty=2)
          })

levelplot(t(toy_ind_1[n:1,]), 
          xlab="", ylab="", main="", 
          col.regions=colorRampPalette(brewer.pal(9, "RdBu"))(50), colorkey=F,
          scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
          panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = c(40,80),lty=2)
          })

levelplot(t(toy_ind_2[n:1,]), 
          xlab="", ylab="", main="", 
          col.regions=colorRampPalette(brewer.pal(9, "RdBu"))(50), colorkey=F,
          scales=list(x=list(at=c(40,80),label=c(40,80)),y=list(at=c(20,40,60,80),label=c(80,60,40,20))),
          panel = function(...){
                  panel.levelplot(...)
                  panel.abline(v = c(40,80),lty=2)
          })

# Signal Heatmaps - Blocks
# heatmap(Signal1,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
#         xlab=NULL, ylab=NULL, margins=c(4,4), main="Signal-1", 
#         col= colorRampPalette(brewer.pal(9, "RdBu"))(50))
# heatmap(Signal2,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
#         xlab=NULL, ylab=NULL, margins=c(4,4), main="Signal-2", 
#         col= colorRampPalette(brewer.pal(9, "RdBu"))(50))
# heatmap(Signal3,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
#         xlab=NULL, ylab=NULL, margins=c(4,4), main="Signal-3", 
#         col= colorRampPalette(brewer.pal(9, "RdBu"))(50))

# U and V maps
heatmap(V,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
        xlab=NULL, ylab=NULL, margins=c(4,4), main="", 
        col= colorRampPalette(brewer.pal(9, "RdBu"))(50))

heatmap(U,scale = "none",Rowv = NA,Colv = NA, cexRow = 0.5,cexCol = 0.5,
        xlab=NULL, ylab=NULL, margins=c(4,4), main="", 
        col= colorRampPalette(brewer.pal(9, "RdBu"))(50))





