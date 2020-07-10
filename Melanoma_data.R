

rm(list = ls()) 

# Source Routines ---------------------------------------------------------

files.sources = list.files(path = "./routines", full.names = TRUE)
sapply(files.sources, source)

 
# Melanoma experiment -----------------------------------------------------
Celltypes=get(load("./Melanoma/Mela_Data_Types.rdata"))
Mela_Data_List=get(load("./Melanoma/Mela_Data.rdata"))

dataY1=scale(as.matrix(Mela_Data_List[[1]]),center=T,scale=T) #mean color
dataY2=scale(as.matrix(Mela_Data_List[[2]]),center=T,scale=T) #sd color 
dataY3=scale(as.matrix(Mela_Data_List[[3]]),center=T,scale=T) #mean shape
dataY4=scale(as.matrix(Mela_Data_List[[4]]),center=T,scale=T) #sd shape

dataY=cbind(dataY1,dataY2,dataY3,dataY4)
dataX1=Celltypes[[1]] #Only 2 gps
dataX=dataX1; gpb=c(1,1)
#dataX=as.matrix(scale(dataX1,center=T,scale=F)[,-2]); gpb=1

gp=c(rep(1,18),rep(2,18),rep(3,13),rep(4,13));  gv=4 

out1 <- cobs(dataY = dataY, 
             dataX = dataX1,  
             Y.bk.idx = gp, 
             X.bk.idx = gpb,
             n.comp = 1,
             alpha.v = 0,
             lambda.v = 0.3)

summary.cobs(out1, Y = dataY, X = dataX1)



####### Tune As in the paper ---------------------------------------------------------


Mela_out <- cobs.tune(dataY = dataY, 
                 dataX = dataX,  
                 Y.bk.idx = gp, 
                 X.bk.idx = gpb,
                 alpha.b = 1,
                 lambda.b = 0,
                 n.comp = 13,
                 alpha.v = 0.5,
                 method = "BIC-Low",
                 M.lambda.v = 40,
                 orth=F)
# Summary of cobs factorization 
summary.cobs(Mela_out,Y=dataY,X=dataX)

cobs_V=Mela_out$V
cobs_B=Mela_out$B
cobs_F=Mela_out$F.est

Score=dataY%*%cobs_V

# Inspection of block and variable selctions ------------------------------

out=Mela_out
# Block-wise segmentation (V)
cobs.out <- out
Y.idx <- cobs.out$Y.bk.idx
uniq.Y.idx <- unique(Y.idx)
K <- length(uniq.Y.idx)
r <- ncol(cobs.out$V) 

V.str <- matrix(nrow = K, ncol = r)
for( j in 1:K){
  if (r > 1){
    V.str[j,] <- colSums ( (cobs.out$V[Y.idx == uniq.Y.idx[j], ])^2  > 0 )
  }
  else{
    V.str[j] <- sum( (cobs.out$V[Y.idx == uniq.Y.idx[j], ])^2 > 0  )}
}
colnames(V.str) <- paste("Comp",1:r)
rownames(V.str) <- uniq.Y.idx
V.str
uniq.Y.idx




# Orth = TRUE option ------------------------------------------------------

# t(out$V) %*% out$V 
#  
# 
# out <- cobs.tune(dataY = dataY, 
#                   dataX = dataX,  
#                   Y.bk.idx = gp, 
#                   X.bk.idx = gpb,
#                   alpha.b = 1,
#                   lambda.b = 0,
#                   n.comp = 13,
#                   alpha.v = 0.5,
#                   method = "BIC-Low",
#                   M.lambda.v = 20,
#                   orth=T)
# # Summary of cobs factorization 
# summary.cobs(out)
# 
# t(out$V) %*% out$V

##### CHI plot, compare with PCA & SupSVD & AJIVE & SIFA & SLIDE & GFA#####
Mela_Cluster=as.vector(as.integer(Celltypes[[1]][,2]+1)) # Mela=1, Nevi=2
names(Mela_Cluster)=c(1:348)

#PCA
library(rARPACK)
Mela_PCA=svds(dataY,k=13)
Mela_PCA_Score=dataY%*%Mela_PCA$v

#SupSVD
Mela_SupSVD=SuperPCA::SupPCA(Y=as.matrix(scale(dataX1,center=T,scale=F)[,-2]), X=dataY, r=13)
Mela_SupSVD_Proj_Score=dataY%*%Mela_SupSVD$V
Mela_SupSVD_EST_Score=Mela_SupSVD$U

#AJIVE
# par(mfrow=c(2,2))
# plot(prcomp(dataY1)$sdev,type="b",xlab="component",ylab="",main="Y1")  #No very clear elbow
# plot(prcomp(dataY2)$sdev,type="b",xlab="component",ylab="",main="Y2")
# plot(prcomp(dataY3)$sdev,type="b",xlab="component",ylab="",main="Y3")
# plot(prcomp(dataY4)$sdev,type="b",xlab="component",ylab="",main="Y4")
min(which(cumsum((prcomp(dataY1)$sdev)^2)/sum((prcomp(dataY1)$sdev)^2)*100>70)) #4 #>75 5 #>80 6
min(which(cumsum((prcomp(dataY2)$sdev)^2)/sum((prcomp(dataY2)$sdev)^2)*100>70)) #5 >75 6 #>80 9
min(which(cumsum((prcomp(dataY3)$sdev)^2)/sum((prcomp(dataY3)$sdev)^2)*100>70)) #2 >75 3 #>80 4
min(which(cumsum((prcomp(dataY4)$sdev)^2)/sum((prcomp(dataY4)$sdev)^2)*100>70)) #4 >75 5 #>80 6
devtools::install_github("idc9/r_jive")
library(ajive)
library(cowplot)

Mela_AJIVE=ajive(blocks=list(dataY1,dataY2,dataY3,dataY4), initial_signal_ranks=c(4,5,2,4),
                 full=F,joint_rank=NA)
I1=as.matrix(Mela_AJIVE$block_decomps[[1]][['individual']][['v']])
I2=as.matrix(Mela_AJIVE$block_decomps[[2]][['individual']][['v']])
I3=as.matrix(Mela_AJIVE$block_decomps[[3]][['individual']][['v']])
I4=as.matrix(Mela_AJIVE$block_decomps[[4]][['individual']][['v']])
J1=as.matrix(Mela_AJIVE$block_decomps[[1]][['joint']][['v']])
J2=as.matrix(Mela_AJIVE$block_decomps[[2]][['joint']][['v']])
J3=as.matrix(Mela_AJIVE$block_decomps[[3]][['joint']][['v']])
J4=as.matrix(Mela_AJIVE$block_decomps[[4]][['joint']][['v']])

J_AJIVE=rbind(J1,J2,J3,J4)
dim(J_AJIVE);dim(I1);dim(I2);dim(I3);dim(I4) #Use these ranks to apply SIFA in Matlab
library(pracma)
AJIVE_V=as.matrix(cbind(J_AJIVE,blkdiag(I1,I2,I3,I4)))
AJIVE_V=AJIVE_V[,apply(AJIVE_V,2,sum)!=0]
Mela_AJIVE_V=apply(AJIVE_V,2,function(x) x/norm(x,"2")); dim(Mela_AJIVE_V)
apply(Mela_AJIVE_V,2,function(x) norm(x,"2"))

Mela_AJIVE_Score=dataY%*%Mela_AJIVE_V

# SIFA
#Mela_SIFA_B <- as.matrix(read.csv("~/Desktop/New Experiment/New_RealData_1028/Mela_SIFA_B.csv", header=FALSE, stringsAsFactors=FALSE))
#Mela_SIFA_V <- as.matrix(read.csv("~/Desktop/New Experiment/New_RealData_1028/Mela_SIFA_V.csv", header=FALSE, stringsAsFactors=FALSE))
Mela_SIFA_V <- as.matrix(Mela_SIFA_V)
Mela_SIFA_Score=dataY%*%Mela_SIFA_V

# COBS Score
Mela_SCARF_Score=Score


# SLIDE
source("./slide-paper-master/AuxillaryFunctions.R")
source("./slide-paper-master/SLIDEfunctions.R")

pvec= as.numeric(table(gp))
out_s <- standardizeX(dataY, pvec, center = T)
SLIDE_X <- out_s$X
svec <- out_s$svec
Mela_fit_seq<- 
  solve_optim1_seq_restarts(X=SLIDE_X, lambda_seq = NULL, pvec = pvec, 
                            k_max = 1000, eps = 1e-8, reduced = F, rank_total = 13,
                            n_lambda=10 ,lambda_max = max(svec), lambda_min = 0.1)
out_struct <- get_structure_v2(Mela_fit_seq, pvec)

# Select the best structure from the list
outbcv <- bcv_optim1_structure_centering(X=SLIDE_X, pvec = pvec, structure_list = out_struct$Slist, n_fold = 3, p_fold=3, k_max = 2000, eps = 1e-8, center = F)
SLIDE_structure <- outbcv$structure_min

Mela_SLIDE_param <- est_givenranks_v4(X = SLIDE_X, pvec = pvec, pattern = SLIDE_structure, k_max = 1000, eps = 1e-7)

slide_V <- Mela_SLIDE_param$V
slide_V=apply(slide_V,2,function(x) x/norm(x,"2"))

Mela_Slide_Score <- dataY%*%slide_V
#Mela_Slide_Score <- Mela_SLIDE_param$U

# GFA
gfa_opts = CCAGFA::getDefaultOpts()
gfa_opts$verbose=0 # do not print iteration details in GFA
fit_gfa = CCAGFA::GFAexperiment(Y = list(dataY1,dataY2,dataY3,dataY4),
                                      K = 13, opts = gfa_opts, Nrep=10)
fit_gfa_trim = CCAGFA::GFAtrim(fit_gfa)
GFA_V=do.call("rbind", fit_gfa_trim$W)
GFA_V=apply(GFA_V,2,function(x) x/norm(x,"2")) # make norm 1

Mela_GFA_Score <- dataY%*%GFA_V

# Calinski-Harabasz Index: SS_B/SS_W*(N-k)/(k-1), max is the best
library(fpc)

CHI=list()
CHI_Single=list()
for (i in 1:13){
  CHI[['SCARF_CHI']][i]=calinhara(Mela_SCARF_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['SCARF']][i]=calinhara(Mela_SCARF_Score[,i],Mela_Cluster,2)
  CHI[['PCA_CHI']][i]=calinhara(Mela_PCA_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['PCA']][i]=calinhara(Mela_PCA_Score[,i],Mela_Cluster,2)
  CHI[['SupSVD_Proj_CHI']][i]=calinhara(Mela_SupSVD_Proj_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['SupSVD_Proj']][i]=calinhara(Mela_SupSVD_Proj_Score[,i],Mela_Cluster,2)
  CHI[['SupSVD_Est_CHI']][i]=calinhara(Mela_SupSVD_EST_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['SupSVD_Est']][i]=calinhara(Mela_SupSVD_EST_Score[,i],Mela_Cluster,2)
  CHI[['AJIVE_CHI']][i]=calinhara(Mela_AJIVE_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['AJIVE']][i]=calinhara(Mela_AJIVE_Score[,i],Mela_Cluster,2)
  CHI[['SIFA_CHI']][i]=calinhara(Mela_SIFA_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['SIFA']][i]=calinhara(Mela_SIFA_Score[,i],Mela_Cluster,2)
  CHI[['SLIDE_CHI']][i]=calinhara(Mela_Slide_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['SLIDE']][i]=calinhara(Mela_Slide_Score[,i],Mela_Cluster,2)
  CHI[['GFA_CHI']][i]=calinhara(Mela_GFA_Score[,1:i],Mela_Cluster,2)
  CHI_Single[['GFA']][i]=calinhara(Mela_GFA_Score[,i],Mela_Cluster,2)
}

par(mfrow=c(1,1))
plot(CHI_Single[['SLIDE']],pch="",ylab="Calinski Harabasz Index",
     xlab="Component",bty="n",main="Single Component CHI - Sorted",ylim=c(0,250))
points(sort(CHI_Single[['SCARF']],decreasing = T),pch=19,type="b",col="red",lwd=3)
#points(sort(CHI_Single[['PCA']],decreasing = T),pch=22,type="b",col="blue")
points(sort(CHI_Single[['SLIDE']],decreasing = T),pch=22,type="b",col="blue",lty=2)
points(sort(CHI_Single[['GFA']],decreasing = T),pch=24,type="b",col="darkgreen",lty=3)
points(sort(CHI_Single[['AJIVE']],decreasing = T),pch=21,type="b",col="gold3",lty=4)
points(sort(CHI_Single[['SupSVD_Est']],decreasing = T),pch=25,type="b",col="darkred",lty=5)
points(sort(CHI_Single[['SIFA']],decreasing = T),pch=23,type="b",col="darkorchid",lty=6)

legend("topright", legend=c("COBS","SLIDE","GFA","AJIVE","SupSVD","SIFA"),
       col=c("red", "blue","darkgreen","gold3","darkred","darkorchid"), lty=c(1:6),pch=c(19,22,24,21,25,23), cex=0.8)

#quartz.save(file="COBS_Mela_CHI_Sorted.png", type = "png", device = dev.cur(),
#            dpi = 300, width=10, height=7)

par(mfrow=c(1,1))
plot(CHI[['SLIDE_CHI']],pch="",ylab="Calinski Harabasz Index",
     xlab="Number of Components",bty="n",main="Cumulative Components CHI",ylim=c(0,180))
points(CHI[['SCARF_CHI']],pch=19,type="b",col="red",lwd=3)
points(CHI[['SLIDE_CHI']],pch=22,type="b",col="blue",lty=2)
points(CHI[['GFA_CHI']],pch=24,type="b",col="darkgreen",lty=3)
points(CHI[['AJIVE_CHI']],pch=21,type="b",col="gold3",lty=4)
points(CHI[['SupSVD_Est_CHI']],pch=25,type="b",col="darkred",lty=5)
points(CHI[['SIFA_CHI']],pch=23,type="b",col="darkorchid",lty=6)

legend("topright", legend=c("COBS","SLIDE","GFA","AJIVE","SupSVD","SIFA"),
       col=c("red", "blue","darkgreen","gold3","darkred","darkorchid"), lty=c(1:6),pch=c(19,22,24,21,25,23), cex=0.8)
#quartz.save(file="COBS_Mela_CHI.png", type = "png", device = dev.cur(),
#            dpi = 300, width=10, height=7)

#### Variance Explained #####

# Overall Variation explained proportion  
SSTO=norm(dataY,"F")^2
SSR_Supervision=norm(dataX%*%cobs_B%*%t(cobs_V),"F")^2
SSR_Unknown=norm(cobs_F%*%t(cobs_V),"F")^2
SSE=SSTO-(SSR_Supervision+SSR_Unknown)

SSR_Supervision/SSTO # 0.06252173
SSR_Unknown/SSTO # [1] 0.4286798
SSE/SSTO # [1] 0.5087985

# Mean Color variation explained proportion
SSTO_MC=norm(dataY1,"F")^2
SSX_MC=norm(dataX%*%cobs_B%*%t(cobs_V[1:18,]),"F")^2
SSF_MC=norm(cobs_F%*%t(cobs_V[1:18,]),"F")^2
SSE_MC=SSTO_MC-SSX_MC-SSF_MC

SSX_MC/SSTO_MC # [1] 0.04553838
SSF_MC/SSTO_MC # [1] 0.4867047
SSE_MC/SSTO_MC # [1] 0.4677569

# S.d. Color variation explained proportion
SSTO_SdC=norm(dataY2,"F")^2
SSX_SdC=norm(dataX%*%cobs_B%*%t(cobs_V[19:36,]),"F")^2
SSF_SdC=norm(cobs_F%*%t(cobs_V[19:36,]),"F")^2
SSE_SdC=SSTO_SdC-SSX_SdC-SSF_SdC

SSX_SdC/SSTO_SdC # [1] 0.02460716
SSF_SdC/SSTO_SdC # [1] 0.3487126
SSE_SdC/SSTO_SdC # [1] 0.6266802

# Mean Shape variation explained proportion
SSTO_MS=norm(dataY3,"F")^2
SSX_MS=norm(dataX%*%cobs_B%*%t(cobs_V[37:49,]),"F")^2
SSF_MS=norm(cobs_F%*%t(cobs_V[37:49,]),"F")^2
SSE_MS=SSTO_MS-SSX_MS-SSF_MS

SSX_MS/SSTO_MS # [1] 0.1784489
SSF_MS/SSTO_MS # [1] 0.4987097
SSE_MS/SSTO_MS # [1] 0.3228414

# Sd Shape variation explained proportion
SSTO_SdS=norm(dataY4,"F")^2 
SSX_SdS=norm(dataX%*%cobs_B%*%t(cobs_V[50:62,]),"F")^2
SSF_SdS=norm(cobs_F%*%t(cobs_V[50:62,]),"F")^2
SSE_SdS=SSTO_SdS-SSX_SdS-SSF_SdS

SSX_SdS/SSTO_SdS # [1] 0.02260705
SSF_SdS/SSTO_SdS # [1] 0.3890315
SSE_SdS/SSTO_SdS # [1] 0.5883615

