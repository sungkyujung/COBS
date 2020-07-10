##### BRCA New experient: Train and Test #######

rm(list = ls()) 

# Source Routines ---------------------------------------------------------

files.sources = list.files(path = "./routines", full.names = TRUE)
sapply(files.sources, source)


# Read the BRCA data ------------------------------------------------------

source("Read_BRCA.R")
dataY1 <-dataY[,1:200] # Gene Expression  
dataY2 <-dataY[,201:400] # Methylation
dataY3 <-dataY[,401:600] # CNV

#  Half training and half testing ####
table(Pam50)
set.seed(12345)
Basal_train <- sample(which(Pam50=="Basal"),ceiling(133*0.5),replace=F)
Her2_train <- sample(which(Pam50=="Her2"),ceiling(51*0.5),replace=F)
LumA_train <- sample(which(Pam50=="LumA"),ceiling(394*0.5),replace=F)
LumB_train <- sample(which(Pam50=="LumB"),ceiling(161*0.5),replace=F)
Normal_train <- sample(which(Pam50=="Normal"),ceiling(31*0.5),replace=F)
train_samples <- sort(c(Basal_train,Her2_train,LumA_train,LumB_train,Normal_train))
dataY_train <- scale(dataY[train_samples,],center = T,scale = F)
dataX_train <- dataX[train_samples,]
dataY_test <- scale(dataY[-train_samples,],center=T,scale=F)
dataX_test <- dataX[-train_samples,]
Subtypes_train <- Pam50[train_samples,]
Subtypes_test <- Pam50[-train_samples,]

##### COBS #####
BRCA_train <- cobs.tune(dataY = dataY_train, 
                      dataX = dataX_train,  
                      Y.bk.idx = gp, 
                      X.bk.idx = gpb,
                      n.comp = 2,  # n.comp=62 to calculate the variation explained table
                      alpha.b = 1,
                      lambda.b = 0.0025,
                      alpha.v = 0,
                      method = "BIC",
                      M.lambda.v = 20,
                      orth=F)
summary_train <- summary.cobs(BRCA_train, Y = dataY_train, X = dataX_train)

cobs_V_train <- BRCA_train$V
cobs_B_train <- BRCA_train$B
cobs_F_train <- BRCA_train$F.est

Score_test_cobs <- dataY_test%*%cobs_V_train
round(diag(var(Score_test_cobs)),2)

Score_test_cobs_scatter <- data.frame(Score_test_cobs[,1:2],Subtypes_test)
colnames(Score_test_cobs_scatter)=c("Score1", "Score2","Subtypes")
colors <- Subtypes_test
library(ggpubr)
ScoreScatter_cobs <- ggscatter(Score_test_cobs_scatter, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",
                               #xlim=c(-15,25),ylim=c(-15,40)
                               ) + border()
ScoreScatter_cobs

# Cluster of the First 2 component Score:
d_cobs <- dist(Score_test_cobs[,1:2], method = "euclidean")
hc_cobs <- hclust(d_cobs, method = "ward.D2")
plot(hc_cobs,labels=F)
sub_grp_cobs <- cutree(hc_cobs, k = 5)
rect.hclust(hc_cobs, k = 5, border = 2:6)
tab_cobs=list
for(i in 1:5){
print(table(Score_test_cobs_scatter$Subtypes[which(sub_grp_cobs==i)]))
}

mclust::adjustedRandIndex(Subtypes_test, sub_grp_cobs)
#fpc::calinhara(Score_test_cobs[,5],as.numeric(Subtypes_test),5)

##### SLIDE #####
source("~/Desktop/slide-paper-master/AuxillaryFunctions.R")
source("~/Desktop/slide-paper-master/SLIDEfunctions.R")

pvec=c(200,200,200)
out_s <- standardizeX(dataY_train, pvec, center = T)
SLIDE_X <- out_s$X
svec <- out_s$svec
BRCA_fit_seq<- 
  solve_optim1_seq_restarts(X=SLIDE_X, lambda_seq = NULL, pvec = pvec, 
                            k_max = 1000, eps = 1e-8, reduced = F, rank_total = 5,
                            n_lambda=10 ,lambda_max = max(svec), lambda_min = 0.1)
out_struct <- get_structure_v2(BRCA_fit_seq, pvec)

# Select the best structure from the list
outbcv <- bcv_optim1_structure_centering(X=SLIDE_X, pvec = pvec, structure_list = out_struct$Slist, n_fold = 3, p_fold=3, k_max = 2000, eps = 1e-8, center = F)
SLIDE_structure <- outbcv$structure_min

BRCA_SLIDE_param <- est_givenranks_v4(X = SLIDE_X, pvec = pvec, pattern = SLIDE_structure, k_max = 1000, eps = 1e-7)

slide_V_train <- BRCA_SLIDE_param$V
slide_V_train=apply(slide_V_train,2,function(x) x/norm(x,"2"))
apply(slide_V_train,2,function(x) norm(x,"2"))

Score_test_slide <- dataY_test%*%slide_V_train
round(diag(var(Score_test_slide)),2)

Score_test_slide_scatter <- data.frame(Score_test_slide[,1:2],Subtypes_test)
colnames(Score_test_slide_scatter)=c("Score1", "Score2","Subtypes")
colors <- Subtypes_test
library(ggpubr)
ScoreScatter_slide <- ggscatter(Score_test_slide_scatter, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",
                                #xlim=c(-10,20),ylim=c(-15,10)
                                ) + border()
ScoreScatter_slide

# Cluster of the First 10 component Score:
d_slide<- dist(Score_test_slide[,1:2], method = "euclidean")
hc_slide <- hclust(d_slide, method = "ward.D2")
plot(hc_slide,labels=F)
sub_grp_slide <- cutree(hc_slide, k = 5)
rect.hclust(hc_slide, k = 5, border = 2:6)
for(i in 1:5){
  print(table(Score_test_slide_scatter$Subtypes[which(sub_grp_slide==i)]))
}

mclust::adjustedRandIndex(Subtypes_test, sub_grp_slide)
#fpc::calinhara(Score_test_slide[,1],as.numeric(Subtypes_test),5)

##### GFA #####
gfa_opts = CCAGFA::getDefaultOpts()
gfa_opts$verbose=0 # do not print iteration details in GFA
fit_gfa_train = CCAGFA::GFAexperiment(Y = list(dataY_train[,1:200],dataY_train[,201:400],dataY_train[,401:600]),
                                K = 10, opts = gfa_opts, Nrep=10)
fit_gfa_trim_train = CCAGFA::GFAtrim(fit_gfa_train)
GFA_V_train=do.call("rbind", fit_gfa_trim_train$W)
GFA_V_train=as.matrix(apply(GFA_V_train,2,function(x) x/norm(x,"2"))) # make norm 1

GFA_structure <- matrix(NA, ncol=10, nrow=3)
for (i in 1:10){
  for (j in 1:3){
    GFA_structure[j,i]=ifelse(all(GFA_V_train[which(gp==j),i]==0),0,1)
  }
}

Score_test_gfa <- dataY_test%*%GFA_V_train
diag(var(Score_test_gfa))  #PC9 and #PC10 have largest variances
sort(round(diag(var(Score_test_gfa)),2),decreasing = T)

GFA_variance <- rbind(GFA_structure, round(diag(var(Score_test_gfa)),2))
GFA_variance[,order(-GFA_variance[4,])]

Score_test_gfa_scatter <- data.frame(Score_test_gfa[,1:2],Subtypes_test) 
colnames(Score_test_gfa_scatter)=c("Score1", "Score2","Subtypes")
colors <- Subtypes_test
library(ggpubr)
ScoreScatter_gfa <- ggscatter(Score_test_gfa_scatter, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",
                              #xlim=c(-15,25),ylim=c(-15,30)
                              ) + border()
ScoreScatter_gfa

# Cluster of the First 10 component Score:
d_gfa<- dist(Score_test_gfa[,c(10,9)], method = "euclidean") # The PCs that have the largest variances
hc_gfa <- hclust(d_gfa, method = "ward.D2")
plot(hc_gfa,labels=F)
sub_grp_gfa <- cutree(hc_gfa, k = 5)
rect.hclust(hc_gfa, k = 5, border = 2:6)
for(i in 1:5){
  print(table(Score_test_gfa_scatter$Subtypes[which(sub_grp_gfa==i)]))
}

mclust::adjustedRandIndex(Subtypes_test, sub_grp_gfa)
#fpc::calinhara(Score_test_gfa[,2],as.numeric(Subtypes_test),5)

##### AJIVE #####

min(which(cumsum((prcomp(dataY_train[,1:200])$sdev)^2)/sum((prcomp(dataY_train[,1:200])$sdev)^2)*100>60)) #13
min(which(cumsum((prcomp(dataY_train[,201:400])$sdev)^2)/sum((prcomp(dataY_train[,201:400])$sdev)^2)*100>60)) #12
min(which(cumsum((prcomp(dataY_train[,401:600])$sdev)^2)/sum((prcomp(dataY_train[,401:600])$sdev)^2)*100>65)) #4
devtools::install_github("idc9/r_jive")
library(ajive)
library(cowplot)

BRCA_AJIVE_Train=ajive(blocks=list(dataY_train[,1:200],dataY_train[,201:400],dataY_train[,401:600]), initial_signal_ranks=c(13,12,4),
                 full=F,joint_rank=NA)
I1=as.matrix(BRCA_AJIVE_Train$block_decomps[[1]][['individual']][['v']])
I2=as.matrix(BRCA_AJIVE_Train$block_decomps[[2]][['individual']][['v']])
I3=as.matrix(BRCA_AJIVE_Train$block_decomps[[3]][['individual']][['v']])

J1=as.matrix(BRCA_AJIVE_Train$block_decomps[[1]][['joint']][['v']])
J2=as.matrix(BRCA_AJIVE_Train$block_decomps[[2]][['joint']][['v']])
J3=as.matrix(BRCA_AJIVE_Train$block_decomps[[3]][['joint']][['v']])

J_AJIVE=rbind(J1,J2,J3)
dim(J_AJIVE);dim(I1);dim(I2);dim(I3) #Use these ranks to apply SIFA
library(pracma)
AJIVE_V_train=as.matrix(cbind(J_AJIVE,blkdiag(I1,I2,I3)))
AJIVE_V_train=AJIVE_V_train[,apply(AJIVE_V_train,2,sum)!=0]
ajive_V_train=apply(AJIVE_V_train,2,function(x) x/norm(x,"2")); dim(ajive_V_train)
apply(ajive_V_train,2,function(x) norm(x,"2"))

Score_test_ajive <- dataY_test%*%ajive_V_train
diag(var(Score_test_ajive)) #PC1 and PC24 have the largest variances
sort(round(diag(var(Score_test_ajive)),2),decreasing = T)

AJIVE_structure <- matrix(NA, ncol=26, nrow=3)
for (i in 1:26){
  for (j in 1:3){
    AJIVE_structure[j,i]=ifelse(all(ajive_V_train[which(gp==j),i]==0),0,1)
  }
}

ajive_variance <- rbind(AJIVE_structure, round(diag(var(Score_test_ajive)),2))
ajive_variance[,order(-ajive_variance[4,])[1:5]]

Score_test_ajive_scatter <- data.frame(Score_test_ajive[,1:2],Subtypes_test) # score 4 and 5 have the largest variance
colnames(Score_test_ajive_scatter)=c("Score1", "Score2","Subtypes")
colors <- Subtypes_test
library(ggpubr)
ScoreScatter_ajive <- ggscatter(Score_test_ajive_scatter, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",
                               # xlim=c(-15,25),ylim=c(-20,10)
                                ) + border()
ScoreScatter_ajive

# Cluster of the First 10 component Score:
d_ajive<- dist(Score_test_ajive[,c(1,24)], method = "euclidean")   # Use the PCs that have the largest variances
hc_ajive <- hclust(d_ajive, method = "ward.D2")
plot(hc_ajive,labels=F)
sub_grp_ajive <- cutree(hc_ajive, k = 5)
rect.hclust(hc_ajive, k = 5, border = 2:6)
for(i in 1:5){
  print(table(Score_test_ajive_scatter$Subtypes[which(sub_grp_ajive==i)]))
}

mclust::adjustedRandIndex(Subtypes_test, sub_grp_ajive)
#fpc::calinhara(Score_test_ajive[,1],as.numeric(Subtypes_test),5)

##### SupSVD #####
SupSVD_train <- SuperPCA::SupPCA(Y=as.matrix(scale(dataX_train,center=T,scale=F)[,-5]), 
                                 X=dataY_train, r=10)
SupSVD_V_train <- SupSVD_train$V

Score_test_supsvd <- dataY_test%*%SupSVD_V_train
round(diag(var(Score_test_supsvd)),2)

Score_test_supsvd_scatter <- data.frame(Score_test_supsvd[,1:2],Subtypes_test) # score 4 and 5 have the largest variance
colnames(Score_test_supsvd_scatter)=c("Score1", "Score2","Subtypes")
colors <- Subtypes_test
library(ggpubr)
ScoreScatter_supsvd <- ggscatter(Score_test_supsvd_scatter, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",
                                # xlim=c(-15,25),ylim=c(-20,10)
) + border()
ScoreScatter_supsvd

# Cluster of the First 2 component Score:
d_supsvd<- dist(Score_test_supsvd[,1:2], method = "euclidean")   # Use the first 2 components
hc_supsvd <- hclust(d_supsvd, method = "ward.D2")
plot(hc_supsvd,labels=F)
sub_grp_supsvd <- cutree(hc_supsvd, k = 5)
rect.hclust(hc_supsvd, k = 5, border = 2:6)
for(i in 1:5){
  print(table(Score_test_supsvd_scatter$Subtypes[which(sub_grp_supsvd==i)]))
}

mclust::adjustedRandIndex(Subtypes_test, sub_grp_supsvd)
#fpc::calinhara(Score_test_supsvd[,1],as.numeric(Subtypes_test),5)

##### SIFA - error occured using r package #####
SIFA_train <- SuperPCA::SIFA(X=dataX_train,
                             Y=list(dataY_train[,1:200],dataY_train[,201:400],dataY_train[,401:600]),
                             r0=4,r=c(10,9,3),sparsity=1,type="A")
# Using Matlab
#write.table(dataX_train,file="./For_SIFA_Matlab_BRCA/DataX_train.csv",sep=" ",col.names=F,row.names=F)
#write.table(dataY_train,file="./For_SIFA_Matlab_BRCA/DataY_train.csv",sep=" ",col.names=F,row.names=F)
SIFA_B_train <- as.matrix(read.csv("./For_SIFA_Matlab_BRCA/BRCA_SIFA_B_train.csv", header=FALSE, stringsAsFactors=FALSE))
SIFA_V_train <- as.matrix(read.csv("./For_SIFA_Matlab_BRCA/BRCA_SIFA_V_train.csv", header=FALSE, stringsAsFactors=FALSE))

SIFA_structure <- matrix(NA, ncol=26, nrow=3)
for (i in 1:26){
  for (j in 1:3){
    SIFA_structure[j,i]=ifelse(all(SIFA_V_train[which(gp==j),i]==0),0,1)
  }
}

Score_test_sifa <- dataY_test%*%SIFA_V_train
diag(var(Score_test_sifa)) #PC1 and PC15 have the largest variances
sort(round(diag(var(Score_test_sifa)),2),decreasing = T)

SIFA_variance <- rbind(SIFA_structure, round(diag(var(Score_test_sifa)),2))
SIFA_variance[,order(-SIFA_variance[4,])[1:5]]

Score_test_sifa_scatter <- data.frame(Score_test_sifa[,1:2],Subtypes_test) # score 4 and 5 have the largest variance
colnames(Score_test_sifa_scatter)=c("Score1", "Score2","Subtypes")
colors <- Subtypes_test
library(ggpubr)
ScoreScatter_sifa <- ggscatter(Score_test_sifa_scatter, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",
                                 # xlim=c(-15,25),ylim=c(-20,10)
                                ) + border()
ScoreScatter_sifa

# Cluster of the First 10 component Score:
d_sifa<- dist(Score_test_sifa[,c(1,15)], method = "euclidean")   # Use the first 5 components
hc_sifa <- hclust(d_sifa, method = "ward.D2")
plot(hc_sifa,labels=F)
sub_grp_sifa <- cutree(hc_sifa, k = 5)
rect.hclust(hc_sifa, k = 5, border = 2:6)
for(i in 1:5){
  print(table(Score_test_sifa_scatter$Subtypes[which(sub_grp_sifa==i)]))
}

mclust::adjustedRandIndex(Subtypes_test, sub_grp_sifa)
#fpc::calinhara(Score_test_sifa[,2],as.numeric(Subtypes_test),5)


##### CHI Score of projected test data for above method #####
library(fpc)

CHI=list()
CHI_Single=list()
for (i in 1:10){
  CHI[['COBS_CHI']][i]=calinhara(Score_test_cobs[,1:i],as.numeric(Subtypes_test),5)
  CHI_Single[['COBS']][i]=calinhara(Score_test_cobs[,i],as.numeric(Subtypes_test),5)
  
  CHI[['SLIDE_CHI']][i]=calinhara(Score_test_slide[,1:i],as.numeric(Subtypes_test),5)
  CHI_Single[['SLIDE']][i]=calinhara(Score_test_slide[,i],as.numeric(Subtypes_test),5)
  
  CHI[['GFA_CHI']][i]=calinhara(Score_test_gfa[,1:i],as.numeric(Subtypes_test),5)
  CHI_Single[['GFA']][i]=calinhara(Score_test_gfa[,i],as.numeric(Subtypes_test),5)
  
  CHI[['AJIVE_CHI']][i]=calinhara(Score_test_ajive[,1:i],as.numeric(Subtypes_test),5)
  CHI_Single[['AJIVE']][i]=calinhara(Score_test_ajive[,i],as.numeric(Subtypes_test),5)
  
  CHI[['SupSVD_CHI']][i]=calinhara(Score_test_supsvd[,1:i],as.numeric(Subtypes_test),5)
  CHI_Single[['SupSVD']][i]=calinhara(Score_test_supsvd[,i],as.numeric(Subtypes_test),5)
  
  CHI[['SIFA_CHI']][i]=calinhara(Score_test_sifa[,1:i],as.numeric(Subtypes_test),5)
  CHI_Single[['SIFA']][i]=calinhara(Score_test_sifa[,i],as.numeric(Subtypes_test),5)
}

par(mfrow=c(1,1))
plot(CHI[['COBS_CHI']],pch="",
     ylab="Calinski Harabasz Index",xlab="Component",
     bty="n",main="Cumulative Components CHI",
     ylim=c(-10,700)
     )
points(CHI[['COBS_CHI']],pch=19,type="b",col="red",lwd=3,lty=1)
points(CHI[['SLIDE_CHI']],pch=22,type="b",col="blue",lty=2)
points(CHI[['GFA_CHI']],pch=24,type="b",col="dark green",lty=3)
points(CHI[['AJIVE_CHI']],pch=21,type="b",col="gold3",lty=4)
points(CHI[['SupSVD_CHI']],pch=25,type="b",col="darkred",lty=5)
points(CHI[['SIFA_CHI']],pch=23,type="b",col="darkorchid",lty=6)
legend("topright", legend=c("COBS","SLIDE","GFA","AJIVE","SupSVD","SIFA"),
       col=c("red","blue","dark green","gold3","darkred","darkorchid"), 
       lty=c(1:6),pch=c(19,22,24,21,25,23), cex=0.8)

# quartz.save(file="BRCA_CHI_Cumulative.png", type = "png", device = dev.cur(),
#             dpi = 300, width=10, height=7)

plot(CHI_Single[['COBS']],pch="",
     ylab="Calinski Harabasz Index",xlab="Component",
     bty="n",main="Single Component CHI - Sorted",
     #ylim=c(0,200)
)
points(sort(CHI_Single[['COBS']],decreasing = T),pch=19,type="b",col="red",lwd=3,lty=1)
points(sort(CHI_Single[['SLIDE']],decreasing = T),pch=22,type="b",col="blue",lty=2)
points(sort(CHI_Single[['GFA']],decreasing = T),pch=24,type="b",col="dark green",lty=3)
points(sort(CHI_Single[['AJIVE']],decreasing = T),pch=21,type="b",col="gold3",lty=4)
points(sort(CHI_Single[['SupSVD']],decreasing = T),pch=25,type="b",col="darkred",lty=5)
points(sort(CHI_Single[['SIFA']],decreasing = T),pch=23,type="b",col="darkorchid",lty=6)
legend("topright", legend=c("COBS","SLIDE","GFA","AJIVE","SupSVD","SIFA"),
       col=c("red","blue","dark green","gold3","darkred","darkorchid"), 
       lty=c(1:6),pch=c(19,22,24,21,25,23), cex=0.8)

# quartz.save(file="BRCA_CHI_Single.png", type = "png", device = dev.cur(),
#             dpi = 300, width=10, height=7)


##### ARI #####
ARI=list()
for (i in 1:10){
    d_cobs <- dist(Score_test_cobs[,1:i], method = "euclidean")
    hc_cobs <- hclust(d_cobs, method = "ward.D2")
    sub_grp_cobs <- cutree(hc_cobs, k = 5)
  ARI[['COBS']][i]=mclust::adjustedRandIndex(Subtypes_test, sub_grp_cobs)
    d_slide<- dist(Score_test_slide[,1:i], method = "euclidean")
    hc_slide <- hclust(d_slide, method = "ward.D2")
    sub_grp_slide <- cutree(hc_slide, k = 5)
  ARI[['SLIDE']][i]=mclust::adjustedRandIndex(Subtypes_test, sub_grp_slide)
    d_gfa<- dist(Score_test_gfa[,1:i], method = "euclidean")
    hc_gfa <- hclust(d_gfa, method = "ward.D2")
    sub_grp_gfa <- cutree(hc_gfa, k = 5)
  ARI[['GFA']][i]=mclust::adjustedRandIndex(Subtypes_test, sub_grp_gfa)
    d_ajive<- dist(Score_test_ajive[,1:i], method = "euclidean")
    hc_ajive <- hclust(d_ajive, method = "ward.D2")
    sub_grp_ajive <- cutree(hc_ajive, k = 5)
  ARI[['AJIVE']][i]=mclust::adjustedRandIndex(Subtypes_test, sub_grp_ajive)
    d_supsvd<- dist(Score_test_supsvd[,1:i], method = "euclidean")
    hc_supsvd <- hclust(d_supsvd, method = "ward.D2")
    sub_grp_supsvd <- cutree(hc_supsvd, k = 5)
  ARI[['SupSVD']][i]=mclust::adjustedRandIndex(Subtypes_test, sub_grp_supsvd)
    d_sifa<- dist(Score_test_sifa[,1:i], method = "euclidean")
    hc_sifa <- hclust(d_sifa, method = "ward.D2")
    sub_grp_sifa <- cutree(hc_sifa, k = 5)
  ARI[['SIFA']][i]=mclust::adjustedRandIndex(Subtypes_test, sub_grp_sifa)
}

par(mfrow=c(1,1))
plot(ARI[['COBS']],pch="",
     ylab="Adjusted Rand Index",xlab="Number of Components",
     bty="n",main="Cumulative Components ARI",
     ylim=c(-0.05,0.5)
)
points(ARI[['COBS']],pch=19,type="b",col="red",lwd=3,lty=1)
points(ARI[['SLIDE']],pch=22,type="b",col="blue",lty=2)
points(ARI[['GFA']],pch=24,type="b",col="dark green",lty=3)
points(ARI[['AJIVE']],pch=21,type="b",col="gold3",lty=4)
points(ARI[['SupSVD']],pch=25,type="b",col="darkred",lty=5)
points(ARI[['SIFA']],pch=23,type="b",col="darkorchid",lty=6)
legend("bottomright", legend=c("COBS","SLIDE","GFA","AJIVE","SupSVD","SIFA"),
       col=c("red","blue","dark green","gold3","darkred","darkorchid"), 
       lty=c(1:6),pch=c(19,22,24,21,25,23), cex=0.8)

quartz.save(file="BRCA_ARI.png", type = "png", device = dev.cur(),
             dpi = 300, width=10, height=7)



