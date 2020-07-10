

rm(list = ls()) 

# Source Routines ---------------------------------------------------------

files.sources = list.files(path = "./routines", full.names = TRUE)
sapply(files.sources, source)


# Read the BRCA data ------------------------------------------------------

source("Read_BRCA.R")
dataY1 <-dataY[,1:200] # Gene Expression  
dataY2 <-dataY[,201:400] # Methylation
dataY3 <-dataY[,401:600] # CNV

#alpha.b = 1
#lambda.b = 0.0025
#tune.par.b = list(alpha = alpha.b, lambda = lambda.b)

##### Fixed Parameter #############
BRCA_out1 <- cobs(dataY = dataY, 
             dataX = dataX,  
             Y.bk.idx = gp, 
             X.bk.idx = gpb,
             n.comp = 1,
             alpha.v = 1,
             lambda.v = 0.3)
summary.cobs(BRCA_out1, Y = dataY, X = dataX)

###### Tune ###########

# Proc time data process
# dataX_original <- dataX
# dataY <- cbind(Expr_scarf,Meth_scarf,CNV_scarf)
# gp <- c(rep(1,5125),rep(2,5036),rep(3,6115))
# 
# set.seed(12345)
# Basal_train <- sample(which(Pam50=="Basal"),ceiling(133*0.1),replace=F)
# Her2_train <- sample(which(Pam50=="Her2"),ceiling(51*0.1),replace=F)
# LumA_train <- sample(which(Pam50=="LumA"),ceiling(394*0.1),replace=F)
# LumB_train <- sample(which(Pam50=="LumB"),ceiling(161*0.1),replace=F)
# Normal_train <- sample(which(Pam50=="Normal"),ceiling(31*0.1),replace=F)
# train_samples <- sort(c(Basal_train,Her2_train,LumA_train,LumB_train,Normal_train))
# dataY <- scale(dataY[train_samples,],center = T,scale = F)
# dataX <- dataX_original[train_samples,]
#

#ptm <- proc.time()
BRCA_out <- cobs.tune(dataY = dataY, 
                 dataX = dataX,  
                 Y.bk.idx = gp, 
                 X.bk.idx = gpb,
                 n.comp = 2,  # n.comp=62 to calculate the variation explained table
                 alpha.b = 1,
                 lambda.b = 0.0025,
                 alpha.v = 0,
                 method = "BIC",
                 M.lambda.v = 20,
                 orth=F)
#proc.time()-ptm
# Summary of cobs factorization 
#a<-
  summary.cobs(BRCA_out,X=dataX,Y=dataY)

# Summary 
# 1. number of components, blocks and groups 
# 2. Block-wise and group-wise structure
# 3. variation due to covariate and from unknown source 

 cobs_V=BRCA_out$V
 cobs_B=BRCA_out$B
# cobs_F=BRCA_out$F.est
  
#cobs_V=BRCA_lasso$V
#cobs_B=BRCA_lasso$B
#cobs_F=BRCA_lasso$F.est

# Use the multi-layer estimation for sigma_e^2 and diag Sigma_F
comp_n<-20
se_sq=as.numeric(seEST.full(dataX=dataX,dataY=dataY,B=cobs_B[,1:comp_n], V=cobs_V[,1:comp_n])$sehat^2)
sigma_k_est=0
for (i in 1:comp_n){
  sigma_k_est[i]=(sfEST.full(dataX=dataX,dataY=dataY,b=as.matrix(cobs_B[,i]),v=as.matrix(cobs_V[,i]),se=sqrt(se_sq))$sfhat)^2
}
sf_sq=sigma_k_est
covf=diag(x=sf_sq)

brca_cov_est <- t(as.matrix(cobs_B[,1:comp_n]))%*%cov(dataX)%*%as.matrix(cobs_B[,1:comp_n])+covf
brca_cov_est
diag(brca_cov_est) # in descending order?

##### Score Plots ###########################################
Score=dataY%*%cobs_V
Subtypes=Pam50$V1
library(selectiveInference)

# Use NEW METHOD results above for layer 1 CI
dataY1_CI=dataY
v1_CI=as.matrix(cobs_V[,1])
X1_CI=dataX  #update data for CI
Y1_CI=dataY1_CI%*%v1_CI #update data for CI
data1_CI=list(x=X1_CI,y=Y1_CI)
gfit1=glmnet(x=X1_CI, y=Y1_CI, standardize=FALSE, intercept=FALSE)  
beta1=coef(gfit1, x=X1_CI, y=Y1_CI, s=0.0025/n, exact=TRUE)[-1] #glmnet result
sfit1=SGL(data=data1_CI, index=gpb, type="linear",nlam=1,alpha=1,lambdas=0.0025,standardize=F)  # No intercept
beta1_SGL=sfit1$beta; sfit1$intercept#SGL result, no intercept
beta1; beta1_SGL; cobs_B[,1] #compare glmnet, sgl and COBS original results
Coeff1_CI=fixedLassoInf(x=X1_CI,y=Y1_CI,beta1,lambda=0.0025,family="gaussian",alpha=0.05,intercept = F,type="full") #used glmnet result
Coeff1_CI$ci

B1_CI <- data.frame(Subtypes=levels(Subtypes), B_hat =cobs_B[,1],B_refit=beta1,Lower=c(Coeff1_CI$ci[,1]),Upper =c(Coeff1_CI$ci[,2]))

# Use NEW METHOD results above for layer 2 CI
Temp_est=BRCA_out$est[[1]]
dataY2_CI=dataY-(Temp_est$s[2]*dataY%*% as.matrix(cobs_V[,1])+Temp_est$s[1]*dataX%*%as.matrix(cobs_B[,1]))%*%t(cobs_V[,1])/(sum(Temp_est$s)) 
X2_CI=dataX
v2_CI=as.matrix(cobs_V[,2])
Y2_CI=dataY2_CI%*%v2_CI
data2_CI=list(x=X2_CI,y=Y2_CI) # for SGL calculation
gfit2=glmnet(x=X2_CI, y=Y2_CI, standardize=FALSE, intercept=FALSE,family="gaussian",alpha=1,lambda=0.0025)
beta2=coef(gfit2, x=X2_CI, y=Y2_CI, s=0.0025/n, exact=TRUE,alpha=1)[-1] #glmnet result
sfit2=SGL(data=data2_CI, index=gpb, type="linear",nlam=1,alpha=1,lambdas=0.0025,standardize=F)
beta2_SGL=sfit2$beta; sfit2$intercept #SGL result, no intercept
beta2; beta2_SGL; cobs_B[,2] #compare glmnet, sgl and COBS original results
Coeff2_CI=fixedLassoInf(x=X2_CI,y=Y2_CI,beta2,lambda=0.0025,family="gaussian",alpha=0.05,intercept = F,type="full")
Coeff2_CI$ci

B2_CI <- data.frame(Subtypes=levels(Subtypes), B_hat =cobs_B[,2],B_refit=beta2,Lower=c(Coeff2_CI$ci[,1]),Upper =c(Coeff2_CI$ci[,2]))


# Plot PC1 PC2 Score + CI of B
Score_BRCA=data.frame(Score[,1:2])
colnames(Score_BRCA)=c("Score1", "Score2")
Score_BRCA=cbind(Score_BRCA,Subtypes)
colors <- Pam50$V1
library(ggpubr)
ScoreScatter <- ggscatter(Score_BRCA, x="Score1", y="Score2",color="Subtypes",shape="Subtypes",xlim=c(-15,25),ylim=c(-15,40)) + border()
# Marginal CI plot of b1 (top panel) and b2 (right panel)
xplot <- ggplot(B1_CI, aes(x=B_hat,y=Subtypes))+
  geom_point(aes(shape=Subtypes,color=Subtypes),size=3)+
  #geom_point(aes(x=B_refit),size=1)+
  geom_errorbarh(aes(xmax =Upper, xmin =Lower))+xlim(-15,25)+theme_gray()
yplot <- ggplot(B2_CI, aes(x=Subtypes,y=B_hat))+
  geom_point(aes(shape=Subtypes,color=Subtypes),size=3)+
  #geom_point(aes(y=B_refit),size=1)+
  geom_errorbar(aes(ymax =Upper, ymin =Lower))+ylim(-15,40)+theme_gray()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
corner <- ggplot() + clean_theme()
# Arranging the plot
#BRCA_Score=
ggarrange(xplot+theme(plot.margin=margin(b=.1,r=.1)), corner+theme(plot.margin=margin(l=.1,b=.1)), ScoreScatter+theme(plot.margin=margin(t=.1,r=.1)), yplot+theme(plot.margin=margin(l=.1,t=.1)), 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(3, 1), heights = c(1, 3),
          common.legend = TRUE)

#quartz.save(file="COBS_BRCA_Score.png", type = "png", device = dev.cur(),
#            dpi = 300, width=10, height=10)


#### Variance Explained #####
# Overall Variance explained proportion  
SSTO=norm(dataY,"F")^2
SSR_Supervision=norm(dataX%*%cobs_B%*%t(cobs_V),"F")^2
SSR_Unknown=norm(cobs_F%*%t(cobs_V),"F")^2
SSE=SSTO-SSR_Supervision-SSR_Unknown

SSR_Supervision/SSTO # [1] 0.1566764
SSR_Unknown/SSTO # [1] 0.5598374
SSE/SSTO # [1] 0.2834862

# Gene Variance explained proportion ###################################
Y1=dataY[,1:200]
SST_gene=norm(Y1,"F")^2
SSX_gene=norm(dataX%*%cobs_B%*%t(cobs_V[1:200,]),"F")^2
SSF_gene=norm(cobs_F%*%t(cobs_V[1:200,]),"F")^2
SSE_gene=SST_gene-SSX_gene-SSF_gene

SSX_gene/SST_gene # [1] 0.2684128
SSF_gene/SST_gene # [1] 0.4198128
SSE_gene/SST_gene # [1] 0.3117744

# Among Gene supervision, check each subtype
SS_gene_B=norm(as.matrix(dataX[,1])%*%cobs_B[1,]%*%t(cobs_V[1:200,]),"F")^2
SS_gene_H=norm(as.matrix(dataX[,2])%*%cobs_B[2,]%*%t(cobs_V[1:200,]),"F")^2
SS_gene_L1=norm(as.matrix(dataX[,3])%*%cobs_B[3,]%*%t(cobs_V[1:200,]),"F")^2
SS_gene_L2=norm(as.matrix(dataX[,4])%*%cobs_B[4,]%*%t(cobs_V[1:200,]),"F")^2
SS_gene_N=norm(as.matrix(dataX[,5])%*%cobs_B[5,]%*%t(cobs_V[1:200,]),"F")^2

#SS_gene_B+SS_gene_H+SS_gene_L1+SS_gene_L2+SS_gene_N -SSX_gene <10^-8

SS_gene_B/SSX_gene # [1] 0.6043465
SS_gene_H/SSX_gene # [1] 0.07777413
SS_gene_L1/SSX_gene # [1] 0.1873801
SS_gene_L2/SSX_gene # [1] 0.09800522
SS_gene_N/SSX_gene # [1] 0.03249411

# Methylation Variance explained proportion #########################
Y2=dataY[,201:400]
SST_methy=norm(Y2,"F")^2
SSX_methy=norm(dataX%*%cobs_B%*%t(cobs_V[201:400,]),"F")^2
SSF_methy=norm(cobs_F%*%t(cobs_V[201:400,]),"F")^2
SSE_methy=SST_methy-SSX_methy-SSF_methy

SSX_methy/SST_methy # [1] 0.1219317
SSF_methy/SST_methy # [1] 0.4666488
SSE_methy/SST_methy # [1] 0.4114195

# Among Methy supervision, check each subtype
SS_methy_B=norm(as.matrix(dataX[,1])%*%cobs_B[1,]%*%t(cobs_V[201:400,]),"F")^2
SS_methy_H=norm(as.matrix(dataX[,2])%*%cobs_B[2,]%*%t(cobs_V[201:400,]),"F")^2
SS_methy_L1=norm(as.matrix(dataX[,3])%*%cobs_B[3,]%*%t(cobs_V[201:400,]),"F")^2
SS_methy_L2=norm(as.matrix(dataX[,4])%*%cobs_B[4,]%*%t(cobs_V[201:400,]),"F")^2
SS_methy_N=norm(as.matrix(dataX[,5])%*%cobs_B[5,]%*%t(cobs_V[201:400,]),"F")^2

#SS_methy_B+SS_methy_H+SS_methy_L1+SS_methy_L2+SS_methy_N -SSX_methy <10^-8

SS_methy_B/SSX_methy # [1] 0.5559726
SS_methy_H/SSX_methy  # [1] 0.02977273
SS_methy_L1/SSX_methy  # [1] 0.1066253
SS_methy_L2/SSX_methy # [1] 0.201131
SS_methy_N/SSX_methy # [1] 0.1064984

# CNV Variance explained proportion #################################
Y3=dataY[,401:600]
SST_cnv=norm(Y3,"F")^2
SSX_cnv=norm(dataX%*%cobs_B%*%t(cobs_V[401:600,]),"F")^2
SSF_cnv=norm(cobs_F%*%t(cobs_V[401:600,]),"F")^2
SSE_cnv=SST_cnv-SSX_cnv-SSF_cnv

SSX_cnv/SST_cnv # [1] 0.07968474
SSF_cnv/SST_cnv # [1] 0.7930505
SSE_cnv/SST_cnv # [1] 0.1272647

# Among Gene supervision, check each subtype
SS_cnv_B=norm(as.matrix(dataX[,1])%*%cobs_B[1,]%*%t(cobs_V[401:600,]),"F")^2
SS_cnv_H=norm(as.matrix(dataX[,2])%*%cobs_B[2,]%*%t(cobs_V[401:600,]),"F")^2
SS_cnv_L1=norm(as.matrix(dataX[,3])%*%cobs_B[3,]%*%t(cobs_V[401:600,]),"F")^2
SS_cnv_L2=norm(as.matrix(dataX[,4])%*%cobs_B[4,]%*%t(cobs_V[401:600,]),"F")^2
SS_cnv_N=norm(as.matrix(dataX[,5])%*%cobs_B[5,]%*%t(cobs_V[401:600,]),"F")^2

#SS_cnv_B+SS_cnv_H+SS_cnv_L1+SS_cnv_L2+SS_cnv_N -SSX_cnv <10^-8

SS_cnv_B/SSX_cnv # [1] 0.1728918
SS_cnv_H/SSX_cnv  # [1] 0.3033468
SS_cnv_L1/SSX_cnv  # [1] 0.1683631
SS_cnv_L2/SSX_cnv # [1] 0.2972557
SS_cnv_N/SSX_cnv # [1] 0.05814251

# Correlation Network, alpha_v=0 #########################
name_expr=colnames(dataY1)
name_meth=colnames(dataY2)
name_cnv=colnames(dataY3)

library("qgraph")
Gene_PC1=cobs_V[1:200,1]; names(Gene_PC1)=name_expr
Methy_PC1=cobs_V[201:400,1]; names(Methy_PC1)=name_meth
name_gene_PC1=names(Gene_PC1[rank(-abs(Gene_PC1))<=10])
name_methy_PC1=names(Methy_PC1[rank(-abs(Methy_PC1))<=10])
#Gene_PC1[name_gene_PC1]
#Methy_PC1[name_methy_PC1]
#Corr_top10PC1_GM=rcorr(Y1[,name_gene_PC1],Y2[,name_methy_PC1])$r

CNV_PC2=cobs_V[401:600,2]; names(CNV_PC2)=name_cnv
name_cnv_PC2=names(CNV_PC2[rank(-abs(CNV_PC2))<=10])
Corr_allthree=Hmisc::rcorr(cbind(dataY1[,name_gene_PC1],dataY2[,name_methy_PC1],dataY3[,name_cnv_PC2]))$r
Groups <- c(rep("Gene Expression",10),rep("Methylation",10),rep("CNV",10))
shapes <- c(rep("square",10),rep("circle",10),rep("triangle",10))

par(mar=c(4,4,4,4),mgp=c(2,0.8,0))
par(mfrow = c(1, 1))
Graph_cor <- qgraph(Corr_allthree,graph ="cor", groups=Groups, vTrans=200,label.prop=0.8,vsize=5,
                    borders=F,layout="spring",theme='TeamFortress',palette='pastel', minimum=0.1,
                    aspect=1.5, GLratio=2,
                    shape=shapes,legend=F)
legend(x=-0.9,y=-0.6, legend=c("Gene expression","Methylation","CNV"),horiz=T,
       pch=c(15:17), col=c("palegreen3","skyblue3","lightpink3"),box.col="white")
title("Expression - Methylation - CNV Correlations",line=-4)
# quartz.save(file="Corr_Network.png", type = "png", device = dev.cur(),
#             dpi = 300, width=7, height=7)

# Different choice of alpha_v #########################
BRCA_lasso <- cobs.tune(dataY = dataY, 
                      dataX = dataX,  
                      Y.bk.idx = gp, 
                      X.bk.idx = gpb,
                      n.comp = 5,
                      alpha.b = 1,
                      lambda.b = 0.0025,
                      alpha.v = 1,
                      method = "BIC",
                      M.lambda.v = 20,
                      orth=F)

BRCA_sgl <- cobs.tune(dataY = dataY, 
                        dataX = dataX,  
                        Y.bk.idx = gp, 
                        X.bk.idx = gpb,
                        n.comp = 5,
                        alpha.b = 1,
                        lambda.b = 0.0025,
                        alpha.v = 0.5,
                        method = "BIC",
                        M.lambda.v = 20,
                        orth=F)


# Summary of cobs factorization 
summary.cobs(BRCA_lasso); summary.cobs(BRCA_sgl)

