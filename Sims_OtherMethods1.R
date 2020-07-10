# 100 Simulation Comparisons from SVD, AJIVE, SupSVD, and SIFA #
#rm(list = ls()) 
routines = list.files(path = './routines', full.names = TRUE)
sapply(routines, source)


# Data Process ---------------------------------------------------------
# Import dataY, dataX, true_B, and true_V first;
# Read and separate data from enviroment
# Input n=sample size

# Large n, Small p
load("./SimulatedData/Sim_COBS.Rdata") # In this rdata, still called SCARF so far.
load("./SimulatedData/Sim_AJIVE.Rdata")
load("./SimulatedData/Sim_SIFA.Rdata")

# Small n, Large p
load("./SimulatedData/Sim_COBS_small_n.Rdata")
load("./SimulatedData/Sim_AJIVE_small_n.Rdata")
load("./SimulatedData/Sim_SIFA_small_n.Rdata")


#Change Data according to different settings, in the COBS data, the name is still SCARF so far.
Sim_dataX=SCARF_LSF_X #LSF = large SF, SSF= Small SF
Sim_dataY=SCARF_LSF_Y
true_B=SCARF_True_B
true_V=SCARF_True_V

# Setting of the groups, according to true_B and true_V
n=dim(Sim_dataY[[1]])[1]
p=dim(Sim_dataY[[1]])[2]
q=dim(Sim_dataX[[1]])[2]
gv=4; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b
gpb=ceiling(1:q / size.groupb)
gp=ceiling(1:p/ size.groupv)

##### SVD #####
PA_SVD=rep(NA,100) #Principal Angle
Grass_SVD=rep(NA,100) #Grassmannian Distance
Frob_SVD=rep(NA,100) #Frobenius Norm
loss_SVD=rep(NA,100) #Loss Function
sim_svd=list()  #Store results
rep_n=1
while(rep_n<=100){
  dataY=as.matrix(Sim_dataY[[rep_n]])
  sim_svd[[rep_n]]=svd(dataY,nv=4,nu=4)
  for (j in 1:4){
    if (angle(sim_svd[[rep_n]]$v[,j],true_V[,j])>90){
      sim_svd[[rep_n]]$v[,j]=-sim_svd[[rep_n]]$v[,j]
    }
  }
  PA_SVD[rep_n]=PrinAngle(true_V,sim_svd[[rep_n]]$v)$Max_angle
  Grass_SVD[rep_n]=GrassDist(true_V,sim_svd[[rep_n]]$v)$Grass_Dist
  Frob_SVD[rep_n]=norm(as.matrix(true_V-sim_svd[[rep_n]]$v),"F")^2
  loss_SVD[rep_n]=lossV(M_true=true_V,M_est=sim_svd[[rep_n]]$v,group=4)$Min_loss
  rep_n=rep_n+1
}

#Result Summary
summary_svd=data.frame(matrix(NA,nrow=4,ncol=2),
                         row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm"))
colnames(summary_svd)=c("mean","sd")
summary_svd[1,1]=round(mean(PA_SVD,na.rm=T),digits = 2)
summary_svd[1,2]=round(sd(PA_SVD,na.rm=T),digits=2)
summary_svd[2,1]=round(mean(Grass_SVD,na.rm=T),digits = 2)
summary_svd[2,2]=round(sd(Grass_SVD,na.rm=T),digits=2)
summary_svd[3,1]=round(mean(loss_SVD,na.rm=T),digits = 2)
summary_svd[3,2]=round(sd(loss_SVD,na.rm=T),digits=2)
summary_svd[4,1]=round(mean(Frob_SVD,na.rm=T),digits = 2)
summary_svd[4,2]=round(sd(Frob_SVD,na.rm=T),digits=2)

#Save Results
#save small n result (Capital words -- Settings, Lower case -- Method)
write.table(summary_svd,file="./Sim_Results_Sn/Summary_AJIVE_S_svd.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_svd,file="./Sim_Results_New/Summary_svd_CBOS_L.csv",sep=" ",col.names=T,row.names=T)
#save(sim_svd,file="sim_svd_AJIVE_S.Rdata")

#Plot estimates of first data set vs. true para
par(mar=c(4,4,4,4),mgp=c(2,0.8,0),oma=c(0,0,2,0))
par(mfrow = c(4, 1))
for ( i in 1:4){
  plot(true_V[,i],ylab=bquote("V"~.(i)),ylim=c(-0.8,0.8),type="l",col="red")
  points(sim_svd[[1]]$v[,i],type="l",lty=2,col="dark green")
  abline(v=cumsum(table(gp)),lty=2,col="gray")
}
title("SVD Estimated V vs. True Settings",outer=T)
#quartz.save(file="SIFA_L_svd.png", type = "png", device = dev.cur(),
#            dpi = 150, width=10, height=10)


##### AJIVE #####
# install.packages('devtools')
devtools::install_github("idc9/r_jive")
library(ajive)
library(cowplot)

PA_AJIVE=rep(NA,100) #Principal Angle
Grass_AJIVE=rep(NA,100) #Grassmannian Distance
Frob_AJIVE=rep(NA,100) #Frobenius Norm
loss_AJIVE=rep(NA,100) #Loss Function
sim_ajive=list()  #Store results
rep_n=1
while(rep_n <= 100){
  dataY=as.matrix(Sim_dataY[[rep_n]])
  data_blocks <- list(dataY[,1:size.groupv],dataY[,(size.groupv+1):(size.groupv*2)],
                      dataY[,(size.groupv*2+1):(size.groupv*3)],dataY[,(size.groupv*3+1):(size.groupv*4)])
  initial_signal_ranks <- c(2,2,2,1) #Scarf(2,2,2,2) #SIFA,AJIVE(2,2,2,1)
  jive_decomp <- ajive(data_blocks, initial_signal_ranks,full=F,joint_rank=NA)
  I1=as.matrix(jive_decomp$block_decomps[[1]][['individual']][['v']])
  I2=as.matrix(jive_decomp$block_decomps[[2]][['individual']][['v']])
  I3=as.matrix(jive_decomp$block_decomps[[3]][['individual']][['v']])
  I4=as.matrix(jive_decomp$block_decomps[[4]][['individual']][['v']])
  J1=as.matrix(jive_decomp$block_decomps[[1]][['joint']][['v']])
  J2=as.matrix(jive_decomp$block_decomps[[2]][['joint']][['v']])
  J3=as.matrix(jive_decomp$block_decomps[[3]][['joint']][['v']])
  J4=as.matrix(jive_decomp$block_decomps[[4]][['joint']][['v']])
  
  J_AJIVE=rbind(J1,J2,J3,J4)
  if(sum(colSums(J_AJIVE)!=0)==length(J_AJIVE[1,])){J_AJIVE=J_AJIVE/apply(J_AJIVE,2,norm.fcn)}
  if(sum(colSums(I1)!=0)==length(I1[1,])){I1=I1/apply(I1,2,norm.fcn)}
  if(sum(colSums(I2)!=0)==length(I2[1,])){I2=I2/apply(I2,2,norm.fcn)}
  if(sum(colSums(I3)!=0)==length(I3[1,])){I3=I3/apply(I3,2,norm.fcn)}
  if(sum(colSums(I4)!=0)==length(I4[1,])){I4=I4/apply(I4,2,norm.fcn)}
  
  AJIVE_V=as.matrix(cbind(J_AJIVE,blkdiag(I1,I2,I3,I4)))  # Use this part for SIFA and JIVE setting
  #AJIVE_V=as.matrix(cbind(blkdiag(I1,I2,I3,I4),J_AJIVE)) # Use this part for SCARF setting
  AJIVE_V=AJIVE_V[,apply(AJIVE_V,2,sum)!=0]
  AJIVE_V=AJIVE_V/apply(AJIVE_V,2,norm.fcn)
  
  sim_ajive[[rep_n]]=AJIVE_V
  
  PA_AJIVE[rep_n]=PrinAngle(true_V,AJIVE_V)$Max_angle
  Grass_AJIVE[rep_n]=GrassDist(true_V,AJIVE_V)$Grass_Dist
  loss_AJIVE[rep_n]=lossV(M_true=true_V,M_est=AJIVE_V,group=4)$Min_loss
  
  if(length(AJIVE_V[1,])==length(true_V[1,])){ # Adjust this according to the results w.r.t. different initial ranks
    Frob_AJIVE[rep_n]=norm(as.matrix(true_V-AJIVE_V),"F")^2
  } else {Frob_AJIVE[rep_n]=NA}
  cat(rep_n,"simulation done.")
  rep_n=rep_n+1
}

#Result Summary
summary_ajive=data.frame(matrix(NA,nrow=4,ncol=2),
                       row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm"))
colnames(summary_ajive)=c("mean","sd")
summary_ajive[1,1]=round(mean(PA_AJIVE,na.rm=T),digits = 2)
summary_ajive[1,2]=round(sd(PA_AJIVE,na.rm=T),digits=2)
summary_ajive[2,1]=round(mean(Grass_AJIVE,na.rm=T),digits = 2)
summary_ajive[2,2]=round(sd(Grass_AJIVE,na.rm=T),digits=2)
summary_ajive[3,1]=round(mean(loss_AJIVE,na.rm=T),digits = 2)
summary_ajive[3,2]=round(sd(loss_AJIVE,na.rm=T),digits=2)
summary_ajive[4,1]=round(mean(Frob_AJIVE,na.rm=T),digits = 2)
summary_ajive[4,2]=round(sd(Frob_AJIVE,na.rm=T),digits=2)

#Save Results
#save small n result (Capital words -- Settings, Lower case -- Method)
write.table(summary_ajive,file="./Sim_Results_Sn/Summary_AJIVE_S_ajive.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_ajive,file="./Sim_Results_New/Summary_ajive_CBOS_L.csv",sep=" ",col.names=T,row.names=T)
#save(sim_ajive,file="sim_ajive_AJIVE_S.Rdata")

#Plot estimates of first data set vs. true para
par(mar=c(4,4,4,4),mgp=c(2,0.8,0),oma=c(0,0,2,0))
par(mfrow = c(4, 1))
for ( i in 1:4){
  plot(true_V[,i],ylab=bquote("V"~.(i)),ylim=c(-0.8,0.8),type="l",col="red")
  points(sim_ajive[[1]][,i],type="l",lty=2,col="dark green")
  abline(v=cumsum(table(gp)),lty=2,col="gray")
}
title("AJIVE Estimated V vs. True Settings",outer=T)
#quartz.save(file="SIFA_L_ajive.png", type = "png", device = dev.cur(),
#            dpi = 150, width=10, height=10)


##### SupSVD #####
#library(RSpectra)
#library(fBasics)
#library(pracma)
#library(SuperPCA)

PA_supsvd_V=rep(NA,100)
Grass_supsvd_V=rep(NA,100)
loss_supsvd_V=rep(NA,100)
Frob_supsvd_V=rep(NA,100)
Frob_supsvd_B=rep(NA,100)
sim_supsvd=list()
rep_n=1
while(rep_n<=100){
  dataY=as.matrix(Sim_dataY[[rep_n]])
  dataX=as.matrix(Sim_dataX[[rep_n]])
  sim_supsvd[[rep_n]]=SuperPCA::SupPCA(Y=dataX, X=dataY, r=4)
  for (j in 1:4){
    if (angle(sim_supsvd[[rep_n]]$V[,j],true_V[,j])>90) {
      sim_supsvd[[rep_n]]$V[,j]=-sim_supsvd[[rep_n]]$V[,j]
      sim_supsvd[[rep_n]]$B[,j]=-sim_supsvd[[rep_n]]$B[,j]}
  }
  PA_supsvd_V[rep_n]=PrinAngle(true_V,sim_supsvd[[rep_n]]$V)$Max_angle
  Grass_supsvd_V[rep_n]=GrassDist(true_V,sim_supsvd[[rep_n]]$V)$Grass_Dist
  Frob_supsvd_V[rep_n]=norm(as.matrix(true_V-sim_supsvd[[rep_n]]$V),"F")^2
  Frob_supsvd_B[rep_n]=norm(as.matrix(true_B-sim_supsvd[[rep_n]]$B),"F")^2
  loss_supsvd_V[rep_n]=lossV(M_true=true_V,M_est=sim_supsvd[[rep_n]]$V,group=4)$Min_loss
  rep_n=rep_n+1
}

#Result Summary
summary_supsvd=data.frame(matrix(NA,nrow=5,ncol=2),
                         row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm","BFrobNorm"))
colnames(summary_supsvd)=c("mean","sd")
summary_supsvd[1,1]=round(mean(PA_supsvd_V,na.rm=T),digits = 2)
summary_supsvd[1,2]=round(sd(PA_supsvd_V,na.rm=T),digits=2)
summary_supsvd[2,1]=round(mean(Grass_supsvd_V,na.rm=T),digits = 2)
summary_supsvd[2,2]=round(sd(Grass_supsvd_V,na.rm=T),digits=2)
summary_supsvd[3,1]=round(mean(loss_supsvd_V,na.rm=T),digits = 2)
summary_supsvd[3,2]=round(sd(loss_supsvd_V,na.rm=T),digits=2)
summary_supsvd[4,1]=round(mean(Frob_supsvd_V,na.rm=T),digits = 2)
summary_supsvd[4,2]=round(sd(Frob_supsvd_V,na.rm=T),digits=2)
summary_supsvd[5,1]=round(mean(Frob_supsvd_B,na.rm=T),digits = 2)
summary_supsvd[5,2]=round(sd(Frob_supsvd_B,na.rm=T),digits=2)

#Save Results
#save small n result (Capital words -- Settings, Lower case -- Method)
## write.table(summary_supsvd,file="./Sim_Results_Sn/Summary_AJIVE_S_supsvd.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_supsvd,file="./Sim_Results_New/Summary_supsvd_CBOS_L.csv",sep=" ",col.names=T,row.names=T)
#save(sim_supsvd,file="sim_supsvd_AJIVE_S.Rdata")

#Plot estimates of first data set vs. true para
par(mar=c(4,4,4,4),mgp=c(2,0.8,0),oma=c(0,0,2,0))
par(mfrow = c(4, 2))
for ( i in 1:4){
  plot(true_B[,i],ylab=bquote("B"~.(i)),ylim=c(-8,8),type="l",col="gold")
  points(sim_supsvd[[1]]$B[,i],type="l",lty=2,col="dark green")
  abline(v=cumsum(table(gpb)),lty=2,col="gray")
  plot(true_V[,i],ylab=bquote("V"~.(i)),ylim=c(-0.8,0.8),type="l",col="red")
  points(sim_supsvd[[1]]$V[,i],type="l",lty=2,col="dark green")
  abline(v=cumsum(table(gp)),lty=2,col="gray")
}
title("SupSVD Estimated B and V vs. True Settings",outer=T)
#quartz.save(file="SCARF_S_SupSVD.png", type = "png", device = dev.cur(),
#            dpi = 150, width=10, height=10)

##### SIFA - Matlab #####

# For results input from Matlab
# Run Matlab code and get estimates from matlab "Sim_SIFA_SCARF.m"
# read 100 results V and B from csv file; rename as SIFA_V_list and SIFA_B_list

#true_B=AJIVE_True_B
#true_V=AJIVE_True_V

SIFA_V=list()
SIFA_B=list()
PA_SIFA_V=rep(NA,100)
Grass_SIFA_V=rep(NA,100)
loss_SIFA_V=rep(NA,100)
Frob_SIFA_V=rep(NA,100)
Frob_SIFA_B=rep(NA,100)
for (i in 1:100){
  SIFA_V[[i]]=as.matrix(SIFA_V_list_AJIVE_S_Sn[((i-1)*p+1):(i*p),]); # p is the row number of V
  SIFA_B[[i]]=as.matrix(SIFA_B_list_AJIVE_S_Sn[((i-1)*q+1):(i*q),]); # q is the row number of B

  for( j in 1:4){
    if(angle(SIFA_V[[i]][,j],true_V[,j])>90){
      SIFA_V[[i]][,j] = -SIFA_V[[i]][,j]
      SIFA_B[[i]][,j] = -SIFA_B[[i]][,j]
    }
  }

  PA_SIFA_V[i]=PrinAngle(true_V,SIFA_V[[i]])$Max_angle
  Grass_SIFA_V[i]=GrassDist(true_V,SIFA_V[[i]])$Grass_Dist
  loss_SIFA_V[i]=lossV(true_V,SIFA_V[[i]],group = 4)$Min_loss
  Frob_SIFA_V[i]=norm(as.matrix(true_V-SIFA_V[[i]]),"F")^2
  Frob_SIFA_B[i]=norm(as.matrix(true_B-SIFA_B[[i]]),"F")^2
}


##### SIFA - R package #####
# Use the SuperPCA package, but it seems there is some problem:
  # In SIFA code, S2 <- RSpectra::svds(Ycurrent, r[1])
                # U_k_ini <- S2$u
                # D_k_ini <- diag(S2$d)
                # V_k_ini <- S2$v
                # U_k <- U_k_ini %*% D_k_ini
  # If r[1]=1, it will show an error message: Error in U_k_ini %*% D_k_ini : non-conformable arguments
# PA_SIFA_V=rep(NA,100)
# Grass_SIFA_V=rep(NA,100)
# loss_SIFA_V=rep(NA,100)
# Frob_SIFA_V=rep(NA,100)
# Frob_SIFA_B=rep(NA,100)
# sim_sifa=list()
# rep_n=1
# while(rep_n<=100){
#   dataY=as.matrix(Sim_dataY[[rep_n]])
#   dataX=as.matrix(Sim_dataX[[rep_n]])
#   sifa_Y <- list(dataY[,1:size.groupv],dataY[,(size.groupv+1):(size.groupv*2)],
#                  dataY[,(size.groupv*2+1):(size.groupv*3)],dataY[,(size.groupv*3+1):(size.groupv*4)])
#   # Under COBS setting:
#   fit_sifa <- SuperPCA::SIFA(X=dataX,Y=sifa_Y,r0=2,r=c(1,1,0,0),sparsity=1,type="A")
#   SIFA_B <- cbind(B[[1]],B[[2]],B0)
#   SIFA_V <- cbind(pracma::blkdiag(V_ind[[1]],V_ind[[2]],V_ind[[3]],V_ind[[4]]),V_joint)
# 
#   # Under AJIVE and SIFA setting:
#   # fit_sifa <- SuperPCA::SIFA(X=dataX,Y=sifa_Y,r0=1,r=c(1,1,1,0),sparsity=1,type="A")
#   # SIFA_B=cbind(B0,B[[1]],B[[2]],B[[3]]);
#   # SIFA_V=cbind(V_joint, pracma::blkdiag(V_ind{1},V_ind{2},V_ind{3},V_ind{4}));
# 
#   sime_sifa[[rep_n]] <- list(V=SIFA_V, B=SIFA_B)
#   for (j in 1:4){
#     if (angle(sim_sifa[[rep_n]]$V[,j],true_V[,j])>90) {
#       sim_sifa[[rep_n]]$V[,j]=-sim_sifa[[rep_n]]$V[,j]
#       sim_sifa[[rep_n]]$B[,j]=-sim_sifa[[rep_n]]$B[,j]}
#   }
# 
#   PA_SIFA_V[rep_n]=PrinAngle(true_V,sim_sifa[[rep_n]]$V)$Max_angle
#   Grass_SIFA_V[rep_n]=GrassDist(true_V,sim_sifa[[rep_n]]$V)$Grass_Dist
#   Frob_SIFA_V[rep_n]=norm(as.matrix(true_V-sim_sifa[[rep_n]]$V),"F")^2
#   Frob_SIFA_B[rep_n]=norm(as.matrix(true_B-sim_sifa[[rep_n]]$B),"F")^2
#   loss_SIFA_V[rep_n]=lossV(M_true=true_V,M_est=sim_sifa[[rep_n]]$V,group=4)$Min_loss
#   rep_n=rep_n+1
# }

#Result Summary
summary_sifa=data.frame(matrix(NA,nrow=5,ncol=2),
                          row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm","BFrobNorm"))
colnames(summary_sifa)=c("mean","sd")
summary_sifa[1,1]=round(mean(PA_SIFA_V,na.rm=T),digits = 2)
summary_sifa[1,2]=round(sd(PA_SIFA_V,na.rm=T),digits=2)
summary_sifa[2,1]=round(mean(Grass_SIFA_V,na.rm=T),digits = 2)
summary_sifa[2,2]=round(sd(Grass_SIFA_V,na.rm=T),digits=2)
summary_sifa[3,1]=round(mean(loss_SIFA_V,na.rm=T),digits = 2)
summary_sifa[3,2]=round(sd(loss_SIFA_V,na.rm=T),digits=2)
summary_sifa[4,1]=round(mean(Frob_SIFA_V,na.rm=T),digits = 2)
summary_sifa[4,2]=round(sd(Frob_SIFA_V,na.rm=T),digits=2)
summary_sifa[5,1]=round(mean(Frob_SIFA_B,na.rm=T),digits = 2)
summary_sifa[5,2]=round(sd(Frob_SIFA_B,na.rm=T),digits=2)

#Save Results
#save small n result (Capital words -- Settings, Lower case -- Method)
## write.table(summary_sifa,file="./Sim_Results_Sn/Summary_AJIVE_S_sifa.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_sifa,file="./Sim_Results_New/Summary_sifa_CBOS_L.csv",sep=" ",col.names=T,row.names=T)

#Plot estimates of first data set vs. true para
par(mar=c(4,4,4,4),mgp=c(2,0.8,0),oma=c(0,0,2,0))
par(mfrow = c(4, 2))
for ( i in 1:4){
  plot(true_B[,i],ylab=bquote("B"~.(i)),ylim=c(-8,8),type="l",col="gold")
  points(SIFA_B[[1]][,i],type="l",lty=2,col="dark green")
  abline(v=cumsum(table(gpb)),lty=2,col="gray")
  plot(true_V[,i],ylab=bquote("V"~.(i)),ylim=c(-0.8,0.8),type="l",col="red")
  points(SIFA_V[[1]][,i],type="l",lty=2,col="dark green")
  abline(v=cumsum(table(gp)),lty=2,col="gray")
}
title("SIFA Estimated B and V vs. True Settings",outer=T)
#quartz.save(file="SIFA_S_sifa.png", type = "png", device = dev.cur(),
#            dpi = 150, width=10, height=10)


##### COBS Old Method, not used here, use Sims_COBS.R for COBS simulation results #####
# Old Method # 
# --------------------------- Use NewTune_experiment for SCARF #
# time_track=rep(NA,100) #tracking user time for each run
# PA_SCARF_V=rep(NA,100) #Principal Angle
# Grass_SCARF_V=rep(NA,100) #Grassmannian Distance
# Frob_SCARF_V=rep(NA,100) #Frobenius Norm
# Frob_SCARF_B=rep(NA,100) #Frobenius Norm
# loss_SCARF_V=rep(NA,100) #Loss Function
# sim_scarf_tune=list()  #Store results
# rep_n=1
# while(rep_n <= 1){
#   dataY=as.matrix(Sim_dataY[[rep_n]])
#   dataX=as.matrix(Sim_dataX[[rep_n]])
#   Y=dataY
#   X=dataX
#   ptm <- proc.time()
#   #sim_scarf_tune[[rep_n]]=SCARF.tune(Y=dataY,X=dataX,layer=4,Y_index=gp,X_index=gpb,alpha.b=0.2,alpha.v=0.2,lam.b.seq =c(0.15,0.1,0.05,0.01,0.005),lam.v.seq = c(0.27,0.2,0.15,0.07,0.001))
#   #sim_scarf_tune[[rep_n]]=SCARF(Y=dataY,X=dataX,layer=4,Y_index=gp,X_index=gpb,alpha.b=0.2,alpha.v=0.2,lambda.b=0.01,lambda.v=0.07)
#   sim_scarf_tune[[rep_n]]=SCARF.tune(Y=dataY,X=dataX,layer=4,Y_index=gp,X_index=gpb,alpha.b=0.2,alpha.v=0.2,lam.b.seq =0.1,lam.v.seq = 0.1)
#    time_track[rep_n]=(proc.time() - ptm)[1]
#   for (j in 1:4){
#     if (angle(sim_scarf_tune[[rep_n]]$V[,j],true_V[,j])>90) {
#       sim_scarf_tune[[rep_n]]$V[,j]=-sim_scarf_tune[[rep_n]]$V[,j]
#       sim_scarf_tune[[rep_n]]$B[,j]=-sim_scarf_tune[[rep_n]]$B[,j]}
#   }
#   
#   PA_SCARF_V[rep_n]=PrinAngle(true_V,sim_scarf_tune[[rep_n]]$V)$Max_angle
#   Grass_SCARF_V[rep_n]=GrassDist(true_V,sim_scarf_tune[[rep_n]]$V)$Grass_Dist
#   Frob_SCARF_V[rep_n]=norm(as.matrix(true_V-sim_scarf_tune[[rep_n]]$V),"F")
#   Frob_SCARF_B[rep_n]=norm(as.matrix(true_B-sim_scarf_tune[[rep_n]]$B),"F")
#   loss_SCARF_V[rep_n]=lossV(M_true=true_V,M_est=sim_scarf_tune[[rep_n]]$V,group=4)$Min_loss
#   cat(rep_n,"run done.")
#   rep_n=rep_n+1
# }
# 
# #Result Summary
# summary_scarf=data.frame(matrix(NA,nrow=5,ncol=2),
#                            row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm","BFrobNorm"))
# colnames(summary_scarf)=c("mean","sd")
# summary_scarf[1,1]=round(mean(PA_SCARF_V,na.rm=T),digits = 2)
# summary_scarf[1,2]=round(sd(PA_SCARF_V,na.rm=T),digits=2)
# summary_scarf[2,1]=round(mean(Grass_SCARF_V,na.rm=T),digits = 2)
# summary_scarf[2,2]=round(sd(Grass_SCARF_V,na.rm=T),digits=2)
# summary_scarf[3,1]=round(mean(loss_SCARF_V,na.rm=T),digits = 2)
# summary_scarf[3,2]=round(sd(loss_SCARF_V,na.rm=T),digits=2)
# summary_scarf[4,1]=round(mean(Frob_SCARF_V,na.rm=T),digits = 2)
# summary_scarf[4,2]=round(sd(Frob_SCARF_V,na.rm=T),digits=2)
# summary_scarf[5,1]=round(mean(Frob_SCARF_B,na.rm=T),digits = 2)
# summary_scarf[5,2]=round(sd(Frob_SCARF_B,na.rm=T),digits=2)

#Save Results
#write.table(summary_scarf,file="summary_scarf_AJIVE_S.csv",sep=" ",col.names=T,row.names=T)
#save(sim_scarf_tune,file="sim_scarf_AJIVE_S.Rdata")

# ---------------------------------------------------

#Plot estimates of first data set vs. true para
# par(mar=c(4,4,4,4),mgp=c(2,0.8,0),oma=c(0,0,2,0))
# par(mfrow = c(4, 2))
# for ( i in 1:4){
#   plot(true_B[,i],ylab=bquote("B"~.(i)),ylim=c(-8,8),type="l",col="gold")
#   points(sim_scarf_tune[[1]]$B[,i],type="l",lty=2,col="dark green")
#   abline(v=cumsum(table(gpb)),lty=2,col="gray")
#   plot(true_V[,i],ylab=bquote("V"~.(i)),ylim=c(-0.8,0.8),type="l",col="red")
#   points(sim_scarf_tune[[1]]$V[,i],type="l",lty=2,col="dark green")
#   abline(v=cumsum(table(gp)),lty=2,col="gray")
# }
# title("SCARF Estimated B and V vs. True Settings",outer=T)
#quartz.save(file="AJIVE_S_Scarf.png", type = "png", device = dev.cur(),
#            dpi = 150, width=10, height=10)


