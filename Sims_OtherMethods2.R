# 100 Simulation Comparisons with SLIDE, SparsePCA, GFA, RRR #
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
Sim_dataX=AJIVE_SSF_X_Sn #LSF = large SF, SSF= Small SF
Sim_dataY=AJIVE_SSF_Y_Sn
true_B=AJIVE_True_B
true_V=AJIVE_True_V

# Setting of the groups, according to true_B and true_V
n=dim(Sim_dataY[[1]])[1]
p=dim(Sim_dataY[[1]])[2]
q=dim(Sim_dataX[[1]])[2]
gv=4; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b
gpb=ceiling(1:q / size.groupb)
gp=ceiling(1:p/ size.groupv)


##### SLIDE #####
source("./slide-paper-master/AuxillaryFunctions.R")
source("./slide-paper-master/SLIDEfunctions.R")
pvec=as.numeric(table(gp))


PA_slide=rep(NA,100) #Principal Angle
Grass_slide=rep(NA,100) #Grassmannian Distance
Frob_slide=rep(NA,100) #Frobenius Norm
loss_slide=rep(NA,100) #Loss Function
sim_slide=list()  #Store results
rep_n=1
while(rep_n<=100){
  dataY=as.matrix(Sim_dataY[[rep_n]])
  out_s <- standardizeX(dataY, pvec, center = T)
  SLIDE_X <- out_s$X
  svec <- out_s$svec
  fit_seq<- 
    solve_optim1_seq_restarts(X=SLIDE_X, lambda_seq = NULL, pvec = pvec, k_max = 1000, eps = 1e-8, reduced = F, rank_total = 4,
                              n_lambda=10 ,lambda_max = max(svec), lambda_min = 0.1)
  out_struct <- get_structure_v2(fit_seq, pvec)
  
  # Select the best structure from the list
  outbcv <- bcv_optim1_structure_centering(X=SLIDE_X, pvec = pvec, structure_list = out_struct$Slist, n_fold = 3, p_fold=3, k_max = 2000, eps = 1e-8, center = F)
  SLIDE_structure <- outbcv$structure_min
  
  SLIDE_param <- est_givenranks_v4(X = SLIDE_X, pvec = pvec, pattern = SLIDE_structure, k_max = 1000, eps = 1e-7)
  SLIDE_param$V <- apply(SLIDE_param$V,2,function(x) x/norm(x,"2"))
  
  sim_slide[[rep_n]]=SLIDE_param
  
  PA_slide[rep_n]=PrinAngle(true_V,sim_slide[[rep_n]]$V)$Max_angle
  Grass_slide[rep_n]=GrassDist(true_V,sim_slide[[rep_n]]$V)$Grass_Dist
  if (dim(sim_slide[[rep_n]]$V)[2]==dim(true_V)[2]){
    Frob_slide[rep_n]=norm(as.matrix(true_V-sim_slide[[rep_n]]$V),"F")^2
  }
  loss_slide[rep_n]=lossV(M_true=true_V,M_est=sim_slide[[rep_n]]$V,group=4)$Min_loss
  cat(rep_n,"simulation done.")
  rep_n=rep_n+1
}

#Result Summary
summary_slide=data.frame(matrix(NA,nrow=4,ncol=2),
                       row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm"))
colnames(summary_slide)=c("mean","sd")
summary_slide[1,1]=round(mean(PA_slide,na.rm=T),digits = 2)
summary_slide[1,2]=round(sd(PA_slide,na.rm=T),digits=2)
summary_slide[2,1]=round(mean(Grass_slide,na.rm=T),digits = 2)
summary_slide[2,2]=round(sd(Grass_slide,na.rm=T),digits=2)
summary_slide[3,1]=round(mean(loss_slide,na.rm=T),digits = 2)
summary_slide[3,2]=round(sd(loss_slide,na.rm=T),digits=2)
summary_slide[4,1]=round(mean(Frob_slide,na.rm=T),digits = 2)
summary_slide[4,2]=round(sd(Frob_slide,na.rm=T),digits=2)

#save small n result (Capital words -- Settings, Lower case -- Method)
# write.table(summary_slide,file="./Sim_Results_Slide/Small n/Summary_AJIVE_S_slide.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_slide,file="./Sim_Results_Slide/Large n/Summary_slide_AJIVE_S.csv",sep=" ",col.names=T,row.names=T)

##### GFA #####
library(CCAGFA)
pvec=as.numeric(table(gp))
gfa_opts = getDefaultOpts()
gfa_opts$verbose=0 # do not print iteration details in GFA
n_comp=4 # number of components

PA_gfa=rep(NA,100) #Principal Angle
Grass_gfa=rep(NA,100) #Grassmannian Distance
Frob_gfa=rep(NA,100) #Frobenius Norm
loss_gfa=rep(NA,100) #Loss Function
sim_gfa=list()  #Store results
rep_n=1
while(rep_n<=100){
  dataY=as.matrix(Sim_dataY[[rep_n]])
  data_gfa=list()
  for (g in 1:gv){
    data_gfa[[g]]=scale(dataY[,(size.groupv*(g-1)+1):(size.groupv*g)])
  }
  set.seed(12345)
  fit_gfa = GFAexperiment(Y = data_gfa, K = n_comp, opts = gfa_opts, Nrep=10)
  fit_gfa_trim = GFAtrim(fit_gfa)
  GFA_V=do.call("rbind", fit_gfa_trim$W)
  GFA_V=apply(GFA_V,2,function(x) x/norm(x,"2")) # make norm 1

  sim_gfa[[rep_n]]=GFA_V
  
  PA_gfa[rep_n]=PrinAngle(true_V,sim_gfa[[rep_n]])$Max_angle
  Grass_gfa[rep_n]=GrassDist(true_V,sim_gfa[[rep_n]])$Grass_Dist
  if(dim(sim_gfa[[rep_n]])[2]==dim(true_V)[2]){
    Frob_gfa[rep_n]=norm(as.matrix(true_V-sim_gfa[[rep_n]]),"F")^2
  }
  loss_gfa[rep_n]=lossV(M_true=true_V,M_est=sim_gfa[[rep_n]],group=4)$Min_loss
  cat(rep_n,"simulation done.")
  rep_n=rep_n+1
}

#Result Summary
summary_gfa <- data.frame(matrix(NA,nrow=4,ncol=2),
                         row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm"))
colnames(summary_gfa) <- c("mean","sd")
summary_gfa[1,1] <- round(mean(PA_gfa,na.rm=T),digits = 2)
summary_gfa[1,2] <- round(sd(PA_gfa,na.rm=T),digits=2)
summary_gfa[2,1] <- round(mean(Grass_gfa,na.rm=T),digits = 2)
summary_gfa[2,2] <- round(sd(Grass_gfa,na.rm=T),digits=2)
summary_gfa[3,1] <- round(mean(loss_gfa,na.rm=T),digits = 2)
summary_gfa[3,2] <- round(sd(loss_gfa,na.rm=T),digits=2)
summary_gfa[4,1] <- round(mean(Frob_gfa,na.rm=T),digits = 2)
summary_gfa[4,2] <- round(sd(Frob_gfa,na.rm=T),digits=2)

# #save small n result (Capital words -- Settings, Lower case -- Method)
## write.table(summary_gfa,file="./Sim_Results_Sn/Summary_CBOS_L_gfa.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
write.table(summary_gfa,file="./Sim_Results_New/Summary_gfa_AJIVE_S.csv",sep=" ",col.names=T,row.names=T)

##### RRR and SRRR #####
library(rrr)
library(rrpack)

PA_rrr_V=rep(NA,100); PA_srrr_V=rep(NA,100) 
Grass_rrr_V=rep(NA,100); Grass_srrr_V=rep(NA,100)
loss_rrr_V=rep(NA,100); loss_srrr_V=rep(NA,100)
Frob_rrr_V=rep(NA,100); Frob_srrr_V=rep(NA,100)
Frob_rrr_B=rep(NA,100); Frob_srrr_B=rep(NA,100)
sim_rrr=list(); sim_srrr=list()
rep_n=1
while(rep_n<=100){
  dataY <- as.matrix(Sim_dataY[[rep_n]])
  dataX <- as.matrix(Sim_dataX[[rep_n]])
  
  # reduced rank regression given fixed rank
  fit_rrr <- rrpack::rrr.fit(Y=dataY, X=dataX, nrank = 4, weight = NULL, coefSVD = T)
  rrr_B <- as.matrix(fit_rrr$coefSVD$u)%*%diag(fit_rrr$coefSVD$d)
  rrr_V <- as.matrix(fit_rrr$coefSVD$v)
  sim_rrr[[rep_n]]=list(B=rrr_B, V=rrr_V)
  # sparse reduced rank regression tuned by cross validation
  fit_srrr <- rrpack::cv.srrr(Y=dataY, X=dataX, nrank = 4, method = c("glasso"))
  srrr_B <- fit_srrr$U%*%fit_srrr$D
  srrr_V <- fit_srrr$V
  sim_srrr[[rep_n]]=list(B=srrr_B, V=srrr_V)
    
  for (j in 1:4){
    if (angle(sim_rrr[[rep_n]]$V[,j],true_V[,j])>90) {
      sim_rrr[[rep_n]]$V[,j]=-sim_rrr[[rep_n]]$V[,j]
      sim_rrr[[rep_n]]$B[,j]=-sim_rrr[[rep_n]]$B[,j]}
    if (angle(sim_srrr[[rep_n]]$V[,j],true_V[,j])>90) {
      sim_srrr[[rep_n]]$V[,j]=-sim_srrr[[rep_n]]$V[,j]
      sim_srrr[[rep_n]]$B[,j]=-sim_srrr[[rep_n]]$B[,j]}
  }
  PA_rrr_V[rep_n]=PrinAngle(true_V,sim_rrr[[rep_n]]$V)$Max_angle
  Grass_rrr_V[rep_n]=GrassDist(true_V,sim_rrr[[rep_n]]$V)$Grass_Dist
  Frob_rrr_V[rep_n]=norm(as.matrix(true_V-sim_rrr[[rep_n]]$V),"F")^2
  Frob_rrr_B[rep_n]=norm(as.matrix(true_B-sim_rrr[[rep_n]]$B),"F")^2
  loss_rrr_V[rep_n]=lossV(M_true=true_V,M_est=sim_rrr[[rep_n]]$V,group=4)$Min_loss
  
  PA_srrr_V[rep_n]=PrinAngle(true_V,sim_srrr[[rep_n]]$V)$Max_angle
  Grass_srrr_V[rep_n]=GrassDist(true_V,sim_srrr[[rep_n]]$V)$Grass_Dist
  Frob_srrr_V[rep_n]=norm(as.matrix(true_V-sim_srrr[[rep_n]]$V),"F")^2
  Frob_srrr_B[rep_n]=norm(as.matrix(true_B-sim_srrr[[rep_n]]$B),"F")^2
  loss_srrr_V[rep_n]=lossV(M_true=true_V,M_est=sim_srrr[[rep_n]]$V,group=4)$Min_loss
  
  cat(rep_n,"simulation done.")
  rep_n=rep_n+1
}

#Result Summary - rrr
summary_rrr=data.frame(matrix(NA,nrow=5,ncol=2),
                          row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm","BFrobNorm"))
colnames(summary_rrr)=c("mean","sd")
summary_rrr[1,1]=round(mean(PA_rrr_V,na.rm=T),digits = 2)
summary_rrr[1,2]=round(sd(PA_rrr_V,na.rm=T),digits=2)
summary_rrr[2,1]=round(mean(Grass_rrr_V,na.rm=T),digits = 2)
summary_rrr[2,2]=round(sd(Grass_rrr_V,na.rm=T),digits=2)
summary_rrr[3,1]=round(mean(loss_rrr_V,na.rm=T),digits = 2)
summary_rrr[3,2]=round(sd(loss_rrr_V,na.rm=T),digits=2)
summary_rrr[4,1]=round(mean(Frob_rrr_V,na.rm=T),digits = 2)
summary_rrr[4,2]=round(sd(Frob_rrr_V,na.rm=T),digits=2)
summary_rrr[5,1]=round(mean(Frob_rrr_B,na.rm=T),digits = 2)
summary_rrr[5,2]=round(sd(Frob_rrr_B,na.rm=T),digits=2)

#save small n result (Capital words -- Settings, Lower case -- Method)
## write.table(summary_rrr,file="./Sim_Results_Sn/Summary_AJIVE_S_rrr.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_rrr,file="./Sim_Results_New/Summary_rrr_AJIVE_S.csv",sep=" ",col.names=T,row.names=T)

#Result Summary - srrr
summary_srrr=data.frame(matrix(NA,nrow=5,ncol=2),
                       row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm","BFrobNorm"))
colnames(summary_srrr)=c("mean","sd")
summary_srrr[1,1]=round(mean(PA_srrr_V,na.rm=T),digits = 2)
summary_srrr[1,2]=round(sd(PA_srrr_V,na.rm=T),digits=2)
summary_srrr[2,1]=round(mean(Grass_srrr_V,na.rm=T),digits = 2)
summary_srrr[2,2]=round(sd(Grass_srrr_V,na.rm=T),digits=2)
summary_srrr[3,1]=round(mean(loss_srrr_V,na.rm=T),digits = 2)
summary_srrr[3,2]=round(sd(loss_srrr_V,na.rm=T),digits=2)
summary_srrr[4,1]=round(mean(Frob_srrr_V,na.rm=T),digits = 2)
summary_srrr[4,2]=round(sd(Frob_srrr_V,na.rm=T),digits=2)
summary_srrr[5,1]=round(mean(Frob_srrr_B,na.rm=T),digits = 2)
summary_srrr[5,2]=round(sd(Frob_srrr_B,na.rm=T),digits=2)

#save small n result (Capital words -- Settings, Lower case -- Method)
## write.table(summary_srrr,file="./Sim_Results_Sn/Summary_AJIVE_S_srrr.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_srrr,file="./Sim_Results_New/Summary_srrr_AJIVE_S.csv",sep=" ",col.names=T,row.names=T)

#### SPCA #####
library(elasticnet)
library(PMA)

PA_spca=rep(NA,100) #Principal Angle
Grass_spca=rep(NA,100) #Grassmannian Distance
Frob_spca=rep(NA,100) #Frobenius Norm
loss_spca=rep(NA,100) #Loss Function
sim_spca=list()  #Store results
rep_n=1
while(rep_n<=100){
  dataY <- as.matrix(Sim_dataY[[rep_n]])
  bestsumabsv <- PMA::SPC.cv(x=dataY,trace=F)$bestsumabsv # use CV to choose best l1 norm for SPCA
  fit_spca <- elasticnet::spca(x=dataY,K=4,type="predictor",sparse="penalty",para=rep(bestsumabsv,4))
  spca_v <- fit_spca$loadings
  sim_spca[[rep_n]]<- spca_v

  for (j in 1:4){
    if (angle(sim_spca[[rep_n]][,j],true_V[,j])>90){
      sim_spca[[rep_n]][,j]=-sim_spca[[rep_n]][,j]
    }
  }
  PA_spca[rep_n]=PrinAngle(true_V,sim_spca[[rep_n]])$Max_angle
  Grass_spca[rep_n]=GrassDist(true_V,sim_spca[[rep_n]])$Grass_Dist
  Frob_spca[rep_n]=norm(as.matrix(true_V-sim_spca[[rep_n]]),"F")^2
  loss_spca[rep_n]=lossV(M_true=true_V,M_est=sim_spca[[rep_n]],group=4)$Min_loss
  cat(rep_n,"simulation done.")
  rep_n=rep_n+1
}

#Result Summary
summary_spca=data.frame(matrix(NA,nrow=4,ncol=2),
                       row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm"))
colnames(summary_spca)=c("mean","sd")
summary_spca[1,1]=round(mean(PA_spca,na.rm=T),digits = 2)
summary_spca[1,2]=round(sd(PA_spca,na.rm=T),digits=2)
summary_spca[2,1]=round(mean(Grass_spca,na.rm=T),digits = 2)
summary_spca[2,2]=round(sd(Grass_spca,na.rm=T),digits=2)
summary_spca[3,1]=round(mean(loss_spca,na.rm=T),digits = 2)
summary_spca[3,2]=round(sd(loss_spca,na.rm=T),digits=2)
summary_spca[4,1]=round(mean(Frob_spca,na.rm=T),digits = 2)
summary_spca[4,2]=round(sd(Frob_spca,na.rm=T),digits=2)

#save small n result (Capital words -- Settings, Lower case -- Method)
## write.table(summary_spca,file="./Sim_Results_Sn/Summary_AJIVE_S_spca.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_spca,file="./Sim_Results_New/Summary_spca_AJIVE_S.csv",sep=" ",col.names=T,row.names=T)



