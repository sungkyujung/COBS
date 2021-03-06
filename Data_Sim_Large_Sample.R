##### Try settings: #####

# 1. AJIVE (no covariates) vs. 
#    SIFA (no partial) vs. 
#    COBS (covariates & partial);
# 2. Large \Sigma_F vs. Small \Sigma_F

source('~/routines')

##### COBS #####

n=500  # Sample size, signal
p=100  # Dim of Y,V
q=40   # Dim of X,B
g=4    # g is the intrinsic rank of Y

gv=4; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b

# True B
gpb=ceiling(1:q / size.groupb)
B=matrix(rep(0,q*g),ncol=g)
B[11:13,1]=5
B[1:3,2]=-4
B[21:23,3]=-3
B[31:33,4]=2
B_group=data.frame(cbind(gpb,B))
colnames(B_group)=c("gpb","b1","b2","b3","b4")

# True V
gp=ceiling(1:p/ size.groupv)
V=matrix(rep(0,p*g),ncol=g)
V[1:5,1]=1
V[26:30,2]=-1
V[51:55,3]=1;V[76:80,3]=1
V[24:25,4]=1; V[49:50,4]=1; V[74:75,4]=-1; V[99:100,4]=-1
for(i in 1:g){
  V[,i]=V[,i]/norm(as.matrix(V[,i]),"F")
}

# V4.proj=matrix(NA,nrow=p,ncol=g-1)
# for (j in 1:(g-1)){
#   V4.proj[,j]=as.matrix(Proj(vec1=as.matrix(V[,j]),vec2=V[,4]))
# }
# Sum.proj4=as.matrix(apply(V4.proj,1,FUN="sum"))
# V4_proj=V[,4]-Sum.proj4
# V4_new=V4_proj/norm(V4_proj,"F")
# for (i in 1:p){
#   if (abs(V4_new[i])<(10^-9)) {V4_new[i]=0}
# }
# V4_new=V4_new/norm(V4_new,"F")
# V[,4]=V4_new

V_group=data.frame(cbind(gp,V))
colnames(V_group)=c("gp","v1","v2","v3","v4")

true_B=B; SCARF_True_B=true_B
true_V=V; SCARF_True_V=true_V

##### Settings #####
set.seed(1234)
Sigma.X=c(1,1,1,1)
se=1
E=matrix(rnorm(n*p, 0, se),n,p)

#SF=c(2.5,2,1.5,1) #Smaller variance ratio
#SF=c(10,8,6,4) #Larger variance ratio

Sigma.F=diag(SF,ncol=4,nrow=4) # Covariance matrix

true_se=se
true_Sigma.F=Sigma.F

# 100 simulated Data: main data Y and supervision X 
Sim_dataX=list()
Sim_dataY=list()
Sim_Vrank=list()
sim_n=1
try.count=1
while (sim_n <= 100){
  # random normal X
  sim_n=sim_n
  dataX=list()
  for (i in 1:gb){
    dataX[[i]]=matrix(rnorm(n*size.groupb,mean=0,sd=Sigma.X[i]),ncol=size.groupb,nrow=n)
  }
  X=do.call("cbind", dataX)
  #X=matrix(rnorm(n*q,0,1),nrow=n,ncol=q)
  X=scale(X,center = T, scale=F)
  dataX=X
  Sim_dataX[[sim_n]]=dataX
  
  # random normal f & E, with fixed Sigma.F & se
  f=as.matrix(mvrnorm(n=n,mu=matrix(c(0,0,0,0),byrow=F),Sigma=Sigma.F))
  E=matrix(rnorm(n*p, 0, se),n,p)
  # Y generated from the model Y=XBV'+FV'+E  
  Y=X%*%true_B%*%t(true_V)+f%*%t(true_V)+E
  Y=scale(Y,center=T,scale=F)
  dataY=Y
  Sim_dataY[[sim_n]]=dataY
  
  V_rank=rank(-(apply(dataX%*%true_B,2,var)+SF))
  Sim_Vrank[[sim_n]]=V_rank
  if (sum(Sim_Vrank[[sim_n]]==c(1,2,3,4))!=4){sim_n=sim_n} else {sim_n=sim_n+1}
  try.count=try.count+1
}
#Sim_Vrank=do.call("rbind", Sim_Vrank)

# Change this part according to the setting SF; Sigma.F
SCARF_SSF_X=Sim_dataX   # SCARF_LSF_X=Sim_dataX 
SCARF_SSF_Y=Sim_dataY   # SCARF_LSF_Y=Sim_dataY

# save("SCARF_LSF_X","SCARF_LSF_Y",
#      "SCARF_SSF_X","SCARF_SSF_Y",
#      "SCARF_True_B","SCARF_True_V",
#      file="Sim_SCARF.Rdata")

##### SIFA #####
n=500  # Sample size, signal
p=100  # Dim of Y,V
q=40   # Dim of X,B
g=4    # g is the intrinsic rank of Y

gv=4; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b

# True B
gpb=ceiling(1:q / size.groupb)
B=matrix(rep(0,q*g),ncol=g)
B[11:13,1]=5
B[1:3,2]=-4
B[21:23,3]=-3
B[31:33,4]=2
#B_group=data.frame(cbind(gpb,B))
#colnames(B_group)=c("gpb","b1","b2","b3","b4")

B0=as.matrix(B[,1]);B1=as.matrix(B[,2]);B2=as.matrix(B[,3]);B3=as.matrix(B[,4]);B4=matrix(0,nrow=q,ncol=1)
B=as.matrix(cbind(B0,B1,B2,B3))

# True V
gp=ceiling(1:p/ size.groupv)
V=matrix(rep(0,p*g),ncol=g)
V[1:5,1]=1
V[26:30,2]=-1
V[51:55,3]=1
V[24:25,4]=1; V[49:50,4]=1; V[74:75,4]=-1; V[99:100,4]=-1
for(i in 1:g){
  V[,i]=V[,i]/norm(as.matrix(V[,i]),"F")
}

# V4.proj=matrix(NA,nrow=p,ncol=g-1)
# for (j in 1:(g-1)){
#   V4.proj[,j]=as.matrix(Proj(vec1=as.matrix(V[,j]),vec2=V[,4]))
# }
# Sum.proj4=as.matrix(apply(V4.proj,1,FUN="sum"))
# V4_proj=V[,4]-Sum.proj4
# V4_new=V4_proj/norm(V4_proj,"F")
# for (i in 1:p){
#   if (abs(V4_new[i])<(10^-9)) {V4_new[i]=0}
# }
# V4_new=V4_new/norm(V4_new,"F")
# V[,4]=V4_new

#V_group=data.frame(cbind(gp,V))
#colnames(V_group)=c("gp","v1","v2","v3","v4")

library(pracma)
V0=as.matrix(V[,4]);V1=as.matrix(V[1:25,1]);V2=as.matrix(V[26:50,2]);V3=as.matrix(V[51:75,3]);V4=matrix(0,nrow=25,ncol=1)
V=as.matrix(cbind(V0,blkdiag(V1,V2,V3,V4)))[,-5]

true_B=B; SIFA_True_B=true_B
true_V=V; SIFA_True_V=true_V

##### Settings #####
set.seed(1234)
Sigma.X=c(1,1,1,1)
se=1
E=matrix(rnorm(n*p, 0, se),n,p)

# SF0=2.5;SF1=2;SF2=1.5;SF3=1;SF4=0 #Smaller variance ratio
# SF0=10;SF1=8;SF2=6;SF3=4;SF4=0 #Larger variance ratio
SF=c(SF0,SF1,SF2,SF3)
#Sigma.F=diag(SF,ncol=4,nrow=4) # Covariance matrix

true_se=se

# 100 simulated Data: main data Y and supervision X 
Sim_dataX=list()
Sim_dataY=list()
Sim_Vrank=list()
sim_n=1
try.count=1
while (sim_n <= 100){
  # random normal X
  sim_n=sim_n
  dataX=list()
  for (i in 1:gb){
    dataX[[i]]=matrix(rnorm(n*size.groupb,mean=0,sd=Sigma.X[i]),ncol=size.groupb,nrow=n)
  }
  X=do.call("cbind", dataX)
  #X=matrix(rnorm(n*q,0,1),nrow=n,ncol=q)
  X=scale(X,center = T, scale=F)
  dataX=X
  Sim_dataX[[sim_n]]=dataX
  
  # random normal f & E, with fixed Sigma.F & se
  #f=as.matrix(mvrnorm(n=n,mu=matrix(c(0,0,0,0),byrow=F),Sigma=Sigma.F))
  F0=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF0))
  F1=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF1))
  F2=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF2))
  F3=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF3))
  F4=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF4))
  E=matrix(rnorm(n*p, 0, se),n,p)
  E1=as.matrix(E[,1:25]); E2=as.matrix(E[,26:50]); E3=as.matrix(E[,51:75]); E4=as.matrix(E[,76:100]);
  # Y generated from the SIFA model Y_k=J_k+A_K+E_k  
  #Y=X%*%true_B%*%t(true_V)+f%*%t(true_V)+E
  
  U0=dataX%*%B0+F0; U1=dataX%*%B1+F1; U2=dataX%*%B2+F2; U3=dataX%*%B3+F3; U4=dataX%*%B4+F4
  Joint=U0%*%t(V0);
  Jnt1=as.matrix(Joint[,1:25]);Jnt2=as.matrix(Joint[,26:50]);Jnt3=as.matrix(Joint[,51:75]);Jnt4=as.matrix(Joint[,76:100]);
  Ind1=U1%*%t(V1); Ind2=U2%*%t(V2); Ind3=U3%*%t(V3); Ind4=U4%*%t(V4);
  Y1=Jnt1+Ind1+E1;Y2=Jnt2+Ind2+E2;Y3=Jnt3+Ind3+E3;Y4=Jnt4+Ind4+E4;
  
  Y=as.matrix(cbind(Y1,Y2,Y3,Y4))
  Y=scale(Y,center=T,scale=F)
  dataY=Y
  Sim_dataY[[sim_n]]=dataY
  
  V_rank=rank(c(-(var(dataX%*%true_B[,1])+SF0),-(var(dataX%*%true_B[,2])+SF1),-(var(dataX%*%true_B[,3])+SF2),-(var(dataX%*%true_B[,4])+SF3)))
  Sim_Vrank[[sim_n]]=V_rank
  if (sum(Sim_Vrank[[sim_n]]==c(1,2,3,4))!=4){sim_n=sim_n} else {sim_n=sim_n+1}
  try.count=try.count+1
}
#Sim_Vrank=do.call("rbind", Sim_Vrank)

# Change this part according to the setting SF; Sigma.F
SIFA_SSF_X=Sim_dataX   # SIFA_LSF_X=Sim_dataX 
SIFA_SSF_Y=Sim_dataY   # SIFA_LSF_Y=Sim_dataY

# save("SIFA_LSF_X","SIFA_LSF_Y",
#      "SIFA_SSF_X","SIFA_SSF_Y",
#      "SIFA_True_B","SIFA_True_V",
#      file="Sim_SIFA.Rdata")



##### AJIVE #####
n=500  # Sample size, signal
p=100  # Dim of Y,V
q=40   # Dim of X,B
g=4    # g is the intrinsic rank of Y

gv=4; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b

# True B
gpb=ceiling(1:q / size.groupb)
B=matrix(rep(0,q*g),ncol=g)
#B_group=data.frame(cbind(gpb,B))
#colnames(B_group)=c("gpb","b1","b2","b3","b4")

B0=as.matrix(B[,1]);B1=as.matrix(B[,2]);B2=as.matrix(B[,3]);B3=as.matrix(B[,4]);B4=matrix(0,nrow=q,ncol=1)
B=as.matrix(cbind(B0,B1,B2,B3))

# True V
# True V
gp=ceiling(1:p/ size.groupv)
V=matrix(rep(0,p*g),ncol=g)
V[1:5,1]=1
V[26:30,2]=-1
V[51:55,3]=1
V[24:25,4]=1; V[49:50,4]=1; V[74:75,4]=-1; V[99:100,4]=-1
for(i in 1:g){
  V[,i]=V[,i]/norm(as.matrix(V[,i]),"F")
}

# V4.proj=matrix(NA,nrow=p,ncol=g-1)
# for (j in 1:(g-1)){
#   V4.proj[,j]=as.matrix(Proj(vec1=as.matrix(V[,j]),vec2=V[,4]))
# }
# Sum.proj4=as.matrix(apply(V4.proj,1,FUN="sum"))
# V4_proj=V[,4]-Sum.proj4
# V4_new=V4_proj/norm(V4_proj,"F")
# for (i in 1:p){
#   if (abs(V4_new[i])<(10^-9)) {V4_new[i]=0}
# }
# V4_new=V4_new/norm(V4_new,"F")
# V[,4]=V4_new

#V_group=data.frame(cbind(gp,V))
#colnames(V_group)=c("gp","v1","v2","v3","v4")

V0=as.matrix(V[,4]);V1=as.matrix(V[1:25,1]);V2=as.matrix(V[26:50,2]);V3=as.matrix(V[51:75,3]);V4=matrix(0,nrow=25,ncol=1)
V=as.matrix(cbind(V0,blkdiag(V1,V2,V3,V4)))[,-5]

true_B=B; AJIVE_True_B=true_B
true_V=V; AJIVE_True_V=true_V

##### Settings #####
set.seed(1234)
Sigma.X=c(1,1,1,1)
se=1
E=matrix(rnorm(n*p, 0, se),n,p)

# SF0=2.5;SF1=2;SF2=1.5;SF3=1;SF4=0 #Smaller variance ratio
# SF0=10;SF1=8;SF2=6;SF3=4;SF4=0 #Larger variance ratio
SF=c(SF0,SF1,SF2,SF3)
#Sigma.F=diag(SF,ncol=4,nrow=4) # Covariance matrix

true_se=se

# 100 simulated Data: main data Y and supervision X 
Sim_dataX=list()
Sim_dataY=list()
Sim_Vrank=list()
sim_n=1
try.count=1
while (sim_n <= 100){
  # random normal X
  sim_n=sim_n
  #dataX=list()
  #for (i in 1:gb){
  #  dataX[[i]]=matrix(rnorm(n*size.groupb,mean=0,sd=Sigma.X[i]),ncol=size.groupb,nrow=n)
  #}
  #X=do.call("cbind", dataX)
  X=matrix(rnorm(n*q,0,1),nrow=n,ncol=q)
  X=scale(X,center = T, scale=F)
  dataX=X
  Sim_dataX[[sim_n]]=dataX
  
  # random normal f & E, with fixed Sigma.F & se
  #f=as.matrix(mvrnorm(n=n,mu=matrix(c(0,0,0,0),byrow=F),Sigma=Sigma.F))
  F0=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF0))
  F1=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF1))
  F2=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF2))
  F3=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF3))
  F4=as.matrix(mvrnorm(n=n,mu=0,Sigma=SF4))
  E=matrix(rnorm(n*p, 0, se),n,p)
  E1=as.matrix(E[,1:25]); E2=as.matrix(E[,26:50]); E3=as.matrix(E[,51:75]); E4=as.matrix(E[,76:100]);
  # Y generated from the SIFA model Y_k=J_k+A_K+E_k  
  #Y=X%*%true_B%*%t(true_V)+f%*%t(true_V)+E
  
  U0=dataX%*%B0+F0; U1=dataX%*%B1+F1; U2=dataX%*%B2+F2; U3=dataX%*%B3+F3; U4=dataX%*%B4+F4
  Joint=U0%*%t(V0);
  Jnt1=as.matrix(Joint[,1:25]);Jnt2=as.matrix(Joint[,26:50]);Jnt3=as.matrix(Joint[,51:75]);Jnt4=as.matrix(Joint[,76:100]);
  Ind1=U1%*%t(V1); Ind2=U2%*%t(V2); Ind3=U3%*%t(V3); Ind4=U4%*%t(V4);
  Y1=Jnt1+Ind1+E1;Y2=Jnt2+Ind2+E2;Y3=Jnt3+Ind3+E3;Y4=Jnt4+Ind4+E4;
  Y1=scale(Y1,center=T,scale=F)
  Y2=scale(Y2,center=T,scale=F)
  Y3=scale(Y3,center=T,scale=F)
  Y4=scale(Y4,center=T,scale=F)
  
  Y=as.matrix(cbind(Y1,Y2,Y3,Y4))
  dataY=Y
  Sim_dataY[[sim_n]]=dataY
  
  V_rank=rank(c(-(var(dataX%*%true_B[,1])+SF0),-(var(dataX%*%true_B[,2])+SF1),-(var(dataX%*%true_B[,3])+SF2),-(var(dataX%*%true_B[,4])+SF3)))
  Sim_Vrank[[sim_n]]=V_rank
  if (sum(Sim_Vrank[[sim_n]]==c(1,2,3,4))!=4){sim_n=sim_n} else {sim_n=sim_n+1}
  try.count=try.count+1
}
#Sim_Vrank=do.call("rbind", Sim_Vrank)

# Change this part according to the setting SF; Sigma.F
AJIVE_SSF_X=Sim_dataX   # AJIVE_LSF_X=Sim_dataX 
AJIVE_SSF_Y=Sim_dataY   # AJIVE_LSF_Y=Sim_dataY

# save("AJIVE_LSF_X","AJIVE_LSF_Y",
#      "AJIVE_SSF_X","AJIVE_SSF_Y",
#      "AJIVE_True_B","AJIVE_True_V",
#      file="Sim_AJIVE.Rdata")


##### Save data in csv for SIFA in matlab #########################################################
# Sim_dataX_100=do.call("rbind", AJIVE_SSF_X) #every n (sample size) row is one data set
# Sim_dataY_100=do.call("rbind", AJIVE_SSF_Y) #every n (sample size) row is one data set
# 
# write.table(Sim_dataX_100,file="/Users/stargao/Desktop/New Experiment/New_Sims_1028/For_SIFA_Matlab/DataX_Ajive_S.csv",sep=" ",col.names=F,row.names=F)
# write.table(Sim_dataY_100,file="/Users/stargao/Desktop/New Experiment/New_Sims_1028/For_SIFA_Matlab/DataY_Ajive_S.csv",sep=" ",col.names=F,row.names=F)




