##### New Tune Experiment with Simulated Data ######

rm(list = ls()) 
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
Sim_dataX=SIFA_LSF_X_Sn #LSF = large SF, SSF= Small SF
Sim_dataY=SIFA_LSF_Y_Sn
true_B=SIFA_True_B
true_V=SIFA_True_V

# Setting of the groups, according to true_B and true_V
n=dim(Sim_dataY[[1]])[1]
p=dim(Sim_dataY[[1]])[2]
q=dim(Sim_dataX[[1]])[2]
gv=4; size.groupv=p/gv  # group number of Y and v
gb=4; size.groupb=q/gb   # group number of X and b
gpb=ceiling(1:q / size.groupb)
gp=ceiling(1:p/ size.groupv)

# Structured Factorization with Thresholding ------------------------------
# Set tuning parameters and inputs 

# lambda.b > 0 (overall sparsity parameter for coefficient B)
# alpha.b in (0,1) (if 0, group lasso, if 1. lasso. If 1 (or 0), the (group) lasso parameter is lambda.b)
# lambda.v > 0 (overall sparsity parameter for loading vector V)
# alpha.v in (0,1) (if 0, group lasso, if 1. lasso.)
alpha.b=0.2; lambda.b.seq=c(0.15,0.1,0.05,0.01,0.005)
alpha.v=0.2; lambda.v.seq=c(0.3,0.25,0.2,0.15,0.01) #lambda.v.seq=c(0.27,0.2,0.15,0.07,0.001)
n.comp <- 4

# other inputs are dataX, dataY, Y.bk.idx, X.bk.inx 

PA_COBS_V=rep(NA,100) #Principal Angle
Grass_COBS_V=rep(NA,100) #Grassmannian Distance
Frob_COBS_V=rep(NA,100) #Frobenius Norm
Frob_COBS_B=rep(NA,100) #Frobenius Norm
loss_COBS_V=rep(NA,100) #Loss Function
rep_n=1
while(rep_n <= 100){
  # Prepare output
  est <- list()
  
  # Initial value -----------------------------------------------------------
  
  outFlag <- FALSE # A global parameter (if true, stop entirely)
  maxLoop <- 100 # A global parameter (maximum iterates, for each component)
  TOL <- .Machine$double.eps # A global parameter for convergence (~ 2e-16)
  lambda.b.tuned <- NULL
  
  dataX=Sim_dataX[[rep_n]]
  dataY=Sim_dataY[[rep_n]]
  
  X <- as.matrix(dataX) # centered
  Y <- as.matrix(dataY)
  
  data=list()
  data.v=list()
  data.b=list()
  data.s=list()
  bic.list=list()
  
  l_b_track=rep(NA,n.comp)
  l_v_track=rep(NA,n.comp)
  
  for (i.layer in 1:n.comp) {
    
    # Prepare the primary data, subtracting all previous components
    # (done by sequentially subtracting from Y)
    if (i.layer == 1){
      data[[i.layer]] <- as.matrix(dataY) # centered
    }else{ 
      data[[i.layer]] <-  data[[i.layer-1]] - (est[[i.layer-1]]$s[2] *data[[i.layer-1]] %*% est[[i.layer-1]]$v  
                 + est[[i.layer-1]]$s[1] * X %*% est[[i.layer-1]]$b) %*% 
        t(est[[i.layer-1]]$v) / (sum(est[[i.layer-1]]$s)) 
    }
    
    bic_track=matrix(NA,length(lambda.b.seq),length(lambda.v.seq), byrow=T)
    V_store=list()
    B_store=list()
    s_store=list()
    
    for (l1 in 1:length(lambda.b.seq)){
      for (l2 in 1:length(lambda.v.seq)){
        tune.par.b <- list(alpha = alpha.b, lambda = lambda.b.seq[l1])
        tune.par.v <- list(alpha = alpha.v, lambda = lambda.v.seq[l2])
        Y=as.matrix(data[[i.layer]])
        
        # Prepare for iterations --------------------------------------------------
        Y_svd <-svd(Y) 
        init <- get.initial(X = X,Y_svd = Y_svd,Y=Y,param.b=tune.par.b,grp.index.X=gpb) 
        b.now <- init$b.init
        v.now <- init$v.init
        s.now <- init$s.init
        
        convergeFlag <- FALSE # for this component
        cnt <- 1
        
        # The inside loop ---------------------------------------------------------
        # Loop while outFlag == TRUE && convergeFlag == TRUE && cnt < maxLoop
        
        
        while( !outFlag  && !convergeFlag && cnt < maxLoop){
          
          # inside loop
          
          # Update V ----------------------------------------------------------------
          v.next  <- update.v(Y_svd=Y_svd, X=X, Y=Y, b = b.now, s = s.now,param = tune.par.v, grp.index = gp)
          
          # If V is a zero vector (resulting in v = NaN), then stop 
          if( any(is.nan(v.next)) ) {
            outFlag <- TRUE
            break;
          }
          
          v.next <- sign(sum(v.next * v.now)) * v.next
          if (sign(sum(v.next * v.now)) < 0 ) cat("v sign flipped")
          
          # Update B ----------------------------------------------------------------
          b.tmp <- update.b(X = X, Y = Y,v = v.next, param = tune.par.b,grp.index.X = gpb)
          b.next <- b.tmp$bhat 
          
          if (is.null(tune.par.b$lambda)) { 
            lambda.b.tuned <- b.tmp$lambda
          } # only if the tuning parameter for B estimate was not supplied. 
          
          # Update Sigma.est --------------------------------------------------------
          s.next <- unlist(Sigma.est(Y=Y,X=X,b=b.next, v=v.next))
          
          b.diff <- sum( ( b.now - b.next )^2 ) 
          v.diff <- 1- abs(sum( v.now * v.next ) ) 
          if( b.diff < TOL && v.diff < TOL ){
            convergeFlag <- TRUE
          }
          b.now <- b.next 
          v.now <- v.next 
          s.now <- s.next 
          cnt <- cnt + 1 
        }
        
        b_temp=b.now
        v_temp=v.now
        s_temp=s.now
        
        bic=fbic(X=X,Y=as.matrix(data[[i.layer]]),b=b_temp,v=v_temp,se=sqrt(s.now[1]),sf=sqrt(s.now[2]))$bic
        
        bic_track[l1,l2]=bic
        B_store[[(l1-1)*length(lambda.v.seq)+l2]]=b_temp
        V_store[[(l1-1)*length(lambda.v.seq)+l2]]=v_temp
        s_store[[(l1-1)*length(lambda.v.seq)+l2]]=s_temp
      }
    }
    
    min_row=which(bic_track == min(bic_track,na.rm=T), arr.ind = TRUE)[1] # row location of min bic to find lambda.b
    min_col=which(bic_track == min(bic_track,na.rm=T), arr.ind = TRUE)[2] # col location of min bic to find lambda.v
    l_b_track[i.layer]=lambda.b.seq[min_row]
    l_v_track[i.layer]=lambda.v.seq[min_col]
    
    v.now=V_store[[(min_row-1)*length(lambda.v.seq)+min_col]]
    b.now=B_store[[(min_row-1)*length(lambda.v.seq)+min_col]]
    s.now=s_store[[(min_row-1)*length(lambda.v.seq)+min_col]]
  
    # Orthogonalize 
    # v.orig <- v.now
    # if (i.layer > 1){
    #   V.prev <- matrix(nrow = p, ncol = i.layer-1)
    #   for (j in (1:(i.layer-1))) { V.prev[,j] <- est[[j]]$v }
    #   v.orth <- v.now - V.prev %*% t(V.prev) %*%v.now
    #   v.now <- v.orth / sqrt(sum(v.orth*v.orth))
    # }
    
    # Save the result at the "est" list
    est[[i.layer]] <- list(b = b.now, v = v.now, s = s.now, #v.orig = v.orig,
                           lambda.b.tuned = lambda.b.tuned)
    # End of inside loop ------------------------------------------------------
    
    
    # store results for final variance calculation and results comparison
    data.v[[i.layer]]=v.now  # final estimated v_hat for each layer
    data.b[[i.layer]]=b.now
    data.s[[i.layer]]=s.now
#  bic.list[[i.layer]]=bic_track
  }
  
  #Results comparisons
  V_COBS=do.call("cbind",data.v)
  B_COBS=do.call("cbind",data.b)
  # BIC_COBS=bic.list
  
  for (j in 1:n.comp){
    if (angle(V_COBS[,j],true_V[,j])>90) {
      V_COBS[,j]=-V_COBS[,j]
      B_COBS[,j]=-B_COBS[,j]
      }
  }
  
  PA_COBS_V[rep_n]=PrinAngle(true_V,V_COBS)$Max_angle
  Grass_COBS_V[rep_n]=GrassDist(true_V,V_COBS)$Grass_Dist
  Frob_COBS_V[rep_n]=norm(as.matrix(true_V-V_COBS),"F")^2
  Frob_COBS_B[rep_n]=norm(as.matrix(true_B-B_COBS),"F")^2
  loss_COBS_V[rep_n]=lossV(M_true=true_V,M_est=V_COBS,group=4)$Min_loss
  cat(rep_n,"run done.")
  rep_n=rep_n+1
}

summary_COBS=data.frame(matrix(NA,nrow=5,ncol=2),
                         row.names =c("PrinAngle","GrassDist","Min_loss","VFrobNorm","BFrobNorm"))
colnames(summary_COBS)=c("mean","sd")
summary_COBS[1,1]=round(mean(PA_COBS_V,na.rm=T),digits = 2)
summary_COBS[1,2]=round(sd(PA_COBS_V,na.rm=T),digits=2)
summary_COBS[2,1]=round(mean(Grass_COBS_V,na.rm=T),digits = 2)
summary_COBS[2,2]=round(sd(Grass_COBS_V,na.rm=T),digits=2)
summary_COBS[3,1]=round(mean(loss_COBS_V,na.rm=T),digits = 2)
summary_COBS[3,2]=round(sd(loss_COBS_V,na.rm=T),digits=2)
summary_COBS[4,1]=round(mean(Frob_COBS_V,na.rm=T),digits = 2)
summary_COBS[4,2]=round(sd(Frob_COBS_V,na.rm=T),digits=2)
summary_COBS[5,1]=round(mean(Frob_COBS_B,na.rm=T),digits = 2)
summary_COBS[5,2]=round(sd(Frob_COBS_B,na.rm=T),digits=2)

#save small n result (Capital words -- Settings, Lower case -- Method)
# write.table(summary_COBS,file="./Sim_Results_Sn/Summary_AJIVE_S_cobs.csv",sep=" ",col.names=T,row.names=T)

#save large n result (Capital words -- Settings, Lower case -- Method)
#write.table(summary_COBS,file="./Sim_Results_New/Summary_cobs_CBOS_L.csv",sep=" ",col.names=T,row.names=T)

# Use the multi-layer estimation for sigma_e^2 and diag Sigma_F
 se_sq=as.numeric(seEST.full(dataX=dataX,dataY=dataY,B=B_COBS, V=V_COBS)$sehat^2)
 sigma_k_est=0
 for (i in 1:4){
   sigma_k_est[i]=(sfEST.full(dataX=dataX,dataY=dataY,b=as.matrix(B_COBS[,i]),v=as.matrix(V_COBS[,i]),se=sqrt(se_sq))$sfhat)^2
 }
 sf_sq=sigma_k_est
 covf=diag(x=sf_sq)

 cov_est <- t(B_COBS)%*%cov(dataX)%*%B_COBS+covf
 cov_est
 diag(cov_est) # in descending order

