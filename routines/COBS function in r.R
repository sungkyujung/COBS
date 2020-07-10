
norm.fcn=function(a){
  norm(as.matrix(a),"F")
}

##### COBS function with BIC , only tune lambda.b and lambda.v so far#####
# Used in single layer
fbic <- function(X=X,Y=Y,b=b_hat,v=v_hat,se=se_hat,sf=sf_hat){
  kb=sum(b!=0)
  kv=sum(v!=0)
  if(kv==0){bic=NA} else {
    bic=(log(n*p)*(kb+kv)+n*(log(2*pi)+2*(p-1)*log(se)+log(sf^2+se^2))+
           1/(se^2)*(norm(Y-X%*%b%*%t(v),"F")^2)-sf^2/(se^2*(sf^2+se^2))*(norm(Y%*%v-X%*%b,"F")^2))/(n*p)
    #bic=(log(n*p)*(kb+kv)+1/(se^2)*(norm(Y-X%*%b%*%t(v),"F")^2)-sf^2/(se^2*(sf^2+se^2))*(norm(Y%*%v-X%*%b,"F")^2))/(n*p)
  }
  result <- list(bic=bic)
  return(result)
}

# Gram-Schmidt Process Projection function: Project vec2 onto vec1
Proj <- function (vec1, vec2){
  if(all(vec1==0)){
    as.matrix(rep(0,length(vec2)))
  } else{
    vec1=as.matrix(vec1)
    vec2=as.matrix(vec2)
    as.matrix(as.numeric(t(vec1)%*%vec2)/(norm(vec1,"F")^2)*vec1)}
}

##### Evaluation functions #####
# Angle
angle<-function(x,true.value){
  #x=as.matrix(x)
  #true.value=as.matrix(true.value)
  acos(sum(x*true.value)/sqrt(sum(x^2))/sqrt(sum(true.value^2)))/pi*180
  #acos(t(x)%*%true.value/(norm(x,"F")*norm(true.value,"F")))/pi*180
  #round(angle,digits=3)
}

# Principal Angle
PrinAngle <- function (M1 , M2){
  U1=as.matrix(svd(M1)$u)  # can also use Gram-Schmidt to get unitary matrix 
  U2=as.matrix(svd(M2)$u)  # can also use Gram-Schmidt to get unitary matrix
  PA=180/pi*acos(svd(t(U1)%*%U2)$d)
  Max.PA=180/pi*acos(min(svd(t(U1)%*%U2)$d))
  result <- list(Angle=PA,Max_angle=Max.PA)
  return(result)
}

# Grassmannian Metric
GrassDist <- function(M1, M2){
  U1=as.matrix(svd(M1)$u) # can also use Gram-Schmidt to get unitary matrix 
  U2=as.matrix(svd(M2)$u) # can also use Gram-Schmidt to get unitary matrix 
  D1=sqrt(sum(acos(svd(t(U1)%*%U2)$d)^2))
  result<-list(Grass_Dist=D1)
  return(result)
}

# Loss function to evaluate group structure in V
lossV=function(M_true, M_est, group){
  tal=10^(-4)
  M1=as.matrix(M_true); M2=as.matrix(M_est)
  p1=dim(M1)[1]; p2=dim(M2)[1]; p=p1
  if (p1 != p2){ stop("Dimension of matrices do not match! exit...")}
  r1=dim(M1)[2]; r2=dim(M2)[2]; 
  size.group=p/group;  gpM=ceiling(1:p/ size.group)
  K1=matrix(NA,nrow=group,ncol=r1); K2=matrix(NA,nrow=group,ncol=r2)
  for (icol in 1:r1){
    for (irow in 1:group){
      K1[irow,icol]=ifelse(sum(abs(M1[which(gpM==irow),icol]))<=tal,
                           #M1[which(gpM==irow),icol]==0)==size.group,
                           0,1)
    }
  }
  for (icol in 1:r2){
    for (irow in 1:group){
      K2[irow,icol]=ifelse(sum(abs(M2[which(gpM==irow),icol]))<=tal,
                           #M2[which(gpM==irow),icol]==0)==size.group,
                           0,1)
    }
  }
  if (r1 <= r2){ A1 = K1; A2=K2} else {A1=K2; A2=K1} #A1 is the one with smaller dimension
  per=gtools::permutations(n=max(r1,r2),r=max(r1,r2))
  loss=rep(NA, length(per[,1]))
  for (l in 1:length(per[,1])){
    P_A2=A2%*%(diag(max(r1,r2))[,per[l,]])
    loss[l]=sum(abs(A1-P_A2[,1:min(r1,r2)]))+group*(max(r1,r2)-min(r1,r2))
  }
  min.loss=min(loss)
  max.loss=max(loss)
  min.permute=per[which(loss==min.loss),]
  result<-list(Min.permute = min.permute, Min_loss=min.loss, Max_loss=max.loss)
  return(result)
}

