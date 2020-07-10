seEST.full=function(dataX=dataX,dataY=dataY,B,V){
  r=length(B[1,])
  se_sqr=((norm(dataY-dataX%*%B%*%t(V),"F"))^2-sum((apply(dataY%*%V-dataX%*%B,2,norm.fcn))^2))/(n*(p-r))
  se_hat=sqrt(se_sqr)
  result <- list(sehat=se_hat)
  return(result)
}

sfEST.full=function(dataY=dataY,dataX=dataX,b,v,se){ #need to use full data Y in this
  sf_sqr=(norm(dataY%*%v-dataX%*%b,"F"))^2/n-(se^2)
  if (sf_sqr>0) {sf_hat=sqrt(sf_sqr)} else {sf_hat=0}
  result <- list(sfhat=sf_hat)
  return(result)
}

