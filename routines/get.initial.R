
get.initial <- function(X = X, Y_svd = Y_svd, Y = Y, 
                        choosePC = 1, 
                        init.v.given = NA,
                        param.b = list(alpha = 1, lambda = NA), 
                        grp.index.X = NA){
  
  # the optional input choosePC (default 1) enables choosing alternative initial values 
  
  if( is.na(init.v.given) ){
  uY <- Y_svd$u[,choosePC] * Y_svd$d[choosePC]      # initial score vector 
  v.init <- Y_svd$v[,choosePC]                # initial loading vector
  } else{
    v.init <- init.v.given
  }
  # initial values for b and sigmas are given by the routines 
  
  b.init <- update.b(X = X, Y = Y, v = v.init, param = param.b, 
                     grp.index.X = grp.index.X)$bhat 
  s.init <- unlist(Sigma.est(Y = Y, X = X, b.init, v.init))
  
  return(list(v.init = v.init, b.init = b.init, s.init = s.init))
}