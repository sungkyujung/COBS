update.v  <- function(Y_svd = Y_svd, X = X, b = b.now, s = s.now,  
                      param = list(alpha = 0, lambda = 0), 
                      grp.index = NULL, Y=Y) {
   
   
  # lambda > 0 (overall sparsity parameter for loading vector v)
  # alpha in (0,1) controls the degree to which group-sparsity is applied. 
  # if alpha = 0, the strictly group-wise sparse estimate.
  # if alpha = 1, the group structure is ignored. 
  
  
  # Two cases: 
  #  - Case I: sigma_f > 0 
  #  - Case II: sigma_f == 0 
  
  
  # Step 1: Exact max-norm --------------------------------------------------
  
  if (s[2] > 0 ){ 
    #  - Case I: sigma_f > 0 
    d <- Y_svd$d 
    if (norm(b,"2") == 0){ 
      v.cand <- Y_svd$v[,1]
    }else{
      z <- s[1]/s[2] * t(Y_svd $u) %*%  X %*% b
      
      h <- function(t) {
        aa <- vector()
        for (i in 1:length(t)){ 
          aa[i] <- sum ( d^2 * z^2 / (d^2 - t[i])^2  )-1
        }
        return(aa)
      }
      
      t1 <- d[1]^2 + abs(d[1]*z[1])
      t2 <- d[1]^2 + sqrt(sum(d^2 * z^2)) 
      tcand <- bisect.dec(h, l = d[1]^2, t1 = t1, t2 = t2)   
      v.cand <- Y_svd $v %*%  ( d / (tcand - d^2) * z ) 
        
      }
    
      
  }else{ 
    #  - Case II: sigma_f == 0 
    
    z <-    X %*% b 
    zY <- t(z) %*% Y 
    n_phi.cand <-  sum(z^2) + sqrt( zY%*% t(zY) ) 
    v.cand <- list()
    v.cand[[1]] <- zY / as.numeric(sum(z^2) - n_phi.cand )
    
    n_phi.cand <-  sum(z^2) - sqrt( zY%*% t(zY) ) 
    v.cand[[2]] <- zY / as.numeric(sum(z^2) - n_phi.cand )
    
    min.id <- which.min (c(norm(Y - z %*% v.cand[[1]],"F")^2,   norm(Y - z %*% v.cand[[2]],"F")^2)) 
    v.cand <- as.vector(v.cand[[min.id]])
    
  }
  
  
  
  
  # Step 2: Thresholding ---------------------------------------------------- 
  
  thres.group <- (1 - param$alpha) * param$lambda
  thres.variable <- param$alpha * param$lambda
  
  v.cand_thres <- pmax( abs(v.cand) - thres.variable , 0)
  if (sum(v.cand_thres^2) < .Machine$double.eps^2){
    v <- NaN
  }else{
    
    u.index <- unique(grp.index)
    for (k in u.index){
      index.k <- grp.index == u.index[k]
      a_k <- as.matrix(v.cand_thres[index.k])
      norm.a_k <- norm(a_k,"F")
      if (norm.a_k>0){
        v.cand_thres[index.k] <- a_k * max( norm.a_k - thres.group, 0 ) / norm.a_k 
      }
      
    }
    
    v <- sign(v.cand) * v.cand_thres / norm(v.cand_thres, "2")
  }
  if (all(is.na(v))){
    v <- NaN
  }else{
      if (sum(v^2) < .Machine$double.eps^2){
      v <- NaN
    } 
  }
  
  
  return(v)
}












