cobs <- function(dataY, dataX = as.matrix(rep(0,nrow(dataY))), 
                 Y.bk.idx = rep(1,ncol(dataY)),  
                 X.bk.idx = rep(1,ncol(dataX)),
                 n.comp = 1,
                 alpha.b = 0,
                 lambda.b = 0.1,
                 alpha.v = 1,
                 lambda.v = 0.1,
                 tune.par.b = list(alpha = alpha.b, lambda = lambda.b),
                 tune.par.v = list(alpha = alpha.v, lambda = lambda.v),
                 orth = FALSE,
                 v.init = NULL,
                 maxLoop = 100,
                 TOL = .Machine$double.eps,
                 verbose = FALSE,
                 verbose.print.phrase = "cobs"
){
  # input: 
  #  dataY: n x p matrix of concatenated primary data blocks
  #  dataX: n x q matrix of concatenated covariate blocks 
  #         if there is no X, pass the zero vector.
  #  Y.bk.idx: p x 1 vector, indicating block memberships for Y
  #  X.bk.idx: q x 1 vector, indicating block memberships for X 
  #         if there is no group structure for X, treat as one group 
  #  n.comp: maximum number of components to be estimated
  #  lambda.b >= 0 (overall sparsity parameter for coefficient B)
  #  alpha.b in [0,1] (if 0, group lasso, if 1. lasso. If 1 (or 0), the (group) lasso parameter is lambda.b)
  #  lambda.v >= 0 (overall sparsity parameter for loading vector V)
  #  alpha.v in [0,1] (if 0, group-wise sparse, if 1, assume no group.)
  #   if tune.par.b is supplied, then inputs for lambda.b and alpha.b are ignored. 
  #   if tune.par.v is supplied, then inputs for lambda.v and alpha.v are ignored. 
  #  orth = FALSE (Flag for orthogonalization; TRUE or FALSE)
  #  maxLoop = 100 # maximum number of iterates, for each component
  #  TOL = .Machine$double.eps # A global parameter for convergence (~ 2e-16)
  #  verbose = TRUE # Prints out the progress (and prints out convergence inspection plots)
  #  verbose.print.phrase = "scarf" needed when exporting figures 
  # 
  # output:
  # est: a list containing all estimates
  # V  : matrix of loading vectors
  # B  : matrix of coefficients
  # sF2: vector of factor variances
  # se2: error variance 
  
  
  # COBS:COvariate-driven Blockwise Structured Factorization with Thresholding ------------------------------ 
  
  # Prepare output
  est <- list()
  
  # Preprocessing  -----------------------------------------------------------
  
  n <- nrow(dataY)
  p <- ncol(dataY)
  
  X <- as.matrix(dataX)
  q <- ncol(X)
  
  outFlag <- FALSE # A global parameter (if true, stop entirely)
  lambda.b.tuned <- NA 
  
  for (i.layer in 1:n.comp) {
    
    # Prepare the primary data, subtracting all previous components
    # (done by sequentially subtracting from Y)
    if (i.layer == 1){
      Y <- as.matrix(dataY) # centered
    }else{ 
     #  Y <- est[[i.layer - 1]]$Y - est[[i.layer - 1]]$Y %*% est[[i.layer-1]]$v%*% t(est[[i.layer-1]]$v)
       
      Y <-  est[[i.layer - 1]]$Y - ( est[[i.layer-1]]$s[2] * est[[i.layer - 1]]$Y %*% est[[i.layer-1]]$v  
                 + est[[i.layer-1]]$s[1] * X %*% est[[i.layer-1]]$b ) %*% 
        t(est[[i.layer-1]]$v) / (sum(est[[i.layer-1]]$s)) 
    }
    
    
    # Prepare for iterations --------------------------------------------------
    Y_svd <- svd(Y) 
    if (is.null(v.init)) {
      init <- get.initial(X = X, Y_svd = Y_svd, Y = Y,
                          param.b = tune.par.b, 
                          grp.index.X = X.bk.idx)
    }else{
      init <- get.initial(X = X, Y_svd = Y_svd, Y = Y,
                          init.v.given = v.init[,i.layer],
                          param.b = tune.par.b, 
                          grp.index.X = X.bk.idx)
    }
    b.now <- init$b.init
    v.now <- init$v.init
    s.now <- init$s.init
    convergeFlag <- FALSE # for this component
    cnt <- 1
    
    
    if(verbose){
      # This is only for inspection # comment out when not needed
      b.trace <- matrix(  nrow = q, ncol = maxLoop); b.trace[,1] <- b.now
      v.trace <- matrix(  nrow = p, ncol = maxLoop); v.trace[,1] <- v.now
      s.trace <- matrix( nrow = 2, ncol = maxLoop); s.trace[,1] <- s.now
    }
    
    # The inside loop ---------------------------------------------------------
    # Loop while outFlag == TRUE && convergeFlag == TRUE && cnt < maxLoop
     
    while( !outFlag  && !convergeFlag && cnt < maxLoop){
      
      # inside loop
      
      # Update V ----------------------------------------------------------------
      
      v.next  <- update.v(Y_svd = Y_svd, X = X, 
                          b = b.now, s = s.now, 
                          param = tune.par.v, 
                          grp.index = Y.bk.idx, Y = Y)
      
      # If V is a zero vector (resulting in v = NaN), then stop 
      if( any(is.nan(v.next)) ) {
        outFlag <- TRUE
        break;
      }
      
      if (sum(v.next * v.now) != 0 ){
        v.next <- sign(sum(v.next * v.now)) * v.next
      }
      
      
      # Update B ----------------------------------------------------------------
      b.tmp <- update.b(X = X, Y = Y, 
                        v = v.next,  
                        param = tune.par.b, 
                        grp.index.X = X.bk.idx)
      b.next <- b.tmp$bhat 
      
      if (is.na(tune.par.b$lambda)) { 
        lambda.b.tuned <- b.tmp$lambda
      } # only if the tuning parameter for B estimate was not supplied. 
      
      
      
      # Update Sigma.est --------------------------------------------------------
      s.next <- unlist(Sigma.est(Y = Y, X = X, b.next, v.next))
      
      b.diff <- sum( ( b.now - b.next )^2 ) 
      v.diff <- 1- abs(sum( v.now * v.next ) ) 
      if( b.diff < TOL && v.diff < TOL ){
        convergeFlag <- TRUE
      }
      b.now <- b.next 
      v.now <- v.next 
      s.now <- s.next 
      cnt <- cnt + 1 
      
      if(verbose){
        b.trace[,cnt] <- b.next 
        v.trace[,cnt] <- v.next 
        s.trace[,cnt] <- s.next 
      }
    }
    
    
    # Inspect the result
    if (verbose){
      if (outFlag){ 
        cat("Layer", i.layer, ":", 
            "loading vector is estimated zero. Stop execution.\n")
        break; 
      }else{
        cat("Layer", i.layer, ":",
            ifelse(convergeFlag,"Converged",
                   paste("not converged with v-diff", v.diff, "b-diff", b.diff)), 
            "at iteration", cnt, 
            "#zeros(v) = ", sum(v.now == 0) ,"\n")
        
        par( mfrow = c(2,2) )
        iterates <- max(1,cnt-10):cnt  
        plot(iterates, s.trace[1,iterates],"l", main = "se^2")
        plot(iterates, s.trace[2,iterates],"l", main = "sf^2")
        iterates <- max(1,cnt-10):(cnt-1)
        bdiff = colSums( ( b.trace[,1:(cnt-1)] - b.trace[,2:cnt] )^2 )
        plot(iterates, log10(bdiff[iterates]),"l", main = "B difference")
        vdiff = colSums( ( v.trace[,1:(cnt-1)] - v.trace[,2:cnt] )^2 )
        vdiff =  1 - colSums( (v.trace[,1:(cnt-1)] * v.trace[,2:cnt] ) )
        plot(iterates, log10(vdiff[iterates]),"l", main = "V difference")
        dev.copy(png,paste(verbose.print.phrase, '_layer',i.layer,'.png', sep=""))
        dev.off()
        
      } 
    }
    
    # Orthogonalize 
    v.orig <- v.now  
    if (orth){
      if (i.layer > 1){
        V.prev <- matrix(nrow = p, ncol = i.layer-1)
        for (j in (1:(i.layer-1))) { V.prev[,j] <- est[[j]]$v } 
        v.orth <- v.now - V.prev %*% t(V.prev) %*%v.now
        v.now <- v.orth / sqrt(sum(v.orth*v.orth))  
      }
    }
    
    
    # Prediction of F 
    f.now <- s.now[2]/sum(s.now)*(Y%*%v.now-X%*%b.now)
    
    # Save the result at the "est" list
    est[[i.layer]] <- list(b = b.now, v = v.now, s = s.now, 
                           f = f.now, 
                           v.orig = v.orig,
                           lambda.b.tuned = lambda.b.tuned, Y = Y)
    
    
    # End of inside loop ------------------------------------------------------
  }
  
  
  # Factor and error variance re-estimation ---------------------------------
  
  r <- length(est) # final number of components estimated
  V.est <- matrix(nrow = p, ncol = r)
  B.est <- matrix(nrow = q, ncol = r)
  F.est <- matrix(nrow = n, ncol = r)
  for ( j in 1:r){
    V.est[,j] <- est[[j]]$v
    B.est[,j] <- est[[j]]$b
    F.est[,j] <- est[[j]]$f
  }
   
  
  qV <- svd(V.est)$u
  Y <- as.matrix(dataY)
  se2.est <- norm( Y -   Y %*% qV %*%t(qV), "F")^2 / (n * (p-r))
  
  sF2.est <- vector()
  for ( j in 1:r){
    sF2.est[j] <- norm(Y %*% V.est[,j] - X %*% B.est[,j], "2")^2 / n  - se2.est
  }
  
  # Prepare output list -----------------------------------------------------
  return.list <- list(V = V.est, B = B.est, 
                      F.est=F.est,
                      sF2 = sF2.est, 
                      se2 = se2.est,
                      est = est,
                      Y.bk.idx = Y.bk.idx, 
                      X.bk.idx = X.bk.idx)
  return(return.list)
  
}



