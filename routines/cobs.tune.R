cobs.tune <- function(dataY, dataX = as.matrix(rep(0,nrow(dataY))), 
                 Y.bk.idx = rep(1,ncol(dataY)),  
                 X.bk.idx = rep(1,ncol(dataX)),
                 n.comp = 1,
                 alpha.b = 1,
                 lambda.b = 0.0025,
                 tune.par.b = list(alpha = alpha.b, lambda = lambda.b),
                 alpha.v = 0,
                 lambda.v.seq = NULL, 
                 orth = FALSE,
                 v.init = NULL,
                 maxLoop = 100,
                 TOL = .Machine$double.eps,
                 verbose = FALSE,
                 verbose.print.phrase = "cobs",
                 method = c("BIC","manual", "BIC-Low"),
                 M.lambda.v = 20 ){
  # Tuning for lambda, where alpha is fixed. 
  
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
      
      Y <-  est[[i.layer - 1]]$Y - ( est[[i.layer-1]]$s[2] * est[[i.layer - 1]]$Y %*% est[[i.layer-1]]$v+ 
            est[[i.layer-1]]$s[1] * X %*% est[[i.layer-1]]$b )%*% 
            t(est[[i.layer-1]]$v) / (sum(est[[i.layer-1]]$s)) 
    }
    
    
    
    # Find max lambda v -------------------------------------------------------
     if (is.null( lambda.v.seq) ){
    tout <- svd(Y, nv = 1, nu = 0)
    lvmax_variable <- max(abs(tout$v)) / alpha.v
    require(dplyr)
    lvmax_group <- data.frame(v = tout$v, label = Y.bk.idx) %>% group_by(label) %>% 
      summarize(sqnorm = sum(v^2)) %>% select(sqnorm) %>% max()
    lvmax <- min(lvmax_group / (1-alpha.v),lvmax_variable) - 10e-3
    K = length(unique(Y.bk.idx))
   
    lambda.v.seq = lvmax*0.9^(1:M.lambda.v-1)    
     }
    
    # Now tune ----------------------------------------------------------------

    bic = vector(length = length(lambda.v.seq))
    Vmat = matrix(nrow = K, ncol = length(lambda.v.seq))
    cobs.out.list <- list()
    cat(paste("Layer",i.layer,"tuning:"))
    for (l in 1:length(lambda.v.seq)){
      cobs.out <- cobs(dataY = Y, dataX = dataX, 
                       Y.bk.idx = Y.bk.idx, X.bk.idx = X.bk.idx, 
                       tune.par.b = tune.par.b,
                       tune.par.v = list(alpha = alpha.v, lambda = lambda.v.seq[l]),
                       n.comp = 1)
      v.init <- cobs.out$V
      if (method == "BIC-Low"){
        bic.omega <- "LD" 
        } else {bic.omega <- "HD"}
      summary.out <- summary.cobs(cobs.out, verbose = FALSE, 
                                  X = dataX, Y = Y,
                                  BICchoice = bic.omega)
      Vmat[,l] <- summary.out$V.str
      bic[l] <- summary.out$bic
      cobs.out.list[[l]] <- cobs.out
      cat(l)
    }
    
    
    # Draw -------------------------------------------------------------------- 
    {
      require(ggplot2)
      require(tidyr)
      d <- data.frame( log_lambda = log(lambda.v.seq), 
                       bic = bic, 
                       nBlocks = colSums(ifelse(Vmat > 0, 1,0)))
      g <- d %>% 
        ggplot() + 
        geom_bar(mapping = aes(x = log_lambda, y = nBlocks), 
                 stat = "identity", color = gray(0.5), fill = gray(0.8)) + 
        geom_line(aes(x = log_lambda, y = K * (bic -min(bic) ) / diff(range(bic)))) + 
        scale_y_continuous(
          name = expression("# blocks identified"), 
          sec.axis = sec_axis(~ . * diff(range(bic)) / K + min(bic) , name = "BIC")) 
      
      
      uBlocks <- unique(d$nBlocks)
      gid <- uBlocks
      tdata <- NULL
      for (kk in 1:length(uBlocks)){ 
        gid[kk] <- which(min(d$bic[d$nBlocks == uBlocks[kk]]) == d$bic)[1]
        
        g <- g + geom_vline(xintercept = log( lambda.v.seq[ gid[kk] ]), colour="red")
        
        tdata <- rbind(tdata, data.frame(x = log( lambda.v.seq[ gid[kk] ]),
                                         y = d$nBlocks[ gid[kk] ]-0.5,
                                         label  = paste("lambda id = ",  gid[kk] ) )
        )
        
      }
      
      g <- g + geom_text(aes(x = x, y = y, label = label), 
                         angle = 90, vjust = 1, data = tdata) +
        theme_bw() + 
        theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()
        ) + ggtitle(paste("COBS Tuning for Layer ",i.layer))
      
      print(g)
    }
    if (method == "manual"){
    lambda_index <- readline(prompt="Choose your lambda index: ")} else {
      lambda_index <- which.min(bic)
    }
    
    if(verbose){
      dev.copy(png,paste('tune',verbose.print.phrase, '_layer',i.layer,'.png', sep=""))
      dev.off()
    }
    
    current.choice <- cobs.out.list[[as.numeric(lambda_index)]] 
    
    
    # Orthogonalize if needed
    v.now <- current.choice$est[[1]]$v
    v.orig <- v.now  
    if (orth){
      if (i.layer > 1){
        V.prev <- matrix(nrow = p, ncol = i.layer-1)
        for (j in (1:(i.layer-1))) { V.prev[,j] <- est[[j]]$v } 
        v.orth <- v.now - V.prev %*% t(V.prev) %*%v.now
        v.now <- v.orth / sqrt(sum(v.orth*v.orth))  
      }
    }
    b.now <- current.choice$est[[1]]$b
    s.now <- current.choice$est[[1]]$s
    f.now <- s.now[2]/sum(s.now)*(Y%*%v.now-X%*%b.now)
    
    # Save the result at the "est" list
    current.choice$est[[1]]$lambda.v.tuned <- lambda.v.seq[lambda_index] 
    current.choice$est[[1]]$f <- f.now
    current.choice$est[[1]]$v <- v.now
    current.choice$est[[1]]$v.orig <- v.orig
    
    est[[i.layer]] <- current.choice$est[[1]]
    
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
  return.list <- list(V = V.est, B = B.est, F.est=F.est,
                      sF2 = sF2.est, 
                      se2 = se2.est,
                      est = est,
                      Y.bk.idx = Y.bk.idx, 
                      X.bk.idx = X.bk.idx)
  return(return.list)
}



