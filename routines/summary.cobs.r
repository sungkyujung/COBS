summary.cobs <- function(cobs.out,verbose = TRUE, 
                         BICchoice = c("HD","LD"), 
                         X=NULL, # Input X to output summary of covariates, otherwise output is NA
                         Y=Y){
  
  p <- length(cobs.out$Y.bk.idx)
  q <- length(cobs.out$X.bk.idx)
  r <- dim(cobs.out$V)[2]
  
  Y.idx <- cobs.out$Y.bk.idx
  uniq.Y.idx <- unique(Y.idx)
  K <- length(uniq.Y.idx)
  
  X.idx <- cobs.out$X.bk.idx
  uniq.X.idx <- unique(X.idx)
  G <- length(uniq.X.idx)
  
  # Block-wise segmentation (V)
  V.str <- matrix(nrow = K, ncol = r)
  for( j in 1:K){
    if (r > 1){
    V.str[j,] <- colSums( (cobs.out$V[Y.idx == uniq.Y.idx[j], ])^2 ) }
    else{
      V.str[j] <- sum( (cobs.out$V[Y.idx == uniq.Y.idx[j], ])^2 )}
     
  }
  colnames(V.str) <- paste("Comp",1:r)
  rownames(V.str) <- uniq.Y.idx 
  
  # Block-wise segmentation (B)
  B.str <- matrix(nrow = G, ncol = r)
  if (all(X==0)==0 & is.null(X)==0) {
    for( j in 1:G){
      if (r > 1){
        B.str[j,] <- colSums( (cobs.out$B[X.idx == uniq.X.idx[j], ])^2 ) }
      else{
        B.str[j] <- sum( (cobs.out$B[X.idx == uniq.X.idx[j], ])^2 )}
      
    }
    colnames(B.str) <- paste("Comp",1:r)
    rownames(B.str) <- uniq.X.idx 
  }
  
  
  # Variance 
  cov_var=NA
  if (all(X==0)==0 & is.null(X)==0){
    cov_var=diag(t(cobs.out$B) %*% cov(dataX) %*% cobs.out$B)
  }
  variations <- rbind(cov_var , cobs.out$sF2)
  colnames(variations) <- paste("Comp",1:r)
  rownames(variations) <- c("Covariate","Unknown")
  
  bic  <- NA
  if (r == 1){ # if single layer 
    # Used in single layer 
      kb <- sum(cobs.out$B!=0)
      kv <- sum(cobs.out$V!=0)
      if( kv != 0) {
        se <- sqrt( cobs.out$se2 ) 
        sf <- sqrt( cobs.out$sF2 )
        p <- ncol(Y)
        n <- ncol(Y)
        b <- cobs.out$B
        v <- cobs.out$V
        if (BICchoice[1] == "HD"){
        gamma <- 0.1
        omega <- 6 * (1+gamma)
        }else{
          omega <- 1
        }
        bic=(omega*log(n*p)*(kb+kv)+n*(log(2*pi)+2*(p-1)*log(se)+log(sf^2+se^2)) +
               1/(se^2)*(norm(Y-X%*%b%*%t(v),"F")^2) - 
               sf^2/(se^2*(sf^2+se^2))*(norm(Y%*%v-X%*%b,"F")^2))/(n*p) 
        
        
      } 
    }
  
  
  if (verbose){
  cat(    "Covariate-driven Blockwise Structured Factorization \n")
  cat(" has estimated",r,"components.\n")
  cat(" with", K, "blocks of primary data and ",G,"variable groups in covariate\n")
  cat("structure in V:\n")
  print(V.str)
  cat("structure in B:\n")
  print(B.str)
  cat("Variation structure decomposed, with overall var(e) =", cobs.out$se2,",\n")
  print(variations)
  }
  else{
  }
  return(list(V.str = V.str, B.str = B.str, variations = variations, bic = bic))
  
    
}
