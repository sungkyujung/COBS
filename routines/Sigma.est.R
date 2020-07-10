Sigma.est <- function(Y = Y, X = X, b, v){  
  
  n <- nrow(Y)
  p <- ncol(Y)
  error2 <- norm( Y %*% v - X %*% b , "2")^2 
  se2 <- ( norm( Y - X %*% b %*% t(v) , "F")^2 -  error2  ) / ( n * (p-1))
  sf2 <- max(error2/n - se2 ,0)
  list(se2 = se2, sf2  = sf2)
  # returns sigma^2_e and sigma^2_f esimates
} 