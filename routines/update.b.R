update.b <- function(X = X, Y = Y, v = v.now, 
                     param = list(alpha = 1, lambda = NA), 
                     grp.index.X) {
   
  # lambda.b > 0 (overall sparsity parameter for coefficient B)
  # alpha.b in (0,1) 
  # (if 0, group lasso, if 1. lasso. 
  # If 1 (or 0), the (group) lasso parameter is lambda.b)
  # 
  # update b by minimizing || y - xb || + penalty(b). 
  # assuming there is no intercept (when written as a regression y ~ x)
  # if an intercept is needed, explicitly form the x matrix accordingly
  # 
  # Rewrite as a regression problem
  x <- X
  y <- Y %*% v
  
  # Two options, 
  # 1. Lambda is supplied.      (Lambda.choice == FALSE)
  # 2. Lambda is chosen inside. (Lambda.choice == TRUE)
  Lambda.choice <- is.na(param$lambda)
  
  # three cases: 
  # 1. There is no penalty (lambda = 0 or only one covariate in x)
  # 2. Only variable selection is needed (alpha == 1, Gaussian lasso using glmnet)
  # 3. Both variable and group selections are needed (alpha < 1, using SGL)
  
  
  # Case I  
  
  if (NCOL(x) == 1 || identical(param$lambda,0)) {
    bhat <- lm( y ~ x - 1)$coefficient 
    bhat <- ifelse(is.na(bhat),0,bhat)
  } else{
    if ( identical(param$alpha,1)  ){
      # Case II (Gaussian lasso using glmnet)
      # 
      require(glmnet)
      
      if (Lambda.choice){ # option I 
        
        a <- cv.glmnet(x,y, alpha = param$alpha, 
                       standardize = FALSE, intercept = FALSE)
        bhat <- coef(a, s = "lambda.1se")[2:(q+1)]
        param$lambda <- a$lambda.1se
        
      }else{ # option II
        a <-  glmnet(x,y, alpha = param$alpha, lambda = param$lambda,
                     standardize = FALSE, intercept = FALSE)
        bhat <- as.vector(a$beta)
      }
    }else{
      # Case III (using SGL, for the case: alpha < 1)
      require(SGL)
      #if (is.na(grp.index.X)) {
      #   grp.index.X <- rep(1,ncol(X))
      #}
      
      if (Lambda.choice){ # option I 
        a <- cvSGL(data = list(x = x, y = y), alpha = param$alpha,
                   index = grp.index.X, type = "linear", 
                   standardize = FALSE)
        l.min.index <- which.min(a$lldiff)
        l.1se.index <- which( a$lldiff - (a$lldiff[  l.min.index ] + a$llSD[  l.min.index ]) < 0 )
        bhat <- a$fit$beta[,l.1se.index[1]]
        param$lambda <- a$lambdas[l.1se.index[1]]
        }else{ # option II
        a <-  SGL(data = list(x = x, y = y), alpha = param$alpha, lambdas = param$lambda,
                  index = grp.index.X, type = "linear", 
                  standardize = FALSE) 
        bhat <- a$beta
      }
      
    }
  }
  return(list(bhat = bhat, lambda = param$lambda))
}
