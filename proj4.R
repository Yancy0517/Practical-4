newt <- function(theta, func, grad, hess,..., tol=1e-8, fscale=1, maxit=100,
                 max.half=20, eps=1e-6){
  n <-length(theta)
  # mf <- 0
  # mtheta <- 0 
  iter <- 0
  # mg <- 0
  # mHi <-0
  if(iter == 0 && is.infinite(func(theta,...)) && is.infinite(grad(theta,...))){
    stop("The objective or derivatives are not finite at the initial theta")
  }
  
  while(iter <= maxit){
    f <- func(theta,...)
    g <- grad(theta,...)
    
    
    if(is.null(hess)){
      h <- approximate_hess(grad, theta, eps, n)
    }
    else{
      h <- hess(theta,...)
    }
    
    ec <- eigen(h)
    minvalue <- min(ec$values)
    while(minvalue <= 0){
      h <- h + diag(n)
      minvalue <- min(eigen(h)$values)
    }
    
    R <- chol(h)
    delta <- -chol2inv(R) %*% g
    
    half <- 0 
    while(half <= max.half){
      theta1 <- theta + delta/(2 ^ half)
      f1 <- func(theta1,...)
      if(f1 > f){
        half <- half + 1
      }
      else{
        break
      }
    }
    theta <- theta1
    
    if(half > max.half){
      stop("the step fails to reduce the objective despite trying ",max.half,
           " times step halvings")
    }
    
    g1 <- grad(theta,...)
}