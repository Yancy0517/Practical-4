####---------- Group members ---------------------------------------------------
## Linsheng Shu (s2317223)
## Jialong He (s2281875)
## Heyu Nie (s2404675)

####---------- Github Repo -----------------------------------------------------
## https://github.com/Yancy0517/Practical-4.git

####------- Group member contribution ------------------------------------------
##


newt <- function(theta, func, grad, hess,..., tol=1e-8, fscale=1, maxit=100,
                 max.half=20, eps=1e-6){
  n <-length(theta) # number of optimization parameters
  iter <- 0 # iteration times starts from 0

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
    
    if(is.null(hess)){
      h1 <- approximate_hess(grad, theta, eps, n)
    }
    else{
      h1 <- hess(theta,...)
    }
    
    ec <- eigen(h1)
    minvalue <- min(ec$values)
    if(max(abs(g1)) < tol * (abs(f1) + fscale) ){
      if(minvalue <= 0){
        stop("The Hessian is not positive definite at convergence")
      }
      else{
        final = list('f' = f1, 'theta' = theta, 'iter' = iter, 'g' = g1, 
                     'Hi' = chol2inv(chol(h1)))
        return(final)
      }
    }

    iter <- iter + 1 
  }
  
  if(iter > maxit){
    stop("The maximum iteration(", maxit, ") is reached without convergence")
  }
  
  final = list('f' = mf, 'theta' = mtheta, 'iter' = iter, 
               'g' = mg, 'Hi' = hi)
  return(final)
}

# If the Hessian matrix is not provided, we need to manually obtain an
# approximation to it by finite differencing of gradient vector which is
# returned by function "grad".
approximate_hess <- function(grad, theta, eps, n){
  grad0 <- grad(theta) # original gradient 
  h <- matrix(0, n, n) # an empty matrix for Hessian
  
  # loop over parameters
  for(i in 1:n){
    # increase i-th parameter by "eps" which is the finite difference interval
    th1 <- theta; th1[i] <- th1[i] + eps  
    grad1 <- grad(th1) # obtain new gradient from updated parameter
    h[i,] <- (grad1 - grad0)/eps # approximate second derivative
  }
  (t(h) + h) / 2 # make it exactly symmetric
}