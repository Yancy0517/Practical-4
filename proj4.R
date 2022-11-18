####---------- Group members ---------------------------------------------------
## Linsheng Shu (s2317223)
## Jialong He (s2281875)
## Heyu Nie (s2404675)

####---------- Github Repo -----------------------------------------------------
## https://github.com/Yancy0517/Practical-4.git

####------- Group member contribution ------------------------------------------
## Linsheng finishes the major construction of the optimization function "newt".
## Heyu addes errors and warnings part and some comments. (30%)
## Jialong writes the function of the hessian matrix when it is not provided
## and some comments. (30%)

################################################################################
# The "newt" function below is basically for optimization problems 
# or to minimize a given objective function.
# The idea is to implement  Newton’s method. Increase our initial values from
# "theta" by some tiny step and see if the objective function decreases.
# Repeat the process several times and finally reach the minimum.
################################################################################

# Function "newt" has several arguments.
# "theta" is initial values for the optimization parameters; 
# "func" is the objective function; "grad" is the gradient function
# "hess" is the Hessian matrix function or 2nd derivatives for each parameter
# "..." is any arguments of "func"; "grad"; "hess" after the first argument
# "tol" is the convergence tolerance; "fscale" used in convergence testing;
# "maxit" is the maximum number of iterations; "max.half" is the maximum number
# of times a step should be halved; "eps" is the finite difference interval
# to use when "hess" is not provided.

# "newt" finally return a list including minimal value of "func", the parameters
# the number of iterations, the gradient vector and the inverse matrix of "hess" 
# at that minimum.
newt <- function(theta, func, grad, hess,..., tol=1e-8, fscale=1, maxit=100,
                 max.half=20, eps=1e-6){
  n <-length(theta) # number of parameters
  iter <- 0 # iteration starts from 0
  
  # if the objective or derivatives are not finite at the initial "theta",
  # stop everything and raise errors
  if(is.infinite(func(theta,...)) | any(is.infinite(grad(theta,...)))){
    stop("The objective or derivatives are not finite at the initial theta")
  }
  
  # keep looking for the minimal value as if the number of iterations 
  # did not reach the maximum number 100
  while(iter <= maxit){
    f <- func(theta, ...) # the objective value
    g <- grad(theta, ...) # the gradient
    
    # if "hess" not provided, we obtain one through a new function
    if(is.null(hess)){
      h <- approximate_hess(grad, theta, eps, n) # see the second function
    }
    else{
      h <- hess(theta, ...) # if "hess" provided 
    }
    
    ec <- eigen(h) # eigen decomposition
    minvalue <- min(ec$values) # the minimal eigenvalues
    
    # if the minimal eigenvalues are less than or equal to 0, the hessian matrix
    # is not positive definite and we need to perturb it to be so
    while(minvalue <= 0){
      h <- h + diag(n) # add an identity matrix to it
      minvalue <- min(eigen(h)$values) # update the minimal eigenvalue
    }
    
    R <- chol(h) # cholesky decomposition
    delta <- -chol2inv(R) %*% g # the tiny step which decreases the objective
    
    half <- 0 # initialize the number of times the step should be halved
    while(half <= max.half){
      theta1 <- theta + delta/(2 ^ half) # update parameters
      f1 <- func(theta1,...) # new objective
      
      # compare this two objectives 
      if(f1 > f){ # if not decreasing the objective, halve the step, "half"+1
        half <- half + 1  
      }
      else{ # if objective decreases, no need to halve the step
        break
      }
    }
    # theta <- theta1
    
    # if the number of times the step is halved exceed its maximum, 
    # stop the process and raise error
    if(half > max.half){ 
      stop("the step fails to reduce the objective despite trying ",max.half,
           " times step halvings")
    }
    
    g1 <- grad(theta1, ...) # if not exceed, update gradient with new "theta"
    
    if(is.null(hess)){ # if no hessian provided
      # update "hess" with our defined function
      h1 <- approximate_hess(grad, ..., theta1, eps, n) 
    }
    else{
      h1 <- hess(theta1, ...) # or use function "hess" to update
    }
    
    ec <- eigen(h1) # eigen decomposition for our new hessian matrix
    minvalue <- min(ec$values) # the minimal eigenvalue
    # if the convergence reaches
    if(max(abs(g1)) < tol * (abs(f1) + fscale) ){
      if(minvalue <= 0){ # also if hessian is not +ve definite, stop procedure
        stop("The Hessian is not positive definite at convergence")
      }
      else{ # if hessian is +ve definite, return required list finally
        final = list('f' = f1, 'theta' = theta, 'iter' = iter, 'g' = g1, 
                     'Hi' = chol2inv(chol(h1)))
        return(final)
      }
    }

    iter <- iter + 1 # at the end of each loop, + 1 iterations
  }
  
  # when iterations greater than its maximum, raise errors.
  if(iter > maxit){
    stop("The maximum iteration(", maxit, ") is reached without convergence")
  }
}
################################################################################
# If the Hessian matrix is not provided, we need to manually obtain an
# approximation to it by finite differencing of gradient vector which is
# returned by function "grad".
# Function "approximate_hess" below is created for this case.
approximate_hess <- function(grad,..., theta, eps, n){
  grad0 <- grad(theta, ...) # original gradient 
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
################################################################################