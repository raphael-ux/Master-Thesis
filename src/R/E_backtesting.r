#rewrite e_Q in our case 1-p and -r and -z
e_Q <- function(L,r,z = NULL,p){
    I <- as.numeric(L > r)
    out <- (1/(1-p))*I
    return(out)
}

e_ES <- function(L,r,z,p){
    out <- pmax(0,L-z)/((1-p)*(r-z))
    return(out)
}

M_process <- function(L,r,z,e,lambda,p){
    T <- length(L)
    out <- numeric(T)
    X <- 1 - lambda + lambda*e(L,r,z,p)
    out <- cumprod(X)  
    return(out)
}

lambda_GREE <- function(L, r, z, e, p) {
  T <- length(L)
  lambda_star <- numeric(T)
  lambda_star[1] <- 0
  val <- e(L, r, z, p)
  for (t in 2:T) {
    g <- function(lambda) {
      mean(log(1 - lambda + lambda * val[1:(t - 1)]))
    }
    opt <- optimize(g, interval = c(0, 0.5), maximum = TRUE)
    lambda_star[t] <- opt$maximum
  }
  return(lambda_star)
}

lambda_GREL <- function(L, r, z, e, p) {
  T <- length(L)
  lambda_star <- numeric(T)
  lambda_star[1] <- 0
  
  for (t in 2:T) {
    val <- e(L[1:(t-1)], r[t], z[t], p)
    g <- function(lambda) mean(log(1 - lambda + lambda * val))
    opt <- optimize(g, interval = c(0, 0.5), maximum = TRUE)
    lambda_star[t] <- opt$maximum
  }
  
  return(lambda_star)
}

lambda_GREM <- function(L,r,z,e,p){
    T <- length(L)

    lambda_GE <- lambda_GREE(L,r,z,e,p)
    lambda_GL <- lambda_GREL(L,r,z,e,p)

    lambda_star <- numeric(T)
    lambda_star[1] <- 0

    M_GREE <- M_process(L,r,z,e,lambda_GE,p)
    M_GREL <- M_process(L,r,z,e,lambda_GL,p)

    for(t in 2:T){
        lambda_star[t] <- (M_GREE[t-1]*lambda_GE[t] + M_GREL[t-1]*lambda_GL[t])/(M_GREE[t-1] + M_GREL[t-1])
    }
    return(lambda_star)
}

e_backtest <- function(L, r, z, e, lambda = "GREM", threshold, p){
    
    if(lambda == "GREE"){
        lambda <- lambda_GREE(L, r, z, e, p)
        M <- M_process(L, r, z, e, lambda, p)
    } else if(lambda == "GREL"){
        lambda <- lambda_GREL(L, r, z, e, p)
        M <- M_process(L, r, z, e, lambda, p)
    } else if(lambda == "GREM"){
        print("test1")
        lambda <- lambda_GREM(L, r, z, e, p)
        M <- M_process(L, r, z, e, lambda, p)
    } else {
        stop("lambda must be 'GREE', 'GREL', or 'GREM'")
    }
    idx <- which(M >= 1/threshold)
    return(length(idx)/length(M))
}
