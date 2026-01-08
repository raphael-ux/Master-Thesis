library("GAS")
library(DEoptim)
data("sp500ret")
library(nloptr)
library(data.table)
library(parallel)
library(numDeriv)
library(pbapply)


UniFilter_hybridmodel <- function(vPw, in_sample_data, out_sample_data = NULL, alphas = c(0.01, 0.05), forecast = FALSE) {
  in_sample_data <- as.numeric(in_sample_data)
  if(!is.null(out_sample_data)) out_sample_data <- as.numeric(out_sample_data)

  data <- if(forecast) {
    stopifnot(!is.null(out_sample_data))
    c(in_sample_data,out_sample_data)
  } else {
    in_sample_data
  }
  n <- length(alphas)
  iK <- 2 * n
  if(forecast){
    iH <- length(out_sample_data)
    iT <- length(data)
  }else{
    iT <- length(data)
  }

  A <- vPw[1]
  B <- vPw[2]

  a <- vPw[seq(3, by = 2, length.out = n)]
  b <- vPw[seq(4, by = 2, length.out = n)]

  delta <- vPw[length(vPw)]
  
  states <- matrix(0, nrow = iK, ncol = iT)
  
  kt <- numeric(iT)
  kt[1] <- 0
  states[,1] <- c(rbind(a, b)) * exp(kt[1])   

  alphas_mat <- matrix(alphas, nrow = n, ncol = 1)

  v_idx <- seq(1, iK, by=2)
  e_idx <- seq(2, iK, by=2)
  
  for(t in 2:iT) {
    y <- data[t-1]
    
    prev_state <-  states[, t-1]
    v <- prev_state[v_idx]
    e <- prev_state[e_idx]
    
    I <- as.numeric(y <= v)
    lambda <- (1 / e) * ((1 / alphas) * I * y - e)  
    
    current_kt <- kt[t-1]
    kt[t] <- B * current_kt + A * sum(lambda) + delta * log(abs(y) + 1e-8)
    
    states[v_idx, t] <- a * exp(kt[t])
    states[e_idx, t] <- b * exp(kt[t])
  }
  if(forecast){
    states_list <- split(states[,(iT-iH + 1):iT], row(states[,(iT-iH + 1):iT]))
    return(states_list)
  }else{
    states_list <- split(states, row(states))
    return(states_list)
  }
}

UniGASOptimiser_FLZ_hybridmodel <- function(vPw, data, alphas = c(0.01,0.05)) {
    iK <-2*length(alphas)
    iT <- length(data)
    vY <- as.numeric(data)
    res <- UniFilter_hybridmodel(vPw = vPw , in_sample_data = data  , alphas = alphas) 
    lossMat <- matrix(NA, nrow = iT, ncol = length(alphas))
    for (j in 1:length(alphas)){
        v <- res[[2*j-1]]
        e <- res[[2*j]]
        lossMat[, j] <-  FZLoss(data = data, VaR = v, ES = e, alpha = alphas[j])
  }
  overallLoss <- mean(rowSums(lossMat), na.rm = TRUE)
  return(overallLoss)
}



UniFit_hybridmodel <- function(data, alphas = c(0.01,0.05), fn.optimizer = fn.nloptr_hybrid.optim   ) {
    vY <- data
    iT <- length(vY)
    iK <- 2*length(alphas)
    
    va <-quantile(data, alphas) 
    vb <- sapply(seq_along(alphas), function(i) {
        mean(data[data <= va[i]], na.rm = TRUE)
    })
    vab <- as.vector(rbind(va, vb))  
    vPw <- c(0.01,0.98,vab,0.001)
    optimiser <- fn.optimizer(par0 = vPw,data = data,alphas = alphas ,FUN = UniGASOptimiser_FLZ_hybridmodel)

    if (optimiser$convergence != 0) {
        warning(paste("solver failed to converge. Convergence flag:", optimiser$convergence))
    }
    return(optimiser)
}



UniRoll_hybridmodel <- function(data, alphas = c(0.01), ForecastLength = 50L, Nstart = NULL, RefitEvery = 2L, RefitWindow = c("moving","recursive"), cluster = NULL, Compute.SE = FALSE){
    vY <- data
    iT <- length(vY)
    iK <- 2*length(alphas)

    lFits <- list()
    lForecasts <- vector("list", iK)
    lpara <- list()
    lData <- list()
    lOut <- list()

    iStart <- iT - ForecastLength

    FitIndex <- seq(iStart, iT, RefitEvery)
    if (tail(FitIndex, 1L) == iT) {
        FitIndex = FitIndex[-length(FitIndex)]
    }

    if (RefitWindow[1L] == "recursive") {
        for (i in 1:length(FitIndex)) {
            lData[[i]] = vY[1:FitIndex[i]]
        }
    }
    if (RefitWindow[1L] == "moving") {
        for (i in 1:length(FitIndex)) {
            lData[[i]] = vY[(FitIndex[i] - iStart + 1L):FitIndex[i]]
        }
    }


    for (i in seq_along(FitIndex)) {
      if (i != length(FitIndex)) {
        lOut[[i]] <- vY[(FitIndex[i] + 1L):FitIndex[i + 1L]]
      } else {
        lOut[[i]] <- vY[(FitIndex[i] + 1L):length(vY)]
      }
    }

  if (is.null(cluster)) {
    lFits <- lapply(seq_along(lData), function(i) {message("Running fit ", i, " of ", length(lData));UniFit_hybridmodel(lData[[i]], alphas = alphas)})
  } else {
    clusterEvalQ(cluster, {library(GAS);library(nloptr)})
    clusterExport(cluster, c("UniFit_hybridmodel", "fn.nloptr_hybrid.optim","UniGASOptimiser_FLZ_hybridmodel","UniFilter_hybridmodel"))
    lFits <- parLapply(cluster,lData,UniFit_hybridmodel,alphas = alphas)
  }
  
  lpara  <- lapply(lFits, function(f) f$pars)

  For <- lapply(seq_along(lOut), function(i, lpara, lData, lOut) {
    UniFilter_hybridmodel(
      vPw = lpara[[i]],
      in_sample_data = lData[[i]],
      out_sample_data = lOut[[i]],
      alphas = alphas,
      forecast = TRUE
    )
  }, lpara = lpara, lData = lData, lOut = lOut)

  for (k in 1:iK) {
    lForecasts[[k]] <- unlist(lapply(For, function(f) f[[k]]))
  }
  return(list(Forecasts = lForecasts, Parameters = lpara))
}


fn.nloptr_hybrid.optim <- function(par0, data, alphas, FUN) {
  
  obj_fun <- function(par) {
    return(FUN(par, data = data, alphas = alphas))
  }

  ineq_constraint <- function(par) {
      n <- length(alphas)
      a <- par[seq(3, by = 2, length.out = n)]
      b <- par[seq(4, by = 2, length.out = n)]
      return(b-a)
  }

  lower <- c(1e-6, 0.5, rep(-Inf, 2*length(alphas)),0)
  upper <- c(0.01, 0.99, rep(Inf, 2*length(alphas)),0.01)
  
  res <- nloptr(x0 = par0,
                eval_f = obj_fun,
                lb = lower,
                ub = upper,
                eval_g_ineq = ineq_constraint,
                opts = list(algorithm = "NLOPT_LN_COBYLA",  
                            xtol_rel = 1e-6,
                            maxeval = 1000))
  
  out <- list(pars = res$solution,
              value = res$objective,
              convergence = res$status,
              message = res$message)
  return(out)
}

