library("GAS")
library(data.table)
library(rugarch)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(dplyr)
library(cubature)
library(future.apply)
library(memoise)
library(pracma)
sourceCpp("/Users/raphaelbreaud/Downloads/GAS/src/Filters.cpp")

FF5_data <- fread(
  "/Users/raphaelbreaud/Downloads/GAS/data/F-F_Research_Data_5_Factors_2x3_daily.csv",
  skip = 3,       
  sep = ",",      
  header = FALSE  
)
setnames(FF5_data, c("date", "Mkt_RF", "SMB", "HML", "RMW", "CMA", "RF"))
FF5_data[, c("Mkt_RF","SMB","HML","RMW","CMA","RF") := lapply(.SD, as.numeric), .SDcols = c("Mkt_RF","SMB","HML","RMW","CMA","RF")]
FF5_data <- FF5_data[!is.na(date)]
FF5_data[, date := as.Date(as.character(date), format = "%Y%m%d")]

crsp_data <- fread("/Users/raphaelbreaud/Downloads/GAS/data/CRSP_1926_2024_daily.csv")
crsp_data[, PRC := abs(PRC)]
crsp_data[, RET_numeric := as.numeric(RET)]
crsp_data <- crsp_data[PRC > 5]
crsp_data[, market_cap := PRC * SHROUT * 1000]
crsp <- crsp_data[, .(PERMNO, date, market_cap)]
setorder(crsp, date, -market_cap)
setkey(crsp, date, PERMNO)
all_dates <- sort(unique(crsp_data[!is.na(RET_numeric)]$date))
all_dates_dt <- data.table(date = all_dates)
permnos <- unique(crsp_data$PERMNO)
n_splits <- 5
permno_splits <- split(permnos, cut(seq_along(permnos), n_splits, labels = FALSE))

wide_tables <- vector("list", length = n_splits)
for (i in 1:n_splits) {
  dt_subset <- crsp_data[PERMNO %in% permno_splits[[i]]]
  wide_dt <- dcast(dt_subset, date ~ PERMNO, value.var = "RET_numeric")
  wide_dt <- merge(all_dates_dt, wide_dt, by = "date", all.x = TRUE)
  wide_tables[[i]] <- wide_dt
}

get_asset <- function(permno, permno_splits) {
  idx <- which(sapply(permno_splits, function(x) permno %in% x))
  if(length(idx) == 0){
    print(idx)
  }
  return(idx)
}

presence_maps <- list()

for (i in 1:length(wide_tables)) {
  presence_maps[[i]] <- !is.na(wide_tables[[i]][, -1]) 
}


available_in_split <- function(split_index, idx,F) {
  presence <- presence_maps[[split_index]] 
  if (idx < F){ return(rep(FALSE, ncol(presence)))}
  past_ok   <- colMeans(presence[(idx-F + 1):(idx),]) == 1
  return(past_ok)
}


available_full <- function(current_date, F) {
  eligible_permnos <- c()

  for (i in 1:length(wide_tables)) {
    dates_split <- wide_tables[[i]]$date
    idx <- which(current_date == dates_split)
    if (length(idx) != 0) {
      available_flags <- available_in_split(i, idx, F) 
      permnos_split <- colnames(wide_tables[[i]])[-1] 
      eligible_permnos <- c(eligible_permnos, permnos_split[available_flags])
    }
  }
  eligible_permnos <- as.integer(eligible_permnos)
  eligible_dt <- data.table(PERMNO = eligible_permnos)
  setkey(eligible_dt, PERMNO)
  return(head(crsp[.(current_date), nomatch = 0][eligible_dt, on = "PERMNO"]$PERMNO,500))
}


uniform_portfolio <- function(group){
  n <- length(group)
  return(rep(1/n,n))
}

get_data_asset <- function(asset = asset, current_date = current_date){
  split_idx <- get_asset(permno = asset, permno_splits = permno_splits)
  dt_wide <- wide_tables[[split_idx]]
  data_asset <- as.numeric(dt_wide[date %in% current_date, ..asset][[1]])
  return(as.numeric(data_asset))
}

uGARCHpec_std  <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)),
  distribution.model = "std"
)

rolling_portfolio <- function(initial_date,dates,H,F,type = c("volatility", "ES" , "VaR") ,cluster = NULL){
    n_quantiles <- 5
    portfolio_returns <- list()
    available_assets <- list()
    lData <- list()
    sorting <- list()
    ES <- list()
    VaR <- list()
    vola <- list()
    iT <- length(dates)  

    initial_date_idx <- match(initial_date, dates)
    if (length(initial_date_idx) == 0 || initial_date_idx < F) {
        stop("incorrect initial date")
    }
    cDates <- dates[seq(from = initial_date_idx, to = iT, by = H)] 
    for(i in 1:(length(cDates)-1)){
        date <- cDates[i]
        print(date)
        date_idx <- match(date,dates)
        lData[[paste(date)]] <- dates[(date_idx - F + 1):(date_idx)]
        available_assets[[paste(date)]] <- available_full(current_date = date  , F = F)
    }
    print("fitting part")

    for(i in 1:(length(cDates)-1)){
        current_date <- cDates[i]
        past_window <- lData[[paste(current_date)]]
        available <- available_assets[[paste(current_date)]]

        print(paste("i =", i, "/", (length(cDates)-1), " - current_date:", current_date))
        print(paste("available assets" ,length(available)))

        sorting[[paste(current_date)]] <- list()
        ES[[paste(current_date)]] <- list()
        VaR[[paste(current_date)]] <- list()
        vola[[paste(current_date)]] <- list()

        for(asset in available){

            data_asset <- get_data_asset(asset = as.character(asset), current_date = past_window)
            log_data_asset <- as.numeric(log(1 + data_asset))
            if(anyNA(data_asset)){
              print(paste("NA in the asset" , asset , "in the date " , current_date))
            }
            fit_try <- tryCatch(
                  ugarchfit(spec = uGARCHspec_std, data = log_data_asset),
                  error = function(e) NULL
                 )
            
            if (is.null(fit_try) || fit_try@fit$convergence != 0) {
                  available <- setdiff(available, asset)
                  cat("GARCH failed – skipping:", asset, "\n")
                next
                }

            forecast <- ugarchforecast(fitORspec = fit_try, data = NULL, n.ahead = H)
            
            sorting[[paste(current_date)]][[as.character(asset)]] <- mean(as.numeric(sigma(forecast)))
            vola[[paste(current_date)]][[as.character(asset)]] <- mean(as.numeric(sigma(forecast)))
            VaR[[paste(current_date)]][[as.character(asset)]] <- min(quantile(forecast, probs = 0.01))

            alpha  <- 0.01
            u_grid <- seq(1e-7, alpha, length.out = 1000)
            q_mat <- sapply(u_grid, function(u)as.numeric(quantile(forecast, probs = u)))
            q_mat <- t(q_mat)   
            ES_h <- apply(q_mat, 2, function(q_vals) {trapz(u_grid, q_vals) / alpha})
            ES[[paste(current_date)]][[as.character(asset)]] <- min(ES_h)

        }
    }
    
    

    for(i in 1:(length(cDates)-1)){
        current_date <- cDates[i]
        sorting_vec <- unlist(sorting[[paste(current_date)]])
        sorting_vec <- sorting_vec[!is.na(sorting_vec)]
        print(paste("i =", i, "/", (length(cDates)-1), " - current_date:", current_date))
        if(type == "volatility"){
          qtiles <- ntile(sorting_vec, n_quantiles) 
        }else{
          qtiles <- ntile(-sorting_vec, n_quantiles)
        }
        quintile_groups <- lapply(1:n_quantiles, function(q) {names(sorting_vec[qtiles == q])})
        if(is.null(portfolio_returns[[paste(current_date)]])){
            portfolio_returns[[paste(current_date)]] <- vector("list", n_quantiles)
        }

        for(j in 1:n_quantiles){
            group <- quintile_groups[[j]]
            if(length(group) == 0){
                portfolio_returns[[paste(current_date)]][[j]] <- NaN
            }else{
                idx <- match(current_date,dates)
                re <-  sapply(group, function(asset){
                  asset_data <- get_data_asset(asset = as.character(asset), current_date = dates[(idx+1):(idx+H)])
                  if(anyNA(asset_data)){
                    print("OUCH")
                    return(-1)
                  }else{
                    return(prod(asset_data + 1) - 1)
                  }                   
                })
                portfolio_returns[[paste(current_date)]][[j]] <- mean(re)
            }
        }
    }
  return(list(portfolio = portfolio_returns , dates = cDates , ES = ES , VaR = VaR , vola = vola))
}



####################################################################################



dates <- all_dates 
initial_date <- dates[[10570]]

H <- 22

test <- rolling_portfolio(initial_date = initial_date ,dates = dates,H = 22,F = 1000,type = "volatility" ,cluster = NULL)

portfolio <- test$portfolio
cDates <- as.Date(test$dates)[-length(test$dates)]


Monthly_returns <- list()
cols_FF5 <- c("Mkt_RF","SMB","HML","RMW","CMA","RF")
FF5_data[, (cols_FF5) := lapply(.SD, function(x) x/100), .SDcols = cols_FF5]
for(i in 1:length(cDates)){
  current_date <- cDates[i]
  idx <- match(current_date,dates)
  Monthly_returns[[paste(current_date)]] <- prod(FF5_data[date %in% dates[(idx + 1):(idx + H)]]$RF + 1)-1
}

return <- list()

for(i in 1:5){
  return[[i]] <- sapply(portfolio, function(x) x[[i]])
}

Cumulative_returns <- list()
for (i in 1:length(return)) {
  Cumulative_returns[[i]] <- cumsum(log(return[[i]] + 1))
}

volatility <- list()
for(i in 1:length(return)){
  volatility[[i]] <- sd(return[[i]])
}

sharpe_ratio <- list()
for(i in 1 : length(return)){
  sharpe_ratio[[i]] <- mean(return[[i]] - unlist(Monthly_returns))/(volatility[[i]])
}

ES <- test$ES
VaR <- test$VaR
vola <- test$vola

saveRDS(ES, "/Users/raphaelbreaud/Downloads/GAS/results/ES_GARCH.rds")
saveRDS(VaR, "/Users/raphaelbreaud/Downloads/GAS/results/VaR_GARCH.rds")
saveRDS(vola, "/Users/raphaelbreaud/Downloads/GAS/results/vola_GARCH.rds")