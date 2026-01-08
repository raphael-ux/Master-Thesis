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
    print(permno)
    print(permno_splits)
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
  past_ok   <- colMeans(presence[(idx-F + 1):(idx), ]) == 1
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

get_sorting <- function(type, filter = FALSE,filt = NULL , forecast = FALSE , forr = NULL){
  if(type == "volatility"){
    if(filter){
      last_col <- ncol(filt$mTheta)
      return(filt$mTheta[2, last_col])
    }else if(forecast){
      return(mean(as.numeric(forr@Forecast$PointForecast[,"scale"])))
    }
  }else if(type == "VaR"){
    if(filter){
      last_col <- ncol(filt$mTheta)
      mTheta <- filt$mTheta[, last_col , drop = FALSE]   
      mQuantile = Quantiles(mTheta, "std", 0.01)
      return(mQuantile)
    }else if(forecast){
      VaR <- quantile(forr,0.01)
      percent_VaR <- min(VaR)
      return(percent_VaR)
    }
  }else if(type == "ES"){
    if(filter){
      last_col <- ncol(filt$mTheta)
      mTheta <- filt$mTheta[, last_col , drop = FALSE]   
      mQuantile = Quantiles(mTheta, "std", 0.01)
      ES <- adaptIntegrate(Quantiles, lowerLimit = 1e-7, upperLimit = 0.01,mTheta = mTheta[, , drop = FALSE], Dist = "std")$integral
      ES <- t(t(ES) / 0.01)
      return(ES)
    }else if(forecast){
      ES <- ES(forr,0.01)
      percent_ES <- min(ES)
      return(percent_ES)
    }
  }
}



uGASSpec_std <- UniGASSpec(Dist = "std", ScalingType = "Identity",GASPar = list(location = FALSE, scale = TRUE, shape = TRUE))



rolling_portfolio <- function(initial_date,dates,H,F,type = c("volatility", "ES" , "VaR") ,cluster = NULL){

    n_quantiles <- 5
    iT <- length(dates) 
    portfolio_returns <- list()
    ES <- list()
    VaR <- list()
    vola <- list()
    lData <- list()
    available_assets <- list()
    initial_date_idx <- match(initial_date, dates)
    if (length(initial_date_idx) == 0 || initial_date_idx < F) {
        stop("incorrect initial date")
    }
    cDates <- dates[seq(from = initial_date_idx, to = iT, by = H)] 

    for(i in 1:(length(cDates)-1)){
        date <- cDates[i]
        date_idx <- match(date,dates)
        lData[[paste(date)]] <- dates[(date_idx - F + 1):(date_idx)]
    }

    first_date <- cDates[1]
    last_date <- cDates[length(cDates)]
    idx_first <- match(first_date,dates)
    idx_last <- match(last_date , dates)
    dates_ <- dates[idx_first:(idx_last-1)]

    for(i in 1:length(dates_)){
        date <- dates_[i]
        available_assets[[paste(date)]] <- available_full(current_date = date  , F = F)
    }  

    for(i in 1:length(dates_)){
        sorting <- list()
        
        current_date <- dates_[i]
        available_ <- available_assets[[paste(current_date)]]

        ES[[paste(current_date)]] <- list()
        VaR[[paste(current_date)]] <- list()
        vola[[paste(current_date)]] <- list()

        if(current_date %in% cDates){
            fit <- list()
            path <- list()
            available <- available_
            past_window <- lData[[paste(current_date)]]
            print(paste("i =", i, "/", (length(dates_)), " - current_date:", current_date))
            print(paste("available assets" ,length(available)))
            for(asset in available){
                data_asset <- get_data_asset(asset = as.character(asset), current_date = past_window)
                log_data_asset <- as.numeric(log(1 + data_asset))
                fit[[as.character(asset)]] <- UniGASFit(GASSpec = uGASSpec_std, data = log_data_asset)
                forecast <-  UniGASFor(uGASFit = fit[[as.character(asset)]], H = 1, Roll = FALSE, B = 10000L, Bands = c(0.1,0.15,0.85,0.9), ReturnDraws = TRUE)
                sorting[[as.character(asset)]] <- get_sorting(type = "volatility", filter = FALSE,filt = NULL , forecast = TRUE , forr = forecast)#getting the forecasted volatility 
                ES[[paste(current_date)]][[as.character(asset)]] <- get_sorting(type = "ES", filter = FALSE,filt = NULL , forecast = TRUE , forr = forecast)#getting the forecasted volatility 
                VaR[[paste(current_date)]][[as.character(asset)]] <- get_sorting(type = "VaR", filter = FALSE,filt = NULL , forecast = TRUE , forr = forecast)#getting the forecasted volatility 
                vola[[paste(current_date)]][[as.character(asset)]] <- get_sorting(type = "volatility", filter = FALSE,filt = NULL , forecast = TRUE , forr = forecast)#getting the forecasted volatility 
                path[[as.character(asset)]] <- log_data_asset
            }
        }else{
            available__ <- intersect(available_, available)
            print(paste("i =", i, "/", (length(dates_)), " - current_date:", current_date))
            print(paste("available assets" ,length(available__)))
            for(asset in available__){
                fit_asset <- fit[[as.character(asset)]]
                vKappa <- fit_asset@Estimates$lParList$vKappa 
                mA <- fit_asset@Estimates$lParList$mA
                mB <- fit_asset@Estimates$lParList$mB 

                data_asset <- log(get_data_asset(asset = as.character(asset) , current_date = current_date) + 1) #getting the return for this date 
                if(anyNA(data_asset)){
                    print("NA in data")
                    print(asset)
                }
                path[[as.character(asset)]] <- c(path[[as.character(asset)]],data_asset) #adding this return to my in_sample_data
                forecast <- GASFilter_univ(path[[as.character(asset)]], vKappa, mA, mB, iT = length(path[[as.character(asset)]]), iK = 3,Dist = "std", ScalingType = "Identity")#forecast next day 
                last_col <- ncol(forecast$mTheta)
                sorting[[as.character(asset)]] <- get_sorting(type = "volatility", filter = TRUE,filt = forecast , forecast = FALSE , forr = NULL)#getting the forecasted volatility 
                ES[[paste(current_date)]][[as.character(asset)]] <- get_sorting(type = "ES", filter = TRUE,filt = forecast , forecast = FALSE , forr = NULL)#getting the forecasted volatility 
                VaR[[paste(current_date)]][[as.character(asset)]] <- get_sorting(type = "VaR", filter = TRUE,filt = forecast , forecast = FALSE , forr = NULL)#getting the forecasted volatility 
                vola[[paste(current_date)]][[as.character(asset)]] <- get_sorting(type = "volatility", filter = TRUE,filt = forecast , forecast = FALSE , forr = NULL)#getting the forecasted volatility 
            }
        }
        
        sorting_vec <- unlist(sorting)
        sorting_vec <- sorting_vec[!is.na(sorting_vec)]

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
                print("Group null")
                portfolio_returns[[paste(current_date)]][[j]] <- NaN
            }else{
                idx <- match(current_date,dates)
                re <-  sapply(group, function(asset){
                  asset_data <- get_data_asset(asset = as.character(asset), current_date = dates[(idx+1)])
                  if(anyNA(asset_data)){
                    return(-1)
                  }else{
                    return(asset_data)
                  }                   
                })
                portfolio_returns[[paste(current_date)]][[j]] <- mean(re)
            }
        }
    }
    return(list(portfolio = portfolio_returns , ES = ES , VaR = VaR , vola = vola))
}



dates <- all_dates 
initial_date <- dates[[10570]]

H <- 22

test1 <- rolling_portfolio(initial_date = initial_date ,dates = dates,H = 22,F = 1000,type = "volatility" ,cluster = NULL)

test <- test1$portfolio

cDates <- dates[10570:length(dates)]


Monthly_returns <- list()
cols_FF5 <- c("Mkt_RF","SMB","HML","RMW","CMA","RF")
FF5_data[, (cols_FF5) := lapply(.SD, function(x) x/100), .SDcols = cols_FF5]
for(i in 1:length(cDates)){
  current_date <- cDates[i]
  Monthly_returns[[paste(current_date)]] <- FF5_data[date %in% current_date]$RF
}




return <- list()
for(i in 1:5){
  return[[i]] <- sapply(test, function(x) x[[i]])
}

Cumulative_returns <- list()
for (i in 1:length(return)) {
  Cumulative_returns[[i]] <- cumsum(log(return[[i]] + 1))
}


volatility <- list()
for (i in 1:length(return)) {
  volatility[[i]] <- sd(return[[i]])
}

Sharpe_ratio <- list()

for(i in 1:length(return)){
    Sharpe_ratio[[i]] <- mean(return[[i]] - unlist(Monthly_returns))/(volatility[[i]])
}


ES <- test1$ES
VaR <- test1$VaR
vola <- test1$vola

saveRDS(ES, "/Users/raphaelbreaud/Downloads/GAS/results/ES_daily.rds")
saveRDS(VaR, "/Users/raphaelbreaud/Downloads/GAS/results/VaR_daily.rds")
saveRDS(vola, "/Users/raphaelbreaud/Downloads/GAS/results/vola_daily.rds")
