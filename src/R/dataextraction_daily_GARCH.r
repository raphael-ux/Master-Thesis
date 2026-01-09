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
  return(head(crsp[.(current_date), nomatch = 0][eligible_dt, on = "PERMNO"]$PERMNO,100))
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

extraction <- function(sorting,initial_date,dates,H,F,type = c("volatility", "ES" , "VaR") ,cluster = NULL){
    n_quantiles <- 5
    portfolio_returns <- list()
    available_assets <- list()
    lData <- list()
    ES <- list()
    vola <- list()
    iT <- length(dates)  

    initial_date_idx <- match(initial_date, dates)
    if (length(initial_date_idx) == 0 || initial_date_idx < F) {
        stop("incorrect initial date")
    }
    cDates <- dates[seq(from = initial_date_idx, to = iT, by = H)]  

    first_date <- cDates[1]
    last_date <- cDates[length(cDates)]
    idx_first <- match(first_date,dates)
    idx_last <- match(last_date , dates)
    cDates <- dates[idx_first:idx_last]

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
                    asset_data <- get_data_asset(asset = as.character(asset), current_date = dates[(idx+1):(idx+1)])
                    if(anyNA(asset_data)){
                      return(-1)
                    }else{
                      return(prod(asset_data + 1) - 1)
                      
                    }                   
                  })
                  portfolio_returns[[paste(current_date)]][[j]] <- mean(re)
              }
          }
      }
  return(list(portfolio = portfolio_returns , dates = cDates))
   
}



dates <- all_dates
initial_date <- dates[[10570]]
H <- 22

ES <- readRDS("../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/ES_daily_GARCH.rds")
VaR <- readRDS("../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/VaR_daily_GARCH.rds")
vola <- readRDS("../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/vola_daily_GARCH.rds")

test_vola <- extraction(sorting = vola,initial_date = initial_date ,dates = dates,H = 22,F = 1000,type = "volatility" ,cluster = NULL)
portfolio_vola <- test_vola$portfolio

test_VaR <- extraction(sorting = VaR,initial_date = initial_date ,dates = dates,H = 22,F = 1000,type = "VaR" ,cluster = NULL)
portfolio_VaR <- test_VaR$portfolio

test_ES <- extraction(sorting = ES,initial_date = initial_date ,dates = dates,H = 22,F = 1000,type = "ES" ,cluster = NULL)
portfolio_ES <- test_ES$portfolio

cDates <- as.Date(test_ES$dates)[-length(test_ES$dates)]

return_ES <- list()


for(i in 1:5){
  return_ES[[i]] <- sapply(portfolio_ES, function(x) x[[i]])
}

Cumulative_returns_ES <- list()
for (i in 1:length(return_ES)) {
  Cumulative_returns_ES[[i]] <- cumsum(log(return_ES[[i]] + 1))
}

volatility_ES <- list()
for(i in 1:length(return_ES)){
  volatility_ES[[i]] <- sd(return_ES[[i]])
}

return_VaR <- list()


for(i in 1:5){
  return_VaR[[i]] <- sapply(portfolio_VaR, function(x) x[[i]])
}

Cumulative_returns_VaR <- list()
for (i in 1:length(return_VaR)) {
  Cumulative_returns_VaR[[i]] <- cumsum(log(return_VaR[[i]] + 1))
}

volatility_VaR <- list()
for(i in 1:length(return_ES)){
  volatility_VaR[[i]] <- sd(return_VaR[[i]])
}

return_vola <- list()

for(i in 1:5){
  return_vola[[i]] <- sapply(portfolio_vola, function(x) x[[i]])
}

Cumulative_returns_vola <- list()
for (i in 1:length(return_vola)) {
  Cumulative_returns_vola[[i]] <- cumsum(log(return_vola[[i]] + 1))
}

volatility_vola <- list()
for(i in 1:length(return_vola)){
  volatility_vola[[i]] <- sd(return_vola[[i]])
}

H <- 1

Monthly_returns <- list()
cols_FF5 <- c("Mkt_RF","SMB","HML","RMW","CMA","RF")
FF5_data[, (cols_FF5) := lapply(.SD, function(x) x/100), .SDcols = cols_FF5]
for(i in 1:length(cDates)){
  current_date <- cDates[i]
  idx <- match(current_date,dates)
  Monthly_returns[[paste(current_date)]] <- prod(FF5_data[date %in% dates[(idx + 1):(idx + H)]]$RF + 1)-1
}



Sharpe_ratio_ES <- list()

for(i in 1:length(return_ES)){
    Sharpe_ratio_ES[[i]] <- mean(return_ES[[i]] - unlist(Monthly_returns))/(volatility_ES[[i]])
}

Sharpe_ratio_VaR <- list()

for(i in 1:length(return_VaR)){
    Sharpe_ratio_VaR[[i]] <- mean(return_VaR[[i]] - unlist(Monthly_returns))/(volatility_VaR[[i]])
}

Sharpe_ratio_vola <- list()

for(i in 1:length(return_vola)){
    Sharpe_ratio_vola[[i]] <- mean(return_vola[[i]] - unlist(Monthly_returns))/(volatility_vola[[i]])
}

Sharpe_ratio_ES_df <- as.data.frame(Sharpe_ratio_ES)
Sharpe_ratio_VaR_df <- as.data.frame(Sharpe_ratio_VaR)
Sharpe_ratio_vola_df <- as.data.frame(Sharpe_ratio_vola)

write.csv(Sharpe_ratio_ES_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/Sharpe_ratio_ES_GAS.csv",row.names = FALSE)
write.csv(Sharpe_ratio_VaR_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/Sharpe_ratio_VaR_GAS.csv",row.names = FALSE)
write.csv(Sharpe_ratio_vola_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/Sharpe_ratio_vola_GAS.csv",row.names = FALSE)


returns_ES_df <- as.data.frame(return_ES)
returns_VaR_df <- as.data.frame(return_VaR)
returns_vola_df <- as.data.frame(return_vola)


write.csv(returns_ES_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/return_ES_GAS.csv",row.names = FALSE)
write.csv(returns_VaR_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/return_VaR_GAS.csv",row.names = FALSE)
write.csv(returns_vola_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/return_vola_GAS.csv",row.names = FALSE)


Cumulative_returns_ES_df <- as.data.frame(Cumulative_returns_ES)
Cumulative_returns_VaR_df <- as.data.frame(Cumulative_returns_VaR)
Cumulative_returns_vola_df <- as.data.frame(Cumulative_returns_vola)


write.csv(Cumulative_returns_ES_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/cumulatives_returns_ES_GAS.csv",row.names = FALSE)
write.csv(Cumulative_returns_VaR_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/cumulatives_returns_VaR_GAS.csv",row.names = FALSE)
write.csv(Cumulative_returns_vola_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/cumulatives_returns_vola_GAS.csv",row.names = FALSE)


volatility_ES_df <- as.data.frame(volatility_ES)
volatility_VaR_df <-  as.data.frame(volatility_VaR)
volatility_vola_df <- as.data.frame(volatility_vola)

write.csv(volatility_ES_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/volatility_ES_GAS.csv",row.names = FALSE)
write.csv(volatility_VaR_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/volatility_VaR_GAS.csv",row.names = FALSE)
write.csv(volatility_vola_df,"../results/Portfolios_Backtesting/GARCH_Sorting/Daily_Sorting/volatility_vola_GAS.csv",row.names = FALSE)


names(Cumulative_returns_vola) <- as.character(1:length(Cumulative_returns_vola))

all_values <- unlist(Cumulative_returns_vola)
ylim_range <- range(all_values, na.rm = TRUE) 

colors <- rainbow(length(Cumulative_returns_vola))  

plot(cDates, Cumulative_returns_vola[[1]], type = "l", col = 1,
     xlab = "Date", ylab = "Cumulative Return",
     main = "Cumulative Returns Over Time",
     ylim = ylim_range , xaxt = "n")

for (i in 1:length(Cumulative_returns_vola)) {
  lines(cDates, Cumulative_returns_vola[[i]], col = colors[i], lwd = 1)
}

axis.Date(1, at = seq(min(cDates), max(cDates), by = "1 year"), format = "%Y", las = 2, cex.axis = 0.8)
legend("topleft", legend = names(Cumulative_returns_vola), col = colors, lwd = 1)

