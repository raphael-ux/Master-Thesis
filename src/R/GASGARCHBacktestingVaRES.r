library("GAS")
library(forecast)
data("sp500ret")
library(tidyr)
library(dplyr)
library(knitr)
source("/Users/raphaelbreaud/Downloads/GAS/R/hybridmodelmultiplealpha.r")
source("/Users/raphaelbreaud/Downloads/GAS/R/OneFactorGAS.r")
source("/Users/raphaelbreaud/Downloads/GAS/R/GARCHFZprof.r")
source("/Users/raphaelbreaud/Downloads/GAS/R/E_backtesting.r")
source("/Users/raphaelbreaud/Downloads/GAS/R/hybridmodel.r")
source("/Users/raphaelbreaud/Downloads/GAS/R/hybridmodel2.r")
library(rugarch)
library(data.table)
library(parallel)
library(esback)
library(pracma)


#DATA################################################################################################
crsp_data <- fread("/Users/raphaelbreaud/Downloads/GAS/data/CRSP_1926_2024_daily.csv")
crsp_data[, PRC := abs(PRC)] #putting the price positive
crsp_data[, log_ret := log(1 + as.numeric(RET)), by = PERMNO]
crsp_data[, market_cap := PRC * SHROUT * 1000] #compute market_cap (adding a new col)
latest_market_cap <- crsp_data[, .SD[which.max(date)], by = PERMNO] #get the last dates for each stock 
latest_market_cap <- latest_market_cap[, .(PERMNO, market_cap)] # get the last market_cap
top100_permnos <- latest_market_cap[order(-market_cap)][1:100, PERMNO] #get the last 100 to stocks (just PERMNO)
crsp_top500_data <- crsp_data[PERMNO %in% top500_permnos] #get the 100 hundread last stock 
######################################################################################################
crsp_top100_data[, log_ret_numeric := as.numeric(log_ret)] # the return can sometime be not numeric, have to hanlde it 
crsp_top100_data_sub <- crsp_top100_data[, .(PERMNO, date, log_ret_numeric)] #keep PERMNO, date and RET_numeric 
ret_wide_inv <- dcast(crsp_top500_data_sub, date ~ PERMNO, value.var = "log_ret_numeric") # col asset, row the price 
######################################################################################################
data <- ret_wide_inv[,1:100]
vAsset_not_clean <- colnames(data)[-1]
vAlpha <- c(0.01, 0.05) 
########### MODEL SPECIFICATION ######################################################################

#GAS-std
uGASSpec_std <- UniGASSpec(Dist = "std", ScalingType = "Identity",
                            GASPar = list(location = FALSE, scale = TRUE, shape = TRUE))

#GARCH-std
uGARCHpec_std  <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)),
  distribution.model = "std"
)


########### Forecast parameters #####################################################################

RefitEvery  <- 22 # refit every 22 days of traiding 
RefitWindow <- "moving" # keep a lenght of 1500 data

########### Univariate Rolling ######################################################################
luGASRoll_std  <- list() #GAS-std
luGASRoll_GG <- list() #Hybrid GAS-GARCH
#garch###############################################################################################
luGASRoll_GH <- list() #GARCH-std
Mu_garch_std <- list() #Mu GARCH-std
Sigma_garch_std <- list() #Sigma GARCH-std
Shape_garch_std <- list() #Shape GARCH-std
########################cluster#######################################################################
cluster <- makeCluster(4)
#######################

clean_data <- list()
data_clean <- list()
vAsset <- list()

iH <- list()
for(sAsset in vAsset_not_clean){
  clean_data[[sAsset]] <- data[[sAsset]]
  clean_data[[sAsset]] <- clean_data[[sAsset]][!is.na(clean_data[[sAsset]])]
  if(length(clean_data[[sAsset]]) >= 2500 & sAsset!="14542"){
    vAsset <- c(vAsset, sAsset)
    data_clean[[sAsset]] <- clean_data[[sAsset]][(length(clean_data[[sAsset]]) - 2500 + 1):length(clean_data[[sAsset]])]
    iH[[sAsset]] <- round(0.4*2500)
  }
}

vAsset <- vAsset[1:50]
i = 0
for(sAsset in vAsset){
  i <- i + 1 
}
iN <- length(vAsset)

for (i in 1:iN) {
  cat(i, "\n") 
  data_cl <- data_clean[[vAsset[[i]]]]
  H <- iH[[vAsset[[i]]]]
  print(mean(sapply(data_cl, function(x) sum(!is.finite(x)))))

  #GAS-std and GAS-norm##########################################################################################
  luGASRoll_std[[i]]  <- UniGASRoll(data = data_cl, GASSpec = uGASSpec_std, ForecastLength = H,
                                   RefitEvery = RefitEvery, RefitWindow = RefitWindow , cluster = cluster )
}
  ###############################################################################################################
stopCluster(cluster)
#GARCH###########################################################################################################
garch_cluster <- makeCluster(2) 
clusterEvalQ(garch_cluster, { library(rugarch) })
clusterExport(garch_cluster, c("uGARCHpec_std" , "uGARCHpec_norm" , "data_clean", "vAsset", "iH", "RefitEvery", "RefitWindow", "vAlpha"))

  #GARCH-std

  luGASRoll_GH <- parLapply(
    garch_cluster,
    vAsset,
    function(sAsset) {
      roll <- ugarchroll(
        spec          = uGARCHpec_std,
        data          = data_clean[[sAsset]],
        n.ahead       = 1,
        n.start       = length(data_clean[[sAsset]]) - iH[[sAsset]],
        refit.every   = RefitEvery,
        refit.window  = RefitWindow,
        solver        = "hybrid",
        calculate.VaR = TRUE,
        VaR.alpha     = vAlpha
      )
      list(
        roll  = roll,
        Mu    = roll@forecast$density[, "Mu"],
        Sigma = roll@forecast$density[, "Sigma"],
        Shape = roll@forecast$density[, "Shape"]
      )
    }
  )

stopCluster(garch_cluster)

###############################################################################

Mu_garch_std    <- lapply(luGASRoll_GH, `[[`, "Mu")
Sigma_garch_std <- lapply(luGASRoll_GH, `[[`, "Sigma")
Shape_garch_std <- lapply(luGASRoll_GH, `[[`, "Shape")
luGASRoll_GH <- lapply(luGASRoll_GH, `[[`, "roll")

Shape_gas_std <- lapply(luGASRoll_std, function(x) {x@Forecast$PointForecast[, 3]})
Mu_gas_std <- lapply(luGASRoll_std, function(x) {x@Forecast$PointForecast[, 1]})

#GAS-GARCH-std #################################################################

luGASRoll_GG <- luGASRoll_GH

for(i in 1:length(luGASRoll_GH)){
  luGASRoll_GG[[i]]@forecast$density$Shape <- as.numeric(Shape_gas_std[[i]])
}

###############################################################################



names(luGASRoll_std)  <- vAsset 
names(luGASRoll_GG) <- vAsset
names(luGASRoll_GH) <- vAsset
names(Mu_garch_std)<- vAsset
names(Sigma_garch_std)<- vAsset
names(Shape_garch_std)<- vAsset



################## EMPIRICAL ##################


lVaR <- list()
lES  <- list()

vDist <- c("GG")

iD <- length(vDist)

#obtain the forecasted Value-at-Risk for each distribution in vDist 

lVaR <- lapply( vAsset, function(sAsset) {

  iH_asset <- iH[[sAsset]]
  mOut <- array(NA, dim = c(iH_asset, length(vDist), length(vAlpha)),
                dimnames = list(1:iH_asset, vDist, vAlpha))

  
  for (sDist in vDist) {
    obj_GAS <- get(paste("luGASRoll_", sDist, sep = "")) 
    
    for (iAlpha in seq_along(vAlpha)) {  
      dAlpha <- vAlpha[iAlpha]
      mOut[, sDist, paste(dAlpha)] <- quantile(luGASRoll_GG[[sAsset]],dAlpha)
    }
  }
  return(mOut)
})

names(lVaR) <- vAsset

#obtain the forecasted Expected-Shortfall for each distribution in vDist 

lES <- lapply( vAsset, function(sAsset) {

  iH_asset <- iH[[sAsset]]
  
  mOut <- array(NA, dim = c(iH_asset, length(vDist), length(vAlpha)),
                dimnames = list(1:iH_asset, vDist, vAlpha))

  for (sDist in vDist) {
    obj_GAS <- get(paste("luGASRoll_", sDist, sep = ""))
    
    for (iAlpha in seq_along(vAlpha)) { 
      dAlpha <- vAlpha[iAlpha]
      alpha <- dAlpha
      u_grid <- seq(1e-7, alpha, length.out = 1000)

      q_mat <- sapply(u_grid, function(u) as.numeric(quantile(luGASRoll_GG[[sAsset]], u)))
      q_mat <- t(q_mat)  # now rows 

      mOut[, sDist, paste(dAlpha)] <- apply(q_mat, 2, function(q_vals) {trapz(u_grid, q_vals) / alpha})
    }
  }
  return(mOut)
})

names(lES) <- vAsset #3 dimension lenghtforecast-models-alpha


################## Var Backtesting ##################

lVaRBacktest <- list()

for (sAsset in vAsset) {
  lVaRBacktest[[sAsset]] <- list()
  data_cl <- data_clean[[sAsset]]
  H <- iH[[sAsset]]
  vRealised <- tail(data_cl,H)

  for (sDist in vDist) {
    lVaRBacktest[[sAsset]][[sDist]] = list()
    for (dAlpha in vAlpha) {
      if(anyNA(lVaR[[sAsset]][, sDist, paste(dAlpha)])){
        print("NA in ")
        print(sDist)
        print(paste(dAlpha))
        print(sAsset)
      }
      if(anyNA(vRealised)){
        print("NA is data")
      }
      lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]] <- BacktestVaR(data = vRealised, VaR = lVaR[[sAsset]][, sDist, paste(dAlpha)], alpha = dAlpha)
    }
  }
}

###################ES backtesting ###############################

lESCCBacktest <- list()

#Conditional Calibration test####################################

for (sAsset in vAsset) {
  lESCCBacktest[[sAsset]] <- list()
  data_cl <- data_clean[[sAsset]]
  H <- iH[[sAsset]]
  vRealised <- tail(data_cl,H)


  for (sDist in vDist) {
    lESCCBacktest[[sAsset]][[sDist]] = list()
    for (dAlpha in vAlpha) {
      if(anyNA(vRealised)){
        print("NA in returns")
      }
      if(anyNA(unlist(lVaR[[sAsset]][, sDist, paste(dAlpha)]))){
        print("NA in VaR for asset")
        print(paste(sAsset))
      }
      if(anyNA(unlist(lES[[sAsset]][, sDist, paste(dAlpha)]))){
        print("NA in ES for asset")
        print(paste(sAsset))
      }
      lESCCBacktest[[sAsset]][[sDist]][[paste(dAlpha)]] <- cc_backtest(r = vRealised, q = lVaR[[sAsset]][, sDist, paste(dAlpha)], e = lES[[sAsset]][, sDist, paste(dAlpha)] , s = NULL, alpha = dAlpha, hommel = TRUE)
    }
  }
}

#Expected Shortfall Regression test#####################################################

lESERBacktest <- list()

for (sAsset in vAsset) {
  lESERBacktest[[sAsset]] <- list()
  data_cl <- data_clean[[sAsset]]
  H <- iH[[sAsset]]
  vRealised <- tail(data_cl,H)

  for (sDist in vDist) {
    lESERBacktest[[sAsset]][[sDist]] = list()
    for (dAlpha in vAlpha) {
      lESERBacktest[[sAsset]][[sDist]][[paste(dAlpha)]] <- esr_backtest(r = vRealised ,q = lVaR[[sAsset]][, sDist, paste(dAlpha)],e = lES[[sAsset]][, sDist, paste(dAlpha)] ,alpha = dAlpha, version =2 , B=0 )
    }
  }
}

#e-backtesting test #####################################################################

lESETBacktest <- list()

for (sAsset in vAsset) {
  lESETBacktest[[sAsset]] <- list()
  data_cl <- data_clean[[sAsset]]
  H <- iH[[sAsset]]
  vRealised <- tail(data_cl,H)

  for (sDist in vDist) {
    lESETBacktest[[sAsset]][[sDist]] = list()
    for (dAlpha in vAlpha) {
      print(1-dAlpha)
      lESETBacktest[[sAsset]][[sDist]][[paste(dAlpha)]] <- e_backtest(L = -vRealised, r = -lES[[sAsset]][, sDist, paste(dAlpha)], z = -lVaR[[sAsset]][, sDist, paste(dAlpha)], e = e_ES, lambda = "GREM", threshold = 0.05, p = (1-dAlpha))
    }
  }
}

################## Var and ES joint Backtesting ###########################################

lVaRESBacktest <- list()

for (sAsset in vAsset) {
  lVaRESBacktest[[sAsset]] <- list()
  data_cl <- data_clean[[sAsset]]
  H <- iH[[sAsset]]
  vRealised <- tail(data_cl,H)

  for (sDist in vDist) {
    lVaRESBacktest[[sAsset]][[sDist]] = list()
    for (dAlpha in vAlpha) {

      lVaRESBacktest[[sAsset]][[sDist]][[paste(dAlpha)]] <- FZLoss(data = vRealised,
                                                                   VaR = lVaR[[sAsset]][, sDist, paste(dAlpha)],
                                                                   ES  = lES[[sAsset]][, sDist, paste(dAlpha)],
                                                                   alpha = dAlpha) 
    }
  }
}

#lVaRESBacktest 3 dimensions : assets-models-alpha like lVaRESBacktest[["AA"]][["std"]][["0.01"]] 
# can do mean(lVaRESBacktest[["AA"]][["std"]][["0.01"]]) and obtain the mean of the FZ loss.
#print(mean(lVaRESBacktest[["AA"]][["std"]][["0.01"]], na.rm = TRUE)) #not taking the NAN to acount, if ES >=0 then will obtain a NAN 



################ TABLES ########################

vModels <- paste("GAS", vDist, sep = "-")
vTestNames <- c("LRuc", "LRcc", "DQ", "Loss-QL", "Loss-FZL", "AE", "ADmean", "ADmax" )
vESCCTestNames <- c("twosided_simple" , "onesided_simple")
vInfoCoef <- c("mean","variance")
cBackTest <- array(NA, dim = c(length(vAsset), length(vTestNames),length(vModels), length(vAlpha)),
                   dimnames = list(vAsset, vTestNames, vModels, vAlpha)) #VaR backtest 
#Conditional Calibration test
cBackTestCC <- array(NA , dim = c(length(vAsset) , length(vESCCTestNames) ,length(vModels) , length(vAlpha)) , dimnames = list(vAsset , vESCCTestNames  , vModels , vAlpha))

#Expected Shortfall Regression test 
cBackTestESER <- array(NA , dim = c(length(vAsset)  , length(vModels) , length(vAlpha)) , dimnames = list(vAsset  , vModels , vAlpha))

#e-backtesting 
cBackTestESET <- array(NA , dim = c(length(vAsset)  , length(vModels) , length(vAlpha)) , dimnames = list(vAsset  , vModels , vAlpha))

cRoll <- data.frame()
for (sAsset in vAsset) {
  for (dAlpha in vAlpha) {
    for (sDist in vDist) {
      #p-value LRuc test 
      cBackTest[sAsset, "LRuc", paste("GAS", sDist, sep = "-"), paste(dAlpha)]     <- lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$LRuc["Pvalue"]

      #p-values LRcc test 
      cBackTest[sAsset, "LRcc", paste("GAS", sDist, sep = "-"), paste(dAlpha)]     <- lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$LRcc["Pvalue"]

      #p-value DQ test 
      cBackTest[sAsset, "DQ", paste("GAS", sDist, sep = "-"), paste(dAlpha)]       <- as.numeric(lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$DQ$pvalue)

      #AE values
      cBackTest[sAsset, "AE", paste("GAS", sDist, sep = "-"), paste(dAlpha)]       <- lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$AE

      #other metrics
      cBackTest[sAsset, "ADmean", paste("GAS", sDist, sep = "-"), paste(dAlpha)]   <- lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$AD["ADmean"]
      cBackTest[sAsset, "ADmax", paste("GAS", sDist, sep = "-"), paste(dAlpha)]    <- lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$AD["ADmax"]
      cBackTest[sAsset, "Loss-QL", paste("GAS", sDist, sep = "-"), paste(dAlpha)]  <- lVaRBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$Loss$Loss
      cBackTest[sAsset, "Loss-FZL", paste("GAS", sDist, sep = "-"), paste(dAlpha)] <- mean(lVaRESBacktest[[sAsset]][[sDist]][[paste(dAlpha)]], na.rm = TRUE)

      #two sided p-value Conditional Calibration test
      cBackTestCC[sAsset , "twosided_simple" , paste("GAS", sDist, sep = "-"), paste(dAlpha)] <- lESCCBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$pvalue_twosided_simple

      #one sided p-vaue Conditional Calibration test
      cBackTestCC[sAsset , "onesided_simple" , paste("GAS", sDist, sep = "-"), paste(dAlpha)] <- lESCCBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$pvalue_onesided_simple

      #p-value Expected Shortfall Regression test
      cBackTestESER[sAsset, paste("GAS", sDist, sep = "-"), paste(dAlpha)] <- lESERBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]$pvalue_twosided_asymptotic

      #e-backtesting metric 
      cBackTestESET[sAsset, paste("GAS", sDist, sep = "-"), paste(dAlpha)] <- lESETBacktest[[sAsset]][[sDist]][[paste(dAlpha)]]
    }
  }
}


#tables #######################################################################################################

sDir_Tables <- "/Users/raphaelbreaud/Downloads/GAS/tables"  
alpha1 <- paste(vAlpha[1])
alpha2 <- paste(vAlpha[2])

#############GG###################################
if("GG" %in% vDist){
  mTab_Var_GG <- format(round(cbind(cBackTest[, "LRuc", "GAS-GG", alpha1], cBackTest[, "LRcc", "GAS-GG", alpha1] , cBackTest[, "DQ", "GAS-GG", alpha1], cBackTest[, "AE", "GAS-GG", alpha1]),2), digits = 4, scientific = FALSE) #
  mTab_Var_GG_numeric <- cbind(cBackTest[, "LRuc", "GAS-GG", alpha1], cBackTest[, "LRcc", "GAS-GG", alpha1] , cBackTest[, "DQ", "GAS-GG", alpha1], cBackTest[, "AE", "GAS-GG", alpha1])
  mTab_ES_GG <- format(round(cbind(cBackTestCC[, "twosided_simple", "GAS-GG", alpha1] ,   cBackTestESER[ , "GAS-GG" , alpha1]    ,   cBackTestESET[ , "GAS-GG" , alpha1]) , 2), digits = 4 , scientific = FALSE)
  mTab_ES_GG_numeric <- cbind(cBackTestCC[, "twosided_simple", "GAS-GG", alpha1] ,   cBackTestESER[ , "GAS-GG" , alpha1]    ,   cBackTestESET[ , "GAS-GG" , alpha1])


  for (i in 1:3) {
  mTab_Var_GG[mTab_Var_GG_numeric[, i] >= 0.05, i] <-
    paste0("\\textbf{", mTab_Var_GG[mTab_Var_GG_numeric[, i] >= 0.05, i], "}")
  }

   for (i in 1:2) {
    mTab_ES_GG[mTab_ES_GG_numeric[, i] >= 0.05, i] <-
    paste0("\\textbf{", mTab_ES_GG[mTab_ES_GG_numeric[, i] >= 0.05, i], "}")
  }
  
  write.table(rbind(col_labels, cbind(rownames(mTab_Var_GG), mTab_Var_GG)),
            file = paste(sDir_Tables, "/mTab_Var_GG.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

  write.table(rbind(col_labels, cbind(rownames(mTab_ES_GG), mTab_ES_GG)),
            file = paste(sDir_Tables, "/mTab_ES_GG.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

}