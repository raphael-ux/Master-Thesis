library("GAS")
library(forecast)
data("sp500ret")
library(tidyr)
library(dplyr)
library(knitr)
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
crsp_top100_data <- crsp_data[PERMNO %in% top100_permnos] #get the 100 hundread last stock 
######################################################################################################
crsp_top100_data[, log_ret_numeric := as.numeric(log_ret)] # the return can sometime be not numeric, have to hanlde it 
crsp_top100_data_sub <- crsp_top100_data[, .(PERMNO, date, log_ret_numeric)] #keep PERMNO, date and RET_numeric 
ret_wide_inv <- dcast(crsp_top100_data_sub, date ~ PERMNO, value.var = "log_ret_numeric") # col asset, row the price 
######################################################################################################
data <- ret_wide_inv[,1:100]
vAsset_not_clean <- colnames(data)[-1]
vAlpha <- c(0.01, 0.05) 
########### MODEL SPECIFICATION ######################################################################

#GARCH-std
uGARCHpec_std  <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)),
  distribution.model = "std"
)

#GARCH-norm
uGARCHpec_norm  <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)),
  distribution.model = "norm"
)

########### Forecast parameters #####################################################################

RefitEvery  <- 22 # refit every 22 days of traiding 
RefitWindow <- "moving" # keep a lenght of 1500 data

#garch###############################################################################################
luGASRoll_GH <- list() #GARCH-std
Mu_garch_std <- list() #Mu GARCH-std
Sigma_garch_std <- list() #Sigma GARCH-std
Shape_garch_std <- list() #Shape GARCH-std

luGASRoll_GHN <- list() #GARCH-norm
Mu_garch_norm <- list() #Mu GARCH-norm
Sigma_garch_norm <- list() #Sigma GARCH-norm
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

  #GARCH-norm

  luGASRoll_GHN <- parLapply(
    garch_cluster,
    vAsset,
    function(sAsset) {
      roll <- ugarchroll(
        spec          = uGARCHpec_norm,
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
        Sigma = roll@forecast$density[, "Sigma"]
     )
    }
  )

stopCluster(garch_cluster)

###############################################################################

Mu_garch_std    <- lapply(luGASRoll_GH, `[[`, "Mu")
Sigma_garch_std <- lapply(luGASRoll_GH, `[[`, "Sigma")
Shape_garch_std <- lapply(luGASRoll_GH, `[[`, "Shape")
luGASRoll_GH <- lapply(luGASRoll_GH, `[[`, "roll")

Mu_garch_norm    <- lapply(luGASRoll_GHN, `[[`, "Mu")
Sigma_garch_norm <- lapply(luGASRoll_GHN, `[[`, "Sigma")
luGASRoll_GHN <- lapply(luGASRoll_GHN, `[[`, "roll")



#GARCH#######

names(luGASRoll_GH) <- vAsset
names(Mu_garch_std)<- vAsset
names(Sigma_garch_std)<- vAsset
names(Shape_garch_std)<- vAsset

names(luGASRoll_GHN) <- vAsset
names(Mu_garch_norm)<- vAsset
names(Sigma_garch_norm)<- vAsset

#############

################## parameters GARCH ###########


################## EMPIRICAL ##################


lVaR <- list()
lES  <- list()

vDist <- c("GH","GHN")

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
      colname <- paste0("alpha(", dAlpha*100, "%)")
      mOut[, sDist, paste(dAlpha)] <- obj_GAS[[sAsset]]@forecast$VaR[, colname]
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
        
       if(sDist == "GH") {

        alpha <- dAlpha
        u_grid <- seq(1e-7, alpha, length.out = 1000)

        q_mat <- sapply(u_grid, function(u) as.numeric(quantile(luGASRoll_GH[[sAsset]], u)))
        q_mat <- t(q_mat)  # now rows 

        mOut[, sDist, paste(dAlpha)] <- apply(q_mat, 2, function(q_vals) {
          trapz(u_grid, q_vals) / alpha
        })

      }else if(sDist == "GHN"){
        alpha <- dAlpha
        u_grid <- seq(1e-7, alpha, length.out = 1000)

        q_mat <- sapply(u_grid, function(u) as.numeric(quantile(luGASRoll_GHN[[sAsset]], u)))
        q_mat <- t(q_mat)  # now rows 

        mOut[, sDist, paste(dAlpha)] <- apply(q_mat, 2, function(q_vals) {
          trapz(u_grid, q_vals) / alpha
        })
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


################## Var and ES joint Backtesting, DM test ##################
vDistPairs <- combn(vDist, 2, simplify = FALSE)
vDistPairsvs <- sapply(vDistPairs, function(pair) paste(pair[1], pair[2], sep = "_vs_"))

dmResults <- list()
dmResults2 <- list()


for (sAsset in vAsset) {
  dmResults[[sAsset]] <- list()
  dmResults2[[sAsset]] <- list()
  
  for (dAlpha in vAlpha) {
    dmResults[[sAsset]][[paste(dAlpha)]] <- list()
    dmResults2[[sAsset]][[paste(dAlpha)]] <- list()
    
    for (pair in vDistPairs) {
      model1 <- pair[1]
      model2 <- pair[2]
      
      loss1 <- lVaRESBacktest[[sAsset]][[model1]][[paste(dAlpha)]]
      loss2 <- lVaRESBacktest[[sAsset]][[model2]][[paste(dAlpha)]]
      
      # remove NA/NaN
      valid_idx <- !is.na(loss1) & !is.nan(loss1) & !is.na(loss2) & !is.nan(loss2)

      # keep only valid points in both series
      loss1_clean <- loss1[valid_idx]
      loss2_clean <- loss2[valid_idx]
      
      dmResults[[sAsset]][[paste(dAlpha)]][[paste(model1, model2, sep = "_vs_")]] <- 
        dm.test(loss1_clean, loss2_clean,
        alternative = "less",
        h = 1, power = 2) #test if model 1 is better than model 2.
      dmResults2[[sAsset]][[paste(dAlpha)]][[paste(model1, model2, sep = "_vs_")]] <- 
        dm.test(loss1_clean, loss2_clean,
        alternative = "greater",
        h = 1, power = 2) #test if model 1 is worst than model 2.
    }
  }
}

#dmResults[["AA"]][["0.01"]][["norm_vs_std"]]$statistic or dmResults[["AA"]][["0.01"]][["norm_vs_std"]]$p.value


################## Extract Losses ##################



################ TABLES ########################

vModels <- paste("GAS", vDist, sep = "-")
vTestNames <- c("LRuc", "LRcc", "DQ", "Loss-QL", "Loss-FZL", "AE", "ADmean", "ADmax" )
vESCCTestNames <- c("twosided_simple" , "onesided_simple")
vInfoCoef <- c("mean","variance")
cBackTest <- array(NA, dim = c(length(vAsset), length(vTestNames),length(vModels), length(vAlpha)),
                   dimnames = list(vAsset, vTestNames, vModels, vAlpha)) #VaR backtest 

cBackTestDM <- array(NA , dim = c(length(vAsset) , length(vAlpha) , length(vDistPairsvs)) ,
                   dimnames = list(vAsset , vAlpha , vDistPairsvs)) #DM FZloss test 1

cBackTestDM2 <- array(NA , dim = c(length(vAsset) , length(vAlpha) , length(vDistPairsvs)) ,
                   dimnames = list(vAsset , vAlpha , vDistPairsvs)) #DM FZloss test 2

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
    for (pair in vDistPairsvs){
        #DM test 1
        cBackTestDM[sAsset,  paste(dAlpha) , pair] <- dmResults[[sAsset]][[paste(dAlpha)]][[pair]]$p.value
        #DM test 2
        cBackTestDM2[sAsset,  paste(dAlpha) , pair] <- dmResults2[[sAsset]][[paste(dAlpha)]][[pair]]$p.value
    }
  }
}


#tables #######################################################################################################

sDir_Tables <- "../results/VaR_ES_Backtesting/GARCH_Tables"  
alpha1 <- paste(vAlpha[1])
alpha2 <- paste(vAlpha[2])

#LRcc p-value  ################################################################################################


mTab_LRcc <- format(round(cbind(cBackTest[, "LRcc", , alpha1], cBackTest[, "LRcc", , alpha2]), 2), digits = 4, scientific = FALSE)
mTab_LRcc_numeric <- cbind(cBackTest[, "LRcc", , alpha1], cBackTest[, "LRcc", , alpha2]) 
mTab_LRcc[mTab_LRcc_numeric >= 0.05] <- paste("\\grb{",mTab_LRcc[mTab_LRcc_numeric >= 0.05], "}", sep = "")
mTab_LRcc[, ncol(mTab_LRcc)] <- paste(mTab_LRcc[, ncol(mTab_LRcc)], "\\\\")
mTab_LRcc <- mTab_LRcc[sort(rownames(mTab_LRcc)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_LRcc), mTab_LRcc)),
            file = paste(sDir_Tables, "/mTab_LRcc.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)




#LRuc p-value ################################################################################################


mTab_LRuc <- format(round(cbind(cBackTest[, "LRuc", , alpha1], cBackTest[, "LRuc", , alpha2]), 2), digits = 4, scientific = FALSE)
mTab_LRuc_numeric <- cbind(cBackTest[, "LRuc", , alpha1], cBackTest[, "LRuc", , alpha2])
mTab_LRuc[mTab_LRuc_numeric >= 0.05] <- paste("\\grb{",mTab_LRuc[mTab_LRuc_numeric >= 0.05], "}", sep = "")
mTab_LRuc[, ncol(mTab_LRuc)] <- paste(mTab_LRuc[, ncol(mTab_LRuc)], "\\\\")
mTab_LRuc <- mTab_LRuc[sort(rownames(mTab_LRuc)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_LRuc), mTab_LRuc)),
            file = paste(sDir_Tables, "/mTab_LRuc.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

###QL loss table #################################################


mTab_QL <- format(round(cbind(cBackTest[, "Loss-QL", , alpha1], cBackTest[, "Loss-QL", , alpha2]), 4), digits = 4, scientific = FALSE) #table as text values 
mTab_QL_numeric <- cbind(cBackTest[, "Loss-QL", , alpha1], cBackTest[, "Loss-QL", , alpha2]) #table as numeric values

for (i in 1:nrow(mTab_QL)) {
  # find the index of the minimum (numeric values!)


  min_val1 <- min(round(mTab_QL_numeric[i, 1:iD], 4))
  min_val2 <- min(round(mTab_QL_numeric[i, (iD + 1):(2*iD)], 4))

# find all indices where the rounded value equals the rounded minimum
  min_indices1 <- which(round(mTab_QL_numeric[i, 1:iD], 4) == min_val1)
  min_indices2 <- which(round(mTab_QL_numeric[i, (iD + 1):(2*iD)], 4) == min_val2) + iD


  # wrap the corresponding text value with LaTeX bold command
  mTab_QL[i, min_indices1] <- paste("\\textbf{", mTab_QL[i, min_indices1], "}", sep = "")
  mTab_QL[i, min_indices2] <- paste("\\textbf{", mTab_QL[i, min_indices2], "}", sep = "")
}

mTab_QL[, ncol(mTab_QL)] <- paste(mTab_QL[, ncol(mTab_QL)], "\\\\")

mTab_QL <- mTab_QL[sort(rownames(mTab_QL)), ]


# ---- Add custom column names ----
# Extract distribution names (vDist) and alphas (vAlpha)
col_labels <- c(outer(vDist, c("0.01", "0.05"), 
                      FUN = function(d, a) paste(d, "-", a)))

# Now prepend an empty string for the "Asset" column
col_labels <- c("Asset", col_labels)

# Write file with header
write.table(rbind(col_labels, cbind(rownames(mTab_QL), mTab_QL)),
            file = paste(sDir_Tables, "/mTab_QLLoss.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)






#DQ table ###########################################################################################################
  

mTab_DQ <- format(round(cbind(cBackTest[, "DQ", , alpha1], cBackTest[, "DQ", , alpha2]), 2), digits = 4, scientific = FALSE) 
mTab_DQ_numeric <- cbind(cBackTest[, "DQ", , alpha1], cBackTest[, "DQ", , alpha2]) 
mTab_DQ[mTab_DQ_numeric >= 0.05] <- paste("\\grb{",mTab_DQ[mTab_DQ_numeric >= 0.05], "}", sep = "")
mTab_DQ[, ncol(mTab_DQ)] <- paste(mTab_DQ[, ncol(mTab_DQ)], "\\\\")
mTab_DQ <- mTab_DQ[sort(rownames(mTab_DQ)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"),FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_DQ), mTab_DQ)),
            file = paste(sDir_Tables, "/mTab_DQ.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

## DM p-values #####################################################################################################

mTab_DM <- format(round(cbind(cBackTestDM[ , alpha1 , ], cBackTestDM[ , alpha2 , ]), 2), digits = 4, scientific = FALSE) 
mTab_DM_numeric <- cbind(cBackTestDM[ , alpha1 , ],cBackTestDM[ , alpha2 , ])
mTab_DM[mTab_DM_numeric <= 0.05] <- paste("\\grb{",mTab_DM[mTab_DM_numeric <= 0.05], "}", sep = "")
mTab_DM[, ncol(mTab_DM)] <- paste(mTab_DM[, ncol(mTab_DM)], "\\\\")
mTab_DM <- mTab_DM[sort(rownames(mTab_DM)), ]

col_labels <- c(outer(vDistPairsvs, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_DM), mTab_DM)),
            file = paste(sDir_Tables, "/mTab_DM.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)




mTab_DM2 <- format(round(cbind(cBackTestDM2[ , alpha1 , ], cBackTestDM2[ , alpha2 , ]), 2), digits = 4, scientific = FALSE) 
mTab_DM_numeric2 <- cbind(cBackTestDM2[ , alpha1 , ],cBackTestDM2[ , alpha2 , ])
mTab_DM2[mTab_DM_numeric2 <= 0.05] <- paste("\\grb{",mTab_DM2[mTab_DM_numeric2 <= 0.05], "}", sep = "")
mTab_DM2[, ncol(mTab_DM2)] <- paste(mTab_DM2[, ncol(mTab_DM2)], "\\\\")
mTab_DM2 <- mTab_DM2[sort(rownames(mTab_DM2)), ]

col_labels <- c(outer(vDistPairsvs, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_DM2), mTab_DM2)),
            file = paste(sDir_Tables, "/mTab_DM2.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)



# average FZ loss ###################################################################################################




mTab_FZloss <- format(round(cbind(cBackTest[, "Loss-FZL", , alpha1], cBackTest[, "Loss-FZL", , alpha2]), 2), digits = 4, scientific = FALSE) #table as text values 
mTab_FZloss_numeric <- cbind(cBackTest[, "Loss-FZL", , alpha1], cBackTest[, "Loss-FZL", , alpha2]) #table as numeric values

for (i in 1:nrow(mTab_FZloss)) {
  # find the minimum value for each block
  min_val1 <- min(mTab_FZloss_numeric[i, 1:iD])
  min_val2 <- min(mTab_FZloss_numeric[i, (iD + 1):(2*iD)])
  
  # find all indices that equal the minimum
  min_indices1 <- which(mTab_FZloss_numeric[i, 1:iD] == min_val1)
  min_indices2 <- which(mTab_FZloss_numeric[i, (iD + 1):(2*iD)] == min_val2) + iD
  
  # wrap the corresponding text values with LaTeX bold
  mTab_FZloss[i, min_indices1] <- paste0("\\textbf{", mTab_FZloss[i, min_indices1], "}")
  mTab_FZloss[i, min_indices2] <- paste0("\\textbf{", mTab_FZloss[i, min_indices2], "}")
}

mTab_FZloss[, ncol(mTab_FZloss)] <- paste(mTab_FZloss[, ncol(mTab_FZloss)], "\\\\")
mTab_FZloss <- mTab_FZloss[sort(rownames(mTab_FZloss)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_FZloss), mTab_FZloss)),
            file = paste(sDir_Tables, "/mTab_FZloss.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)



###AD table ##########################################################


mTab_AD <- format(round(cbind(cBackTest[, "ADmean", , alpha1], cBackTest[, "ADmean", , alpha2]), 2), digits = 4, scientific = FALSE) 
mTab_AD[, ncol(mTab_AD)] <- paste(mTab_AD[, ncol(mTab_AD)], "\\\\")
mTab_AD <- mTab_AD[sort(rownames(mTab_AD)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"),FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_AD), mTab_AD)),
            file = paste(sDir_Tables, "/mTab_AD.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

###AE table ##########################################################


mTab_AE <- format(round(cbind(cBackTest[, "AE", , alpha1], cBackTest[, "AE", , alpha2]), 2), digits = 4, scientific = FALSE) #table as text values 
mTab_AE[, ncol(mTab_AE)] <- paste(mTab_AE[, ncol(mTab_AE)], "\\\\")
mTab_AE <- mTab_AE[sort(rownames(mTab_AE)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"),FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_AE), mTab_AE)),
            file = paste(sDir_Tables, "/mTab_AE.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

######CC ES backtesting #####################################################

mTab_EScc <- format(round(cbind(cBackTestCC[, "twosided_simple", , alpha1], cBackTestCC[, "twosided_simple", , alpha2]), 2), digits = 4, scientific = FALSE) 
mTab_EScc_numeric <- cbind(cBackTestCC[, "twosided_simple", , alpha1], cBackTestCC[, "twosided_simple", , alpha2]) 
mTab_EScc[mTab_EScc_numeric >= 0.05] <- paste("\\grb{",mTab_EScc[mTab_EScc_numeric >= 0.05], "}", sep = "")
mTab_EScc[, ncol(mTab_EScc)] <- paste(mTab_EScc[, ncol(mTab_EScc)], "\\\\")
mTab_EScc <- mTab_EScc[sort(rownames(mTab_EScc)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_EScc), mTab_EScc)),
            file = paste(sDir_Tables, "/mTab_EScc.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)


##############ERS ES backtesting ####################################################

mTab_ESers <- format(round(cbind(cBackTestESER[, , alpha1], cBackTestESER[, , alpha2]), 2), digits = 4, scientific = FALSE) 
mTab_ESers_numeric <- cbind(cBackTestESER[, , alpha1], cBackTestESER[, , alpha2])
mTab_ESers[mTab_ESers_numeric >= 0.05] <- paste("\\grb{",mTab_ESers[mTab_ESers_numeric >= 0.05], "}", sep = "")
mTab_ESers[, ncol(mTab_ESers)] <- paste(mTab_ESers[, ncol(mTab_ESers)], "\\\\")
mTab_ESers <- mTab_ESers[sort(rownames(mTab_ESers)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_ESers), mTab_ESers)),
            file = paste(sDir_Tables, "/mTab_ESers.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ES E-backtesting ######################################################
mTab_ESet <- format(round(cbind(cBackTestESET[, , alpha1], cBackTestESET[, , alpha2]), 2), digits = 4, scientific = FALSE) 
mTab_ESet[, ncol(mTab_ESet)] <- paste(mTab_ESet[, ncol(mTab_ESet)], "\\\\")
mTab_ESet <- mTab_ESet[sort(rownames(mTab_ESet)), ]

col_labels <- c(outer(vDist, c("0.01", "0.05"), FUN = function(d, a) paste(d, "-", a)))
col_labels <- c("Asset", col_labels)

write.table(rbind(col_labels, cbind(rownames(mTab_ESet), mTab_ESet)),
            file = paste(sDir_Tables, "/mTab_ESet.txt", sep = ""),
            sep = " & ", quote = FALSE, row.names = FALSE, col.names = FALSE)

