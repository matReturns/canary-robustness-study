



# ============================================================
# Canary Robustness Study
# Testing canary pairs across multiple canary momentum windows
# ============================================================

library(quantmod)
library(PerformanceAnalytics)
library(xts)
library(zoo)

stratStats <- function(rets, digits = 4) {
  
  rets <- na.omit(rets)
  annual <- apply.yearly(rets, Return.cumulative)
  
  stats <- rbind(
    Return.annualized(rets),
    StdDev.annualized(rets),
    SharpeRatio.annualized(rets, Rf = 0),
    maxDrawdown(rets),
    CalmarRatio(rets),
    apply(annual, 2, min, na.rm = TRUE)
  )
  
  rownames(stats) <- c(
    "Annualized Return",
    "Annualized Volatility",
    "Annualized Sharpe",
    "Worst Drawdown",
    "Calmar Ratio",
    "Worst Calendar Year"
  )
  
  return(round(stats, digits))
}




run_canary_robustness_test <- function(dataStartDate = "2006-01-01",
                                       analysisStartDate = "2007-07-01",
                                       endDate = Sys.Date(),
                                       investableAssets = c("SPY", "EFA", "EEM", "VNQ", "DBC"),
                                       canaryPairs = list(
                                         VWO_BND = c("VWO", "BND"),
                                         EFA_BND = c("EFA", "BND"),
                                         EEM_BND = c("EEM", "BND"),
                                         HYG_BND = c("HYG", "BND"),
                                         DBC_BND = c("DBC", "BND"),
                                         SPY_BND = c("SPY", "BND")
                                       ),
                                       canaryWindowSets = list(
                                         W10_20 = c(10, 20),
                                         W20_60 = c(20, 60),
                                         W21_63 = c(21, 63),
                                         W63_126 = c(63, 126)
                                       ),
                                       nAssets = 2,
                                       assetMomWindows = c(63, 126),
                                       riskOffMode = c("conditional_defensive", "defensive", "cash"),
                                       defensiveAsset = "IEF",
                                       defensiveMomWindows = c(20, 60),
                                       cashAsset = NULL,
                                       cashReturn = 0,
                                       rebalanceOn = "months",
                                       verbose = TRUE) {
  
  riskOffMode <- match.arg(riskOffMode)
  
  # -----------------------------
  # Helper functions
  # -----------------------------
  cumulative_return <- function(x) {
    prod(1 + as.numeric(x), na.rm = FALSE) - 1
  }
  
  momentum_score <- function(R, symbol, endRow, windows) {
    
    scores <- sapply(windows, function(w) {
      
      startRow <- endRow - w + 1
      
      if (startRow < 1) {
        return(NA_real_)
      }
      
      cumulative_return(R[startRow:endRow, symbol])
    })
    
    sum(scores, na.rm = FALSE)
  }
  
  get_year_return <- function(rets, year = "2022") {
    
    yearRange <- paste0(year, "-01-01/", year, "-12-31")
    yr <- rets[yearRange]
    
    if (NROW(yr) == 0) {
      return(NA_real_)
    }
    
    return(as.numeric(Return.cumulative(yr)))
  }
  
  # -----------------------------
  # Symbols needed
  # -----------------------------
  canarySymbols <- unique(unlist(canaryPairs))
  
  downloadSymbols <- unique(c(
    investableAssets,
    canarySymbols,
    defensiveAsset,
    cashAsset
  ))
  
  downloadSymbols <- downloadSymbols[!is.na(downloadSymbols)]
  
  if (is.null(cashAsset)) {
    cashName <- "CASH"
  } else {
    cashName <- cashAsset
  }
  
  tradeAssets <- unique(c(
    investableAssets,
    defensiveAsset,
    cashName
  ))
  
  # -----------------------------
  # Download adjusted prices
  # -----------------------------
  priceList <- list()
  
  for (sym in downloadSymbols) {
    
    if (verbose) {
      message("Downloading: ", sym)
    }
    
    tmp <- getSymbols(
      sym,
      from = dataStartDate,
      to = endDate,
      auto.assign = FALSE,
      warnings = FALSE
    )
    
    px <- Ad(tmp)
    colnames(px) <- sym
    priceList[[sym]] <- px
  }
  
  prices <- do.call(merge, priceList)
  prices <- na.omit(prices)
  
  assetReturns <- Return.calculate(prices)
  assetReturns <- na.omit(assetReturns)
  
  if (is.null(cashAsset)) {
    
    cashReturns <- xts(
      rep(cashReturn / 252, NROW(assetReturns)),
      order.by = index(assetReturns)
    )
    
    colnames(cashReturns) <- "CASH"
    
    allReturns <- merge(assetReturns, cashReturns)
    
  } else {
    
    allReturns <- assetReturns
  }
  
  allReturns <- na.omit(allReturns)
  
  if (verbose) {
    message("Return data begins: ", as.character(first(index(allReturns))))
    message("Return data ends:   ", as.character(last(index(allReturns))))
  }
  
  # -----------------------------
  # Rebalance dates
  # -----------------------------
  rebalanceRows <- endpoints(allReturns, on = rebalanceOn)
  rebalanceRows <- rebalanceRows[rebalanceRows > 0]
  rebalanceRows <- rebalanceRows[rebalanceRows <= NROW(allReturns)]
  
  maxLookback <- max(
    assetMomWindows,
    unlist(canaryWindowSets),
    defensiveMomWindows
  )
  
  signalRows <- rebalanceRows[rebalanceRows > maxLookback]
  
  # ============================================================
  # Baseline 1: Equal Weight
  # ============================================================
  equalWeight <- Return.portfolio(
    R = allReturns[, investableAssets],
    weights = rep(1 / length(investableAssets), length(investableAssets)),
    rebalance_on = rebalanceOn
  )
  
  colnames(equalWeight) <- "EqualWeight"
  
  # ============================================================
  # Baseline 2: Top-N Momentum without Canary
  # ============================================================
  run_topn_momentum <- function() {
    
    weightsMat <- matrix(
      NA_real_,
      nrow = NROW(allReturns),
      ncol = length(investableAssets)
    )
    
    colnames(weightsMat) <- investableAssets
    
    for (signalRow in signalRows) {
      
      assetScores <- sapply(investableAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, assetMomWindows)
      })
      
      rankedAssets <- names(sort(assetScores, decreasing = TRUE))
      selectedAssets <- rankedAssets[1:nAssets]
      
      w <- rep(0, length(investableAssets))
      names(w) <- investableAssets
      w[selectedAssets] <- 1 / nAssets
      
      weightsMat[signalRow, ] <- w[investableAssets]
    }
    
    weights <- xts(weightsMat, order.by = index(allReturns))
    weights <- na.locf(weights, na.rm = FALSE)
    weightsLag <- lag(weights, k = 1)
    
    rets <- xts(
      rowSums(weightsLag * allReturns[, investableAssets], na.rm = FALSE),
      order.by = index(allReturns)
    )
    
    colnames(rets) <- paste0("Top", nAssets, "_Momentum")
    
    return(list(
      returns = rets,
      weights = weights,
      weightsLag = weightsLag
    ))
  }
  
  topN <- run_topn_momentum()
  topNReturns <- topN$returns
  
  # ============================================================
  # Run one canary pair / one window set
  # ============================================================
  run_one_canary <- function(pairName, canaryAssets, windowName, canaryMomWindows) {
    
    strategyName <- paste(pairName, windowName, sep = "_")
    
    weightsMat <- matrix(
      NA_real_,
      nrow = NROW(allReturns),
      ncol = length(tradeAssets)
    )
    
    colnames(weightsMat) <- tradeAssets
    
    riskBudgetVec <- rep(NA_real_, NROW(allReturns))
    
    signalRecords <- vector("list", length(signalRows))
    recordCounter <- 1
    
    for (signalRow in signalRows) {
      
      signalDate <- index(allReturns)[signalRow]
      
      # -----------------------------
      # Rank investable assets by momentum
      # -----------------------------
      assetScores <- sapply(investableAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, assetMomWindows)
      })
      
      rankedAssets <- names(sort(assetScores, decreasing = TRUE))
      selectedAssets <- rankedAssets[1:nAssets]
      
      # -----------------------------
      # Canary risk budget
      # -----------------------------
      canaryScores <- sapply(canaryAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, canaryMomWindows)
      })
      
      riskBudget <- mean(canaryScores > 0)
      riskOffWeight <- 1 - riskBudget
      
      # -----------------------------
      # Start weights
      # -----------------------------
      w <- rep(0, length(tradeAssets))
      names(w) <- tradeAssets
      
      # Risk-on allocation
      w[selectedAssets] <- riskBudget / nAssets
      
      # Risk-off allocation
      if (riskOffWeight > 0) {
        
        if (riskOffMode == "cash") {
          
          w[cashName] <- w[cashName] + riskOffWeight
          
        } else if (riskOffMode == "defensive") {
          
          w[defensiveAsset] <- w[defensiveAsset] + riskOffWeight
          
        } else if (riskOffMode == "conditional_defensive") {
          
          defensiveScore <- momentum_score(
            allReturns,
            defensiveAsset,
            signalRow,
            defensiveMomWindows
          )
          
          if (!is.na(defensiveScore) && defensiveScore > 0) {
            w[defensiveAsset] <- w[defensiveAsset] + riskOffWeight
          } else {
            w[cashName] <- w[cashName] + riskOffWeight
          }
        }
      }
      
      weightsMat[signalRow, ] <- w[tradeAssets]
      riskBudgetVec[signalRow] <- riskBudget
      
      signalRecords[[recordCounter]] <- data.frame(
        date = as.Date(signalDate),
        strategy = strategyName,
        canaryPair = pairName,
        windowSet = windowName,
        canaryAssets = paste(canaryAssets, collapse = ", "),
        canaryWindows = paste(canaryMomWindows, collapse = ", "),
        selectedAssets = paste(selectedAssets, collapse = ", "),
        riskBudget = riskBudget,
        riskOffWeight = riskOffWeight,
        stringsAsFactors = FALSE
      )
      
      recordCounter <- recordCounter + 1
    }
    
    # -----------------------------
    # Convert signal weights to daily weights
    # -----------------------------
    weights <- xts(weightsMat, order.by = index(allReturns))
    weights <- na.locf(weights, na.rm = FALSE)
    weightsLag <- lag(weights, k = 1)
    
    # Strategy returns
    rets <- xts(
      rowSums(weightsLag * allReturns[, tradeAssets], na.rm = FALSE),
      order.by = index(allReturns)
    )
    
    colnames(rets) <- strategyName
    
    retsLive <- rets[paste0(analysisStartDate, "/")]
    retsLive <- na.omit(retsLive)
    
    liveWeights <- weightsLag[index(retsLive), tradeAssets]
    
    riskBudget <- xts(riskBudgetVec, order.by = index(allReturns))
    riskBudget <- na.locf(riskBudget, na.rm = FALSE)
    riskBudgetLag <- lag(riskBudget, k = 1)
    liveRiskBudget <- riskBudgetLag[index(retsLive)]
    
    signalLog <- do.call(rbind, signalRecords)
    
    # -----------------------------
    # Exposure stats
    # -----------------------------
    exposure <- data.frame(
      Strategy = strategyName,
      CanaryPair = pairName,
      WindowSet = windowName,
      CanaryWindows = paste(canaryMomWindows, collapse = "_"),
      avgRiskBudget = mean(as.numeric(liveRiskBudget), na.rm = TRUE),
      pctFullRiskOn = mean(as.numeric(liveRiskBudget) == 1, na.rm = TRUE),
      pctHalfRiskOn = mean(as.numeric(liveRiskBudget) == 0.5, na.rm = TRUE),
      pctRiskOff = mean(as.numeric(liveRiskBudget) == 0, na.rm = TRUE),
      avgRiskAssetWeight = mean(rowSums(liveWeights[, investableAssets], na.rm = TRUE), na.rm = TRUE),
      avgDefensiveWeight = mean(as.numeric(liveWeights[, defensiveAsset]), na.rm = TRUE),
      avgCashWeight = mean(as.numeric(liveWeights[, cashName]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    exposure[, c(
      "avgRiskBudget",
      "pctFullRiskOn",
      "pctHalfRiskOn",
      "pctRiskOff",
      "avgRiskAssetWeight",
      "avgDefensiveWeight",
      "avgCashWeight"
    )] <- round(exposure[, c(
      "avgRiskBudget",
      "pctFullRiskOn",
      "pctHalfRiskOn",
      "pctRiskOff",
      "avgRiskAssetWeight",
      "avgDefensiveWeight",
      "avgCashWeight"
    )], 4)
    
    return(list(
      returns = retsLive,
      weights = weights,
      weightsLag = weightsLag,
      riskBudget = riskBudget,
      signalLog = signalLog,
      exposure = exposure
    ))
  }
  
  # ============================================================
  # Run all canary pairs across all window sets
  # ============================================================
  results <- list()
  returnsList <- list()
  exposureList <- list()
  signalLogList <- list()
  
  for (pairName in names(canaryPairs)) {
    
    for (windowName in names(canaryWindowSets)) {
      
      if (verbose) {
        message("Running: ", pairName, " | ", windowName)
      }
      
      res <- run_one_canary(
        pairName = pairName,
        canaryAssets = canaryPairs[[pairName]],
        windowName = windowName,
        canaryMomWindows = canaryWindowSets[[windowName]]
      )
      
      strategyName <- paste(pairName, windowName, sep = "_")
      
      results[[strategyName]] <- res
      returnsList[[strategyName]] <- res$returns
      exposureList[[strategyName]] <- res$exposure
      signalLogList[[strategyName]] <- res$signalLog
    }
  }
  
  robustnessReturns <- do.call(merge, returnsList)
  exposureStats <- do.call(rbind, exposureList)
  signalLog <- do.call(rbind, signalLogList)
  
  rownames(exposureStats) <- NULL
  rownames(signalLog) <- NULL
  
  # -----------------------------
  # Combine baselines and robustness returns
  # -----------------------------
  mainReturns <- merge(
    equalWeight,
    topNReturns,
    robustnessReturns
  )
  
  mainReturns <- mainReturns[paste0(analysisStartDate, "/")]
  mainReturns <- na.omit(mainReturns)
  
  robustnessReturns <- robustnessReturns[index(mainReturns)]
  robustnessReturns <- na.omit(robustnessReturns)
  
  # -----------------------------
  # Summary tables
  # -----------------------------
  summary <- stratStats(mainReturns)
  robustnessSummary <- stratStats(robustnessReturns)
  
  summaryTable <- data.frame(
    Strategy = colnames(mainReturns),
    AnnualizedReturn = as.numeric(summary["Annualized Return", ]),
    AnnualizedVolatility = as.numeric(summary["Annualized Volatility", ]),
    Sharpe = as.numeric(summary["Annualized Sharpe", ]),
    WorstDrawdown = as.numeric(summary["Worst Drawdown", ]),
    Calmar = as.numeric(summary["Calmar Ratio", ]),
    WorstCalendarYear = as.numeric(summary["Worst Calendar Year", ]),
    stringsAsFactors = FALSE
  )
  
  summaryTable$Return2022 <- sapply(
    summaryTable$Strategy,
    function(s) get_year_return(mainReturns[, s], "2022")
  )
  
  summaryTable$Return2022 <- round(summaryTable$Return2022, 4)
  
  # Add parsed canary pair and window set
  summaryTable$CanaryPair <- NA_character_
  summaryTable$WindowSet <- NA_character_
  
  for (i in seq_len(nrow(summaryTable))) {
    
    strategyName <- summaryTable$Strategy[i]
    
    if (strategyName %in% names(results)) {
      summaryTable$CanaryPair[i] <- results[[strategyName]]$exposure$CanaryPair
      summaryTable$WindowSet[i] <- results[[strategyName]]$exposure$WindowSet
    }
  }
  
  # Put columns in a nicer order
  summaryTable <- summaryTable[, c(
    "Strategy",
    "CanaryPair",
    "WindowSet",
    "AnnualizedReturn",
    "AnnualizedVolatility",
    "Sharpe",
    "WorstDrawdown",
    "Calmar",
    "WorstCalendarYear",
    "Return2022"
  )]
  
  rankedBySharpe <- summaryTable[order(-summaryTable$Sharpe), ]
  rankedByCalmar <- summaryTable[order(-summaryTable$Calmar), ]
  rankedByDrawdown <- summaryTable[order(summaryTable$WorstDrawdown), ]
  rankedByReturn <- summaryTable[order(-summaryTable$AnnualizedReturn), ]
  
  rownames(rankedBySharpe) <- NULL
  rownames(rankedByCalmar) <- NULL
  rownames(rankedByDrawdown) <- NULL
  rownames(rankedByReturn) <- NULL
  
  annualReturns <- apply.yearly(mainReturns, Return.cumulative)
  
  # -----------------------------
  # Canary-pair averages across windows
  # -----------------------------
  canaryOnly <- subset(summaryTable, !is.na(CanaryPair))
  
  pairSummary <- aggregate(
    cbind(
      AnnualizedReturn,
      AnnualizedVolatility,
      Sharpe,
      WorstDrawdown,
      Calmar,
      WorstCalendarYear,
      Return2022
    ) ~ CanaryPair,
    data = canaryOnly,
    FUN = mean
  )
  
  pairSummary <- pairSummary[order(-pairSummary$Sharpe), ]
  rownames(pairSummary) <- NULL
  
  # -----------------------------
  # Window-set averages across canary pairs
  # -----------------------------
  windowSummary <- aggregate(
    cbind(
      AnnualizedReturn,
      AnnualizedVolatility,
      Sharpe,
      WorstDrawdown,
      Calmar,
      WorstCalendarYear,
      Return2022
    ) ~ WindowSet,
    data = canaryOnly,
    FUN = mean
  )
  
  windowSummary <- windowSummary[order(-windowSummary$Sharpe), ]
  rownames(windowSummary) <- NULL
  
  # -----------------------------
  # Heatmap-friendly tables
  # -----------------------------
  sharpeMatrix <- xtabs(
    Sharpe ~ CanaryPair + WindowSet,
    data = canaryOnly
  )
  
  calmarMatrix <- xtabs(
    Calmar ~ CanaryPair + WindowSet,
    data = canaryOnly
  )
  
  drawdownMatrix <- xtabs(
    WorstDrawdown ~ CanaryPair + WindowSet,
    data = canaryOnly
  )
  
  
  # -----------------------------
  # Round numeric columns only
  # -----------------------------
  round_numeric_df <- function(df, digits = 4) {
    
    numericCols <- sapply(df, is.numeric)
    df[numericCols] <- lapply(df[numericCols], round, digits)
    
    return(df)
  }
  
  
  summaryTable <- round_numeric_df(summaryTable, 4)
  rankedBySharpe <- round_numeric_df(rankedBySharpe, 4)
  rankedByCalmar <- round_numeric_df(rankedByCalmar, 4)
  rankedByDrawdown <- round_numeric_df(rankedByDrawdown, 4)
  rankedByReturn <- round_numeric_df(rankedByReturn, 4)
  pairSummary <- round_numeric_df(pairSummary, 4)
  windowSummary <- round_numeric_df(windowSummary, 4)
  exposureStats <- round_numeric_df(exposureStats, 4)
  return(list(
    mainReturns = mainReturns,
    robustnessReturns = robustnessReturns,
    summary = summary,
    robustnessSummary = robustnessSummary,
    summaryTable = summaryTable,
    rankedBySharpe = rankedBySharpe,
    rankedByCalmar = rankedByCalmar,
    rankedByDrawdown = rankedByDrawdown,
    rankedByReturn = rankedByReturn,
    pairSummary = pairSummary,
    windowSummary = windowSummary,
    exposureStats = exposureStats,
    annualReturns = round(annualReturns, 4),
    sharpeMatrix = round(sharpeMatrix, 4),
    calmarMatrix = round(calmarMatrix, 4),
    drawdownMatrix = round(drawdownMatrix, 4),
    results = results,
    signalLog = signalLog,
    prices = prices,
    allReturns = allReturns,
    settings = list(
      dataStartDate = dataStartDate,
      analysisStartDate = as.character(first(index(mainReturns))),
      endDate = as.character(last(index(mainReturns))),
      investableAssets = investableAssets,
      canaryPairs = canaryPairs,
      canaryWindowSets = canaryWindowSets,
      nAssets = nAssets,
      assetMomWindows = assetMomWindows,
      riskOffMode = riskOffMode,
      defensiveAsset = defensiveAsset,
      defensiveMomWindows = defensiveMomWindows,
      cashAsset = cashAsset,
      cashReturn = cashReturn,
      rebalanceOn = rebalanceOn
    )
  ))
}




canaryRobustness <- run_canary_robustness_test(
  dataStartDate = "2006-01-01",
  analysisStartDate = "2007-07-01",
  endDate = "2026-07-09",
  investableAssets = c("SPY", "EFA", "EEM", "VNQ", "DBC"),
  canaryPairs = list(
    VWO_BND = c("VWO", "BND"),
    EFA_BND = c("EFA", "BND"),
    EEM_BND = c("EEM", "BND"),
    HYG_BND = c("HYG", "BND"),
    DBC_BND = c("DBC", "BND"),
    SPY_BND = c("SPY", "BND")
  ),
  canaryWindowSets = list(
    W10_20 = c(10, 20),
    W20_60 = c(20, 60),
    W21_63 = c(21, 63),
    W63_126 = c(63, 126)
  ),
  nAssets = 2,
  assetMomWindows = c(63, 126),
  riskOffMode = "conditional_defensive",
  defensiveAsset = "IEF",
  defensiveMomWindows = c(20, 60),
  cashAsset = NULL,
  cashReturn = 0,
  rebalanceOn = "months",
  verbose = TRUE
)




canaryRobustness$summaryTable
canaryRobustness$rankedBySharpe
canaryRobustness$rankedByCalmar
canaryRobustness$rankedByDrawdown
canaryRobustness$rankedByReturn

canaryRobustness$pairSummary
canaryRobustness$windowSummary

canaryRobustness$sharpeMatrix
canaryRobustness$calmarMatrix
canaryRobustness$drawdownMatrix

canaryRobustness$exposureStats
canaryRobustness$annualReturns
canaryRobustness$settings




top6SharpeNames <- canaryRobustness$rankedBySharpe$Strategy[1:6]

charts.PerformanceSummary(
  canaryRobustness$mainReturns[, top6SharpeNames],
  main = "Top Canary Signal Robustness Results",
  wealth.index = T
)




top6CalmarNames <- canaryRobustness$rankedByCalmar$Strategy[1:6]

charts.PerformanceSummary(
  canaryRobustness$mainReturns[, top6CalmarNames],
  main = "Top Canary Robustness Results by Calmar Ratio",
  wealth.index = T
)




barplot(
  canaryRobustness$pairSummary$Sharpe,
  names.arg = canaryRobustness$pairSummary$CanaryPair,
  las = 2,
  main = "Average Sharpe by Canary Pair",
  ylab = "Average Sharpe Across Window Sets"
)




barplot(
  canaryRobustness$windowSummary$Sharpe,
  names.arg = canaryRobustness$windowSummary$WindowSet,
  las = 2,
  main = "Average Sharpe by Canary Momentum Window",
  ylab = "Average Sharpe Across Canary Pairs"
)




barplot(
  canaryRobustness$pairSummary$WorstDrawdown,
  names.arg = canaryRobustness$pairSummary$CanaryPair,
  las = 2,
  main = "Average Worst Drawdown by Canary Pair",
  ylab = "Average Worst Drawdown"
)




canaryRobustness$pairSummary
canaryRobustness$windowSummary
canaryRobustness$rankedBySharpe
canaryRobustness$rankedByCalmar
canaryRobustness$sharpeMatrix
canaryRobustness$calmarMatrix































