library(dplyr)
library(tidystats)
library(knitr)
library(kableExtra)

source("experiments.R")


##
# Arguments:
#  studiesAnalysisDfPath (string): path where exact analysis was done using web app
#    to determine whether a result was a false positive or not.
#  
makePlausibleFPTable <- function(studiesAnalysisDfPath = "data/StudiesAnalysis.csv", 
                                 useSynced = FALSE) 
{
  plausibleFPTable <-
    read.csv(studiesAnalysisDfPath) %>%
    filter(IncludeInt == 1) %>% 
    group_by(ArticleTag) %>%
    summarize(PlausibleFPs = sum(PlausibleInt), ConditionCount = n(), PlausibleFPRate = mean(PlausibleInt), )
  
  write.csv(plausibleFPTable, file = "data/output/PlausibleFPTable.csv", row.names = FALSE)
  
  return (plausibleFPTable)
}


##
# Generate synthetic data based using parameters from experimental conditions
# that yielded plausibly false positives to test how frequently simulated simple
# consensus data yields false positives.
#
# Arguments:
#  studiesAnalysisDfPath (string): path where exact analysis was done using web app
#    to determine whether a result was a false positive or not.
#  ntrials (int): Number of synthetic datasets to generate and t-test for each of the
#    treatments yielding plausibly false positive results.
#  limit (int): Limit the number of treatments to test for initial development.
#
makeTtestSimulationsTable <- function(studiesAnalysisDfPath = "data/StudiesAnalysis.csv",
                                      outputPath = "data/output/TtestFitTable.csv",
                                      diagnosticSavePath = "data/diagnostic/tTestFits.RDS",
                                      ntrials = 2, limit = 2)
{
  # Load analysis of studies from web interface stored as CSV.
  studiesDf <- read.csv(studiesAnalysisDfPath)

  # Filter out implausible false positives.
  studiesDf <- filter(studiesDf, (IncludeInt == 1) & (PlausibleInt == 1))
  
  # Limit number of plausible false positives for preliminary development.
  if (limit > 0)
  {
    studiesDf <- head(studiesDf, limit)
  }
  
  ## Apply t-tests to generated data over n trials.
  # Initialize results dataframe with same length of studiesDf and appropriate
  # additional columns.
  resultsDf <- data.frame(ArticleTag = c(), 
                          TreatmentTag = c(),
                          TrialIndex = c(),
                          tTestPvalue = c(),
                          ExpectedPower = c())
  
  allTestsDiagnostic <- list()
  for (rowIdx in 1:nrow(studiesDf))
  {
    # Get treatment row of interest.
    row <- studiesDf[rowIdx, ]
    
    # Calculate number of bins from min and max opinion bin value.
    nBins <- row$MaxBinValue - row$MinBinValue + 1
    
    # Assemble trials dataframe.
    trialsDf <- data.frame(ArticleTag = c(), TreatmentTag = c(), 
                           TrialIndex = c(), LatentMean = c(), LatentSDPre = c(),
                           LatentSDPost = c(), tTestPvalue = c(), ObservedShift = c(),
                           ExpectedPower = c())
    
    expectedPower = power.t.test(delta = row$ObservedShift, n = row$N, 
                                 sd = (row$LatentSDPost + row$LatentSDPre) / 2.0,
                                 sig.level = 0.1)$power
                                 
    for (trialIdx in 1:ntrials)
    {
      # Run t-test ntrials experiments for treatment row and trial index.
      tTestResult <- tTestExperiment(row$N, row$MinBinValue, nBins, 
                                     row$LatentMean, row$LatentSDPre, 
                                     row$LatentSDPost)
      # return (tTestResult$p.value)
      trial_row = data.frame(ArticleTag = row$ArticleTag, TreatmentTag = row$TreatmentTag,
                             TrialIndex = trialIdx, LatentMean = row$LatentMean,
                             ObservedShift = row$ObservedShift,
                             LatentSDPre = row$LatentSDPre, LatentSDPost = row$LatentSDPost,
                             tTestPvalue = tTestResult$p.value,
                             ExpectedPower = expectedPower)
      
      resultsDf <- bind_rows(resultsDf, trial_row)
      
      # Create a test ID for this treatment and add the test output to diagnostic set.
      testId <- paste(row$ArticleTag, row$TreatmentTag, trialIdx, sep = "_")
      allTestsDiagnostic <- allTestsDiagnostic %>% add_stats(tTestResult, identifier = testId)
    }
  }
  
  write.csv(resultsDf, outputPath)
  saveRDS(allTestsDiagnostic, file = diagnosticSavePath)
  
  # return (resultsDf)
}

summarizeTTestFitTable <- function(fitTablePath = "data/output/TtestFitTable.csv",
                                   summaryTTablePath = "data/output/TtestSummaryTable.csv",
                                   significanceVal = 0.1)
{
  fitTableDf <- 
    read.csv(fitTablePath) %>%
    group_by(ArticleTag, TreatmentTag, ObservedShift, LatentSDPre, LatentSDPost, ExpectedPower) %>%
    summarize(FractionSignificant = mean(tTestPvalue < significanceVal), .groups='keep')
  
  write.csv(fitTableDf, summaryTTablePath, row.names = FALSE)
  
  return (fitTableDf)
}


latexifyPlausibleFPTable <- function(plausibleFPTablePath = "data/output/PlausibleFPTable.csv",
                                     latexifiedPath = "papers/plausibleFP-table.tex")
{
    plausibleFP_table <- read.csv(plausibleFPTablePath)
    plausibleFP_latex <- 
        kable(plausibleFP_table, "latex", longtable = T, booktabs = T, digits=2) %>%
        kable_styling(latex_options = c("repeat_header"))

    writeLines(plausibleFP_latex, latexifiedPath) 
}


latexifyTTestExperiment <- function(tTestTablePath = "data/output/TtestSummaryTable.csv",
                                    latexifiedPath = "paper/t-test-table.tex")
{
    tTestTable <- read.csv(tTestTablePath)
    tTestTable$Article_Condition <- paste(tTestTable$ArticleTag, tTestTable$TreatmentTag, sep = "_")
    
    # tTestTable <- subset(tTestTable, select = -c(ExpectedPower, ArticleTag, TreatmentTag))
    tTestTable <- subset(tTestTable, select = -c(ArticleTag, TreatmentTag))
    
    col_order <- c("Article_Condition", "LatentMean", "LatentSDPre", 
                   "LatentSDPost", "FractionSignificant")
    col_order <- c("Article_Condition", "ObservedShift", "LatentSDPre", 
                   "LatentSDPost", "ExpectedPower", "FractionSignificant")
    
    
    tTestTable <- tTestTable[, col_order]
    
    tTestTable <- tTestTable %>% rename(FPRate = FractionSignificant,
                                        Est.FPRate = ExpectedPower)
    
    ttestTable_latex <- 
        kable(tTestTable, "latex", align = "lccccc", longtable = T, booktabs = T, digits=2) %>%
        kable_styling(latex_options = c("repeat_header"))

    writeLines(ttestTable_latex, latexifiedPath)
}


doAnalyses <- function(studiesAnalysisDfPath = "data/StudiesAnalysis.csv")
{
  
}
