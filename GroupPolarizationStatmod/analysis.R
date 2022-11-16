library(dplyr)

##
# Arguments:
#  studiesAnalysisDfPath (string): path where exact analysis was done using web app
#    to determine whether a result was a false positive or not.
#  
makePlausibleFPTable <- function(studiesAnalysisDfPath = "data/StudiesAnalysis.csv", useSynced = FALSE) 
{
  plausibleFPTable <-
    read.csv(studiesAnalysisDfPath) %>%
    filter(IncludeInt == 1) %>% 
    group_by(ArticleTag) %>%
    summarize(PlausibleFPRate = mean(PlausibleInt))
  
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
makeTtestReliabilityTable <- function(studiesAnalysisDfPath = "data/StudiesAnalysis.csv",
                                      outputPath = "data/output/TtestFitTable.csv",
                                      diagnosticSavePath = "data/diagnostic/tTestFits.RDS",
                                      ntrials = 2, limit = 2)
{
  # Load analysis of studies from web interface stored as CSV.
  studiesDf <- read.csv(studiesAnalysisDfPath)

  # Filter out implausible false positives.
  studiesDf <- filter(studiesDf, (IncludeInt == 1) && (PlausibleInt == 1))
  
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
                           TrialIndex = c(), tTestPvalue = c(), 
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
                             TrialIndex = trialIdx, tTestPvalue = tTestResult$p.value,
                             ExpectedPower = expectedPower)
      
      resultsDf <- bind_rows(resultsDf, trial_row)
      
      # Create a test ID for this treatment and add the test output to diagnostic set.
      testId <- paste(row$ArticleTag, row$TreatmentTag, trialIdx, sep = "_")
      allTestsDiagnostic <- allTestsDiagnostic %>% add_stats(tTestResult, identifier = testId)
    }
  }
  
  write.csv(resultsDf, outputPath)
  saveRDS(allTestsDiagnostic, file = diagnosticSavePath)
  
  return (resultsDf)
}

doAnalyses <- function(studiesAnalysisDfPath = "data/StudiesAnalysis.csv")
{
  
}