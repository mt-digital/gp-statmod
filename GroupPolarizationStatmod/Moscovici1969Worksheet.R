source("model.R")
source("experiments.R")

# Number of responses is the number per group (4) x the number of groups 
# (10 - DGaul or 5 - US) x the number of items (12 - DG, 11 - US).
nResponsesDeGaulle <- 4 * 10 * 12
nResponsesAmericans <- 4 * 5 * 11

# See Table 4, p. 132.
deGaullePreMean <- 0.9
deGaullePostMean <- 1.19

americansPreMean <- -0.61
americansPostMean <- -1.04

####   EXAMPLE SUCCESSFUL PARAMETER COMBINATIONS FOR BOTH CASES   ####

###  DeGaulle  ###
latentMean <- 1.2
preSD <- 6.0
postSD <- 0.15

dgPreData <- simulatedObservation(nResponsesDeGaulle, -3, 7, latentMean, preSD)
dgPostData <- simulatedObservation(nResponsesDeGaulle, -3, 7, latentMean, postSD)

### Americans  ###
latentMean <- -1.05
preSD <- 12.0
postSD <- 0.85

usPreData <- simulatedObservation(nResponsesAmericans, -3, 7, latentMean, preSD)
usPostData <- simulatedObservation(nResponsesAmericans, -3, 7, latentMean, postSD)


####   NOW TO GENERATE FALSE POSITIVES FOR THE PAPER   ####
moscoviciFalsePositives <- function(nTrials = 10, tol = 0.05)
{
    deGaullePreSDs <- c(5.0, 5.5, 6.0, 6.5, 7.0)
    deGaullePostSDs <- c(0.05, 0.10, 0.15, 0.2, 0.25)

    americansPreSDs <- c(10.0, 11.0, 12.0, 13.0, 14.0)
    americansPostSDs <- c(0.75, 0.8, 0.85, 0.9, 0.95)
    
    for (latentMean in c(1.1, 1.2, 1.3))
    {
        res <- oneConditionExperiment(deGaullePreMean, deGaullePostMean,
                                      latentMean,
                                      preDelibSDs = deGaullePreSDs,
                                      postDelibSDs = deGaullePostSDs,
                                      nTrials = nTrials, tol = tol,
                                      firstBinValue = -3, nBins = 7)

        write.csv(res,
                  paste("data/MoscoviciDG_latentMean=", latentMean,
                        "_tol=", tol, ".csv", sep = "")
                  )
        print(paste("Finished with latentMean =", latentMean))
    }

    for (latentMean in c(-1.0, -1.05, -1.1))
    {
        res <- oneConditionExperiment(americansPreMean, americansPostMean,
                                      latentMean, 
                                      preDelibSDs = americansPreSDs,
                                      postDelibSDs = americansPostSDs,
                                      nTrials = nTrials, tol = tol,
                                      firstBinValue = -3, nBins = 7)

        write.csv(res,
                  paste("data/MoscoviciAmericans_latentMean=", latentMean,
                        "_tol=", tol, ".csv", sep = "")
                  )

        print(paste("Finished with latentMean =", latentMean))
    }
}
