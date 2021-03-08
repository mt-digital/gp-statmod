##
# Script to see if shifts in Schkade, Sunstein, and Hastie (2010) may plausibly
# be due to a combination of narrowing opinion variance and the use of an
# ordinal scale for measuring opinions.
#
# "Magic" values for latent means and std devs came from an informal 
# inverse search via trial and error in the R console. 

# Means for different groups/questions come from Table 1 (??) in Schkade;
# standard deviations are estimated by reading off the plot in Figure 2.
require(testit)
source("model.R")


# Reported liberal mean opinions.
globWarmLibPreMean <- 9.19; affirmActLibPreMean <- 5.81; civUnionsLibPreMean <- 9.22;
globWarmLibPostMean <- 9.44; affirmActLibPostMean <- 6.38; civUnionsLibPostMean <- 9.69;
combinedLibPreMean <- 8.07; combinedLibPostMean <- 8.50;

# Reported conservative mean opinions.
globWarmConsPreMean <- 5.13; affirmActConsPreMean <- 2.84; civUnionsConsPreMean <- 2.48;
globWarmConsPostMean <- 2.97; affirmActConsPostMean <- 1.61; civUnionsConsPostMean <- 2.19;
combinedConsPreMean <- 3.48; combinedConsPostMean <- 2.26;


# Runnable demonstration that the simulated values using the following latent
# means and standard deviations reproduce the observed shifts and estimates of
# variance. Use just one set of SDs for now, probably want to use 3-4 for
# sensitivity analyses later.

# Parameters for all cases for Schkade 2010.
N = 30;
latPreSd <- 6.0;
latPostSd <- 1.5;
firstBinValue <- 1;
nBins <- 10;


basicFalseEffectDemo <- function(verbose = TRUE)
{

    if (verbose) 
    {
        print("*****  BASIC DEMONSTRATION OF GENERATING FALSE GROUP POLARIZATION EFFECTS *****")
    }

    print("*** Global warming question --- Liberals (FALSE POSITIVE) ***")
    latentMean <- 9.9; 

    simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPreSd);
    while(abs(globWarmLibPreMean - mean(simObsPre)) > 0.05)
    {
        simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPreSd);
    }
    if (verbose)
    {
        print(paste("Simulated observed pre-deliberation mean: ", mean(simObsPre)))
        print(paste("Simulated observed pre-deliberation SD: ", sd(simObsPre)))
    }
    assert(abs(globWarmLibPreMean - mean(simObsPre)) <= 0.05)

    simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPostSd);
    # assert(FALSE)
    while(abs(globWarmLibPostMean - mean(simObsPost)) > 0.05)
    {
        simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPostSd);
    }
    if (verbose) 
    {
        print(paste("Simulated observed post-deliberation mean: ", mean(simObsPost)))
        print(paste("Simulated observed post-deliberation SD: ", sd(simObsPost)))
        print("These simulated observations match reported shifts, so this case is suspect for a false positive.")
        print("")
        print("Because it is a FP candidate, we run the t-test for detecting whether two distributions have identical means.")
        modelInput <- makeModelInput(simObsPre, simObsPost);
        modelFit <- frequentistModel(modelInput)
        print(modelFit)
    }
    assert(abs(globWarmLibPostMean - mean(simObsPost)) <= 0.05)


    print("*** Civil unions question --- Conservatives (FALSE POSITIVE) ***")
    latentMean <- 2.0; 

    simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPreSd);
    while(abs(civUnionsConsPreMean - mean(simObsPre)) > 0.05)
    {
        simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPreSd);
    }
    if (verbose)
    {
        print(paste("Simulated observed pre-deliberation mean: ", mean(simObsPre)))
        print(paste("Simulated observed pre-deliberation SD: ", sd(simObsPre)))
    }
    assert(abs(civUnionsConsPreMean - mean(simObsPre)) <= 0.05)

    simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPostSd);
    while(abs(civUnionsConsPostMean - mean(simObsPost)) > 0.05)
    {
        simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPostSd);
    }
    if (verbose)
    {
        print(paste("Simulated observed post-deliberation mean: ", mean(simObsPost)))
        print(paste("Simulated observed post-deliberation SD: ", sd(simObsPost)))
        print("These simulated observations match reported shifts, so this case is suspect for a false positive.")
        print("")
    }
    assert(abs(civUnionsConsPostMean - mean(simObsPost)) <= 0.05)


    print("*** Civil unions question --- Liberals (FALSE POSITIVE) ***")
    latentMean <- 10; 

    simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPreSd);
    while(abs(civUnionsLibPreMean - mean(simObsPre)) > 0.05)
    {
        simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPreSd);
    }
    if (verbose) 
    {
        print(paste("Simulated observed pre-deliberation mean: ", mean(simObsPre)))
        print(paste("Simulated observed pre-deliberation SD: ", sd(simObsPre)))
    }
    assert(abs(civUnionsLibPreMean - mean(simObsPre)) <= 0.05)

    simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPostSd);
    while(abs(civUnionsLibPostMean - mean(simObsPost)) > 0.05)
    {
        simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latPostSd);
    }

    if (verbose) 
    {
        print(paste("Simulated observed post-deliberation mean: ", mean(simObsPost)))
        print(paste("Simulated observed post-deliberation SD: ", sd(simObsPost)))
        print("These simulated observations match reported shifts, so this case is suspect for a false positive.")
        print("")
    }
    assert(abs(civUnionsLibPostMean - mean(simObsPost)) <= 0.05)
}


##
# Initial sketch of data to be used to demonstrate the range of latent
# variances that can give rise to reported opinion shifts.
#
oneConditionExperiment <- 
    function(
        expectedPreObsMean, expectedPostObsMean, latentMean, 
        # preDelibSDs = c(2.5, 3.0, 3.5, 4.0, 4.5), 
        preDelibSDs = c(3.0, 3.5, 4.0, 4.5, 5.0), 
        postDelibSDs = c(1.0, 1.5, 2.0, 2.5, 3.0),
        nTrials = 10, 
        N = 30,
        tol = 0.05
    )
{
    # Convenience function for generating simulated data for this 
    # experiment.
    getSimData <- function(latentSD)
    {
        return (simulatedObservation(N, 1, 10, latentMean, latentSD))
    }

    # Initialize results dataframe. No need to keep individual trials, only
    # need the success rate for a given parameter combination.
    result <- data.frame(matrix(ncol=4, 
                                nrow=length(preDelibSDs) * length(postDelibSDs)
                         ))

    colnames(result) <- c(
        'latentMean', 'preDelibSD', 'postDelibSD', 'successRate'
    )
    # Use this to index (preSD, postSD) pairs and add results to DF.
    paramSetIdx = 1
    for (preSD in preDelibSDs)
    {
        for (postSD in postDelibSDs)
        {
            trialResults <- c()
            for (trialIdx in seq(1, nTrials))
            {
                simObsPre <- getSimData(preSD)
                simObsPost <- getSimData(postSD)
                if (approxEqual(expectedPreObsMean, mean(simObsPre), tol = tol) &&
                    approxEqual(expectedPostObsMean, mean(simObsPost), tol = tol))
                {
                    success = 1
                }
                else
                {
                    success = 0
                }
                trialResults <- append(trialResults, success)
            }
            
            result[paramSetIdx, ] = c(latentMean, preSD, postSD, mean(trialResults))
            paramSetIdx = paramSetIdx + 1
        }
    }

    return (result)
}


##
# Run all conditions, saved to CSV files with names indicating the condition,
# e.g., BoulderGlobalWarming_mean=9.8.csv or COSpringsCivilUnions_mean=
allFalsePositiveConditions <- function(nTrials = 10, tol = 0.05)
{
    # Boulder, Global Warming and Civil Unions can use same sets of params.
    civilUnionPreSDs = c(3.5, 4.0, 4.5, 5.0, 5.5)
    civilUnionPostSDs = c(0.0, 0.5, 1.0, 1.5, 2.0)
    for (latentMean in c(9.8, 10.0, 10.2))
    {
        gwResult <- oneConditionExperiment(globWarmLibPreMean, 
                                           globWarmLibPostMean, latentMean,
                                           nTrials = nTrials, tol = tol)
        write.csv(gwResult, 
                  paste("data/BoulderGlobalWarming_latentMean=", latentMean, 
                        "_tol=", tol, ".csv", 
                        sep = ""))

        cuResult <- oneConditionExperiment(civUnionsLibPreMean,
                                           civUnionsLibPostMean, latentMean,
                                           preDelibSDs = civilUnionPreSDs,
                                           postDelibSDs = civilUnionPostSDs,
                                           nTrials = nTrials, tol = tol)
        write.csv(cuResult, 
                  paste("data/BoulderCivilUnions_latentMean=", latentMean,
                        "_tol=", tol, ".csv", 
                        sep = ""))
        print(paste("Finished latentMean =", latentMean))
    }

    # CO Springs, Civil Unions.
    preSDs = c(4.0, 4.5, 5.0, 5.5, 6.0)
    for (latentMean in c(1.9, 2.0, 2.1, 2.2))
    {
        cuResult <- oneConditionExperiment(civUnionsConsPreMean,
                                           civUnionsConsPostMean, latentMean,
                                           preDelibSDs = preSDs,
                                           nTrials = nTrials, tol = tol)
        write.csv(cuResult, 
                  paste("data/COSpringsCivilUnions_latentMean=", latentMean, 
                        "_tol=", tol, ".csv", 
                        sep = ""))

        print(paste("Finished latentMean =", latentMean))
    }
}



approxEqual <- function(val1, val2, tol = 0.05)
{
    return (abs(val1 - val2) <= tol)
}
