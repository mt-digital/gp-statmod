##
# Experiments to understand the interpretation of existing evidence for 
# group polarization and to design new experiments.
#

# require(glue)
require(MASS)

source("JAGSModels.R")

## INITIALIZE JAGS PARAMETERS (Following Kruschke (2015) DBDA2E-utilities.R)

# Check that required packages are installed:
want = c("parallel","rjags","runjags","compute.es")
have = want %in% rownames(installed.packages())
if ( any(!have) ) { install.packages( want[!have] ) }

# Load rjags. Assumes JAGS is already installed.
try( library(rjags) )

# Load runjags. Assumes JAGS is already installed.
try( library(runjags) )
try( runjags.options( inits.warning=FALSE , rng.warning=FALSE ) )

# set default number of chains and parallelness for MCMC:
library(parallel) # for detectCores().
nCores = detectCores() 

if ( !is.finite(nCores) ) { nCores = 1 } 

if ( nCores > 4 ) 
{ 
  nChainsDefault = 4  # because JAGS has only 4 rng's.
  runjagsMethodDefault = "rjparallel"
} else if ( nCores == 4 ) { 
  nChainsDefault = 3  # save 1 core for other processes.
  runjagsMethodDefault = "rjparallel"
} else { 
  nChainsDefault = 3 
  runjagsMethodDefault = "rjags" # NOT parallel
}  

fileNameRoot = "JAGSOutput-"

calculateBayesian <- function(N, firstBinValue, nBins, latentMean, latentSd)
{
    data <- makeJAGSModelData(N, firstBinValue, nBins, latentMean, latentSd)

    parameters = c( "mu" , "sigma" , "thresh" )
    numSavedSteps = 20000
    thinSteps = 5 
    adaptSteps = 500  # Number of steps to "tune" the samplers
    burnInSteps = 1000
    saveName = fileNameRoot 
    runjagsMethod = runjagsMethodDefault
    nChains = nChainsDefault
    
    ordRunJagsOut <- run.jags(method="parallel", 
                              model="OrdAsOrdAndMet-OrdModel.txt", 
                              monitor=parameters, 
                              data=ordDataList,  
                              n.chains=nChains,
                              adapt=adaptSteps,
                              burnin=burnInSteps, 
                              sample=ceiling(numSavedSteps/nChains),
                              thin=thinSteps,
                              summarise=FALSE,
                              plots=FALSE)
    ordCodaSamples = as.mcmc.list(ordRunJagsOut)
    # resulting codaSamples object has these indices: 
    #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
    if (!is.null(saveName)) {
        save(ordCodaSamples , file=paste(saveName,"-Ord-Mcmc.Rdata",sep=""))
    }

    return (ordRunJagsOut)
}


tTestExperiment <- function(N, firstBinValue, nBins, latentMean, latentSdPre, 
                            latentSdPost, paired = TRUE, var.equal = FALSE)
{
    simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latentSdPre)
    # print(simObsPre)
    print(mean(simObsPre))
    # print("")
    simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latentSdPost)
    # print(simObsPost)
    print(mean(simObsPost))

    # df <- data.frame("pre" = simObsPre, "post" = simObsPost)

    return (t.test(simObsPre, simObsPost, paired = paired, var.equal = var.equal))
    # return (t.test(pre, post, data=df, paired = paired, var.equal = var.equal))
}


orderedProbitExperiment <- function(N, firstBinValue, nBins, latentMean,
                                        latentSdPre, latentSdPost, paired = TRUE,
                                        var.equal = FALSE)
{
    simObsPre <- simulatedObservation(N, firstBinValue, nBins, latentMean, latentSdPre)
    simObsPost <- simulatedObservation(N, firstBinValue, nBins, latentMean, latentSdPost)

    # Call helper function to run Bayesian analysis adapted from L&K 2018.
}


simulatedObservation <- function(N, firstBinValue, nBins, latentMean, latentSd)
{
    # Draw latent opinion data for each participant; parameters are set
    # to empirical, as opposed to population, values.
    
    # mvrnorm uses variance, not standard deviation, for its spread parameter.
    # latentVariance <- latentSd^2
    # latentData <- mvrnorm(N, latentMean, latentVariance, empirical=TRUE)
    latentData <- rnorm(N, latentMean, latentSd)

    # Create the bins and thresholds to be used to bin latentData.
    bins <- seq(from=firstBinValue, length.out=nBins)
    thresholds <- bins[1:length(bins)-1] + 0.5 

    # Bin latent data transforming to ordinal scale.
    ordinalData <- 0 * latentData
    for (sIdx in 1:length(ordinalData)) 
    {
        # Identify bin index from multiple comparisons of latent data with
        # threshold values. c(-Inf, thresholds) is a 1-D vector of 
        # -Inf followed by the threshold values. The `which` function gets 
        # the indexes where the greater than statement is TRUE. I.e., if the
        # opinion is not greater than any of the thresholds, and only -Inf,
        # then it will go in the first bin. If it is only greater than the
        # first threshold (second position in the -Inf, thresholds vector)
        # then the opinion goes in the second bin, and so on.
        binIdx <- max(which(latentData[sIdx] > c(-Inf, thresholds)))

        # Assign this subject's ordinal data to be the binned value.
        ordinalData[sIdx] <- bins[binIdx]
    }
    
    return (ordinalData)
}


##
# figureDir is where to save figures, if the argument is given.
#
simulatedObservations <- function(N = 10000, firstBinValue = 1, nBins = 10,
                                  latentMean = 6.4, preSD = 100, postSD = 1.0,
                                  figureDir = FALSE)
{
    pre <- simulatedObservation(N, firstBinValue, nBins, latentMean, preSD)
    post <- simulatedObservation(N, firstBinValue, nBins, latentMean, postSD)

    print(paste("Pre mean =", mean(pre), "; sd =", sd(pre)))
    print(paste("Post mean =", mean(post), "; sd =", sd(post)))

    if (is.character(figureDir))
    {
        breaks = (firstBinValue:(firstBinValue + nBins + 1)) - 0.5

        pdf(paste(figureDir, "/", "N=", N, "_FBV=", firstBinValue, "_latentMean=",
                  latentMean, "_SD=", preSD, ".pdf", sep=""))
        hist(pre, breaks = breaks, main = "Pre-deliberation opinions")
        dev.off()

        pdf(paste(figureDir, "/", "N=", N, "_FBV=", firstBinValue, "_latentMean=",
                  latentMean, "_SD=", postSD, ".pdf", sep=""))
        hist(post, breaks = breaks, main = "Post-deliberation opinions")
        dev.off()
    }

    return (cbind(pre, post));
}


##
# Initial sketch of data to be used to demonstrate the range of latent
# variances that can give rise to reported opinion shifts.
#
oneConditionExperiment <- 
    function(
        expectedPreObsMean, expectedPostObsMean, latentMean, 
        preDelibSDs = c(3.0, 3.5, 4.0, 4.5, 5.0), 
        postDelibSDs = c(1.0, 1.5, 2.0, 2.5, 3.0),
        nTrials = 10, 
        N = 30,
        tol = 0.05,
        firstBinValue = 1,
        nBins = 10
    )
{
    # Convenience function for generating simulated data for this 
    # experiment.
    getSimData <- function(latentSD)
    {
        return (simulatedObservation(N, firstBinValue, nBins, latentMean, latentSD))
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


approxEqual <- function(val1, val2, tol = 0.05)
{
    return (abs(val1 - val2) <= tol)
}


# Calculate the ranges of latent means that successfully generate a false
# detection of group polarization for several case studies.
#
# Arguments:
#     caseStudyDataFile (string): location of case studies CSV
#     step (numeric): Resolution of successful region over which solutions for
#                     latent standard deviations could be found.
latentMeanRanges <- function(caseStudyDataFile = "caseStudies.csv",
                             step = 0.1, guess = 1.0)
{
    csData <- read.csv(caseStudyDataFile)
    names(csData) <- csData$CaseStudyTag

    nCS <- nrow(csData)

    # Record the latent means that successfully generated a false
    # group polarization detection, plus the one that resulted in least error.
    csData$minLatentMean <- numeric(nCS)
    csData$maxLatentMean <- numeric(nCS)
    csData$bestLatentMean <- numeric(nCS)

    for (case in csData)
    {
        # If no solutions are found for a given case study 
        # report the means as NULL.
        if (is.null(knownSuccessMean))
        {
            case$minLatentMean <- NULL
            case$maxLatentMean <- NULL
            case$bestLatentMean <- NULL
        }
        else
        {
            resUp <- stepThroughSolutions(case, "up")
            resDown <- stepThroughSolutions(case, "down")

            latentMeans <- append(resUp[1], resDown[1])
            case$minLatentMean <- min(latentMeans)
            case$maxLatentMean <- max(latentMeans)

            case$bestLatentMean = min(resUp[2], resDown[2])

            case$bestSSE = min(resUp[3], resDown[3])
        }

        # Replace existing row with processed data.
        csData[case$CaseStudyTag, ] <- case
    }

    # Save processed dataframe to disk.
    write.csv(csData, "caseStudiesProcessed.csv")
}


##
# Arguments:
#     case (data.frame row): Single experimental condition case from case
#                            studies.
#
# stepThroughSolutions <- function(case, stepDir = "up", step = step, guess = guess)
# {
#     latentMean <- case$KnownSuccessfulMean
#     bestMean <- latentMean

#     kVec <- case$minBinValue:case$maxBinValue

#     # Solve for pre- and post-SDs using hillclimbing.
#     latentPreSDResult <- solveForLatentSD(
#         kVec, latentMean, cs$ObservedMeanPre, guess
#     )
#     latentPostSDResult <- solveForLatentSD(
#         kVec, latentMean, cs$ObservedMeanPost, guess
#     )

#     # Step latent mean up from known solution, tracking best solution.
#     successful <- TRUE
#     # Want to make sure the first SSE is less than the init val, so set to Inf.
#     currentSSE <- Inf
#     while (successful) 
#     {
#         # Append successful latent mean to list.
#         successfulMeans <- append(successfulMeans, latentMean)

#         # See if the SSE of this solution is 
#         currentSSE <- SSE(latentPreSDResult, latentPostSDResult)
#         if (currentSSE < bestSSE)
#         {
#             bestSSE <- currentSSE
#             bestMean <- latentMean
#         }

#         # Increment latentMean for the next loop if it happens.
#         if (stepDir == "up")
#             latentMean <- latentMean + step
#         else if (stepDir == "down")
#             latentMean <- latentMean - step
#         else
#             stop("Provided stepDir must be 'up' or 'down'")

#         # Solve for pre- and post-SDs using hillclimbing.
#         latentPreSDResult <- solveForLatentSD(
#             kVec, latentMean, cs$ObservedMeanPre, guess
#         )
#         latentPostSDResult <- solveForLatentSD(
#             kVec, latentMean, cs$ObservedMeanPost, guess
#         )
        
#         # See if results are within tol.
#         successful <- hillclimbSuccess(latentPreSDResult, 
#                                        latentPostSDResult)
#     }

#     return (c(successfulMeans, bestMean, bestSSE))
# }


function(latentPreSDResult, latentPostSDResult)
{
    return (latentPreSDResult[2] + latentPostSDResult[2])
}


hillclimbSuccess <- function(latentPreSDResult, latentPostSDResult, tol = 1e-3) 
{
    # The second element in the hillclimb results is the squared error.
    return (latentPreSDResult[2] < tol && latentPostSDResult[2] < tol)
}
