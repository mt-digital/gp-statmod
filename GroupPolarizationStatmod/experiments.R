##
# Experiments to understand the interpretation of existing evidence for 
# group polarization and to design new experiments.
#

# require(glue)
require(MASS)

# source('plot.R')
# source('model.R')


simulatedObservation <- function(N, firstBinValue, nBins, latentMean, latentStd)
{
    # Draw latent opinion data for each participant; parameters are set
    # to empirical, as opposed to population, values.
    
    # mvrnorm uses variance, not standard deviation, for its spread parameter.
    latentVariance <- latentStd^2
    latentData <- mvrnorm(N, latentMean, latentVariance, empirical=TRUE)

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


# ##
# # Experiment modeling experimental procedure used by 
# # Schkade, Sunstein, and Hastie (2010)
# # "When deliberation produces extremism".
# #
# schkade2010 <- function(nParticipants=63,  # 34 women
#                         nConservativeGroups=5, nDemocratGroups=5,
#                         minInGroup=5, maxInGroup=7, 
#                         nBins=10, 
#                         nQuestions=3,  # gbl warming, affirm action, civil unions
#                         preLibSd=2.3^2,  # Multiply by 2 since given in terms
#                         postLibSd=1.0^2, # of standard deviation in 
#                         preConsSd=2.8^2, # Schkade, et al.
#                         postConsSd=0.75^2,
#                         plot=FALSE,
#                         saveFig=FALSE
# )
# {
#     # Create groups according to paper. See p. 229 for start
#     # of "Procedures and Results of Study.

#     # Need to triple check that, indeed, the mean across all groups is
#     # being used. Here I just choose one super-group (conservative/liberal)
#     # to have 31 members and the other to have 32. 
#     nLiberals <- 31
#     nConservatives <- 32
    
#     # Participants gave their opinions on three political issues.
#     questions <- c('Global warming', 'Affirmative action', 'Civil unions')

#     # Schkade, et al., claim to have found significant effects. We want to know
#     # if these are potentially spurious, i.e., there is no significant effect
#     # of deliberation. These means I assign below are close to the initial
#     # pre-value given in Table 1.
#     libMeans = c(
#         'Global warming' = 9.2,
#         'Affirmative action' = 5.8,
#         'Civil unions' = 9.2
#     )
#     consMeans = c(
#         'Global warming' = 5.13,
#         'Affirmative action' = 2.8,
#         'Civil unions' = 2.4
#     )
    
#     # Schkade, et al., report variances in Figure 2. I selected two points from
#     # the scatterplot in Fig 2 to get the pre- and post-discussion variances
#     # for defaults for the liberals and conservatives. 
#     # A future sensitivity check should examine different pre-and post-
#     # discussion variances selected from that figure.
#     #
#     # With the method explained, generate simulated data.
#     # Start for now with just civil unions since this is the most
#     # polarized issue.
#     libMean <- libMeans[['Civil unions']]

#     # See 
#     libPreData <- simulatedObservation(nLiberals, 1, nBins, 
#                                        libMean, preLibSd)
#     libPostData <- simulatedObservation(nLiberals, 1, nBins, 
#                                         libMean, postLibSd)
     
#     consMean <- consMeans[['Civil unions']]
#     consPreData <- simulatedObservation(nConservatives, 1, 10, 
#                                         consMean, preConsSd)

#     consPostData <- simulatedObservation(nConservatives, 1, 10, 
#                                          consMean, postConsSd)

#     # Run model on each simulated dataset. 
#     libData <- makeModelInput(libPreData, libPostData)
#     consData <- makeModelInput(consPreData, consPostData)

#     libFitted <- frequentistModel(libData)
#     consFitted <- frequentistModel(consData)

#     if (plot)
#     {
#         for (party in c("Liberal", "Conservative"))
#         {
#             # Each result is a list with the simulated "input data" and
#             # frequentist model fit results. 
#             if (party == "Liberal")
#             {
#                 plotFreq(libData, libFitted, nBins)
#                 latentMean = libMean
#                 initialLatentSd = preLibSd
#                 finalLatentSd = postLibSd
#             }
#             else
#             {
#                 plotFreq(consData, consFitted, nBins)
#                 latentMean = consMean
#                 initialLatentSd = preConsSd
#                 finalLatentSd = postConsSd
#             }
            
#             if (saveFig)
#             {
#                 fileName <- glue(
#                     '~/workspace/Presentations/gp-statmod/Figures/nBins={nBins}',
#                     '_party={party}',
#                     '_lm1={latentMean}_lm2={latentMean}',
#                     '_sd1={initialLatentSd}',
#                     '_sd2={finalLatentSd}')
                
#                 saveGraph(fileName)
#             }
#         }
#     }

#     return (c(libData=libData, consData=consData, 
#              libFitted=libFitted, consFitted=consFitted))
# }


# # Use this function to better understand how false negatives occur.
# compareNoMeanChange <- function(nBins=10, firstBinValue=1, 
#                                 latentMeans=c(5.0, 7.0, 9.0),
#                                 initialLatentSd=4.0, finalLatentSd=1.0, N=100,
#                                 plot=TRUE, saveFig=FALSE)
# {
#     results <- c()
#     for (latentMean in latentMeans)
#     {
#         result <- twoGroupComparison(N, firstBinValue, nBins, 
#                                      latentMeans=c(latentMean, latentMean),
#                                      latentSds=c(initialLatentSd, 
#                                                  finalLatentSd),
#         ) # rngSeed=42
        
#         results <- append(results, result)
        
#         if (plot)
#         {
#             # Each result is a list with the simulated "input data" and
#             # frequentist model fit results. 
#             plotFreq(result$inputData, result$fittedFreq, nBins)
            
#             if (saveFig)
#             {
#                 fileName <- glue('report/Figures/nBins={nBins}',
#                                  '_lm1={latentMean}_lm2={latentMean}',
#                                  '_sd1={initialLatentSd}',
#                                  '_sd2={finalLatentSd}')
                
#                 saveGraph(fileName)
#             }
#         }
#     }
    
#     return (results)
# }


# # Use this function to better understand how false negatives occur.
# compareChangedMean <- function(nBins=10, initialLatentMeans=c(2.5, 3.5, 4.5),
#                                meanChange=1.0, initialLatentSd=1.0, 
#                                finalLatentSd=4.0, firstBinValue=1, N=100,
#                                plot=TRUE, saveFig=FALSE)
# {
#     results <- c()
#     for (latentMean in initialLatentMeans)
#     {
#         initialLatentMean <- latentMean
#         finalLatentMean <- latentMean + meanChange
#         result <- twoGroupComparison(N, firstBinValue, nBins, 
#                                      latentMeans=c(initialLatentMean, 
#                                                    finalLatentMean),
#                                      latentSds=c(initialLatentSd, 
#                                                  finalLatentSd),
#                                      )
        
#         results <- append(results, result)
        
#         if (plot)
#         {
#             # Each result is a list with the simulated "input data" and
#             # frequentist model fit results. 
#             plotFreq(result$inputData, result$fittedFreq, nBins)
            
#             if (saveFig)
#             {
#                 fileName <- glue('report/Figures/nBins={nBins}',
#                                  '_lm1={initialLatentMean}',
#                                  '_lm2={finalLatentMean}',
#                                  '_sd1={initialLatentSd}',
#                                  '_sd2={finalLatentSd}')
                
#                 saveGraph(fileName)
#             }
#         }
#     }
    
#     return (results)
# }

