##
# Experiments to understand the interpretation of existing evidence for 
# group polarization and to design new experiments.
#

require(glue)

source('plot.R')
source('model.R')

##
# Experiment modeling experimental procedure used by 
# Schkade, Sunstein, and Hastie (2010)
# "When deliberation produces extremism".
#
schkade2010 <- function(nParticipants=63,  # 34 women
                        nConservativeGroups=5, nDemocratGroups=5,
                        minInGroup=5, maxInGroup=7, 
                        nBins=10, 
                        nQuestions=3,  # gbl warming, affirm action, civil unions
                        preLibSd=2.3^2,  # Multiply by 2 since given in terms
                        postLibSd=1.0^2, # of standard deviation in 
                        preConsSd=2.8^2, # Schkade, et al.
                        postConsSd=0.75^2,
                        plot=FALSE,
                        saveFig=FALSE
)
{
    # Create groups according to paper. See p. 229 for start
    # of "Procedures and Results of Study.

    # Need to triple check that, indeed, the mean across all groups is
    # being used. Here I just choose one super-group (conservative/liberal)
    # to have 31 members and the other to have 32. 
    nLiberals <- 31
    nConservatives <- 32
    
    # Participants gave their opinions on three political issues.
    questions <- c('Global warming', 'Affirmative action', 'Civil unions')

    # Schkade, et al., claim to have found significant effects. We want to know
    # if these are potentially spurious, i.e., there is no significant effect
    # of deliberation. These means I assign below are close to the initial
    # pre-value given in Table 1.
    libMeans = c(
        'Global warming' = 9.2,
        'Affirmative action' = 5.8,
        'Civil unions' = 9.2
    )
    consMeans = c(
        'Global warming' = 5.13,
        'Affirmative action' = 2.8,
        'Civil unions' = 2.4
    )
    
    # Schkade, et al., report variances in Figure 2. I selected two points from
    # the scatterplot in Fig 2 to get the pre- and post-discussion variances
    # for defaults for the liberals and conservatives. 
    # A future sensitivity check should examine different pre-and post-
    # discussion variances selected from that figure.
    #
    # With the method explained, generate simulated data.
    # Start for now with just civil unions since this is the most
    # polarized issue.
    libMean <- libMeans[['Civil unions']]

    # See 
    libPreData <- simulatedObservation(nLiberals, 1, nBins, 
                                       libMean, preLibSd)
    libPostData <- simulatedObservation(nLiberals, 1, nBins, 
                                        libMean, postLibSd)
     
    consMean <- consMeans[['Civil unions']]
    consPreData <- simulatedObservation(nConservatives, 1, 10, 
                                        consMean, preConsSd)

    consPostData <- simulatedObservation(nConservatives, 1, 10, 
                                         consMean, postConsSd)

    # Run model on each simulated dataset. 
    libData <- makeModelInput(libPreData, libPostData)
    consData <- makeModelInput(consPreData, consPostData)

    libFitted <- frequentistModel(libData)
    consFitted <- frequentistModel(consData)

    if (plot)
    {
        for (party in c("Liberal", "Conservative"))
        {
            # Each result is a list with the simulated "input data" and
            # frequentist model fit results. 
            if (party == "Liberal")
            {
                plotFreq(libData, libFitted, nBins)
                latentMean = libMean
                initialLatentSd = preLibSd
                finalLatentSd = postLibSd
            }
            else
            {
                plotFreq(consData, consFitted, nBins)
                latentMean = consMean
                initialLatentSd = preConsSd
                finalLatentSd = postConsSd
            }
            
            if (saveFig)
            {
                fileName <- glue(
                    '~/workspace/Presentations/gp-statmod/Figures/nBins={nBins}',
                    '_party={party}',
                    '_lm1={latentMean}_lm2={latentMean}',
                    '_sd1={initialLatentSd}',
                    '_sd2={finalLatentSd}')
                
                saveGraph(fileName)
            }
        }
    }

    return (c(libData=libData, consData=consData, 
             libFitted=libFitted, consFitted=consFitted))
}


# Use this function to better understand how false negatives occur.
compareNoMeanChange <- function(nBins=10, firstBinValue=1, 
                                latentMeans=c(5.0, 7.0, 9.0),
                                initialLatentSd=4.0, finalLatentSd=1.0, N=100,
                                plot=TRUE, saveFig=FALSE)
{
    results <- c()
    for (latentMean in latentMeans)
    {
        result <- twoGroupComparison(N, firstBinValue, nBins, 
                                     latentMeans=c(latentMean, latentMean),
                                     latentSds=c(initialLatentSd, 
                                                 finalLatentSd),
        ) # rngSeed=42
        
        results <- append(results, result)
        
        if (plot)
        {
            # Each result is a list with the simulated "input data" and
            # frequentist model fit results. 
            plotFreq(result$inputData, result$fittedFreq, nBins)
            
            if (saveFig)
            {
                fileName <- glue('report/Figures/nBins={nBins}',
                                 '_lm1={latentMean}_lm2={latentMean}',
                                 '_sd1={initialLatentSd}',
                                 '_sd2={finalLatentSd}')
                
                saveGraph(fileName)
            }
        }
    }
    
    return (results)
}


# Use this function to better understand how false negatives occur.
compareChangedMean <- function(nBins=10, initialLatentMeans=c(2.5, 3.5, 4.5),
                               meanChange=1.0, initialLatentSd=1.0, 
                               finalLatentSd=4.0, firstBinValue=1, N=100,
                               plot=TRUE, saveFig=FALSE)
{
    results <- c()
    for (latentMean in initialLatentMeans)
    {
        initialLatentMean <- latentMean
        finalLatentMean <- latentMean + meanChange
        result <- twoGroupComparison(N, firstBinValue, nBins, 
                                     latentMeans=c(initialLatentMean, 
                                                   finalLatentMean),
                                     latentSds=c(initialLatentSd, 
                                                 finalLatentSd),
                                     )
        
        results <- append(results, result)
        
        if (plot)
        {
            # Each result is a list with the simulated "input data" and
            # frequentist model fit results. 
            plotFreq(result$inputData, result$fittedFreq, nBins)
            
            if (saveFig)
            {
                fileName <- glue('report/Figures/nBins={nBins}',
                                 '_lm1={initialLatentMean}',
                                 '_lm2={finalLatentMean}',
                                 '_sd1={initialLatentSd}',
                                 '_sd2={finalLatentSd}')
                
                saveGraph(fileName)
            }
        }
    }
    
    return (results)
}
