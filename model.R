##
# model.R
#
# Refactored version of code from OSF page for Liddell & Krushke (2018).
# Working on customizing for group polarization modeling.
#
# Author: Matthew A. Turner
# Date: 2020-09-06
# Based on Liddell & Kruschke code available at https://osf.io/53ce9/
#
# See also Liddell, T. M., & Kruschke, J. K. (2018). Analyzing ordinal data with 
# metric models: What could possibly go wrong? Journal of Experimental 
# Social Psychology, 79. https://doi.org/10.1016/j.jesp.2018.08.009

require(MASS)

source("DBDA2E-utilities.R")


##
# Create a two-group model with two means and standard deviations, fit a 
# metric and ordinal model to each, and compare the measured values of 
# Cohen's d. I think of group 1 as the pre-discussion group and 
# group 2 as the post-discussion group.
#
# Arguments:
#     N (int): number of participants
#     firstBinValue (int): e.g. -3 on a seven-point Likert scale, or 
#         1 in Schkade, Sunstein, & Hastie (2007/2010)
#     nBins (int): number of ordinal bins to use
#     latentMeans (float, float): latent means of group 1 and group 2
#     latentSds (float, float): latent std devs of group 1 and group 2
# 
twoGroupComparison = function(N, firstBinValue, nBins, 
                              latentMeans=c(5, 5), latentSds=c(4, 1),
                              rngSeed=FALSE)
{
    if (rngSeed)
    {
        set.seed(rngSeed)
    }
    observation1 = simulatedObservation(N, firstBinValue, nBins, 
                                        latentMeans[1], latentSds[1])

    observation2 = simulatedObservation(N, firstBinValue, nBins,
                                        latentMeans[2], latentSds[2])

    ## 
    # Run metric and ordinal JAGS analysis on simulations. 
    inputData = makeModelInput(observation1, observation2)
    jagsData = assembleJAGSData(inputData, nBins)

    # Metric model fitting with t-test and frequentist hypothesis testing.
    fittedFrequentistModel = frequentistModel(inputData)
    # Metric model fitting with JAGS and Bayesian significance measures.
    fittedMetricModel = metricModel(jagsData$metDataList)
    # Ordinal model fitting with JAGS and Bayesian significance measures.
    # fittedOrdinalModel = ordinalModel(jagsData$ordDataList)

    # plotComparison(inputData, fittedFrequentistModel,
    #                fittedMetricModel, fittedOrdinalModel)

    return (list(inputData=inputData,
                 fittedFreq=fittedFrequentistModel,
                 fittedMetric=fittedMetricModel))
                 # fittedOrd=fittedOrdinalModel))
}


##
# Show the frequentist and Bayesian fit statistics and simulations overlaid
# on simulated data. Adapted from code generating Figure 2 in 
# L&K (2018).
# 
plotComparison = function(inputData, modelResults) 
{
    source('plot.R')
    plotFreq(modelResults$fittedFreq)
    plotMetric(modelResults$fittedFreq)
    plotOrdinal(modelResults$fittedFreq)
}



##
# Runs the metric frequentist (t-test) model on input dataframe.
#
frequentistModel = function(data) 
{
    preDiscussion = data$resp[data$cond == 1]
    postDiscussion = data$resp[data$cond == 2]

    tInfo = t.test(preDiscussion, postDiscussion)
    effSz = effectSize(preDiscussion, postDiscussion)

    return (list(tInfo=tInfo, effSz=effSz))
}

##
# Cohen's d.
#
effectSize = function(x, y)
{
    sdN = function(x) { sqrt(mean((x - mean(x))^2)) }

    # See p. 331 in L&K (2018).
    numer = (mean(x) - mean(y)) 
    denom = sqrt((sdN(x)^2 + sdN(y)^2)) / 2.0

    return (numer / denom)
}

##
# Runs the metric JAGS model (defined within) on input dataframe.
#
metricModel = function(metDataList) 
{
    # Define JAGS model. Here g is not "group" but instead "condition", as in
    # pre- and post-discussion opinion for the single group. Then pre- and
    # post-discussion each have their own mean and standard deviation, 
    # defined in the second for loop. The first for loop indicates that
    # individual opinions are drawn from a normal distribution with pre-
    # or post-discussion mean and standard deviation. 1/sigma[g[i]]^2 is
    # used because of the convention in Bayesian analysis to have the second
    # normal distribution parameter be precision, the reciprocal of variance,
    # instead of variance itself.
    metModelString = "
        model {
          for ( i in 1:Ntotal ) {
            meanY[i] ~ dnorm( mu[g[i]] , 1/sigma[g[i]]^2  )
          }
          for ( gIdx in 1:nG ) { 
            mu[gIdx] ~ dnorm( (1+max(nYlevels))/2 , 1/(max(nYlevels))^2 )
            sigma[gIdx] ~ dunif( 0.01 , max(nYlevels)*10 )
          }
        }
    "
    # Need to persist to disk, apparently.
    writeLines( metModelString , con="OrdAsOrdAndMet-MetModel.txt" )

    # RUN THE CHAINS FOR METRIC MODEL
    parameters = c( "mu" , "sigma" )
    numSavedSteps = 20000 # 20000 
    thinSteps = 5 # 5 
    adaptSteps = 500               # Number of steps to "tune" the samplers
    burnInSteps = 1000
    # saveName=fileNameRoot 
    runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
    nChains=nChainsDefault # from DBDA2E-utilities
    metRunJagsOut <- run.jags(method="parallel", # runjagsMethod ,
                              model="OrdAsOrdAndMet-MetModel.txt", 
                              monitor=parameters, 
                              data=metDataList,  
                              #inits=initsList, 
                              n.chains=nChains,
                              adapt=adaptSteps,
                              burnin=burnInSteps, 
                              sample=ceiling(numSavedSteps/nChains),
                              thin=thinSteps,
                              summarise=FALSE,
                              plots=FALSE)
    
    return (metRunJagsOut)
}


##
# Runs the ordinal JAGS model (defined within) on input dataframe.
#
ordinalModel = function(ordDataList) {

    ordModelString = "
        model {
            for ( i in 1:Ntotal ) {
              y[i] ~ dcat( pr[i,1:nYlevels[q[i]]] )
              pr[i,1] <- pnorm( thresh[q[i],1] ,
                                mu[g[i]] , 1/sigma[g[i]]^2 )
              for ( k in 2:(nYlevels[q[i]]-1) ) {
                pr[i,k] <- max( 0 ,  pnorm( thresh[q[i], k ]  ,
                                            mu[g[i]] , 1/sigma[g[i]]^2 )
                                   - pnorm( thresh[q[i],k-1] ,
                                            mu[g[i]] , 1/sigma[g[i]]^2 ) )
              }
              pr[i,nYlevels[q[i]]] <- (
                 1 - pnorm( thresh[q[i],nYlevels[q[i]]-1] ,
                            mu[g[i]] , 1/sigma[g[i]]^2 ) )
            }
            for ( gIdx in 1:nG ) { 
              mu[gIdx] ~ dnorm( (1+max(nYlevels))/2 , 1/(max(nYlevels))^2 )
              sigma[gIdx] ~ dunif( 0.01 , max(nYlevels)*10 )
            }
            # Prior on thresh[q,k]. Stochastic for all except thresh[1,1] and thresh[1,last].
            for ( qIdx in 1 ) {
              for ( kIdx in 2:(nYlevels[qIdx]-2) ) { # 1 and nYlevels-1 are fixed
                thresh[qIdx,kIdx] ~ dnorm( kIdx+0.5 , 1/2^2 )
              }
            }
            for ( qIdx in 2:nQ ) {
              for ( kIdx in 1:(nYlevels[qIdx]-1) ) { 
                thresh[qIdx,kIdx] ~ dnorm( kIdx+0.5 , 1/2^2 )
              }
            }
        }
      " # close quote for ordModelString

    writeLines( ordModelString , con="OrdAsOrdAndMet-OrdModel.txt" )
    #-----------------------------------------------------------------------------
    # INTIALIZE THE CHAINS.
    # Let JAGS do it...
    #-----------------------------------------------------------------------------
    # RUN THE CHAINS FOR ORDERED-PROBIT MODEL.
    parameters = c( "mu" , "sigma" , "thresh" )
    numSavedSteps = 20000 # 20000 
    thinSteps = 5 # 5 
    adaptSteps = 500               # Number of steps to "tune" the samplers
    burnInSteps = 1000
    # saveName=fileNameRoot 
    runjagsMethod=runjagsMethodDefault # from DBDA2E-utilities
    nChains=nChainsDefault # from DBDA2E-utilities

    ordRunJagsOut <- run.jags(method="parallel", # runjagsMethod ,
                              model="OrdAsOrdAndMet-OrdModel.txt", 
                              monitor=parameters, 
                              data=ordDataList,  
                              #inits=initsList, 
                              n.chains=nChains,
                              adapt=adaptSteps,
                              burnin=burnInSteps, 
                              sample=ceiling(numSavedSteps/nChains),
                              thin=thinSteps,
                              summarise=FALSE,
                              plots=FALSE)
    
    return (ordRunJagsOut)
}


##
# The JAGS program requires the data be in a particular format in a dataframe.
# This function combines two simulated observations into a properly formatted
# dataframe.
#
makeModelInput = function(observationData1, observationData2) 
{
    # Create columns for input dataframe as vectors. XXX Not sure which of 
    # these are actually required, just copying L&K for now.
    resp = item = subID = cond = 
        rep(0, prod(dim(observationData1)) + prod(dim(observationData2)))
    dfIdx = 1

    # Infer number of participants, XXX currently "subjects" which will change.
    N = length(observationData1)
    if (N != length(observationData2)) 
    { 
        stop("pre- and post-discussion observations not same length")
    }

    # Notes: 
    #   - item always 1 because for now just using one-question survey.
    #   - subID is the same in both for loops
    #   - cond indicates whether pre-discussion (1) or post-discussion (2)
    dfIdx = 1
    for (sIdx in 1:N)
    {
       resp[dfIdx] = observationData1[sIdx]
       item[dfIdx] = 1
       subID[dfIdx] = sIdx
       cond[dfIdx] = 1

       dfIdx = dfIdx + 1
    }
    for (sIdx in 1:N)
    {
       resp[dfIdx] = observationData2[sIdx]
       item[dfIdx] = 1
       subID[dfIdx] = sIdx
       cond[dfIdx] = 2

       dfIdx = dfIdx + 1
    }

    return (data.frame(resp=resp, item=item, subID=subID, cond=cond))
}


##
# JAGS apparently requires data be "assembled" a certain way according to the
# L&K code at l:669 in their OrdinalScaleGroupJags.R. 
#
assembleJAGSData = function(df, nBins)
{
    # XXX Most of the following is yanked from L&K, of questionable necessity.

    # Rename and reclass y values for convenience:
    y = as.numeric(df[, "resp"])
    # Do some checking that data make sense:
    if ( any( y!=round(y) ) ) { stop("All y values must be integers (whole numbers).") }
    if ( any( y < 1 ) ) { stop("All y values must be 1 or larger.") }

    totalNObs = length(y)

    qName = "item"
    sName = "subID"
    gName = "cond"

    # Question, subject, and group data vectors.
    q = as.numeric(as.factor(df[,qName]))
    s = as.numeric(as.factor(df[,sName])) # not used
    g = as.numeric(as.factor(df[,gName]))
    qLevels = levels(as.factor(df[,qName]))
    sLevels = levels(as.factor(df[,sName])) # not used
    gLevels = levels(as.factor(df[,gName]))
    nQ = max(q)
    nS = max(s) # not used
    nG = max(g)

    # Initialize all thresholds to be NA, which means they will be fit unless
    # specified.
    # thresh = rep(NA, nBins)
    thresh = matrix(NA, nrow=nQ, ncol=nBins - 1) # default to NA

    # Fix low and high thresholds. 
    thresh[1] = 1 + 0.5
    thresh[nBins - 1] = nBins - 0.5

    # Differs from L&K script because we have only one question, no need to
    # average over questions.
    metDataList = list(
        y=y, g=g, nYlevels=nBins,
        nG=nG, Ntotal=totalNObs
    )

    # Identical to L&K except list value variable names are my own.
    ordDataList = list(
        y=y, q=q, g=g, thresh=thresh, nYlevels=nBins,
        nQ=nQ, nG=nG, Ntotal=totalNObs
    )

    return (list(
        metDataList = metDataList,
        ordDataList = ordDataList
    ))
}


##
# Function for generating simulated data for a single-question group 
# polarization simulation. Thresholds will be evenly spaced by 1.
#
# In this preliminary implementation all groups have same stdandard deviation.
# TODO: add option for std deviation to be a random variable at the group 
# level. 
#
# Arguments:
#      N (int): number of participants
#      firstBinValue (int): e.g. -3 on a seven-point Likert scale, or 
#          1 in Schkade, Sunstein, & Hastie (2007/2010)
#      nBins (int): number of ordinal bins to use
#      latentMean (float): opinion on the continuous, latent psychological scale
#      latentStd (float): standard deviation of opinion, latent psychological scale
#      
simulatedObservation = function(N, firstBinValue, nBins, latentMean, latentStd)
{
    # Draw latent opinion data for each participant; parameters are set
    # to empirical, as opposed to population, values.
    latentData = mvrnorm(N, latentMean, latentStd, empirical=TRUE)

    # Create the bins and thresholds to be used to bin latentData.
    bins = seq(from=firstBinValue, length.out=nBins)
    thresholds = bins[1:length(bins)-1] + 0.5 

    # Bin latent data transforming to ordinal scale.
    ordinalData = 0 * latentData
    for (sIdx in 1:length(ordinalData)) 
    {
        # Identify bin index from multiple comparisons of latent data with
        # threshold values.
        binIdx = max(which(latentData[sIdx] > c(-Inf, thresholds)))

        # Assign this subject's ordinal data to be the binned value.
        ordinalData[sIdx] = bins[binIdx]
    }
    
    return (ordinalData)
}
