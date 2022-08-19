##
# Large-N exact calculations.
##

# Probability that a participant reports opinion in bin k, i.e., between
# theta_{k-1} and theta_k.
binProb <- function(thisBin, bin1, K, latentMean, latentSD)
{
    # Thresholds in theta vector start at bin1 + 0.5 (see Eq. 6 in L&K supplement).
    interiorThresholds <- c(bin1:(bin1 + K - 1) + 0.5)

    # Get this by solving thisbin = (bin1 - 1) + binIdx, where binIdx starts from 1.
    binIdx <- thisBin - bin1 + 1
    # But it only works if k is not the first...
    if (thisBin == bin1)
    {
        theta_km1 = -Inf
    }
    else
    {
        theta_km1 <- interiorThresholds[binIdx - 1]
    }
    # ... or the last values.
    if (thisBin == bin1 + K - 1)
    {
        theta_k = Inf
    }
    else
    {
        theta_k <- interiorThresholds[binIdx]
    }
    # From Equation 2 in Liddell & Kruschke (2018), p. 340.
    return (pnorm((theta_k - latentMean) / latentSD) - 
            pnorm((theta_km1 - latentMean) / latentSD))
}


makeProbVec <- function(kVec, latentMean, latentSD)
{
    L <- kVec[1]
    K <- length(kVec)

    probVec <- numeric(K)

    for (ii in 1:K)
    {
        res <- binProb(kVec[ii], L, K, latentMean, latentSD)
        probVec[ii] <- binProb(kVec[ii], L, K, latentMean, latentSD)
    }

    return (probVec);
}


# Simulated observed mean denoted <o_t> in paper.
meanObs <- function(kVec, probVec)
{
    return (sum(kVec * probVec));
}


# Simulated observed standard deviation denoted S_t in paper.
sdObs <- function(kVec, probVec) {
    expectation <- meanObs(kVec, probVec)
    variVec <- (kVec - expectation)^2

    return ( sum(variVec * probVec)^(0.5) );
}
