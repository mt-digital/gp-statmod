

##
# Large-N exact calculations.
##

# Probability that a participant reports opinion in bin k, i.e., between
# theta_{k-1} and theta_k.
binProb <- function(k, L, K, latentMean, latentSD)
{
    # Thresholds in theta vector start at L + 0.5 (see Eq. 6 in L&K supplement).
    thetaVec <- c(L:(L + K - 1) + 0.5)

    # Get this by solving k = (L - 1) + binIdx, where binIdx starts from 1.
    binIdx <- k - L + 1

    # But it only works if k is not the first...
    if (k == L)
    {
        theta_km1 = -Inf
    }
    else
    {
        theta_km1 <- thetaVec[binIdx - 1]
    }

    # ... or the last values.
    if (k == L + K - 1)
    {
        theta_k = Inf
    }
    else
    {
        theta_k <- thetaVec[binIdx]
    }
    print(theta_k)
    print(latentMean)
    print(latentSD)
    print("")
    ret <- pnorm((theta_k - latentMean) / latentSD) - pnorm((theta_km1 - latentMean) / latentSD)
    print("binProb ret:")
    print(ret)
    print("")
    return (
        pnorm((theta_k - latentMean) / latentSD) - 
        pnorm((theta_km1 - latentMean) / latentSD)
    );
}


makeProbVec <- function(kVec, latentMean, latentSD)
{
    L <- kVec[1]
    print(kVec)
    K <- length(kVec)
    print(K)

    probVec <- numeric(K)
    # for (k in kVec)
    for (ii in 1:K)
    {
        print(ii)
        res <- binProb(kVec[ii], L, K, latentMean, latentSD)
        print(res)
        probVec[ii] <- binProb(kVec[ii], L, K, latentMean, latentSD)
        # probVec <- append(probVec, c(binProb(k, L, K, latentMean, latentSD)))
        # print(probVec)
    }

    # print("in model:")
    # print(probVec)

    return (probVec);
}


meanObs <- function(kVec, probVec)
{
    return (sum(kVec * probVec));
}


sdObs <- function(kVec, probVec) {
    expectation <- meanObs(kVec, probVec)
    # print(expectation)
    variVec <- (kVec - expectation)^2
    # print(variVec)

    return ( sum(variVec * probVec)^(0.5) );
}
