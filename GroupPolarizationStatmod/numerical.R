##
# This module implements the numerical procedure to find what latent standard
# deviation results in the observed mean in the large N limit, given latent
# mean µ, the bin vectors, and thresholds.
#
# This estimation from the large N limit can then be compared with finite N
# simulations with N specified by referring to studies we review in our
# paper. 
#
# library(cmna)
source('model.R')


solveForLatentSD <- function(kVec, latentMean, observedMean, guess, 
                             step = 0.01, maxIts = 10000, tol = 1e-6,
                             verbose = FALSE)
{
    sqErr <- function(latentSD)
    {
        simulatedObservedMean <- meanObs(
            kVec,
            makeProbVec(kVec, latentMean, latentSD)
        )

        return ((observedMean - simulatedObservedMean)^2);
    }

    ret <- hillclimbing(sqErr, guess, step, maxIts, tol)

    if (verbose)
        return (ret)
    else
        return (ret[1])
}


##
# Called hill climbing, but it's really hill tumbling to find the minimum.
#
hillclimbing <- function(f, x, stepSize = 0.1, maxIts = 1e5, tol = 1e-4)
{
    n <- length(x)

    xcurr <- x
    ycurr <- f(x)
    
    its <- 0
    change <- 0.0

    while(its < maxIts)
    {
        xnext <- rnorm(1, xcurr, stepSize * sqrt(0.5 * ycurr))
        ynext <- f(xnext)

        change <- ynext - ycurr
        if (abs(change) < tol)
        {
            return (c(xcurr, ycurr, its));
        }
        
        # A negative change means ynext is less that ycurr, and closer/better
        # to minimum solution.
        if (change < 0)
        {
            xcurr <- xnext
            ycurr <- ynext
        }
    }

    return (c(xcurr, ycurr, its));
}
