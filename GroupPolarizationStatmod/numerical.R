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
                             step = 0.01, maxIts = 10000, tol = 1e-6)
{
    # Keep the squared error function to be minimized through hillclimbing
    # separate from hillclimbing itself.
    sqErr <- function(latentSD)
    {
        simulatedObservedMean <- meanObs(
            kVec,
            makeProbVec(kVec, latentMean, latentSD)
        )

        return ((observedMean - simulatedObservedMean)^2);
    }

    return (hillclimbing(sqErr, guess, step, maxIts, tol))
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
    totalIts <- 0
    change <- 0.0

    while(its < maxIts)
    {
        its <- its + 1
        totalIts <- its

        xnext <- rnorm(1, xcurr, stepSize)
        ynext <- f(xnext)

        change <- ynext - ycurr

        if (abs(change) < tol)
        {
            return (c(xcurr, ycurr, totalIts));
        }
        
        # A negative change means ynext is less that ycurr, and closer/better
        # to minimum solution.
        if (change < 0)
        {
            xcurr <- xnext
            ycurr <- ynext
        }
    }

    return (c(xcurr, ycurr, totalIts));
}
