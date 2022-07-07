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
    errfunc <- function(latentSD)
    {
        simulatedObservedMean <- meanObs(
            kVec,
            makeProbVec(kVec, latentMean, latentSD)
        )
        
        # return ((observedMean - simulatedObservedMean)^2);
        return (simulatedObservedMean - observedMean);
    }

    return (hillclimbing(errfunc, guess, step, maxIts, tol))
}


##
# Called hill climbing, but it's really hill tumbling to find the minimum.
#
hillclimbing <- function(errfunc, sd_guess, stepSize = 0.1, maxIts = 1e5, tol = 1e-4)
{
    sd_curr <- sd_guess
    # sqerr_curr <- errfunc(sd_curr)[2]
    # err_curr <-  errfunc(sd_curr)[2]
    # sqerr_curr <- err_curr^2
    err_curr <- errfunc(sd_curr)
    sqerr_curr <- err_curr^2
    
    its <- 0
    change <- 0.0

    while(its < maxIts)
    {
        its <- its + 1

        # Undirected, random search for a better standard deviation to 
        # generate pre or post observed mean.
        sd_next <- rnorm(1, sd_curr, stepSize)
        err_next <- errfunc(sd_next)
        sqerr_next <- err_next^2

        change <- sqerr_next - sqerr_curr

        if (abs(change) < tol)
        {
            return (c(sd_curr, err_curr, its));
        }
        
        # A negative change means err_next is less that sqerr_curr, and closer/better
        # to minimum solution.
        if (change < 0)
        {
            sd_curr <- sd_next
            err_curr <- err_next
            sqerr_curr <- sqerr_next
        }
    }
    return (c(sd_curr, err_curr, its));
}
